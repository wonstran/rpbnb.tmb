#' Fit a bivariate random-parameter negative binomial model (TMB)
#'
#' Maximum simulated likelihood via TMB (Template Model Builder) with
#' automatic differentiation. Supports Famoye/Sarmanov and discrete-copula
#' (Frank, Gaussian, Clayton) dependence.
#'
#' @param formula_1,formula_2 Model formulas for the two count outcomes.
#' @param data A data frame.
#' @param random_1,random_2 Random coefficient specs per equation. NULL (all
#'   fixed), character vector (all Normal), or named list with per-variable
#'   distribution specs (\code{"normal"}, \code{"lognormal"}, \code{"uniform"},
#'   \code{"triangular"}).
#' @param draws Number of Halton simulation draws.
#' @param seed Random seed for draws.
#' @param start Optional starting parameter vector (named or positional).
#' @param dependence Dependence structure: "famoye", "independence", or a
#'   \code{copula()} object for copula dependence.
#' @param control An object from \code{rpbnb_tmb_control()}.
#' @param inference Inference storage: \code{"full"} for a full covariance,
#'   \code{"diag"} for standard errors only, or \code{"none"} to skip Hessian
#'   calculations. In diagonal mode, \code{vcov()} returns \code{NA} for
#'   covariance cross-terms.
#' @param keep Retained fit state: \code{"postfit"} keeps data needed for
#'   marginal effects, \code{"compact"} keeps estimates only, and
#'   \code{"full"} also retains the TMB objective and responses. Low-level
#'   diagnostics that access \code{fit$obj} require \code{"full"}.
#' @param poisson_1,poisson_2 Fit the corresponding margin at its exact Poisson
#'   limit (NB2 dispersion m = 0, held fixed).
#' @return An object of class \code{rpbnb_tmb_fit}. The \code{sdreport} field
#'   is a compact package-owned summary and does not retain a second TMB tape.
#'
#'   \code{boundary_report} is a character vector naming the reported
#'   quantities whose estimates are pinned against an implementation bound --
#'   the frozen Famoye interval, the Frank overflow guard, or an \code{exp()}
#'   clamp -- or whose delta-method derivative has collapsed numerically. Those
#'   estimates are set by the implementation rather than identified by the
#'   data, so their standard errors and covariances are reported as \code{NA}
#'   and a warning is raised. An empty vector means no such constraint bound.
#'
#'   \code{boundary_sides} names, for each entry of \code{boundary_report},
#'   which end was reached: \code{"lower"}, \code{"upper"}, or
#'   \code{"degenerate"} when the derivative collapsed rather than a bound
#'   being hit. The two sides call for opposite remedies -- a lower dispersion
#'   clamp means the margin is effectively Poisson, an upper one means it is
#'   degenerately over-dispersed -- so this is retained on the fit rather than
#'   only mentioned in the warning.
#'
#'   \code{lambda_bounds} is a named numeric vector \code{c(lower =, upper =)}
#'   giving the admissible Famoye dependence interval, and \code{NULL} for
#'   every other dependence structure. The bounds are computed once at the
#'   starting values and held fixed for the whole fit; they are \emph{not}
#'   recomputed at the optimum. A \code{lam} estimate at either end is
#'   therefore an artefact of the starting values rather than a property of the
#'   data, which is why the field is exposed.
#'   \code{\link{rpbnb_tmb_dependence_profile}} reports a likelihood-based
#'   interval where the delta-method standard error collapses to \code{NA}, but
#'   for Famoye that interval is mapped through this same frozen box and so
#'   cannot escape it: widen the box by refitting from better starting values
#'   before treating such an interval as informative.
#' @export
#' @examples
#' \dontrun{
#' sim <- simulate_rpbnb_tmb(n = 300,
#'   beta1 = c("(Intercept)" = 0.2, x1 = 0.4),
#'   beta2 = c("(Intercept)" = 0.1, x1 = -0.3),
#'   dispersion = c(m1 = 0.4, m2 = 0.5), seed = 1)
#' fit <- fit_rpbnb_tmb(y1 ~ x1, y2 ~ x1, data = sim$data, draws = 100)
#' coef(fit)
#' }
fit_rpbnb_tmb <- function(formula_1, formula_2, data,
                          random_1 = NULL, random_2 = NULL,
                          draws = 400L, seed = 1234L, start = NULL,
                          dependence = "famoye",
                          control = rpbnb_tmb_control(),
                          inference = c("full", "diag", "none"),
                          keep = c("postfit", "compact", "full"),
                          poisson_1 = FALSE, poisson_2 = FALSE) {
  had_random_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_random_seed) {
    saved_random_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }
  on.exit({
    if (had_random_seed) {
      assign(".Random.seed", saved_random_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)

  stopifnot(is.data.frame(data))
  inference <- match.arg(inference)
  keep <- match.arg(keep)
  .chk_poisson_flag(poisson_1, "poisson_1")
  .chk_poisson_flag(poisson_2, "poisson_2")
  if (length(draws) != 1L || !is.numeric(draws) || is.na(draws) ||
      !is.finite(draws) || draws < 1 || draws != floor(draws) ||
      draws > .Machine$integer.max) {
    stop("draws must be one whole number greater than or equal to 1.",
         call. = FALSE)
  }
  draws <- as.integer(draws)

  # Resolve dependence
  family_code <- .resolve_family_code(dependence)
  # `copula(par = )` is a simulation argument. The fitting path reads only the
  # family, and the dependence parameter is always estimated -- so accepting a
  # `par` here would validate it, ignore it, and return a fit indistinguishable
  # from one that never supplied it. Reject it instead of silently discarding
  # a value the caller had every reason to think was doing something.
  if (inherits(dependence, "rpbnb_copula") && !is.null(dependence$par)) {
    stop("`copula(par = )` is a simulation argument and has no effect on ",
         "fitting: the dependence parameter is always estimated. Drop `par` ",
         "to fit, or use `start` to set its working-scale starting value.",
         call. = FALSE)
  }
  # Codes 1-3 are exactly the copula families.
  if (family_code >= 1L && (isTRUE(poisson_1) || isTRUE(poisson_2))) {
    stop("Poisson-limit margins not supported with copula dependence.")
  }

  # Prepare data
  prep <- .prepare_bnb_data(formula_1, formula_2, data)
  Y1 <- prep$Y1; Y2 <- prep$Y2
  X1 <- prep$X1; X2 <- prep$X2
  n <- prep$n

  # Parse random specs
  spec1 <- parse_rand_spec(random_1)
  spec2 <- parse_rand_spec(random_2)
  rand_idx1 <- match(spec1$names, colnames(X1))
  rand_idx2 <- match(spec2$names, colnames(X2))
  if (any(is.na(rand_idx1))) stop("Random names not found in formula_1.")
  if (any(is.na(rand_idx2))) stop("Random names not found in formula_2.")
  dist1 <- match(spec1$dist, c("normal", "lognormal", "uniform", "triangular")) - 1L
  dist2 <- match(spec2$dist, c("normal", "lognormal", "uniform", "triangular")) - 1L
  sign1 <- as.integer(spec1$sign)
  sign2 <- as.integer(spec2$sign)
  q1 <- length(rand_idx1); q2 <- length(rand_idx2)
  k1 <- ncol(X1); k2 <- ncol(X2)
  total_rand <- q1 + q2
  effective_draws <- if (total_rand > 0L) draws else 1L
  configured_threads <- .configure_tmb_threads(
    control$n_cores,
    max_threads = control$max_threads,
    parallel_tape = control$parallel_tape
  )
  .check_tmb_workload(
    n = n,
    draws = effective_draws,
    family_code = family_code,
    max_workload = control$max_workload,
    n_threads = configured_threads,
    parallel_tape = control$parallel_tape
  )

  # Generate Halton draws
  set.seed(seed)
  if (total_rand > 0) {
    Z <- halton_uniform(draws, total_rand, burn = control$halton_burn)
    Z1 <- if (q1 > 0) Z[, 1:q1, drop = FALSE] else matrix(0, draws, 0)
    Z2 <- if (q2 > 0) Z[, (q1 + 1):total_rand, drop = FALSE] else matrix(0, draws, 0)
  } else {
    Z1 <- matrix(0, 1, 0); Z2 <- matrix(0, 1, 0)
  }

  # Parameter names
  scale_lab <- function(dist, cols) {
    vapply(seq_along(dist), function(j) {
      paste0(rand_dist_registry[[dist[j]]]$scale_label, cols[j])
    }, character(1))
  }
  par_names <- c(paste0("b1:", colnames(X1)), paste0("b2:", colnames(X2)),
                 if (q1 > 0) scale_lab(spec1$dist, paste0("1:", colnames(X1)[rand_idx1])),
                 if (q2 > 0) scale_lab(spec2$dist, paste0("2:", colnames(X2)[rand_idx2])),
                 "log_m1", "log_m2",
                 if (family_code >= 0L) "z_dep" else NULL)
  # Default start values
  default_start <- c(rep(0, k1 + k2),
                     if (q1 > 0) rep(log(0.2), q1),
                     if (q2 > 0) rep(log(0.2), q2),
                     log(0.5), log(0.5),
                     if (family_code >= 0L) {
                       # Starting Frank exactly at independence causes its
                       # removable singularity to erase the dependence score
                       # from the recorded AD tape.
                       if (family_code == 1L) 0.1 else 0
                     } else NULL)
  start <- .resolve_start(start, default_start, par_names)
  names(start) <- par_names

  # Poisson margins: pin log_m at log(POISSON_M)
  fixed_names <- c(if (isTRUE(poisson_1)) "log_m1",
                   if (isTRUE(poisson_2)) "log_m2")
  if (length(fixed_names)) {
    start[fixed_names] <- log(POISSON_M)
  }
  free <- !(par_names %in% fixed_names)
  npar <- sum(free)

  # Famoye: compute frozen lambda bounds at start
  lamLo <- 0; lamHi <- 0
  if (family_code == 0L) {
    # Extract beta/sd/dispersion from start
    i1 <- 1:k1; i2 <- (k1 + 1):(k1 + k2)
    s1 <- if (q1 > 0) exp(start[(k1 + k2 + 1):(k1 + k2 + q1)]) else numeric(0)
    s2 <- if (q2 > 0) exp(start[(k1 + k2 + q1 + 1):(k1 + k2 + q1 + q2)]) else numeric(0)
    idx_end <- k1 + k2 + q1 + q2
    m1_start <- if (poisson_1) 0 else exp(start[idx_end + 1])
    m2_start <- if (poisson_2) 0 else exp(start[idx_end + 2])
    xb1 <- as.vector(X1 %*% start[i1]); xb2 <- as.vector(X2 %*% start[i2])
    Rloc <- if (total_rand > 0) draws else 1L
    lamLo <- -Inf; lamHi <- Inf
    for (r in seq_len(Rloc)) {
      dev1 <- if (q1 > 0) rand_realize(matrix(Z1[r, ], 1, q1), spec1$dist, sign1, start[i1][rand_idx1], s1)$dev[1, ] else numeric(0)
      dev2 <- if (q2 > 0) rand_realize(matrix(Z2[r, ], 1, q2), spec2$dist, sign2, start[i2][rand_idx2], s2)$dev[1, ] else numeric(0)
      eta1 <- if (q1 > 0) xb1 + as.vector(X1[, rand_idx1, drop = FALSE] %*% dev1) else xb1
      eta2 <- if (q2 > 0) xb2 + as.vector(X2[, rand_idx2, drop = FALSE] %*% dev2) else xb2
      mu1_r <- pmin(pmax(exp(eta1), 1e-300), 1e15)
      mu2_r <- pmin(pmax(exp(eta2), 1e-300), 1e15)
      b <- lambda_bounds_vec(c_val(mu1_r, m1_start), c_val(mu2_r, m2_start))
      lamLo <- max(lamLo, b[1]); lamHi <- min(lamHi, b[2])
    }
    if (!(lamLo < lamHi)) {
      stop("Invalid lambda bounds at starting values.")
    }
  }

  # Build TMB object
  tmb_data <- .build_tmb_data(Y1, Y2, X1, X2, rand_idx1, rand_idx2,
                              Z1, Z2, dist1, dist2, sign1, sign2,
                              family_code, poisson_1, poisson_2,
                              lamLo, lamHi)
  # Map fixed parameters (e.g., pinned log_m for Poisson margins)
  map <- list()
  if (length(fixed_names)) {
    for (fn in fixed_names) {
      map[[fn]] <- factor(NA)  # fix to starting value
    }
  }
  # Independence: z_dep is in the C++ template but not estimated
  if (family_code < 0L) {
    map[["z_dep"]] <- factor(NA)
  }

  i1 <- 1:k1; i2 <- (k1 + 1):(k1 + k2)
  idx_end <- k1 + k2 + q1 + q2

  configured <- .make_rpbnb_tmb_object(
    data = tmb_data,
    parameters = list(
      beta1 = start[i1],
      beta2 = start[i2],
      log_sd1 = if (q1 > 0) start[(k1 + k2 + 1):(k1 + k2 + q1)] else numeric(0),
      log_sd2 = if (q2 > 0) start[(k1 + k2 + q1 + 1):(k1 + k2 + q1 + q2)] else numeric(0),
      log_m1 = start[idx_end + 1],
      log_m2 = start[idx_end + 2],
      z_dep = if (family_code >= 0L) start[idx_end + 3] else 0
    ),
    map = if (length(map) > 0) map else NULL,
    silent = control$print_level == 0L,
    n_cores = configured_threads,
    max_threads = configured_threads,
    parallel_tape = control$parallel_tape
  )
  obj <- configured$obj

  # Optimize with nlminb
  opt <- stats::nlminb(
    start = obj$par,
    objective = obj$fn,
    gradient = obj$gr,
    control = list(
      iter.max = control$iterlim,
      eval.max = control$iterlim * 2,
      rel.tol = control$reltol,
      trace = max(0, control$print_level - 1)
    )
  )

  # Build named coefficient vector matching par_names
  coef_vec <- start
  coef_vec[free] <- opt$par
  names(coef_vec) <- par_names

  # Construct result
  value <- opt$objective
  ll_hat <- -value  # nll -> logLik

  inference_result <- .rpbnb_inference(
    obj = obj,
    par = opt$par,
    coef = coef_vec,
    par_names = par_names,
    free = free,
    mode = inference,
    family_code = family_code,
    lamLo = lamLo,
    lamHi = lamHi
  )
  m1_hat <- as.numeric(inference_result$report["m1"])
  m2_hat <- as.numeric(inference_result$report["m2"])
  .warn_boundary_report(
    inference_result$boundary_report, family_code,
    sides = inference_result$boundary_sides
  )

  rp_meta <- list(
    Z1 = Z1, Z2 = Z2,
    dist1 = spec1$dist, dist2 = spec2$dist,
    sign1 = sign1, sign2 = sign2
  )

  result <- list(
    coef = coef_vec,
    vcov = inference_result$vcov,
    vcov_diag = inference_result$vcov_diag,
    se = inference_result$se,
    logLik = ll_hat,
    nobs = n,
    npar = npar,
    m1 = m1_hat,
    m2 = m2_hat,
    dependence = dependence,
    X1 = if (keep == "compact") NULL else X1,
    X2 = if (keep == "compact") NULL else X2,
    Y1 = if (keep == "full") Y1 else NULL,
    Y2 = if (keep == "full") Y2 else NULL,
    rand_idx1 = rand_idx1,
    rand_idx2 = rand_idx2,
    rp_meta = if (keep == "compact") NULL else rp_meta,
    model_meta = prep$model_meta,
    optimizer = opt,
    sdreport = inference_result$sdreport,
    obj = if (keep == "full") obj else NULL,
    inference = inference,
    boundary_report = inference_result$boundary_report,
    # Kept on the fit, not only in the warning: the two sides of a clamp call
    # for opposite remedies, and a suppressed or long-since-scrolled-past
    # warning must not be the only place that distinction survives.
    boundary_sides = inference_result$boundary_sides,
    # Frozen at the starting values, not recomputed at the optimum, so a
    # Famoye lambda estimate can be pinned by a bound the data never implied.
    # Retained so that artefact is inspectable rather than invisible.
    lambda_bounds = if (family_code == 0L) c(lower = lamLo, upper = lamHi),
    keep = keep,
    parallel = list(requested = control$n_cores,
                    realized = configured$n_cores),
    call = match.call()
  )
  class(result) <- "rpbnb_tmb_fit"
  result
}
