#' @keywords internal
#' @noRd
POISSON_M <- 1e-6

#' @keywords internal
#' @noRd
.chk_poisson_flag <- function(x, arg) {
  if (!(is.logical(x) && length(x) == 1L && !is.na(x))) {
    stop("`", arg, "` must be a single non-missing logical.")
  }
}

#' @keywords internal
#' @noRd
.check_counts <- function(y, label) {
  y <- as.numeric(y)
  if (any(!is.finite(y))) stop("Response ", label, " has non-finite values.")
  if (any(y < 0)) stop("Response ", label, " has negative values.")
  if (any(abs(y - round(y)) > 1e-8)) stop("Response ", label, " has non-integer values.")
  as.integer(round(y))
}

#' @keywords internal
#' @noRd
.prepare_bnb_data <- function(formula_1, formula_2, data) {
  vars <- unique(c(all.vars(stats::terms(formula_1, data = data)),
                   all.vars(stats::terms(formula_2, data = data))))
  missing_vars <- vars[!vars %in% names(data)]
  if (length(missing_vars)) stop("Variables not found: ", paste(missing_vars, collapse = ", "))
  ok <- stats::complete.cases(data[, vars, drop = FALSE])
  if (!any(ok)) stop("No complete cases.")
  data_cc <- data[ok, , drop = FALSE]
  mf1 <- stats::model.frame(formula_1, data = data_cc)
  mf2 <- stats::model.frame(formula_2, data = data_cc)
  Y1 <- .check_counts(stats::model.response(mf1), "1")
  Y2 <- .check_counts(stats::model.response(mf2), "2")
  X1 <- stats::model.matrix(formula_1, mf1)
  X2 <- stats::model.matrix(formula_2, mf2)
  list(Y1 = Y1, Y2 = Y2, X1 = X1, X2 = X2,
       cn1 = colnames(X1), cn2 = colnames(X2),
       n = length(Y1), data = data_cc)
}

#' @keywords internal
#' @noRd
.bound_mu <- function(X, beta) {
  pmin(pmax(as.vector(exp(X %*% beta)), 1e-300), 1e15)
}

#' @keywords internal
#' @noRd
rand_dist_registry <- list(
  normal = list(
    base = "normal",
    u_to_base = function(u) stats::qnorm(u),
    coef = function(b, s, base, sign) b + s * base,
    dev = function(b, s, base, sign) s * base,
    dloc_factor = function(b, s, base, coef) rep(1, length(base)),
    dscale = function(b, s, base, coef) s * base,
    scale_label = "log_sd"
  ),
  lognormal = list(
    base = "normal",
    u_to_base = function(u) stats::qnorm(u),
    coef = function(b, s, base, sign) sign * exp(b + s * base),
    dev = function(b, s, base, sign) sign * exp(b + s * base) - b,
    dloc_factor = function(b, s, base, coef) coef,
    dscale = function(b, s, base, coef) coef * base * s,
    scale_label = "log_s"
  ),
  uniform = list(
    base = "uniform",
    u_to_base = function(u) u,
    coef = function(b, s, base, sign) b + s * (2 * base - 1),
    dev = function(b, s, base, sign) s * (2 * base - 1),
    dloc_factor = function(b, s, base, coef) rep(1, length(base)),
    dscale = function(b, s, base, coef) s * (2 * base - 1),
    scale_label = "log_w"
  ),
  triangular = list(
    base = "uniform",
    u_to_base = function(u) tri_icdf(u),
    coef = function(b, s, base, sign) b + s * base,
    dev = function(b, s, base, sign) s * base,
    dloc_factor = function(b, s, base, coef) rep(1, length(base)),
    dscale = function(b, s, base, coef) s * base,
    scale_label = "log_w"
  )
)

tri_icdf <- function(u) {
  ifelse(u < 0.5, -1 + sqrt(2 * u), 1 - sqrt(2 * (1 - u)))
}

parse_rand_spec <- function(spec) {
  if (is.null(spec) || length(spec) == 0) {
    return(list(names = character(0), dist = character(0),
                sign = numeric(0), scale = numeric(0)))
  }
  valid <- names(rand_dist_registry)
  if (is.character(spec) && is.null(names(spec))) {
    return(list(names = spec, dist = rep("normal", length(spec)),
                sign = rep(1, length(spec)), scale = rep(NA_real_, length(spec))))
  }
  if (!is.list(spec) || is.null(names(spec)) || any(!nzchar(names(spec)))) {
    stop("random spec must be a character vector or named list.")
  }
  nm <- names(spec)
  dist <- character(length(nm))
  sgn <- numeric(length(nm))
  scale <- numeric(length(nm))
  for (i in seq_along(nm)) {
    v <- spec[[i]]
    if (is.character(v) && length(v) == 1L) {
      d <- v; this_sign <- 1; this_scale <- NA_real_
    } else if (is.list(v)) {
      d <- if (is.null(v$dist)) "normal" else v$dist
      this_sign <- if (is.null(v$sign)) 1 else v$sign
      this_scale <- if (!is.null(v$scale)) v$scale else if (!is.null(v$sd)) v$sd else NA_real_
      if (!is.null(v$sign) && d != "lognormal") {
        stop("`sign` is only for lognormal.")
      }
    } else {
      stop("Invalid random spec for '", nm[i], "'.")
    }
    if (!d %in% valid) stop("Unknown distribution '", d, "'.")
    if (!this_sign %in% c(-1, 1)) stop("`sign` must be -1 or 1.")
    dist[i] <- d; sgn[i] <- this_sign; scale[i] <- this_scale
  }
  list(names = nm, dist = dist, sign = sgn, scale = scale)
}

chk_rand_spec <- function(spec, bv, lbl) {
  miss <- spec$names[!spec$names %in% names(bv)]
  if (length(miss)) stop("random name(s) ", paste(miss, collapse = ", "), " not in ", lbl, ".")
  if (length(spec$names)) {
    if (any(is.na(spec$scale))) stop("Each random coefficient needs a `scale`/`sd` in ", lbl, ".")
    if (any(!is.finite(spec$scale)) || any(spec$scale <= 0)) {
      stop("Random scales must be finite and positive in ", lbl, ".")
    }
  }
}

chk_dispersion <- function(dispersion) {
  if (!all(c("m1", "m2") %in% names(dispersion))) {
    stop("`dispersion` must be c(m1 = ., m2 = .).")
  }
  m <- dispersion[c("m1", "m2")]
  if (any(!is.finite(m)) || any(m <= 0)) stop("Dispersions must be finite and positive.")
}

rand_realize <- function(U, dist, sign, b, s) {
  R <- nrow(U); q <- ncol(U)
  base <- coef <- dev <- dloc <- dscale <- matrix(0, nrow = R, ncol = q)
  for (j in seq_len(q)) {
    reg <- rand_dist_registry[[dist[j]]]
    bj <- reg$u_to_base(U[, j])
    cj <- reg$coef(b[j], s[j], bj, sign[j])
    base[, j] <- bj
    coef[, j] <- cj
    dev[, j] <- reg$dev(b[j], s[j], bj, sign[j])
    dloc[, j] <- reg$dloc_factor(b[j], s[j], bj, cj)
    dscale[, j] <- reg$dscale(b[j], s[j], bj, cj)
  }
  list(base = base, coef = coef, dev = dev, dloc = dloc, dscale = dscale)
}

d_const <- function() 1 - exp(-1)

c_val <- function(mu, m) {
  if (m == 0) return(exp(-d_const() * mu))
  (1 + d_const() * m * mu)^(-1 / m)
}

lambda_bounds_vec <- function(c1, c2) {
  lam_min <- -1 / pmax((1 - c1) * (1 - c2), c1 * c2)
  lam_max <-  1 / pmax(c1 * (1 - c2), c2 * (1 - c1))
  c(max(lam_min), min(lam_max))
}

.resolve_start <- function(start, default, par_names, label = "start") {
  if (is.null(start)) { names(default) <- par_names; return(default) }
  if (!is.numeric(start)) stop("`", label, "` must be numeric.")
  nm <- names(start)
  if (is.null(nm) || all(!nzchar(nm))) {
    if (length(start) != length(par_names)) {
      stop("`", label, "` must have length ", length(par_names), ".")
    }
    names(start) <- par_names
    return(start)
  }
  if (any(!nzchar(nm))) stop("Mix of named/unnamed in `", label, "`.")
  if (anyDuplicated(nm)) stop("Duplicate names in `", label, "`.")
  unknown <- setdiff(nm, par_names)
  if (length(unknown)) stop("Unknown names in `", label, "`: ", paste(unknown, collapse = ", "))
  out <- default; names(out) <- par_names
  out[nm] <- start
  out
}

.observed_info_vcov <- function(info, par_names, label = "model") {
  info <- (info + t(info)) / 2
  ev <- try(eigen(info, symmetric = TRUE, only.values = TRUE), silent = TRUE)
  bad_eig <- inherits(ev, "try-error") || any(!is.finite(ev$values))
  min_eig <- if (bad_eig) NA_real_ else min(ev$values)
  max_eig <- if (bad_eig) NA_real_ else max(ev$values)
  pd <- isTRUE(!bad_eig && min_eig > 0)
  ridge <- 0
  if (!pd) {
    ridge <- if (bad_eig) 1e-2 else 1e-8 - min_eig
    info <- info + diag(ridge, nrow(info))
  }
  vc <- try(solve(info), silent = TRUE)
  inv_method <- "solve"
  if (inherits(vc, "try-error")) { vc <- .pseudo_inv(info); inv_method <- "pseudo_inv" }
  dimnames(vc) <- list(par_names, par_names)
  se <- sqrt(pmax(diag(vc), 0)); names(se) <- par_names
  cond <- if (pd) max_eig / min_eig else NA_real_
  list(vcov = vc, se = se,
       diag = list(min_eigenvalue = min_eig, max_eigenvalue = max_eig,
                   condition = cond, ridge = ridge, inversion = inv_method,
                   positive_definite = pd, repaired = !pd, label = label))
}

# Moore-Penrose pseudoinverse via SVD (base R, no MASS dependency)
.pseudo_inv <- function(x, tol = sqrt(.Machine$double.eps)) {
  s <- svd(x)
  d <- s$d
  d[d > tol] <- 1 / d[d > tol]
  d[d <= tol] <- 0
  s$v %*% diag(d, nrow = length(d)) %*% t(s$u)
}

.free_index_vcov <- function(info, par_names, free, label = "model") {
  if (all(free)) return(.observed_info_vcov(info, par_names, label = label))
  sub <- .observed_info_vcov(info[free, free, drop = FALSE], par_names[free], label = label)
  p <- length(par_names)
  vc <- matrix(NA_real_, p, p, dimnames = list(par_names, par_names))
  vc[free, free] <- sub$vcov
  se <- rep(NA_real_, p); names(se) <- par_names
  se[free] <- sub$se
  list(vcov = vc, se = se, diag = sub$diag)
}

#' @export
copula <- function(family, par = NULL) {
  family <- match.arg(family, c("frank", "normal", "kimeldorf"))
  if (!is.null(par)) {
    if (family == "normal" && (abs(par) >= 1)) stop("|rho| must be < 1.")
    if (family == "kimeldorf" && par <= 0) stop("Clayton theta must be > 0.")
    if (family == "frank" && !is.finite(par)) stop("Frank theta must be finite.")
  }
  structure(list(family = family, par = par), class = "rpbnb_copula")
}

#' Control parameters for rpbnb_tmb estimators
#' @param iterlim Maximum optimizer iterations.
#' @param reltol Relative convergence tolerance.
#' @param print_level Verbosity for nlminb.
#' @param n_cores Number of OpenMP threads for TMB.
#' @param max_threads Maximum OpenMP threads permitted for one fit. Increase
#'   explicitly to opt into a larger memory footprint.
#' @param max_workload Maximum weighted observation-draw evaluations permitted
#'   before TMB tape construction. Use \code{Inf} to disable the guard.
#' @param halton_burn Number of leading Halton points to discard.
#' @export
rpbnb_tmb_control <- function(iterlim = 500L,
                              reltol = 1e-8,
                              print_level = 0L,
                              n_cores = 1L,
                              max_threads = 4L,
                              max_workload = 2e6,
                              halton_burn = 300L) {
  if (length(n_cores) != 1L || !is.numeric(n_cores) ||
      is.na(n_cores) || !is.finite(n_cores) ||
      n_cores < 1 || n_cores != floor(n_cores) ||
      n_cores > .Machine$integer.max) {
    stop("n_cores must be one whole number greater than or equal to 1.",
         call. = FALSE)
  }
  if (length(max_threads) != 1L || !is.numeric(max_threads) ||
      is.na(max_threads) || !is.finite(max_threads) ||
      max_threads < 1 || max_threads != floor(max_threads) ||
      max_threads > .Machine$integer.max) {
    stop("max_threads must be one whole number greater than or equal to 1.",
         call. = FALSE)
  }
  if (length(max_workload) != 1L || !is.numeric(max_workload) ||
      is.na(max_workload) || max_workload <= 0) {
    stop("max_workload must be one positive number or Inf.", call. = FALSE)
  }
  structure(list(iterlim = as.integer(iterlim), reltol = reltol,
                  print_level = as.integer(print_level),
                  n_cores = as.integer(n_cores),
                  max_threads = as.integer(max_threads),
                  max_workload = as.numeric(max_workload),
                  halton_burn = as.integer(halton_burn)),
            class = "rpbnb_tmb_control")
}
