#' Draw one count from the Famoye conditional pmf without truncating its tail
#'
#' Replaces "enumerate y2 to the 0.999 NB quantile, renormalize, sample".
#' That discarded only about 0.1% of the mass, but for an over-dispersed NB
#' that 0.1% sits at the largest counts and carries a real share of the mean:
#' at mu = 5 the sampled second margin came out 4.1% low at dispersion 10 and
#' 13.8% low at dispersion 50, so `simulate_rpbnb_tmb()` did not reproduce the
#' NB2 margin it was asked for.
#'
#' The Famoye conditional pmf sums to exactly 1 over the full support, because
#' `c2 = E[exp(-Y2)]` makes the tilt integrate to zero. Inverse-CDF sampling
#' against that known total is therefore exact: accumulate probabilities from
#' zero upward in geometrically growing blocks and return the first count whose
#' cumulative probability reaches `u * total_mass`. Nothing is discarded and
#' nothing is renormalized.
#'
#' Famoye only. The copula path used to run through here too, with
#' `total_mass = P(Y1 = y1)`, but that required evaluating the copula CDF at
#' two arguments one marginal step apart -- which loses all precision once the
#' CDF saturates. It now samples copula uniforms directly; see
#' [.copula_conditional_inverse()].
#'
#' @param u One uniform draw.
#' @param total_mass Analytic total of `block_pmf` over the full support.
#' @param block_pmf Vectorized conditional pmf, evaluated on a run of counts.
#' @param first_block Size of the first block; later blocks double. Sized from
#'   the NB quantile so the common case still resolves in one vectorized call.
#' @param fallback Used only when `total_mass` is not usably positive, which
#'   means the conditioning event itself underflowed; falls back to the
#'   marginal inverse CDF rather than silently returning zero.
#' @return One count, always a double -- `qnbinom()` returns double and counts
#'   here can exceed integer range at extreme dispersions, so every path
#'   agrees on double rather than varying by which branch was taken.
#' @keywords internal
#' @noRd
.sample_conditional_count <- function(u, total_mass, block_pmf, first_block,
                                      fallback) {
  if (!is.finite(total_mass) || total_mass <= .Machine$double.eps) {
    return(as.numeric(fallback(u)))
  }
  target <- u * total_mass
  cumulative <- 0
  lo <- 0L
  size <- max(as.integer(first_block), 1L)
  # Terminates because the accumulated mass converges to total_mass while the
  # block size doubles; the ceiling is a backstop against a malformed pmf.
  while (lo <= 1e7) {
    p <- block_pmf(lo:(lo + size - 1L))
    p[!is.finite(p) | p < 0] <- 0
    cs <- cumulative + cumsum(p)
    hit <- which(cs >= target)
    if (length(hit)) return(as.numeric(lo + hit[1L] - 1L))
    cumulative <- cs[length(cs)]
    # Floating-point shortfall: the mass is spent but rounding kept the
    # cumulative sum just under the target. The last count is the answer.
    if (cumulative >= total_mass * (1 - 1e-12)) {
      return(as.numeric(lo + size - 1L))
    }
    lo <- lo + size
    size <- size * 2L
  }
  stop("Conditional count sampling did not converge; the conditional pmf ",
       "does not integrate to its stated total.", call. = FALSE)
}

#' Invert a copula's conditional distribution: V given U and a uniform W
#'
#' Sampling copula uniforms by conditional inversion and then applying
#' `qnbinom()` reproduces the intended joint pmf exactly, without ever
#' evaluating the copula CDF at two nearly equal arguments. That differencing
#' is what made the previous enumeration approach fail: for Frank at theta
#' near its supported ceiling the CDF saturates, so a rectangle spanning one
#' marginal probability step is computed as the difference of numbers that
#' agree to full precision, and the results ranged from identically zero to
#' 70x the correct mass.
#'
#' @param family One of "frank", "normal", "kimeldorf".
#' @param theta Natural-scale dependence parameter; 0 means independence.
#' @param u,w Equal-length uniform vectors. `u` is the first margin's uniform,
#'   `w` the independent uniform driving the conditional draw.
#' @return A uniform vector the same length as `u`.
#' @keywords internal
#' @noRd
.copula_conditional_inverse <- function(family, theta, u, w) {
  # runif() excludes the endpoints, but qnorm() is infinite there and a clamp
  # is cheaper than reasoning about the generator. EPS also bounds the output:
  # v is mathematically open in (0,1) but attains its limits in floating point
  # at extreme theta, and qnbinom(1, ...) is Inf.
  eps <- 1e-15
  u <- pmin(pmax(u, eps), 1 - eps)
  w <- pmin(pmax(w, eps), 1 - eps)

  v <- switch(family,
    frank = {
      if (abs(theta) < 1e-10) return(w)
      # Everything in log space. The naive form,
      #   v = (log D - log N)/theta,
      #   D = w + (1-w) e^{-theta u},  N = (1-w) e^{-theta u} + w e^{-theta},
      # is algebraically right but computes e^{-theta u} directly, which
      # overflows past theta = -710 and underflows past theta = +745 --
      # inside the domain copula() accepts, where it returned v identically 1
      # and hence Inf counts. D and N are each a sum of two positive terms, so
      # log-sum-exp evaluates them at any finite theta.
      log_w <- log(w)
      log_1mw <- log1p(-w)
      shared <- log_1mw - theta * u          # log of the (1-w) e^{-theta u} term
      (.log_sum_exp2(log_w, shared) -
         .log_sum_exp2(shared, log_w - theta)) / theta
    },
    normal = {
      # No bivariate-normal CDF needed here -- conditional inversion of a
      # Gaussian copula is one univariate normal draw, which is why this path
      # no longer depends on the pbivnorm package. |rho| < 1 is enforced by
      # copula(), so the scale below is strictly positive.
      stats::pnorm(theta * stats::qnorm(u) +
                     sqrt(1 - theta^2) * stats::qnorm(w))
    },
    kimeldorf = {
      if (theta < 1e-10) return(w)
      # v = (A*B + 1)^(-1/theta) with A = w^(-theta/(1+theta)) - 1 and
      # B = u^(-theta). B alone overflows around theta = 700 for mid-range u,
      # and the fitting link admits theta up to exp(20) ~ 4.9e8, at which the
      # direct form returned a point mass at zero rather than a uniform. A is
      # bounded because -theta/(1+theta) -> -1, so only B needs log space.
      log_a <- log(expm1(-theta / (1 + theta) * log(w)))
      log_ab <- log_a - theta * log(u)
      # v = exp(-softplus(log_ab) / theta), but forming softplus(log_ab) first
      # and dividing afterwards throws the answer away when theta * log(u)
      # overflows -- which happens for every u < exp(-1) once theta approaches
      # double.xmax, sending v to 0 for 37% of draws. Divide term by term:
      # for the upper branch
      #   softplus(L)/theta = log_a/theta + (-log u) + log1p(exp(-L))/theta,
      # so v = u * exp(-(log_a + log1p(exp(-L))) / theta). L then survives only
      # inside exp(-L), where an overflow to Inf yields the correct limit of 0
      # instead of destroying the result.
      #
      # Kept in that multiplicative form rather than exp(log(u) - c/theta):
      # it avoids an exp(log(u)) round-trip costing about an ulp, and it makes
      # the comonotonic limit exact instead of merely close. At theta near
      # double.xmax, c/theta lands around 1e-308 -- NOT underflowed to zero,
      # which is sixteen orders further down -- but far below half an ulp at
      # 1 (1.1e-16), so exp(-c/theta) rounds to exactly 1 and u * 1 preserves
      # u bitwise.
      ifelse(
        log_ab > 0,
        u * exp(-(log_a + log1p(exp(-log_ab))) / theta),
        exp(-log1p(exp(log_ab)) / theta)
      )
    },
    stop("Unknown copula family: ", family, call. = FALSE)
  )

  pmin(pmax(v, eps), 1 - eps)
}

#' Log of exp(a) + exp(b), without forming either exponential
#' @keywords internal
#' @noRd
.log_sum_exp2 <- function(a, b) {
  m <- pmax(a, b)
  out <- m + log1p(exp(-abs(a - b)))
  # Both terms zero: m - m is NaN above, but the sum really is zero.
  out[is.infinite(m) & m < 0] <- -Inf
  out
}

#' Simulate data from a bivariate NB process
#'
#' Generates count data from a bivariate negative binomial model with
#' random coefficients and Famoye/Sarmanov or copula dependence.
#'
#' @param n Number of observations.
#' @param beta1,beta2 Named coefficient vectors (must include "(Intercept)").
#' @param random_1,random_2 Random coefficient specs (same format as fit_rpbnb_tmb).
#' @param dispersion Named vector c(m1 = , m2 = ) of NB2 dispersions.
#' @param dependence "famoye", "independence", or a copula() object.
#' @param lambda Famoye dependence parameter (0 = independence).
#' @param covariates Optional data frame. If NULL, standard-normal covariates generated.
#' @param seed Optional integer seed.
#' @return List with $data (data.frame), $truth (parameters), $settings.
#' @export
#' @examples
#' sim <- simulate_rpbnb_tmb(n = 200,
#'   beta1 = c("(Intercept)" = 0.2, x1 = 0.4),
#'   beta2 = c("(Intercept)" = 0.1, x1 = -0.3),
#'   dispersion = c(m1 = 0.4, m2 = 0.5), seed = 1)
#' head(sim$data)
simulate_rpbnb_tmb <- function(n, beta1, beta2,
                               random_1 = NULL, random_2 = NULL,
                               dispersion = c(m1 = 0.5, m2 = 0.5),
                               dependence = "famoye",
                               lambda = 0, covariates = NULL, seed = NULL) {
  # Validity
  stopifnot("(Intercept)" %in% names(beta1), "(Intercept)" %in% names(beta2))
  chk_dispersion(dispersion)
  spec1 <- parse_rand_spec(random_1)
  spec2 <- parse_rand_spec(random_2)
  chk_rand_spec(spec1, beta1, "beta1")
  chk_rand_spec(spec2, beta2, "beta2")
  m1 <- dispersion["m1"]; m2 <- dispersion["m2"]
  r1 <- 1/m1; r2 <- 1/m2
  valid_dependence <- identical(dependence, "famoye") ||
    identical(dependence, "independence") ||
    inherits(dependence, "rpbnb_copula")
  if (!valid_dependence) {
    stop('`dependence` must be "famoye", "independence", or copula().',
         call. = FALSE)
  }

  if (!is.null(seed)) set.seed(seed)
  if (is.null(covariates)) {
    vars <- setdiff(union(names(beta1), names(beta2)), "(Intercept)")
    # An intercept-only model has no covariates, and `lapply()` over an empty
    # character vector yields a 0x0 frame -- which then collides with the
    # n-row outcomes in the final data.frame(). Carry n rows and no columns.
    covariates <- if (length(vars)) {
      stats::setNames(as.data.frame(lapply(vars, function(v) rnorm(n))), vars)
    } else {
      data.frame(row.names = seq_len(n))
    }
  }
  # Build design matrices
  X1 <- cbind(`(Intercept)` = rep(1, n))
  X2 <- cbind(`(Intercept)` = rep(1, n))
  for (nm in setdiff(names(beta1), "(Intercept)"))
    X1 <- cbind(X1, covariates[[nm]])
  for (nm in setdiff(names(beta2), "(Intercept)"))
    X2 <- cbind(X2, covariates[[nm]])
  colnames(X1) <- names(beta1); colnames(X2) <- names(beta2)
  beta1 <- beta1[colnames(X1)]; beta2 <- beta2[colnames(X2)]

  # Realized coefficients
  coef1 <- matrix(rep(beta1, each = n), n, length(beta1))
  colnames(coef1) <- names(beta1)
  coef2 <- matrix(rep(beta2, each = n), n, length(beta2))
  colnames(coef2) <- names(beta2)

  if (length(spec1$names)) {
    U1 <- matrix(runif(n * length(spec1$names)), n, length(spec1$names))
    r1z <- rand_realize(U1, spec1$dist, spec1$sign,
                        beta1[spec1$names], spec1$scale)
    for (j in seq_along(spec1$names)) {
      nm <- spec1$names[j]
      coef1[, nm] <- r1z$coef[, j]
    }
  }
  if (length(spec2$names)) {
    U2 <- matrix(runif(n * length(spec2$names)), n, length(spec2$names))
    r2z <- rand_realize(U2, spec2$dist, spec2$sign,
                        beta2[spec2$names], spec2$scale)
    for (j in seq_along(spec2$names)) {
      nm <- spec2$names[j]
      coef2[, nm] <- r2z$coef[, j]
    }
  }

  mu1 <- exp(rowSums(X1 * coef1))
  mu2 <- exp(rowSums(X2 * coef2))

  # Generate counts from joint distribution
  if (identical(dependence, "independence")) {
    y1 <- rnbinom(n, mu = mu1, size = r1)
    y2 <- rnbinom(n, mu = mu2, size = r2)
  } else if (identical(dependence, "famoye")) {
    if (identical(lambda, 0) || isTRUE(lambda == 0)) {
      y1 <- rnbinom(n, mu = mu1, size = r1)
      y2 <- rnbinom(n, mu = mu2, size = r2)
    } else {
    # Use conditional enumeration for Famoye
    c1 <- c_val(mu1, m1); c2 <- c_val(mu2, m2)
    bnds <- lambda_bounds_vec(c1, c2)
    if (lambda < bnds[1] || lambda > bnds[2])
      stop("lambda outside admissible bounds.")
    # Simulate via the exact conditional P(y2 | y1). Because c2 = E[exp(-Y2)]
    # under NB2, the Famoye tilt integrates to zero and the conditional pmf
    # sums to exactly 1 over the FULL support -- so inverse-CDF sampling
    # against a known total of 1 is exact and needs no renormalization.
    y1 <- rnbinom(n, mu = mu1, size = r1)
    y2 <- numeric(n)
    for (i in seq_len(n)) {
      tilt1 <- exp(-y1[i]) - c1[i]
      y2[i] <- .sample_conditional_count(
        u = runif(1),
        total_mass = 1,
        block_pmf = local({
          mu_i <- mu2[i]; c2_i <- c2[i]
          function(yv) {
            dnbinom(yv, mu = mu_i, size = r2) *
              (1 + lambda * tilt1 * (exp(-yv) - c2_i))
          }
        }),
        first_block = qnbinom(0.999, mu = mu2[i], size = r2) + 1L,
        fallback = function(u) qnbinom(u, mu = mu2[i], size = r2)
      )
    }
    }
  } else if (inherits(dependence, "rpbnb_copula")) {
    # Copula uniforms first, NB quantiles second -- not conditional-pmf
    # enumeration. Inverse-transforming copula uniforms gives exactly the
    # intended joint pmf, P(Y1=y1, Y2=y2) = the C-rectangle over the two
    # marginal CDF steps, but never evaluates that rectangle. Enumerating it
    # instead meant differencing a saturated copula CDF at two arguments a
    # marginal probability apart: at theta = 30 the evaluated rectangles
    # summed to anywhere from 0 to 70x the analytic marginal, which is a
    # cancellation artefact rather than a distributional statement. It is also
    # fully vectorized, where enumeration looped over observations.
    theta <- if (!is.null(dependence$par)) dependence$par else 0
    u1 <- runif(n)
    u2 <- .copula_conditional_inverse(dependence$family, theta, u1, runif(n))
    y1 <- qnbinom(u1, mu = mu1, size = r1)
    y2 <- qnbinom(u2, mu = mu2, size = r2)
  }

  data <- data.frame(y1 = y1, y2 = y2, covariates)
  truth <- list(beta1 = beta1, beta2 = beta2, random_1 = spec1,
                random_2 = spec2, dispersion = dispersion,
                lambda = lambda)
  list(data = data, truth = truth, coef1 = coef1, coef2 = coef2,
       mu1 = mu1, mu2 = mu2, settings = list(n = n, seed = seed))
}
