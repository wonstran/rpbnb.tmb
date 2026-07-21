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

  if (!is.null(seed)) set.seed(seed)
  if (is.null(covariates)) {
    vars <- setdiff(union(names(beta1), names(beta2)), "(Intercept)")
    covariates <- as.data.frame(lapply(vars, function(v) rnorm(n)))
    names(covariates) <- vars
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
  cop_family <- if (inherits(dependence, "rpbnb_copula")) dependence$family else NULL
  if (identical(dependence, "independence") || lambda == 0) {
    y1 <- rnbinom(n, mu = mu1, size = r1)
    y2 <- rnbinom(n, mu = mu2, size = r2)
  } else if (identical(dependence, "famoye")) {
    # Use accept-reject or conditional method for Famoye
    d <- d_const()
    c1 <- c_val(mu1, m1); c2 <- c_val(mu2, m2)
    bnds <- lambda_bounds_vec(c1, c2)
    if (lambda < bnds[1] || lambda > bnds[2])
      stop("lambda outside admissible bounds.")
    # Simulate via joint pmf enumeration for small counts + conditional
    y1 <- rnbinom(n, mu = mu1, size = r1)
    y2 <- numeric(n)
    for (i in seq_len(n)) {
      # Conditional pmf P(y2 | y1) proportional to P(y1,y2)
      # Enumerate a reasonable range
      max_y2 <- qnbinom(0.999, mu = mu2[i], size = r2)
      y2_vals <- 0:max_y2
      nb2 <- dnbinom(y2_vals, mu = mu2[i], size = r2)
      dep <- 1 + lambda * (exp(-y1[i]) - c1[i]) * (exp(-y2_vals) - c2[i])
      dep <- pmax(dep, 0)
      probs <- nb2 * dep
      probs <- probs / sum(probs)
      y2[i] <- sample(y2_vals, 1, prob = probs)
    }
  } else if (inherits(dependence, "rpbnb_copula")) {
    y1 <- rnbinom(n, mu = mu1, size = r1)
    y2 <- numeric(n)
    for (i in seq_len(n)) {
      max_y2 <- qnbinom(0.999, mu = mu2[i], size = r2)
      y2_vals <- 0:max_y2
      # Joint pmf via copula CDF rectangle
      theta <- if (!is.null(dependence$par)) dependence$par else 0
      pmf <- function(yv) {
        a <- pnbinom(y1[i], mu = mu1[i], size = r1)
        am <- if (y1[i] > 0) pnbinom(y1[i] - 1, mu = mu1[i], size = r1) else 0
        b <- pnbinom(yv, mu = mu2[i], size = r2)
        bm <- if (yv > 0) pnbinom(yv - 1, mu = mu2[i], size = r2) else 0
        switch(cop_family,
          frank = { C <- function(u,v) {
            if (abs(theta) < 1e-10) return(u*v)
            et <- exp(-theta)
            -log(1 + (exp(-theta*u)-1)*(exp(-theta*v)-1)/(et-1))/theta
          }},
          normal = { C <- function(u,v) {
            pbivnorm::pbivnorm(qnorm(pmax(pmin(u,1-1e-15),1e-15)),
                               qnorm(pmax(pmin(v,1-1e-15),1e-15)), theta)
          }},
          kimeldorf = { C <- function(u,v) {
            if (theta < 1e-10) return(u*v)
            pmax(u^(-theta) + v^(-theta) - 1, 0)^(-1/theta)
          }}
        )
        C(a,b) - C(am,b) - C(a,bm) + C(am,bm)
      }
      probs <- vapply(y2_vals, pmf, numeric(1))
      probs <- pmax(probs, 0)
      probs <- probs / sum(probs)
      y2[i] <- sample(y2_vals, 1, prob = probs)
    }
  }

  data <- data.frame(y1 = y1, y2 = y2, covariates)
  truth <- list(beta1 = beta1, beta2 = beta2, random_1 = spec1,
                random_2 = spec2, dispersion = dispersion,
                lambda = lambda)
  list(data = data, truth = truth, coef1 = coef1, coef2 = coef2,
       mu1 = mu1, mu2 = mu2, settings = list(n = n, seed = seed))
}
