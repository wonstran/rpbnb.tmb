# The SML path must be bit-identical before and after the Laplace feature.
# BASELINE was captured from the pre-change code; see the implementation plan.
BASELINE_SML_LOGLIK <- -949.6478422037

sml_baseline_fit <- function() {
  set.seed(1)
  n <- 300
  x1 <- rnorm(n)
  X <- cbind(`(Intercept)` = 1, x1 = x1)
  mu <- exp(X %*% c(0.3, 0.5))
  dat <- data.frame(
    y1 = rnbinom(n, mu = mu, size = 1 / 0.4),
    y2 = rnbinom(n, mu = mu, size = 1 / 0.5),
    x1 = x1
  )
  fit_rpbnb_tmb(y1 ~ x1, y2 ~ x1, data = dat, random_1 = "x1",
                dependence = "famoye", draws = 50, seed = 7)
}

test_that("adding latent parameters leaves the SML tape unchanged", {
  skip_on_cran()
  fit <- sml_baseline_fit()
  expect_equal(fit$logLik, BASELINE_SML_LOGLIK, tolerance = 1e-10)
})

test_that("laplace rejects unsupported random-coefficient distributions", {
  skip_on_cran()
  set.seed(3)
  n <- 120
  x1 <- rnorm(n)
  dat <- data.frame(y1 = rpois(n, 2), y2 = rpois(n, 2), x1 = x1)
  expect_error(
    fit_rpbnb_tmb(y1 ~ x1, y2 ~ x1, data = dat,
                  random_1 = list(x1 = "uniform"),
                  dependence = "independence", draws = 20,
                  method = "laplace"),
    "uniform"
  )
})

test_that("laplace requires at least one random coefficient", {
  skip_on_cran()
  set.seed(4)
  n <- 120
  x1 <- rnorm(n)
  dat <- data.frame(y1 = rpois(n, 2), y2 = rpois(n, 2), x1 = x1)
  expect_error(
    fit_rpbnb_tmb(y1 ~ x1, y2 ~ x1, data = dat,
                  dependence = "independence", draws = 20,
                  method = "laplace"),
    "random coefficient"
  )
})

test_that("laplace keeps draws meaningful for post-estimation averaging", {
  skip_on_cran()
  set.seed(5)
  n <- 200
  x1 <- rnorm(n)
  dat <- data.frame(y1 = rpois(n, 2), y2 = rpois(n, 2), x1 = x1)
  expect_no_warning(
    fit <- fit_rpbnb_tmb(y1 ~ x1, y2 ~ x1, data = dat, random_1 = "x1",
                         dependence = "independence", draws = 25,
                         method = "laplace")
  )
  expect_identical(nrow(fit$rp_meta$Z1), 25L)
  expect_identical(fit$method, "laplace")
})
