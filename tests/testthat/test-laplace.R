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
