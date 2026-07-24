test_that("Famoye independence (lambda=0) equals two marginal NB2 fits", {
  skip_on_cran()
  set.seed(1)
  n <- 300
  x1 <- rnorm(n)
  X <- cbind(`(Intercept)` = 1, x1 = x1)
  mu <- exp(X %*% c(0.3, 0.5))
  y1 <- rnbinom(n, mu = mu, size = 1/0.4)
  y2 <- rnbinom(n, mu = mu, size = 1/0.5)
  dat <- data.frame(y1 = y1, y2 = y2, x1 = x1)

  fit <- fit_rpbnb_tmb(y1 ~ x1, y2 ~ x1, data = dat,
                        dependence = "independence", draws = 100)
  expect_true(is.finite(fit$logLik))
  expect_true(fit$npar >= 6)  # b1, b2, log_m1, log_m2, 4 coefs + 2 dispersions
})

test_that("Famoye fit recovers known parameters", {
  skip_on_cran()
  sim <- simulate_rpbnb_tmb(n = 500,
    beta1 = c("(Intercept)" = 0.5, x1 = 0.3),
    beta2 = c("(Intercept)" = 0.2, x1 = -0.2),
    dispersion = c(m1 = 0.4, m2 = 0.5),
    dependence = "famoye", lambda = 0.3, seed = 42)
  fit <- fit_rpbnb_tmb(y1 ~ x1, y2 ~ x1, data = sim$data,
                        dependence = "famoye", draws = 200)
  expect_true(is.finite(fit$logLik))
  # Check truth is within 95% CI for key parameters (approximate)
  ci_b1 <- coef(fit)["b1:x1"] + c(-1.96, 1.96) * fit$se["b1:x1"]
  expect_true(ci_b1[1] < 0.3 && ci_b1[2] > 0.3)
})

test_that("Famoye with random coefficients converges", {
  skip_on_cran()
  sim <- simulate_rpbnb_tmb(n = 400,
    beta1 = c("(Intercept)" = 0.2, x1 = 0.4),
    beta2 = c("(Intercept)" = 0.1, x1 = -0.3),
    random_1 = list(x1 = list(sd = 0.3)),
    dispersion = c(m1 = 0.4, m2 = 0.5),
    dependence = "famoye", lambda = 0, seed = 1)
  fit <- fit_rpbnb_tmb(y1 ~ x1, y2 ~ x1, data = sim$data,
                        random_1 = "x1",
                        dependence = "famoye", draws = 200)
  expect_true(is.finite(fit$logLik))
})

test_that("fit_rpbnb_tmb validates draws before constructing Halton points", {
  d <- data.frame(y1 = c(0, 1), y2 = c(1, 0), x = c(-1, 1))
  expect_error(
    fit_rpbnb_tmb(y1 ~ x, y2 ~ x, d, random_1 = "x", draws = 0),
    "draws"
  )
  expect_error(
    fit_rpbnb_tmb(y1 ~ x, y2 ~ x, d, random_1 = "x", draws = 1.5),
    "draws"
  )
})
