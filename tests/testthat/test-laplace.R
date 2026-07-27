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
  # The validation runs over the union of both equations' distributions, so a
  # bad distribution supplied via random_2 alone must be caught too -- and
  # triangular (untested above) must be rejected in either equation.
  expect_error(
    fit_rpbnb_tmb(y1 ~ x1, y2 ~ x1, data = dat,
                  random_2 = list(x1 = "uniform"),
                  dependence = "independence", draws = 20,
                  method = "laplace"),
    "uniform"
  )
  expect_error(
    fit_rpbnb_tmb(y1 ~ x1, y2 ~ x1, data = dat,
                  random_1 = list(x1 = "triangular"),
                  dependence = "independence", draws = 20,
                  method = "laplace"),
    "triangular"
  )
  expect_error(
    fit_rpbnb_tmb(y1 ~ x1, y2 ~ x1, data = dat,
                  random_2 = list(x1 = "triangular"),
                  dependence = "independence", draws = 20,
                  method = "laplace"),
    "triangular"
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
  # The defining behavior of method = "laplace" is that u1 (n x q1 = 200 x 1)
  # is integrated out as a TMB random effect rather than optimized directly.
  # If `random` were dropped, those 200 latents would join the optimized
  # ("outer") parameter vector, so length(fit$optimizer$par) would be
  # 7 + 200 = 207 instead of 7. The exact count below (beta1: 2, beta2: 2,
  # log_sd1: 1, log_m1: 1, log_m2: 1 -- z_dep is mapped off under
  # dependence = "independence", and log_sd2 doesn't exist since random_2
  # is unused) pins that no latent leaked into the optimizer and that no
  # fixed parameter was accidentally freed or mapped off.
  expect_identical(length(fit$optimizer$par), 7L)
})

test_that("laplace engages TMB's sparse random-effect machinery", {
  skip_on_cran()
  sim <- simulate_rpbnb_tmb(
    n = 200,
    beta1 = c("(Intercept)" = 0.2, x1 = 0.4),
    beta2 = c("(Intercept)" = 0.1, x1 = -0.3),
    random_1 = list(x1 = list(sd = 0.3)),
    dispersion = c(m1 = 0.4, m2 = 0.5),
    dependence = "famoye", lambda = 0, seed = 11
  )
  fit <- fit_rpbnb_tmb(y1 ~ x1, y2 ~ x1, data = sim$data,
                       random_1 = "x1", dependence = "famoye",
                       draws = 50, method = "laplace", keep = "full")
  # One latent per observation per random coefficient. If this is 0 the fit
  # silently ran a fixed-effect model at u = 0 instead of integrating.
  expect_identical(length(fit$obj$env$random), 200L)
  expect_true(is.finite(fit$logLik))
})

test_that("laplace and sml agree on data simulated from the true process", {
  skip_on_cran()
  sim <- simulate_rpbnb_tmb(
    n = 500,
    beta1 = c("(Intercept)" = 0.5, x1 = 0.3),
    beta2 = c("(Intercept)" = 0.2, x1 = -0.2),
    random_1 = list(x1 = list(sd = 0.4)),
    dispersion = c(m1 = 0.4, m2 = 0.5),
    dependence = "famoye", lambda = 0.2, seed = 21
  )
  common <- list(formula_1 = y1 ~ x1, formula_2 = y2 ~ x1,
                 data = sim$data, random_1 = "x1",
                 dependence = "famoye", draws = 200, seed = 99)

  fit_sml <- do.call(fit_rpbnb_tmb, c(common, list(method = "sml")))
  fit_lap <- do.call(fit_rpbnb_tmb, c(common, list(method = "laplace")))

  expect_true(is.finite(fit_lap$logLik))
  expect_true(isTRUE(fit_lap$sdreport$pdHess))

  # Same parameter vector, same names, so a direct comparison is meaningful.
  expect_identical(names(coef(fit_lap)), names(coef(fit_sml)))

  # Agreement judged against sampling noise rather than an absolute tolerance:
  # these are two approximations to the same integral, not two computations of
  # the same number.
  key <- c("b1:x1", "b2:x1")
  for (nm in key) {
    tol <- 3 * max(fit_sml$se[[nm]], fit_lap$se[[nm]])
    expect_lt(abs(coef(fit_sml)[[nm]] - coef(fit_lap)[[nm]]), tol)
  }
})

test_that("laplace integrates random effects from both equations simultaneously", {
  skip_on_cran()
  n <- 250
  sim <- simulate_rpbnb_tmb(
    n = n,
    beta1 = c("(Intercept)" = 0.3, x1 = 0.3),
    beta2 = c("(Intercept)" = 0.1, x2 = -0.2),
    random_1 = list(x1 = list(sd = 0.3)),
    random_2 = list(x2 = list(sd = 0.25)),
    dispersion = c(m1 = 0.4, m2 = 0.5),
    dependence = "famoye", lambda = 0, seed = 41
  )
  fit <- fit_rpbnb_tmb(y1 ~ x1, y2 ~ x2, data = sim$data,
                       random_1 = "x1", random_2 = "x2",
                       dependence = "famoye",
                       draws = 50, method = "laplace", keep = "full")
  # q1 = q2 = 1 (one random slope in each equation). Task 6's acceptance fit
  # uses four random slopes in each equation, and every test up to this one
  # exercises random_1 alone (q2 == 0). If only the first equation's latents
  # reached MakeADFun -- e.g. random_names dropped "u2", or q2 silently
  # collapsed to 0 -- this would report 250 (n * q1) instead of 500
  # (n * (q1 + q2)), leaving u2 fixed at zero rather than integrated.
  expect_identical(length(fit$obj$env$random), as.integer(n * (1L + 1L)))
  expect_true(is.finite(fit$logLik))
})

test_that("laplace gives the same answer single- and multi-threaded", {
  skip_on_cran()
  sim <- simulate_rpbnb_tmb(
    n = 300,
    beta1 = c("(Intercept)" = 0.3, x1 = 0.35),
    beta2 = c("(Intercept)" = 0.1, x1 = -0.25),
    random_1 = list(x1 = list(sd = 0.35)),
    dispersion = c(m1 = 0.4, m2 = 0.5),
    dependence = "famoye", lambda = 0.1, seed = 31
  )
  common <- list(formula_1 = y1 ~ x1, formula_2 = y2 ~ x1,
                 data = sim$data, random_1 = "x1",
                 dependence = "famoye", draws = 50, seed = 5,
                 method = "laplace")

  fit1 <- do.call(fit_rpbnb_tmb, c(common, list(
    control = rpbnb_tmb_control(n_cores = 1L))))
  fit4 <- do.call(fit_rpbnb_tmb, c(common, list(
    control = rpbnb_tmb_control(n_cores = 4L, max_threads = 4L))))

  # Threading partitions the accumulation; it must not change the objective.
  expect_equal(fit1$logLik, fit4$logLik, tolerance = 1e-6)
  expect_equal(unname(coef(fit1)), unname(coef(fit4)), tolerance = 1e-4)
})
