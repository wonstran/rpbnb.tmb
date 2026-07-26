test_that("Frank random-parameter fit has finite curvature from default start", {
  skip_on_cran()
  data_path <- system.file(
    "extdata", "rwm1984_bnb.csv",
    package = "rpbnb.tmb", mustWork = TRUE
  )
  d <- utils::read.csv(data_path)[1:120, ]

  ctrl <- rpbnb_tmb_control(iterlim = 100, n_cores = 1)
  fit <- expect_no_warning(fit_rpbnb_tmb(
    docvis ~ age + hhninc + educ + female + married + kids,
    hospvis ~ age + educ + outwork + female + self,
    data = d,
    random_1 = "hhninc",
    random_2 = "educ",
    dependence = copula("frank"),
    draws = 20,
    seed = 20240712,
    control = ctrl
  ))

  expect_equal(fit$optimizer$convergence, 0)
  expect_true(isTRUE(fit$sdreport$pdHess))
  expect_true(all(is.finite(fit$se)))
  expect_gt(abs(unname(coef(fit)["z_dep"]) - 0.1), 1e-6)
})

test_that("Frank copula fit converges", {
  skip_on_cran()
  set.seed(1)
  n <- 300
  x1 <- rnorm(n)
  dat <- data.frame(y1 = rnbinom(n, mu = exp(0.3 + 0.5 * x1), size = 1/0.4),
                    y2 = rnbinom(n, mu = exp(0.1 - 0.2 * x1), size = 1/0.5),
                    x1 = x1)
  fit <- fit_rpbnb_tmb(y1 ~ x1, y2 ~ x1, data = dat,
                        dependence = copula("frank"), draws = 150)
  expect_true(is.finite(fit$logLik))
})

test_that("Gaussian copula fit converges", {
  skip_on_cran()
  set.seed(1)
  n <- 300
  x1 <- rnorm(n)
  dat <- data.frame(y1 = rnbinom(n, mu = exp(0.3 + 0.5 * x1), size = 1/0.4),
                    y2 = rnbinom(n, mu = exp(0.1 - 0.2 * x1), size = 1/0.5),
                    x1 = x1)
  fit <- fit_rpbnb_tmb(y1 ~ x1, y2 ~ x1, data = dat,
                        dependence = copula("normal"), draws = 150)
  expect_true(is.finite(fit$logLik))
})

test_that("Clayton copula fit converges", {
  skip_on_cran()
  set.seed(1)
  n <- 300
  x1 <- rnorm(n)
  dat <- data.frame(y1 = rnbinom(n, mu = exp(0.3 + 0.5 * x1), size = 1/0.4),
                    y2 = rnbinom(n, mu = exp(0.1 - 0.2 * x1), size = 1/0.5),
                    x1 = x1)
  fit <- fit_rpbnb_tmb(y1 ~ x1, y2 ~ x1, data = dat,
                        dependence = copula("kimeldorf"), draws = 150)
  expect_true(is.finite(fit$logLik))
})
