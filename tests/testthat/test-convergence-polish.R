# Covers the two changes made after the truck Gaussian/Laplace fit returned
# nlminb code 0, a reassuring "relative convergence (4)", and NA standard
# errors on five equation-2 coefficients: the convergence-polish loop, and
# Poisson-limit margins under copula dependence.

polish_data <- function(n = 250, seed = 21) {
  set.seed(seed)
  x1 <- rnorm(n)
  mu1 <- exp(0.4 + 0.5 * x1)
  mu2 <- exp(0.2 - 0.3 * x1)
  data.frame(
    y1 = rpois(n, mu1),
    y2 = rpois(n, mu2),
    x1 = x1
  )
}

test_that("rpbnb_tmb_control validates gradtol and restarts", {
  expect_error(rpbnb_tmb_control(gradtol = 0), "gradtol")
  expect_error(rpbnb_tmb_control(gradtol = -1), "gradtol")
  expect_error(rpbnb_tmb_control(gradtol = c(1e-5, 1e-5)), "gradtol")
  expect_error(rpbnb_tmb_control(gradtol = NA_real_), "gradtol")
  expect_error(rpbnb_tmb_control(restarts = -1L), "restarts")
  expect_error(rpbnb_tmb_control(restarts = 1.5), "restarts")
  expect_error(rpbnb_tmb_control(restarts = NA_integer_), "restarts")
  # Zero is the documented way back to the single-nlminb() path and must not
  # be rejected alongside the negatives.
  expect_identical(rpbnb_tmb_control(restarts = 0L)$restarts, 0L)
  expect_identical(rpbnb_tmb_control()$restarts, 10L)
  expect_identical(rpbnb_tmb_control()$gradtol, 1e-5)
})

test_that("the polish loop never returns a worse objective than one nlminb call", {
  skip_on_cran()
  d <- polish_data()
  ctrl0 <- rpbnb_tmb_control(n_cores = 1, restarts = 0L)
  ctrl <- rpbnb_tmb_control(n_cores = 1)

  # The data is Poisson, so a free dispersion collapses and the Hessian does
  # not factor -- that is the situation under test, not an accident, and both
  # fits warn about it.
  plain <- suppressWarnings(
    fit_rpbnb_tmb(y1 ~ x1, y2 ~ x1, data = d,
                  dependence = copula("normal"), draws = 20,
                  seed = 5, control = ctrl0)
  )
  polished <- suppressWarnings(
    fit_rpbnb_tmb(y1 ~ x1, y2 ~ x1, data = d,
                  dependence = copula("normal"), draws = 20,
                  seed = 5, control = ctrl)
  )

  expect_identical(plain$optimizer$restarts, 0L)
  # Monotonicity is the whole contract: a restart that lands above where it
  # started is discarded, so polishing can only ever move the objective down.
  expect_lte(polished$optimizer$objective, plain$optimizer$objective)
  expect_gte(logLik(polished), logLik(plain))
  expect_true(is.finite(polished$optimizer$max_abs_gradient))
  expect_true(polished$optimizer$gradient_tolerance > 0)
})

test_that("Poisson-limit margins fit under every copula family", {
  skip_on_cran()
  d <- polish_data()
  # "kimeldorf" is Clayton's name in copula().
  for (family in c("normal", "frank", "kimeldorf")) {
    fit <- fit_rpbnb_tmb(
      y1 ~ x1, y2 ~ x1, data = d,
      dependence = copula(family), draws = 20, seed = 5,
      poisson_1 = TRUE, poisson_2 = TRUE,
      control = rpbnb_tmb_control(n_cores = 1)
    )
    # log_m1/log_m2 are mapped out, so their standard errors are NA by
    # construction; every parameter that is actually estimated must have one.
    estimated <- setdiff(names(fit$se), c("log_m1", "log_m2"))
    expect_true(all(is.finite(fit$se[estimated])),
                info = paste("family:", family))
    expect_true(isTRUE(fit$sdreport$pdHess), info = paste("family:", family))
    expect_true(all(is.na(fit$se[c("log_m1", "log_m2")])),
                info = paste("family:", family))
  }
})

test_that("Poisson-limit margins recover a Poisson DGP better than a free dispersion", {
  skip_on_cran()
  # The point of allowing the limit under a copula is not tidiness. Data with
  # no residual over-dispersion drives m to the floor, and at m ~ 1e-6 the NB2
  # size 1/m is ~1e6, so the marginal CDF is a difference of logarithms that
  # agree to eleven digits: the value survives and the curvature does not.
  # Pinning the limit removes the parameter rather than estimating it into
  # that regime.
  d <- polish_data(n = 400, seed = 3)
  ctrl <- rpbnb_tmb_control(n_cores = 1)

  free <- suppressWarnings(
    fit_rpbnb_tmb(y1 ~ x1, y2 ~ x1, data = d,
                  dependence = copula("normal"), draws = 20,
                  seed = 5, control = ctrl)
  )
  pinned <- fit_rpbnb_tmb(y1 ~ x1, y2 ~ x1, data = d,
                          dependence = copula("normal"), draws = 20,
                          seed = 5, poisson_1 = TRUE, poisson_2 = TRUE,
                          control = ctrl)

  # The free fit should agree that there is nothing to estimate -- and it pays
  # for that agreement with a Hessian that will not factor, which is the whole
  # reason the pinned form had to become reachable under a copula.
  expect_lt(free$m1, 1e-3)
  expect_lt(free$m2, 1e-3)
  expect_false(isTRUE(free$sdreport$pdHess))
  expect_true(isTRUE(pinned$sdreport$pdHess))
  # Same likelihood to within the two collapsed parameters, but the pinned fit
  # pays for two fewer of them.
  expect_equal(as.numeric(logLik(pinned)), as.numeric(logLik(free)),
               tolerance = 1e-3)
  expect_identical(pinned$npar, free$npar - 2L)
  expect_equal(unname(coef(pinned)[1:4]), unname(coef(free)[1:4]),
               tolerance = 1e-2)
})

test_that("a non-positive-definite Hessian is warned about, not just left as NA", {
  skip_on_cran()
  # Previously the only signal that inference had failed was NA appearing in a
  # printed table; nlminb reported code 0 and a success message either way.
  set.seed(9)
  n <- 60
  x1 <- rnorm(n)
  # Two counts that are almost entirely zero leave the dependence and the
  # dispersions with nothing to identify them.
  d <- data.frame(y1 = rbinom(n, 1, 0.03), y2 = rbinom(n, 1, 0.03), x1 = x1)
  fit <- suppressWarnings(
    fit_rpbnb_tmb(y1 ~ x1, y2 ~ x1, data = d, dependence = copula("normal"),
                  draws = 10, seed = 2,
                  control = rpbnb_tmb_control(n_cores = 1))
  )
  if (isFALSE(fit$sdreport$pdHess)) {
    expect_warning(
      fit_rpbnb_tmb(y1 ~ x1, y2 ~ x1, data = d, dependence = copula("normal"),
                    draws = 10, seed = 2,
                    control = rpbnb_tmb_control(n_cores = 1)),
      "not positive definite"
    )
  } else {
    succeed("Hessian factored on this platform; nothing to warn about.")
  }
})
