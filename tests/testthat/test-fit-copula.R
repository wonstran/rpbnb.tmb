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

test_that("copula margins reproduce the exact NB2 mass function", {
  skip_on_cran()
  # The copula likelihood needs the NB2 CDF at y and y - 1, and their
  # difference has to be the NB2 pmf to full precision -- p_obs is a second
  # difference of copula CDFs, so any error in a margin is amplified, not
  # averaged away.
  #
  # Frank at z_dep = 0 is exactly the independence copula (the template's
  # near-independence branch returns u * v verbatim), so the joint density
  # collapses to the product of the two marginal masses. That makes this an
  # exact check of the template's CDF against R's dnbinom rather than a
  # tolerance-tuned approximation.
  # Both a small-count and a large-count regime. The second matters on its own:
  # the CDF is accumulated term by term, so rounding grows with the count, and
  # the truck workload this package is built for reaches counts of 266.
  scenarios <- list(
    small = list(b1 = c(0.3, 0.5), b2 = c(0.1, -0.2), m1 = 0.4, m2 = 0.5),
    large = list(b1 = c(3.0, 0.5), b2 = c(2.5, -0.2), m1 = 0.4, m2 = 0.5)
  )
  for (nm in names(scenarios)) {
    s <- scenarios[[nm]]
    set.seed(11)
    n <- 200
    x1 <- rnorm(n)
    mu1 <- exp(s$b1[1] + s$b1[2] * x1)
    mu2 <- exp(s$b2[1] + s$b2[2] * x1)
    dat <- data.frame(
      y1 = rnbinom(n, mu = mu1, size = 1 / s$m1),
      y2 = rnbinom(n, mu = mu2, size = 1 / s$m2),
      x1 = x1
    )
    fit <- fit_rpbnb_tmb(y1 ~ x1, y2 ~ x1, data = dat,
                         dependence = copula("frank"), draws = 2,
                         keep = "full",
                         control = rpbnb_tmb_control(iterlim = 1, n_cores = 1))

    # beta1 (2), beta2 (2), log_m1, log_m2, z_dep -- no random coefficients.
    par <- c(s$b1, s$b2, log(s$m1), log(s$m2), 0)
    expect_identical(length(fit$obj$par), length(par))

    expected <-
      -sum(stats::dnbinom(dat$y1, mu = mu1, size = 1 / s$m1, log = TRUE)) -
      sum(stats::dnbinom(dat$y2, mu = mu2, size = 1 / s$m2, log = TRUE))

    expect_equal(as.numeric(fit$obj$fn(par)), expected,
                 tolerance = 1e-10, label = nm)
    # Not vacuous only if the "large" scenario really does exercise long
    # summations; 266 is the truck data's maximum.
    if (nm == "large") expect_gt(max(dat$y1), 100)
  }
})

test_that("Frank's closed-form cell probability matches the second difference", {
  skip_on_cran()
  # frank_cell_prob() replaces C(a,b) - C(a',b) - C(a,b') + C(a',b') with a
  # single log1p so the tail stops cancelling. That rearrangement is exact, so
  # in a WELL-CONDITIONED regime -- small counts, where the naive second
  # difference is still accurate -- the two must agree to near machine
  # precision. This is the check on the algebra; the tail behaviour it was
  # written for is what test-laplace.R exercises.
  set.seed(7)
  n <- 150
  x1 <- rnorm(n)
  b1 <- c(0.2, 0.3); b2 <- c(0.1, -0.2); m1 <- 0.4; m2 <- 0.5
  mu1 <- exp(b1[1] + b1[2] * x1)
  mu2 <- exp(b2[1] + b2[2] * x1)
  dat <- data.frame(y1 = rnbinom(n, mu = mu1, size = 1 / m1),
                    y2 = rnbinom(n, mu = mu2, size = 1 / m2),
                    x1 = x1)
  fit <- fit_rpbnb_tmb(y1 ~ x1, y2 ~ x1, data = dat,
                       dependence = copula("frank"), draws = 2, keep = "full",
                       control = rpbnb_tmb_control(iterlim = 1, n_cores = 1))

  z_dep <- 2
  theta <- FRANK_THETA_MAX * tanh(z_dep / FRANK_THETA_MAX)
  frank_cdf <- function(u, v, th) {
    -log1p(expm1(-th * u) * expm1(-th * v) / expm1(-th)) / th
  }
  a <- stats::pnbinom(dat$y1, size = 1 / m1, mu = mu1)
  am <- stats::pnbinom(dat$y1 - 1, size = 1 / m1, mu = mu1)
  b <- stats::pnbinom(dat$y2, size = 1 / m2, mu = mu2)
  bm <- stats::pnbinom(dat$y2 - 1, size = 1 / m2, mu = mu2)
  p <- frank_cdf(a, b, theta) - frank_cdf(am, b, theta) -
    frank_cdf(a, bm, theta) + frank_cdf(am, bm, theta)
  # Guard the premise: if any cell probability were near the rounding floor the
  # reference itself would be noise and the comparison would prove nothing.
  expect_gt(min(p), 1e-8)

  par <- c(b1, b2, log(m1), log(m2), z_dep)
  expect_equal(as.numeric(fit$obj$fn(par)), -sum(log(p)), tolerance = 1e-10)
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
