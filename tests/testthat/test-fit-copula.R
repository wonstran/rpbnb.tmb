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

test_that("NB2 margin recursion survives P(Y = 0) underflow near the mode", {
  skip_on_cran()
  # nb2_cdf_pair() used to seed its recursion with the linear-space P(Y = 0)
  # and multiply forward from there. At mu = size = 2000, P(Y = 0) = 0.5^2000
  # underflows to exact 0 in double precision, and every later term stayed
  # zero because each step only ever multiplies the previous one -- even
  # though P(Y = 2000), right at the mode, is an ordinary representable
  # probability (0.0063). That handed the optimizer a flat objective pinned
  # at the 1e-300 cell-probability floor instead of the real likelihood.
  # Frank at z_dep = 0 is exactly independence, so the joint density is the
  # product of the two marginal masses -- an exact check against dnbinom().
  n <- 5
  dat <- data.frame(y1 = rep(2000, n), y2 = rep(2000, n))
  fit <- fit_rpbnb_tmb(y1 ~ 1, y2 ~ 1, data = dat,
                       dependence = copula("frank"), draws = 2, keep = "full",
                       control = rpbnb_tmb_control(iterlim = 1, n_cores = 1))

  par <- c(log(2000), log(2000), log(1 / 2000), log(1 / 2000), 0)
  expected <- -n * 2 * stats::dnbinom(2000, mu = 2000, size = 2000, log = TRUE)
  expect_equal(as.numeric(fit$obj$fn(par)), expected, tolerance = 1e-8)
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

test_that("Clayton's closed-form cell probability matches the second difference", {
  skip_on_cran()
  # Companion to the Frank check above. clayton_cell_prob() replaces
  # C(a,b) - C(a',b) - C(a,b') + C(a',b') with a factored form whose O(xy)
  # second-order term is written in closed form, so the tail stops cancelling.
  # The rearrangement is exact, so in a WELL-CONDITIONED regime the two must
  # agree to near machine precision.
  set.seed(5)
  n <- 200
  x1 <- rnorm(n)
  b1 <- c(0.9, 0.4); b2 <- c(1.1, -0.3); m1 <- 0.4; m2 <- 0.5
  mu1 <- exp(b1[1] + b1[2] * x1)
  mu2 <- exp(b2[1] + b2[2] * x1)
  dat <- data.frame(y1 = rnbinom(n, mu = mu1, size = 1 / m1),
                    y2 = rnbinom(n, mu = mu2, size = 1 / m2),
                    x1 = x1)
  fit <- fit_rpbnb_tmb(y1 ~ x1, y2 ~ x1, data = dat,
                       dependence = copula("kimeldorf"), draws = 2,
                       keep = "full",
                       control = rpbnb_tmb_control(iterlim = 1, n_cores = 1))

  clayton_cdf <- function(u, v, th) {
    ifelse(u == 0 | v == 0, 0,
           pmax(u^(-th) + v^(-th) - 1, 1e-300)^(-1 / th))
  }
  a  <- stats::pnbinom(dat$y1,     size = 1 / m1, mu = mu1)
  am <- stats::pnbinom(dat$y1 - 1, size = 1 / m1, mu = mu1)
  b  <- stats::pnbinom(dat$y2,     size = 1 / m2, mu = mu2)
  bm <- stats::pnbinom(dat$y2 - 1, size = 1 / m2, mu = mu2)

  # Both signs of dependence, since theta = exp(z_dep) is one-sided.
  for (z_dep in c(-1.5, 0, 1.2)) {
    th <- exp(z_dep)
    p <- clayton_cdf(a, b, th) - clayton_cdf(am, b, th) -
      clayton_cdf(a, bm, th) + clayton_cdf(am, bm, th)
    # Guard the premise: near the rounding floor the reference is itself noise
    # and the comparison would prove nothing.
    expect_gt(min(p), 1e-8)
    par <- c(b1, b2, log(m1), log(m2), z_dep)
    expect_equal(as.numeric(fit$obj$fn(par)), -sum(log(p)),
                 tolerance = 1e-10, label = paste("z_dep =", z_dep))
  }
})

test_that("Clayton cell probabilities stay positive in the far tail", {
  skip_on_cran()
  # The regime the closed form exists for. At the truck fit's own starting
  # values -- all slopes zero, so mu = 1 against counts running to 266 -- the
  # naive second difference returns a NON-POSITIVE probability for 243 of the
  # 3,487 observations and a strictly NEGATIVE one for 11, where the true cell
  # probabilities run down to 1e-136. Under SML that only corrupts the
  # objective; under Laplace the negative ones put negative curvature into the
  # inner Hessian and TMB's inner Newton fails outright. The closed form is
  # positive by construction, so the objective must be finite and far below
  # the value the 1e-300 floor would produce.
  path <- system.file("extdata", "export_dense_all.csv",
                      package = "rpbnb.tmb", mustWork = TRUE)
  d <- utils::read.csv(path)
  dat <- data.frame(y1 = d$ALL_3, y2 = d$C_DISTR)
  expect_gt(max(dat$y1), 200)  # not vacuous: the tail must actually be long

  fit <- fit_rpbnb_tmb(y1 ~ 1, y2 ~ 1, data = dat,
                       dependence = copula("kimeldorf"), draws = 2,
                       keep = "full",
                       control = rpbnb_tmb_control(iterlim = 1, n_cores = 1))

  par <- c(0, 0, log(0.5), log(0.5), 0)  # mu = 1, m = 0.5, theta = 1
  nll <- as.numeric(fit$obj$fn(par))
  expect_true(is.finite(nll))
  expect_true(all(is.finite(fit$obj$gr(par))))
  # 243 floored cells would add 243 * -log(1e-300) = about 167,860 nats. The
  # true value is roughly 42,900, so this separates the two unambiguously.
  expect_lt(nll, 1e5)
})

test_that("Gaussian's strip integral matches adaptive quadrature", {
  skip_on_cran()
  # gaussian_cell_prob() replaces the second difference of four bivariate
  # normal CDFs with the single conditional strip integral
  #   int_{q(a')}^{q(a)} phi(z) [Phi((q(b) - rho z)/s) - Phi((q(b') - rho z)/s)] dz
  # evaluated on a fixed 20-point Gauss-Legendre rule. The reference here is
  # stats::integrate() on the same integral -- adaptive Gauss-Kronrod, a
  # different algorithm entirely -- so this checks the rule is fine enough,
  # not merely that two copies of the same code agree.
  set.seed(9)
  n <- 150
  x1 <- rnorm(n)
  b1 <- c(0.8, 0.4); b2 <- c(1.0, -0.3); m1 <- 0.4; m2 <- 0.5
  mu1 <- exp(b1[1] + b1[2] * x1)
  mu2 <- exp(b2[1] + b2[2] * x1)
  dat <- data.frame(y1 = rnbinom(n, mu = mu1, size = 1 / m1),
                    y2 = rnbinom(n, mu = mu2, size = 1 / m2),
                    x1 = x1)
  fit <- fit_rpbnb_tmb(y1 ~ x1, y2 ~ x1, data = dat,
                       dependence = copula("normal"), draws = 2, keep = "full",
                       control = rpbnb_tmb_control(iterlim = 1, n_cores = 1))

  safe_qnorm <- function(p) {
    stats::qnorm(pmin(pmax(p, 1e-15), 1 - 1e-15))
  }
  cell <- function(qa, qam, qb, qbm, rho) {
    s <- sqrt(1 - rho^2)
    if (qa <= qam) return(0)
    stats::integrate(
      function(z) stats::dnorm(z) *
        (stats::pnorm((qb - rho * z) / s) - stats::pnorm((qbm - rho * z) / s)),
      lower = qam, upper = qa, rel.tol = 1e-12
    )$value
  }

  qa  <- safe_qnorm(stats::pnbinom(dat$y1,     size = 1 / m1, mu = mu1))
  qam <- safe_qnorm(stats::pnbinom(dat$y1 - 1, size = 1 / m1, mu = mu1))
  qb  <- safe_qnorm(stats::pnbinom(dat$y2,     size = 1 / m2, mu = mu2))
  qbm <- safe_qnorm(stats::pnbinom(dat$y2 - 1, size = 1 / m2, mu = mu2))

  # z_dep is capped at 1.0 on purpose. By z_dep = 1.4 (rho = 0.885) the
  # smallest cell here is 7e-17, and at that size stats::integrate() and a
  # 2000-point Gauss-Legendre rule disagree with each other by 6e-3 -- the
  # reference stops being a reference, so a comparison there would measure
  # nothing.
  for (z_dep in c(-0.8, 0, 0.5, 1.0)) {
    rho <- tanh(z_dep)
    p <- mapply(cell, qa, qam, qb, qbm, MoreArgs = list(rho = rho))
    # Guard the premise: a reference near the rounding floor proves nothing.
    expect_gt(min(p), 1e-10)
    par <- c(b1, b2, log(m1), log(m2), z_dep)
    expect_equal(as.numeric(fit$obj$fn(par)), -sum(log(p)),
                 tolerance = 1e-8, label = paste("z_dep =", z_dep))
  }
})

test_that("Gaussian cell probabilities stay non-negative in the far tail", {
  skip_on_cran()
  # Companion to the Clayton tail test. At the truck starting values (mu = 1
  # against counts running to 266) the four-corner second difference returns a
  # non-positive cell probability for 457 to 600 of the 3,487 observations
  # depending on rho, of which 212 to 355 are strictly NEGATIVE -- the negative
  # curvature that stops the inner Newton. The strip integral is an integral of
  # a non-negative integrand over an ordered interval, so it cannot go negative.
  #
  # Note this asserts LESS than the Clayton counterpart: 245 of these cells hit
  # the qnorm clamp, because Gaussian must pass through qnorm(), which is
  # singular where the NB2 CDF saturates. Those still floor. The claim is that
  # nothing goes negative and the objective stays differentiable, not that
  # every cell is recovered.
  #
  # The threshold below is not arbitrary and not vacuous: 245 floored cells
  # contribute 245 * -log(1e-300) = 169,240 nats, putting the strip integral at
  # 196,000-203,000, while the four-corner form lands at 334,000-432,000. Any
  # value under 2.5e5 could only come from the strip form.
  path <- system.file("extdata", "export_dense_all.csv",
                      package = "rpbnb.tmb", mustWork = TRUE)
  d <- utils::read.csv(path)
  dat <- data.frame(y1 = d$ALL_3, y2 = d$C_DISTR)
  expect_gt(max(dat$y1), 200)  # not vacuous: the tail must actually be long

  fit <- fit_rpbnb_tmb(y1 ~ 1, y2 ~ 1, data = dat,
                       dependence = copula("normal"), draws = 2, keep = "full",
                       control = rpbnb_tmb_control(iterlim = 1, n_cores = 1))

  for (z_dep in c(-0.55, 0.2, 0.7)) {
    par <- c(0, 0, log(0.5), log(0.5), z_dep)
    nll <- as.numeric(fit$obj$fn(par))
    expect_true(is.finite(nll), label = paste("z_dep =", z_dep))
    expect_true(all(is.finite(fit$obj$gr(par))),
                label = paste("z_dep =", z_dep))
    expect_lt(nll, 2.5e5)
  }
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
