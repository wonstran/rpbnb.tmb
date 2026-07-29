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
  # evaluated on a five-panel 16-point Gauss-Legendre rule. The reference here
  # is stats::integrate() on the same integral -- adaptive Gauss-Kronrod, a
  # different algorithm entirely -- so this checks the rule is fine enough,
  # not merely that two copies of the same code agree.
  #
  # This test covers moderate dependence only; see the high-correlation test
  # below for the regime where the panel layout is what earns its keep.
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

test_that("Gaussian stays accurate at strong correlation", {
  skip_on_cran()
  # Regression for the defect that one fixed rule across the whole quantile
  # interval could not see. The integrand is a smoothed indicator of
  # [q(b')/rho, q(b)/rho] whose edges have width sqrt(1-rho^2)/|rho|; as
  # |rho| -> 1 that collapses while the interval does not, so the nodes step
  # over the plateau entirely. The previous rule returned 1.7e-25 for a cell
  # whose true probability is 0.0401 and lost 53.8 log-likelihood units on a
  # single observation.
  #
  # rho = 0.9999 is INSIDE the supported domain: .classify_boundary() in
  # R/inference.R marks a Gaussian estimate as a boundary only once tanh() has
  # actually saturated, so this is reported as an ordinary interior estimate.
  #
  # The margins are built comonotonically so the cells stay on the diagonal and
  # keep comfortably representable probabilities even at rho = 0.9999 -- a
  # reference sitting near the rounding floor would prove nothing, which is why
  # the moderate-dependence test above stops at z_dep = 1.
  n <- 60
  m1 <- 0.4; m2 <- 0.5
  u <- seq(0.03, 0.97, length.out = n)
  mu1 <- 3; mu2 <- 6
  dat <- data.frame(
    y1 = stats::qnbinom(u, size = 1 / m1, mu = mu1),
    y2 = stats::qnbinom(u, size = 1 / m2, mu = mu2)
  )
  fit <- fit_rpbnb_tmb(y1 ~ 1, y2 ~ 1, data = dat,
                       dependence = copula("normal"), draws = 2, keep = "full",
                       control = rpbnb_tmb_control(iterlim = 1, n_cores = 1))

  safe_qnorm <- function(p) stats::qnorm(pmin(pmax(p, 1e-15), 1 - 1e-15))
  cell <- function(qa, qam, qb, qbm, rho) {
    if (qa <= qam) return(0)
    s <- sqrt(1 - rho^2)
    stats::integrate(
      function(z) stats::dnorm(z) *
        (stats::pnorm((qb - rho * z) / s) - stats::pnorm((qbm - rho * z) / s)),
      lower = qam, upper = qa, rel.tol = 1e-12, subdivisions = 2000L
    )$value
  }
  qa  <- safe_qnorm(stats::pnbinom(dat$y1,     size = 1 / m1, mu = mu1))
  qam <- safe_qnorm(stats::pnbinom(dat$y1 - 1, size = 1 / m1, mu = mu1))
  qb  <- safe_qnorm(stats::pnbinom(dat$y2,     size = 1 / m2, mu = mu2))
  qbm <- safe_qnorm(stats::pnbinom(dat$y2 - 1, size = 1 / m2, mu = mu2))

  for (z_dep in c(2, 3, 4, 5)) {
    rho <- tanh(z_dep)
    p <- mapply(cell, qa, qam, qb, qbm, MoreArgs = list(rho = rho))
    # Guard the premise: a reference near the rounding floor proves nothing.
    expect_gt(min(p), 1e-8)
    par <- c(log(mu1), log(mu2), log(m1), log(m2), z_dep)
    expect_equal(as.numeric(fit$obj$fn(par)), -sum(log(p)),
                 tolerance = 1e-5,
                 label = sprintf("z_dep = %g (rho = %.6f)", z_dep, rho))
  }
})

test_that("Clayton reaches its comonotonic limit at strong dependence", {
  skip_on_cran()
  # Regression for the 1e15 cap on the factored ratios. The cap was applied
  # before the result raises the ratio to the power -1/theta, so it discarded
  # |k| * (log x - 34.5) in the exponent: at z_dep = 5 the exact symmetric cell
  # is 0.29077 and the capped form returned 0.15036. Clayton tends to the
  # comonotonic copula as theta grows, so this cell must APPROACH P(Y = 1);
  # the cap instead drove it to 2.4e-10 at z_dep = 20, reversing the
  # likelihood's behaviour across a region R/inference.R reports as interior.
  mu <- 1; m <- 0.5; sz <- 1 / m
  dat <- data.frame(y1 = 1L, y2 = 1L)
  fit <- fit_rpbnb_tmb(y1 ~ 1, y2 ~ 1, data = dat,
                       dependence = copula("kimeldorf"), draws = 2,
                       keep = "full",
                       control = rpbnb_tmb_control(iterlim = 1, n_cores = 1))
  A  <- stats::pnbinom(1, size = sz, mu = mu)
  Am <- stats::pnbinom(0, size = sz, mu = mu)
  clayton <- function(u, v, th) (u^(-th) + v^(-th) - 1)^(-1 / th)
  got <- function(z_dep) {
    exp(-as.numeric(fit$obj$fn(c(log(mu), log(mu), log(m), log(m), z_dep))))
  }

  # Exact four-corner comparison where double precision still supports one.
  for (z_dep in c(1, 2, 5)) {
    th <- exp(z_dep)
    exact <- clayton(A, A, th) - 2 * clayton(Am, A, th) + clayton(Am, Am, th)
    expect_gt(exact, 1e-8)          # guard the reference
    expect_equal(got(z_dep), exact, tolerance = 1e-10,
                 label = sprintf("z_dep = %g", z_dep))
  }

  # Beyond z_dep ~ 7 the four-corner reference cancels to exactly 0 in double
  # precision, so the limit itself becomes the only usable check.
  comonotonic <- stats::dnbinom(1, size = sz, mu = mu)
  expect_equal(got(20), comonotonic, tolerance = 1e-7)
  # ... and it must be approached from below, not crossed or collapsed.
  seq_p <- vapply(c(5, 8, 11, 14, 17, 20), got, numeric(1))
  expect_true(all(diff(seq_p) > 0))
  expect_true(all(seq_p < comonotonic))
})

test_that("Clayton's taped objective does not depend on the starting values", {
  skip_on_cran()
  # Regression for axis-branch selection. The branches were chosen by
  # asDouble(am) == 0.0, which resolves when the tape is BUILT, on a
  # parameter-dependent CDF that underflows to exactly zero for positive counts
  # in ordinary regions -- P(Y <= 0) at mu = r = 2000 is 0.5^2000. A tape built
  # at such a start kept the y = 0 axis formula permanently, so the same
  # parameter vector scored 2.256 or 1.592 nats depending only on where its
  # tape had been made. They are now selected from the observed counts.
  mk <- function(mu_tape, m_tape) {
    data <- .build_tmb_data(
      Y1 = 1, Y2 = 1,
      X1 = cbind("(Intercept)" = 1), X2 = cbind("(Intercept)" = 1),
      rand_idx1 = integer(0), rand_idx2 = integer(0),
      Z1 = matrix(0.5, 1L, 1L), Z2 = matrix(0.5, 1L, 1L),
      dist1 = 0L, dist2 = 0L, sign1 = 1L, sign2 = 1L,
      family_code = 3L, pois1 = FALSE, pois2 = FALSE,
      lamLo = 0, lamHi = 0, est_method = 0L)
    pars <- list(beta1 = log(mu_tape), beta2 = log(1),
                 log_sd1 = numeric(0), log_sd2 = numeric(0),
                 log_m1 = log(m_tape), log_m2 = log(0.5), z_dep = 0,
                 u1 = matrix(0, 1L, 0L), u2 = matrix(0, 1L, 0L))
    .make_rpbnb_tmb_object(data = data, parameters = pars,
                           map = list(), n_cores = 1L)$obj
  }
  # The extreme tape point is chosen so the first margin's lower CDF underflows.
  expect_identical(stats::pnbinom(0, size = 2000, mu = 2000), 0)

  target <- c(beta1 = log(1), beta2 = log(1),
              log_m1 = log(0.5), log_m2 = log(0.5), z_dep = 0)
  default_tape <- mk(1, 0.5)$fn(target)
  extreme_tape <- mk(2000, 1 / 2000)$fn(target)
  expect_equal(as.numeric(extreme_tape), as.numeric(default_tape),
               tolerance = 1e-12)
})

test_that("Clayton stays finite and accurate where the CDF saturates", {
  skip_on_cran()
  # Regression for the defect the strong-dependence rewrite introduced. Once
  # the NB2 CDF saturates, log F(y) and log F(y - 1) are equal to the last
  # bit, so building log(u/u') as their difference gives log(0) = -Inf and
  # NaN derivatives. On the truck margins at mu = 1 that is every count from
  # y = 70 up -- 197 of them -- and the difference is already 32% wrong at
  # y = 69 before it collapses. The ratio is therefore built from the marginal
  # mass, which nb2_cdf_pair() carries to full relative precision.
  #
  # This checks the GRADIENT and HESSIAN, not just the value: the -Inf left
  # the objective finite and only poisoned the derivatives, which is why the
  # existing value-only tail test passed while the Laplace fit died at its
  # first inner Newton step.
  counts <- c(0, 1, 5, 40, 69, 70, 100, 180, 266)
  dat <- data.frame(y1 = counts, y2 = rev(counts))   # mixed, incl. deep tail
  fit <- fit_rpbnb_tmb(y1 ~ 1, y2 ~ 1, data = dat,
                       dependence = copula("kimeldorf"), draws = 2,
                       keep = "full",
                       control = rpbnb_tmb_control(iterlim = 1, n_cores = 1))

  for (z_dep in c(-1, 0, 1, 3)) {
    par <- c(0, 0, log(0.5), log(0.5), z_dep)   # mu = 1, m = 0.5
    expect_true(is.finite(as.numeric(fit$obj$fn(par))),
                label = sprintf("value at z_dep = %g", z_dep))
    expect_true(all(is.finite(as.numeric(fit$obj$gr(par)))),
                label = sprintf("gradient at z_dep = %g", z_dep))
    expect_true(all(is.finite(as.numeric(fit$obj$he(par)))),
                label = sprintf("hessian at z_dep = %g", z_dep))
  }

  # Value check in the deep tail against the leading term of the Taylor
  # expansion, which is independent of the bracket algebra being tested:
  #
  #   (1+x+y)^k - (1+x)^k - (1+y)^k + 1 = k(k-1)xy + O(3)
  #
  # This is the check that matters most here, because the two terms of the
  # return are -k*xy and k^2*xy -- the same order. An implementation that
  # loses either one still produces a positive, smooth, finite cell
  # probability, wrong by a factor of k/(k-1). Writing du as
  # log(1+x+y) - log1p(x) - log1p(y) does exactly that: all three are O(x+y)
  # and their difference is O(xy), so it underflows to zero.
  #
  # A whole-formula R reference is deliberately NOT used. Written with x and y
  # in linear space it goes negative from y = 40 on, and written in log space
  # it just restates the implementation. The asymptotic is a different
  # expression arrived at a different way.
  lcdf3 <- function(y, mu, r) {
    lp <- log(r) - log(r + mu); lq <- log(mu) - log(r + mu)
    lt <- r * lp; lcum <- lt; lm1 <- lt
    for (k in seq_len(y)) {
      lm1 <- lcum
      lt <- lt + log(r + k - 1) - log(k) + lq
      m <- max(lcum, lt); lcum <- m + log1p(exp(min(lcum, lt) - m))
    }
    c(lcum, lm1, lt)
  }
  asym_logp <- function(y1, y2, th) {
    A <- lcdf3(y1, 1, 2); B <- lcdf3(y2, 1, 2); k <- -1 / th
    La <- -th * A[1]; Lb <- -th * B[1]; M <- max(La, Lb)
    ls00 <- M + log(exp(La - M) + exp(Lb - M) - exp(-M))
    lx <- function(v) {
      LR <- log1p(exp(v[3] - v[2]))
      -th * v[2] + log(-expm1(-th * LR)) - ls00
    }
    k * ls00 + log(k * (k - 1)) + lx(A) + lx(B)
  }
  for (z_dep in c(-1, 1)) {
    for (y in c(69, 100, 180, 266)) {
      one <- data.frame(y1 = y, y2 = y)
      f1 <- fit_rpbnb_tmb(y1 ~ 1, y2 ~ 1, data = one,
                          dependence = copula("kimeldorf"), draws = 2,
                          keep = "full",
                          control = rpbnb_tmb_control(iterlim = 1,
                                                      n_cores = 1))
      got <- -as.numeric(f1$obj$fn(c(0, 0, log(0.5), log(0.5), z_dep)))
      ref <- asym_logp(y, y, exp(z_dep))
      # Not vacuous only if the cell really is in the saturated tail.
      expect_lt(ref, -100)
      expect_equal(got, ref, tolerance = 1e-9,
                   label = sprintf("log p at y = %d, z_dep = %g", y, z_dep))
    }
  }
})
