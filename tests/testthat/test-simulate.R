simulation_coefficients <- function() {
  list(
    beta1 = c("(Intercept)" = 0.2, x = 0.3),
    beta2 = c("(Intercept)" = 0.1, x = -0.2)
  )
}

test_that("copula simulation ignores the Famoye-only lambda argument", {
  beta <- simulation_coefficients()
  dep <- copula("kimeldorf", par = 5)

  default_lambda <- simulate_rpbnb_tmb(
    n = 300,
    beta1 = beta$beta1,
    beta2 = beta$beta2,
    dependence = dep,
    seed = 2718
  )
  explicit_lambda <- simulate_rpbnb_tmb(
    n = 300,
    beta1 = beta$beta1,
    beta2 = beta$beta2,
    dependence = dep,
    lambda = 1,
    seed = 2718
  )

  expect_identical(default_lambda$data, explicit_lambda$data)
  expect_gt(cor(default_lambda$data$y1, default_lambda$data$y2,
                method = "kendall"), 0.25)
})

test_that("simulation rejects unrecognized dependence values directly", {
  beta <- simulation_coefficients()

  expect_error(
    simulate_rpbnb_tmb(
      n = 10,
      beta1 = beta$beta1,
      beta2 = beta$beta2,
      dependence = "garbage",
      seed = 1
    ),
    '`dependence` must be "famoye", "independence", or copula\\(\\)'
  )
})

test_that("Gaussian copula simulation needs no bivariate-normal CDF", {
  # Conditional inversion of a Gaussian copula is one univariate normal draw,
  # so this path no longer calls pbivnorm::pbivnorm() and pbivnorm was dropped
  # from Suggests. It must work whether or not that package is installed.
  beta <- simulation_coefficients()

  expect_no_error(
    simulate_rpbnb_tmb(
      n = 10,
      beta1 = beta$beta1,
      beta2 = beta$beta2,
      dependence = copula("normal", par = 0.5),
      seed = 1
    )
  )
  expect_false(any(grepl("pbivnorm", deparse(simulate_rpbnb_tmb))))
})

test_that("Famoye conditional sampling inverts the full CDF, tail included", {
  # Deterministic, and it is the assertion that actually discriminates. A
  # Monte Carlo check on the mean does not: at mu = 5, dispersion 10 the old
  # truncated sampler was 0.205 low against a 0.25 five-percent band, so it
  # would have passed. Invert the marginal NB through the same helper the
  # Famoye path uses and require the exact quantile, including past the 0.999
  # cutoff the old implementation could never return.
  mu <- 5
  for (m in c(10, 50)) {
    size <- 1 / m
    cap <- qnbinom(0.999, mu = mu, size = size)
    block <- function(yv) dnbinom(yv, mu = mu, size = size)

    for (u in c(0.01, 0.5, 0.9, 0.999, 0.9999, 0.999999)) {
      drawn <- rpbnb.tmb:::.sample_conditional_count(
        u = u, total_mass = 1, block_pmf = block,
        first_block = 8L,          # forces several block extensions
        fallback = function(u) qnbinom(u, mu = mu, size = size)
      )
      expect_identical(drawn, qnbinom(u, mu = mu, size = size))
    }
    # The far-tail draws must actually exceed the old truncation point,
    # otherwise the test above proves nothing about truncation.
    expect_gt(qnbinom(0.999999, mu = mu, size = size), cap)
  }
})

test_that("an independence-valued copula reproduces the requested NB2 margin", {
  # End-to-end companion to the deterministic test above. The discriminating
  # assertion is the tail one, not the mean: under the old sampler no y2 could
  # exceed qnbinom(0.999), by construction.
  mu <- 5
  m <- 10
  cap <- qnbinom(0.999, mu = mu, size = 1 / m)
  sim <- simulate_rpbnb_tmb(
    n = 20000,
    beta1 = c("(Intercept)" = log(mu), x = 0),
    beta2 = c("(Intercept)" = log(mu), x = 0),
    dispersion = c(m1 = m, m2 = m),
    dependence = copula("frank"),   # par omitted => independence-valued
    seed = 4242
  )

  # ~0.1% of 20000 draws belong above the old cap; seeing any is impossible
  # under truncation, and seeing zero here would be a ~2e-9 fluke.
  expect_gt(sum(sim$data$y2 > cap), 0L)
  # And the margin is unbiased, which truncation also broke (4.795 vs 5).
  expect_equal(mean(sim$data$y2), mu, tolerance = 0.05)
})

test_that("untruncated sampling still reproduces requested dependence", {
  beta <- simulation_coefficients()
  sim <- simulate_rpbnb_tmb(
    n = 4000, beta1 = beta$beta1, beta2 = beta$beta2,
    dependence = copula("kimeldorf", par = 5), seed = 991
  )
  expect_gt(cor(sim$data$y1, sim$data$y2, method = "kendall"), 0.5)
})

test_that("copula conditional inversion round-trips its conditional CDF", {
  # Deterministic and exact -- this is the assertion that would have caught
  # the strong-Frank failure. Each family's conditional CDF h(v|u) = dC/du is
  # written here in a cancellation-free form and applied to the sampler's
  # output; it must return the driving uniform back. Frank is exercised to
  # +/-FRANK_THETA_MAX, the ceiling the fitter supports, because the failure
  # was a precision loss that only appears once the CDF saturates.
  inv <- rpbnb.tmb:::.copula_conditional_inverse
  grid <- expand.grid(u = c(0.01, 0.2, 0.5, 0.8, 0.99),
                      w = c(0.01, 0.25, 0.5, 0.75, 0.99))

  h_frank <- function(theta, u, v) {
    (exp(-theta * (u + v)) - exp(-theta * u)) /
      (exp(-theta) + exp(-theta * (u + v)) -
         exp(-theta * u) - exp(-theta * v))
  }
  for (theta in c(-FRANK_THETA_MAX, -30, -10, -1, 1, 10, 20, 30,
                  FRANK_THETA_MAX)) {
    v <- inv("frank", theta, grid$u, grid$w)
    expect_equal(h_frank(theta, grid$u, v), grid$w, tolerance = 1e-10)
    expect_true(all(v > 0 & v < 1))
  }

  # Only up to theta = 50: past that the reference formula itself overflows on
  # u^(-theta), which is exactly why the implementation is in log space. The
  # domain above 50 is covered by the uniformity test below instead.
  h_clayton <- function(theta, u, v) {
    u^(-theta - 1) * (u^(-theta) + v^(-theta) - 1)^(-1 / theta - 1)
  }
  for (theta in c(0.5, 5, 20, 50)) {
    v <- inv("kimeldorf", theta, grid$u, grid$w)
    expect_equal(h_clayton(theta, grid$u, v), grid$w, tolerance = 1e-10)
  }

  h_normal <- function(rho, u, v) {
    pnorm((qnorm(v) - rho * qnorm(u)) / sqrt(1 - rho^2))
  }
  for (rho in c(-0.9, -0.5, 0.5, 0.9, 0.99)) {
    v <- inv("normal", rho, grid$u, grid$w)
    expect_equal(h_normal(rho, grid$u, v), grid$w, tolerance = 1e-10)
  }

  # theta = 0 must be exactly independence, not merely close to it.
  expect_identical(inv("frank", 0, grid$u, grid$w), grid$w)
  expect_identical(inv("kimeldorf", 0, grid$u, grid$w), grid$w)
})

test_that("the conditional inverse is uniform across the accepted domain", {
  # The bug this guards: the direct forms overflowed/underflowed inside the
  # domain copula() accepts. Frank's exp(-theta*u) dies past theta = -710 and
  # returned v identically 1 (so qnbinom gave Inf counts); Clayton's
  # u^(-theta) dies well before CLAYTON_THETA_MAX, and at the fitting-link
  # ceiling produced a point mass at zero instead of a uniform.
  #
  # These are the PUBLIC limits, not convenient interior values: copula()
  # accepts any finite Frank par and any positive Clayton par, and a fitted
  # Clayton model can report up to exp(20).
  inv <- rpbnb.tmb:::.copula_conditional_inverse
  set.seed(20260727)
  n <- 20000
  u <- runif(n)
  w <- runif(n)

  # The domain ends at +/-.Machine$double.xmax, not at a round number chosen
  # for comfort. Clayton at double.xmax was the last thing to break here:
  # theta * log(u) overflowed for every u < exp(-1), so 37% of draws collapsed
  # to the lower clamp.
  big <- .Machine$double.xmax
  # copula() accepts |rho| < 1, so the Gaussian boundary is the closest
  # double to +/-1 on either side, not a round 1 - 1e-12. Same principle as
  # `big` for the unbounded families: the list has to be the contract, not a
  # number that merely looks extreme.
  near_one <- 1 - .Machine$double.eps / 2
  domain <- list(
    list("frank", c(-big, -1e300, -1e100, -1e6, -1000, -750, -710,
                    -FRANK_THETA_MAX, FRANK_THETA_MAX, 710, 750, 1000,
                    1e6, 1e100, 1e300, big)),
    list("kimeldorf", c(0.5, 20, 500, 1000, 1e5, CLAYTON_THETA_MAX,
                        1e12, 1e100, 1e300, big)),
    list("normal", c(-near_one, -1 + 1e-12, -0.9, 0.9, 1 - 1e-12, near_one))
  )
  # Guard the guard: these must actually be inside what copula() accepts, or
  # the list is testing something the API does not promise.
  expect_no_error(copula("normal", par = near_one))
  expect_no_error(copula("normal", par = -near_one))
  expect_error(copula("normal", par = 1))

  eps <- 1e-15  # the implementation's own clamp bound
  for (spec in domain) {
    for (theta in spec[[2L]]) {
      label <- paste0(spec[[1L]], " theta=", format(theta))
      v <- inv(spec[[1L]], theta, u, w)

      expect_true(all(is.finite(v)), info = label)
      # `> 0 & < 1` alone cannot see a numerical collapse, because the
      # implementation clamps to [eps, 1-eps] immediately before this runs --
      # a degenerate result arrives here already looking interior. Require
      # that NOTHING sits on the clamp bound, so the clamp is provably not
      # load-bearing rather than merely asserted to be.
      expect_false(any(v <= eps | v >= 1 - eps), info = label)
      # And actually uniform -- being finite and interior is not enough.
      expect_equal(mean(v), 0.5, tolerance = 0.02, info = label)
      expect_lt(suppressWarnings(ks.test(v, "punif")$statistic), 0.02)
    }
  }
})

test_that("the conditional inverse approaches its comonotonic limits", {
  # Uniformity alone would also hold for an inverse that ignored u entirely,
  # so pin the dependence too: as theta grows the copula must converge to
  # v = u, and Frank's negative branch to v = 1 - u.
  inv <- rpbnb.tmb:::.copula_conditional_inverse
  set.seed(31)
  u <- runif(5000)
  w <- runif(5000)

  expect_lt(max(abs(inv("frank", 1e6, u, w) - u)), 1e-4)
  expect_lt(max(abs(inv("frank", -1e6, u, w) - (1 - u))), 1e-4)
  expect_lt(max(abs(inv("kimeldorf", CLAYTON_THETA_MAX, u, w) - u)), 1e-4)

  # At the top of the domain, distinguishing bitwise identity from
  # equal-within-tolerance, because they are not the same claim:
  big <- .Machine$double.xmax

  # Bitwise. Clayton's branch is u * exp(-c/theta); c/theta is about 1e-308
  # here -- not zero, the smallest subnormal is sixteen orders lower -- but far
  # below half an ulp at 1, so exp() rounds to exactly 1 and u * 1 is u.
  # Frank's positive branch reduces to u + (log w - log1p(-w))/theta, whose
  # second term is likewise below ulp(u). Both are consequences of IEEE
  # rounding, not luck -- and this is the assertion the earlier Clayton defect
  # failed outright (max|v-u| was 0.37).
  expect_identical(inv("kimeldorf", big, u, w), u)
  expect_identical(inv("frank", big, u, w), u)

  # NOT bitwise: Frank's negative branch differences two quantities of order
  # 1e308 before dividing, so the result carries one ulp of relative error --
  # 3748 of 5000 values differ from 1-u, by at most eps/2. State the bound
  # rather than leaning on expect_equal()'s default tolerance.
  v_neg <- inv("frank", -big, u, w)
  expect_lte(max(abs(v_neg - (1 - u))), .Machine$double.eps)
})

test_that("extreme copula parameters simulate finite counts", {
  # End-to-end companion: theta = -710 used to yield 1000 Inf values out of
  # 1000 draws.
  b <- function(mu) c("(Intercept)" = log(mu), x = 0)
  cases <- list(
    list("frank", -1000), list("frank", -710), list("frank", 1000),
    list("kimeldorf", CLAYTON_THETA_MAX)
  )
  for (case in cases) {
    sim <- simulate_rpbnb_tmb(
      n = 1000, beta1 = b(5), beta2 = b(5),
      dispersion = c(m1 = 1, m2 = 1),
      dependence = copula(case[[1L]], par = case[[2L]]), seed = 2
    )
    label <- paste0(case[[1L]], " theta=", format(case[[2L]]))
    expect_true(all(is.finite(sim$data$y1)), info = label)
    expect_true(all(is.finite(sim$data$y2)), info = label)
    # Near-perfect rank agreement, in the sign the parameter asks for.
    tau <- cor(sim$data$y1, sim$data$y2, method = "kendall")
    expect_gt(sign(case[[2L]]) * tau, 0.9, label = label)
  }
})

test_that("strong Frank dependence simulates rather than aborting", {
  # Regression: enumerating copula rectangles at theta = 30 produced evaluated
  # masses from 0 to 70x the analytic marginal, and the sampler aborted with
  # "did not converge". Both of these calls used to error.
  b <- function(mu) c("(Intercept)" = log(mu), x = 0)
  for (m in c(10, 50)) {
    sim <- expect_no_error(simulate_rpbnb_tmb(
      n = 2000, beta1 = b(5), beta2 = b(5),
      dispersion = c(m1 = m, m2 = m),
      dependence = copula("frank", par = 30),
      seed = if (m == 10) 321 else 123
    ))
    expect_gt(cor(sim$data$y1, sim$data$y2, method = "kendall"), 0.5)
  }

  # And at the ceiling itself, in both directions.
  for (theta in c(-FRANK_THETA_MAX, FRANK_THETA_MAX)) {
    sim <- expect_no_error(simulate_rpbnb_tmb(
      n = 2000, beta1 = b(5), beta2 = b(5),
      dispersion = c(m1 = 1, m2 = 1),
      dependence = copula("frank", par = theta), seed = 8
    ))
    tau <- cor(sim$data$y1, sim$data$y2, method = "kendall")
    expect_gt(sign(theta) * tau, 0.7)
  }
})

test_that("strong dependence leaves both margins exactly NB2", {
  # The first margin is qnbinom() of a uniform, so it is exact by
  # construction; the second is only exact if the conditional inverse is
  # marginally uniform. Check the second directly and at high dependence.
  b <- function(mu) c("(Intercept)" = log(mu), x = 0)
  sim <- simulate_rpbnb_tmb(
    n = 30000, beta1 = b(5), beta2 = b(5),
    dispersion = c(m1 = 1, m2 = 1),
    dependence = copula("frank", par = 20), seed = 606
  )
  # SE of each mean is about 0.032 here, so 0.04 relative is a ~4 SE band.
  expect_equal(mean(sim$data$y1), 5, tolerance = 0.04)
  expect_equal(mean(sim$data$y2), 5, tolerance = 0.04)

  # Goodness of fit on the whole distribution, not just the mean. Chi-square,
  # not ks.test(): KS assumes a continuous reference and on count data it
  # reports D = 0.167, p = 0 for a perfectly correct sample.
  binned <- table(factor(pmin(sim$data$y2, 30L), levels = 0:30))
  expected <- c(dnbinom(0:29, mu = 5, size = 1),
                1 - pnbinom(29, mu = 5, size = 1))
  expect_gt(
    suppressWarnings(chisq.test(as.numeric(binned), p = expected)$p.value),
    0.001
  )
})

test_that("an intercept-only model simulates on every dependence path", {
  # Pre-existing: with no covariates `vars` is empty, and the covariate frame
  # was built 0x0, colliding with the n-row outcomes in the final
  # data.frame(). Covers the direct path and both conditional paths.
  only <- c("(Intercept)" = 0)
  for (dep in list("independence", "famoye", copula("frank", par = 5))) {
    sim <- expect_no_error(simulate_rpbnb_tmb(
      n = 5, beta1 = only, beta2 = only, dependence = dep,
      lambda = if (identical(dep, "famoye")) 0.3 else 0, seed = 1
    ))
    expect_identical(nrow(sim$data), 5L)
    expect_identical(names(sim$data), c("y1", "y2"))
  }
})
