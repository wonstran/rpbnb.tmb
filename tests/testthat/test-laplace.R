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
  # max_workload is set explicitly rather than left to the default. The default
  # auto-detects available memory, and .detect_available_memory_gib() is
  # documented to degrade to NA on platforms where detection is unavailable --
  # a sandbox that blocks subprocesses, an unrecognized OS -- at which point
  # rpbnb_tmb_max_workload() warns and falls back to 8 GiB. That is supported
  # behaviour, so leaving the default here would make expect_no_warning()
  # reject the package's own documented fallback on restricted Windows/CI
  # runners. This test is about draws surviving for post-estimation averaging,
  # not about memory detection; pinning the budget keeps the warning guard
  # meaningful for warnings that would actually indicate a defect.
  expect_no_warning(
    fit <- fit_rpbnb_tmb(y1 ~ x1, y2 ~ x1, data = dat, random_1 = "x1",
                         dependence = "independence", draws = 25,
                         method = "laplace",
                         control = rpbnb_tmb_control(max_workload = Inf))
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

  expect_identical(fit_sml$optimizer$convergence, 0L)
  expect_identical(fit_lap$optimizer$convergence, 0L)

  # Agreement judged against sampling noise rather than an absolute tolerance:
  # these are two approximations to the same integral, not two computations of
  # the same number.
  # log_sd1:x1 is included deliberately: it is the parameter that most
  # directly reflects correctness of the random-effect integration (Jacobian,
  # standard-normal prior, latent-to-deviation scale). b2:x1 belongs to an
  # equation with no random effect and barely discriminates on its own.
  key <- c("b1:x1", "b2:x1", "log_sd1:x1")
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

test_that("sml and laplace degrade together on 8 random slopes, not just 1", {
  skip_on_cran()
  # A review of this branch flagged that the accuracy evidence for Laplace
  # used ONE latent per observation, while the truck acceptance fit uses
  # EIGHT (4 random slopes per equation), and asked whether Laplace might
  # specifically flatten log_sd in that many-latent regime -- i.e. whether
  # the two estimators would disagree once there are enough random slopes for
  # a per-observation sparse-Hessian approximation to behave differently from
  # averaging over Halton draws.
  #
  # They do not disagree: both estimators locate the same scales, and both
  # collapse on x4 (whose true sd, 0.30, is the smallest of the four and is
  # the least identified from n = 600 observations). That collapse is a
  # property of the likelihood/data, not a flaw of either approximation, so
  # this test's job is to show agreement despite it, and the minimum-headroom
  # assertion below exists so the test cannot pass merely because every
  # tolerance is enormous.
  sim <- simulate_rpbnb_tmb(
    n = 600,
    beta1 = c("(Intercept)" = 0.5, x1 = 0.30, x2 = -0.25, x3 = 0.20, x4 = 0.15),
    beta2 = c("(Intercept)" = 0.3, x1 = -0.20, x2 = 0.25, x3 = -0.15, x4 = 0.20),
    random_1 = list(x1 = list(sd = 0.45), x2 = list(sd = 0.35),
                    x3 = list(sd = 0.40), x4 = list(sd = 0.30)),
    random_2 = list(x1 = list(sd = 0.45), x2 = list(sd = 0.35),
                    x3 = list(sd = 0.40), x4 = list(sd = 0.30)),
    dispersion = c(m1 = 0.4, m2 = 0.5),
    dependence = "independence", seed = 4242
  )
  common <- list(formula_1 = y1 ~ x1 + x2 + x3 + x4,
                 formula_2 = y2 ~ x1 + x2 + x3 + x4,
                 data = sim$data,
                 random_1 = c("x1", "x2", "x3", "x4"),
                 random_2 = c("x1", "x2", "x3", "x4"),
                 dependence = "independence",
                 # The acceptance evidence this codifies was run at draws =
                 # 200 (SML took 226s there); halved to 100 here to keep the
                 # suite's runtime down. Agreement was re-checked directly at
                 # draws = 100 before writing this test.
                 draws = 100, seed = 4242)

  fit_sml <- do.call(fit_rpbnb_tmb, c(common, list(method = "sml")))
  fit_lap <- do.call(fit_rpbnb_tmb, c(common, list(method = "laplace")))

  scales <- c("log_sd1:x1", "log_sd1:x2", "log_sd1:x3", "log_sd1:x4",
              "log_sd2:x1", "log_sd2:x2", "log_sd2:x3", "log_sd2:x4")

  # Same tolerance idiom as "laplace and sml agree on data simulated from the
  # true process" above: agreement judged against sampling noise, since these
  # are two approximations to the same integral rather than two computations
  # of the same number.
  for (nm in scales) {
    tol <- 3 * max(fit_sml$se[[nm]], fit_lap$se[[nm]])
    expect_lt(abs(coef(fit_sml)[[nm]] - coef(fit_lap)[[nm]]), tol)
  }

  # Not vacuous: if every scale's tolerance were enormous (as x4's is, both
  # estimators' se(log_sd1:x4) run past 1), the loop above would pass by
  # construction. Require that most scales (x1-x3 in both equations) are
  # actually identified with a reasonably tight standard error, so the
  # agreement demonstrated above is evidence of real correctness rather than
  # an artifact of loose tolerances.
  precise <- vapply(scales, function(nm) fit_sml$se[[nm]] < 0.5, logical(1))
  expect_gte(sum(precise), 4L)
})

test_that("laplace fits every copula family", {
  skip_on_cran()
  # Every Laplace test above uses famoye or independence, whose margins are
  # dnbinom2. The copula families instead need the NB2 *CDF*, and the outer
  # Laplace gradient differentiates the joint negative log-likelihood three
  # times. Taking that CDF from TMB's pbeta() made all three copula families
  # fail here -- nlminb stopped one or two steps in with "inner newton
  # optimization failed during gradient calculation" and an NA/NaN gradient --
  # because pbeta()'s third derivatives are NaN over wide regions of ordinary
  # parameter values. See tests/testthat/test-fit-copula.R for the companion
  # check that the replacement reproduces the NB2 mass function exactly.
  sim <- simulate_rpbnb_tmb(
    n = 300,
    beta1 = c("(Intercept)" = 0.3, x1 = 0.35),
    beta2 = c("(Intercept)" = 0.1, x1 = -0.25),
    random_1 = list(x1 = list(sd = 0.35)),
    dispersion = c(m1 = 0.4, m2 = 0.5),
    dependence = copula("frank", par = 3), seed = 31
  )
  for (fam in c("frank", "normal", "kimeldorf")) {
    fit <- fit_rpbnb_tmb(y1 ~ x1, y2 ~ x1, data = sim$data, random_1 = "x1",
                         dependence = copula(fam), draws = 50, seed = 5,
                         method = "laplace",
                         control = rpbnb_tmb_control(n_cores = 1))
    expect_true(is.finite(fit$logLik), label = fam)
    expect_identical(fit$optimizer$convergence, 0L, label = fam)
    expect_true(all(is.finite(fit$se)), label = fam)
  }
})

test_that("laplace gives the same answer single- and multi-threaded", {
  skip_on_cran()
  # Mirrors tests/testthat/test-parallel.R: on a runtime whose TMB build
  # supports only one thread, requesting 4 cores below silently realizes 1,
  # which would make this a red test rather than a skip.
  previous <- TMB::openmp(DLL = "rpbnb.tmb")
  on.exit(TMB::openmp(n = previous, DLL = "rpbnb.tmb"), add = TRUE)
  supported <- as.integer(TMB::openmp(max = TRUE, DLL = "rpbnb.tmb")[[1L]])
  if (supported < 2L) skip("TMB runtime supports only one thread")

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

  # Without this the test is a tautology: if OpenMP is unavailable, fit4 runs
  # single-threaded and the comparisons below compare a fit to itself.
  expect_gt(fit4$parallel$realized, 1L)

  # Threading partitions the accumulation; it must not change the objective.
  expect_equal(fit1$logLik, fit4$logLik, tolerance = 1e-6)
  expect_equal(unname(coef(fit1)), unname(coef(fit4)), tolerance = 1e-4)
})
