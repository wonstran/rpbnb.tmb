test_that("the shared dependence link reproduces the reported natural values", {
  z <- 0.37
  lo <- -0.4
  hi <- 0.9
  coefs <- c(log_m1 = log(0.5), log_m2 = log(0.4), z_dep = z)

  for (fc in 0:3) {
    link <- .rpbnb_dependence_link(fc, lamLo = lo, lamHi = hi)
    report <- .rpbnb_natural_report(coefs, family_code = fc,
                                    lamLo = lo, lamHi = hi)
    values <- link$map(z)

    expect_identical(names(values), link$names, info = paste("family", fc))
    for (nm in link$names) {
      expect_equal(report[[nm]]$value, unname(values[[nm]]),
                   info = paste("family", fc, nm))
    }
  }
})

test_that("independence has no dependence link", {
  expect_null(.rpbnb_dependence_link(-1L))
})

dp_test_data <- function(n = 200L) {
  set.seed(250725)
  x <- stats::runif(n, -1, 1)
  shared <- stats::rnorm(n, 0, 0.3)
  data.frame(
    y1 = stats::rpois(n, exp(0.4 + 0.3 * x + shared)),
    y2 = stats::rpois(n, exp(0.3 - 0.2 * x + shared)),
    x  = x
  )
}

test_that("the working-scale profile brackets the working estimate", {
  fit <- fit_rpbnb_tmb(
    y1 ~ x, y2 ~ x, dp_test_data(),
    dependence = copula("normal"), draws = 1L, keep = "full"
  )

  ci <- .dependence_profile_ci(fit, level = 0.95)

  expect_false(is.null(ci))
  expect_length(as.numeric(ci), 2L)
  expect_false(is.null(attr(ci, "profile")))
  z_hat <- unname(coef(fit)[["z_dep"]])
  expect_lt(as.numeric(ci)[1], z_hat)
  expect_gt(as.numeric(ci)[2], z_hat)
})

test_that("a fit without the TMB objective yields NULL, not an error", {
  fit <- fit_rpbnb_tmb(
    y1 ~ x, y2 ~ x, dp_test_data(),
    dependence = copula("normal"), draws = 1L, keep = "postfit"
  )

  expect_null(fit$obj)
  expect_null(.dependence_profile_ci(fit, level = 0.95))
})

test_that("a caller-supplied trace does not silently degrade the profile to NULL", {
  fit <- fit_rpbnb_tmb(
    y1 ~ x, y2 ~ x, dp_test_data(),
    dependence = copula("normal"), draws = 1L, keep = "full"
  )

  # trace is a real, documented TMB::tmbprofile() argument, and this
  # function's own docs promise `...` is forwarded to it. Before the fix,
  # supplying it collided with the hardcoded `trace = FALSE` inside
  # do.call(), which threw "formal argument matched by multiple actual
  # arguments" -- caught by the inner try() and converted to NULL, so the
  # caller silently fell back to a Wald interval with no diagnostic.
  ci_trace_true <- .dependence_profile_ci(fit, level = 0.95, trace = TRUE)
  expect_false(is.null(ci_trace_true))
  expect_length(as.numeric(ci_trace_true), 2L)

  ci_trace_false <- .dependence_profile_ci(fit, level = 0.95, trace = FALSE)
  expect_false(is.null(ci_trace_false))
  expect_length(as.numeric(ci_trace_false), 2L)
})

test_that("the pre-call last.par.best is restored after profiling, even when it differs from optimizer$par", {
  fit <- fit_rpbnb_tmb(
    y1 ~ x, y2 ~ x, dp_test_data(),
    dependence = copula("normal"), draws = 1L, keep = "full"
  )

  # obj is fit$obj's environment, not a copy: .dependence_profile_ci()
  # overwrites obj$env$last.par.best with fit$optimizer$par before profiling.
  # TMB::tmbprofile() restores that back to fit$optimizer$par on its own
  # on.exit, but nothing restores it to whatever it held *before*
  # .dependence_profile_ci() was called -- so any caller-visible state
  # predating this call is lost. A freshly-fitted `fit$obj$env$last.par.best`
  # already equals fit$optimizer$par, which would let a missing fix pass this
  # test by coincidence -- force a sentinel value that differs from
  # fit$optimizer$par so the assertion actually exercises the restore path.
  sentinel <- fit$optimizer$par + 1
  fit$obj$env$last.par.best <- sentinel

  ci <- .dependence_profile_ci(fit, level = 0.95)
  expect_false(is.null(ci))
  expect_identical(fit$obj$env$last.par.best, sentinel)
})
