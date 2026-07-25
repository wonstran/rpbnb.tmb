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
  # trace = TRUE makes tmbprofile() print a "Profile value: ..." line per
  # step; capture.output() swallows that progress spew without touching what
  # is actually being asserted (that a caller-supplied trace still profiles).
  invisible(capture.output(
    ci_trace_true <- .dependence_profile_ci(fit, level = 0.95, trace = TRUE)
  ))
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

test_that("the profile interval brackets the estimate on the natural scale", {
  fit <- fit_rpbnb_tmb(
    y1 ~ x, y2 ~ x, dp_test_data(),
    dependence = copula("normal"), draws = 1L, keep = "full"
  )

  pr <- rpbnb_tmb_dependence_profile(fit)

  expect_identical(pr$parameter, c("rho", "tau"))
  expect_identical(unique(pr$method), "profile")
  expect_true(all(pr$lower < pr$estimate))
  expect_true(all(pr$estimate < pr$upper))
  expect_true(all(pr$level == 0.95))
})

test_that("a wider level nests the narrower one", {
  fit <- fit_rpbnb_tmb(
    y1 ~ x, y2 ~ x, dp_test_data(),
    dependence = copula("normal"), draws = 1L, keep = "full"
  )

  narrow <- rpbnb_tmb_dependence_profile(fit, level = 0.95)
  wide   <- rpbnb_tmb_dependence_profile(fit, level = 0.99)

  expect_true(all(wide$lower <= narrow$lower))
  expect_true(all(wide$upper >= narrow$upper))
})

test_that("Famoye endpoints stay inside the frozen lambda bounds", {
  fit <- fit_rpbnb_tmb(
    y1 ~ x, y2 ~ x, dp_test_data(),
    dependence = "famoye", draws = 1L, keep = "full"
  )

  pr <- suppressWarnings(rpbnb_tmb_dependence_profile(fit))
  bounds <- fit$lambda_bounds

  expect_identical(pr$parameter, "lam")
  # Holds whether or not the optimum pinned: a pinned fit returns NA endpoints,
  # and every finite endpoint must lie in the box.
  ends <- c(pr$lower, pr$upper)
  ends <- ends[is.finite(ends)]
  expect_true(all(ends >= bounds[["lower"]]))
  expect_true(all(ends <= bounds[["upper"]]))
})

test_that("a fit without the TMB objective warns and degrades to Wald", {
  fit <- fit_rpbnb_tmb(
    y1 ~ x, y2 ~ x, dp_test_data(),
    dependence = copula("normal"), draws = 1L, keep = "postfit"
  )

  expect_warning(
    pr <- rpbnb_tmb_dependence_profile(fit),
    "falling back to a Wald interval"
  )
  expect_identical(unique(pr$method), "wald")
  expect_true(all(is.finite(pr$lower)))
  expect_true(all(pr$lower < pr$estimate))
  expect_null(attr(pr, "profile"))
})

test_that("a reloaded fit still profiles because TMB retapes from stored data", {
  fit <- fit_rpbnb_tmb(
    y1 ~ x, y2 ~ x, dp_test_data(),
    dependence = copula("normal"), draws = 1L, keep = "full"
  )
  path <- tempfile(fileext = ".rds")
  saveRDS(fit, path)
  reloaded <- readRDS(path)
  on.exit(unlink(path), add = TRUE)

  # The reloaded external pointer really is dead (reloaded$obj$env$ADFun$ptr
  # is a nil pointer) but that does not make the objective unusable: TMB
  # detects the nil pointer and retapes from the data and parameters it kept
  # in obj$env, so reloaded$obj$fn() returns the same value it did before
  # saving and the liveness probe in .dependence_profile_ci() passes. The
  # earlier version of this test assumed retaping doesn't happen and expected
  # a Wald downgrade; that was verified false by comparing
  # reloaded$obj$fn(reloaded$optimizer$par) against the pre-save value
  # (bit-identical, ruling out a dangling-pointer read of freed memory).
  # Retaping needs the package's compiled DLL loaded in the session, which
  # always holds while the package is in use.
  expect_no_warning(pr <- rpbnb_tmb_dependence_profile(reloaded))
  expect_identical(unique(pr$method), "profile")
})

test_that("independence has nothing to profile", {
  fit <- fit_rpbnb_tmb(
    y1 ~ x, y2 ~ x, dp_test_data(),
    dependence = "independence", draws = 1L, keep = "full"
  )

  expect_error(rpbnb_tmb_dependence_profile(fit), "nothing to profile")
})

test_that("level and fit are validated", {
  fit <- fit_rpbnb_tmb(
    y1 ~ x, y2 ~ x, dp_test_data(),
    dependence = copula("normal"), draws = 1L, keep = "full"
  )

  expect_error(rpbnb_tmb_dependence_profile(fit, level = 1), "between 0 and 1")
  expect_error(rpbnb_tmb_dependence_profile(fit, level = c(0.9, 0.95)),
               "between 0 and 1")
  expect_error(rpbnb_tmb_dependence_profile(list()), "rpbnb_tmb_fit")
})

test_that('inference = "none" returns the estimate with NA endpoints', {
  fit <- fit_rpbnb_tmb(
    y1 ~ x, y2 ~ x, dp_test_data(),
    dependence = copula("normal"), draws = 1L,
    keep = "postfit", inference = "none"
  )

  pr <- suppressWarnings(rpbnb_tmb_dependence_profile(fit))

  expect_true(all(is.finite(pr$estimate)))
  expect_true(all(is.na(pr$lower)))
  expect_true(all(is.na(pr$upper)))
})

test_that("the reported estimate matches the fit's own natural value", {
  for (dep in list("famoye", copula("frank"))) {
    fit <- fit_rpbnb_tmb(
      y1 ~ x, y2 ~ x, dp_test_data(),
      dependence = dep, draws = 1L, keep = "full"
    )
    pr <- suppressWarnings(rpbnb_tmb_dependence_profile(fit))
    reported <- summary(fit$sdreport, "report")

    for (nm in pr$parameter) {
      expect_equal(
        pr$estimate[pr$parameter == nm],
        unname(reported[nm, "Estimate"]),
        info = nm
      )
    }
  }
})
