# RUNTIME: this is the slowest file in the suite -- about 6.2 minutes for 59
# assertions on a 2026-era laptop (R 4.5.1, TMB 1.9.21, Windows), against ~12
# seconds for the next slowest. TMB::tmbprofile() refits the model at every
# point along the profile grid, and no amount of fixture sharing removes that.
#
# It is SLOW, NOT HUNG. Recorded here because that has now cost two reviewers
# real time: one abandoned a run at ~23 minutes, another timed out at 2 and
# then at 6 minutes -- the second missing completion by about twelve seconds
# and reporting the result as inconclusive rather than passing. Budget 10+
# minutes, or run it with a filter and no wall-clock limit.
#
# If you are shortening the feedback loop, `devtools::test(filter = "...")`
# on any other file returns in seconds; leave this one to a full run.

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

# Shared across most `.dependence_profile_ci()` / `rpbnb_tmb_dependence_profile()`
# tests below: fitting `copula("normal")` on `dp_test_data()` with
# `keep = "full"` is deterministic (dp_test_data() seeds internally) and was
# previously repeated by seven separate tests -- refitting was the single
# largest contributor to this file's runtime. `.dependence_profile_ci()`
# restores every mutation it makes to `fit$obj$env$last.par.best` on exit
# (see R/inference.R), so tests that only call it or
# `rpbnb_tmb_dependence_profile()` and read the result can safely share this
# one fit. The exception is the "last.par.best is restored" test below: it
# deliberately plants a sentinel value in `fit$obj$env$last.par.best` *before*
# calling `.dependence_profile_ci()` and asserts that the sentinel -- not
# `fit$optimizer$par` -- is what comes back afterwards. Sharing that test's
# fit would leave the sentinel sitting in the shared object's environment for
# every test that runs after it, silently changing what "before" means for
# any of them; that test keeps its own dedicated fit.
dp_shared_fit <- fit_rpbnb_tmb(
  y1 ~ x, y2 ~ x, dp_test_data(),
  dependence = copula("normal"), draws = 1L, keep = "full"
)

test_that("the working-scale profile brackets the working estimate", {
  fit <- dp_shared_fit

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
  fit <- dp_shared_fit

  # trace is a real, documented TMB::tmbprofile() argument, and this
  # function's own docs promise `...` is forwarded to it. Before the fix,
  # supplying it collided with the hardcoded `trace = FALSE` inside
  # do.call(), which threw "formal argument matched by multiple actual
  # arguments" -- caught by the inner try() and converted to NULL, so the
  # caller silently fell back to a Wald interval with no diagnostic.
  # trace = TRUE makes tmbprofile() print a "Profile value: ..." line per
  # step; capture.output() swallows that progress spew without touching what
  # is actually being asserted (that a caller-supplied trace still profiles).
  trace_output <- capture.output(
    ci_trace_true <- .dependence_profile_ci(fit, level = 0.95, trace = TRUE)
  )
  expect_false(is.null(ci_trace_true))
  expect_length(as.numeric(ci_trace_true), 2L)
  # Proves trace = TRUE actually reached tmbprofile() and produced output,
  # not just that the call didn't error.
  expect_true(any(grepl("Profile value", trace_output)))

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
  fit <- dp_shared_fit

  pr <- rpbnb_tmb_dependence_profile(fit)

  expect_identical(pr$parameter, c("rho", "tau"))
  expect_identical(unique(pr$method), "profile")
  expect_true(all(pr$lower < pr$estimate))
  expect_true(all(pr$estimate < pr$upper))
  expect_true(all(pr$level == 0.95))
})

test_that("a wider level nests the narrower one", {
  fit <- dp_shared_fit

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
  # If the fit fully pins, this filter can yield a zero-length vector and
  # all() over it is vacuously TRUE, which would let the two checks below
  # pass without ever comparing anything. Guard against that so the test can
  # actually fail.
  expect_true(length(ends) > 0)
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
  fit <- dp_shared_fit
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
  fit <- dp_shared_fit

  expect_error(rpbnb_tmb_dependence_profile(fit, level = 1), "between 0 and 1")
  expect_error(rpbnb_tmb_dependence_profile(fit, level = c(0.9, 0.95)),
               "between 0 and 1")
  expect_error(rpbnb_tmb_dependence_profile(list()), "rpbnb_tmb_fit")
})

test_that("lincomb and slice in ... are rejected rather than silently ignored", {
  # `.dependence_profile_ci()` checks `dots` for these before it ever touches
  # `fit$obj`, so no model fit is needed to exercise the guard: a caller
  # supplying `lincomb` would otherwise have tmbprofile() silently profile a
  # different linear combination of parameters than z_dep, and `slice = TRUE`
  # would silently swap a likelihood profile for a likelihood slice -- both
  # produce a plausible-looking wrong answer with `method` still saying
  # "profile", so they must error instead of being dropped like `obj`/`name`.
  expect_error(
    .dependence_profile_ci(list(), level = 0.95, lincomb = c(1, 0, 0)),
    "lincomb"
  )
  expect_error(
    .dependence_profile_ci(list(), level = 0.95, slice = TRUE),
    "slice"
  )
})

test_that("a liveness-probe failure warns with the real cause, not the keep = full advice", {
  # Fix 1 regression: before the fix, every NULL from .dependence_profile_ci()
  # -- whether fit$obj was NULL, the liveness probe failed, or tmbprofile()
  # itself threw -- was reported with the same "refit with keep = \"full\""
  # advice. That is only correct for the first case. Here the objective is
  # present but broken, which .dependence_profile_ci() must distinguish and
  # report as such, quoting the real underlying error.
  fit <- dp_shared_fit
  fit$obj <- list(
    fn = function(par) stop("simulated objective failure"),
    env = new.env()
  )

  expect_warning(
    pr <- rpbnb_tmb_dependence_profile(fit),
    "simulated objective failure"
  )
  expect_identical(unique(pr$method), "wald")
  warning_text <- tryCatch(
    { rpbnb_tmb_dependence_profile(fit); NA_character_ },
    warning = function(w) conditionMessage(w)
  )
  expect_false(grepl("keep = \"full\"", warning_text, fixed = TRUE))
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
