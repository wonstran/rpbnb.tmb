# The dense truck report was published with continuous elasticities of 1e-16
# through 1e-18 -- printed as 0.0000, which reads as "no effect" rather than
# "this number means nothing".  The cause was centring: an elasticity carries a
# leading x-bar factor, and a centred regressor has x-bar = 0.  `scaling`
# restates both diagnostics in the covariates' original units.

scale_data <- function(n = 300, seed = 11) {
  set.seed(seed)
  # x1 is deliberately far from zero and widely scaled -- the SR40_MI3 shape
  # that made this matter in the first place.
  x1 <- rnorm(n, mean = 45, sd = 6)
  x2 <- rbinom(n, 1, 0.4)
  mu1 <- exp(-1.2 + 0.03 * x1 + 0.4 * x2)
  mu2 <- exp(-0.9 + 0.02 * x1 - 0.3 * x2)
  data.frame(y1 = rpois(n, mu1), y2 = rpois(n, mu2), x1 = x1, x2 = x2)
}

test_that(".scaling_vec defaults to identity and validates its input", {
  cn <- c("(Intercept)", "x1", "x2")
  sv <- rpbnb.tmb:::.scaling_vec(NULL, cn)
  expect_identical(unname(sv$center), c(0, 0, 0))
  expect_identical(unname(sv$scale), c(1, 1, 1))

  sv <- rpbnb.tmb:::.scaling_vec(list(x1 = c(center = 45, scale = 6)), cn)
  expect_identical(unname(sv$center), c(0, 45, 0))
  expect_identical(unname(sv$scale), c(1, 6, 1))

  # A name the design does not carry is ignored, not an error: callers pass one
  # scaling list for two equations that need not share covariates.
  expect_silent(rpbnb.tmb:::.scaling_vec(list(zz = c(center = 1, scale = 2)), cn))

  expect_error(rpbnb.tmb:::.scaling_vec(list(x1 = c(center = 0)), cn), "center")
  expect_error(
    rpbnb.tmb:::.scaling_vec(list(x1 = c(center = 0, scale = 0)), cn),
    "positive"
  )
  expect_error(rpbnb.tmb:::.scaling_vec(c(x1 = 1), cn), "named list")
})

test_that("scaling restates a standardized fit in the covariates' own units", {
  skip_on_cran()
  d <- scale_data()
  ctr <- mean(d$x1); scl <- stats::sd(d$x1)
  ds <- d; ds$x1 <- (d$x1 - ctr) / scl

  ctrl <- rpbnb_tmb_control(n_cores = 1)
  # No random coefficients, so the raw and standardized designs really are the
  # same model reparameterized -- which is what makes this an invariance test
  # rather than a comparison of two different fits.
  common <- list(formula_1 = y1 ~ x1 + x2, formula_2 = y2 ~ x1 + x2,
                 dependence = copula("normal"), method = "sml",
                 draws = 20, seed = 3, control = ctrl)
  fit_raw <- suppressWarnings(do.call(fit_rpbnb_tmb, c(common, list(data = d))))
  fit_std <- suppressWarnings(do.call(fit_rpbnb_tmb, c(common, list(data = ds))))

  sc <- list(x1 = c(center = ctr, scale = scl))
  quiet <- function(e) invisible(capture.output(r <- e, type = "output"))

  quiet(me_raw <- rpbnb_tmb_marginal_effects(fit_raw, which = "y1"))
  quiet(me_std <- rpbnb_tmb_marginal_effects(fit_std, which = "y1"))
  quiet(me_bak <- rpbnb_tmb_marginal_effects(fit_std, which = "y1", scaling = sc))
  quiet(el_raw <- rpbnb_tmb_elasticities(fit_raw, which = "y1"))
  quiet(el_bak <- rpbnb_tmb_elasticities(fit_std, which = "y1", scaling = sc))

  i <- match("x1", me_std$Name)
  # The chain rule, exactly: a standardized AME is the raw one times the scale.
  expect_equal(me_bak$Estimate[i], me_std$Estimate[i] / scl)
  expect_equal(me_bak$`Std. Error`[i], me_std$`Std. Error`[i] / scl)
  # And it agrees with the fit that never standardized at all.
  expect_equal(me_bak$Estimate[i], me_raw$Estimate[i], tolerance = 1e-5)
  expect_equal(el_bak$Estimate[i], el_raw$Estimate[i], tolerance = 1e-5)

  # Centring is what destroys the elasticity, so it must be what `scaling`
  # restores: without it the number is ~0, with it it is the real elasticity.
  quiet(el_std <- rpbnb_tmb_elasticities(fit_std, which = "y1"))
  expect_lt(abs(el_std$Estimate[i]), 1e-8)
  expect_gt(abs(el_bak$Estimate[i]), 0.1)

  # x2 is binary and was never transformed, so nothing about it may move.
  j <- match("x2", me_std$Name)
  expect_identical(me_bak$Estimate[j], me_std$Estimate[j])
  expect_identical(me_bak$`Std. Error`[j], me_std$`Std. Error`[j])
  expect_identical(el_bak$Estimate[j], el_std$Estimate[j])

  # Omitting `scaling` must leave the existing behaviour untouched.
  quiet(me_none <- rpbnb_tmb_marginal_effects(fit_std, which = "y1"))
  expect_identical(me_none$Estimate, me_std$Estimate)
})

test_that("raw_scale adds a section and is absent by default", {
  results_dir <- file.path(tempdir(), "truck-results-rawscale")
  unlink(results_dir, recursive = TRUE)
  args <- list(model_summary = "Summary: fitted model",
               marginal_effects = "x  1.25",
               elasticities = "x  0.50",
               results_dir = results_dir,
               timestamp = as.POSIXct("2026-08-03 10:53:08", tz = "UTC"))

  plain <- readLines(do.call(rpbnb.tmb:::.write_truck_results_markdown, args),
                     warn = FALSE)
  expect_false(any(grepl("original units", plain)))

  withr <- do.call(rpbnb.tmb:::.write_truck_results_markdown,
                   c(args, list(raw_scale = "x  2.50",
                                scaling = list(x = c(center = 45, scale = 6)))))
  report <- readLines(withr, warn = FALSE)
  expect_true(any(grepl("## Marginal effects and elasticities in original units",
                        report, fixed = TRUE)))
  expect_true(any(grepl("center =      45.0000   scale =       6.0000", report)))
  expect_true(any(grepl("x  2.50", report, fixed = TRUE)))
  # Additive: the three original sections survive unchanged.
  expect_true(all(c("## Model fit summary",
                    "## Average marginal effects (AME)",
                    "## Elasticities / semi-elasticities (AME)") %in% report))
})

# `LNAADT_3` is log(AADT).  Fed to the elasticity formula as an ordinary
# regressor it returns x-bar * b = 8.59 -- the elasticity with respect to the
# LOG of traffic.  The elasticity with respect to traffic is b itself, 0.871.
# Nothing in the output distinguishes the two, which is what makes it dangerous.

test_that(".log_vars_flag validates against the design", {
  cn <- c("(Intercept)", "lv", "x2"); sel <- 2:3; is_bin <- c(FALSE, TRUE)
  expect_identical(rpbnb.tmb:::.log_vars_flag(NULL, cn, sel, is_bin),
                   c(FALSE, FALSE))
  expect_identical(rpbnb.tmb:::.log_vars_flag("lv", cn, sel, is_bin),
                   c(TRUE, FALSE))
  expect_error(rpbnb.tmb:::.log_vars_flag("nope", cn, sel, is_bin), "not found")
  expect_error(rpbnb.tmb:::.log_vars_flag("x2", cn, sel, is_bin), "binary")
  expect_error(rpbnb.tmb:::.log_vars_flag(1L, cn, sel, is_bin), "character")
})

test_that("log_vars reports per unit of the underlying variable", {
  skip_on_cran()
  set.seed(19)
  n <- 300
  v  <- exp(rnorm(n, mean = 9.9, sd = 0.75))     # the raw variable (e.g. AADT)
  lv <- log(v)
  x2 <- rbinom(n, 1, 0.4)
  mu1 <- exp(-6 + 0.8 * lv + 0.3 * x2)
  mu2 <- exp(-5 + 0.5 * lv - 0.2 * x2)
  d <- data.frame(y1 = rpois(n, mu1), y2 = rpois(n, mu2), lv = lv, x2 = x2)

  fit <- suppressWarnings(fit_rpbnb_tmb(
    y1 ~ lv + x2, y2 ~ lv + x2, data = d, dependence = copula("normal"),
    method = "sml", draws = 20, seed = 4,
    control = rpbnb_tmb_control(n_cores = 1)))
  quiet <- function(e) invisible(capture.output(r <- e, type = "output"))

  quiet(el_plain <- rpbnb_tmb_elasticities(fit, which = "y1"))
  quiet(el_log   <- rpbnb_tmb_elasticities(fit, which = "y1", log_vars = "lv"))
  i <- match("lv", el_plain$Name)
  b <- unname(fit$coef["b1:lv"])

  # No random coefficient on lv, so d(log mu)/d(log v) is exactly the slope.
  expect_equal(el_log$Estimate[i], b, tolerance = 1e-6)
  # The untreated version is x-bar * b -- much larger, and the whole problem.
  expect_equal(el_plain$Estimate[i], mean(lv) * b, tolerance = 1e-6)
  expect_gt(abs(el_plain$Estimate[i]), 5 * abs(el_log$Estimate[i]))
  expect_identical(el_log$Type[i], "log-continuous")

  # AME per unit of v: dmu/dv = (dmu/dlv)/v.
  quiet(me_log <- rpbnb_tmb_marginal_effects(fit, which = "y1", log_vars = "lv"))
  quiet(me_pl  <- rpbnb_tmb_marginal_effects(fit, which = "y1"))
  mu <- as.vector(predict(fit)[, 1])
  expect_equal(me_log$Estimate[i], mean(b * mu / v), tolerance = 1e-5)
  expect_identical(me_log$Type[i], "log-continuous")

  # The binary column must be untouched by any of it.
  j <- match("x2", me_pl$Name)
  expect_identical(me_log$Estimate[j], me_pl$Estimate[j])
  expect_identical(me_log$`Std. Error`[j], me_pl$`Std. Error`[j])

  # log_vars composes with scaling: standardizing the logged column must not
  # change the per-unit-of-v answer, since s cancels out of both formulas.
  ctr <- mean(lv); scl <- stats::sd(lv)
  ds <- d; ds$lv <- (lv - ctr) / scl
  fit_s <- suppressWarnings(fit_rpbnb_tmb(
    y1 ~ lv + x2, y2 ~ lv + x2, data = ds, dependence = copula("normal"),
    method = "sml", draws = 20, seed = 4,
    control = rpbnb_tmb_control(n_cores = 1)))
  quiet(el_both <- rpbnb_tmb_elasticities(
    fit_s, which = "y1", log_vars = "lv",
    scaling = list(lv = c(center = ctr, scale = scl))))
  expect_equal(el_both$Estimate[i], el_log$Estimate[i], tolerance = 1e-4)
})
