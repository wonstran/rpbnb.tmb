inference_test_data <- function(n = 120L) {
  set.seed(240724)
  x <- seq(-1, 1, length.out = n)
  data.frame(
    y1 = stats::rnbinom(n, size = 2.5, mu = exp(0.3 + 0.2 * x)),
    y2 = stats::rnbinom(n, size = 3.0, mu = exp(0.1 - 0.25 * x)),
    x = x
  )
}

fit_inference_mode <- function(inference = "full", keep = "postfit") {
  fit_rpbnb_tmb(
    y1 ~ x, y2 ~ x,
    data = inference_test_data(),
    dependence = "independence",
    inference = inference,
    keep = keep,
    control = rpbnb_tmb_control(iterlim = 150L, n_cores = 1L)
  )
}

test_that("full inference supplies compatible fixed and report summaries", {
  fit <- fit_inference_mode("full", "full")

  expect_s3_class(fit$sdreport, "rpbnb_sdreport")
  fixed <- summary(fit$sdreport, "fixed")
  report <- summary(fit$sdreport, "report")
  expect_identical(rownames(fixed), names(coef(fit)))
  expect_true(all(c("Estimate", "Std. Error") %in% colnames(fixed)))
  expect_true(all(c("m1", "m2") %in% rownames(report)))
  expect_equal(unname(report["m1", "Estimate"]), fit$m1)
  expect_equal(unname(report["m2", "Estimate"]), fit$m2)
  expect_true(is.matrix(vcov(fit)))
  expect_identical(dim(vcov(fit)), c(fit$npar, fit$npar))

  legacy <- TMB::sdreport(fit$obj)
  expect_equal(
    unname(fixed[, "Std. Error"]),
    unname(summary(legacy, "fixed")[, "Std. Error"]),
    tolerance = 1e-6
  )
  expect_equal(
    unname(report[, c("Estimate", "Std. Error")]),
    unname(summary(legacy, "report")[, c("Estimate", "Std. Error")]),
    tolerance = 1e-6
  )
})

test_that("diagonal and no-inference modes materialize covariance lazily", {
  diagonal <- fit_inference_mode("diag", "postfit")
  expect_null(diagonal$vcov)
  expect_length(diagonal$vcov_diag, diagonal$npar)
  diagonal_vcov <- vcov(diagonal)
  expect_true(all(is.finite(diag(diagonal_vcov))))
  expect_true(all(is.na(diagonal_vcov[row(diagonal_vcov) != col(diagonal_vcov)])))

  none <- fit_inference_mode("none", "postfit")
  expect_null(none$vcov)
  expect_true(all(is.na(none$se)))
  expect_true(all(is.na(vcov(none))))
  expect_true(is.na(none$sdreport$pdHess))
})

test_that("retention modes keep only requested post-fit state", {
  postfit <- fit_inference_mode("none", "postfit")
  expect_null(postfit$obj)
  expect_null(postfit$Y1)
  expect_null(postfit$Y2)
  expect_true(is.matrix(postfit$X1))
  expect_true(is.list(postfit$rp_meta))

  compact <- fit_inference_mode("none", "compact")
  expect_null(compact$obj)
  expect_null(compact$X1)
  expect_null(compact$X2)
  expect_null(compact$rp_meta)

  full <- fit_inference_mode("none", "full")
  expect_true(is.list(full$obj))
  expect_true(is.numeric(full$Y1))
  expect_true(is.numeric(full$Y2))
  expect_true(is.matrix(full$X1))
  expect_true(is.list(full$rp_meta))
})

test_that("inference and retention arguments are matched", {
  d <- inference_test_data(10L)
  expect_error(
    fit_rpbnb_tmb(y1 ~ x, y2 ~ x, d, inference = "partial"),
    "arg"
  )
  expect_error(
    fit_rpbnb_tmb(y1 ~ x, y2 ~ x, d, keep = "everything"),
    "arg"
  )
})

test_that("Frank reporting is smooth at independence and bounded at extremes", {
  at_independence <- .rpbnb_natural_report(
    c(log_m1 = 0, log_m2 = 0, z_dep = 0),
    family_code = 1L, lamLo = 0, lamHi = 0
  )
  at_extreme <- .rpbnb_natural_report(
    c(log_m1 = 0, log_m2 = 0, z_dep = -500),
    family_code = 1L, lamLo = 0, lamHi = 0
  )

  expect_equal(at_independence$theta$value, 0)
  expect_equal(at_independence$theta$derivative, 1)
  expect_equal(at_independence$tau$value, 0)
  expect_true(abs(at_extreme$theta$value) <= 35)
  expect_true(is.finite(at_extreme$tau$value))
})

test_that("Famoye saturation is marked as boundary inference", {
  report <- .rpbnb_natural_report(
    c(log_m1 = 0, log_m2 = 0, z_dep = 14),
    family_code = 0L, lamLo = -3, lamHi = 4
  )

  expect_true(report$lam$boundary)
})

# Every dependence link is bounded, so every family can saturate.  Flagging
# only Famoye would leave the same silent-boundary defect in the three copulas.
natural_report <- function(z, family_code, lamLo = 0, lamHi = 0) {
  .rpbnb_natural_report(
    c(log_m1 = 0, log_m2 = 0, z_dep = z),
    family_code = family_code, lamLo = lamLo, lamHi = lamHi
  )
}

test_that("Frank saturation against its bounded link is marked as boundary", {
  # 35 * tanh(z / 35) saturates; z = 87.6 reproduces an observed fit that
  # capped theta at 34.5 while reporting an interior-looking standard error.
  saturated <- natural_report(87.6, 1L)
  expect_true(saturated$theta$boundary)
  expect_true(saturated$tau$boundary)

  interior <- natural_report(5, 1L)
  expect_false(interior$theta$boundary)
  expect_false(interior$tau$boundary)
})

test_that("a strong but genuine Gaussian correlation keeps its standard error", {
  # |rho| < 1 is the model's real parameter domain, not an implementation cap.
  # rho = 0.964 with drho/dz = 0.07 is a well-identified estimate; suppressing
  # its standard error would be a false negative, not caution.
  strong <- natural_report(2, 2L)
  expect_gt(abs(strong$rho$value), 0.96)
  expect_gt(strong$rho$derivative, 0)
  expect_false(strong$rho$boundary)
  expect_false(strong$tau$boundary)
})

test_that("Gaussian is flagged only once tanh has actually saturated", {
  # tanh(z) reaches 1 in double precision near z = 19.1; past that the
  # delta-method derivative is exactly zero and the standard error is a lie.
  saturated <- natural_report(25, 2L)
  expect_equal(saturated$rho$derivative, 0)
  expect_true(saturated$rho$boundary)
  expect_true(saturated$tau$boundary)
  expect_true(is.finite(saturated$tau$derivative))
})

test_that("a hard clamp flags once reached, but not before", {
  # pmin/pmax clamps are the identity strictly inside (-20, 20): at z = 19.3
  # the value, objective and derivative are exactly what an unclamped
  # implementation would give, so flagging there would be a false positive.
  expect_false(natural_report(1, 3L)$theta$boundary)
  expect_false(natural_report(19.3, 3L)$theta$boundary)

  # Exact equality is a kink, not an interior point: the left derivative is
  # exp(20) and the right is 0, so the delta method has no unique derivative.
  # The template's CondExp picks the pass-through branch and clamp_deriv()
  # matches it, but that is a mechanical detail, not identified inference.
  expect_true(natural_report(20, 3L)$theta$boundary)
  expect_identical(natural_report(20, 3L)$theta$side, "upper")
  expect_true(natural_report(-20, 3L)$theta$boundary)
  expect_identical(natural_report(-20, 3L)$theta$side, "lower")

  expect_true(natural_report(20.5, 3L)$theta$boundary)
  expect_true(natural_report(25, 3L)$theta$boundary)
  expect_true(natural_report(25, 3L)$tau$boundary)
})

test_that("a clamped dispersion parameter is flagged with its side", {
  dispersion_report <- function(log_m1) {
    .rpbnb_natural_report(
      c(log_m1 = log_m1, log_m2 = 0, z_dep = 0),
      family_code = -1L, lamLo = 0, lamHi = 0
    )
  }
  expect_false(dispersion_report(19.3)$m1$boundary)
  expect_false(dispersion_report(0)$m1$boundary)

  low <- dispersion_report(-25)
  expect_true(low$m1$boundary)
  expect_identical(low$m1$side, "lower")
  expect_false(low$m2$boundary)

  high <- dispersion_report(25)
  expect_true(high$m1$boundary)
  expect_identical(high$m1$side, "upper")

  # The kink at exact equality is flagged too.
  expect_true(dispersion_report(-20)$m1$boundary)
  expect_identical(dispersion_report(-20)$m1$side, "lower")
})

test_that("boundary side metadata survives on the fitted object", {
  # A suppressed or scrolled-past warning must not be the only place the side
  # is recorded: lower and upper clamps call for opposite remedies.
  fit <- fit_rpbnb_tmb(
    y1 ~ x, y2 ~ x,
    data = inference_test_data(),
    dependence = copula("frank"),
    start = c(z_dep = 200),
    control = rpbnb_tmb_control(iterlim = 1L, n_cores = 1L)
  ) |> suppressWarnings()

  expect_true(length(fit$boundary_report) > 0)
  expect_named(fit$boundary_sides)
  expect_setequal(names(fit$boundary_sides), fit$boundary_report)
  expect_true(all(
    fit$boundary_sides %in% c("lower", "upper", "degenerate")
  ))
})

test_that("independence and interior dispersion reports are not flagged", {
  report <- natural_report(0, -1L)
  expect_false(report$m1$boundary)
  expect_false(report$m2$boundary)
})

test_that("the Famoye boundary threshold is pinned on both sides", {
  # The frozen admissible interval is an artefact of the starting values, so
  # proximity to it -- not just derivative collapse -- invalidates inference.
  expect_false(natural_report(3.8, 0L, -3, 4)$lam$boundary)
  expect_true(natural_report(4.0, 0L, -3, 4)$lam$boundary)
})

test_that("a clamped dispersion warns about dispersion, not dependence", {
  messages <- function(...) {
    collected <- character(0)
    withCallingHandlers(
      .warn_boundary_report(...),
      warning = function(cond) {
        collected <<- c(collected, conditionMessage(cond))
        invokeRestart("muffleWarning")
      }
    )
    paste(collected, collapse = " | ")
  }

  dispersion_only <- messages("m1", family_code = 1L,
                              sides = c(m1 = "lower"))
  expect_match(dispersion_only, "[Dd]ispersion")
  # Naming the Frank link here would be a false diagnosis: m1 has nothing to
  # do with the dependence family, and "change dependence family" is bad advice.
  expect_false(grepl("Frank", dispersion_only, fixed = TRUE))
  # Poisson advice belongs to the lower clamp, where m -> 0 ...
  expect_match(dispersion_only, "poisson_1")

  # ... and is the opposite model at the upper clamp, where m is enormous.
  upper <- messages("m1", family_code = 1L, sides = c(m1 = "upper"))
  expect_match(upper, "over-dispersed")
  expect_false(grepl("use poisson_1", upper, fixed = TRUE))

  dependence_only <- messages("theta", family_code = 1L)
  expect_match(dependence_only, "Frank")
  expect_false(grepl("[Dd]ispersion", dependence_only))

  both <- messages(c("m1", "theta"), family_code = 1L)
  expect_match(both, "[Dd]ispersion")
  expect_match(both, "Frank")

  expect_silent(.warn_boundary_report(character(0), family_code = 1L))
})

# R re-implements every natural-scale transform that the template ADREPORTs.
# Comparing estimates alone would miss a derivative mismatch, and it is the
# derivative that sets the reported standard errors.
cpp_report_value <- function(family_code, log_m1 = log(0.6), z_dep = 0.5) {
  x <- seq(-0.8, 0.8, length.out = 8L)
  data <- .build_tmb_data(
    Y1 = c(0, 1, 2, 0, 1, 3, 0, 2), Y2 = c(1, 0, 1, 2, 0, 1, 3, 0),
    X1 = cbind(1, x), X2 = cbind(1, x),
    rand_idx1 = integer(0), rand_idx2 = integer(0),
    Z1 = matrix(numeric(0), 1, 0), Z2 = matrix(numeric(0), 1, 0),
    dist1 = integer(0), dist2 = integer(0),
    sign1 = integer(0), sign2 = integer(0),
    family_code = family_code, pois1 = FALSE, pois2 = FALSE,
    lamLo = -1, lamHi = 1, est_method = 0L
  )
  parameters <- list(
    beta1 = c(0.1, 0.2), beta2 = c(-0.1, -0.15),
    log_sd1 = numeric(0), log_sd2 = numeric(0),
    log_m1 = log_m1, log_m2 = log(0.7), z_dep = z_dep,
    u1 = matrix(0, 8L, 0L), u2 = matrix(0, 8L, 0L)
  )
  obj <- suppressWarnings(
    .make_rpbnb_tmb_object(
      data = data, parameters = parameters, n_cores = 1L
    )$obj
  )
  invisible(obj$fn(obj$par))
  # Supplying an identity fixed-parameter Hessian makes the reported standard
  # error equal |d report / d par| exactly, so this extracts the template's own
  # AD Jacobian.  A finite difference cannot be used: the log_m clamp is a kink,
  # and a central difference there straddles both branches and returns half the
  # true one-sided derivative.
  identity_hessian <- diag(length(obj$par))
  report <- suppressWarnings(
    summary(TMB::sdreport(obj, hessian.fixed = identity_hessian), "report")
  )
  list(value = report[, "Estimate"], derivative = report[, "Std. Error"])
}

expect_transforms_agree <- function(family_code, z_dep = 0.5,
                                    log_m1 = log(0.6), label = "") {
  cpp <- cpp_report_value(family_code, log_m1 = log_m1, z_dep = z_dep)
  r <- .rpbnb_natural_report(
    c(log_m1 = log_m1, log_m2 = log(0.7), z_dep = z_dep),
    family_code = family_code, lamLo = -1, lamHi = 1
  )
  shared <- intersect(names(cpp$value), names(r))
  expect_gt(length(shared), 0)
  for (item in shared) {
    # Relative comparison: these quantities span 1e-09 (a clamped dispersion)
    # to 1e+08 (a saturating Clayton theta), so an absolute tolerance would
    # accept literally any answer at the small end -- including the zero the
    # old clamp_deriv() returned.
    scale <- max(abs(r[[item]]$derivative), abs(cpp$derivative[[item]]), 1e-300)
    expect_equal(
      unname(cpp$value[[item]]), r[[item]]$value,
      tolerance = 1e-8,
      info = sprintf("%s %s: value", label, item)
    )
    expect_lt(
      abs(abs(r[[item]]$derivative) - cpp$derivative[[item]]) / scale,
      1e-6
    )
  }
}

test_that("R and C++ transforms agree in value and derivative", {
  # Interior, and close enough to each family's cap that the transforms are in
  # their nonlinear regime.
  for (family_code in c(0L, 1L, 2L, 3L)) {
    expect_transforms_agree(family_code, z_dep = 0.5,
                            label = sprintf("family %d interior", family_code))
  }
  expect_transforms_agree(0L, z_dep = 8, label = "Famoye near cap")
  expect_transforms_agree(1L, z_dep = 60, label = "Frank near cap")
  expect_transforms_agree(2L, z_dep = 6, label = "Gaussian near saturation")
  expect_transforms_agree(3L, z_dep = 19, label = "Clayton near clamp")
})

test_that("R and C++ agree on the dispersion clamp, including at equality", {
  for (log_m1 in c(-1, -19, -20, log(0.6))) {
    expect_transforms_agree(
      -1L, log_m1 = log_m1, label = sprintf("log_m1 = %g", log_m1)
    )
  }
})

test_that("the dispersion clamp contract fails if clamp_deriv regresses", {
  # Proving the regression test can actually fail.  The defect it guards
  # against returned exactly zero at |log_m| == 20, where the template passes
  # the input through and the true derivative is exp(-20).
  cpp <- cpp_report_value(-1L, log_m1 = -20)
  expect_equal(unname(cpp$derivative[["m1"]]), exp(-20), tolerance = 1e-12)

  regressed <- 0
  scale <- max(abs(regressed), cpp$derivative[["m1"]], 1e-300)
  expect_gt(abs(abs(regressed) - cpp$derivative[["m1"]]) / scale, 1e-6)
})

test_that("report covariance stays consistent with report standard errors", {
  fit <- fit_rpbnb_tmb(
    y1 ~ x, y2 ~ x,
    data = inference_test_data(),
    dependence = copula("frank"),
    start = c(z_dep = 200),
    inference = "full",
    control = rpbnb_tmb_control(iterlim = 1L, n_cores = 1L)
  ) |> suppressWarnings()

  flagged <- fit$boundary_report
  expect_true(length(flagged) > 0)

  # Withdrawing a standard error has to withdraw the matching covariances too,
  # otherwise a caller can reconstruct the quantity we declared unidentified.
  expect_true(all(is.na(fit$sdreport$sd[flagged])))
  for (item in flagged) {
    expect_true(all(is.na(fit$sdreport$cov[item, ])))
    expect_true(all(is.na(fit$sdreport$cov[, item])))
  }

  # Where inference did survive, the diagonal must still be the squared
  # standard error.  (This fit uses iterlim = 1, so the Hessian need not be
  # positive definite and unflagged entries may legitimately be NaN.)
  unflagged <- setdiff(names(fit$sdreport$sd), flagged)
  usable <- unflagged[is.finite(fit$sdreport$sd[unflagged])]
  expect_equal(
    unname(diag(fit$sdreport$cov)[usable]),
    unname(fit$sdreport$sd[usable]^2),
    tolerance = 1e-8
  )
})

test_that("boundary dependence parameters report NA standard errors", {
  fit <- fit_rpbnb_tmb(
    y1 ~ x, y2 ~ x,
    data = inference_test_data(),
    dependence = copula("frank"),
    start = c(z_dep = 200),
    inference = "full",
    control = rpbnb_tmb_control(iterlim = 1L, n_cores = 1L)
  ) |> suppressWarnings()

  expect_true("theta" %in% fit$boundary_report)
  report <- summary(fit$sdreport, "report")
  expect_true(is.na(report["theta", "Std. Error"]))
})

test_that("a saturated copula fit warns instead of failing silently", {
  expect_warning(
    fit_rpbnb_tmb(
      y1 ~ x, y2 ~ x,
      data = inference_test_data(),
      dependence = copula("frank"),
      start = c(z_dep = 200),
      control = rpbnb_tmb_control(iterlim = 1L, n_cores = 1L)
    ),
    "boundary"
  )
})

test_that("draw means stream with numerical equivalence", {
  xb <- c(-0.4, 0.2, 0.7)
  xr <- matrix(c(-1, 0.5, 2), ncol = 1L)
  dev <- matrix(c(-0.3, 0.1, 0.4, 0.8), ncol = 1L)
  weights <- c(0.7, 1.1, 1.4, 1.8)

  dense <- vapply(seq_len(nrow(dev)), function(r) {
    exp(xb + as.vector(xr %*% dev[r, ]))
  }, numeric(length(xb)))
  expect_equal(.draw_mean_exp(xb, xr, dev), rowMeans(dense))
  expect_equal(
    .draw_mean_weighted_exp(xb, xr, dev, weights),
    rowMeans(sweep(dense, 2L, weights, `*`))
  )
})

test_that("compact fits reject marginal effects before accessing fit state", {
  compact <- fit_inference_mode("none", "compact")
  expect_error(
    rpbnb_tmb_marginal_effects(compact, which = "y1"),
    'keep = "postfit"'
  )
})
