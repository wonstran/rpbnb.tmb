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

test_that("Gaussian and Clayton saturation are marked as boundary", {
  expect_true(natural_report(4, 2L)$rho$boundary)
  expect_true(natural_report(4, 2L)$tau$boundary)
  expect_false(natural_report(0.5, 2L)$rho$boundary)

  expect_true(natural_report(25, 3L)$theta$boundary)
  expect_true(natural_report(25, 3L)$tau$boundary)
  expect_false(natural_report(1, 3L)$theta$boundary)
})

test_that("independence and dispersion reports are never spuriously flagged", {
  report <- natural_report(0, -1L)
  expect_false(report$m1$boundary)
  expect_false(report$m2$boundary)
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
