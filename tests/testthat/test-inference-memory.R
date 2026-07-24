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

test_that("default inference avoids TMB sdreport retaping", {
  fit_source <- paste(
    readLines(
      testthat::test_path("..", "..", "R", "fit_rpbnb_tmb.R"),
      warn = FALSE
    ),
    collapse = "\n"
  )
  expect_false(grepl("TMB::sdreport(", fit_source, fixed = TRUE))
})

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

test_that("marginal effects do not construct observation-by-draw matrices", {
  source <- paste(
    readLines(
      testthat::test_path("..", "..", "R", "marginal_effects.R"),
      warn = FALSE
    ),
    collapse = "\n"
  )
  expect_false(grepl("mu_mat <- matrix", source, fixed = TRUE))
  expect_false(grepl("mu_mat_t <- matrix", source, fixed = TRUE))
  expect_false(grepl("vapply(seq_len(R)", source, fixed = TRUE))
  expect_false(grepl("vapply(seq_len(R_t)", source, fixed = TRUE))
})

test_that("compact fits reject marginal effects before accessing fit state", {
  compact <- fit_inference_mode("none", "compact")
  expect_error(
    rpbnb_tmb_marginal_effects(compact, which = "y1"),
    'keep = "postfit"'
  )
})
