api_test_data <- function(n = 80L) {
  set.seed(240724)
  x <- seq(-1, 1, length.out = n)
  data.frame(
    y1 = stats::rpois(n, exp(0.2 + 0.25 * x)),
    y2 = stats::rpois(n, exp(0.1 - 0.15 * x)),
    x = x
  )
}

test_that("mapped Poisson margins retain coefficients and count only free df", {
  fit <- fit_rpbnb_tmb(
    y1 ~ x, y2 ~ x, api_test_data(),
    dependence = "independence",
    poisson_1 = TRUE, poisson_2 = TRUE,
    inference = "none"
  )

  expect_true(all(is.finite(coef(fit))))
  expect_equal(unname(coef(fit)[c("log_m1", "log_m2")]), rep(log(POISSON_M), 2))
  expect_equal(fit$npar, 4L)
  expect_equal(attr(logLik(fit), "df"), 4L)
})

test_that("base information criteria honor their generic arguments", {
  fit <- fit_rpbnb_tmb(
    y1 ~ x, y2 ~ x, api_test_data(),
    dependence = "independence", inference = "none"
  )

  expect_equal(AIC(fit, k = 0), -2 * fit$logLik)
  expect_equal(
    BIC(fit),
    -2 * fit$logLik + log(fit$nobs) * fit$npar
  )
})

test_that("fitting preserves an existing caller RNG stream", {
  d <- api_test_data(20L)
  set.seed(9182)
  seed_before <- .Random.seed

  fit_rpbnb_tmb(
    y1 ~ x, y2 ~ x, d,
    random_1 = "x", draws = 5L,
    dependence = "independence", inference = "none"
  )

  expect_identical(.Random.seed, seed_before)
})

test_that("fitting does not create a caller RNG stream when none existed", {
  d <- api_test_data(20L)
  old_seed_exists <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (old_seed_exists) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (old_seed_exists) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    rm(".Random.seed", envir = .GlobalEnv)
  }

  fit_rpbnb_tmb(
    y1 ~ x, y2 ~ x, d,
    random_1 = "x", draws = 5L,
    dependence = "independence", inference = "none"
  )

  expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
})

test_that("predict returns integrated responses and link-scale predictions", {
  data <- api_test_data(80L)
  fit <- fit_rpbnb_tmb(
    y1 ~ x, y2 ~ x, data = data,
    dependence = "independence",
    inference = "none",
    control = rpbnb_tmb_control(iterlim = 100L)
  )

  response <- predict(fit)
  link <- predict(fit, type = "link")
  expect_identical(dim(response), c(fit$nobs, 2L))
  expect_identical(colnames(response), c("y1", "y2"))
  expect_equal(exp(link), response)

  newdata <- data.frame(x = c(-0.5, 0.5))
  new_response <- predict(fit, newdata = newdata)
  expect_identical(dim(new_response), c(2L, 2L))
  expect_true(all(is.finite(new_response)))
})

test_that("compact fits can predict new data and reject unavailable rows", {
  data <- api_test_data(60L)
  fit <- fit_rpbnb_tmb(
    y1 ~ x, y2 ~ x, data = data,
    dependence = "independence",
    inference = "none",
    keep = "compact",
    control = rpbnb_tmb_control(iterlim = 100L)
  )

  expect_error(predict(fit), "newdata")
  expect_identical(
    dim(predict(fit, newdata = data.frame(x = 0))),
    c(1L, 2L)
  )
})

test_that("diagonal inference warns before marginal-effect standard errors", {
  data <- api_test_data(80L)
  fit <- fit_rpbnb_tmb(
    y1 ~ x, y2 ~ x, data = data,
    dependence = "independence",
    inference = "diag",
    control = rpbnb_tmb_control(iterlim = 100L)
  )

  expect_warning(
    visible <- withVisible(
      rpbnb_tmb_marginal_effects(fit, which = "y1")
    ),
    "full covariance"
  )
  expect_false(visible$visible)
  expect_true(all(is.na(visible$value$`Std. Error`)))
})

test_that("dependence specifications resolve to one shared family code", {
  expect_identical(.resolve_family_code("independence"), -1L)
  expect_identical(.resolve_family_code("famoye"), 0L)
  expect_identical(.resolve_family_code(copula("frank")), 1L)
  expect_identical(.resolve_family_code(copula("normal")), 2L)
  expect_identical(.resolve_family_code(copula("kimeldorf")), 3L)
  # "gaussian" is not a valid dependence string; only copula("normal") is.
  expect_error(.resolve_family_code("gaussian"), "must be")
})

test_that("a transformed formula cannot misalign the two equations", {
  # complete.cases() on the raw variables passes here -- the NaN is created by
  # log() inside formula_1 only. Before both frames shared one row mask, eq 1
  # silently used original rows 2-6 while eq 2 used rows 1-6, and the template
  # (which takes n from Y1 and indexes Y2/X2 with it) paired observation i of
  # one outcome with observation i+1 of the other. It fitted without complaint.
  d <- data.frame(
    y1 = c(0, 1, 1, 2, 2, 3),
    y2 = c(0, 0, 1, 1, 2, 2),
    x  = c(-1, 1, 2, 3, 4, 5)
  )

  prep <- suppressWarnings(.prepare_bnb_data(y1 ~ log(x), y2 ~ x, d))

  expect_identical(length(prep$Y1), length(prep$Y2))
  expect_identical(nrow(prep$X1), nrow(prep$X2))
  expect_identical(prep$n, nrow(prep$X1))
  # Row 1 (x = -1) is the only one dropped, from BOTH equations.
  expect_identical(prep$Y1, as.integer(d$y1[-1L]))
  expect_identical(prep$Y2, as.integer(d$y2[-1L]))
  expect_identical(prep$n, 5L)
})

test_that("a transformation dropping rows from either side drops both", {
  d <- data.frame(
    y1 = c(0, 1, 1, 2, 2, 3),
    y2 = c(0, 0, 1, 1, 2, 2),
    x  = c(-1, 1, 2, 3, 4, 5),
    w  = c(1, 2, 3, 4, 5, -1)
  )

  # One row lost to each equation, at opposite ends: four survive in both.
  prep <- suppressWarnings(.prepare_bnb_data(y1 ~ log(x), y2 ~ log(w), d))

  expect_identical(prep$n, 4L)
  expect_identical(prep$Y1, as.integer(d$y1[2:5]))
  expect_identical(prep$Y2, as.integer(d$y2[2:5]))
  expect_identical(nrow(prep$data), 4L)
})

test_that("copula(par = ) is rejected by the fitter rather than ignored", {
  # It is validated by copula() and then read only by simulate_rpbnb_tmb().
  # Accepting it here would return a fit bit-identical to one that never
  # passed it, which is the failure mode this rejection exists to prevent.
  d <- api_test_data(40L)

  expect_error(
    fit_rpbnb_tmb(y1 ~ x, y2 ~ x, d, dependence = copula("normal", par = 0.75),
                  inference = "none"),
    "simulation argument"
  )
  expect_no_error(
    fit_rpbnb_tmb(y1 ~ x, y2 ~ x, d, dependence = copula("normal"),
                  inference = "none",
                  control = rpbnb_tmb_control(iterlim = 5L))
  )
})
