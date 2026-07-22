parallel_tmb_fixture <- function() {
  empty_draws <- matrix(numeric(0), nrow = 1L, ncol = 0L)
  data <- .build_tmb_data(
    Y1 = c(0, 1), Y2 = c(1, 0),
    X1 = matrix(1, nrow = 2L), X2 = matrix(1, nrow = 2L),
    rand_idx1 = integer(0), rand_idx2 = integer(0),
    Z1 = empty_draws, Z2 = empty_draws,
    dist1 = integer(0), dist2 = integer(0),
    sign1 = integer(0), sign2 = integer(0),
    family_code = -1L, pois1 = FALSE, pois2 = FALSE,
    lamLo = -1, lamHi = 1, n_cores = 1L
  )
  parameters <- list(
    beta1 = 0, beta2 = 0,
    log_sd1 = numeric(0), log_sd2 = numeric(0),
    log_m1 = log(0.5), log_m2 = log(0.5), z_dep = 0
  )
  list(data = data, parameters = parameters,
       map = list(z_dep = factor(NA)))
}

test_that("rpbnb_tmb_control validates n_cores", {
  expect_error(rpbnb_tmb_control(n_cores = 0), "n_cores")
  expect_error(rpbnb_tmb_control(n_cores = -1), "n_cores")
  expect_error(rpbnb_tmb_control(n_cores = NA_integer_), "n_cores")
  expect_error(rpbnb_tmb_control(n_cores = 1.5), "n_cores")
  expect_error(rpbnb_tmb_control(n_cores = "2"), "n_cores")
  expect_error(rpbnb_tmb_control(n_cores = c(1, 2)), "n_cores")
  expect_equal(rpbnb_tmb_control(n_cores = 2)$n_cores, 2L)
})

test_that("TMB thread configuration realizes a valid request", {
  previous <- TMB::openmp(DLL = "rpbnb.tmb")
  on.exit(TMB::openmp(n = previous, DLL = "rpbnb.tmb"), add = TRUE)

  expect_identical(
    expect_no_error(.configure_tmb_threads(1L, DLL = "rpbnb.tmb")),
    1L
  )
})

test_that("TMB object construction records the realized thread count", {
  fixture <- parallel_tmb_fixture()
  previous <- TMB::openmp(DLL = "rpbnb.tmb")
  on.exit(TMB::openmp(n = previous, DLL = "rpbnb.tmb"), add = TRUE)

  expect_no_error({
    configured <- .make_rpbnb_tmb_object(
      data = fixture$data,
      parameters = fixture$parameters,
      map = fixture$map,
      n_cores = 1L
    )
    expect_named(configured, c("obj", "n_cores"))
    expect_identical(configured$n_cores, 1L)
    expect_true(is.function(configured$obj$fn))
    expect_true(is.finite(configured$obj$fn(configured$obj$par)))
  })
})

test_that("fitted models record requested and realized threads", {
  set.seed(240722)
  x <- seq(-1, 1, length.out = 100)
  d <- data.frame(
    y1 = stats::rnbinom(100, size = 2, mu = exp(0.2 + 0.3 * x)),
    y2 = stats::rnbinom(100, size = 3, mu = exp(-0.1 - 0.2 * x)),
    x = x
  )

  fit <- fit_rpbnb_tmb(
    y1 ~ x, y2 ~ x, data = d,
    dependence = "independence",
    control = rpbnb_tmb_control(iterlim = 100L, n_cores = 1L)
  )

  expect_identical(fit$parallel, list(requested = 1L, realized = 1L))
})
