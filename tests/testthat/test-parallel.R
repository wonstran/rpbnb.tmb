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
    lamLo = -1, lamHi = 1
  )
  parameters <- list(
    beta1 = 0, beta2 = 0,
    log_sd1 = numeric(0), log_sd2 = numeric(0),
    log_m1 = log(0.5), log_m2 = log(0.5), z_dep = 0
  )
  list(data = data, parameters = parameters,
       map = list(z_dep = factor(NA)))
}

copula_parallel_fixture <- function(family_code) {
  x <- seq(-0.8, 0.8, length.out = 8L)
  z1 <- matrix(c(0.15, 0.35, 0.65, 0.85), ncol = 1L)
  z2 <- matrix(c(0.75, 0.25, 0.55, 0.45), ncol = 1L)
  data <- .build_tmb_data(
    Y1 = c(0, 1, 2, 0, 1, 3, 0, 2),
    Y2 = c(1, 0, 1, 2, 0, 1, 3, 0),
    X1 = cbind("(Intercept)" = 1, x = x),
    X2 = cbind("(Intercept)" = 1, x = x),
    rand_idx1 = 2L, rand_idx2 = 2L,
    Z1 = z1, Z2 = z2,
    dist1 = 0L, dist2 = 0L,
    sign1 = 1L, sign2 = 1L,
    family_code = family_code, pois1 = FALSE, pois2 = FALSE,
    lamLo = 0, lamHi = 0
  )
  parameters <- list(
    beta1 = c(0.1, 0.2), beta2 = c(-0.1, -0.15),
    log_sd1 = log(0.2), log_sd2 = log(0.25),
    log_m1 = log(0.6), log_m2 = log(0.7),
    z_dep = if (family_code == 1L) 1.5 else atanh(0.25)
  )
  list(data = data, parameters = parameters, map = NULL)
}

parallel_fit_data <- function() {
  set.seed(240722)
  x <- seq(-1, 1, length.out = 100)
  data.frame(
    y1 = stats::rnbinom(100, size = 2, mu = exp(0.2 + 0.3 * x)),
    y2 = stats::rnbinom(100, size = 3, mu = exp(-0.1 - 0.2 * x)),
    x = x
  )
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
  d <- parallel_fit_data()

  fit <- fit_rpbnb_tmb(
    y1 ~ x, y2 ~ x, data = d,
    dependence = "independence",
    control = rpbnb_tmb_control(iterlim = 100L, n_cores = 1L)
  )

  expect_identical(fit$parallel, list(requested = 1L, realized = 1L))
})

test_that("serial and parallel fitted models are numerically equivalent", {
  previous <- TMB::openmp(DLL = "rpbnb.tmb")
  on.exit(TMB::openmp(n = previous, DLL = "rpbnb.tmb"), add = TRUE)
  supported <- as.integer(TMB::openmp(max = TRUE, DLL = "rpbnb.tmb")[[1L]])
  if (supported < 2L) skip("TMB runtime supports only one thread")

  d <- parallel_fit_data()
  fit_serial <- fit_rpbnb_tmb(
    y1 ~ x, y2 ~ x, data = d,
    dependence = "independence",
    control = rpbnb_tmb_control(iterlim = 100L, n_cores = 1L)
  )
  fit_parallel <- fit_rpbnb_tmb(
    y1 ~ x, y2 ~ x, data = d,
    dependence = "independence",
    control = rpbnb_tmb_control(iterlim = 100L, n_cores = 2L)
  )

  expect_equal(coef(fit_parallel), coef(fit_serial), tolerance = 1e-6)
  expect_equal(fit_parallel$se, fit_serial$se, tolerance = 1e-5)
  expect_equal(fit_parallel$optimizer$objective,
               fit_serial$optimizer$objective, tolerance = 1e-7)
  expect_true(all(is.finite(coef(fit_parallel))))
  expect_true(all(is.finite(fit_parallel$se)))
  expect_identical(fit_parallel$parallel,
                   list(requested = 2L, realized = 2L))
})

test_that("print reports the realized TMB thread count", {
  x <- structure(
    list(
      logLik = -10, nobs = 5L, npar = 1L,
      dependence = "independence", coef = c(beta = 0),
      parallel = list(requested = 4L, realized = 2L)
    ),
    class = "rpbnb_tmb_fit"
  )
  output <- capture.output(print(x))
  expect_true(any(output == "  TMB threads: 2"))

  x$parallel <- NULL
  legacy_output <- expect_no_error(capture.output(print(x)))
  expect_false(any(grepl("TMB threads", legacy_output, fixed = TRUE)))
})

test_that("model DLL reports its compile-time OpenMP capability", {
  makeconf <- file.path(R.home("etc"), .Platform$r_arch, "Makeconf")
  if (!file.exists(makeconf)) makeconf <- file.path(R.home("etc"), "Makeconf")
  openmp_line <- grep(
    "^SHLIB_OPENMP_CXXFLAGS[[:space:]]*=",
    readLines(makeconf, warn = FALSE), value = TRUE
  )
  toolchain_openmp <- length(openmp_line) == 1L &&
    nzchar(trimws(sub("^[^=]*=", "", openmp_line)))

  fixture <- parallel_tmb_fixture()
  previous <- TMB::openmp(DLL = "rpbnb.tmb")
  on.exit(TMB::openmp(n = previous, DLL = "rpbnb.tmb"), add = TRUE)
  supported <- as.integer(TMB::openmp(max = TRUE, DLL = "rpbnb.tmb")[[1L]])
  requested <- if (toolchain_openmp && supported >= 2L) 2L else 1L
  configured <- .make_rpbnb_tmb_object(
    data = fixture$data,
    parameters = fixture$parameters,
    map = fixture$map,
    n_cores = requested
  )

  expect_identical(configured$n_cores, requested)
  expect_identical(
    configured$obj$report()$openmp_compiled,
    as.integer(toolchain_openmp)
  )
})

test_that("C++ likelihood declares a TMB parallel accumulator", {
  cpp <- readLines(
    testthat::test_path("..", "..", "src", "rpbnb.tmb.cpp"),
    warn = FALSE
  )
  expect_true(any(grepl("parallel_accumulator<Type> nll", cpp, fixed = TRUE)))
})

test_that("serial and parallel copula objectives and gradients agree", {
  previous <- TMB::openmp(DLL = "rpbnb.tmb")
  on.exit(TMB::openmp(n = previous, DLL = "rpbnb.tmb"), add = TRUE)
  supported <- as.integer(TMB::openmp(max = TRUE, DLL = "rpbnb.tmb")[[1L]])
  if (supported < 2L) skip("TMB runtime supports only one thread")

  reference_fn <- c(`1` = 23.1997247775494, `2` = 23.0779112497932)
  reference_gr <- list(
    `1` = c(-0.234213555458111, -0.644306804887814,
            -0.305971489013748, -0.377991335683169,
            0.0028883314805063, -0.00572813869692664,
            0.468873423307357, 0.365165705931843,
            0.739199683181961),
    `2` = c(-0.0442292322068991, -0.604280050235007,
            -0.668739747888821, -0.551604808939026,
            0.00734542308301755, -0.00647710160389995,
            0.443987374428227, 0.355502072666529,
            3.79146327117669)
  )

  for (family_code in c(1L, 2L)) {
    fixture <- copula_parallel_fixture(family_code)
    serial <- .make_rpbnb_tmb_object(
      data = fixture$data, parameters = fixture$parameters,
      map = fixture$map, n_cores = 1L
    )
    par <- serial$obj$par
    serial_fn <- serial$obj$fn(par)
    serial_gr <- serial$obj$gr(par)
    key <- as.character(family_code)
    expect_equal(serial_fn, unname(reference_fn[[key]]), tolerance = 1e-10)
    expect_equal(as.numeric(serial_gr), reference_gr[[key]], tolerance = 1e-8)

    parallel <- .make_rpbnb_tmb_object(
      data = fixture$data, parameters = fixture$parameters,
      map = fixture$map, n_cores = 2L
    )
    parallel_fn <- parallel$obj$fn(par)
    parallel_gr <- parallel$obj$gr(par)

    expect_equal(parallel_fn, serial_fn, tolerance = 1e-10)
    expect_equal(parallel_gr, serial_gr, tolerance = 1e-8)
    expect_true(all(is.finite(parallel_gr)))
  }
})

test_that("parallel benchmark reports timing and numerical agreement", {
  benchmark_path <- testthat::test_path(
    "..", "..", "inst", "benchmark_parallel.R"
  )
  expect_true(file.exists(benchmark_path))
  if (!file.exists(benchmark_path)) return(invisible())

  benchmark <- paste(readLines(benchmark_path, warn = FALSE), collapse = "\n")
  expect_match(benchmark, "realized", fixed = TRUE)
  expect_match(benchmark, "speedup", fixed = TRUE)
  expect_match(benchmark, "coef_diff", fixed = TRUE)
  expect_match(benchmark, "se_diff", fixed = TRUE)
  expect_match(benchmark, "objective_diff", fixed = TRUE)
  expect_match(benchmark, "gradient_diff", fixed = TRUE)
})

test_that("demo scripts use physical-core counts", {
  demo_paths <- testthat::test_path(
    "..", "..", "inst",
    c("fit_rpbnb_diff_copula.R", "fit_rpbnb_diff_famoye.R")
  )
  copula_demo <- paste(
    readLines(demo_paths[[1L]], warn = FALSE),
    collapse = "\n"
  )
  famoye_demo <- paste(
    readLines(demo_paths[[2L]], warn = FALSE),
    collapse = "\n"
  )

  expect_match(copula_demo, "detectCores(logical = FALSE)", fixed = TRUE)
  expect_match(copula_demo, "max(1L, min(4L", fixed = TRUE)
  expect_match(
    copula_demo,
    '.example_positive_integer("RPBNB_N_CORES", default_cores)',
    fixed = TRUE
  )
  expect_match(copula_demo, "n_cores     = example_cores", fixed = TRUE)

  expect_match(
    famoye_demo,
    "n_cores <- parallel::detectCores(logical = FALSE)",
    fixed = TRUE
  )
  expect_match(
    famoye_demo,
    "n_cores     = n_cores",
    fixed = TRUE
  )
  expect_false(grepl("min(2L", famoye_demo, fixed = TRUE))
})
