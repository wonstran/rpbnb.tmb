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
    lamLo = -1, lamHi = 1, est_method = 0L
  )
  parameters <- list(
    beta1 = 0, beta2 = 0,
    log_sd1 = numeric(0), log_sd2 = numeric(0),
    log_m1 = log(0.5), log_m2 = log(0.5), z_dep = 0,
    u1 = matrix(0, 2L, 0L), u2 = matrix(0, 2L, 0L)
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
    lamLo = 0, lamHi = 0, est_method = 0L
  )
  parameters <- list(
    beta1 = c(0.1, 0.2), beta2 = c(-0.1, -0.15),
    log_sd1 = log(0.2), log_sd2 = log(0.25),
    log_m1 = log(0.6), log_m2 = log(0.7),
    z_dep = if (family_code == 1L) 1.5 else atanh(0.25),
    u1 = matrix(0, 8L, 1L), u2 = matrix(0, 8L, 1L)
  )
  # SML fixture: latents are tape constants at zero, matching the map applied
  # in fit_rpbnb_tmb() for est_method == 0 (see R/fit_rpbnb_tmb.R). Fixed here
  # so obj$par and the hardcoded reference_gr values below are unaffected.
  list(data = data, parameters = parameters,
       map = list(u1 = factor(rep(NA_integer_, 8L)),
                  u2 = factor(rep(NA_integer_, 8L))))
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

test_that("rpbnb_tmb_control validates memory guardrails", {
  expect_error(rpbnb_tmb_control(max_threads = 0), "max_threads")
  expect_error(rpbnb_tmb_control(max_threads = 1.5), "max_threads")
  expect_error(rpbnb_tmb_control(max_threads = "2"), "max_threads")
  expect_error(rpbnb_tmb_control(max_workload = 0), "max_workload")
  expect_error(rpbnb_tmb_control(max_workload = NA_real_), "max_workload")
  expect_error(rpbnb_tmb_control(max_workload = "large"), "max_workload")

  control <- rpbnb_tmb_control(max_threads = 3, max_workload = Inf)
  expect_identical(control$max_threads, 3L)
  expect_identical(control$max_workload, Inf)
  expect_error(rpbnb_tmb_control(parallel_tape = NA), "parallel_tape")
  expect_error(rpbnb_tmb_control(parallel_tape = 1), "parallel_tape")
  expect_false(rpbnb_tmb_control()$parallel_tape)
})

test_that("thread configuration respects the memory-aware cap", {
  previous <- TMB::openmp(DLL = "rpbnb.tmb")
  on.exit(TMB::openmp(n = previous, DLL = "rpbnb.tmb"), add = TRUE)
  supported <- as.integer(TMB::openmp(max = TRUE, DLL = "rpbnb.tmb")[[1L]])
  if (supported < 2L) skip("TMB runtime supports only one thread")

  expect_warning(
    realized <- .configure_tmb_threads(
      n_cores = min(4L, supported),
      max_threads = 1L,
      DLL = "rpbnb.tmb"
    ),
    "max_threads"
  )
  expect_identical(realized, 1L)
})

test_that("TMB tape construction is sequential unless explicitly enabled", {
  previous <- TMB::config(DLL = "rpbnb.tmb")$tape.parallel
  on.exit(
    TMB::config(tape.parallel = previous, DLL = "rpbnb.tmb"),
    add = TRUE
  )

  .configure_tmb_threads(
    n_cores = 1L, max_threads = 1L,
    parallel_tape = FALSE, DLL = "rpbnb.tmb"
  )
  expect_identical(TMB::config(DLL = "rpbnb.tmb")$tape.parallel, 0L)

  .configure_tmb_threads(
    n_cores = 1L, max_threads = 1L,
    parallel_tape = TRUE, DLL = "rpbnb.tmb"
  )
  expect_identical(TMB::config(DLL = "rpbnb.tmb")$tape.parallel, 1L)
})

test_that("the workload guard rejects exactly at its limit and opts out", {
  limit <- 100 * 10 * TAPE_CALIBRATION$family_weight[["gaussian"]]
  expect_error(
    .check_tmb_workload(
      n = 100L, draws = 10L, family_code = 2L,
      max_workload = limit - 1
    ),
    "max_workload"
  )
  expect_no_error(
    .check_tmb_workload(
      n = 100L, draws = 10L, family_code = 2L, max_workload = limit
    )
  )
  expect_no_error(
    .check_tmb_workload(
      n = .Machine$integer.max, draws = .Machine$integer.max,
      family_code = 2L, max_workload = Inf
    )
  )
})

test_that("concurrent taping, and only that, multiplies the budget", {
  weights <- TAPE_CALIBRATION$family_weight
  families <- c(independence = -1L, famoye = 0L, frank = 1L,
                gaussian = 2L, clayton = 3L)

  # Sequential taping: the thread count is irrelevant to the guard.
  for (nm in names(families)) {
    limit <- 1000 * 400 * weights[[nm]]
    expect_no_error(.check_tmb_workload(
      n = 1000L, draws = 400L, family_code = families[[nm]],
      max_workload = limit, n_threads = 8L, parallel_tape = FALSE
    ))
    expect_error(.check_tmb_workload(
      n = 1000L, draws = 400L, family_code = families[[nm]],
      max_workload = limit - 1, n_threads = 8L, parallel_tape = FALSE
    ), "max_workload")
  }

  # Concurrent taping: n_threads full-size recordings are alive at once.
  concurrent <- 1000 * 400 * weights[["famoye"]] * 8
  expect_error(.check_tmb_workload(
    n = 1000L, draws = 400L, family_code = 0L,
    max_workload = concurrent - 1, n_threads = 8L, parallel_tape = TRUE
  ), "max_workload")
  expect_no_error(.check_tmb_workload(
    n = 1000L, draws = 400L, family_code = 0L,
    max_workload = concurrent, n_threads = 8L, parallel_tape = TRUE
  ))
})

# Asserting that one particular dataset size lands under the limit locks in a
# coincidence: it passes at 96.8% of budget and fails on a slightly larger
# dataset that costs no more memory to fit.  Assert the guard's properties and
# a minimum headroom for realistic work instead.
test_that("the workload guard is linear in its inputs", {
  # A passing check returns the weighted workload it computed.
  cost <- function(n, draws, family_code, ...) {
    .check_tmb_workload(
      n = n, draws = draws, family_code = family_code,
      max_workload = 1e12, ...
    )
  }
  base <- cost(1000L, 100L, 0L)
  expect_equal(base, 1e5)
  expect_equal(cost(2000L, 100L, 0L), 2 * base)
  expect_equal(cost(1000L, 200L, 0L), 2 * base)
  # Family weights are measured, not assumed, and they differ: Frank costs
  # about 2.9x Famoye per unit.  The guard must read them from the single
  # calibration source rather than hard-coding any of them.
  weights <- TAPE_CALIBRATION$family_weight
  expect_equal(cost(1000L, 100L, -1L), weights[["independence"]] * base)
  expect_equal(cost(1000L, 100L, 1L), weights[["frank"]] * base)
  expect_equal(cost(1000L, 100L, 2L), weights[["gaussian"]] * base)
  expect_equal(cost(1000L, 100L, 3L), weights[["clayton"]] * base)
  expect_gt(weights[["frank"]], 2)
})

test_that("threads multiply the workload only under concurrent taping", {
  cost <- function(...) {
    .check_tmb_workload(
      n = 1000L, draws = 100L, family_code = 0L,
      max_workload = 1e12, ...
    )
  }
  expect_equal(cost(n_threads = 8L, parallel_tape = FALSE), 1e5)
  expect_equal(cost(n_threads = 8L, parallel_tape = TRUE), 8e5)
})

test_that("the control default calls through to rpbnb_tmb_max_workload()", {
  # An earlier version of this test asserted the default equals the fixed
  # calibration constant directly. That stopped being the general case once
  # the default started auto-detecting memory; this test now asserts the
  # actual current contract instead, on both paths that default can take.
  testthat::local_mocked_bindings(
    .detect_available_memory_gib = function() 10
  )
  expect_identical(
    rpbnb_tmb_control()$max_workload,
    rpbnb_tmb_max_workload()
  )
})

test_that("the control default falls back to the calibration constant when detection fails", {
  testthat::local_mocked_bindings(
    .detect_available_memory_gib = function() NA_real_
  )
  expect_warning(
    result <- rpbnb_tmb_control()$max_workload,
    "Could not detect"
  )
  expect_identical(result, .calibration_default_workload())
})

test_that("the documentation is generated from the calibration", {
  # `@eval .calibration_doc()` means the Rd text cannot drift from the
  # constants: it is produced from them at roxygenise time.  Assert the
  # generator actually reports the shipped values, so a silently edited
  # constant shows up as changed documentation rather than a stale claim.
  doc <- .calibration_doc()
  expect_match(doc, as.character(TAPE_CALIBRATION$peak_bytes_per_unit),
               fixed = TRUE)
  expect_match(doc, as.character(TAPE_CALIBRATION$tape_bytes_per_unit),
               fixed = TRUE)
  expect_match(doc, format(.calibration_default_workload(),
                           scientific = FALSE, big.mark = ","), fixed = TRUE)
  expect_match(doc, "rpbnb_tmb_max_workload", fixed = TRUE)
  expect_match(doc, "auto-detect", fixed = TRUE)
  expect_match(doc, "frank 3.6", fixed = TRUE)
  # Every family the guard can be asked about must be a measured one.
  expect_setequal(
    names(TAPE_CALIBRATION$family_weight), TAPE_CALIBRATION$measured_families
  )
})

test_that("the committed calibration data backs the published constants", {
  path <- system.file("extdata", "memory_calibration.csv",
                      package = "rpbnb.tmb")
  skip_if(path == "", "calibration data not installed")
  raw <- utils::read.csv(path)

  # Every family carrying a weight was actually measured.
  expect_setequal(
    sort(unique(raw$family)), sort(TAPE_CALIBRATION$measured_families)
  )
  # And the published weights are no smaller than the observed ratios.
  # PEAK, not tape: the guard budgets peak working set, so weights validated
  # against retained tape would not be validated against the thing they scale.
  # The two differ enough to matter -- Frank retains 2.87x Famoye but peaks at
  # 3.53x, which is exactly how a tape-derived 2.9 came to under-budget it.
  famoye <- aggregate(peak_mib ~ units, subset(raw, family == "famoye"), mean)
  for (fam in setdiff(unique(raw$family), "famoye")) {
    cells <- subset(raw, family == fam)
    ratios <- vapply(seq_len(nrow(cells)), function(i) {
      base <- famoye$peak_mib[famoye$units == cells$units[i]]
      if (length(base)) cells$peak_mib[i] / base else NA_real_
    }, numeric(1))
    expect_gte(
      TAPE_CALIBRATION$family_weight[[fam]], max(ratios, na.rm = TRUE) - 1e-8
    )
  }
})

test_that("the README's restated calibration figures match the calibration", {
  # README.md hand-copies the weights, the per-unit slope, the default
  # workload and the budget. Nothing read it, so changing its Frank weight to
  # 9.9 left every test and the manual verifier green -- the same class of
  # silent drift the reference-manual verifier was extended to catch.
  readme <- testthat::test_path("..", "..", "README.md")
  skip_if_not(file.exists(readme), "README.md not available (installed pkg)")
  text <- paste(readLines(readme, warn = FALSE), collapse = "\n")

  # Compare NUMBERS, not formatted strings: the README writes `7e5` where
  # format() gives "7e+05", and `1.0` where format() gives "1". A string
  # comparison fails on those without any drift having occurred.
  extract <- function(pattern) {
    hit <- regmatches(text, regexpr(pattern, text, perl = TRUE))
    expect_length(hit, 1L)
    as.numeric(hit)
  }

  # Weights, from the table row: | frank | **3.6** | 3.530 | 2.867 |
  for (fam in names(TAPE_CALIBRATION$family_weight)) {
    row <- regmatches(
      text, regexpr(paste0("\\|\\s*", fam, "\\s*\\|[^|]*\\|"), text)
    )
    expect_length(row, 1L)
    expect_equal(as.numeric(gsub("[^0-9.]", "", row)),
                 TAPE_CALIBRATION$family_weight[[fam]],
                 info = paste("README weight for", fam))
  }

  expect_equal(extract("(?<=The default of `)[^`]+(?=` units)"),
               .calibration_default_workload())
  expect_equal(extract("(?<=about )[0-9.]+(?= GiB of peak memory)"),
               TAPE_CALIBRATION$budget_gib)
  # Prose rounds the ratio to one decimal, so compare at that precision.
  expect_equal(extract("(?<=about )[0-9.]+(?=x the retained tape)"),
               TAPE_CALIBRATION$peak_over_tape_and_eval, tolerance = 0.02)
})

test_that("the peak calibration constants are reproducible from the data", {
  # peak_bytes_per_unit scales every budget the guard enforces, so it must be
  # recoverable from the committed measurements rather than asserted. Same
  # large-workload slope the tape figure uses, applied to peak_mib.
  path <- system.file("extdata", "memory_calibration.csv",
                      package = "rpbnb.tmb")
  skip_if(path == "", "calibration data not installed")
  raw <- utils::read.csv(path)
  famoye <- aggregate(cbind(tape_mib, peak_mib) ~ units,
                      subset(raw, family == "famoye"), mean)
  top <- utils::tail(famoye[order(famoye$units), ], 3L)

  peak_slope <- coef(stats::lm(peak_mib ~ units, top))[["units"]] * 1024^2
  expect_equal(TAPE_CALIBRATION$peak_bytes_per_unit, peak_slope,
               tolerance = 0.02)

  tape_slope <- coef(stats::lm(tape_mib ~ units, top))[["units"]] * 1024^2
  expect_equal(TAPE_CALIBRATION$tape_bytes_per_unit, tape_slope,
               tolerance = 0.02)

  # The published ratio divides by tape + eval; naming it "retained" alone
  # would be off by more than 40%.
  expect_equal(
    TAPE_CALIBRATION$peak_over_tape_and_eval,
    median(raw$peak_mib / (raw$tape_mib + raw$eval_mib)),
    tolerance = 0.01
  )
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
  fixture <- parallel_tmb_fixture()
  previous <- TMB::openmp(DLL = "rpbnb.tmb")
  on.exit(TMB::openmp(n = previous, DLL = "rpbnb.tmb"), add = TRUE)
  supported <- as.integer(TMB::openmp(max = TRUE, DLL = "rpbnb.tmb")[[1L]])
  requested <- if (supported >= 2L) 2L else 1L
  configured <- .make_rpbnb_tmb_object(
    data = fixture$data,
    parameters = fixture$parameters,
    map = fixture$map,
    n_cores = requested
  )

  expect_identical(configured$n_cores, requested)
  reported <- configured$obj$report()$openmp_compiled
  expect_true(reported %in% c(0L, 1L))
  if (supported >= 2L) expect_identical(reported, 1L)
})

test_that("serial and parallel copula objectives and gradients agree", {
  previous <- TMB::openmp(DLL = "rpbnb.tmb")
  on.exit(TMB::openmp(n = previous, DLL = "rpbnb.tmb"), add = TRUE)
  supported <- as.integer(TMB::openmp(max = TRUE, DLL = "rpbnb.tmb")[[1L]])
  if (supported < 2L) skip("TMB runtime supports only one thread")

  reference_fn <- c(`1` = 23.1990464774657, `2` = 23.0779112497932)
  reference_gr <- list(
    `1` = c(-0.234009598804351, -0.644262309598928,
            -0.30606552853166, -0.378131107639618,
            0.00289799326474705, -0.00572295641355119,
            0.468916002251157, 0.365254841128196,
            0.737716614784872),
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
    serial_he <- serial$obj$he(par)
    key <- as.character(family_code)
    expect_equal(serial_fn, unname(reference_fn[[key]]), tolerance = 1e-10)
    expect_equal(as.numeric(serial_gr), reference_gr[[key]], tolerance = 1e-8)
    if (family_code == 2L) {
      expect_equal(
        as.numeric(serial_gr),
        numDeriv::grad(serial$obj$fn, par),
        tolerance = 1e-6
      )
    }

    parallel <- .make_rpbnb_tmb_object(
      data = fixture$data, parameters = fixture$parameters,
      map = fixture$map, n_cores = 2L
    )
    parallel_fn <- parallel$obj$fn(par)
    parallel_gr <- parallel$obj$gr(par)
    parallel_he <- parallel$obj$he(par)

    expect_equal(parallel_fn, serial_fn, tolerance = 1e-10)
    expect_equal(parallel_gr, serial_gr, tolerance = 1e-8)
    expect_equal(parallel_he, serial_he, tolerance = 1e-6)
    expect_true(all(is.finite(parallel_gr)))
    expect_true(all(is.finite(parallel_he)))
  }
})

test_that("Frank and NB kernels stay finite at extreme working parameters", {
  frank <- copula_parallel_fixture(1L)
  frank_obj <- .make_rpbnb_tmb_object(
    data = frank$data, parameters = frank$parameters,
    map = frank$map, n_cores = 1L
  )$obj
  frank_par <- frank_obj$par
  frank_par["z_dep"] <- -500
  expect_true(is.finite(frank_obj$fn(frank_par)))
  expect_true(all(is.finite(frank_obj$gr(frank_par))))

  nb <- parallel_tmb_fixture()
  nb_obj <- .make_rpbnb_tmb_object(
    data = nb$data, parameters = nb$parameters,
    map = nb$map, n_cores = 1L
  )$obj
  nb_par <- nb_obj$par
  nb_par[c("beta1", "beta2")] <- -100
  expect_true(is.finite(nb_obj$fn(nb_par)))
  expect_true(all(is.finite(nb_obj$gr(nb_par))))
})
