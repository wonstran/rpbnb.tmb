load_famoye_example_helpers <- function() {
  path <- testthat::test_path(
    "..", "..", "inst", "fit_rpbnb_diff_famoye.R"
  )
  expressions <- parse(path)
  helper_names <- c(
    ".famoye_example_positive_integer",
    ".famoye_example_observation_count",
    ".famoye_example_fit_diagnostics",
    ".famoye_example_with_memory_guidance"
  )
  helper_env <- new.env(parent = baseenv())

  for (expression in expressions) {
    if (is.call(expression) &&
        identical(expression[[1L]], as.name("<-")) &&
        as.character(expression[[2L]]) %in% helper_names) {
      eval(expression, envir = helper_env)
    }
  }
  helper_env
}

test_that("Famoye example reads positive integer configuration", {
  helpers <- load_famoye_example_helpers()
  old <- Sys.getenv("RPBNB_FAMOYE_DRAWS", unset = NA_character_)
  on.exit(
    if (is.na(old)) {
      Sys.unsetenv("RPBNB_FAMOYE_DRAWS")
    } else {
      Sys.setenv(RPBNB_FAMOYE_DRAWS = old)
    },
    add = TRUE
  )

  Sys.unsetenv("RPBNB_FAMOYE_DRAWS")
  expect_identical(
    helpers$.famoye_example_positive_integer(
      "RPBNB_FAMOYE_DRAWS", 100L
    ),
    100L
  )

  Sys.setenv(RPBNB_FAMOYE_DRAWS = "40")
  expect_identical(
    helpers$.famoye_example_positive_integer(
      "RPBNB_FAMOYE_DRAWS", 100L
    ),
    40L
  )
})

test_that("Famoye example rejects invalid integer configuration", {
  helpers <- load_famoye_example_helpers()
  old <- Sys.getenv("RPBNB_FAMOYE_N_CORES", unset = NA_character_)
  on.exit(
    if (is.na(old)) {
      Sys.unsetenv("RPBNB_FAMOYE_N_CORES")
    } else {
      Sys.setenv(RPBNB_FAMOYE_N_CORES = old)
    },
    add = TRUE
  )

  for (value in c("abc", "1.5", "0", "-1", "Inf", "1,2")) {
    Sys.setenv(RPBNB_FAMOYE_N_CORES = value)
    expect_error(
      helpers$.famoye_example_positive_integer(
        "RPBNB_FAMOYE_N_CORES", 1L
      ),
      "RPBNB_FAMOYE_N_CORES"
    )
  }
})

test_that("Famoye example caps observations to available rows", {
  helpers <- load_famoye_example_helpers()

  expect_identical(
    helpers$.famoye_example_observation_count(50L, 100L),
    50L
  )
  expect_message(
    actual <- helpers$.famoye_example_observation_count(150L, 100L),
    "150.*100"
  )
  expect_identical(actual, 100L)
})

test_that("Famoye example reports fit diagnostics", {
  helpers <- load_famoye_example_helpers()
  fit <- list(
    parallel = list(requested = 2L, realized = 2L),
    optimizer = list(convergence = 0L, message = "relative convergence"),
    sdreport = list(pdHess = TRUE)
  )
  output <- paste(
    capture.output(
      helpers$.famoye_example_fit_diagnostics(fit, 1.25)
    ),
    collapse = "\n"
  )

  expect_match(output, "1.25 s")
  expect_match(output, "requested=2")
  expect_match(output, "realized=2")
  expect_match(output, "code=0")
  expect_match(output, "positive-definite Hessian: yes")
})

test_that("Famoye example translates bad_alloc into workload guidance", {
  helpers <- load_famoye_example_helpers()
  expect_error(
    helpers$.famoye_example_with_memory_guidance(
      stop("Caught exception 'std::bad_alloc'", call. = FALSE)
    ),
    paste0(
      "Restart R.*RPBNB_FAMOYE_N_OBS.*",
      "RPBNB_FAMOYE_DRAWS.*RPBNB_FAMOYE_N_CORES"
    )
  )
})

test_that("Famoye example preserves unrelated error conditions", {
  helpers <- load_famoye_example_helpers()
  original <- structure(
    list(message = "unrelated failure", call = NULL),
    class = c("famoye_example_test_error", "error", "condition")
  )
  caught <- tryCatch(
    helpers$.famoye_example_with_memory_guidance(stop(original)),
    error = identity
  )

  expect_s3_class(caught, "famoye_example_test_error")
  expect_identical(conditionMessage(caught), "unrelated failure")
})
