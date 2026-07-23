load_copula_example_helpers <- function() {
  path <- testthat::test_path(
    "..", "..", "inst", "fit_rpbnb_diff_copula.R"
  )
  expressions <- parse(path)
  helper_names <- c(
    ".example_positive_integer",
    ".example_observation_count",
    ".example_fit_diagnostics"
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

test_that("copula example reads positive integer configuration", {
  helpers <- load_copula_example_helpers()
  old <- Sys.getenv("RPBNB_DRAWS", unset = NA_character_)
  on.exit(
    if (is.na(old)) {
      Sys.unsetenv("RPBNB_DRAWS")
    } else {
      Sys.setenv(RPBNB_DRAWS = old)
    },
    add = TRUE
  )

  Sys.unsetenv("RPBNB_DRAWS")
  expect_identical(
    helpers$.example_positive_integer("RPBNB_DRAWS", 20L),
    20L
  )

  Sys.setenv(RPBNB_DRAWS = "40")
  expect_identical(
    helpers$.example_positive_integer("RPBNB_DRAWS", 20L),
    40L
  )
})

test_that("copula example rejects invalid integer configuration", {
  helpers <- load_copula_example_helpers()
  old <- Sys.getenv("RPBNB_N_CORES", unset = NA_character_)
  on.exit(
    if (is.na(old)) {
      Sys.unsetenv("RPBNB_N_CORES")
    } else {
      Sys.setenv(RPBNB_N_CORES = old)
    },
    add = TRUE
  )

  for (value in c("abc", "1.5", "0", "-1", "Inf", "1,2")) {
    Sys.setenv(RPBNB_N_CORES = value)
    expect_error(
      helpers$.example_positive_integer("RPBNB_N_CORES", 1L),
      "RPBNB_N_CORES"
    )
  }
})

test_that("copula example caps observations to available rows", {
  helpers <- load_copula_example_helpers()

  expect_identical(
    helpers$.example_observation_count(50L, 100L),
    50L
  )
  expect_message(
    actual <- helpers$.example_observation_count(150L, 100L),
    "150.*100"
  )
  expect_identical(actual, 100L)
})

test_that("copula example reports parallel and convergence diagnostics", {
  helpers <- load_copula_example_helpers()
  fit <- list(
    parallel = list(requested = 4L, realized = 2L),
    optimizer = list(convergence = 0L, message = "relative convergence"),
    sdreport = list(pdHess = TRUE)
  )

  output <- paste(
    capture.output(helpers$.example_fit_diagnostics(fit, 1.25)),
    collapse = "\n"
  )

  expect_match(output, "1.25 s")
  expect_match(output, "requested=4")
  expect_match(output, "realized=2")
  expect_match(output, "code=0")
  expect_match(output, "relative convergence")
  expect_match(output, "positive-definite Hessian: yes")
})
