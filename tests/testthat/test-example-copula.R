copula_example_text <- function() {
  path <- testthat::test_path(
    "..", "..", "inst", "fit_rpbnb_diff_copula.R"
  )
  paste(readLines(path, warn = FALSE), collapse = "\n")
}

test_that("copula example has three editable workload settings", {
  script <- copula_example_text()

  expect_match(script, "n_obs <- 500L", fixed = TRUE)
  expect_match(script, "draws <- 20L", fixed = TRUE)
  expect_match(
    script,
    "n_cores <- parallel::detectCores(logical = FALSE)",
    fixed = TRUE
  )
  expect_match(
    script,
    "if (is.na(n_cores)) n_cores <- 1L",
    fixed = TRUE
  )
})

test_that("copula example uses settings directly", {
  script <- copula_example_text()

  expect_match(script, "min(n_obs, nrow(data))", fixed = TRUE)
  expect_match(script, "draws      = draws", fixed = TRUE)
  expect_match(script, "n_cores     = n_cores", fixed = TRUE)
  expect_match(script, "fit$parallel$realized", fixed = TRUE)
  expect_match(script, "fit$optimizer$convergence", fixed = TRUE)
  expect_match(script, "fit$sdreport$pdHess", fixed = TRUE)
})

test_that("copula example has no configuration helpers", {
  script <- copula_example_text()

  expect_false(grepl(".example_", script, fixed = TRUE))
  expect_false(grepl("RPBNB_", script, fixed = TRUE))
})
