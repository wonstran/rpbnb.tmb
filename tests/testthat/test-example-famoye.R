famoye_example_text <- function() {
  path <- testthat::test_path(
    "..", "..", "inst", "fit_rpbnb_diff_famoye.R"
  )
  paste(readLines(path, warn = FALSE), collapse = "\n")
}

test_that("Famoye example has three editable workload settings", {
  script <- famoye_example_text()

  expect_match(script, "n_obs <- 5000L", fixed = TRUE)
  expect_match(script, "draws <- 400L", fixed = TRUE)
  expect_match(
    script,
    "n_cores <- 8",
    fixed = TRUE
  )
  expect_match(
    script,
    "if (is.na(n_cores)) n_cores <- 1L",
    fixed = TRUE
  )
})

test_that("Famoye example uses settings directly", {
  script <- famoye_example_text()

  expect_match(script, "min(n_obs, nrow(data))", fixed = TRUE)
  expect_match(script, "draws      = draws", fixed = TRUE)
  expect_match(script, "n_cores     = n_cores", fixed = TRUE)
  expect_match(script, "fit$parallel$realized", fixed = TRUE)
  expect_match(script, "fit$optimizer$convergence", fixed = TRUE)
  expect_match(script, "fit$sdreport$pdHess", fixed = TRUE)
})

test_that("Famoye example has no configuration helpers", {
  script <- famoye_example_text()

  expect_false(grepl(".famoye_example_", script, fixed = TRUE))
  expect_false(grepl("RPBNB_FAMOYE_", script, fixed = TRUE))
})
