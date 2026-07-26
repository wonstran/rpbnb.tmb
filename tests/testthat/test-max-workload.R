test_that(".detect_available_memory_gib() never errors and returns a sane value", {
  result <- .detect_available_memory_gib()

  expect_length(result, 1L)
  expect_true(is.numeric(result))
  # Either a real reading or an explicit "could not detect" signal -- never
  # NaN, never negative, never an error propagating out of this function.
  expect_true(is.na(result) || (is.finite(result) && result >= 0))
})

test_that("an explicit budget_gib is pure arithmetic with no fraction applied", {
  expect_equal(
    rpbnb_tmb_max_workload(budget_gib = 8),
    .calibration_default_workload()
  )
  expect_equal(
    rpbnb_tmb_max_workload(budget_gib = 16),
    signif(16 * 1024^3 / TAPE_CALIBRATION$peak_bytes_per_unit, 1)
  )
  # fraction is ignored entirely on this path -- not even applied silently.
  expect_equal(
    rpbnb_tmb_max_workload(budget_gib = 16, fraction = 0.1),
    rpbnb_tmb_max_workload(budget_gib = 16, fraction = 0.9)
  )
})

test_that("detected memory is discounted by fraction", {
  testthat::local_mocked_bindings(
    .detect_available_memory_gib = function() 10
  )
  expect_equal(
    rpbnb_tmb_max_workload(),
    signif(0.8 * 10 * 1024^3 / TAPE_CALIBRATION$peak_bytes_per_unit, 1)
  )
  expect_equal(
    rpbnb_tmb_max_workload(fraction = 0.5),
    signif(0.5 * 10 * 1024^3 / TAPE_CALIBRATION$peak_bytes_per_unit, 1)
  )
})

test_that("failed detection warns and falls back to the calibration default", {
  testthat::local_mocked_bindings(
    .detect_available_memory_gib = function() NA_real_
  )
  expect_warning(
    result <- rpbnb_tmb_max_workload(),
    "Could not detect"
  )
  expect_equal(result, .calibration_default_workload())
})

test_that("fraction is validated", {
  expect_error(rpbnb_tmb_max_workload(fraction = 0), "fraction")
  expect_error(rpbnb_tmb_max_workload(fraction = 1.5), "fraction")
  expect_error(rpbnb_tmb_max_workload(fraction = NA_real_), "fraction")
  expect_error(rpbnb_tmb_max_workload(fraction = c(0.5, 0.5)), "fraction")
})

test_that("budget_gib is validated when supplied", {
  expect_error(rpbnb_tmb_max_workload(budget_gib = "8"), "budget_gib")
  expect_error(rpbnb_tmb_max_workload(budget_gib = -1), "budget_gib")
  expect_error(rpbnb_tmb_max_workload(budget_gib = 0), "budget_gib")
  expect_error(rpbnb_tmb_max_workload(budget_gib = Inf), "budget_gib")
  expect_error(rpbnb_tmb_max_workload(budget_gib = c(8, 16)), "budget_gib")
})
