test_that(".detect_available_memory_gib() never errors and returns a sane value", {
  result <- .detect_available_memory_gib()

  expect_length(result, 1L)
  expect_true(is.numeric(result))
  # Either a real reading or an explicit "could not detect" signal -- never
  # NaN, never negative, never an error propagating out of this function.
  expect_true(is.na(result) || (is.finite(result) && result >= 0))
})
