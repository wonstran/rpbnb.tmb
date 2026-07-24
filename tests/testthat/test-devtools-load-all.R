test_that("devtools can load the source package", {
  skip_if_not_installed("devtools")

  package_root <- normalizePath(test_path("..", ".."))
  expect_no_error(devtools::load_all(package_root, quiet = TRUE))
})
