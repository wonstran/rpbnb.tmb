test_that("truck results writer creates the timestamped Markdown report", {
  results_dir <- file.path(tempdir(), "truck-results")
  unlink(results_dir, recursive = TRUE)

  output_path <- rpbnb.tmb:::.write_truck_results_markdown(
    model_summary = c("Summary: fitted model", "Log-likelihood: -12.5"),
    marginal_effects = c(
      "--- Marginal effects (equation 1, AME) ---",
      "x  1.25"
    ),
    elasticities = c(
      "--- Elasticities (equation 1, AME) ---",
      "x  0.50"
    ),
    results_dir = results_dir,
    timestamp = as.POSIXct("2026-07-27 14:35:09", tz = "UTC")
  )

  expect_identical(
    output_path,
    file.path(results_dir, "results_2026-07-27-143509.md")
  )
  expect_true(file.exists(output_path))

  report <- readLines(output_path, warn = FALSE)
  expect_identical(report, c(
    "# RP-BNB truck model results",
    "",
    "Generated: 2026-07-27 14:35:09 UTC",
    "",
    "## Model fit summary",
    "",
    "```text",
    "Summary: fitted model",
    "Log-likelihood: -12.5",
    "```",
    "",
    "## Average marginal effects (AME)",
    "",
    "```text",
    "--- Marginal effects (equation 1, AME) ---",
    "x  1.25",
    "```",
    "",
    "## Elasticities / semi-elasticities (AME)",
    "",
    "```text",
    "--- Elasticities (equation 1, AME) ---",
    "x  0.50",
    "```"
  ))
})

test_that("shipped truck examples remain valid R syntax", {
  # system.file() rather than test_path("..", "..", "inst", ...): installation
  # flattens inst/, so the relative path resolves only in a source tree and
  # turns this into a release-check failure. system.file() resolves the source
  # inst/ under devtools::load_all() as well, so the assertion stays live in
  # both development and installed-package testing.
  # Discovered rather than listed. A hand-written list silently shrinks its own
  # coverage: this test previously named two of the four shipped truck scripts
  # while claiming in its title to cover them all, so a syntax regression in
  # truck_rpbnb_diff_famoye_laplace.R or truck_rpbnb_diff_frank_laplace.R would
  # have gone unnoticed. The pattern is anchored and narrow so it cannot start
  # sweeping up unrelated files.
  root <- system.file(package = "rpbnb.tmb", mustWork = TRUE)
  scripts <- list.files(root, pattern = "^truck_[A-Za-z0-9_]+\\.R$")
  # Guard the discovery itself: if the pattern or the layout ever stops
  # matching, this fails loudly instead of vacuously passing over zero files.
  expect_gte(length(scripts), 4L)
  expect_true(all(c("truck_rpbnb_diff_famoye_dense.R",
                    "truck_rpbnb_diff_famoye_laplace.R",
                    "truck_rpbnb_diff_frank_laplace.R",
                    "truck_rpbnb_diff_kimeldorf_laplace.R") %in% scripts))

  for (nm in scripts) {
    script <- system.file(nm, package = "rpbnb.tmb", mustWork = TRUE)
    # parse() raises on invalid syntax, and its own message names the file, so
    # a syntax error identifies itself; `label` additionally names the file on
    # the (unlikely) path where parse succeeds but returns something odd.
    # expect_no_error() is not used here because it forbids extra arguments.
    expect_true(is.expression(parse(file = script)), label = nm)
  }
})
