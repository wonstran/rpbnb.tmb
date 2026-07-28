# Write the dense truck example's captured reports as Markdown.
# Internal helper used by inst/truck_rpbnb_diff_famoye_dense.R.
.write_truck_results_markdown <- function(
    model_summary,
    marginal_effects,
    elasticities,
    results_dir = "results",
    timestamp = Sys.time()) {
  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

  file_stamp <- format(timestamp, "%Y-%m-%d-%H%M%S")
  generated_at <- format(timestamp, "%Y-%m-%d %H:%M:%S %Z")
  output_path <- file.path(
    results_dir,
    sprintf("results_%s.md", file_stamp)
  )

  report <- c(
    "# RP-BNB truck model results",
    "",
    sprintf("Generated: %s", generated_at),
    "",
    "## Model fit summary",
    "",
    "```text",
    model_summary,
    "```",
    "",
    "## Average marginal effects (AME)",
    "",
    "```text",
    marginal_effects,
    "```",
    "",
    "## Elasticities / semi-elasticities (AME)",
    "",
    "```text",
    elasticities,
    "```"
  )

  writeLines(report, output_path, useBytes = TRUE)
  invisible(output_path)
}
