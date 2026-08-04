# Write the dense truck example's captured reports as Markdown.
# Internal helper used by inst/truck_rpbnb_diff_famoye_dense.R.
#
# `dependence` and `method` are optional and default to NULL, which reproduces
# the original three-section report byte for byte -- the four Famoye/Frank/
# Clayton scripts in inst/ pass neither, and test-results-export.R pins the
# exact lines.  Supplying them adds the "Model information" header, which
# matters once several dependence structures are fitted to the same data into
# the same results/ directory: the reports are named only by timestamp, so
# without it nothing inside the file records which model produced it.
.write_truck_results_markdown <- function(
    model_summary,
    marginal_effects,
    elasticities,
    dependence = NULL,
    method = NULL,
    raw_scale = NULL,
    scaling = NULL,
    coef_orig_units = NULL,
    results_dir = "results",
    timestamp = Sys.time()) {
  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

  file_stamp <- format(timestamp, "%Y-%m-%d-%H%M%S")
  generated_at <- format(timestamp, "%Y-%m-%d %H:%M:%S %Z")
  output_path <- file.path(
    results_dir,
    sprintf("results_%s.md", file_stamp)
  )

  model_information <- character(0)
  if (!is.null(dependence) || !is.null(method)) {
    model_information <- c(
      "## Model information",
      "",
      if (!is.null(dependence)) {
        sprintf("- Dependence: %s", .dependence_label(dependence))
      },
      if (!is.null(method)) sprintf("- Method: %s", method),
      ""
    )
  }

  # The coefficient table in "Model fit summary" is per standard deviation of
  # each continuous predictor once the design was standardized.  When the
  # caller supplies the original-unit restatement (an affine transform of the
  # standardized estimates computed by the fit script), it is shown right
  # after the summary so the two are read together.
  coef_orig_units_section <- character(0)
  if (!is.null(coef_orig_units)) {
    coef_orig_units_section <- c(
      "",
      "## Coefficients in original units",
      "",
      paste(
        "The design was standardized before fitting, so the coefficient table",
        "above is per standard deviation of each continuous predictor.  The",
        "table below restates the same fit in the covariates' original units:",
        "continuous slopes divide by the scale, the intercept absorbs the",
        "centring shift, and binary 0/1 coefficients are unchanged.  Standard",
        "errors are delta-method on the full covariance.  These tables are for",
        "display; the fitted design itself remains standardized."
      ),
      "",
      "```text",
      .scaling_note(scaling),
      coef_orig_units,
      "```"
    )
  }

  # A standardized design makes the two sections above unreadable in opposite
  # ways, and silently: a scaled coefficient's effect is per standard deviation
  # rather than per unit, and a CENTRED regressor has sample mean zero, which
  # drives its elasticity's leading x-bar factor to zero.  The dense truck
  # report carried continuous elasticities of 1e-16 through 1e-18 for exactly
  # that reason -- printed to four decimals as 0.0000, which reads as "no
  # effect" rather than "this number means nothing".  So when the caller
  # standardized, the raw-unit tables are appended and the transform is
  # recorded, rather than leaving the reader to notice.
  raw_scale_section <- character(0)
  if (!is.null(raw_scale)) {
    scaling_note <- .scaling_note(scaling)
    raw_scale_section <- c(
      "",
      "## Marginal effects and elasticities in original units",
      "",
      paste(
        "The design was standardized before fitting, so the two sections",
        "above are per standard deviation, and the elasticities of the",
        "centred regressors are identically zero (their sample mean is zero)."
      ),
      paste(
        "The tables below restate the same fit in the covariates' original",
        "units. Binary regressors were not transformed, so their effects and",
        "semi-elasticities are unchanged."
      ),
      "",
      "```text",
      scaling_note,
      raw_scale,
      "```"
    )
  }

  report <- c(
    "# RP-BNB truck model results",
    "",
    sprintf("Generated: %s", generated_at),
    "",
    model_information,
    "## Model fit summary",
    "",
    "```text",
    model_summary,
    "```",
    coef_orig_units_section,
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
    "```",
    raw_scale_section
  )

  writeLines(report, output_path, useBytes = TRUE)
  invisible(output_path)
}

# Standardization note shared by the two "original units" sections: the raw
# center/scale table for the continuous predictors.  NULL when the transform
# was not supplied (c() drops it, so the surrounding section still renders).
.scaling_note <- function(scaling) {
  if (is.null(scaling)) return(NULL)
  tab <- do.call(rbind, lapply(names(scaling), function(v) {
    sprintf("  %-12s center = %12.4f   scale = %12.4f",
            v, scaling[[v]][["center"]], scaling[[v]][["scale"]])
  }))
  c("Standardization applied before fitting:", "", as.vector(tab), "")
}

# A human-readable name for whatever `fit$dependence` holds.  deparse() is what
# print.rpbnb_tmb_fit() uses, and on a copula object it renders the whole
# structure() call -- correct, unreadable, and three lines long in a bulleted
# list.  Anything unrecognised still falls through to deparse() rather than
# being dropped, since a wrong label is worse than an ugly one.
.dependence_label <- function(dependence) {
  if (inherits(dependence, "rpbnb_copula")) {
    family <- dependence$family
    pretty <- switch(family,
      normal = "Gaussian copula",
      frank = "Frank copula",
      kimeldorf = "Clayton copula",
      paste(family, "copula")
    )
    return(sprintf("%s (%s)", pretty, family))
  }
  if (is.character(dependence) && length(dependence) == 1L) {
    return(switch(dependence,
      famoye = "Famoye/Sarmanov",
      independence = "independence",
      dependence
    ))
  }
  paste(deparse(dependence), collapse = " ")
}
