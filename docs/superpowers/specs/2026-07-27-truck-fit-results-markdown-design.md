# Truck Fit Results Markdown Export Design

## Goal

Update `inst/truck_rpbnb_diff_famoye_dense.R` so every successful run exports
the model-fit summary, average marginal effects, and elasticities to:

`results/results_YYYY-MM-DD-HHMMSS.md`

The existing console report must remain available.

## Design

Capture each requested report section when it is already evaluated:

1. Capture the printed output of `summary(fit)`.
2. Capture the printed output and returned value of
   `rpbnb_tmb_marginal_effects(fit, which = "both")`.
3. Capture the printed output and returned value of
   `rpbnb_tmb_elasticities(fit, which = "both")`.

Echo the captured text back to the console so the script retains its current
interactive behavior. Do not call either diagnostic function a second time,
because both perform nontrivial post-estimation calculations.

After all three sections succeed, create `results/` if necessary, generate a
timestamp with `format(Sys.time(), "%Y-%m-%d-%H%M%S")`, and write one Markdown
document. The document contains:

- a title;
- the generation timestamp;
- a `Model fit summary` section;
- an `Average marginal effects (AME)` section;
- an `Elasticities / semi-elasticities (AME)` section.

The captured console tables are placed in fenced `text` blocks. This preserves
all existing labels and formatting without adding a Markdown-table dependency
or reshaping heterogeneous result tables.

## Error Handling

The results directory is created with `recursive = TRUE`. If directory creation
or file writing fails, the script stops with the base R error rather than
silently claiming an export was produced. Export occurs only after all requested
report sections have completed, so a successful file is not partially populated.

## Verification

Add a focused static test for the example script that verifies:

- the three report calls are captured rather than recomputed for export;
- captured output is echoed to the console;
- the results directory is created;
- the filename uses the required date-time pattern;
- all three Markdown headings and fenced output blocks are written.

Run that focused test, parse the updated R script, and inspect the final diff.
The full truck model is intentionally not executed because it requires the
large bundled dataset, 300 simulation draws, 12 cores, and a high memory budget.
