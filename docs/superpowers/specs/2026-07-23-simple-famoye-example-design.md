# Simple Famoye Example Design

## Goal

Make `inst/fit_rpbnb_diff_famoye.R` easy to read and edit while using every
detected physical CPU core by default.

## Settings

The script exposes three ordinary variables near the top:

```r
n_obs <- 5000L
draws <- 400L
n_cores <- parallel::detectCores(logical = FALSE)
if (is.na(n_cores)) n_cores <- 1L
```

The available dataset has fewer than 5,000 rows, so the script uses
`min(n_obs, nrow(data))` when selecting observations. Users change workload
values by editing these three lines. The script does not read environment
variables.

## Fit and Output

The existing formulas, random coefficients, Famoye dependence, seed, model
summary, fitted-mean section, dependence report, marginal effects, and
elasticities remain unchanged.

The model call uses `draws` and `n_cores` directly. Before fitting, the script
prints the selected observation count, draw count, and requested physical-core
count. After fitting, it prints elapsed time, requested and realized TMB
threads, optimizer convergence, and Hessian status.

The four Famoye-specific configuration, diagnostics, and allocation-guidance
helper functions are removed. The equivalent diagnostics are printed directly
after fitting.

## Memory Trade-off

The requested `5,000` observations, `400` draws, and all physical cores can
exhaust memory while TMB constructs its AD tape. The script intentionally
keeps these user-selected settings and does not silently reduce them. Users
must lower `n_obs`, `draws`, or `n_cores` if `std::bad_alloc` occurs.

## Verification

Static tests verify that the script has the three editable settings, detects
physical rather than logical cores, caps observations to available rows, uses
the variables in the model call, and no longer contains the removed helper
functions. A small temporary copy of the script is run with reduced settings
to verify real fitting without changing the requested defaults.
