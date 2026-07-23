# Simple Copula Example Design

## Goal

Make `inst/fit_rpbnb_diff_copula.R` follow the simplified Famoye example's
structure while retaining a safer default workload for the more expensive
discrete-copula likelihood.

## Settings

The script exposes three ordinary variables near the top:

```r
n_obs <- 500L
draws <- 20L
n_cores <- parallel::detectCores(logical = FALSE)
if (is.na(n_cores)) n_cores <- 1L
```

Users change workload values by editing these lines. The script does not read
environment variables. It selects at most the requested number of rows with
`min(n_obs, nrow(data))`.

The lower observation and draw defaults are intentional. The discrete Frank
copula evaluates negative-binomial CDF corners for each simulation draw and is
substantially more expensive than the Famoye likelihood.

## Fit and Output

The existing formulas, random coefficients, Frank copula dependence, seed,
model summary, fitted-mean section, dependence report, and marginal effects
remain unchanged. The user's current commented-out elasticity section also
remains commented out.

The model call uses `draws` and `n_cores` directly. Before fitting, the script
prints the selected observation count, draw count, and requested core count.
After fitting, it directly prints:

- elapsed fitting time;
- requested and realized TMB thread counts;
- optimizer convergence code and message;
- whether `sdreport` reports a positive-definite Hessian.

## Removed Configuration Layer

The `.example_positive_integer`, `.example_observation_count`, and
`.example_fit_diagnostics` helpers are removed. The `RPBNB_N_OBS`,
`RPBNB_DRAWS`, and `RPBNB_N_CORES` environment-variable interface and its
documentation are also removed.

## Verification

Static tests verify that the script:

- contains the three editable settings and physical-core detection fallback;
- caps observations to available rows;
- uses `draws` and `n_cores` directly in the fit;
- retains fit diagnostics;
- no longer contains the removed helpers or environment-variable names.

A temporary copy runs with a small workload to verify a real copula fit without
changing the source defaults. Existing focused copula and parallel regression
tests are also run.

User changes in the copula elasticity section,
`inst/fit_rpbnb_diff_famoye.R`, `.gitignore`, and `.Rbuildignore` remain
outside this change.
