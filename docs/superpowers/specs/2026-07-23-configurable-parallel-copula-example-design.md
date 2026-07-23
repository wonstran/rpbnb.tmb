# Configurable Parallel Copula Example Design

## Goal

Update `inst/fit_rpbnb_diff_copula.R` so it remains an easy sourceable
example while allowing users to configure the sample size, simulation draws,
and TMB thread count without editing the file.

## Interface

The script reads three optional environment variables:

- `RPBNB_N_OBS`, default `500`
- `RPBNB_DRAWS`, default `20`
- `RPBNB_N_CORES`, default to the smaller of four and the detected number of
  physical cores

Each value must be a single positive integer. Missing or empty variables use
their defaults. Invalid values stop immediately with an error that names the
variable.

The script remains runnable from the package root with:

```r
source("inst/fit_rpbnb_diff_copula.R")
```

Users can override settings before sourcing:

```r
Sys.setenv(
  RPBNB_N_OBS = 500,
  RPBNB_DRAWS = 20,
  RPBNB_N_CORES = 4
)
source("inst/fit_rpbnb_diff_copula.R")
```

## Runtime Behavior

The existing Frank-copula model, formulas, random coefficients, seed, summary,
fitted means, dependence report, marginal effects, and elasticities remain
unchanged.

Before fitting, the script prints the selected observation count, draw count,
and requested core count. The observation count is capped at the number of
available data rows, with an informative message when capping occurs.

After fitting, the script prints:

- elapsed fitting time;
- requested and realized TMB thread counts from `fit$parallel`;
- optimizer convergence code and message;
- whether `sdreport` reports a positive-definite Hessian.

## Error Handling

Configuration parsing fails before model fitting when an environment variable
is non-numeric, non-integral, zero, negative, non-finite, or contains multiple
values. Dataset and package-loading behavior otherwise remain unchanged.

## Verification

Static tests source the configuration helper portion in an isolated
environment and verify defaults, overrides, invalid-value errors, and
observation capping without running the expensive model. Existing installed
package copula and parallel regression tests continue to verify fitting and
OpenMP behavior.
