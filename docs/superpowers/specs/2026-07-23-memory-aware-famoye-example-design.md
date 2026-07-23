# Memory-Aware Famoye Example Design

## Goal

Update `inst/fit_rpbnb_diff_famoye.R` so users can control its workload
without editing the file and receive useful diagnostics when TMB cannot
allocate the automatic-differentiation tape.

## Configuration

The script reads three optional, script-specific environment variables:

- `RPBNB_FAMOYE_N_OBS`, default `500`
- `RPBNB_FAMOYE_DRAWS`, default `100`
- `RPBNB_FAMOYE_N_CORES`, default to the smaller of two and the detected
  number of physical cores

Script-specific names prevent configuration intended for the copula example
from changing the Famoye workload. Each value must be one finite, positive
whole number no larger than `.Machine$integer.max`. Missing or empty values
use their defaults.

The requested observation count is capped at the number of available rows,
with a message that reports the requested and available counts.

## Runtime Behavior

The existing formulas, random coefficients, Famoye dependence, seed, summary,
fitted-mean section, dependence report, marginal effects, and elasticities
remain unchanged.

Before fitting, the script prints the selected observation count, draw count,
and requested core count. After fitting, it prints:

- elapsed fitting time;
- requested and realized TMB thread counts;
- optimizer convergence code and message;
- whether `sdreport` reports a positive-definite Hessian.

The script remains runnable from the package root:

```r
source("inst/fit_rpbnb_diff_famoye.R")
```

Users can configure it before sourcing:

```r
Sys.setenv(
  RPBNB_FAMOYE_N_OBS = 500,
  RPBNB_FAMOYE_DRAWS = 100,
  RPBNB_FAMOYE_N_CORES = 2
)
source("inst/fit_rpbnb_diff_famoye.R")
```

## Allocation Failure Guidance

The model-fit call is wrapped with error handling. Errors that contain
`std::bad_alloc` are rethrown with a concise explanation that TMB exhausted
memory while constructing its AD tape and advice to restart R or reduce
`RPBNB_FAMOYE_N_OBS`, `RPBNB_FAMOYE_DRAWS`, or
`RPBNB_FAMOYE_N_CORES`. Other errors are rethrown without changing their
message or class.

## Verification

Tests parse and evaluate only the example helper definitions in an isolated
environment. They verify defaults and overrides, invalid values, observation
capping, fit diagnostics, enhanced `std::bad_alloc` guidance, and unchanged
handling of unrelated errors. Installed-package Famoye and parallel tests
remain green, and a small configured Famoye fit verifies the complete script.
