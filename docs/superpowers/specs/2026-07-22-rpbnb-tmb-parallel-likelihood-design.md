# AD-Safe Parallel Likelihood Design

## Goal

Accelerate a single `fit_rpbnb_tmb()` call with `control$n_cores` while preserving the finite gradients, Hessians, coefficients, and standard errors of the corrected serial TMB implementation.

## Selected approach

Use TMB's `parallel_accumulator<Type>` to partition the log-likelihood by observation. Each observation is statistically independent conditional on the parameters and simulation draws, so its complete simulated log-likelihood contribution is a valid parallel accumulation unit.

Raw `#pragma omp parallel for` loops are not permitted inside the objective. They previously shared one AD tape and mutable likelihood storage across threads, producing NaN gradients. TMB automatic parallelization is not selected because it is more version-dependent and less explicit. Process-level R parallelism remains appropriate for multiple independent fits but is outside this feature's single-fit scope.

## C++ architecture

The objective retains the current parameter transforms, AD-aware `CppAD::CondExp*` clamps, copula formulas, and reported quantities.

Before parallel accumulation, it computes shared parameter-dependent values once:

- fixed linear predictors `xb1` and `xb2`;
- transformed dispersions and dependence parameters;
- random-coefficient deviations for every draw and random coefficient in `dev1` and `dev2` matrices;
- data-only Famoye values such as `exp(-Y1)` and `exp(-Y2)`.

The likelihood loop is reorganized from draw-first to observation-first:

1. Construct `parallel_accumulator<Type> nll(this)`.
2. For each observation `i`, allocate an `R`-element vector of draw-specific log-likelihoods.
3. For each draw `r`, construct scalar `eta1` and `eta2` from the fixed predictor and precomputed deviations.
4. Evaluate the selected marginal and dependence likelihood using the existing numerically stable formulas.
5. Reduce the draw likelihoods with the current dynamic log-sum-exp calculation.
6. Add the negative observation contribution with `nll -= log_contrib`.
7. Return `nll` after the existing `ADREPORT` calls.

This removes the `n x R` likelihood matrix. It replaces it with one `R`-element local vector per observation/tape, reducing shared mutable state and peak likelihood-storage memory.

## R thread configuration

`rpbnb_tmb_control()` continues to expose `n_cores`, with validation that it is one finite positive integer.

Immediately before `TMB::MakeADFun()`, `fit_rpbnb_tmb()` configures the loaded model DLL:

```r
TMB::openmp(n = control$n_cores, DLL = "rpbnb.tmb")
```

Every fit sets the requested value explicitly so a previous fit cannot leak its thread count into the next fit. `n_cores = 1L` selects the same parallel-accumulator code path with one tape and serves as the serial reference.

If the DLL was built without OpenMP and the user requests more than one core, the fit continues with one core and emits one clear warning. The fitted object records both `requested_n_cores` and `actual_n_cores` for diagnostics.

## Build configuration

Both `src/Makevars` and `src/Makevars.win` add R's portable OpenMP variables:

```make
PKG_CXXFLAGS = $(TMB_CXXFLAGS) $(SHLIB_OPENMP_CXXFLAGS) ...
PKG_LIBS = $(TMB_LIBS) $(SHLIB_OPENMP_CXXFLAGS)
```

`Makevars.win` retains `-Wa,-mbig-obj`, which is required for this nested-AD template on Windows. If `SHLIB_OPENMP_CXXFLAGS` is empty on a platform, the package builds without OpenMP and uses the documented single-core fallback.

## Numerical and behavioral requirements

- One-core and multi-core objective values agree within `1e-8` at identical parameters and draws.
- One-core and multi-core gradients agree within `1e-7` at identical parameters and draws.
- Fitted coefficients and standard errors agree within `1e-5` for a deterministic regression fixture.
- Frank, Gaussian, Clayton, Famoye, and independence likelihoods remain finite.
- The Frank random-parameter regression retains a positive-definite Hessian and finite standard errors.
- Poisson flags and all four random-coefficient distributions retain their existing semantics.
- Parallel reduction may change the last floating-point digits, but must not change convergence status or statistical conclusions.

## Testing

Add an internal helper that builds an objective without optimizing so tests can compare `fn` and `gr` at exactly the same parameter vector under one and multiple threads.

Tests cover:

1. `rpbnb_tmb_control()` rejects zero, negative, missing, non-finite, and non-scalar `n_cores` values.
2. A Frank random-parameter objective produces equivalent one-core and two-core values and gradients.
3. One-core and two-core Frank fits produce equivalent coefficients and finite standard errors.
4. Existing forced copula and Famoye test files continue to pass.
5. A clean Windows source build succeeds with OpenMP and the big-object flag.

A separate benchmark script compares elapsed time for one core and up to four cores on the 500-observation, 20-draw demonstration. Timing is reported for manual verification, not asserted in unit tests because shared CI and desktop load make timing thresholds unreliable.

## User-facing behavior

Existing calls remain valid. Users enable parallel fitting through the existing interface:

```r
control <- rpbnb_tmb_control(n_cores = 4L)
fit <- fit_rpbnb_tmb(..., control = control)
```

The fit object exposes the realized configuration in `fit$parallel`. Demo scripts stop using `parallel::detectCores()` without a cap; they request `min(4L, parallel::detectCores())` to avoid oversubscription.

## Acceptance criteria

The feature is accepted when a clean rebuilt package passes the parallel equivalence tests, the existing family regression tests, and the 500-observation demonstration with finite inference. A manual benchmark must show more than one active worker and report elapsed times, but a specific speedup is not required because the workload and hardware determine scaling.
