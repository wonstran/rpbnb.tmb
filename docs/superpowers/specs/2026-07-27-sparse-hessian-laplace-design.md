# Sparse-Hessian Laplace Estimation Path Design

## Goal

Add a second estimator to `fit_rpbnb_tmb()` that integrates the random coefficients with TMB's
Laplace approximation instead of Halton simulation, so that fits whose tape size is currently
bounded by `n x draws` become bounded by `n` alone.

The concrete target is `inst/truck_rpbnb_diff_famoye_dense.R`: 12,083 observations, four random
slopes per equation, Famoye dependence, 300 draws. That configuration exhausts memory during
`MakeADFun()` on a 31 GB machine. Under the new path it must fit.

## Background: why the current path runs out of memory

The simulated maximum likelihood (SML) objective tapes every one of the `n x R` conditional
density evaluations. `src/rpbnb.tmb.cpp:357` loops over observations, `:360` over draws, and the
full Famoye or copula density is recorded inside both. Measured tape cost is roughly 853 bytes per
observation-draw for Famoye, so the truck workload records about 12,083 x 300 x 853 = 3.1 GB of
tape, and peak working set during construction runs several times that.

No amount of thread or workload tuning changes the scaling. The memory investigation in
`comments/review_2026-07-24-1928.md` established that `parallel_tape = FALSE` removes a thread
multiplier on the peak, but the `n x R` term itself is intrinsic to SML.

The Laplace approximation removes `R` from the tape. TMB records the *joint* negative log
likelihood over data and latent variables once per observation, then integrates the latents
numerically using the sparse Hessian of that joint objective. Tape size becomes O(n).

## Selected approach

One template, one DLL, with a `method` switch. The SML and Laplace branches share every density
function.

This is possible because the existing per-draw block already *is* the conditional likelihood. At
`src/rpbnb.tmb.cpp:311`, `compute_dev(b, s, base, dist_code, sign_code)` maps a standard-normal
`base` to a coefficient deviation. In SML the `base` is a Halton-derived constant
(`u_to_base()` applies `qnorm` to a uniform draw); in Laplace the `base` is a latent parameter that
is already standard normal. Nothing else about the deviation, the linear predictor, the negative
binomial kernels, the Famoye factor, the copula terms, or the numerical clamps differs between the
two estimators.

Two alternatives were considered and rejected. A separate template sharing extracted headers would
require refactoring validated numerical C++ before adding any feature, and would double the compile
time and DLL management for no structural gain. Expressing the integral through TMB's newer
`intern` / `integrate` facilities would avoid template changes but those features are thinly
documented, interact poorly with `parallel_accumulator`, and would not reuse the existing family
code paths cleanly.

## Scope

In scope: normal and lognormal random-coefficient distributions under Laplace. All four dependence
families. All existing inference, prediction, and marginal-effect methods.

Out of scope: uniform and triangular distributions under Laplace. Their bounded support and kinked
densities are not twice differentiable in a way the inner Newton solver can use, and smoothing them
would silently change the assumed distribution. These remain SML-only and are rejected with a clear
error.

Also out of scope: exposing TMB `inner.control` tuning knobs, changing the SML estimator in any
way, and adaptive quadrature refinements beyond the plain Laplace approximation.

## C++ architecture

### New inputs

```cpp
DATA_INTEGER(est_method);    // 0 = SML, 1 = Laplace
PARAMETER_MATRIX(u1);        // n x q1 latent standard normals
PARAMETER_MATRIX(u2);        // n x q2 latent standard normals
```

The data field is named `est_method` rather than `method` to keep it visibly distinct from
`MakeADFun()`'s own `method` argument, which selects the outer optimizer and is unrelated.

`u1` and `u2` carry one row per observation. Each observation is an independent unit with its own
random-coefficient vector; the SML draws are a shared quadrature grid approximating that
observation-level integral, so an observation-level latent is the faithful Laplace counterpart.

### Branch structure

The precomputation of `dev1` / `dev2` at `src/rpbnb.tmb.cpp:336` is guarded to run only when
`est_method == 0`, since it is indexed by draw and has no meaning under Laplace.

Inside the observation loop, the draw loop runs `R` times under SML and exactly once under Laplace.
Only two things change in the Laplace branch:

1. Deviations are read from the latents rather than the precomputed draw matrix:
   `compute_dev(beta1(col), sd1(j), u1(i, j), dist1(j), sign1(j))`. Note `u_to_base()` is *not*
   applied — the latent is already the standard-normal base that `u_to_base()` would have produced.
2. The reduction changes. SML keeps its dynamic log-sum-exp minus `logR`. Laplace adds the single
   conditional contribution plus the latent prior:

```cpp
nll -= log_draw(0);
for (int j = 0; j < q1; j++) nll -= dnorm(u1(i, j), Type(0), Type(1), true);
for (int j = 0; j < q2; j++) nll -= dnorm(u2(i, j), Type(0), Type(1), true);
```

Everything between — `mu1`, `mu2`, the family dispatch, the Sarmanov validity penalty at
`src/rpbnb.tmb.cpp:402`, the copula branches, the trailing `ADREPORT` calls — is untouched and
shared.

### Parallel accumulation

`parallel_accumulator<Type>` is retained for both branches. Each observation's contribution under
Laplace is self-contained (its conditional density plus its own latent priors), which is a valid
accumulation unit. TMB supports `parallel_accumulator` together with `random=`, applying
`reorder_random` to the tape.

This combination is a verification point, not an assumption. If the Laplace path produces wrong
gradients or fails to converge under multiple threads, the fallback is a plain `Type nll = 0`
accumulator in the Laplace branch; the memory goal does not depend on threading.

## R architecture

### Public interface

`fit_rpbnb_tmb()` gains `method = c("sml", "laplace")`, matched with `match.arg()`, defaulting to
`"sml"`. Existing calls are unaffected.

Validation under `method = "laplace"`:

- any `dist1` or `dist2` code outside `{0, 1}` (normal, lognormal) is an error naming the offending
  distribution and directing the user to `method = "sml"`;
- `total_rand == 0` is an error, since there is nothing to integrate;
- a `draws` value supplied by the user warns that it is ignored. Detection uses `missing(draws)`,
  so the default never warns.

### Object construction

The parameter list at `R/fit_rpbnb_tmb.R:250` gains `u1 = matrix(0, n, q1)` and
`u2 = matrix(0, n, q2)`.

Under SML both are mapped off:

```r
map[["u1"]] <- factor(rep(NA_integer_, n * q1))
map[["u2"]] <- factor(rep(NA_integer_, n * q2))
```

A zero-length factor is the correct value when `q1` or `q2` is zero. Mapped parameters become tape
constants, and the SML branch never reads them, so the SML tape is bit-identical to today's.

Under Laplace, `.make_rpbnb_tmb_object()` is called with `random = c("u1", "u2")`, dropping any name
whose dimension is zero. The `random` argument already exists at `R/tmb_helpers.R:45` and is
currently always `NULL`; this is the first caller to use it.

`Z1` and `Z2` become 1 x q dummy matrices under Laplace and no Halton sequence is generated, so
`control$halton_burn` and `seed` have no effect on that path.

### Guardrails

`.check_tmb_workload()` is called with `draws = 1L` under Laplace. The tape genuinely is O(n), so
the guard measures the right quantity rather than being disabled. The truck workload evaluates to
12,083 x 1 x 1 = 1.2e4 against a default `max_workload` of 2e6, which passes without the caller
needing to override anything.

`.configure_tmb_threads()` is unchanged, including the `parallel_tape = FALSE` default that keeps
tape construction sequential.

Expected truck-fit footprint: roughly 10 MB of likelihood tape, an ADGrad tape of similar order, and
a sparse Hessian over 12,083 x 8 = 96,664 latents. That Hessian is block diagonal with 12,083
independent 8 x 8 blocks, because no latent is shared across observations, so its Cholesky factor
carries no fill beyond the blocks themselves.

### Inference and methods

Unchanged. With `random=` supplied, `obj$fn` and `obj$gr` are the Laplace-approximated marginal
likelihood and its gradient over the same fixed-parameter vector the optimizer already handles. The
`stats::optimHess()` call at `R/inference.R:506` operates on that marginal objective and remains
correct. `predict()`, `rpbnb_tmb_marginal_effects()`, and `rpbnb_tmb_elasticities()` read fitted
coefficients and need no change.

The fit object records `method` so `print()` and `summary()` can label which estimator produced the
result. Reporting the estimator is required — two fits of the same model under different estimators
must not be indistinguishable in their printed output.

## Error handling

Distribution, zero-random-coefficient, and ignored-`draws` conditions are checked in R before any
tape is built, so they cost nothing and produce R-level messages rather than C++ exceptions.

The inner Newton solver may fail to converge on badly scaled data. TMB signals this through
non-finite `obj$fn` values, which the existing optimizer path already surfaces via
`fit$optimizer$convergence` and the `pdHess` flag on the sdreport. No new machinery is added; the
Laplace path reuses the existing convergence reporting.

The Sarmanov validity penalty at `src/rpbnb.tmb.cpp:402` is a value-only barrier with zero gradient,
as documented in that comment block. Under Laplace it sits inside the inner optimization as well as
the outer, so a fit that lands on the penalty is invalid in the same way and for the same reason it
is under SML. This is pre-existing behavior, not a regression, and is out of scope here.

## Testing

**SML regression.** An existing fixture is fit before and after the change and must return an
identical log likelihood. This is the direct evidence that the mapped `u1` / `u2` parameters never
enter the SML tape. The full existing suite must pass unmodified.

**Laplace mechanism.** On a small synthetic fit, assert that `length(obj$env$random)` equals
`n * (q1 + q2)`, confirming the sparse path is actually engaged rather than silently falling back.

**Laplace correctness.** Simulate from a known RP-BNB process with n around 500 and one or two
normal random slopes, then fit the same data under both estimators. Laplace and SML coefficient
estimates must agree within a few standard errors. This is the test that decides whether the feature
is usable, and it is written before the truck script is attempted.

**Error paths.** Uniform or triangular distribution with `method = "laplace"` errors; no random
coefficients with `method = "laplace"` errors; explicitly supplied `draws` with
`method = "laplace"` warns.

**Acceptance.** `inst/truck_rpbnb_diff_famoye_laplace.R`, a clone of the dense script with
`method = "laplace"`, completes on this machine and reports a positive-definite Hessian. This is the
fit that currently fails.

## Risks

**Laplace accuracy with few observations per latent.** Each unit contributes two counts but carries
up to eight latents, so the conditional posterior of the latents is dominated by the standard normal
prior. That makes the inner problem well conditioned and close to Gaussian, which favors the
approximation, but the likelihood curvature contribution is still being approximated and the
resulting bias is not known in advance for this model. The Laplace-versus-SML agreement test above
is the gate. If the two estimators disagree materially on data simulated from the true process, that
is a finding to report and act on, not a result to ship quietly.

**`parallel_accumulator` with `random=`.** Supported by TMB but unverified for this template. Fallback
is a serial accumulator in the Laplace branch, as described above.

**Windows compile size.** The template already needs `-Wa,-mbig-obj` in `src/Makevars.win`. Adding a
second branch grows the object file further. If the build breaks, the mitigation is factoring the
family dispatch into a function template rather than reverting the design.

## Documentation

Roxygen for the new `method` argument states plainly that SML and Laplace are different
approximations to the same integral, that they agree asymptotically but need not agree closely on a
given dataset, and that Laplace removes the `draws` dimension from the memory cost. The README memory
section names Laplace as the answer for large `n` rather than only describing thread and workload
tuning.
