# Design: `rpbnb_tmb_dependence_profile()`

Date: 2026-07-25

## Problem

`summary()` reports `NA` for the standard error of a dependence parameter whose
link derivative has collapsed. For Famoye this is common, because the admissible
lambda interval is computed once at the starting values and held fixed for the
whole fit (`fit_rpbnb_tmb.R`, "Famoye: compute frozen lambda bounds at start"),
so lambda is mapped through a logistic into a box that the optimizer cannot
leave.

The suppression is correct. `SE(lam) = |dlam/dz| * SE(z_dep)`
(`inference.R`, `.rpbnb_inference`) and
`dlam/dz = span * (1 - 2*eps) * sigma(z) * (1 - sigma(z))`
(`inference.R`, `.rpbnb_natural_report`). At a pinned optimum `sigma(z)`
approaches 0 or 1, so the derivative approaches 0 while `SE(z_dep)` is typically
diverging. The product is `0 * Inf`: either absurdly small, claiming precision
that does not exist, or numerically unstable. `report_sd[boundary_report] <- NA_real_`
is applied after the number is computed, deliberately.

There is therefore no valid symmetric standard error to recover on the natural
scale. What users need instead is an *interval*, obtained without dividing by a
vanishing derivative.

## Goal

One exported function returning a confidence interval for the dependence
parameter of any fitted `rpbnb_tmb_fit`, valid at or near a boundary.

Out of scope: fixing lambda at a supplied value (a `par` slot for Famoye
analogous to `copula(family, par =)`). Recorded as possible later work; the
constrained-refit machinery it needs is not built here.

## Key facts established by code reading

These constrain the design and were each verified against the installed sources
rather than assumed.

1. `TMB::tmbprofile()` profiles from `obj$env$last.par.best`, **not** `obj$par`.
2. `fit$optimizer` (the whole `nlminb` return) is retained in every `keep` mode.
   `nlminb` is called with `start = obj$par`, so `fit$optimizer$par` holds the
   free parameters at the optimum, carrying TMB-style names in `obj`'s order and
   directly assignable to `obj$env$last.par.best`.
3. `fit$obj` is retained only under `keep = "full"`. Rebuilding it from a
   lesser fit is impossible: it needs `Y1`/`Y2`, which are themselves retained
   only under `keep = "full"`.
4. `z_dep` is a scalar `PARAMETER`, so `tmbprofile()`'s
   `sum(names(par) == name) != 1` uniqueness check passes.
5. `TMB:::confint.tmbprofile(object, parm, level, ...)` accepts `level` and
   locates endpoints with `approx()`. An uncrossed profile therefore yields
   **`NA`**, not `+/-Inf`.
6. Every dependence link is monotone increasing in `z`, so mapping interval
   endpoints through it is valid for every family:
   - Famoye: `lamLo + span * (eps + (1 - 2*eps) * plogis(z))`, `eps = 1e-6`
   - Frank: `FRANK_THETA_MAX * tanh(z / FRANK_THETA_MAX)`
   - Gaussian: `tanh(z)`
   - Clayton: `exp(pmin(pmax(z, -20), 20))`
   Each family's Kendall tau is monotone in its own parameter as well
   (`.frank_tau`; `2/pi * asin(rho)`; `theta / (theta + 2)`), so tau rows are
   equally valid.
7. `fit$se` is named by `par_names` and retained in every `keep` mode and both
   non-`none` `inference` modes. It is a better source for `SE(z_dep)` than
   `vcov(fit)`, which is `NULL` under `inference = "diag"` and reconstructed by
   the `vcov` method.
8. Independence maps `z_dep` out via `map[["z_dep"]] <- factor(NA)`, so it has
   no profile and no standard error.

## API

```r
rpbnb_tmb_dependence_profile(fit,
                             level  = 0.95,
                             method = c("profile", "wald"),
                             ...)
```

- `fit` — an `rpbnb_tmb_fit`.
- `level` — one number strictly in `(0, 1)`.
- `method` — `"profile"` (default) or `"wald"`.
- `...` — passed to `TMB::tmbprofile()` (e.g. `ytol`, `parm.range`).

### Return value

A data frame, one row per reported dependence quantity, with columns:

| column | meaning |
|---|---|
| `parameter` | `"lam"`, or `"theta"`/`"rho"` plus `"tau"` |
| `estimate` | natural-scale point estimate, from the shared link at `z_hat` |
| `lower`, `upper` | interval endpoints on the natural scale |
| `level` | the `level` argument, recycled |
| `method` | `"profile"` or `"wald"` — the method actually used, not the one requested |

`attr(result, "profile")` holds the raw `tmbprofile` object when the profile
path ran, and is `NULL` otherwise, so `plot()` on the profile stays available.

Rows are ordered with the family's own parameter first, then `tau`.

## Architecture

Three internal pieces plus the exported function.

### `.resolve_family_code(dependence)` — new, in `utilities.R`

Maps a `dependence` value to a family code, returning `-1` for
`"independence"`, `0` for `"famoye"`, and `1`/`2`/`3` for the `rpbnb_copula`
families. `fit_rpbnb_tmb()` currently inlines this `switch`; it is changed to
call the helper so the profile function cannot drift from it. Invalid input
raises the same error `fit_rpbnb_tmb()` raises today.

### `.rpbnb_dependence_link(family_code, lamLo, lamHi)` — new, in `inference.R`

Returns a list with:

- `names` — character vector of reported quantity names for the family, in
  output order (`"lam"`; or `"theta"`/`"rho"` then `"tau"`).
- `map` — function of `z` returning a named numeric vector of the natural-scale
  values, matching `names`.

`.rpbnb_natural_report()` is refactored to obtain its `value` entries from this
`map`, keeping its own derivative expressions and boundary logic unchanged. This
is the anti-drift measure: there is exactly one definition of each link, and a
test asserts the profile function and `summary()` agree at `z_hat`.

Returns `NULL` for `family_code < 0`.

### `.dependence_profile_ci(fit, level, ...)` — new, in `inference.R`

The profile path. Returns a length-2 unnamed numeric of **working-scale**
endpoints for `z_dep` (either may be `NA`), carrying the raw profile in
`attr(, "profile")`; or `NULL` to signal "fall back to Wald". `NULL` and
`c(NA, NA)` are distinct: `NULL` means the profile could not be attempted and
Wald should run, `c(NA, NA)` means it ran and did not bracket.

1. Return `NULL` if `fit$obj` is `NULL`.
2. Probe liveness: `try(fit$obj$fn(fit$optimizer$par), silent = TRUE)`, and
   require a finite scalar result. Return `NULL` on error or non-finite result.
   (See the corrected note below on what this does and does not catch — a
   `saveRDS()` round-trip is **not** one of the cases it rejects.)
3. Restore the optimum: `fit$obj$env$last.par.best <- fit$optimizer$par`.
4. `pr <- TMB::tmbprofile(fit$obj, "z_dep", trace = FALSE, ...)`.
5. `ci <- confint(pr, level = level)`.
6. If either endpoint is `NA` **and** `"ytol"` is not among `names(list(...))`,
   retry once with `ytol = 10` (against `tmbprofile`'s default of `2`). Still
   `NA` — keep the `NA`.
7. Attach `pr` as an attribute and return.

**Corrected during implementation: reloaded fits still profile.** This section
originally claimed that a fit round-tripped through `saveRDS()` loses its
external pointer and must therefore degrade to a Wald interval. That was
verified false on TMB 1.9.21:

- `reloaded$obj$env$ADFun` prints `$ptr <pointer: (nil)>` — the pointer really
  is dead.
- `reloaded$obj$fn(reloaded$optimizer$par)` nevertheless returns a value
  bit-identical to the pre-save one (`366.8121741` in the diagnostic run).

TMB detects the nil pointer and retapes from the data and parameters it keeps in
`obj$env`, so the objective is fully usable and the resulting profile is
correct. The bit-identical value rules out a dangling-pointer read of freed
memory. Retaping requires the package's compiled DLL to be loaded in the
session, which always holds when the package is in use.

The probe in step 2 is therefore still correct and worth keeping — it rejects an
objective that genuinely cannot be evaluated — but a `saveRDS()` round-trip is
not such a case, and forcing Wald there would discard a valid profile interval
for no benefit. Test 5 asserts the profile path succeeds after a round-trip; the
Wald fallback is covered instead by the `keep = "postfit"` test, where `fit$obj`
is genuinely `NULL`.

**Residual caution.** The retaping behaviour above is observed on TMB 1.9.21,
not a contract TMB documents as stable. If a future TMB stops retaping, test 5
fails loudly rather than silently returning a wrong interval, which is the
failure mode to prefer. The separate assumption that TMB raises an R-level error
rather than crashing on an unusable pointer remains unenforceable by this design,
since there is no portable way to inspect pointer validity without calling into
it.

### `rpbnb_tmb_dependence_profile()` — new, in `inference.R`

Orchestrates:

1. Validate `fit` class and `level`.
2. `family_code <- .resolve_family_code(fit$dependence)`; error if `< 0`.
3. If `method = "profile"`, call `.dependence_profile_ci()`. On `NULL`, warn and
   set `method <- "wald"`.
4. Wald path: `z_hat +/- qnorm(1 - (1 - level)/2) * fit$se[["z_dep"]]`.
5. Map `z_hat` and both endpoints through `link$map`, assemble the data frame.
6. Warn, naming the side, if any endpoint is `NA`.

## Error handling

| condition | behaviour |
|---|---|
| `fit` not an `rpbnb_tmb_fit` | error |
| `level` not one number in `(0, 1)` | error |
| independence dependence | error: no dependence parameter is estimated |
| `fit$obj` absent or pointer stale | warning, degrade to Wald |
| `inference = "none"` on the Wald path | warning, `NA` endpoints, `estimate` still returned |
| profile endpoint uncrossed after retry | warning naming the side, `NA` for that endpoint |
| Clayton `z_hat` at the `exp` clamp | no special case; the link is the identity inside `(-20, 20)` and flat outside, so a mapped endpoint is correct either way |

Degrading rather than erroring on a missing `obj` follows the precedent in
`rpbnb_tmb_marginal_effects()`, which warns about a missing full covariance and
still returns point estimates.

## Testing

New file `tests/testthat/test-dependence-profile.R`. Fits are small
(`n` around 300, `draws` 50) to keep the suite fast; the existing suite's peak
memory contract is unaffected because no fit here approaches `max_workload`.

1. **Brackets the estimate.** Famoye fit: `lower < estimate < upper`.
2. **Respects the frozen box.** Both endpoints lie within `fit$lambda_bounds`.
3. **Nesting in `level`.** The `0.99` interval contains the `0.95` interval.
4. **Wald fallback warns and works.** A `keep = "postfit"` fit warns and returns
   a finite interval with `method == "wald"`.
5. **Reloaded fits still profile.** A fit round-tripped through
   `saveRDS()`/`readRDS()` returns `method == "profile"` with no warning,
   because TMB retapes from the data it keeps in `obj$env`. (Originally
   specified as a Wald degradation; corrected during implementation — see
   "Corrected during implementation" above.)
6. **Independence errors.**
7. **Frank returns both rows.** `parameter` is exactly `c("theta", "tau")`, in
   that order.
8. **Anti-drift guard.** `estimate` equals the value `summary(fit)` reports for
   the same quantity, for Famoye and for one copula family. This is the test
   that makes the shared-link refactor worth doing.
9. **`inference = "none"`.** Warns, `NA` endpoints, non-`NA` estimate.

## Documentation

- Roxygen block with `@export`, and `@examples` in `\dontrun{}` matching the
  style of `fit_rpbnb_tmb()`.
- Cross-reference from the `lambda_bounds` documentation in `fit_rpbnb_tmb()`
  and from the boundary warning's remedy text, so a user who hits the warning is
  pointed at the function that answers it.
- `NAMESPACE` and `man/` regenerated with roxygen.
