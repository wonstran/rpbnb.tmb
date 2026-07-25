# Address Review 2026-07-24 19:28 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Resolve the confirmed peak-memory, correctness, API, test-harness, and packaging defects in `comments/review_2026-07-24-1928.md`, and record the disposition of every review item.

**Architecture:** Keep simulated likelihood semantics stable while making TMB tape construction sequential by default, repairing mapped-parameter reconstruction, and hardening dependence/numerical boundaries. Replace source-text tests with installed-package behavioral tests, make predictions and simulations honest, and use `R CMD check` as the final integration gate. Extreme-tail Gaussian quadrature redesign and broad demo/data restructuring are documented as follow-up work where a safe atomic derivative redesign is larger than this patch.

**Tech Stack:** R 4.5, TMB 1.9.21, C++17, OpenMP, testthat 3, roxygen2.

## Global Constraints

- Work only in `.worktrees/codex-address-review-20260724` on `codex/address-review-20260724`.
- Preserve the user's uncommitted `inst/fit_rpbnb_diff_famoye.R`, `comments/`, and `dist/` state in the main checkout.
- Write a behavioral regression before every production behavior change and observe the expected failure.
- Keep `optimHess()` rather than `obj$he()` because the exact Hessian tape has a measured multi-gigabyte peak.
- Do not change ordinary-fit Gaussian copula values or derivatives while addressing the extreme-tail review separately.

---

### Task 1: Sequential TMB taping and calibrated workload policy

**Files:**
- Modify: `R/utilities.R`
- Modify: `R/tmb_helpers.R`
- Modify: `R/fit_rpbnb_tmb.R`
- Modify: `README.md`
- Test: `tests/testthat/test-parallel.R`

**Interfaces:**
- Produces: `rpbnb_tmb_control(..., parallel_tape = FALSE)`
- Produces: `.configure_tmb_threads(..., parallel_tape, DLL)` returning the realized evaluation thread count.
- Produces: `.check_tmb_workload(..., n_threads, parallel_tape)` using measured family weights and a thread multiplier only when tapes are constructed concurrently.

- [ ] Add behavioral tests that query `TMB::config(DLL = "rpbnb.tmb")` after configuration, verify `tape.parallel == 0` by default and `1` on explicit opt-in, and verify Frank/Gaussian workload weights at their measured ratios.
- [ ] Run the focused parallel test and confirm failures for the missing `parallel_tape` control and old weight of four.
- [ ] Configure `TMB::config(tape.parallel = as.integer(parallel_tape), DLL = DLL)`, resolve threads before the workload check, and use weights 1.0 for independence/Famoye/Frank and 1.25 for Gaussian/Clayton.
- [ ] Re-run the focused test and confirm serial/parallel objectives, gradients, and Hessians remain equal.

### Task 2: Mapped parameters, model degrees of freedom, RNG, and information criteria

**Files:**
- Modify: `R/fit_rpbnb_tmb.R`
- Modify: `R/methods.R`
- Modify: `NAMESPACE`
- Test: `tests/testthat/test-api-correctness.R`

**Interfaces:**
- `fit$coef` retains named pinned parameters while `fit$npar` counts only free parameters.
- `logLik(fit)` exposes the correct free-parameter `df`.
- Base `AIC.default()` and `BIC.default()` consume `logLik()`; package-specific AIC/BIC methods are removed.

- [ ] Add failing tests for both Poisson-pinned margins, correct `df`, `AIC(fit, k = 0)`, and preservation of an existing and absent `.Random.seed`.
- [ ] Run the new test file and confirm the mapped-position crash, inflated `df`, ignored `k`, and RNG mutation.
- [ ] Reconstruct coefficients as `coef_vec <- start; coef_vec[free] <- opt$par`, count `sum(free)`, scope Halton seeding with save/restore logic, and remove the AIC/BIC S3 methods.
- [ ] Re-run the new tests and existing inference tests.

### Task 3: Dependence simulation correctness

**Files:**
- Modify: `R/simulate_rpbnb_tmb.R`
- Test: `tests/testthat/test-simulation.R`

**Interfaces:**
- Copula simulation dispatch depends only on the `rpbnb_copula` object.
- Copula simulation requires a finite family-valid `par`.
- Gaussian simulation gives an explicit `pbivnorm` installation error when unavailable.

- [ ] Add failing tests showing non-zero association from Clayton simulation with default `lambda`, clear rejection of unknown dependence, and validation of missing copula parameters.
- [ ] Run the tests and confirm the default-lambda independence bug.
- [ ] Dispatch independence, Famoye, and copulas exclusively by `dependence`; add a final validation error; guard `pbivnorm`; include copula truth in the returned object.
- [ ] Re-run simulation and fitting regressions.

### Task 4: Famoye-boundary and C++ numerical hardening

**Files:**
- Modify: `R/inference.R`
- Modify: `R/fit_rpbnb_tmb.R`
- Modify: `src/rpbnb.tmb.cpp`
- Test: `tests/testthat/test-numerical-guards.R`
- Test: `tests/testthat/test-inference-memory.R`

**Interfaces:**
- Famoye saturation returns the estimate with `NA` uncertainty and emits a fit warning.
- Frank uses the smooth bounded transform `theta = 35 * tanh(z_dep / 35)` and a removable-singularity series near zero.
- Negative-binomial means use an objective-safe lower log-mean bound of `-35`.

- [ ] Add failing tests for saturated Famoye reports, finite Frank objectives at extreme working parameters, finite NB objectives at intercept `-40`, and finite Frank gradients at independence.
- [ ] Run the focused tests and confirm each existing failure.
- [ ] Mark saturated Famoye delta-method uncertainty unavailable, warn from the fit, mirror the bounded Frank transform in R/C++, use `expm1`/`log1p` plus a local independence series, tighten the log-mean floor, and replace unsafe Debye zero handling.
- [ ] Rebuild the DLL and verify objective/gradient/Hessian checks for all copula families.

### Task 5: Prediction and diagnostic API behavior

**Files:**
- Modify: `R/utilities.R`
- Modify: `R/fit_rpbnb_tmb.R`
- Modify: `R/methods.R`
- Modify: `R/marginal_effects.R`
- Modify: `inst/fit_rpbnb_diff_famoye.R`
- Modify: `inst/fit_rpbnb_diff_copula.R`
- Test: `tests/testthat/test-predict-diagnostics.R`

**Interfaces:**
- `predict.rpbnb_tmb_fit(object, newdata = NULL, type = c("response", "link"))` returns an `n × 2` data frame of draw-integrated marginal means or their logs.
- Marginal-effect exports return invisibly after printing once.
- Diagonal/no-inference fits warn that marginal-effect covariance is unavailable.

- [ ] Add failing tests for fitted/new-data prediction dimensions and values, one-print diagnostic behavior, and the diagonal-inference warning.
- [ ] Run the tests and confirm the current coefficient-vector stub, duplicate auto-print, and silent `NA` uncertainty.
- [ ] Retain formula terms/contrasts with post-fit state, implement draw-integrated prediction using `.draw_mean_exp()`, return exported diagnostics invisibly, warn on missing covariance, and use the distribution registry for scale labels.
- [ ] Update demos to call `predict(fit)` and re-run the diagnostics tests.

### Task 6: Installed-package behavioral tests and check-safe fixtures

**Files:**
- Delete: `tests/testthat/test-devtools-load-all.R`
- Delete: `tests/testthat/test-example-copula.R`
- Delete: `tests/testthat/test-example-famoye.R`
- Modify: `tests/testthat/test-fit-copula.R`
- Modify: `tests/testthat/test-inference-memory.R`
- Modify: `tests/testthat/test-parallel.R`
- Modify: `tests/testthat/test-against-rpbnb.R`
- Modify: `inst/benchmark_parallel.R`

**Interfaces:**
- Tests run against an installed package without reading unavailable `R/` or `src/` sources.
- Numerical invariants replace frozen source strings and golden gradients.

- [ ] Run `R CMD check --no-manual` on the pre-change archive and retain its expected source-path failures as RED evidence.
- [ ] Remove runtime `devtools::load_all()` and source-grep assertions; locate installed scripts with `system.file()` only where script parsing is useful.
- [ ] Cross-check `obj$gr` against `numDeriv::grad` for Frank, Gaussian, and Clayton; assert Hessian symmetry and serial/parallel equality without frozen digits.
- [ ] Remove `library(rpbnb)` cross-test contamination and skip only the external malformed fixture tests.
- [ ] Run installed-package tests and confirm no source-path failures.

### Task 7: Package metadata, documentation, and build hygiene

**Files:**
- Modify: `DESCRIPTION`
- Modify: `R/rpbnb.tmb-package.R`
- Modify: `R/utilities.R`
- Modify: `NAMESPACE`
- Modify: `src/Makevars`
- Modify: `src/Makevars.win`
- Modify: `.Rbuildignore`
- Modify: `LICENSE`
- Create: `LICENSE.md`
- Generate: `man/*.Rd`

**Interfaces:**
- Roxygen owns `NAMESPACE`.
- `copula()` has public documentation.
- `RcppEigen` is declared in `LinkingTo` and Makevars contains no shell probe.

- [ ] Add documentation/check expectations by running `devtools::document()` and `R CMD check --no-manual` to capture current warnings.
- [ ] Document `copula()`, declare all used stats imports, add `RcppEigen` to `LinkingTo`, remove redundant author fields and Makevars shell probes, remove superseded dead helpers, and expand `.Rbuildignore`.
- [ ] Convert `LICENSE` to the required DCF stub and retain the full MIT text in `LICENSE.md`.
- [ ] Regenerate `NAMESPACE`/Rd files and run `R CMD build` plus `R CMD check --no-manual`.

### Task 8: Review response and final verification

**Files:**
- Create: `comments/response_2026-07-24-2236.md`
- Modify: `docs/superpowers/plans/2026-07-24-address-review-1928.md`

**Interfaces:**
- Response maps H1–H6, M1–M10, C1–C8, O1–O4, low-severity items, and coverage gaps to fixed, deferred, rejected, or already-correct status with validation evidence.

- [ ] Write the response with exact files changed, behavioral evidence, and explicit reasons for deferring adaptive extreme-tail BVN quadrature, bulk fixture relocation, and demo deduplication.
- [ ] Run a clean source install and all installed-package tests.
- [ ] Run `R CMD build` and `R CMD check --no-manual`; record ERROR/WARNING/NOTE counts.
- [ ] Run `git diff --check`, inspect status/diff scope, and ensure the main checkout's user-owned files are unchanged.
