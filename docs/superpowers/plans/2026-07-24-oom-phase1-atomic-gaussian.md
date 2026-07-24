# OOM Phase 1 and Atomic Gaussian Kernel Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Reduce fit-time and post-fit memory while preserving the current simulated likelihood and compressing Gaussian copula quadrature into an atomic TMB operation.

**Architecture:** Add workload, inference, and retention policies at the R boundary; replace unconditional `sdreport()` with inference from the existing tape; stream post-fit draw reductions; and wrap the existing Gaussian quadrature with a known-derivative TMB atomic. Keep numerical and OpenMP compatibility under regression tests.

**Tech Stack:** R 4.5, TMB 1.9.21, C++17, OpenMP, testthat 3.

## Global Constraints

- Work only in `.worktrees/codex-oom-phase1` on `codex/oom-phase1`.
- Preserve the simulated-likelihood estimator and existing coefficient names.
- Do not call `TMB::sdreport()` from `fit_rpbnb_tmb()`.
- Preserve `summary(fit$sdreport, "fixed")` and `"report"` through a compact compatibility class.
- Treat the two malformed external-fixture failures in `test-against-rpbnb.R` as baseline failures.
- Do not modify the user's main-checkout Famoye example or untracked `dist/`.

---

### Task 1: Workload and thread policies

**Files:**
- Modify: `R/utilities.R`
- Modify: `R/tmb_helpers.R`
- Modify: `R/fit_rpbnb_tmb.R`
- Test: `tests/testthat/test-parallel.R`
- Test: `tests/testthat/test-fit-famoye.R`

**Interfaces:**
- Produces: `rpbnb_tmb_control(..., max_threads = 4L, max_workload = 2e6)`
- Produces: `.check_tmb_workload(n, draws, family_code, max_workload)`
- Produces: `.make_rpbnb_tmb_object(..., max_threads)`

- [x] **Step 1: Write failing tests**

Add tests asserting that `max_threads` caps a larger request, invalid caps fail,
oversized weighted workloads fail before object construction, and
`max_workload = Inf` opts out.

- [x] **Step 2: Verify RED**

Run:

```powershell
Rscript -e ".libPaths(c(normalizePath('.rlib'), .libPaths())); devtools::load_all('.'); testthat::test_file('tests/testthat/test-parallel.R')"
```

Expected: failures for missing control fields and helpers.

- [x] **Step 3: Implement the policy**

Validate positive integer `max_threads`, positive numeric or infinite
`max_workload`, calculate `n * max(1, draws) * c(1, 4, 4, 4)[family_code]`,
and pass the capped thread count to `TMB::openmp()`.

- [x] **Step 4: Verify GREEN**

Run the focused parallel and Famoye tests; expect zero new failures.

### Task 2: Inference and retention modes

**Files:**
- Create: `R/inference.R`
- Modify: `R/fit_rpbnb_tmb.R`
- Modify: `R/methods.R`
- Modify: `NAMESPACE`
- Test: `tests/testthat/test-inference-memory.R`
- Modify: `tests/testthat/test-parallel.R`
- Modify: `inst/benchmark_parallel.R`

**Interfaces:**
- Produces: `.rpbnb_inference(obj, par, par_names, mode, report_context)`
- Produces: `.rpbnb_natural_report(par, report_context)`
- Produces: `summary.rpbnb_sdreport(object, select = c("all", "fixed", "report"), ...)`
- `fit$vcov` is a full matrix only in full mode; `fit$vcov_diag` stores diagonal mode.

- [x] **Step 1: Write failing mode tests**

Test all three inference modes, `vcov()` materialization, compact sdreport
summaries, and all three retention modes. Assert that default postfit drops
`obj`, full retains it, compact drops `X1`, `X2`, and `rp_meta`, and no mode
calls `TMB::sdreport()`.

- [x] **Step 2: Verify RED**

Run the new test file against the baseline package and confirm failures are
caused by missing arguments and compatibility class.

- [x] **Step 3: Implement inference helpers**

Use `optimHess`, a symmetric Cholesky check, full or columnwise solves, and
sparse one-parameter delta gradients. Construct an `rpbnb_sdreport` list with
fixed/report estimates and standard errors.

- [x] **Step 4: Integrate retention modes**

Use `match.arg()` in `fit_rpbnb_tmb()`, build the result from shared compact
fields, and conditionally add postfit or full fields. Update `vcov()` and
marginal-effect covariance access. Request `keep = "full"` in benchmark code
that evaluates stored objectives.

- [x] **Step 5: Verify GREEN**

Reinstall and run inference, fit, method, and parallel tests.

### Task 3: Streaming marginal effects

**Files:**
- Modify: `R/marginal_effects.R`
- Test: `tests/testthat/test-inference-memory.R`
- Test: `tests/testthat/test-fit-famoye.R`

**Interfaces:**
- Produces: `.draw_mean_mu(xb, xr, dev)` returning an `n`-vector
- Produces: `.draw_mean_binary(xb, xr, dev)` returning an `n`-vector

- [x] **Step 1: Write failing allocation-structure tests**

Add source-structure assertions forbidding `matrix(0, xu, R)` and
`vapply(seq_len(R), ...)` in the marginal-effect implementation. Add numerical
equivalence tests against a small explicit `n * draws` reference.

- [x] **Step 2: Verify RED**

Run the new focused tests and confirm the source-structure assertions fail.

- [x] **Step 3: Implement streaming sums**

Accumulate `mu_sum <- mu_sum + exp(eta)` and derivative sums inside draw loops,
then divide by `R`. Use the same helpers from `.estimand_t()` so numerical
Jacobian evaluations remain streaming.

- [x] **Step 4: Verify GREEN**

Run marginal-effect, inference, and Famoye tests; expect matching estimates and
no dense draw matrices.

### Task 4: Atomic Gaussian copula quadrature

**Files:**
- Modify: `src/rpbnb.tmb.cpp`
- Modify: `tests/testthat/test-parallel.R`

**Interfaces:**
- Produces: `gaussian_bvn_quadrature(h, k, rho)`
- Produces: `gaussian_bvn_quadrature_gradient(h, k, rho)`
- Produces: `TMB_ATOMIC_VECTOR_FUNCTION(gaussian_bvn_atomic, ...)`

- [x] **Step 1: Write failing atomic tests**

Require a known-derivative atomic in the source. Extend the Gaussian fixture
test to compare serial and two-thread objective, gradient, and Hessian,
including serial-first construction in the same R process.

- [x] **Step 2: Verify RED**

Run the focused source test and confirm the missing atomic registration fails.

- [x] **Step 3: Implement the atomic wrapper**

Move the existing quadrature unchanged to a top-level function, encode its
first derivatives in a thread-safe atomic reverse rule, and replace four direct
quadrature calls with atomic calls using inputs `(h, k, rho)`.

- [x] **Step 4: Rebuild and verify GREEN**

Run `R CMD INSTALL --preclean`, then the complete parallel test file. Objective,
gradient, and Hessian tolerances are `1e-10`, `1e-8`, and `1e-6`.

### Task 5: Documentation and complete verification

**Files:**
- Modify: `README.md`
- Modify: `R/fit_rpbnb_tmb.R`
- Modify: `R/utilities.R`

- [x] **Step 1: Document controls**

Document `inference`, `keep`, `max_threads`, `max_workload`, compact-fit
limitations, and `keep = "full"` for low-level TMB diagnostics.

- [x] **Step 2: Rebuild documentation**

Run `devtools::document()` and confirm `NAMESPACE` and Rd files reflect the
public arguments.

- [x] **Step 3: Install from clean source**

Run:

```powershell
R CMD INSTALL --preclean --library=.rlib .
```

Expected: successful C++17/OpenMP build and package load.

- [x] **Step 4: Run complete tests**

Run the full installed-package test suite. Expect all tests except the two
recorded `test-against-rpbnb.R` external-fixture failures to pass.

- [x] **Step 5: Inspect changes**

Run `git diff --check`, `git status --short`, and a focused diff review. Do not
commit unless the user requests a commit.
