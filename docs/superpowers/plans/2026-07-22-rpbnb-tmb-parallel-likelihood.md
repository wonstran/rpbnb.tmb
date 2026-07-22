# AD-Safe Parallel Likelihood Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add genuine multi-core execution to a single `fit_rpbnb_tmb()` fit while preserving the finite, AD-safe likelihood and producing numerically equivalent objectives, gradients, coefficients, and standard errors.

**Architecture:** TMB owns the OpenMP thread pool through `TMB::openmp()`. The R layer validates and records the realized thread count. The C++ objective precomputes draw-level random deviations, then evaluates independent observations through `parallel_accumulator<Type>` with the Monte Carlo integration retained inside each observation. This avoids shared mutable AD state and leaves one-thread execution as the deterministic fallback.

**Tech Stack:** R, testthat, TMB, C++17, CppAD, OpenMP, GNU Make variables supplied by R (`SHLIB_OPENMP_CXXFLAGS`).

## Global Constraints

- Preserve all existing AD-safety fixes, including conditional expressions, finite copula limits, and dynamic log-sum-exp scaling.
- Parallelize across observations only. Do not parallelize Monte Carlo draws independently because each observation's log-sum-exp reduction must stay local and stable.
- Use `parallel_accumulator<Type>` only through `+=` or `-=` operations.
- Treat `n_cores = 1L` as the reference execution path and supported fallback.
- Do not use elapsed-time thresholds in automated tests; timing is too variable for CI.
- Keep unrelated pre-existing working-tree changes untouched.

Before running commands in a fresh PowerShell session on the current machine, make R available without changing the persistent system configuration:

```powershell
$env:Path = 'C:\Program Files\R\R-4.5.1\bin;' + $env:Path
```

---

## Task 1: Validate and Configure the Requested Thread Count

**Files:**

- Modify: `R/utilities.R`
- Modify: `R/tmb_helpers.R`
- Modify: `R/fit_rpbnb_tmb.R`
- Create: `tests/testthat/test-parallel.R`

- [ ] **Step 1: Write failing validation tests**

Add focused tests showing that `rpbnb_tmb_control()` rejects invalid core counts and normalizes valid scalar integers.

```r
test_that("rpbnb_tmb_control validates n_cores", {
  expect_error(rpbnb_tmb_control(n_cores = 0), "n_cores")
  expect_error(rpbnb_tmb_control(n_cores = -1), "n_cores")
  expect_error(rpbnb_tmb_control(n_cores = NA_integer_), "n_cores")
  expect_error(rpbnb_tmb_control(n_cores = 1.5), "n_cores")
  expect_error(rpbnb_tmb_control(n_cores = "2"), "n_cores")
  expect_error(rpbnb_tmb_control(n_cores = c(1, 2)), "n_cores")
  expect_equal(rpbnb_tmb_control(n_cores = 2)$n_cores, 2L)
})
```

- [ ] **Step 2: Run the focused test and confirm it fails**

```powershell
Rscript -e "devtools::load_all(quiet = TRUE); testthat::test_file('tests/testthat/test-parallel.R')"
```

Expected: at least one invalid value is silently coerced or accepted.

- [ ] **Step 3: Implement strict control validation**

In `rpbnb_tmb_control()`, require `n_cores` to be one finite, whole number greater than or equal to one before coercing it to integer. Reject values larger than `.Machine$integer.max`.

- [ ] **Step 4: Add a small internal TMB thread configurator**

In `R/tmb_helpers.R`, add `.configure_tmb_threads(n_cores, DLL = "rpbnb.tmb")` that:

1. obtains the supported maximum with `TMB::openmp(max = TRUE, DLL = DLL)`;
2. treats a missing, non-finite, or less-than-one maximum as one;
3. chooses `realized <- min(n_cores, supported)`;
4. calls `TMB::openmp(n = realized, DLL = DLL)`;
5. warns when the request is reduced, including requested and realized values; and
6. returns the realized integer directly for storage.

Do not alter global OpenMP state before the package DLL is loaded.

- [ ] **Step 5: Extract construction of the TMB object**

Move the direct `TMB::MakeADFun()` call from `fit_rpbnb_tmb()` into `.make_rpbnb_tmb_object(data, parameters, map = NULL, random = NULL, silent = TRUE, n_cores = 1L)`. Configure threads immediately before calling `MakeADFun()`. Return exactly `list(obj = obj, n_cores = realized)`.

This helper is intentionally internal so tests can compare identical fixtures at different thread counts without running the optimizer.

- [ ] **Step 6: Store the realized count in fitted objects**

Set `fit$parallel <- list(requested = control$n_cores, realized = configured$n_cores)` without changing existing result fields. The stored value must describe the thread setting used to construct and evaluate that fit.

- [ ] **Step 7: Run the focused test until green**

```powershell
Rscript -e "devtools::load_all(quiet = TRUE); testthat::test_file('tests/testthat/test-parallel.R')"
```

- [ ] **Step 8: Commit Task 1**

```powershell
git add R/utilities.R R/tmb_helpers.R R/fit_rpbnb_tmb.R tests/testthat/test-parallel.R
git commit -m "feat: validate and configure TMB threads"
```

---

## Task 2: Compile and Link the Package with OpenMP

**Files:**

- Modify: `src/Makevars`
- Modify: `src/Makevars.win`
- Modify: `src/rpbnb.tmb.cpp`
- Modify: `tests/testthat/test-parallel.R`

- [ ] **Step 1: Add a failing compiled-capability test**

Build a deterministic internal TMB fixture and compare `obj$report()$openmp_compiled` with whether `R CMD config SHLIB_OPENMP_CXXFLAGS` returns a non-empty flag. On the current toolchain the expected value is one. This field does not exist yet, so the test records the missing build guarantee rather than trusting the configurable thread count alone.

Restore the prior TMB thread setting with `on.exit()` so this test cannot leak state into later tests or add a new test dependency.

- [ ] **Step 2: Confirm the test fails before the build contract exists**

Run the focused test against a clean installed package, not only `load_all()`.

```powershell
$lib = 'C:/tmp/rpbnb-par-lib'
New-Item -ItemType Directory -Force -Path $lib | Out-Null
R CMD INSTALL --preclean --library=$lib .
Rscript -e ".libPaths(c('C:/tmp/rpbnb-par-lib', .libPaths())); library(rpbnb.tmb); testthat::test_file('tests/testthat/test-parallel.R')"
```

Expected: `openmp_compiled` is absent.

- [ ] **Step 3: Add portable OpenMP build flags**

In both make files, append the R toolchain-provided flag to compilation and linking:

```make
PKG_CXXFLAGS += $(SHLIB_OPENMP_CXXFLAGS)
PKG_LIBS += $(SHLIB_OPENMP_CXXFLAGS)
```

Retain the Windows assembler workaround `-Wa,-mbig-obj` already needed for this large TMB translation unit. Do not hard-code `-fopenmp`; R supplies the platform-appropriate value.

- [ ] **Step 4: Report the target DLL's compile-time OpenMP capability**

In `src/rpbnb.tmb.cpp`, report a scalar integer without changing the likelihood:

```cpp
int openmp_compiled = 0;
#ifdef _OPENMP
openmp_compiled = 1;
#endif
REPORT(openmp_compiled);
```

This distinguishes the target model DLL's compile mode from `TMB::openmp(max = TRUE)`, whose maximum comes from the installed TMB package.

- [ ] **Step 5: Perform a clean rebuild and rerun the test**

```powershell
$lib = 'C:/tmp/rpbnb-par-lib'
R CMD INSTALL --preclean --library=$lib .
Rscript -e ".libPaths(c('C:/tmp/rpbnb-par-lib', .libPaths())); library(rpbnb.tmb); testthat::test_file('tests/testthat/test-parallel.R')"
```

Expected on this OpenMP-capable R toolchain: `openmp_compiled` equals one, and requesting two realizes two when the machine maximum is at least two. When `SHLIB_OPENMP_CXXFLAGS` is empty, expect `openmp_compiled == 0L` and skip only assertions that require multiple threads.

- [ ] **Step 6: Commit Task 2**

```powershell
git add src/Makevars src/Makevars.win src/rpbnb.tmb.cpp tests/testthat/test-parallel.R
git commit -m "build: enable OpenMP for TMB objective"
```

---

## Task 3: Restructure the Objective for Observation-Level Parallelism

**Files:**

- Modify: `src/rpbnb.tmb.cpp`
- Modify: `R/fit_rpbnb_tmb.R`
- Modify: `tests/testthat/test-parallel.R`

- [ ] **Step 1: Add structural and numerical parallelization tests**

Build one deterministic, small fixture with fixed Monte Carlo draws and no optimization. Construct the serial TMB object through `.make_rpbnb_tmb_object()`, immediately capture its objective and gradient, then construct the two-thread object and capture its values. Do not retain both objects and evaluate them after reconfiguration because TMB's thread setting is DLL-global. Compare the captured values at the same parameter vector:

```r
expect_equal(obj_parallel$fn(par), obj_serial$fn(par), tolerance = 1e-10)
expect_equal(obj_parallel$gr(par), obj_serial$gr(par), tolerance = 1e-8)
expect_true(all(is.finite(obj_parallel$gr(par))))
```

Cover at least the Gaussian and Frank copulas because they exercise different dependence transformations and the previously repaired Frank singularity. Skip only the two-thread comparison when the runtime truly lacks OpenMP.

Also add a source-level contract asserting that `src/rpbnb.tmb.cpp` declares `parallel_accumulator<Type>`. Numerical equality alone cannot prove that multiple likelihood regions exist, because a serial objective would also return equal numbers under a two-thread configuration.

- [ ] **Step 2: Run the test and establish the red state**

```powershell
Rscript -e "devtools::load_all(quiet = TRUE); testthat::test_file('tests/testthat/test-parallel.R')"
```

Expected: the structural contract fails because the objective still lacks `parallel_accumulator<Type>`.

- [ ] **Step 3: Remove the unused per-fit thread value from TMB data**

Delete `DATA_INTEGER(n_cores)` and the `_OPENMP`/`omp_set_num_threads(n_cores)` block from `src/rpbnb.tmb.cpp`, then remove the corresponding entry and argument from `.build_tmb_data()`. Thread configuration is owned by `TMB::openmp()`, not objective data, so retaining both would create two competing sources of truth.

- [ ] **Step 4: Precompute draw-level random deviations**

Before entering the observation loop, transform the supplied draws into `matrix<Type> dev1(R, q1)` and `matrix<Type> dev2(R, q2)`. Preserve the existing distribution transforms and draw order exactly. These matrices are read-only inside the parallel region.

Verify that this refactor alone gives the same one-thread objective and gradient as the pre-refactor code using a temporary saved reference value in the test fixture.

- [ ] **Step 5: Reorganize the likelihood from draw-first to observation-first**

Replace the shared `matrix<Type> ll_all(n, R)` workflow with this ownership pattern:

```cpp
parallel_accumulator<Type> nll(this);

for (int i = 0; i < n; ++i) {
  vector<Type> log_draw(R);

  for (int r = 0; r < R; ++r) {
    // Compute this observation/draw's marginal parameters and copula mass.
    log_draw(r) = log_joint_pmf;
  }

  Type max_log = log_draw(0);
  for (int r = 1; r < R; ++r) {
    max_log = CppAD::CondExpGt(log_draw(r), max_log, log_draw(r), max_log);
  }

  Type scaled_sum = Type(0);
  for (int r = 0; r < R; ++r) {
    scaled_sum += exp(log_draw(r) - max_log);
  }

  Type log_contribution = max_log + log(scaled_sum) - log(Type(R));
  nll -= log_contribution;
}

return nll;
```

Keep all working values local to the observation iteration. Do not write to shared matrices, counters, or diagnostics in the loop.

- [ ] **Step 6: Preserve every family and copula branch exactly**

Move the existing marginal PMF and copula rectangle-mass calculations into the new loop without algebraic simplification. In particular, retain:

- `CppAD::CondExp*` clamps and comparisons;
- Gaussian dependence computed on the AD tape;
- Frank's nonzero default and near-zero limit;
- Clayton boundary handling;
- finite probability floors before `log()`; and
- dynamic log-sum-exp maximum selection.

- [ ] **Step 7: Run focused numerical-equivalence tests**

```powershell
Rscript -e "devtools::load_all(quiet = TRUE); testthat::test_file('tests/testthat/test-parallel.R')"
```

Expected: serial and parallel objective/gradient comparisons pass for every fixture.

- [ ] **Step 8: Run existing family and copula regressions**

```powershell
Rscript -e "devtools::load_all(quiet = TRUE); testthat::test_file('tests/testthat/test-fit-copula.R'); testthat::test_file('tests/testthat/test-fit-famoye.R')"
```

Expected: no reappearance of NA curvature or standard errors.

- [ ] **Step 9: Commit Task 3**

```powershell
git add src/rpbnb.tmb.cpp R/fit_rpbnb_tmb.R tests/testthat/test-parallel.R
git commit -m "feat: parallelize TMB likelihood by observation"
```

---

## Task 4: Verify Fitted-Model Equivalence and Expose Diagnostics

**Files:**

- Modify: `R/fit_rpbnb_tmb.R`
- Modify: `R/methods.R`
- Modify: `tests/testthat/test-parallel.R`

- [ ] **Step 1: Add a fitted-model equivalence test**

Fit the same seeded, compact data set with `n_cores = 1L` and `n_cores = 2L`. Use the same draws and starting values. Assert:

- coefficient vectors agree with absolute tolerance `1e-6`;
- reported standard errors agree within `1e-5`;
- objective values agree within `1e-7`;
- all coefficients and standard errors are finite; and
- the parallel metadata records the actual thread count.

If the runtime supports only one thread, keep the serial assertions and skip only the two-thread equivalence portion.

- [ ] **Step 2: Add a failing print-method test**

Capture the printed fitted object and expect the exact line `TMB threads: 2` using the realized value, not merely the request.

- [ ] **Step 3: Show the realized thread count in output**

Update the relevant print method to display `fit$parallel$realized`. Keep backward compatibility for older fitted objects that do not have a `parallel` field by omitting the line or treating the value as one.

- [ ] **Step 4: Run focused and existing tests**

```powershell
Rscript -e "devtools::load_all(quiet = TRUE); testthat::test_file('tests/testthat/test-parallel.R'); testthat::test_file('tests/testthat/test-fit-copula.R'); testthat::test_file('tests/testthat/test-fit-famoye.R')"
```

- [ ] **Step 5: Commit Task 4**

```powershell
git add R/fit_rpbnb_tmb.R R/methods.R tests/testthat/test-parallel.R
git commit -m "test: verify parallel fit equivalence"
```

---

## Task 5: Add a Reproducible Parallel Benchmark and Update Demos

**Files:**

- Create: `inst/benchmark_parallel.R`
- Modify: `inst/fit_rpbnb_diff_copula.R`
- Modify: `inst/fit_rpbnb_diff_famoye.R`

- [ ] **Step 1: Create a benchmark script**

Use a seeded 500-observation, 20-draw Frank model, matching the NA regression workload. Fit once with one core and once with a detected count capped at four, while reusing identical data, draws, and starting values. If core detection returns `NA`, use one.

The script must print:

- requested and realized cores for each fit;
- elapsed seconds;
- speedup as `serial_elapsed / parallel_elapsed`;
- maximum absolute coefficient difference; and
- maximum absolute standard-error difference.

Stop with an error if coefficients or standard errors are non-finite or differ beyond the tolerances established in Task 4. Report speedup but do not fail solely because speedup is below one on a small or busy machine.

- [ ] **Step 2: Update example scripts to exercise available cores safely**

Define a bounded example setting:

```r
detected_cores <- parallel::detectCores(logical = FALSE)
demo_cores <- if (is.na(detected_cores)) 1L else max(1L, min(4L, detected_cores))
```

Pass it through `rpbnb_tmb_control(n_cores = demo_cores)`. Keep the examples runnable on one-core systems.

- [ ] **Step 3: Run the benchmark once**

```powershell
Rscript inst/benchmark_parallel.R
```

Expected: finite results, numerical equivalence, realized multi-core use when available, and an informational timing comparison.

- [ ] **Step 4: Commit Task 5**

```powershell
git add inst/benchmark_parallel.R inst/fit_rpbnb_diff_copula.R inst/fit_rpbnb_diff_famoye.R
git commit -m "perf: add parallel likelihood benchmark"
```

---

## Task 6: Clean-Build Verification and Handoff

**Files:**

- Verify all files changed in Tasks 1-5
- Modify only files required by failures uncovered here

- [ ] **Step 1: Install from a clean precompiled state**

```powershell
$lib = 'C:/tmp/rpbnb-par-lib'
New-Item -ItemType Directory -Force -Path $lib | Out-Null
R CMD INSTALL --preclean --library=$lib .
```

Expected: the large Windows TMB object compiles with the big-object assembler flag and links with the OpenMP runtime.

- [ ] **Step 2: Run the installed-package focused suites**

```powershell
Rscript -e ".libPaths(c('C:/tmp/rpbnb-par-lib', .libPaths())); library(rpbnb.tmb); testthat::test_file('tests/testthat/test-parallel.R'); testthat::test_file('tests/testthat/test-fit-copula.R'); testthat::test_file('tests/testthat/test-fit-famoye.R')"
```

- [ ] **Step 3: Run package-level checks that are reliable in the current repository**

```powershell
Rscript -e "devtools::test()"
git diff --check
git status --short
```

Record any pre-existing package-check failures separately from regressions introduced by this work.

- [ ] **Step 4: Run the installed-package benchmark**

```powershell
Rscript inst/benchmark_parallel.R
```

Record the machine's logical/physical core count, realized TMB thread count, serial time, parallel time, and numerical differences in the handoff.

- [ ] **Step 5: Review the final diff for scope and AD safety**

Confirm that:

- no raw `#pragma omp` was added around AD code;
- no parameter-dependent ordinary C++ `if` was introduced;
- no shared mutable buffer is written by the observation loop;
- every requested thread count is validated and the realized count is visible;
- one-thread execution remains supported; and
- unrelated user changes remain intact.

- [ ] **Step 6: Commit only verification-driven adjustments, if any**

Stage each verification-adjusted file by its explicit path, review `git diff --cached`, then commit with `git commit -m "fix: address parallel verification findings"`.

Skip this commit when verification required no source changes.

---

## Completion Evidence

Implementation is complete only when the handoff includes:

1. the clean install command and successful result;
2. focused serial/parallel objective and gradient equivalence results;
3. fitted coefficient and standard-error equivalence results;
4. the existing Frank/NA regression result with finite curvature and SEs;
5. the realized TMB thread count for a request above one;
6. benchmark timings without claiming a guaranteed speedup; and
7. `git diff --check` output plus a scoped working-tree summary.
