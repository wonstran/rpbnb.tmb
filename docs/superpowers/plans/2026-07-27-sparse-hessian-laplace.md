# Sparse-Hessian Laplace Estimation Path Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add `method = "laplace"` to `fit_rpbnb_tmb()`, integrating random coefficients with TMB's sparse-Hessian Laplace approximation instead of Halton simulation, so tape size scales with `n` rather than `n x draws`.

**Architecture:** One C++ template, one DLL, switched by a new `est_method` data field. New `PARAMETER_MATRIX(u1)`/`u2` hold one latent standard-normal row per observation. The Laplace branch reuses every existing density function and changes only two things inside the observation loop: deviations read from the latents rather than the precomputed draw matrix, and the log-sum-exp reduction becomes a single contribution plus `dnorm` priors. Under SML the latents are `map`-fixed and never read, so that tape is unchanged.

**Tech Stack:** R 4.5.1, TMB 1.9.21, CppAD/TMBad, testthat 3e, devtools, Rtools (Windows).

## Global Constraints

- Spec: `docs/superpowers/specs/2026-07-27-sparse-hessian-laplace-design.md`. Read it before Task 1.
- `Rscript` is at `C:\Program Files\R\R-4.5.1\bin\Rscript.exe`. Every command below uses `$RS` as shorthand; in PowerShell set it once per session:
  ```powershell
  $RS = "C:\Program Files\R\R-4.5.1\bin\Rscript.exe"
  ```
- Package name / DLL name is `rpbnb.tmb` (dot, not underscore). The repo directory is `rpbnb_tmb`.
- TMB matches parameters by name but `checkParameterOrder = TRUE` is on by default: the R `parameters` list order MUST match the `PARAMETER*` declaration order in the template. New parameters go last in both.
- Default `method` is `"sml"`. No existing call site, test, demo script, or documented behavior may change.
- Laplace supports `"normal"` and `"lognormal"` random-coefficient distributions only.
- `draws` is NOT ignored or warned about under Laplace. It still sizes the Famoye lambda-bound loop and the `rp_meta$Z1`/`Z2` grids that `predict()` and marginal effects average over.
- Compilation of this template needs `-Wa,-mbig-obj` on Windows; it is already in `src/Makevars.win`. Do not remove it.
- Commit after every task.

---

### Task 1: Add latent parameters and the `est_method` switch without changing SML behavior

This task adds the new template inputs and wires them through R, but implements **no** Laplace logic. Its entire deliverable is: the template accepts the new inputs, and every existing fit produces a bit-identical result. This is the safety net that makes Task 2 reviewable.

**Files:**
- Modify: `src/rpbnb.tmb.cpp:216` (after `DATA_SCALAR(lamHi)`), `src/rpbnb.tmb.cpp:225` (after `PARAMETER(z_dep)`)
- Modify: `R/tmb_helpers.R:102-118` (`.build_tmb_data`)
- Modify: `R/fit_rpbnb_tmb.R:229-264` (data build, map, parameters list)
- Test: `tests/testthat/test-laplace.R` (create)

**Interfaces:**
- Consumes: nothing from earlier tasks.
- Produces: `.build_tmb_data(..., lamLo, lamHi, est_method)` — `est_method` is a new trailing required argument, integer `0L` or `1L`, placed in the returned list as `est_method`. Template data field `est_method`; template parameters `u1` (n x q1 matrix) and `u2` (n x q2 matrix), declared after `z_dep`.

- [ ] **Step 1: Capture the SML baseline log-likelihood**

This is a regression test, so the assertion value must come from the code as it exists *now*, before any edit. Run:

```powershell
& $RS -e "devtools::load_all('.'); set.seed(1); n <- 300; x1 <- rnorm(n); X <- cbind(`(Intercept)`=1, x1=x1); mu <- exp(X %*% c(0.3,0.5)); dat <- data.frame(y1=rnbinom(n, mu=mu, size=1/0.4), y2=rnbinom(n, mu=mu, size=1/0.5), x1=x1); f <- fit_rpbnb_tmb(y1~x1, y2~x1, data=dat, random_1='x1', dependence='famoye', draws=50, seed=7); cat(sprintf('%.10f', f\$logLik), '\n')"
```

Record the printed number. It is referred to below as `<BASELINE>`. Paste the literal value into the test in Step 2 — do not round it and do not substitute a value from any other run.

- [ ] **Step 2: Write the SML regression test**

Create `tests/testthat/test-laplace.R`:

```r
# The SML path must be bit-identical before and after the Laplace feature.
# BASELINE was captured from the pre-change code; see the implementation plan.
BASELINE_SML_LOGLIK <- <BASELINE>

sml_baseline_fit <- function() {
  set.seed(1)
  n <- 300
  x1 <- rnorm(n)
  X <- cbind(`(Intercept)` = 1, x1 = x1)
  mu <- exp(X %*% c(0.3, 0.5))
  dat <- data.frame(
    y1 = rnbinom(n, mu = mu, size = 1 / 0.4),
    y2 = rnbinom(n, mu = mu, size = 1 / 0.5),
    x1 = x1
  )
  fit_rpbnb_tmb(y1 ~ x1, y2 ~ x1, data = dat, random_1 = "x1",
                dependence = "famoye", draws = 50, seed = 7)
}

test_that("adding latent parameters leaves the SML tape unchanged", {
  skip_on_cran()
  fit <- sml_baseline_fit()
  expect_equal(fit$logLik, BASELINE_SML_LOGLIK, tolerance = 1e-10)
})
```

- [ ] **Step 3: Run the regression test against unmodified code**

```powershell
& $RS -e "devtools::load_all('.'); testthat::test_file('tests/testthat/test-laplace.R')"
```

Expected: PASS. If it fails here the baseline was mis-captured — fix that before editing anything.

- [ ] **Step 4: Add the data field to the template**

In `src/rpbnb.tmb.cpp`, immediately after `DATA_SCALAR(lamHi);` (line 216), add:

```cpp
  // 0 = simulated maximum likelihood (Halton draws)
  // 1 = Laplace approximation (latent u1/u2 integrated by TMB)
  DATA_INTEGER(est_method);
```

The field is named `est_method` rather than `method` to stay visibly distinct from `MakeADFun()`'s own `method` argument, which selects the outer optimizer and is unrelated.

- [ ] **Step 5: Add the latent parameters to the template**

In `src/rpbnb.tmb.cpp`, immediately after `PARAMETER(z_dep);` (line 225), add:

```cpp
  // Latent standard normals, one row per observation. Under est_method == 0
  // these are map-fixed at zero in R and never read; under est_method == 1
  // they are TMB random effects.
  PARAMETER_MATRIX(u1);
  PARAMETER_MATRIX(u2);
```

Then, so an unused-variable warning does not appear while Task 2 is still pending, add immediately after the `int R = ...` line (line 233):

```cpp
  (void)u1; (void)u2; (void)est_method;
```

This line is deleted in Task 2.

- [ ] **Step 6: Add `est_method` to the data builder**

In `R/tmb_helpers.R`, change the `.build_tmb_data` signature (line 102-105) to take a trailing `est_method`:

```r
.build_tmb_data <- function(Y1, Y2, X1, X2, rand_idx1, rand_idx2,
                            Z1, Z2, dist1, dist2, sign1, sign2,
                            family_code, pois1, pois2,
                            lamLo, lamHi, est_method) {
```

and add the field to the returned list, after `lamHi` (line 116):

```r
    lamLo = as.numeric(lamLo), lamHi = as.numeric(lamHi),
    est_method = as.integer(est_method)
```

- [ ] **Step 7: Pass `est_method` and the latent parameters from the fitting function**

In `R/fit_rpbnb_tmb.R`, change the `.build_tmb_data()` call (lines 229-232) to:

```r
  tmb_data <- .build_tmb_data(Y1, Y2, X1, X2, rand_idx1, rand_idx2,
                              Z1, Z2, dist1, dist2, sign1, sign2,
                              family_code, poisson_1, poisson_2,
                              lamLo, lamHi, est_method = 0L)
```

After the independence `map[["z_dep"]]` block (line 243), add the latent map entries:

```r
  # SML: the latents are tape constants at zero and are never read by the
  # template. Guarded on non-zero length because TMB rejects a map entry for a
  # zero-length parameter.
  if (n * q1 > 0) map[["u1"]] <- factor(rep(NA_integer_, n * q1))
  if (n * q2 > 0) map[["u2"]] <- factor(rep(NA_integer_, n * q2))
```

In the `parameters` list (lines 250-258), add `u1` and `u2` **last**, after `z_dep`, to match the template declaration order:

```r
    parameters = list(
      beta1 = start[i1],
      beta2 = start[i2],
      log_sd1 = if (q1 > 0) start[(k1 + k2 + 1):(k1 + k2 + q1)] else numeric(0),
      log_sd2 = if (q2 > 0) start[(k1 + k2 + q1 + 1):(k1 + k2 + q1 + q2)] else numeric(0),
      log_m1 = start[idx_end + 1],
      log_m2 = start[idx_end + 2],
      z_dep = if (family_code >= 0L) start[idx_end + 3] else 0,
      u1 = matrix(0, n, q1),
      u2 = matrix(0, n, q2)
    ),
```

- [ ] **Step 8: Recompile and re-run the regression test**

```powershell
& $RS -e "devtools::load_all('.', recompile = TRUE); testthat::test_file('tests/testthat/test-laplace.R')"
```

Expected: PASS, with the same log-likelihood to 1e-10. A failure here means the new parameters reached the SML tape, which they must not.

- [ ] **Step 9: Run the full existing suite**

```powershell
& $RS -e "devtools::load_all('.'); testthat::test_dir('tests/testthat')"
```

Expected: all pass. This is also what exercises the `q1 == q2 == 0` case (`test-fit-famoye.R:12` fits with no random coefficients), confirming TMB accepts the zero-column `u1`/`u2` matrices.

- [ ] **Step 10: Commit**

```bash
git add src/rpbnb.tmb.cpp R/tmb_helpers.R R/fit_rpbnb_tmb.R tests/testthat/test-laplace.R
git commit -m "feat: add est_method switch and latent parameters, SML unchanged"
```

---

### Task 2: Implement the Laplace branch in the C++ template

**Files:**
- Modify: `src/rpbnb.tmb.cpp:233` (drop the `(void)` line, force `R = 1`), `:336-350` (guard `dev1`/`dev2` precomputation), `:363-368` (deviation source), `:504-514` (reduction)

**Interfaces:**
- Consumes: `est_method`, `u1`, `u2` from Task 1.
- Produces: template behavior for `est_method == 1` — the joint negative log-likelihood over data and latents, suitable for `MakeADFun(random = c("u1","u2"))`. No R-visible signature change.

- [ ] **Step 1: Force a single "draw" under Laplace**

In `src/rpbnb.tmb.cpp`, replace the `(void)u1; (void)u2; (void)est_method;` line added in Task 1 with:

```cpp
  // Laplace evaluates the conditional density once at the current latent
  // values; there is no draw dimension to average over.
  if (est_method == 1) R = 1;
```

`R` is declared on the preceding line (`int R = (q1 + q2 > 0) ? Z1.rows() : 1;`) and is already non-const, so this is a plain reassignment.

- [ ] **Step 2: Guard the per-draw deviation precomputation**

The `dev1`/`dev2` matrices at lines 336-350 are indexed by draw and read `Z1`/`Z2`, which are dummy matrices under Laplace. Wrap the loop:

```cpp
  matrix<Type> dev1(R, q1), dev2(R, q2);
  if (est_method == 0) {
    for (int r = 0; r < R; r++) {
      for (int j = 0; j < q1; j++) {
        Type base = u_to_base(Z1(r, j), dist1(j));
        int col = rand_idx1(j);
        dev1(r, j) = compute_dev(beta1(col), sd1(j), base,
                                 dist1(j), sign1(j));
      }
      for (int j = 0; j < q2; j++) {
        Type base = u_to_base(Z2(r, j), dist2(j));
        int col = rand_idx2(j);
        dev2(r, j) = compute_dev(beta2(col), sd2(j), base,
                                 dist2(j), sign2(j));
      }
    }
  }
```

- [ ] **Step 3: Read deviations from the latents under Laplace**

Replace the two inner `eta` accumulation loops (lines 363-368) with:

```cpp
      for (int j = 0; j < q1; j++) {
        int col = rand_idx1(j);
        // u_to_base() is deliberately NOT applied to the latent: it maps a
        // uniform draw to a standard-normal base, and the latent already IS
        // that base.
        Type d = (est_method == 1)
          ? compute_dev(beta1(col), sd1(j), u1(i, j), dist1(j), sign1(j))
          : dev1(r, j);
        eta1 += X1(i, col) * d;
      }
      for (int j = 0; j < q2; j++) {
        int col = rand_idx2(j);
        Type d = (est_method == 1)
          ? compute_dev(beta2(col), sd2(j), u2(i, j), dist2(j), sign2(j))
          : dev2(r, j);
        eta2 += X2(i, col) * d;
      }
```

`est_method` is a compile-time-invariant int, so the branch hoists out of the loop; keeping one loop avoids duplicating the linear-predictor code.

- [ ] **Step 4: Branch the reduction**

Replace the log-sum-exp block (lines 504-514) with:

```cpp
    if (est_method == 1) {
      nll -= log_draw(0);
      for (int j = 0; j < q1; j++)
        nll -= dnorm(u1(i, j), Type(0), Type(1), true);
      for (int j = 0; j < q2; j++)
        nll -= dnorm(u2(i, j), Type(0), Type(1), true);
    } else {
      Type max_log = log_draw(0);
      for (int r = 1; r < R; r++) {
        max_log = CppAD::CondExpGt(log_draw(r), max_log,
                                   log_draw(r), max_log);
      }
      Type scaled_sum = Type(0);
      for (int r = 0; r < R; r++) {
        scaled_sum += exp(log_draw(r) - max_log);
      }
      Type log_contribution = max_log + log(scaled_sum) - logR;
      nll -= log_contribution;
    }
```

Everything after this — the `ADREPORT` calls and the Kendall's tau block — is unchanged and shared by both estimators.

- [ ] **Step 5: Recompile and confirm SML is still unchanged**

```powershell
& $RS -e "devtools::load_all('.', recompile = TRUE); testthat::test_file('tests/testthat/test-laplace.R')"
```

Expected: PASS at the same baseline. The Laplace branch is now present but unreachable from R until Task 3.

- [ ] **Step 6: Run the full suite**

```powershell
& $RS -e "devtools::load_all('.'); testthat::test_dir('tests/testthat')"
```

Expected: all pass.

- [ ] **Step 7: Commit**

```bash
git add src/rpbnb.tmb.cpp
git commit -m "feat: implement Laplace branch in TMB template"
```

---

### Task 3: Wire `method = \"laplace\"` through the R fitting function

**Files:**
- Modify: `R/fit_rpbnb_tmb.R:70-77` (signature), `:91` (match.arg), `:138-154` (validation and workload), `:229-232` (dummy Z), `:234-243` (map guard), `:248-264` (random argument), `:313-349` (result field)

**Interfaces:**
- Consumes: `est_method` data field and Laplace template branch from Tasks 1-2; `.make_rpbnb_tmb_object(random =)` at `R/tmb_helpers.R:45`, which already exists and has never had a non-`NULL` caller.
- Produces: `fit_rpbnb_tmb(..., method = c("sml", "laplace"))`. The returned `rpbnb_tmb_fit` gains a `method` element holding the matched string.

- [ ] **Step 1: Write the failing validation tests**

Append to `tests/testthat/test-laplace.R`:

```r
test_that("laplace rejects unsupported random-coefficient distributions", {
  skip_on_cran()
  set.seed(3)
  n <- 120
  x1 <- rnorm(n)
  dat <- data.frame(y1 = rpois(n, 2), y2 = rpois(n, 2), x1 = x1)
  expect_error(
    fit_rpbnb_tmb(y1 ~ x1, y2 ~ x1, data = dat,
                  random_1 = list(x1 = "uniform"),
                  dependence = "independence", draws = 20,
                  method = "laplace"),
    "uniform"
  )
})

test_that("laplace requires at least one random coefficient", {
  skip_on_cran()
  set.seed(4)
  n <- 120
  x1 <- rnorm(n)
  dat <- data.frame(y1 = rpois(n, 2), y2 = rpois(n, 2), x1 = x1)
  expect_error(
    fit_rpbnb_tmb(y1 ~ x1, y2 ~ x1, data = dat,
                  dependence = "independence", draws = 20,
                  method = "laplace"),
    "random coefficient"
  )
})

test_that("laplace keeps draws meaningful for post-estimation averaging", {
  skip_on_cran()
  set.seed(5)
  n <- 200
  x1 <- rnorm(n)
  dat <- data.frame(y1 = rpois(n, 2), y2 = rpois(n, 2), x1 = x1)
  expect_no_warning(
    fit <- fit_rpbnb_tmb(y1 ~ x1, y2 ~ x1, data = dat, random_1 = "x1",
                         dependence = "independence", draws = 25,
                         method = "laplace")
  )
  expect_identical(nrow(fit$rp_meta$Z1), 25L)
  expect_identical(fit$method, "laplace")
})
```

- [ ] **Step 2: Run to verify they fail**

```powershell
& $RS -e "devtools::load_all('.'); testthat::test_file('tests/testthat/test-laplace.R')"
```

Expected: the three new tests FAIL with an unused-argument error for `method`. The Task 1 regression test still passes.

- [ ] **Step 3: Add the argument and match it**

In `R/fit_rpbnb_tmb.R`, change the signature (lines 70-77) to add `method` **last**, after `poisson_2`:

```r
fit_rpbnb_tmb <- function(formula_1, formula_2, data,
                          random_1 = NULL, random_2 = NULL,
                          draws = 400L, seed = 1234L, start = NULL,
                          dependence = "famoye",
                          control = rpbnb_tmb_control(),
                          inference = c("full", "diag", "none"),
                          keep = c("postfit", "compact", "full"),
                          poisson_1 = FALSE, poisson_2 = FALSE,
                          method = c("sml", "laplace")) {
```

Appending rather than inserting is deliberate: inserting `method` before `poisson_1` would shift the positional index of `poisson_1` and `poisson_2`, silently breaking any positional caller. The global constraint that no existing call site may change binds here.

Next to the existing `keep <- match.arg(keep)` (line 92), add:

```r
  method <- match.arg(method)
```

- [ ] **Step 4: Add Laplace validation**

Immediately after `total_rand <- q1 + q2` (line 140), insert:

```r
  if (identical(method, "laplace")) {
    if (total_rand == 0L) {
      stop("method = \"laplace\" needs at least one random coefficient to ",
           "integrate. Use method = \"sml\", or specify random_1/random_2.",
           call. = FALSE)
    }
    bad_dist <- unique(c(spec1$dist, spec2$dist))
    bad_dist <- bad_dist[!bad_dist %in% c("normal", "lognormal")]
    if (length(bad_dist)) {
      stop("method = \"laplace\" supports only normal and lognormal random ",
           "coefficients; got ", paste(bad_dist, collapse = ", "),
           ". Use method = \"sml\" for these distributions.", call. = FALSE)
    }
  }
```

- [ ] **Step 5: Make the workload guard measure the Laplace tape**

Replace `effective_draws` (line 141, now shifted down by the validation block) with:

```r
  # Laplace tapes one conditional evaluation per observation, so the draw
  # dimension is genuinely absent.  It is not free, though: the cost that
  # replaces it is the latent dimension n * (q1 + q2), which sizes the random
  # -effect vector and its sparse Hessian.  Budgeting `1L` here would leave
  # max_workload bounding nothing at all on the very path introduced to solve
  # a memory problem, so the multiplier is the per-observation latent count.
  effective_draws <- if (identical(method, "laplace")) {
    total_rand
  } else if (total_rand > 0L) {
    draws
  } else {
    1L
  }
```

- [ ] **Step 6: Pass dummy draws to the template only**

Halton generation at line 158 stays exactly as it is — its output still feeds the Famoye lambda-bound loop and `rp_meta`. Change only what reaches the template. Replace the `.build_tmb_data()` call from Task 1 with:

```r
  # Under Laplace the template never indexes a draw dimension; Z1/Z2 stay real
  # in R because the lambda-bound loop and rp_meta still consume them.
  Z1_tmb <- if (identical(method, "laplace")) matrix(0, 1, q1) else Z1
  Z2_tmb <- if (identical(method, "laplace")) matrix(0, 1, q2) else Z2
  tmb_data <- .build_tmb_data(Y1, Y2, X1, X2, rand_idx1, rand_idx2,
                              Z1_tmb, Z2_tmb, dist1, dist2, sign1, sign2,
                              family_code, poisson_1, poisson_2,
                              lamLo, lamHi,
                              est_method = if (identical(method, "laplace")) 1L else 0L)
```

- [ ] **Step 7: Map the latents only under SML**

Replace the map block added in Task 1 Step 7 with:

```r
  # SML: the latents are tape constants at zero and are never read by the
  # template. Laplace: they are random effects and must stay free.
  if (!identical(method, "laplace")) {
    if (n * q1 > 0) map[["u1"]] <- factor(rep(NA_integer_, n * q1))
    if (n * q2 > 0) map[["u2"]] <- factor(rep(NA_integer_, n * q2))
  }
```

- [ ] **Step 8: Declare the random effects**

Immediately before the `.make_rpbnb_tmb_object()` call (line 248), add:

```r
  random_names <- if (identical(method, "laplace")) {
    c(if (q1 > 0) "u1", if (q2 > 0) "u2")
  } else {
    NULL
  }
```

and add the argument to the call, after `map`:

```r
    map = if (length(map) > 0) map else NULL,
    random = random_names,
```

- [ ] **Step 9: Record the estimator on the fit**

In the `result` list, immediately after `inference = inference,` (line 335), add:

```r
    method = method,
```

Two fits of the same model under different estimators must not be indistinguishable in the object they return.

- [ ] **Step 10: Run the validation tests**

```powershell
& $RS -e "devtools::load_all('.'); testthat::test_file('tests/testthat/test-laplace.R')"
```

Expected: all four tests PASS, including the Task 1 SML baseline.

- [ ] **Step 11: Run the full suite**

```powershell
& $RS -e "devtools::load_all('.'); testthat::test_dir('tests/testthat')"
```

Expected: all pass.

- [ ] **Step 12: Commit**

```bash
git add R/fit_rpbnb_tmb.R tests/testthat/test-laplace.R
git commit -m "feat: add method = laplace to fit_rpbnb_tmb"
```

---

### Task 4: Verify the sparse path engages and agrees with SML

This task produces the evidence that decides whether the feature is usable. Step 5 is the gate described in the spec's risk section.

**Files:**
- Modify: `tests/testthat/test-laplace.R`

**Interfaces:**
- Consumes: `fit_rpbnb_tmb(method = "laplace")` from Task 3.
- Produces: no code interface; produces the correctness evidence later tasks depend on.

- [ ] **Step 1: Write the sparse-mechanism test**

Append to `tests/testthat/test-laplace.R`:

```r
test_that("laplace engages TMB's sparse random-effect machinery", {
  skip_on_cran()
  sim <- simulate_rpbnb_tmb(
    n = 200,
    beta1 = c("(Intercept)" = 0.2, x1 = 0.4),
    beta2 = c("(Intercept)" = 0.1, x1 = -0.3),
    random_1 = list(x1 = list(sd = 0.3)),
    dispersion = c(m1 = 0.4, m2 = 0.5),
    dependence = "famoye", lambda = 0, seed = 11
  )
  fit <- fit_rpbnb_tmb(y1 ~ x1, y2 ~ x1, data = sim$data,
                       random_1 = "x1", dependence = "famoye",
                       draws = 50, method = "laplace", keep = "full")
  # One latent per observation per random coefficient. If this is 0 the fit
  # silently ran a fixed-effect model at u = 0 instead of integrating.
  expect_identical(length(fit$obj$env$random), 200L)
  expect_true(is.finite(fit$logLik))
})
```

- [ ] **Step 2: Run it**

```powershell
& $RS -e "devtools::load_all('.'); testthat::test_file('tests/testthat/test-laplace.R')"
```

Expected: PASS. A `length()` of `0` means `random_names` never reached `MakeADFun`; revisit Task 3 Step 8.

- [ ] **Step 3: Write the SML-versus-Laplace agreement test**

Append to `tests/testthat/test-laplace.R`:

```r
test_that("laplace and sml agree on data simulated from the true process", {
  skip_on_cran()
  sim <- simulate_rpbnb_tmb(
    n = 500,
    beta1 = c("(Intercept)" = 0.5, x1 = 0.3),
    beta2 = c("(Intercept)" = 0.2, x1 = -0.2),
    random_1 = list(x1 = list(sd = 0.4)),
    dispersion = c(m1 = 0.4, m2 = 0.5),
    dependence = "famoye", lambda = 0.2, seed = 21
  )
  common <- list(formula_1 = y1 ~ x1, formula_2 = y2 ~ x1,
                 data = sim$data, random_1 = "x1",
                 dependence = "famoye", draws = 200, seed = 99)

  fit_sml <- do.call(fit_rpbnb_tmb, c(common, list(method = "sml")))
  fit_lap <- do.call(fit_rpbnb_tmb, c(common, list(method = "laplace")))

  expect_true(is.finite(fit_lap$logLik))
  expect_true(isTRUE(fit_lap$sdreport$pdHess))

  # Same parameter vector, same names, so a direct comparison is meaningful.
  expect_identical(names(coef(fit_lap)), names(coef(fit_sml)))

  # Agreement judged against sampling noise rather than an absolute tolerance:
  # these are two approximations to the same integral, not two computations of
  # the same number.
  key <- c("b1:x1", "b2:x1")
  for (nm in key) {
    tol <- 3 * max(fit_sml$se[[nm]], fit_lap$se[[nm]])
    expect_lt(abs(coef(fit_sml)[[nm]] - coef(fit_lap)[[nm]]), tol)
  }
})
```

- [ ] **Step 4: Run it**

```powershell
& $RS -e "devtools::load_all('.'); testthat::test_file('tests/testthat/test-laplace.R')"
```

Expected: PASS.

- [ ] **Step 5: Write the multithreaded-Laplace test**

The spec flags `parallel_accumulator` combined with `random=` as supported by TMB but unverified for this template. Every test so far runs single-threaded, and the Task 6 acceptance script runs on 12 cores — so without this step, the first exercise of the combination would be the largest, slowest, hardest-to-diagnose fit in the plan. Append to `tests/testthat/test-laplace.R`:

```r
test_that("laplace gives the same answer single- and multi-threaded", {
  skip_on_cran()
  sim <- simulate_rpbnb_tmb(
    n = 300,
    beta1 = c("(Intercept)" = 0.3, x1 = 0.35),
    beta2 = c("(Intercept)" = 0.1, x1 = -0.25),
    random_1 = list(x1 = list(sd = 0.35)),
    dispersion = c(m1 = 0.4, m2 = 0.5),
    dependence = "famoye", lambda = 0.1, seed = 31
  )
  common <- list(formula_1 = y1 ~ x1, formula_2 = y2 ~ x1,
                 data = sim$data, random_1 = "x1",
                 dependence = "famoye", draws = 50, seed = 5,
                 method = "laplace")

  fit1 <- do.call(fit_rpbnb_tmb, c(common, list(
    control = rpbnb_tmb_control(n_cores = 1L))))
  fit4 <- do.call(fit_rpbnb_tmb, c(common, list(
    control = rpbnb_tmb_control(n_cores = 4L, max_threads = 4L))))

  # Threading partitions the accumulation; it must not change the objective.
  expect_equal(fit1$logLik, fit4$logLik, tolerance = 1e-6)
  expect_equal(unname(coef(fit1)), unname(coef(fit4)), tolerance = 1e-4)
})
```

- [ ] **Step 6: Run it**

```powershell
& $RS -e "devtools::load_all('.'); testthat::test_file('tests/testthat/test-laplace.R')"
```

Expected: PASS. If the multithreaded fit produces a different objective, non-finite gradients, or fails to converge, the cause is `parallel_accumulator` interacting with `random=`. The spec's named fallback is to use a plain accumulator in the Laplace branch: in `src/rpbnb.tmb.cpp`, replace `parallel_accumulator<Type> nll(this);` with a conditional that uses `Type nll = 0;` when `est_method == 1`. The memory goal does not depend on threading, so taking the fallback is an acceptable outcome — but report it rather than applying it silently, since it costs Laplace users their parallelism.

- [ ] **Step 7: If the agreement test fails, stop and report**

A failure in Step 4 is a finding about the estimator, not a test to loosen. Do **not** widen the tolerance to make it pass. Record the two coefficient vectors, both standard-error vectors, both log-likelihoods, and `fit_lap$optimizer$convergence`, then report before continuing. The spec names this outcome explicitly as something to surface rather than ship quietly.

- [ ] **Step 8: Commit**

```bash
git add tests/testthat/test-laplace.R
git commit -m "test: verify laplace engages sparse path and agrees with sml"
```

---

### Task 5: Document the estimator

**Files:**
- Modify: `R/fit_rpbnb_tmb.R:13` (`@param draws`), and the roxygen block above `fit_rpbnb_tmb` (add `@param method`)
- Modify: `README.md:57-62` (memory section)

**Interfaces:**
- Consumes: the `method` argument from Task 3.
- Produces: regenerated `man/fit_rpbnb_tmb.Rd`.

- [ ] **Step 1: Document `method` and amend `draws`**

In `R/fit_rpbnb_tmb.R`, replace the `@param draws` line (line 13) with:

```r
#' @param draws Number of Halton simulation draws. Under
#'   \code{method = "sml"} this sets the simulation grid the likelihood is
#'   averaged over, and tape size scales with \code{nrow(data) * draws}. Under
#'   \code{method = "laplace"} it does not affect the likelihood or the tape,
#'   but still sizes the Halton grid used for the frozen Famoye lambda bounds
#'   and for the post-estimation averaging in \code{predict()} and the
#'   marginal-effect functions.
```

and add, immediately after the `@param poisson_2` block (roxygen `@param` order is independent of signature order, but keeping them aligned aids review):

```r
#' @param method Estimator for the random-coefficient integral. \code{"sml"}
#'   (default) uses simulated maximum likelihood over Halton draws.
#'   \code{"laplace"} uses TMB's Laplace approximation with a sparse Hessian
#'   over one latent vector per observation, which removes \code{draws} from
#'   the memory cost and makes tape size scale with \code{nrow(data)} alone.
#'
#'   The two are different approximations to the same integral. They agree
#'   asymptotically but need not agree closely on any given dataset, so a
#'   Laplace fit is not a drop-in reproduction of an SML fit. \code{"laplace"}
#'   supports \code{"normal"} and \code{"lognormal"} random coefficients only,
#'   and requires at least one random coefficient.
```

- [ ] **Step 2: Regenerate documentation**

```powershell
& $RS -e "devtools::document()"
```

Expected: `man/fit_rpbnb_tmb.Rd` updated, no errors.

- [ ] **Step 3: Amend the README memory section**

In `README.md`, after the paragraph ending "about 8 GiB of peak memory for one fit." (line 62), insert:

```markdown
All of the above bounds the *simulated* likelihood, whose tape scales with
`nrow(data) * draws`. `method = "laplace"` integrates the random coefficients
with TMB's Laplace approximation instead, taping one conditional evaluation per
observation and integrating the latents through a sparse Hessian. Tape size then
scales with `nrow(data)` alone and the draw budget stops binding. That is the
option to reach for when a fit exhausts memory during `MakeADFun()`, rather than
raising `max_workload` against RAM you do not have. It supports normal and
lognormal random coefficients, and it is a different approximation to the same
integral — see `?fit_rpbnb_tmb`.
```

- [ ] **Step 4: Confirm the package still loads and documents cleanly**

```powershell
& $RS -e "devtools::load_all('.'); testthat::test_dir('tests/testthat')"
```

Expected: all pass.

- [ ] **Step 5: Commit**

```bash
git add R/fit_rpbnb_tmb.R man/ README.md
git commit -m "docs: document method = laplace and its effect on draws"
```

---

### Task 6: Acceptance — fit the truck workload that currently exhausts memory

**Files:**
- Create: `inst/truck_rpbnb_diff_famoye_laplace.R`

**Interfaces:**
- Consumes: everything above.
- Produces: the acceptance evidence for the feature.

- [ ] **Step 1: Create the acceptance script**

Create `inst/truck_rpbnb_diff_famoye_laplace.R`. It mirrors `inst/truck_rpbnb_diff_famoye_dense.R` with `method = "laplace"`, the default `max_workload` restored (the override exists only because the SML tape needed it), and reporting of peak memory:

```r
#!/usr/bin/env Rscript
# =============================================================================
# rpbnb.tmb -- truck crash RP-BNB, Famoye dependence, Laplace estimator.
#
# Same model and data as truck_rpbnb_diff_famoye_dense.R, which exhausts memory
# during MakeADFun() because the SML tape scales with n * draws. Laplace tapes
# one conditional evaluation per observation instead, so the draw budget no
# longer binds and the default max_workload is sufficient.
#
# Run from the package root:
#   library(rpbnb.tmb); source("inst/truck_rpbnb_diff_famoye_laplace.R")
# =============================================================================

library(here)
devtools::load_all()

sep <- function() cat("\n", paste(rep("=", 72), collapse = ""), "\n", sep = "")

draws   <- 300L   # post-estimation averaging only under Laplace
n_cores <- 12L

data <- read.csv(here("inst\\extdata", "export_dense_all.csv"))

cat("=== RP-BNB on truck all crashes (Laplace) ===\n")
cat("Observations :", nrow(data), "\n")
cat("Cores asked  :", n_cores, "\n")

f1 <- ALL_3  ~ LNAADT_3+SR40_MI3+MPD_ME+MPD_STD+IRI_ME+RUT_L+SP50GE+ACCPNTS+SIGNAL1+NEAR_SIG+AUXLNUM++DP10_ME
f2 <- C_DISTR ~ LNAADT_3+SR40_MI3+MPD_ME+MPD_STD+RUT_9+ACCPNTS+SIGNAL1+NEAR_SIG+AUXLNUM+DP01_ME+CS_MINAB

t_fit <- system.time(
  fit <- fit_rpbnb_tmb(
    formula_1  = f1,
    formula_2  = f2,
    data       = data,
    random_1   = c("SR40_MI3", "AUXLNUM", "MPD_ME", "MPD_STD"),
    random_2   = c("SR40_MI3", "RUT_9", "MPD_ME", "MPD_STD"),
    dependence = "famoye",
    draws      = draws,
    seed       = 20240712,
    method     = "laplace",
    control    = rpbnb_tmb_control(
      print_level = 1,
      n_cores     = n_cores,
      max_threads = n_cores
    )
  )
)[["elapsed"]]

cat(sprintf("\nEstimation finished in %.2f s\n", t_fit))
# Deliberately not reporting gc() figures: the TMB tape lives on the C++ heap
# and is invisible to R's garbage collector, so gc() would understate exactly
# the quantity this script exists to test. Watch the process working set
# externally (Task Manager, or inst/benchmark_memory.R) if a number is needed.
cat(sprintf("Optimizer: code=%d, message=%s\n",
            fit$optimizer$convergence, fit$optimizer$message))
cat("sdreport positive-definite Hessian:",
    if (isTRUE(fit$sdreport$pdHess)) "yes" else "no", "\n")

sep(); cat("MODEL SUMMARY\n"); sep()
summary(fit)

sep(); cat("DEPENDENCE (sdreport)\n"); sep()
if (!is.null(fit$sdreport)) print(summary(fit$sdreport, "report"))
```

- [ ] **Step 2: Run the acceptance script**

```powershell
& $RS -e "source('inst/truck_rpbnb_diff_famoye_laplace.R')"
```

Expected: the fit completes without `std::bad_alloc`, prints an optimizer code and a `pdHess` line, and prints a summary. This is the fit that fails today.

- [ ] **Step 3: Record the outcome**

Note the elapsed time, the peak working set observed externally, the convergence code, and whether `pdHess` is `yes`. If the fit completes but `pdHess` is `no` or the optimizer reports non-convergence, that is a real result to report — the memory goal is met but the fit is not usable, and that distinction must reach the user rather than being buried under "it ran".

- [ ] **Step 4: Run the full suite one final time**

```powershell
& $RS -e "devtools::load_all('.'); testthat::test_dir('tests/testthat')"
```

Expected: all pass.

- [ ] **Step 5: Commit**

```bash
git add inst/truck_rpbnb_diff_famoye_laplace.R
git commit -m "example: truck Famoye fit via Laplace estimator"
```

---

## Spec Coverage

| Spec section | Task |
|---|---|
| C++ new inputs (`est_method`, `u1`, `u2`) | 1 |
| C++ branch structure (dev guard, latent deviations, reduction) | 2 |
| `parallel_accumulator` retained; fallback if it misbehaves | 2 (no change to the accumulator); verified by Task 4 Steps 5-6, which also spell out the fallback |
| R public interface (`method`, validation) | 3 |
| Object construction (parameters, map, dummy Z, `random=`) | 3 |
| `draws` dual role under Laplace | 3 (Step 6 comment, test in Step 1), 5 (roxygen) |
| Workload guard with `draws = 1L` | 3 Step 5 |
| Inference/methods unchanged | 3 Step 9 (`method` recorded); no other change needed |
| SML regression | 1 |
| Laplace mechanism | 4 |
| Laplace correctness | 4 |
| Error paths | 3 |
| Acceptance (truck fit) | 6 |
| Documentation | 5 |
