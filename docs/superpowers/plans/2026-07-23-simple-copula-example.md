# Simple Copula Example Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Simplify the Frank-copula example so it follows the Famoye example's plain editable-settings structure.

**Architecture:** Replace the environment-variable parsing helpers with three top-level workload variables and use them directly in row selection and the fit call. Keep the model and report sections unchanged, including the user's currently commented-out elasticity section, and retain concise diagnostics inline after fitting.

**Tech Stack:** R 4.5, `rpbnb.tmb`, TMB/OpenMP, `testthat`, Git

## Global Constraints

- Default to `n_obs <- 500L` and `draws <- 20L`.
- Default `n_cores` to every detected physical core, with `1L` as the `NA` fallback.
- Do not read `RPBNB_N_OBS`, `RPBNB_DRAWS`, or `RPBNB_N_CORES`.
- Preserve the formulas, random coefficients, Frank copula, seed, summaries, dependence report, and marginal effects.
- Preserve the user's commented-out copula elasticity section.
- Do not modify or stage the user's changes in `inst/fit_rpbnb_diff_famoye.R`, `.gitignore`, or `.Rbuildignore`.

## File Structure

- Modify `tests/testthat/test-example-copula.R`: replace helper-unit tests with static tests for the simple script contract.
- Modify `tests/testthat/test-parallel.R`: update copula example assertions to expect all physical cores and direct `n_cores` use.
- Modify `inst/fit_rpbnb_diff_copula.R`: remove configuration helpers and use plain settings and inline diagnostics.

---

### Task 1: Simplify the Copula Example

**Files:**

- Modify: `tests/testthat/test-example-copula.R`
- Modify: `tests/testthat/test-parallel.R:255-285`
- Modify: `inst/fit_rpbnb_diff_copula.R:1-130`

**Interfaces:**

- Consumes: `parallel::detectCores(logical = FALSE)`, `fit_rpbnb_tmb()`, `rpbnb_tmb_control()`
- Produces: editable `n_obs`, `draws`, and `n_cores` script variables used directly by the example

- [ ] **Step 1: Replace the copula helper tests with failing contract tests**

Replace `tests/testthat/test-example-copula.R` with:

```r
copula_example_text <- function() {
  path <- testthat::test_path(
    "..", "..", "inst", "fit_rpbnb_diff_copula.R"
  )
  paste(readLines(path, warn = FALSE), collapse = "\n")
}

test_that("copula example has three editable workload settings", {
  script <- copula_example_text()

  expect_match(script, "n_obs <- 500L", fixed = TRUE)
  expect_match(script, "draws <- 20L", fixed = TRUE)
  expect_match(
    script,
    "n_cores <- parallel::detectCores(logical = FALSE)",
    fixed = TRUE
  )
  expect_match(
    script,
    "if (is.na(n_cores)) n_cores <- 1L",
    fixed = TRUE
  )
})

test_that("copula example uses settings directly", {
  script <- copula_example_text()

  expect_match(script, "min(n_obs, nrow(data))", fixed = TRUE)
  expect_match(script, "draws      = draws", fixed = TRUE)
  expect_match(script, "n_cores     = n_cores", fixed = TRUE)
  expect_match(script, "fit$parallel$realized", fixed = TRUE)
  expect_match(script, "fit$optimizer$convergence", fixed = TRUE)
  expect_match(script, "fit$sdreport$pdHess", fixed = TRUE)
})

test_that("copula example has no configuration helpers", {
  script <- copula_example_text()

  expect_false(grepl(".example_", script, fixed = TRUE))
  expect_false(grepl("RPBNB_", script, fixed = TRUE))
})
```

In `tests/testthat/test-parallel.R`, replace the four copula assertions in
`"demo scripts use physical-core counts"` with:

```r
expect_match(
  copula_demo,
  "n_cores <- parallel::detectCores(logical = FALSE)",
  fixed = TRUE
)
expect_match(
  copula_demo,
  "n_cores     = n_cores",
  fixed = TRUE
)
expect_false(grepl("min(4L", copula_demo, fixed = TRUE))
```

- [ ] **Step 2: Run the focused tests and verify the RED state**

Run:

```powershell
$env:NOT_CRAN='true'
& 'C:\Program Files\R\R-4.5.1\bin\Rscript.exe' -e ".libPaths(c('C:/Users/litabook/repos/rpbnb_tmb/.worktrees/.rlib', .libPaths())); testthat::test_dir('tests/testthat', filter='(example-copula|parallel)', reporter='summary', package='rpbnb.tmb', load_package='installed')"
```

Expected: `test-example-copula.R` and the copula assertions in
`test-parallel.R` fail because the script still has `example_*` variables,
environment-variable helpers, and a four-core cap.

- [ ] **Step 3: Remove the configuration layer and add plain settings**

In `inst/fit_rpbnb_diff_copula.R`, remove the optional `Sys.setenv()` comment,
the `.example_positive_integer`, `.example_observation_count`, and
`.example_fit_diagnostics` functions, and the `detected_cores`,
`default_cores`, `example_n_obs`, `example_draws`, and `example_cores`
assignments.

Immediately after `sep`, add:

```r
# ---- Settings ---------------------------------------------------------------
n_obs <- 500L
draws <- 20L
n_cores <- parallel::detectCores(logical = FALSE)
if (is.na(n_cores)) n_cores <- 1L
```

Replace data selection and workload output with:

```r
data <- read.csv(file.path("data", "rwm1984_bnb.csv"))
data <- data[seq_len(min(n_obs, nrow(data))), , drop = FALSE]
cat("=== RP-BNB (different formulas, Frank copula) on rwm1984 health counts ===\n")
cat("Observations :", nrow(data), "\n")
cat("Draws        :", draws, "\n")
cat("Cores asked  :", n_cores, "\n")
```

Use the settings directly in the fit:

```r
    draws      = draws,
    seed       = 20240712,
    control    = rpbnb_tmb_control(
      print_level = 1,
      n_cores     = n_cores
    )
```

Replace `.example_fit_diagnostics(fit, t_fit)` with:

```r
cat(sprintf("\nEstimation finished in %.2f s\n", t_fit))
cat(sprintf(
  "TMB threads: requested=%d, realized=%d\n",
  fit$parallel$requested, fit$parallel$realized
))
cat(sprintf(
  "Optimizer: code=%d, message=%s\n",
  fit$optimizer$convergence, fit$optimizer$message
))
cat(
  "sdreport positive-definite Hessian:",
  if (isTRUE(fit$sdreport$pdHess)) "yes" else "no",
  "\n"
)
```

Do not change the commented lines:

```r
#sep(); cat("ELASTICITIES / SEMI-ELASTICITIES (AME)\n"); sep()
#rpbnb_tmb_elasticities(fit, which = "both")
```

- [ ] **Step 4: Run the focused tests and verify the GREEN state**

Run:

```powershell
$env:NOT_CRAN='true'
& 'C:\Program Files\R\R-4.5.1\bin\Rscript.exe' -e ".libPaths(c('C:/Users/litabook/repos/rpbnb_tmb/.worktrees/.rlib', .libPaths())); invisible(parse(file='inst/fit_rpbnb_diff_copula.R')); testthat::test_dir('tests/testthat', filter='(example-copula|parallel)', reporter='summary', package='rpbnb.tmb', load_package='installed')"
```

Expected: the script parses and both focused test files pass.

- [ ] **Step 5: Review and commit only the simplification**

Run:

```powershell
git diff --check
git diff -- inst/fit_rpbnb_diff_copula.R tests/testthat/test-example-copula.R tests/testthat/test-parallel.R
git add -p -- inst/fit_rpbnb_diff_copula.R
git add -- tests/testthat/test-example-copula.R tests/testthat/test-parallel.R
git diff --cached --check
git diff --cached
```

At the `git add -p` prompts, stage the settings, data, fit, and diagnostics
hunks. Do not stage the pre-existing elasticity-comment hunk. Verify the
cached diff contains the simplification but leaves the elasticity section
unchanged from `HEAD`.

Commit:

```powershell
git commit -m "refactor: simplify copula example settings"
```

Expected: the commit contains only the copula simplification and its tests.
The user's elasticity, Famoye, `.gitignore`, and `.Rbuildignore` changes remain
unstaged.

---

### Task 2: Verify a Real Copula Fit

**Files:**

- Verify: `inst/fit_rpbnb_diff_copula.R`

**Interfaces:**

- Consumes: the simplified source script from Task 1
- Produces: evidence that direct settings, OpenMP threads, optimization, and Hessian reporting work together

- [ ] **Step 1: Run a reduced temporary copy**

Run:

```powershell
& 'C:\Program Files\R\R-4.5.1\bin\Rscript.exe' -e ".libPaths(c('C:/Users/litabook/repos/rpbnb_tmb/.worktrees/.rlib', .libPaths())); x <- readLines('inst/fit_rpbnb_diff_copula.R'); x <- sub('^n_obs <- 500L$', 'n_obs <- 120L', x); x <- sub('^n_cores <- parallel::detectCores\\(logical = FALSE\\)$', 'n_cores <- 2L', x); x <- sub('print_level = 1', 'print_level = 0', x, fixed=TRUE); p <- tempfile(fileext='.R'); writeLines(x, p); e <- new.env(parent=globalenv()); out <- capture.output(source(p, local=e)); cat(grep('^(Observations|Draws|Cores asked|Estimation finished|TMB threads|Optimizer|sdreport)', out, value=TRUE), sep='\n'); fit_check <- e[['fit']]; stopifnot(fit_check[['parallel']][['requested']] == 2L, fit_check[['parallel']][['realized']] == 2L, fit_check[['optimizer']][['convergence']] == 0L, isTRUE(fit_check[['sdreport']][['pdHess']]))"
```

Expected output includes:

```text
Observations : 120
Draws        : 20
Cores asked  : 2
TMB threads: requested=2, realized=2
Optimizer: code=0
sdreport positive-definite Hessian: yes
```

- [ ] **Step 2: Confirm repository state**

Run:

```powershell
git status --short
git log -1 --oneline
git diff --check
```

Expected: the simplification commit is at `HEAD`; only the user's preserved
copula elasticity, Famoye, `.gitignore`, and `.Rbuildignore` changes remain.
