# Simple Famoye Example Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the Famoye example's helper-based configuration with three editable settings and use every detected physical core by default.

**Architecture:** Keep the example as one linear sourceable script. Static tests verify the settings and data flow; a reduced temporary copy verifies real execution without running the high-memory defaults.

**Tech Stack:** R 4.5, testthat, TMB, base R.

## Global Constraints

- Set `n_obs <- 5000L`.
- Set `draws <- 400L`.
- Set `n_cores <- parallel::detectCores(logical = FALSE)` with an `NA` fallback of one.
- Remove all `.famoye_example_*` helpers and all `RPBNB_FAMOYE_*` environment-variable handling.
- Preserve formulas, Famoye dependence, seed, summaries, marginal effects, and elasticities.
- Preserve unrelated `.gitignore` and `.Rbuildignore` changes.

---

### Task 1: Simplify and verify the Famoye example

**Files:**
- Modify: `inst/fit_rpbnb_diff_famoye.R`
- Modify: `tests/testthat/test-example-famoye.R`
- Modify: `tests/testthat/test-parallel.R`

**Interfaces:**
- Produces: editable `n_obs`, `draws`, and `n_cores` script variables.
- Consumes: `parallel::detectCores(logical = FALSE)` and the available data row count.

- [ ] **Step 1: Replace helper tests with failing simple-script tests**

Replace `tests/testthat/test-example-famoye.R` with:

```r
famoye_example_text <- function() {
  path <- testthat::test_path(
    "..", "..", "inst", "fit_rpbnb_diff_famoye.R"
  )
  paste(readLines(path, warn = FALSE), collapse = "\n")
}

test_that("Famoye example has three editable workload settings", {
  script <- famoye_example_text()

  expect_match(script, "n_obs <- 5000L", fixed = TRUE)
  expect_match(script, "draws <- 400L", fixed = TRUE)
  expect_match(
    script,
    "n_cores <- parallel::detectCores(logical = FALSE)",
    fixed = TRUE
  )
  expect_match(script, "if (is.na(n_cores)) n_cores <- 1L", fixed = TRUE)
})

test_that("Famoye example uses settings directly", {
  script <- famoye_example_text()

  expect_match(script, "min(n_obs, nrow(data))", fixed = TRUE)
  expect_match(script, "draws      = draws", fixed = TRUE)
  expect_match(script, "n_cores     = n_cores", fixed = TRUE)
  expect_match(script, "fit$parallel$realized", fixed = TRUE)
  expect_match(script, "fit$optimizer$convergence", fixed = TRUE)
  expect_match(script, "fit$sdreport$pdHess", fixed = TRUE)
})

test_that("Famoye example has no configuration helpers", {
  script <- famoye_example_text()

  expect_false(grepl(".famoye_example_", script, fixed = TRUE))
  expect_false(grepl("RPBNB_FAMOYE_", script, fixed = TRUE))
})
```

Update the Famoye assertions in `tests/testthat/test-parallel.R` to require:

```r
expect_match(
  famoye_demo,
  "n_cores <- parallel::detectCores(logical = FALSE)",
  fixed = TRUE
)
expect_match(famoye_demo, "n_cores     = n_cores", fixed = TRUE)
expect_false(grepl("min(2L", famoye_demo, fixed = TRUE))
```

- [ ] **Step 2: Run tests and verify RED**

Run:

```powershell
& 'C:\Program Files\R\R-4.5.1\bin\Rscript.exe' -e ".libPaths(c('C:/Users/litabook/repos/rpbnb_tmb/.worktrees/.rlib', .libPaths())); testthat::test_file('tests/testthat/test-example-famoye.R', reporter='summary')"
```

Expected: failures because the script still contains Famoye helper functions
and environment-variable configuration.

- [ ] **Step 3: Simplify the script**

After `library(rpbnb.tmb)`, use:

```r
n_obs <- 5000L
draws <- 400L
n_cores <- parallel::detectCores(logical = FALSE)
if (is.na(n_cores)) n_cores <- 1L
```

Select the available rows and print settings:

```r
data <- read.csv(file.path("data", "rwm1984_bnb.csv"))
data <- data[seq_len(min(n_obs, nrow(data))), , drop = FALSE]
cat("=== RP-BNB (different formulas, Famoye) on rwm1984 health counts ===\n")
cat("Observations :", nrow(data), "\n")
cat("Draws        :", draws, "\n")
cat("Cores asked  :", n_cores, "\n")
```

Use `draws` and `n_cores` directly in `fit_rpbnb_tmb()`. After the fit, print:

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

- [ ] **Step 4: Run focused and regression tests**

Run:

```powershell
$env:NOT_CRAN='true'
& 'C:\Program Files\R\R-4.5.1\bin\Rscript.exe' -e ".libPaths(c('C:/Users/litabook/repos/rpbnb_tmb/.worktrees/.rlib', .libPaths())); testthat::test_dir('tests/testthat', filter='(example-famoye|parallel)', reporter='summary', package='rpbnb.tmb', load_package='installed')"
```

Expected: all selected tests pass.

- [ ] **Step 5: Verify a reduced real run**

Create a temporary copy in R, replace only the three setting lines with
`120L`, `20L`, and `2L`, then source that temporary copy:

```powershell
& 'C:\Program Files\R\R-4.5.1\bin\Rscript.exe' -e ".libPaths(c('C:/Users/litabook/repos/rpbnb_tmb/.worktrees/.rlib', .libPaths())); x <- readLines('inst/fit_rpbnb_diff_famoye.R'); x <- sub('^n_obs <- 5000L$', 'n_obs <- 120L', x); x <- sub('^draws <- 400L$', 'draws <- 20L', x); x <- sub('^n_cores <- parallel::detectCores\\\\(logical = FALSE\\\\)$', 'n_cores <- 2L', x); p <- tempfile(fileext='.R'); writeLines(x, p); source(p)"
```

Expected: 120 observations, 20 draws, two requested and realized threads,
optimizer code zero, and no error.

- [ ] **Step 6: Review and commit**

Run:

```powershell
git diff --check
git add -- inst/fit_rpbnb_diff_famoye.R tests/testthat/test-example-famoye.R tests/testthat/test-parallel.R
git commit -m "refactor: simplify Famoye example settings"
```

Expected: the simplified script and its tests are committed; unrelated local
changes remain unstaged.
