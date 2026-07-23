# Configurable Parallel Copula Example Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the Frank-copula example configurable through environment variables and report its actual parallel execution and convergence status.

**Architecture:** Keep the example as one sourceable script. Add three small, side-effect-light helper functions at the top of the script for positive-integer parsing, observation-count capping, and post-fit diagnostics; the main script consumes those helpers before and after the existing fit. Tests parse and evaluate only those helper assignments in an isolated environment, so verification does not run the expensive model.

**Tech Stack:** R 4.5, testthat, TMB, base R environment variables.

## Global Constraints

- Preserve the current Frank-copula formulas, random coefficients, seed, summary, fitted means, dependence report, marginal effects, and elasticities.
- Support `RPBNB_N_OBS` with default `500`, `RPBNB_DRAWS` with default `20`, and `RPBNB_N_CORES` with a default capped at four detected physical cores.
- Accept only one finite, positive, whole-number value no larger than `.Machine$integer.max`.
- Cap the requested observation count at the available row count and explain the cap.
- Keep the script runnable with `source("inst/fit_rpbnb_diff_copula.R")`.
- Preserve unrelated unstaged changes in the working tree.

---

## File Structure

- Modify `inst/fit_rpbnb_diff_copula.R`: configuration helpers, environment-variable consumption, configuration output, and fit diagnostics.
- Create `tests/testthat/test-example-copula.R`: isolated behavioral tests for the example helpers.

### Task 1: Configuration parsing and observation capping

**Files:**
- Modify: `inst/fit_rpbnb_diff_copula.R`
- Create: `tests/testthat/test-example-copula.R`

**Interfaces:**
- Produces: `.example_positive_integer(name, default)` returning one positive integer.
- Produces: `.example_observation_count(requested, available)` returning the usable integer row count.
- Consumes: `Sys.getenv()`, `.Machine$integer.max`, and the loaded data row count.

- [ ] **Step 1: Write the failing helper-loading and configuration tests**

```r
load_copula_example_helpers <- function() {
  path <- testthat::test_path("..", "..", "inst",
                             "fit_rpbnb_diff_copula.R")
  expressions <- parse(path)
  helper_names <- c(
    ".example_positive_integer",
    ".example_observation_count",
    ".example_fit_diagnostics"
  )
  helper_env <- new.env(parent = baseenv())

  for (expression in expressions) {
    if (is.call(expression) &&
        identical(expression[[1L]], as.name("<-")) &&
        as.character(expression[[2L]]) %in% helper_names) {
      eval(expression, envir = helper_env)
    }
  }

  helper_env
}

test_that("copula example reads positive integer configuration", {
  helpers <- load_copula_example_helpers()
  old <- Sys.getenv("RPBNB_DRAWS", unset = NA_character_)
  on.exit(if (is.na(old)) Sys.unsetenv("RPBNB_DRAWS") else
    Sys.setenv(RPBNB_DRAWS = old), add = TRUE)

  Sys.unsetenv("RPBNB_DRAWS")
  expect_identical(
    helpers$.example_positive_integer("RPBNB_DRAWS", 20L),
    20L
  )

  Sys.setenv(RPBNB_DRAWS = "40")
  expect_identical(
    helpers$.example_positive_integer("RPBNB_DRAWS", 20L),
    40L
  )
})

test_that("copula example rejects invalid integer configuration", {
  helpers <- load_copula_example_helpers()
  old <- Sys.getenv("RPBNB_N_CORES", unset = NA_character_)
  on.exit(if (is.na(old)) Sys.unsetenv("RPBNB_N_CORES") else
    Sys.setenv(RPBNB_N_CORES = old), add = TRUE)

  for (value in c("abc", "1.5", "0", "-1", "Inf", "1,2")) {
    Sys.setenv(RPBNB_N_CORES = value)
    expect_error(
      helpers$.example_positive_integer("RPBNB_N_CORES", 1L),
      "RPBNB_N_CORES"
    )
  }
})

test_that("copula example caps observations to available rows", {
  helpers <- load_copula_example_helpers()

  expect_identical(helpers$.example_observation_count(50L, 100L), 50L)
  expect_message(
    actual <- helpers$.example_observation_count(150L, 100L),
    "150.*100"
  )
  expect_identical(actual, 100L)
})
```

- [ ] **Step 2: Run the focused test and verify RED**

Run:

```powershell
& 'C:\Program Files\R\R-4.5.1\bin\Rscript.exe' -e ".libPaths(c('C:/Users/litabook/repos/rpbnb_tmb/.worktrees/.rlib', .libPaths())); testthat::test_file('tests/testthat/test-example-copula.R', reporter='summary')"
```

Expected: FAIL because `.example_positive_integer()` and
`.example_observation_count()` are not present in the script.

- [ ] **Step 3: Add the minimal configuration helpers**

Add after `library(rpbnb.tmb)`:

```r
.example_positive_integer <- function(name, default) {
  raw <- Sys.getenv(name, unset = "")
  if (!nzchar(raw)) return(as.integer(default))

  value <- suppressWarnings(as.numeric(raw))
  valid <- length(value) == 1L &&
    is.finite(value) &&
    value >= 1 &&
    value == floor(value) &&
    value <= .Machine$integer.max
  if (!valid) {
    stop(name, " must be a single positive integer.", call. = FALSE)
  }
  as.integer(value)
}

.example_observation_count <- function(requested, available) {
  if (requested > available) {
    message(
      "RPBNB_N_OBS requested ", requested,
      " rows; using all ", available, " available rows."
    )
  }
  as.integer(min(requested, available))
}
```

Replace fixed demo values with:

```r
detected_cores <- parallel::detectCores(logical = FALSE)
default_cores <- if (is.na(detected_cores)) {
  1L
} else {
  max(1L, min(4L, as.integer(detected_cores)))
}
example_n_obs <- .example_positive_integer("RPBNB_N_OBS", 500L)
example_draws <- .example_positive_integer("RPBNB_DRAWS", 20L)
example_cores <- .example_positive_integer("RPBNB_N_CORES", default_cores)
```

After reading the CSV, select rows with:

```r
example_n_obs <- .example_observation_count(example_n_obs, nrow(data))
data <- data[seq_len(example_n_obs), , drop = FALSE]
cat("Observations :", nrow(data), "\n")
cat("Draws        :", example_draws, "\n")
cat("Cores asked  :", example_cores, "\n")
```

Pass `draws = example_draws` and `n_cores = example_cores` to
`fit_rpbnb_tmb()`.

- [ ] **Step 4: Run the focused test and verify GREEN**

Run the command from Step 2.

Expected: the three configuration tests pass.

### Task 2: Parallel and convergence diagnostics

**Files:**
- Modify: `inst/fit_rpbnb_diff_copula.R`
- Modify: `tests/testthat/test-example-copula.R`

**Interfaces:**
- Produces: `.example_fit_diagnostics(fit, elapsed)` printing elapsed time,
  requested and realized threads, optimizer status, and Hessian status.
- Consumes: `fit$parallel`, `fit$optimizer`, and `fit$sdreport$pdHess`.

- [ ] **Step 1: Write the failing diagnostics test**

```r
test_that("copula example reports parallel and convergence diagnostics", {
  helpers <- load_copula_example_helpers()
  fit <- list(
    parallel = list(requested = 4L, realized = 2L),
    optimizer = list(convergence = 0L, message = "relative convergence"),
    sdreport = list(pdHess = TRUE)
  )

  output <- capture.output(helpers$.example_fit_diagnostics(fit, 1.25))

  expect_match(output, "1.25 s")
  expect_match(output, "requested=4")
  expect_match(output, "realized=2")
  expect_match(output, "code=0")
  expect_match(output, "relative convergence")
  expect_match(output, "positive-definite Hessian: yes")
})
```

- [ ] **Step 2: Run the focused test and verify RED**

Run the Task 1 Step 2 command.

Expected: FAIL because `.example_fit_diagnostics()` is absent.

- [ ] **Step 3: Add the diagnostics helper and call it**

Add after the configuration helpers:

```r
.example_fit_diagnostics <- function(fit, elapsed) {
  requested <- fit$parallel$requested
  realized <- fit$parallel$realized
  convergence <- fit$optimizer$convergence
  optimizer_message <- fit$optimizer$message
  pd_hessian <- isTRUE(fit$sdreport$pdHess)

  cat(sprintf("\nEstimation finished in %.2f s\n", elapsed))
  cat(sprintf(
    "TMB threads: requested=%d, realized=%d\n",
    requested, realized
  ))
  cat(sprintf(
    "Optimizer: code=%d, message=%s\n",
    convergence, optimizer_message
  ))
  cat(
    "sdreport positive-definite Hessian:",
    if (pd_hessian) "yes" else "no",
    "\n"
  )
}
```

Replace the existing elapsed-time `cat()` call with:

```r
.example_fit_diagnostics(fit, t_fit)
```

- [ ] **Step 4: Run the focused test and verify GREEN**

Run the Task 1 Step 2 command.

Expected: all four example tests pass.

### Task 3: Full verification and focused commit

**Files:**
- Verify: `inst/fit_rpbnb_diff_copula.R`
- Verify: `tests/testthat/test-example-copula.R`

**Interfaces:**
- Consumes all helpers and installed-package test infrastructure.
- Produces a sourceable configurable parallel example with regression evidence.

- [ ] **Step 1: Reinstall the current source package**

Run:

```powershell
& 'C:\Program Files\R\R-4.5.1\bin\R.exe' CMD INSTALL --preclean --library='C:\Users\litabook\repos\rpbnb_tmb\.worktrees\.rlib' .
```

Expected: installation succeeds and compilation includes `-fopenmp`.

- [ ] **Step 2: Run focused, copula, and parallel tests**

Run:

```powershell
$env:NOT_CRAN='true'
& 'C:\Program Files\R\R-4.5.1\bin\Rscript.exe' -e ".libPaths(c('C:/Users/litabook/repos/rpbnb_tmb/.worktrees/.rlib', .libPaths())); testthat::test_dir('tests/testthat', filter='(example-copula|fit-copula|parallel)', reporter='summary', package='rpbnb.tmb', load_package='installed')"
```

Expected: all selected tests pass with no failures.

- [ ] **Step 3: Run the example with a small explicit configuration**

Run:

```powershell
$env:RPBNB_N_OBS='120'
$env:RPBNB_DRAWS='10'
$env:RPBNB_N_CORES='2'
& 'C:\Program Files\R\R-4.5.1\bin\Rscript.exe' -e ".libPaths(c('C:/Users/litabook/repos/rpbnb_tmb/.worktrees/.rlib', .libPaths())); source('inst/fit_rpbnb_diff_copula.R')"
```

Expected: output reports 120 observations, 10 draws, requested and realized
thread counts, optimizer status, and Hessian status; the script completes
without an error.

- [ ] **Step 4: Review the diff and commit only the example update**

Run:

```powershell
git diff --check
git diff -- inst/fit_rpbnb_diff_copula.R tests/testthat/test-example-copula.R
git add -- inst/fit_rpbnb_diff_copula.R tests/testthat/test-example-copula.R
git commit -m "feat: make parallel copula example configurable"
```

Expected: the commit contains only the example script and its test. Existing
unrelated `.gitignore`, Famoye-example, and `.Rbuildignore` changes remain
unstaged.
