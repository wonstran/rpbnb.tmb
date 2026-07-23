# Memory-Aware Famoye Example Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the Famoye example safely configurable and provide actionable guidance when TMB exhausts memory while constructing its AD tape.

**Architecture:** Keep the example sourceable and self-contained. Add four Famoye-prefixed helper functions for integer configuration, row-count capping, fit diagnostics, and allocation-error translation; tests parse only those function assignments into an isolated environment so they do not run the model.

**Tech Stack:** R 4.5, testthat, TMB, base R condition handling and environment variables.

## Global Constraints

- Use `RPBNB_FAMOYE_N_OBS` with default `500`.
- Use `RPBNB_FAMOYE_DRAWS` with default `100`.
- Use `RPBNB_FAMOYE_N_CORES` with a default capped at two detected physical cores.
- Accept only one finite, positive whole number no larger than `.Machine$integer.max`.
- Preserve the existing formulas, random coefficients, Famoye dependence, seed, summary, fitted-mean section, dependence report, marginal effects, and elasticities.
- Translate only errors containing `std::bad_alloc`; preserve unrelated error messages and classes.
- Preserve unrelated unstaged `.gitignore` and `.Rbuildignore` changes.

---

## File Structure

- Modify `inst/fit_rpbnb_diff_famoye.R`: Famoye-specific configuration, safer defaults, diagnostics, and allocation guidance.
- Create `tests/testthat/test-example-famoye.R`: isolated behavioral tests for every new helper.
- Modify `tests/testthat/test-parallel.R`: reflect the Famoye example's two-core default bound and validated override.

### Task 1: Famoye configuration and bounded workload

**Files:**
- Create: `tests/testthat/test-example-famoye.R`
- Modify: `tests/testthat/test-parallel.R`
- Modify: `inst/fit_rpbnb_diff_famoye.R`

**Interfaces:**
- Produces: `.famoye_example_positive_integer(name, default)` returning one positive integer.
- Produces: `.famoye_example_observation_count(requested, available)` returning the usable row count.
- Consumes: the three `RPBNB_FAMOYE_*` variables and detected physical cores.

- [ ] **Step 1: Write failing configuration tests**

Create `tests/testthat/test-example-famoye.R`:

```r
load_famoye_example_helpers <- function() {
  path <- testthat::test_path(
    "..", "..", "inst", "fit_rpbnb_diff_famoye.R"
  )
  expressions <- parse(path)
  helper_names <- c(
    ".famoye_example_positive_integer",
    ".famoye_example_observation_count",
    ".famoye_example_fit_diagnostics",
    ".famoye_example_with_memory_guidance"
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

test_that("Famoye example reads positive integer configuration", {
  helpers <- load_famoye_example_helpers()
  old <- Sys.getenv("RPBNB_FAMOYE_DRAWS", unset = NA_character_)
  on.exit(
    if (is.na(old)) Sys.unsetenv("RPBNB_FAMOYE_DRAWS") else
      Sys.setenv(RPBNB_FAMOYE_DRAWS = old),
    add = TRUE
  )

  Sys.unsetenv("RPBNB_FAMOYE_DRAWS")
  expect_identical(
    helpers$.famoye_example_positive_integer(
      "RPBNB_FAMOYE_DRAWS", 100L
    ),
    100L
  )
  Sys.setenv(RPBNB_FAMOYE_DRAWS = "40")
  expect_identical(
    helpers$.famoye_example_positive_integer(
      "RPBNB_FAMOYE_DRAWS", 100L
    ),
    40L
  )
})

test_that("Famoye example rejects invalid integer configuration", {
  helpers <- load_famoye_example_helpers()
  old <- Sys.getenv("RPBNB_FAMOYE_N_CORES", unset = NA_character_)
  on.exit(
    if (is.na(old)) Sys.unsetenv("RPBNB_FAMOYE_N_CORES") else
      Sys.setenv(RPBNB_FAMOYE_N_CORES = old),
    add = TRUE
  )

  for (value in c("abc", "1.5", "0", "-1", "Inf", "1,2")) {
    Sys.setenv(RPBNB_FAMOYE_N_CORES = value)
    expect_error(
      helpers$.famoye_example_positive_integer(
        "RPBNB_FAMOYE_N_CORES", 1L
      ),
      "RPBNB_FAMOYE_N_CORES"
    )
  }
})

test_that("Famoye example caps observations to available rows", {
  helpers <- load_famoye_example_helpers()
  expect_identical(
    helpers$.famoye_example_observation_count(50L, 100L),
    50L
  )
  expect_message(
    actual <- helpers$.famoye_example_observation_count(150L, 100L),
    "150.*100"
  )
  expect_identical(actual, 100L)
})
```

Replace the demo assertion in `tests/testthat/test-parallel.R` with:

```r
test_that("demo scripts use bounded physical-core counts", {
  demo_paths <- testthat::test_path(
    "..", "..", "inst",
    c("fit_rpbnb_diff_copula.R", "fit_rpbnb_diff_famoye.R")
  )
  copula_demo <- paste(
    readLines(demo_paths[[1L]], warn = FALSE), collapse = "\n"
  )
  famoye_demo <- paste(
    readLines(demo_paths[[2L]], warn = FALSE), collapse = "\n"
  )

  expect_match(copula_demo, "detectCores(logical = FALSE)", fixed = TRUE)
  expect_match(copula_demo, "max(1L, min(4L", fixed = TRUE)
  expect_match(
    copula_demo,
    '.example_positive_integer("RPBNB_N_CORES", default_cores)',
    fixed = TRUE
  )
  expect_match(copula_demo, "n_cores     = example_cores", fixed = TRUE)

  expect_match(famoye_demo, "detectCores(logical = FALSE)", fixed = TRUE)
  expect_match(famoye_demo, "max(1L, min(2L", fixed = TRUE)
  expect_match(
    famoye_demo,
    paste0(
      '.famoye_example_positive_integer(',
      '"RPBNB_FAMOYE_N_CORES", famoye_default_cores)'
    ),
    fixed = TRUE
  )
  expect_match(
    famoye_demo,
    "n_cores     = famoye_example_cores",
    fixed = TRUE
  )
})
```

- [ ] **Step 2: Run tests and verify RED**

Run:

```powershell
& 'C:\Program Files\R\R-4.5.1\bin\Rscript.exe' -e ".libPaths(c('C:/Users/litabook/repos/rpbnb_tmb/.worktrees/.rlib', .libPaths())); testthat::test_file('tests/testthat/test-example-famoye.R', reporter='summary'); testthat::test_file('tests/testthat/test-parallel.R', reporter='summary')"
```

Expected: failures because the Famoye helpers and new bounded configuration
are absent.

- [ ] **Step 3: Add configuration helpers and consume them**

Add after `library(rpbnb.tmb)`:

```r
.famoye_example_positive_integer <- function(name, default) {
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

.famoye_example_observation_count <- function(requested, available) {
  if (requested > available) {
    message(
      "RPBNB_FAMOYE_N_OBS requested ", requested,
      " rows; using all ", available, " available rows."
    )
  }
  as.integer(min(requested, available))
}
```

Replace the fixed settings with:

```r
detected_cores <- parallel::detectCores(logical = FALSE)
famoye_default_cores <- if (is.na(detected_cores)) {
  1L
} else {
  max(1L, min(2L, as.integer(detected_cores)))
}
famoye_example_n_obs <- .famoye_example_positive_integer(
  "RPBNB_FAMOYE_N_OBS", 500L
)
famoye_example_draws <- .famoye_example_positive_integer(
  "RPBNB_FAMOYE_DRAWS", 100L
)
famoye_example_cores <- .famoye_example_positive_integer(
  "RPBNB_FAMOYE_N_CORES", famoye_default_cores
)
```

Select data and print configuration:

```r
famoye_example_n_obs <- .famoye_example_observation_count(
  famoye_example_n_obs, nrow(data)
)
data <- data[seq_len(famoye_example_n_obs), , drop = FALSE]
cat("Observations :", nrow(data), "\n")
cat("Draws        :", famoye_example_draws, "\n")
cat("Cores asked  :", famoye_example_cores, "\n")
```

Pass `draws = famoye_example_draws` and
`n_cores = famoye_example_cores` to `fit_rpbnb_tmb()`.

- [ ] **Step 4: Run tests and verify GREEN for Task 1**

Run the Step 2 command.

Expected: configuration tests and the bounded-demo regression pass.

### Task 2: Fit diagnostics and allocation guidance

**Files:**
- Modify: `tests/testthat/test-example-famoye.R`
- Modify: `inst/fit_rpbnb_diff_famoye.R`

**Interfaces:**
- Produces: `.famoye_example_fit_diagnostics(fit, elapsed)` for runtime output.
- Produces: `.famoye_example_with_memory_guidance(value)` preserving unrelated errors and translating `std::bad_alloc`.

- [ ] **Step 1: Write failing diagnostics and condition tests**

Append:

```r
test_that("Famoye example reports fit diagnostics", {
  helpers <- load_famoye_example_helpers()
  fit <- list(
    parallel = list(requested = 2L, realized = 2L),
    optimizer = list(convergence = 0L, message = "relative convergence"),
    sdreport = list(pdHess = TRUE)
  )
  output <- paste(
    capture.output(
      helpers$.famoye_example_fit_diagnostics(fit, 1.25)
    ),
    collapse = "\n"
  )

  expect_match(output, "1.25 s")
  expect_match(output, "requested=2")
  expect_match(output, "realized=2")
  expect_match(output, "code=0")
  expect_match(output, "positive-definite Hessian: yes")
})

test_that("Famoye example translates bad_alloc into workload guidance", {
  helpers <- load_famoye_example_helpers()
  expect_error(
    helpers$.famoye_example_with_memory_guidance(
      stop("Caught exception 'std::bad_alloc'", call. = FALSE)
    ),
    "Restart R.*RPBNB_FAMOYE_N_OBS.*RPBNB_FAMOYE_DRAWS.*RPBNB_FAMOYE_N_CORES"
  )
})

test_that("Famoye example preserves unrelated error conditions", {
  helpers <- load_famoye_example_helpers()
  original <- structure(
    list(message = "unrelated failure", call = NULL),
    class = c("famoye_example_test_error", "error", "condition")
  )
  caught <- tryCatch(
    helpers$.famoye_example_with_memory_guidance(stop(original)),
    error = identity
  )

  expect_s3_class(caught, "famoye_example_test_error")
  expect_identical(conditionMessage(caught), "unrelated failure")
})
```

- [ ] **Step 2: Run the focused test and verify RED**

Run:

```powershell
& 'C:\Program Files\R\R-4.5.1\bin\Rscript.exe' -e ".libPaths(c('C:/Users/litabook/repos/rpbnb_tmb/.worktrees/.rlib', .libPaths())); testthat::test_file('tests/testthat/test-example-famoye.R', reporter='summary')"
```

Expected: failures because the diagnostics and condition helpers are absent.

- [ ] **Step 3: Add diagnostics and condition helpers**

Add:

```r
.famoye_example_fit_diagnostics <- function(fit, elapsed) {
  cat(sprintf("\nEstimation finished in %.2f s\n", elapsed))
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
}

.famoye_example_with_memory_guidance <- function(value) {
  tryCatch(
    value,
    error = function(error) {
      if (grepl(
        "std::bad_alloc", conditionMessage(error), fixed = TRUE
      )) {
        stop(
          paste0(
            "TMB exhausted memory while constructing the AD tape. ",
            "Restart R or reduce RPBNB_FAMOYE_N_OBS, ",
            "RPBNB_FAMOYE_DRAWS, or RPBNB_FAMOYE_N_CORES."
          ),
          call. = FALSE
        )
      }
      stop(error)
    }
  )
}
```

Wrap the fit:

```r
fit <- .famoye_example_with_memory_guidance(
  fit_rpbnb_tmb(
    formula_1 = f1,
    formula_2 = f2,
    data = data,
    random_1 = "hhninc",
    random_2 = "educ",
    dependence = "famoye",
    draws = famoye_example_draws,
    seed = 20240712,
    control = rpbnb_tmb_control(
      print_level = 1,
      n_cores = famoye_example_cores
    )
  )
)
```

Replace the existing elapsed-time line with:

```r
.famoye_example_fit_diagnostics(fit, t_fit)
```

- [ ] **Step 4: Run the focused test and verify GREEN**

Run the Step 2 command.

Expected: all six Famoye example tests pass.

### Task 3: Rebuild, regression verification, and commit

**Files:**
- Verify: `inst/fit_rpbnb_diff_famoye.R`
- Verify: `tests/testthat/test-example-famoye.R`
- Verify: `tests/testthat/test-parallel.R`

**Interfaces:**
- Produces a tested, sourceable memory-aware Famoye example.

- [ ] **Step 1: Reinstall the source package**

Run:

```powershell
& 'C:\Program Files\R\R-4.5.1\bin\R.exe' CMD INSTALL --preclean --library='C:\Users\litabook\repos\rpbnb_tmb\.worktrees\.rlib' .
```

Expected: installation succeeds with `-fopenmp` in compile and link commands.

- [ ] **Step 2: Run Famoye example, fit, and parallel tests**

Run:

```powershell
$env:NOT_CRAN='true'
& 'C:\Program Files\R\R-4.5.1\bin\Rscript.exe' -e ".libPaths(c('C:/Users/litabook/repos/rpbnb_tmb/.worktrees/.rlib', .libPaths())); testthat::test_dir('tests/testthat', filter='(example-famoye|fit-famoye|parallel)', reporter='summary', package='rpbnb.tmb', load_package='installed')"
```

Expected: all selected tests pass.

- [ ] **Step 3: Run a small real Famoye example**

Run:

```powershell
$env:RPBNB_FAMOYE_N_OBS='120'
$env:RPBNB_FAMOYE_DRAWS='20'
$env:RPBNB_FAMOYE_N_CORES='2'
& 'C:\Program Files\R\R-4.5.1\bin\Rscript.exe' -e ".libPaths(c('C:/Users/litabook/repos/rpbnb_tmb/.worktrees/.rlib', .libPaths())); source('inst/fit_rpbnb_diff_famoye.R')"
```

Expected: the script reports 120 observations, 20 draws, requested and
realized thread counts, optimizer status, and Hessian status without error.

- [ ] **Step 4: Review and commit the focused change**

Run:

```powershell
git diff --check
git diff -- inst/fit_rpbnb_diff_famoye.R tests/testthat/test-example-famoye.R tests/testthat/test-parallel.R
git add -- inst/fit_rpbnb_diff_famoye.R tests/testthat/test-example-famoye.R tests/testthat/test-parallel.R
git commit -m "feat: make Famoye example memory-aware"
```

Expected: only the Famoye example and its relevant tests are committed.
Existing `.gitignore` and `.Rbuildignore` changes remain unstaged.
