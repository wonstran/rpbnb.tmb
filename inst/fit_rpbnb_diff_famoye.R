#!/usr/bin/env Rscript
# =============================================================================
# rpbnb.tmb estimation demo -- DIFFERENT formulas, Famoye/Sarmanov dependence.
#
# Fits a random-parameter bivariate NB (RP-BNB) model to the German health-care
# counts in data/rwm1984_bnb.csv using the TMB backend:
#   * y1 = docvis  (doctor visits), y2 = hospvis (hospital visits)
#   * the two equations use DIFFERENT covariate sets (the package does not
#     require the margins to share a design)
#   * a fixed + random coefficient mix: a random slope on one CONTINUOUS
#     covariate per equation (hhninc in eq 1, educ in eq 2); all others fixed
#   * Famoye/Sarmanov dependence between the two counts
#
# Usage demo on real data -- there is no known ground truth, so the script
# reports the fitted model, predictions, and convergence diagnostics.
# Run from the package root:
#   library(rpbnb.tmb); source("inst/fit_rpbnb_diff_famoye.R")
#
# Optional memory-aware configuration before sourcing:
#   Sys.setenv(
#     RPBNB_FAMOYE_N_OBS = 500,
#     RPBNB_FAMOYE_DRAWS = 100,
#     RPBNB_FAMOYE_N_CORES = 2
#   )
# =============================================================================

library(rpbnb.tmb)

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

sep <- function() cat("\n", paste(rep("=", 72), collapse = ""), "\n", sep = "")
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

# ---- 1. Data ----------------------------------------------------------------
data <- read.csv(file.path("data", "rwm1984_bnb.csv"))
famoye_example_n_obs <- .famoye_example_observation_count(
  famoye_example_n_obs, nrow(data)
)
data <- data[seq_len(famoye_example_n_obs), , drop = FALSE]
cat("=== RP-BNB (different formulas, Famoye) on rwm1984 health counts ===\n")
cat("Observations :", nrow(data), "\n")
cat("Draws        :", famoye_example_draws, "\n")
cat("Cores asked  :", famoye_example_cores, "\n")

# ---- 2. Model specification -------------------------------------------------
# The two equations carry DIFFERENT covariates. Random coefficients sit on
# CONTINUOUS regressors only (dummies would be weakly identified: NB dispersion
# vs random-coefficient scale), and each random covariate must appear in its
# own equation's formula.
f1 <- docvis  ~ age + hhninc + educ + female + married + kids
f2 <- hospvis ~ age + educ + outwork + female + self
cat("Equation 1   :", deparse(f1), "\n")
cat("Equation 2   :", deparse(f2), "\n\n")

t_fit <- system.time(
  fit <- .famoye_example_with_memory_guidance(
    fit_rpbnb_tmb(
      formula_1  = f1,
      formula_2  = f2,
      data       = data,
      random_1   = "hhninc", # random slope on hhninc (eq 1), continuous, in f1
      random_2   = "educ",   # random slope on educ (eq 2), continuous, in f2
      dependence = "famoye",
      draws      = famoye_example_draws,
      seed       = 20240712,
      control    = rpbnb_tmb_control(
        print_level = 1,
        n_cores     = famoye_example_cores
      )
    )
  )
)[["elapsed"]]
.famoye_example_fit_diagnostics(fit, t_fit)

# ---- 3. Model summary -------------------------------------------------------
sep(); cat("MODEL SUMMARY\n"); sep()
summary(fit)

# ---- 4. Fitted means (predict) ----------------------------------------------
sep(); cat("FITTED MEANS (predict) -- first 6 observations\n"); sep()
cat("(predict.rpbnb_tmb_fit currently returns the raw coefficients;\n",
    " fitted means are available as fit$mu1 / fit$mu2)\n", sep = "")
print(head(data.frame(mu1 = fit$mu1, mu2 = fit$mu2)))

# ---- 5. Dependence ----------------------------------------------------------
sep(); cat("DEPENDENCE (sdreport)\n"); sep()
sdr <- fit$sdreport
if (!is.null(sdr)) {
  sdr_sum <- summary(sdr, "report")
  print(sdr_sum)
}

# ---- 6. Marginal effects (AME) -----------------------------------------------
sep(); cat("AVERAGE MARGINAL EFFECTS (AME)\n"); sep()
rpbnb_tmb_marginal_effects(fit, which = "both")

# ---- 7. Elasticities / semi-elasticities (AME) --------------------------------
sep(); cat("ELASTICITIES / SEMI-ELASTICITIES (AME)\n"); sep()
rpbnb_tmb_elasticities(fit, which = "both")
