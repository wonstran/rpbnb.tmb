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
# =============================================================================

library(rpbnb.tmb)

sep <- function() cat("\n", paste(rep("=", 72), collapse = ""), "\n", sep = "")

# ---- Settings ---------------------------------------------------------------
# draws is bounded by memory, not by taste: peak working set runs about
# 12 kB per observation-draw (see TAPE_CALIBRATION / inst/benchmark_memory.R).
# Famoye carries weight 1.0, so this configuration is
#   1000 obs x 200 draws x 1.0 = 2.0e5 weighted units  (~2.3 GiB peak),
# comfortably inside the default budget. Raise both deliberately if you have
# the memory; the whole file is 3874 observations.
n_obs <- 1000L
draws <- 200L
n_cores <- 12 #parallel::detectCores(logical = FALSE)
if (is.na(n_cores)) n_cores <- 1L

# ---- 1. Data ----------------------------------------------------------------
data <- read.csv(system.file(
  "extdata", "rwm1984_bnb.csv", package = "rpbnb.tmb", mustWork = TRUE
))
data <- data[seq_len(min(n_obs, nrow(data))), , drop = FALSE]  # comment this line to use all obs
cat("=== RP-BNB (different formulas, Famoye) on rwm1984 health counts ===\n")
cat("Observations :", nrow(data), "\n")
cat("Draws        :", draws, "\n")
cat("Cores asked  :", n_cores, "\n")

# ---- 2. Model specification -------------------------------------------------
# The two equations carry DIFFERENT covariates. Random coefficients sit on
# CONTINUOUS regressors only (dummies would be weakly identified: NB dispersion
# vs random-coefficient scale), and each random covariate must appear in its
# own equation's formula.
f1 <- docvis  ~ age + hhninc + educ + female + married + kids
f2 <- hospvis ~ age + educ + outwork + female + self
cat("Equation 1   :", deparse(f1), "\n")
cat("Equation 2   :", deparse(f2), "\n\n")

# Expect a "pinned at a boundary" dependence warning here: the docvis/hospvis
# association in this real dataset is stronger than Famoye/Sarmanov's range
# can represent (confirmed independent of n_obs and of the z_dep starting
# value -- lam wants to land around 13 against a box of roughly +/-1 either
# way). That is a property of this data under this family, not a fitting
# problem to chase. See fit_rpbnb_diff_copula.R for a dependence family
# (Frank, unbounded theta) whose range covers an association this strong.
t_fit <- system.time(
  fit <- fit_rpbnb_tmb(
    formula_1  = f1,
    formula_2  = f2,
    data       = data,
    random_1   = "hhninc",   # random slope on hhninc (eq 1), continuous, in f1
    random_2   = "educ",     # random slope on educ   (eq 2), continuous, in f2
    dependence = "famoye",
    draws      = draws,
    seed       = 20240712,
    control    = rpbnb_tmb_control(
      print_level = 1,
      n_cores     = n_cores,
      max_threads = n_cores
      # The default budgets 80% of detected available memory, which covers the
      # 2.0e5 weighted units above. To state a budget yourself instead of
      # hand-computing units, pass e.g.
      #   max_workload = rpbnb_tmb_max_workload(budget_gib = 32)
    )
  )
)[["elapsed"]]
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

# ---- 3. Model summary -------------------------------------------------------
sep(); cat("MODEL SUMMARY\n"); sep()
summary(fit)

# ---- 4. Fitted means (predict) ----------------------------------------------
sep(); cat("FITTED MEANS (predict) -- first 6 observations\n"); sep()
print(head(predict(fit)))

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
