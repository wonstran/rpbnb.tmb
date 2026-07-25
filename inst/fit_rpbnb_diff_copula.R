#!/usr/bin/env Rscript
# =============================================================================
# rpbnb.tmb estimation demo -- DIFFERENT formulas, Frank COPULA dependence.
#
# Fits a random-parameter bivariate NB (RP-BNB) model to the German health-care
# counts in data/rwm1984_bnb.csv using the TMB backend:
#   * y1 = docvis  (doctor visits), y2 = hospvis (hospital visits)
#   * the two equations use DIFFERENT covariate sets (the package does not
#     require the margins to share a design)
#   * a fixed + random coefficient mix: a random slope on one CONTINUOUS
#     covariate per equation (hhninc in eq 1, educ in eq 2); all others fixed
#   * discrete Frank-copula dependence (theta) between the two counts
#
# RUNTIME: the copula RP path (discrete-copula pmf + per-draw NB CDF corners) is
# noticeably heavier than the Famoye path. Lower draws or use a data subset for
# a quicker look. Run from the package root:
#   library(rpbnb.tmb); source("inst/fit_rpbnb_diff_copula.R")
# =============================================================================

library(rpbnb.tmb)

sep <- function() cat("\n", paste(rep("=", 72), collapse = ""), "\n", sep = "")

# ---- Settings ---------------------------------------------------------------
# Frank costs about 2.9x Famoye per observation-draw (see TAPE_CALIBRATION /
# inst/benchmark_memory.R), so the same memory budget buys proportionally
# fewer draws. 3874 obs x 100 draws x 2.9 is ~1.1e6 weighted units; raise
# max_workload deliberately if you have the memory for more.
n_obs <- 5000L
draws <- 100L
n_cores <- 8 #parallel::detectCores(logical = FALSE)
if (is.na(n_cores)) n_cores <- 1L

# ---- 1. Data ----------------------------------------------------------------
data <- read.csv(system.file(
  "extdata", "rwm1984_bnb.csv", package = "rpbnb.tmb", mustWork = TRUE
))
data <- data[seq_len(min(n_obs, nrow(data))), , drop = FALSE]
cat("=== RP-BNB (different formulas, Frank copula) on rwm1984 health counts ===\n")
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

t_fit <- system.time(
  fit <- fit_rpbnb_tmb(
    formula_1  = f1,
    formula_2  = f2,
    data       = data,
    random_1   = "hhninc",   # random slope on hhninc (eq 1), continuous, in f1
    random_2   = "educ",     # random slope on educ   (eq 2), continuous, in f2
    dependence = copula("frank"),   # try "normal" or "kimeldorf"
    draws      = draws,
    seed       = 20240712,
    control    = rpbnb_tmb_control(
      print_level = 1,
      n_cores     = n_cores,
      # 3874 x 100 x 2.9 = 1.1e6 weighted units, above the 7e5 default. This
      # is the deliberate opt-in the documentation asks for, not a workaround:
      # at ~12 kB peak per unit it commits to roughly 13 GiB. Lower `draws` if
      # you have less memory than that.
      max_workload = 1.2e6
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

# ---- 5. Dependence (copula report) ------------------------------------------
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
#sep(); cat("ELASTICITIES / SEMI-ELASTICITIES (AME)\n"); sep()
#rpbnb_tmb_elasticities(fit, which = "both")
