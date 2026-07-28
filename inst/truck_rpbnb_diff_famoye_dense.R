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

#library(rpbnb.tmb)
library(here)

devtools::load_all()

sep <- function() cat("\n", paste(rep("=", 72), collapse = ""), "\n", sep = "")


# ---- Settings ---------------------------------------------------------------
# draws is bounded by memory, not by taste: peak working set runs about
# 12 kB per observation-draw (see TAPE_CALIBRATION / inst/benchmark_memory.R),
# so 3874 obs x 200 draws is ~8 GiB peak -- right at the default
# max_workload. Raise both deliberately if you have the memory.
n_obs <- 5000L
draws <- 500L
n_cores <- 12 #parallel::detectCores(logical = FALSE)
if (is.na(n_cores)) n_cores <- 1L

# ---- 1. Data ----------------------------------------------------------------
data <- read.csv(here("inst\\extdata", "export_dense_all.csv"))

data <- data[seq_len(min(n_obs, nrow(data))), , drop = FALSE]  # comment this line to use all obs
cat("=== RP-BNB on truck all crashes ===\n")
cat("Observations :", nrow(data), "\n")
cat("Draws        :", draws, "\n")
cat("Cores asked  :", n_cores, "\n")

# ---- 2. Model specification -------------------------------------------------
# The two equations carry DIFFERENT covariates. Random coefficients sit on
# CONTINUOUS regressors only (dummies would be weakly identified: NB dispersion
# vs random-coefficient scale), and each random covariate must appear in its
# own equation's formula.
f1 <- ALL_3  ~ LNAADT_3+SR40_MI3+MPD_ME+MPD_STD+IRI_ME+RUT_L+SP50GE+ACCPNTS+SIGNAL1+NEAR_SIG+AUXLNUM++DP10_ME
f2 <- C_DISTR ~ LNAADT_3+SR40_MI3+MPD_ME+MPD_STD+RUT_9+ACCPNTS+SIGNAL1+NEAR_SIG+AUXLNUM+DP01_ME+CS_MINAB

#f2 <- C_DISTR ~ SR40_MI3+MPD_ME+LNAADT_3+IRI_ME+SP50LE+ACCPNTS+SIGNAL1+NEAR_SIG+CS_MINAB+DP10_ME

cat("Equation 1   :", deparse(f1), "\n")
cat("Equation 2   :", deparse(f2), "\n\n")

t_fit <- system.time(
  fit <- fit_rpbnb_tmb(
    formula_1  = f1,
    formula_2  = f2,
    data       = data,
    random_1   = c("SR40_MI3","AUXLNUM","MPD_ME","MPD_STD"),   # random slope on hhninc (eq 1), continuous, in f1
    random_2   = c("SR40_MI3","RUT_9","MPD_ME","MPD_STD"),     # random slope on educ   (eq 2), continuous, in f2
    dependence = "famoye",
    draws      = draws,
    seed       = 20240712,
    control    = rpbnb_tmb_control(
      print_level = 1,
      n_cores     = n_cores,
      max_threads = n_cores,
      max_workload = 32*1024^3/12083,
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
model_summary_output <- capture.output(summary(fit))
cat(model_summary_output, sep = "\n")
cat("\n")

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
marginal_effects_output <- capture.output(
  marginal_effects <- rpbnb_tmb_marginal_effects(fit, which = "both")
)
cat(marginal_effects_output, sep = "\n")
cat("\n")

# ---- 7. Elasticities / semi-elasticities (AME) --------------------------------
sep(); cat("ELASTICITIES / SEMI-ELASTICITIES (AME)\n"); sep()
elasticities_output <- capture.output(
  elasticities <- rpbnb_tmb_elasticities(fit, which = "both")
)
cat(elasticities_output, sep = "\n")
cat("\n")

# ---- 8. Export requested results -------------------------------------------
results_path <- .write_truck_results_markdown(
  model_summary = model_summary_output,
  marginal_effects = marginal_effects_output,
  elasticities = elasticities_output
)
cat("Results written to:", results_path, "\n")
