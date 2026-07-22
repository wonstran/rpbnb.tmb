#!/usr/bin/env Rscript
# =============================================================================
# rpbnb.tmb estimation demo -- DIFFERENT formulas, Gaussian COPULA dependence.
#
# Fits a random-parameter bivariate NB (RP-BNB) model to the German health-care
# counts in data/rwm1984_bnb.csv using the TMB backend:
#   * y1 = docvis  (doctor visits), y2 = hospvis (hospital visits)
#   * the two equations use DIFFERENT covariate sets (the package does not
#     require the margins to share a design)
#   * a fixed + random coefficient mix: a random slope on one CONTINUOUS
#     covariate per equation (hhninc in eq 1, educ in eq 2); all others fixed
#   * discrete Gaussian-copula dependence (rho) between the two counts
#
# RUNTIME: the copula RP path (discrete-copula pmf + per-draw NB CDF corners) is
# noticeably heavier than the Famoye path. Lower draws or use a data subset for
# a quicker look. Run from the package root:
#   library(rpbnb.tmb); source("inst/fit_rpbnb_diff_copula.R")
# =============================================================================

library(rpbnb.tmb)

sep <- function() cat("\n", paste(rep("=", 72), collapse = ""), "\n", sep = "")

# ---- 1. Data ----------------------------------------------------------------
data <- read.csv(file.path("data", "rwm1984_bnb.csv"))
# Use first 500 observations (copula path is heavier; fewer obs + fewer draws
# keeps the AD tape manageable and reduces the risk of NaN gradients).
data <- data[1:500, ]
cat("=== RP-BNB (different formulas, Gaussian copula) on rwm1984 health counts ===\n")
cat("Observations :", nrow(data), "\n")

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
    dependence = copula("normal"),   # Gaussian copula; try "frank" or "kimeldorf"
    draws      = 100,        # fewer draws for quicker demo (copula path is heavy)
    seed       = 20240712,
    control    = rpbnb_tmb_control(
      print_level = 1,
      n_cores     = parallel::detectCores()
    )
  )
)[["elapsed"]]
cat(sprintf("\nEstimation finished in %.1f s\n", t_fit))

# ---- 3. Model summary -------------------------------------------------------
sep(); cat("MODEL SUMMARY\n"); sep()
summary(fit)

# ---- 4. Fitted means (predict) ----------------------------------------------
sep(); cat("FITTED MEANS (predict) -- first 6 observations\n"); sep()
cat("(predict.rpbnb_tmb_fit currently returns the raw coefficients;\n",
    " fitted means are available as fit$mu1 / fit$mu2)\n", sep = "")
print(head(data.frame(mu1 = fit$mu1, mu2 = fit$mu2)))

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
sep(); cat("ELASTICITIES / SEMI-ELASTICITIES (AME)\n"); sep()
rpbnb_tmb_elasticities(fit, which = "both")
