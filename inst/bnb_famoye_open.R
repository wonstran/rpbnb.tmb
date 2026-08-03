library(rpbnb)

sep <- function() cat("\n", paste(rep("=", 72), collapse = ""), "\n", sep = "")
setwd("C:\\Users\\litabook\\repos\\truck")

#n_obs <- 5000L
draws <- 500L
n_cores <- 12 #parallel::detectCores(logical = FALSE)
if (is.na(n_cores)) n_cores <- 1L

# ---- 1. Data ----------------------------------------------------------------
data <- read.csv(file.path("data", "export_open_all.csv"))

#data <- data[seq_len(min(n_obs, nrow(data))), , drop = FALSE]  # comment this line to use all obs
#cat("=== RP-BNB on truck all crashes ===\n")
cat("Observations :", nrow(data), "\n")
cat("Draws        :", draws, "\n")
cat("Cores asked  :", n_cores, "\n")

# ---- 2. Model specification -------------------------------------------------
# The two equations carry DIFFERENT covariates. Random coefficients sit on
# CONTINUOUS regressors only (dummies would be weakly identified: NB dispersion
# vs random-coefficient scale), and each random covariate must appear in its
# own equation's formula.
f1 <- ALL_3  ~ SR40_MI3 + MPD_ME + LNAADT_3 + IRI_ME + G_ABG2 + SP50LE + ACCPNTS + SIGNAL1 + NEAR_SIG + CS_MINAB + DP10_ME + RUT_L
f2 <- C_HV ~ SR40_MI3+MPD_ME+LNAADT_3+IRI_ME+SP50LE+ACCPNTS+SIGNAL1+NEAR_SIG+DP10_ME+RUT_890+G_ABG2

cat("Equation 1   :", deparse(f1), "\n")
cat("Equation 2   :", deparse(f2), "\n\n")

t_fit <- system.time(
  fit <- fit_bnb(
    formula_1  = f1,
    formula_2  = f2,
    data       = data,
    dependence = "famoye",
    control    = rpbnb_control(
      print_level = 2
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
#sep(); cat("AVERAGE MARGINAL EFFECTS (AME)\n"); sep()
#marginal_effects_output <- capture.output(
#  marginal_effects <- rpbnb_tmb_marginal_effects(fit, which = "both")
#)
#cat(marginal_effects_output, sep = "\n")
#cat("\n")

# ---- 7. Elasticities / semi-elasticities (AME) --------------------------------
#sep(); cat("ELASTICITIES / SEMI-ELASTICITIES (AME)\n"); sep()
#elasticities_output <- capture.output(
#  elasticities <- rpbnb_tmb_elasticities(fit, which = "both")
#)
#cat(elasticities_output, sep = "\n")
#cat("\n")
