library(rpbnb)

sep <- function() cat("\n", paste(rep("=", 72), collapse = ""), "\n", sep = "")
setwd("C:\\Users\\litabook\\repos\\truck")

#n_obs <- 5000L
draws <- 500L
n_cores <- 12 #parallel::detectCores(logical = FALSE)
if (is.na(n_cores)) n_cores <- 1L

# ---- 1. Data ----------------------------------------------------------------
data <- read.csv(file.path("data", "export_dense_all.csv"))

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
f1 <- ALL_3  ~ LNAADT_3+SR40_MI3+MPD_ME+MPD_STD+IRI_ME+RUT_L+SP50GE+ACCPNTS+SIGNAL1+NEAR_SIG+AUXLNUM++DP10_ME
f2 <- C_HV ~ LNAADT_3+SR40_MI3+MPD_ME+MPD_STD+IRI_ME+RUT_7+RUT_9+SP50GE+ACCPNTS+SIGNAL1+NEAR_SIG+AUXLNUM+DP01_ME+CS_MINAB

cat("Equation 1   :", deparse(f1), "\n")
cat("Equation 2   :", deparse(f2), "\n\n")

t_fit <- system.time(
  fit <- fit_rpbnb(
    formula_1  = f1,
    formula_2  = f2,
    data       = data,
    random_1   = c("SR40_MI3","MPD_ME"),   # random slope on hhninc (eq 1), continuous, in f1
    random_2   = c("SR40_MI3"),     # random slope on educ   (eq 2), continuous, in f2
    dependence = "famoye",
    draws      = draws,
    seed       = 20240712,
    control    = rpbnb_control(
      print_level = 1,
      n_cores     = rpbnb_threads(),
      se_method   = "opg"    # fast BHHH/OPG SEs; the robust default is "numeric"
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

# ---- 3b. Significance of the boundary parameters (LR tests) -----------------
# The natural-scale table shows no Wald z/p for the random-coefficient SDs
# (sd1:hhninc, sd2:educ) or the NB2 dispersions (m1, m2): those are positive
# parameters whose null sits on the boundary of the parameter space (sd = 0, or
# m = 0 = Poisson), where the Wald ratio does not test that null. The correct
# test is a likelihood-ratio test against a properly nested restricted fit, with
# the 50:50 chi-square boundary correction (Self & Liang 1987).
#
# rpbnb_boundary_tests() runs all of them and merges the results: it refits each
# nested restricted model (dropping one random SD at a time -- keeping any other
# random coefficients in the equation -- or pinning a margin at its Poisson
# limit), reusing the full model's draws / seed so the simulated log-likelihoods
# compare on common random numbers. SEs are skipped for the refits (the LR test
# needs only logLik + df).
sep(); cat("BOUNDARY-PARAMETER SIGNIFICANCE (boundary-corrected LR tests)\n"); sep()

bt <- rpbnb_boundary_tests(fit, data,
                           control = rpbnb_control(print_level = 1, n_cores = rpbnb_threads(),
                                                   compute_se = FALSE))
print(bt)

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
