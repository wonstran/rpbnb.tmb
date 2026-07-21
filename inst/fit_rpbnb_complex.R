#!/usr/bin/env Rscript
# =============================================================================
# Fit a random-parameter bivariate NB (RP-BNB) model to the complex sample data
# using rpbnb.tmb (TMB backend). Compare the estimates against known ground-truth.
#
#   * 8 independent variables (3 continuous + 5 dummy)
#   * random coefficients on x_age (eq 1) and x_income (eq 2)
#   * GENUINE Famoye/Sarmanov dependence (non-zero lambda) between y1 and y2
#   * TMB simulated-likelihood estimation with automatic differentiation
# =============================================================================

library(rpbnb.tmb)

# ---- 1. Load data + ground truth -------------------------------------------
data  <- read.csv(file.path("data", "simulated_rpbnb_dependent.csv"))
truth <- readRDS(file.path("data", "simulated_rpbnb_dependent_truth.rds"))

cat("=== Fitting RP-BNB (TMB) to complex sample data ===\n")
cat("Observations :", nrow(data), "\n")
cat(sprintf("True lambda  : %.4f (genuine Famoye/Sarmanov dependence)\n\n", truth$lambda))

# ---- 2. Model specification -------------------------------------------------
# Same 8 covariates in both equations; random coefficients matched to the DGP.
rhs <- ~ x_age + x_income + x_score + d_female + d_urban + d_married + d_college + d_smoker

f1 <- update(rhs, y1 ~ .)
f2 <- update(rhs, y2 ~ .)

t_fit <- system.time(
  fit <- fit_rpbnb_tmb(
    formula_1 = f1,
    formula_2 = f2,
    data      = data,
    random_1  = truth$random_names_1,   # c("x_age")
    random_2  = truth$random_names_2,   # c("x_income")
    dependence = "famoye",
    draws     = 500,
    seed      = 20240712,
    control   = rpbnb_tmb_control(
      print_level = 2,
      n_cores     = parallel::detectCores()
    )
  )
)[["elapsed"]]

cat(sprintf("\nEstimation finished in %.1f s\n\n", t_fit))

# ---- 3. Model summary -------------------------------------------------------
print(fit)

# ---- 4. Compare estimated MEAN coefficients vs truth ------------------------
compare_means <- function(prefix, true_beta) {
  est <- fit$coef[paste0(prefix, ":", names(true_beta))]
  se  <- fit$se[paste0(prefix, ":", names(true_beta))]
  data.frame(
    Variable  = names(true_beta),
    True      = as.numeric(true_beta),
    Estimate  = as.numeric(est),
    StdErr    = as.numeric(se),
    z_vs_true = (as.numeric(est) - as.numeric(true_beta)) / as.numeric(se),
    row.names = NULL
  )
}

cat("\n", paste(rep("=", 72), collapse = ""), "\n", sep = "")
cat("EQUATION 1 (y1) - mean coefficients: estimate vs true\n")
cat(paste(rep("=", 72), collapse = ""), "\n")
print(compare_means("b1", truth$beta1), digits = 4)

cat("\nEQUATION 2 (y2) - mean coefficients: estimate vs true\n")
cat(paste(rep("=", 72), collapse = ""), "\n")
print(compare_means("b2", truth$beta2), digits = 4)

# ---- 5. Compare random-coefficient SDs vs truth -----------------------------
compare_sd <- function(prefix, random_spec) {
  nm  <- names(random_spec)
  # rpbnb.tmb uses log_sd as the scale label (normal distribution)
  est <- exp(fit$coef[paste0("log_sd", prefix, ":", nm)])
  true_sd <- vapply(random_spec, function(z) z$sd, numeric(1))
  data.frame(
    Variable = nm,
    True_SD  = as.numeric(true_sd),
    Est_SD   = as.numeric(est),
    row.names = NULL
  )
}

cat("\n", paste(rep("=", 72), collapse = ""), "\n", sep = "")
cat("RANDOM-COEFFICIENT STANDARD DEVIATIONS: estimate vs true\n")
cat(paste(rep("=", 72), collapse = ""), "\n")
cat("Equation 1:\n")
print(compare_sd("1", truth$random_1), digits = 4)
cat("\nEquation 2:\n")
print(compare_sd("2", truth$random_2), digits = 4)

# ---- 6. Dispersion & dependence vs truth ------------------------------------
cat("\n", paste(rep("=", 72), collapse = ""), "\n", sep = "")
cat("DISPERSION & DEPENDENCE: estimate vs true\n")
cat(paste(rep("=", 72), collapse = ""), "\n")
disp_tbl <- data.frame(
  Parameter = c("m1", "m2", "lambda"),
  True      = c(truth$dispersion[["m1"]], truth$dispersion[["m2"]], truth$lambda),
  Estimate  = c(fit$m1, fit$m2, coef(fit)["z_dep"]),
  row.names = NULL
)
print(disp_tbl, digits = 4)

cat("\nlogLik =", round(fit$logLik, 2),
    "  AIC =", round(AIC(fit), 2),
    "  BIC =", round(BIC(fit), 2),
    "  npar =", fit$npar, "\n")

cat(sprintf("\nNote: TRUE lambda = %.4f (genuine Famoye/Sarmanov dependence)\n",
            truth$lambda))
cat("      The TMB backend should recover a significantly non-zero lambda.\n")
cat("      The raw cor(y1, y2) may be modest -- the Famoye e^{-y} tilt\n")
cat("      concentrates dependence in the low-count region.\n")
