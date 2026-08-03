#library(rpbnb.tmb)
devtools::load_all("/home/wonstran/repos/rpbnb.tmb")

sep <- function() cat("\n", paste(rep("=", 72), collapse = ""), "\n", sep = "")
setwd("/home/wonstran/repos/truck")
# detectCores() reports the hardware CPU count and ignores OMP_NUM_THREADS, so
# on this box it returns 32 while OpenMP -- and therefore TMB -- will only ever
# grant 8.  Asking for 32 just produced a "using 8 supported threads" warning on
# every run; take the binding limit instead.
n_cores <- 24L
draws <- 500L
dependence <- copula("frank")
optimizer <- "laplace"

# Path components stay separate so the same call resolves on POSIX and Windows.
truck_data <- read.csv(file.path("data", "export_dense_all.csv"))

# ---- Standardize the continuous predictors ---------------------------------
# SR40_MI3 carries the random coefficient in BOTH equations and runs over
# [17.7, 68.5] -- it is never near zero.  So SR40_MI3 * (b + sd * u_i) is not a
# random slope, it is a random INTERCEPT in disguise: at the uncentred fit's
# own estimates the latent contributed a SD of 0.90 (eq 1) and 2.04 (eq 2) to
# the log-mean.  That absorbs the overdispersion the NB dispersion is there to
# carry, so log_m1/log_m2 slid to the Poisson limit (m = 4.9e-7, 1.5e-7) where
# the NB curvature vanishes, leaving an indefinite Hessian, seven negative
# variance diagonals and NA standard errors on five equation-2 coefficients.
# Centring turns it back into a slope: the latent variance is now zero at the
# mean of SR40_MI3 and grows toward its tails, which no longer competes with
# the dispersion for the same variance.
#
# Scaling is the second half.  IRI_ME spans 496 next to 0/1 indicators, giving
# kappa(X1) = 3279 and a Hessian condition number near 1e7; PORT's trust region
# collapses on that long before the score does.  Scaled, kappa(X1) = 2.5.
#
# Binary regressors are left as 0/1: centring them only moves the intercept,
# and it costs the "one unit = present" reading the marginal effects rely on.
continuous_vars <- c("LNAADT_3", "SR40_MI3", "MPD_ME", "MPD_STD", "IRI_ME",
                     "ACCPNTS", "AUXLNUM", "DP10_ME", "DP01_ME")
scaling <- lapply(truck_data[continuous_vars], function(x) {
  c(center = mean(x), scale = stats::sd(x))
})
for (v in continuous_vars) {
  truck_data[[v]] <- (truck_data[[v]] - scaling[[v]]["center"]) / scaling[[v]]["scale"]
}
cat("Standardized (centred and scaled) predictors:\n")
print(round(do.call(rbind, scaling), 4))
# Coefficients, marginal effects and elasticities below are therefore on the
# STANDARDIZED scale: a unit is one standard deviation of the raw variable.
# Multiply a coefficient or an AME by 1/scale to return it to raw units.  The
# elasticities for these variables are NOT interpretable as printed -- a
# centred regressor has mean zero, so the X-bar/Y-bar factor the elasticity is
# built from no longer refers to anything; recompute them from the raw means.
# `scaling` is saved next to the fit so that is a lookup, not a rerun.

cat("=== RP-BNB on truck all crashes (Laplace) ===\n")
cat("Dependence   : copula(\"", dependence$family, "\")\n", sep = "")
cat("Cores asked  :", n_cores, "\n")
cat("Optimizer    :", optimizer, "\n")
cat("Draws (if SML):", draws, "\n\n")

f1 <- ALL_3  ~ LNAADT_3+SR40_MI3+MPD_ME+MPD_STD+IRI_ME+RUT_L+SP50GE+ACCPNTS+SIGNAL1+NEAR_SIG+AUXLNUM+DP10_ME
f2 <- C_HV ~ LNAADT_3+SR40_MI3+MPD_ME+MPD_STD+IRI_ME+RUT_7+RUT_9+SP50GE+ACCPNTS+SIGNAL1+NEAR_SIG+AUXLNUM+DP01_ME
r1 <- c("SR40_MI3")
r2 <- c("SR40_MI3")

cat("Equation 1   :", deparse(f1), "\n")
cat("Equation 2   :", deparse(f2), "\n")
cat("Random Parameters 1:", deparse(r1), "\n")
cat("Random Parameters 2:", deparse(r2), "\n\n")

t_fit <- system.time(
  fit <- fit_rpbnb_tmb(
    formula_1  = f1,
    formula_2  = f2,
    data       = truck_data,
    random_1   = r1,
    random_2   = r2,
    dependence = dependence,
    seed       = 20240712,
    method     = optimizer,
    draws = draws,
    # No max_workload override: the header's claim that the default suffices is
    # only true if the guard is actually left on.  This fit's workload is
    # nrow(data) * draws = 1.74e6, comfortably inside the default (which
    # rpbnb_tmb_max_workload() sizes from available memory, so it is not a fixed
    # number), so the guard costs nothing here and still catches an accidentally
    # oversized respecification.
    control    = rpbnb_tmb_control(
      print_level = 1,
      n_cores     = n_cores,
      max_threads = n_cores
    )
  )
)[["elapsed"]]

cat(sprintf("\nEstimation finished in %.2f s\n", t_fit))
# Deliberately not reporting gc() figures: the TMB tape lives on the C++ heap
# and is invisible to R's garbage collector, so gc() would understate exactly
# the quantity this script exists to test.  Watch the process working set
# externally (Task Manager, or inst/benchmark_memory.R) if a number is needed.
cat(sprintf(
  "TMB threads: requested=%d, realized=%d\n",
  fit$parallel$requested, fit$parallel$realized
))
cat(sprintf("Optimizer: code=%d, message=%s\n",
            fit$optimizer$convergence, fit$optimizer$message))
cat("sdreport positive-definite Hessian:",
    if (isTRUE(fit$sdreport$pdHess)) "yes" else "no", "\n")

# ---- Persist the fit ---------------------------------------------------
# This fit costs ~55 min.  Without it on disk every follow-up diagnostic means
# paying that again, so write it before any post-estimation step can fail.
stamp <- format(Sys.time(), "%Y-%m-%d-%H%M%S")
# A distinct prefix: the uncentred fits under fit_normal_dense_<stamp>.rds are
# on a different design and their coefficients are not comparable to these.
fit_path <- file.path("results", paste0("fit_normal_dense_centered_", stamp, ".rds"))
saveRDS(list(fit = fit, scaling = scaling), fit_path)
cat("Fit object saved to:", fit_path, "\n")

# ---- Convergence diagnostics -------------------------------------------
# print_level = 1 buries these under thousands of TMB inner-iteration lines, so
# restate them loudly: a boundary-bound dispersion or a non-PD Hessian makes the
# standard errors below meaningless, and that must not be easy to scroll past.
diagnostics <- character(0)
if (!isTRUE(fit$sdreport$pdHess)) {
  diagnostics <- c(diagnostics,
                   "Hessian is not positive definite: standard errors are unreliable.")
}
if (length(fit$boundary_report)) {
  diagnostics <- c(diagnostics,
                   paste0("Parameters at a constraint bound: ",
                          paste(paste0(fit$boundary_report,
                                       " (", fit$boundary_sides, ")"),
                                collapse = ", ")))
}
if (length(diagnostics)) {
  sep(); cat("CONVERGENCE WARNINGS\n"); sep()
  cat(paste0("  * ", diagnostics, collapse = "\n"), "\n")
} else {
  cat("No boundary or Hessian warnings.\n")
}

# ---- Model summary -----------------------------------------------------
sep(); cat("MODEL SUMMARY\n"); sep()
model_summary_output <- capture.output(summary(fit))
cat(model_summary_output, sep = "\n")
cat("\n")

# ---- Fitted means (predict) ---------------------------------------------
sep(); cat("FITTED MEANS (predict) -- first 6 observations\n"); sep()
print(head(predict(fit)))

# ---- Dependence -----------------------------------------------------------
sep(); cat("DEPENDENCE (sdreport)\n"); sep()
if (!is.null(fit$sdreport)) print(summary(fit$sdreport, "report"))

# ---- Marginal effects (AME) ------------------------------------------------
sep(); cat("AVERAGE MARGINAL EFFECTS (AME)\n"); sep()
marginal_effects_output <- capture.output(
  marginal_effects <- rpbnb_tmb_marginal_effects(fit, which = "both")
)
cat(marginal_effects_output, sep = "\n")
cat("\n")

# ---- Elasticities / semi-elasticities (AME) --------------------------------
sep(); cat("ELASTICITIES / SEMI-ELASTICITIES (AME)\n"); sep()
elasticities_output <- capture.output(
  elasticities <- rpbnb_tmb_elasticities(fit, which = "both")
)
cat(elasticities_output, sep = "\n")
cat("\n")

# ---- Results in the covariates' original units ------------------------------
# The two sections above are per standard deviation, and the elasticities of the
# centred regressors are identically zero because their sample mean is zero.
# `scaling` restates both in raw units: AMEs divide by the scale, elasticities
# recover the x-bar factor.  The fitted design is NOT rebuilt -- doing so would
# reintroduce the (center/scale) * dev random intercept that centring removed.
#
# LNAADT_3 is log(AADT), so without `log_vars` the elasticity formula treats it
# as an ordinary regressor and returns x-bar * b -- the elasticity with respect
# to the LOG of traffic, 8.59, where the elasticity with respect to traffic is
# the coefficient itself, 0.871.  Nothing in the printed table distinguishes
# the two.  Naming it here reports both diagnostics per vehicle per day.
log_vars <- "LNAADT_3"
sep(); cat("MARGINAL EFFECTS AND ELASTICITIES IN ORIGINAL UNITS\n"); sep()
raw_scale_output <- capture.output({
  cat("\n--- Marginal effects (AME, original units) ---\n")
  raw_me <- rpbnb_tmb_marginal_effects(fit, which = "both",
                                       scaling = scaling, log_vars = log_vars)
  cat("\n--- Elasticities (AME, original units) ---\n")
  raw_el <- rpbnb_tmb_elasticities(fit, which = "both",
                                   scaling = scaling, log_vars = log_vars)
})
cat(raw_scale_output, sep = "\n")
cat("\n")

# ---- Export requested results ----------------------------------------------
results_path <- rpbnb.tmb:::.write_truck_results_markdown(
  model_summary = model_summary_output,
  marginal_effects = marginal_effects_output,
  elasticities = elasticities_output,
  dependence = fit$dependence,
  method = fit$method,
  raw_scale = raw_scale_output,
  scaling = scaling
)
cat("Results written to:", results_path, "\n")
