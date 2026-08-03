# NOT library(rpbnb.tmb): the installed 0.3.4 predates the `scaling`/`log_vars`
# arguments used below, and both diagnostic functions take `...`, so under the
# installed build those arguments are SILENTLY SWALLOWED -- the script would
# print standardized numbers under an "original units" heading and only fail
# later, at the exporter.  The source tree carries the same version number, so
# the version string cannot tell the two apart.  Switch back to library() only
# after a version bump and devtools::install().
devtools::load_all("/home/wonstran/repos/rpbnb.tmb")

sep <- function() cat("\n", paste(rep("=", 72), collapse = ""), "\n", sep = "")
setwd("/home/wonstran/repos/truck")
# detectCores() returns 32 on this box, and 32 threads is what killed the
# rsession during the dense Gaussian fit: TMB gives each thread its own tape and
# Hessian workspace, so the process reached 62 GB resident (90 GB virtual)
# against 60 GB of RAM with no swap, and the kernel OOM-killer took it.  Nothing
# warned first -- .check_tmb_workload() multiplies by n_threads only when
# parallel_tape = TRUE (see R/tmb_helpers.R:79), which is not set here, so the
# guard scored this fit at a few thousand against a budget of 4e6.
# 8 threads is the configuration that has actually completed.
n_cores <- 18L
draws <- 500L
dependence <- copula("kimeldorf")

# Path components stay separate: "inst\\extdata" is one directory named
# inst\extdata on POSIX, so the backslash form only resolves on Windows.
data <- read.csv(file.path("data", "export_open_all.csv"))

# ---- Standardize the continuous predictors ---------------------------------
# Both random-coefficient carriers are strictly positive and bounded away from
# zero -- SR40_MI3 over [17.9, 65.1] and MPD_ME over [0.46, 3.62] -- so
# x * (b + sd * u_i) is not a random slope but a random INTERCEPT in disguise,
# with a per-observation SD of sd * x.  That absorbs the overdispersion the NB
# dispersion exists to carry, and the two then compete for the same variance:
# on the dense Gaussian model the same specification drove log_m1/log_m2 to the
# Poisson limit, where the NB curvature vanishes and the Hessian stops being
# positive definite -- seven negative variance diagonals and NA standard errors.
# This data is overdispersed too (var/mean = 32.1 for ALL_3, 2.65 for C_HV), so
# there is real dispersion to lose.  Centring makes the latent variance zero at
# the mean of the carrier and nonzero only in its tails: a slope again.
#
# Scaling is the second half: IRI_ME spans 29-380 next to 0/1 indicators, giving
# kappa(X1) = 3055 and a Hessian condition number near 1e7, on which PORT's
# trust region collapses long before the score does.
#
# Binary regressors stay 0/1 -- centring them only moves the intercept, and it
# costs the "one unit = present" reading the marginal effects rely on.
continuous_vars <- c("SR40_MI3", "MPD_ME", "LNAADT_3", "IRI_ME",
                     "ACCPNTS", "CS_MINAB", "DP10_ME")
scaling <- lapply(data[continuous_vars], function(x) {
  c(center = mean(x), scale = stats::sd(x))
})
for (v in continuous_vars) {
  data[[v]] <- (data[[v]] - scaling[[v]]["center"]) / scaling[[v]]["scale"]
}
cat("Standardized (centred and scaled) predictors:\n")
print(round(do.call(rbind, scaling), 4))

cat("=== RP-BNB on truck all crashes (Laplace) ===\n")
cat("Dependence   : copula(\"", dependence$family, "\")\n", sep = "")
cat("Cores asked  :", n_cores, "\n")

f1 <- ALL_3  ~ SR40_MI3 + MPD_ME + LNAADT_3 + IRI_ME + G_ABG2 + SP50LE + ACCPNTS + SIGNAL1 + NEAR_SIG + CS_MINAB + DP10_ME + RUT_L
f2 <- C_HV ~SR40_MI3+MPD_ME+LNAADT_3+IRI_ME+SP50LE+ACCPNTS+SIGNAL1+NEAR_SIG+CS_MINAB+DP10_ME

cat("Equation 1   :", deparse(f1), "\n")
cat("Equation 2   :", deparse(f2), "\n\n")

t_fit <- system.time(
  fit <- fit_rpbnb_tmb(
    formula_1  = f1,
    formula_2  = f2,
    data       = data,
    random_1   = c("SR40_MI3", "MPD_ME"),
    random_2   = c("SR40_MI3"),
    dependence = dependence,
    seed       = 20240712,
    method     = "laplace",
    draws = draws,
    # No max_workload override: the header's claim that the default suffices is
    # only true if the guard is actually left on.
    control    = rpbnb_tmb_control(
      print_level = 1,
      n_cores     = n_cores,
      max_threads = n_cores,
      max_workload = Inf
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

# ---- Persist the fit and the transform --------------------------------------
# `scaling` has to outlive the session or the raw-unit tables below can never be
# regenerated without paying for the fit again, so the two are saved together.
stamp <- format(Sys.time(), "%Y-%m-%d-%H%M%S")
fit_path <- file.path("results", paste0("fit_frank_open_centered_", stamp, ".rds"))
saveRDS(list(fit = fit, scaling = scaling), fit_path)
cat("Fit object saved to:", fit_path, "\n")

# ---- Convergence diagnostics -------------------------------------------
# print_level = 1 buries these under thousands of TMB inner-iteration lines.  A
# non-PD Hessian makes every standard error below meaningless, and that must not
# be easy to scroll past.
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
# centred regressors are identically zero because their sample mean is zero --
# they print as 0.0000, which reads as "no effect" rather than "this number
# means nothing".  `scaling` restates both: AMEs divide by the scale,
# elasticities recover the x-bar factor.  The fitted design is NOT rebuilt --
# substituting raw x back into the random term would re-add the
# (center/scale) * dev random intercept that centring just removed.
#
# LNAADT_3 is log(AADT).  Without `log_vars` the elasticity formula treats it as
# an ordinary regressor and returns x-bar * b, the elasticity with respect to
# the LOG of traffic -- roughly ten times the elasticity with respect to traffic
# itself, which is the coefficient.  Nothing in the printed table distinguishes
# the two, so it has to be declared.
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
# `:::` because .write_truck_results_markdown() is internal and unexported: the
# unqualified call this replaces could only ever have resolved under
# devtools::load_all(), which exports everything, and would fail under library().
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
