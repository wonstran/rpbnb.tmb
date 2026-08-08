# NOT library(rpbnb.tmb): the installed 0.3.4 predates the `scaling`/`log_vars`
# arguments used below, and both diagnostic functions take `...`, so under the
# installed build those arguments are SILENTLY SWALLOWED -- the script would
# print standardized numbers under an "original units" heading and only fail
# later, at the exporter.  The source tree carries the same version number, so
# the version string cannot tell the two apart.  Switch back to library() only
# after a version bump and devtools::install().
devtools::load_all("C:\\Users\\zwang9\\repos\\rpbnb.tmb")

sep <- function() cat("\n", paste(rep("=", 72), collapse = ""), "\n", sep = "")
setwd("C:\\Users\\zwang9\\repos\\rpbnb.tmb")

# detectCores() returns 24 on this box (31.5 GiB RAM).  8 threads, not more:
# with tape.parallel off TMB still records one tape per parallel region, and
# more regions buy eval speed, not construction speed.
#
# 300 draws is affordable ONLY because gaussian_cell_prob() is checkpointed
# (REGISTER_ATOMIC in src/rpbnb.tmb.cpp): measured on this model at 8 threads,
# peak working set is ~0.4 + 0.050 * draws GiB (0.56 GiB at 10 draws, 2.6 GiB
# at 50), so 300 draws peaks around 13 GiB.  Before the checkpoint the same
# kernel cost ~0.57 GiB of peak PER DRAW -- ~177 GiB at 300 draws, the
# std::bad_alloc this comment replaces -- because the taped quadrature ran
# ~112 pnorm/dnorm nodes per observation-draw and this data also tapes NB CDF
# sums out to counts of ~242.  TAPE_CALIBRATION's per-unit constant is ~22x
# too optimistic for exactly that reason (its fixture uses Poisson(2) counts),
# so do not trust the workload guard's arithmetic on this data; trust the
# probe: run MakeADFun + one fn()/gr() at two small draw counts and
# extrapolate linearly before raising draws further.
n_cores <- 20L
draws <- 500L
# "famoye" and "independence" are plain strings; only frank/normal/kimeldorf
# go through copula().  The label below handles both shapes.
dependence <- "famoye" #copula("frank") #copula("kimeldorf")
est_method <- "laplace"

# setwd() is at the project root; file.path() builds the platform-native
# separator, so this resolves on both Windows and POSIX.
data <- read.csv(file.path("inst", "extdata", "export_open_all.csv"))

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

cat("=== RP-BNB on truck all crashes (", est_method, ") ===\n", sep = "")
dep_label <- if (inherits(dependence, "rpbnb_copula")) {
  paste0("copula(\"", dependence$family, "\")")
} else {
  dependence
}
cat("Dependence   :", dep_label, "\n")
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
    # The 40-draw fit estimated equation 2's SR40_MI3 random SD at exp(-12.4)
    # ~= 0 with SE 156: the data carry no detectable slope heterogeneity in
    # C_HV, and keeping the parameter just burns a Halton dimension and a row
    # of the Hessian on noise.  Dropping it also makes the same draw count
    # integrate a 2-D latent instead of a 3-D one.
    random_2   = c("SR40_MI3"),
    dependence = dependence,
    seed       = 20240712,
    method     = est_method,
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
dir.create("results", recursive = TRUE, showWarnings = FALSE)
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

# ---- Coefficients in original units (display) -------------------------------
# The fit above is on standardized predictors, so each continuous coefficient
# is per standard deviation.  Because the standardization is an affine column
# transform, the original-unit coefficients are exact and need no refit:
# slopes divide by the scale, the intercept absorbs the centring shift
# -sum_j (c_j / s_j) * b_j, and binary 0/1 coefficients are unchanged.  SEs are
# delta-method on the full covariance, so the intercept's cross-covariances
# with the slopes count.  Random-coefficient SDs rescale the same way as
# slopes and are shown per equation in the same format as `summary(fit)`:
# Estimate_log is the log of the original-unit SD (derived, = log_sd - log
# scale), its Std. Error is the fitted log-scale SE (unchanged by the constant
# shift), and Pr(>|z|)/Signif follow.  This section is for display only; the
# fitted design itself stays
# standardized (see the note at the end of the original-units section).  The
# tables use the same fixed-decimal format and significance stars as
# `summary(fit)` (.print_tbl), so the two summaries line up.  The captured
# output is passed to the markdown exporter below as `coef_orig_units`.
sep(); cat("COEFFICIENTS IN ORIGINAL UNITS\n"); sep()
sc_mat <- do.call(rbind, scaling)
coef_orig_units <- function(eq) {
  nm <- names(fit$coef)
  idx <- grep(paste0("^b", eq, ":"), nm)
  var <- sub(paste0("^b", eq, ":"), "", nm[idx])
  A <- diag(length(idx))
  dimnames(A) <- list(nm[idx], nm[idx])
  cont <- which(var %in% continuous_vars & var != "(Intercept)")
  A[cont, cont] <- diag(1 / sc_mat[var[cont], "scale"])
  A[var == "(Intercept)", cont] <-
    -sc_mat[var[cont], "center"] / sc_mat[var[cont], "scale"]
  V <- fit$vcov
  if (is.null(rownames(V))) dimnames(V) <- list(names(fit$coef), names(fit$coef))
  b_orig  <- as.vector(A %*% fit$coef[idx])
  se_orig <- sqrt(pmax(diag(A %*% V[nm[idx], nm[idx]] %*% t(A)), 0))
  z_orig  <- b_orig / se_orig
  p_orig  <- 2 * pnorm(-abs(z_orig))
  data.frame(
    Estimate = b_orig,
    `Std. Error` = se_orig,
    `z value` = z_orig,
    `Pr(>|z|)` = p_orig,
    Signif = rpbnb.tmb:::.signif_stars(p_orig),
    row.names = nm[idx], check.names = FALSE
  )
}
coef_orig_units_output <- capture.output({
  cat("\n--- Equation 1 (y1) ---\n")
  rpbnb.tmb:::.print_tbl(coef_orig_units(1L))
  cat("\n--- Equation 2 (y2) ---\n")
  rpbnb.tmb:::.print_tbl(coef_orig_units(2L))
  # Random-coefficient SDs, per equation and in the model summary's format.
  cat("\n--- Random-coefficient SDs (equation 1) ---\n")
  sd_orig_units <- function(prefix) {
    nm <- names(fit$coef)
    idx <- grep(paste0("^", prefix, ":"), nm)
    var <- sub(paste0("^", prefix, ":"), "", nm[idx])
    log_orig <- fit$coef[idx] - log(sc_mat[var, "scale"])
    pval <- 2 * pnorm(-abs(log_orig / fit$se[idx]))
    data.frame(
      Estimate_log = log_orig,
      Estimate = exp(log_orig),
      `Std. Error` = fit$se[idx],
      `Pr(>|z|)` = pval,
      Signif = rpbnb.tmb:::.signif_stars(pval),
      row.names = sub("^log_", "", nm[idx]), check.names = FALSE
    )
  }
  rpbnb.tmb:::.print_tbl(sd_orig_units("log_sd1"))
  cat("\n--- Random-coefficient SDs (equation 2) ---\n")
  if (any(grepl("^log_sd2:", names(fit$coef)))) {
    rpbnb.tmb:::.print_tbl(sd_orig_units("log_sd2"))
  } else {
    cat("  (no random coefficients in equation 2)\n")
  }
})
cat(coef_orig_units_output, sep = "\n")
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
  scaling = scaling,
  coef_orig_units = coef_orig_units_output
)
cat("Results written to:", results_path, "\n")
