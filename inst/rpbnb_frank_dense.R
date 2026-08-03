library(rpbnb.tmb)

sep <- function() cat("\n", paste(rep("=", 72), collapse = ""), "\n", sep = "")
setwd("/home/wonstran/repos/truck")
n_cores <- parallel::detectCores()
draws <- 500L
dependence <- copula("frank")
optimizer <- "laplace"

# Path components stay separate: "inst\\extdata" is one directory named
# inst\extdata on POSIX, so the backslash form only resolves on Windows.
data <- read.csv(file.path("data", "export_dense_all.csv"))

cat("=== RP-BNB on truck all crashes (Laplace) ===\n")
cat("Dependence   : copula(\"", dependence$family, "\")\n", sep = "")
cat("Cores asked  :", n_cores, "\n")
cat("Optimizer    :", optimizer, "\n")
cat("Draws (if SML):", draws, "\n\n")

f1 <- ALL_3  ~ LNAADT_3+SR40_MI3+MPD_ME+MPD_STD+IRI_ME+RUT_L+SP50GE+ACCPNTS+SIGNAL1+NEAR_SIG+AUXLNUM++DP10_ME
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
    data       = data,
    random_1   = r1,
    random_2   = r2,
    dependence = dependence,
    seed       = 20240712,
    method     = optimizer,
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

# ---- Export requested results ----------------------------------------------
results_path <- rpbnb.tmb:::.write_truck_results_markdown(
  model_summary = model_summary_output,
  marginal_effects = marginal_effects_output,
  elasticities = elasticities_output,
  dependence = fit$dependence,
  method = fit$method
)
cat("Results written to:", results_path, "\n")
