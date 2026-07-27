#!/usr/bin/env Rscript
# =============================================================================
# rpbnb.tmb -- truck crash RP-BNB, Famoye dependence, Laplace estimator.
#
# Same model and data as truck_rpbnb_diff_famoye_dense.R, which exhausts memory
# during MakeADFun() because the SML tape scales with n * draws.  Laplace tapes
# one conditional evaluation per observation and integrates the latents through
# a sparse Hessian, so the draw budget no longer binds and the default
# max_workload is sufficient -- note the absence of the max_workload override
# the SML script needs.
#
# Run from the package root:
#   library(rpbnb.tmb); source("inst/truck_rpbnb_diff_famoye_laplace.R")
# =============================================================================

library(here)
devtools::load_all()

sep <- function() cat("\n", paste(rep("=", 72), collapse = ""), "\n", sep = "")

# Under Laplace `draws` no longer sizes the tape.  It still sizes the Halton
# grid used for the frozen Famoye lambda bounds and for the post-estimation
# averaging in predict() and the marginal-effect functions, so it is not inert.
draws   <- 300L
n_cores <- 12L

data <- read.csv(here("inst\\extdata", "export_dense_all.csv"))

cat("=== RP-BNB on truck all crashes (Laplace) ===\n")
cat("Observations :", nrow(data), "\n")
cat("Cores asked  :", n_cores, "\n")

f1 <- ALL_3  ~ LNAADT_3+SR40_MI3+MPD_ME+MPD_STD+IRI_ME+RUT_L+SP50GE+ACCPNTS+SIGNAL1+NEAR_SIG+AUXLNUM++DP10_ME
f2 <- C_DISTR ~ LNAADT_3+SR40_MI3+MPD_ME+MPD_STD+RUT_9+ACCPNTS+SIGNAL1+NEAR_SIG+AUXLNUM+DP01_ME+CS_MINAB

cat("Equation 1   :", deparse(f1), "\n")
cat("Equation 2   :", deparse(f2), "\n\n")

t_fit <- system.time(
  fit <- fit_rpbnb_tmb(
    formula_1  = f1,
    formula_2  = f2,
    data       = data,
    random_1   = c("SR40_MI3", "AUXLNUM", "MPD_ME", "MPD_STD"),
    random_2   = c("SR40_MI3", "RUT_9", "MPD_ME", "MPD_STD"),
    dependence = "famoye",
    draws      = draws,
    seed       = 20240712,
    method     = "laplace",
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

sep(); cat("MODEL SUMMARY\n"); sep()
summary(fit)

sep(); cat("DEPENDENCE (sdreport)\n"); sep()
if (!is.null(fit$sdreport)) print(summary(fit$sdreport, "report"))
