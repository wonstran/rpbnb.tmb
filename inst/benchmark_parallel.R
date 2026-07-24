#!/usr/bin/env Rscript

# Reproducible serial-versus-parallel benchmark for the AD-safe TMB likelihood.
# Run from the package root after installing the current source tree.

library(rpbnb.tmb)

detected_cores <- parallel::detectCores(logical = FALSE)
parallel_cores <- if (is.na(detected_cores)) {
  1L
} else {
  max(1L, min(4L, detected_cores))
}

data <- utils::read.csv(file.path("data", "rwm1984_bnb.csv"))[1:500, ]
formula_1 <- docvis ~ age + hhninc + educ + female + married + kids
formula_2 <- hospvis ~ age + educ + outwork + female + self

fit_once <- function(n_cores) {
  fit_rpbnb_tmb(
    formula_1,
    formula_2,
    data = data,
    random_1 = "hhninc",
    random_2 = "educ",
    dependence = copula("frank"),
    draws = 20L,
    seed = 20240712L,
    keep = "full",
    control = rpbnb_tmb_control(iterlim = 100L, n_cores = n_cores)
  )
}

serial_elapsed <- system.time(
  fit_serial <- fit_once(1L)
)[["elapsed"]]
parallel_elapsed <- system.time(
  fit_parallel <- fit_once(parallel_cores)
)[["elapsed"]]

coef_diff <- max(abs(coef(fit_parallel) - coef(fit_serial)))
se_diff <- max(abs(fit_parallel$se - fit_serial$se))
speedup <- serial_elapsed / parallel_elapsed

# OpenMP changes floating-point reduction order. Verify the two taped
# likelihoods directly at one common parameter vector before comparing the
# nearby solutions selected by nlminb's stopping rules.
probe_par <- coef(fit_serial)
TMB::openmp(n = fit_serial$parallel$realized, DLL = "rpbnb.tmb")
serial_objective <- fit_serial$obj$fn(probe_par)
serial_gradient <- fit_serial$obj$gr(probe_par)
TMB::openmp(n = fit_parallel$parallel$realized, DLL = "rpbnb.tmb")
parallel_objective <- fit_parallel$obj$fn(probe_par)
parallel_gradient <- fit_parallel$obj$gr(probe_par)
objective_diff <- abs(parallel_objective - serial_objective)
gradient_diff <- max(abs(parallel_gradient - serial_gradient))

if (!all(is.finite(coef(fit_serial))) ||
    !all(is.finite(coef(fit_parallel))) ||
    !all(is.finite(fit_serial$se)) ||
    !all(is.finite(fit_parallel$se))) {
  stop("Benchmark produced non-finite coefficients or standard errors.")
}
if (!is.finite(objective_diff) || objective_diff > 1e-8) {
  stop(sprintf("Objective difference %.3g exceeds tolerance 1e-8.",
               objective_diff))
}
if (!is.finite(gradient_diff) || gradient_diff > 1e-8) {
  stop(sprintf("Gradient difference %.3g exceeds tolerance 1e-8.",
               gradient_diff))
}
# These are optimizer-level comparisons: nlminb may stop at nearby points when
# the reduction order perturbs the objective at machine precision.
if (!is.finite(coef_diff) || coef_diff > 1e-3) {
  stop(sprintf("Coefficient difference %.3g exceeds tolerance 1e-3.", coef_diff))
}
if (!is.finite(se_diff) || se_diff > 1e-3) {
  stop(sprintf("Standard-error difference %.3g exceeds tolerance 1e-3.", se_diff))
}

cat(sprintf(
  "Serial: requested=%d realized=%d elapsed=%.3f s\n",
  fit_serial$parallel$requested,
  fit_serial$parallel$realized,
  serial_elapsed
))
cat(sprintf(
  "Parallel: requested=%d realized=%d elapsed=%.3f s\n",
  fit_parallel$parallel$requested,
  fit_parallel$parallel$realized,
  parallel_elapsed
))
cat(sprintf("speedup=%.3fx\n", speedup))
cat(sprintf("coef_diff=%.12g\n", coef_diff))
cat(sprintf("se_diff=%.12g\n", se_diff))
cat(sprintf("objective_diff=%.12g\n", objective_diff))
cat(sprintf("gradient_diff=%.12g\n", gradient_diff))
