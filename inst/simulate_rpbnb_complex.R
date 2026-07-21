#!/usr/bin/env Rscript
# =============================================================================
# Generate a richer random-parameter bivariate NB (RP-BNB) sample dataset
# using rpbnb.tmb's simulate_rpbnb_tmb() with Famoye/Sarmanov dependence.
#
#   * 5000 observations
#   * 8 independent variables: 3 continuous + 5 dummy (0/1)
#   * a mix of FIXED and RANDOM coefficients in each equation
#   * Genuine Famoye dependence (non-zero lambda) between y1 and y2
# =============================================================================

library(rpbnb.tmb)

n <- 5000L

# ---- 1. Covariates: 3 continuous + 5 dummy ---------------------------------
set.seed(101)
covariates <- data.frame(
  x_age    = rnorm(n, 0, 1),
  x_income = rnorm(n, 0, 1),
  x_score  = rnorm(n, 0, 1),
  d_female  = rbinom(n, 1, 0.50),
  d_urban   = rbinom(n, 1, 0.60),
  d_married = rbinom(n, 1, 0.55),
  d_college = rbinom(n, 1, 0.40),
  d_smoker  = rbinom(n, 1, 0.25)
)

# ---- 2. True coefficient MEANS -----------------------------------------------
beta1 <- c("(Intercept)" = 0.50, x_age = 0.20, x_income = 0.15, x_score = -0.10,
           d_female = 0.30, d_urban = -0.25, d_married = 0.10,
           d_college = 0.20, d_smoker = 0.35)

beta2 <- c("(Intercept)" = 0.30, x_age = -0.10, x_income = 0.25, x_score = 0.15,
           d_female = -0.20, d_urban = 0.30, d_married = -0.15,
           d_college = 0.25, d_smoker = 0.10)

# ---- 3. RANDOM coefficients (continuous -> well identified) -----------------
random_1 <- list(x_age   = list(dist = "normal", sd = 0.30),
                 d_urban = list(dist = "normal", sd = 0.25))
random_2 <- list(x_income = list(dist = "normal", sd = 0.20),
                 d_college = list(dist = "normal", sd = 0.25))

# ---- 4. NB2 dispersions & Famoye dependence (positive lambda) ---------------
dispersion <- c(m1 = 0.50, m2 = 0.60)

# We want a moderate positive dependence. simulate_rpbnb_tmb() draws
# via conditional enumeration and uses a safe non-zero lambda.
# Pick lambda ~ 0.3 * max_valid to stay well inside the admissible range.
# (The exact valid bounds depend on the realized means and dispersions.)
# A pilot call with lambda = 0.3 is safe for these dispersion values.
lambda <- 0.3

# ---- 5. Simulate ------------------------------------------------------------
sim <- simulate_rpbnb_tmb(
  n          = n,
  beta1      = beta1,
  beta2      = beta2,
  random_1   = random_1,
  random_2   = random_2,
  dispersion = dispersion,
  dependence = "famoye",
  lambda     = lambda,
  covariates = covariates,
  seed       = 202
)

data <- sim$data

# ---- 6. Report --------------------------------------------------------------
cat("=== Simulated RP-BNB dataset (TMB, Famoye dependence) ===\n")
cat("Observations :", nrow(data), "\n")
cat("Variables    :", paste(setdiff(names(data), c("y1", "y2")), collapse = ", "), "\n\n")
cat("First 6 rows:\n"); print(head(data))
cat("\nOutcome summaries:\n"); print(summary(data[, c("y1", "y2")]))
cat(sprintf("\ny1: mean=%.3f var=%.3f\n", mean(data$y1), var(data$y1)))
cat(sprintf("y2: mean=%.3f var=%.3f\n", mean(data$y2), var(data$y2)))
cat(sprintf("cor(y1, y2) = %.4f\n", cor(data$y1, data$y2)))

# ---- 7. Persist data + ground-truth parameters ------------------------------
dir.create("data", showWarnings = FALSE)
out_csv <- file.path("data", "simulated_rpbnb_complex.csv")
write.csv(data, out_csv, row.names = FALSE)

truth <- list(
  beta1 = beta1, beta2 = beta2,
  random_1 = random_1, random_2 = random_2,
  dispersion = dispersion, lambda = lambda,
  random_names_1 = names(random_1),
  random_names_2 = names(random_2)
)
out_rds <- file.path("data", "simulated_rpbnb_complex_truth.rds")
saveRDS(truth, out_rds)

cat("\nSaved data  ->", out_csv, "\n")
cat("Saved truth ->", out_rds, "\n")
cat(sprintf("TRUE lambda = %.4f (Famoye dependence)\n", lambda))
