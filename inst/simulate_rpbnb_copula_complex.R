#!/usr/bin/env Rscript
# =============================================================================
# Generate a rich copula random-parameter bivariate NB (RP-BNB) sample dataset
# using rpbnb.tmb's simulate_rpbnb_tmb() with Gaussian copula dependence.
#
#   * 5000 observations
#   * 8 independent variables: 3 continuous + 5 dummy
#   * a mix of FIXED and RANDOM coefficients in each equation
#   * Gaussian copula dependence (rho = 0.5) between y1 and y2
#
# NOTE: random coefficients sit on CONTINUOUS regressors only -- the copula RP
# estimator is weakly identified when a random coefficient sits on a 0/1 dummy
# (NB dispersion vs random-coefficient scale trade-off).
# =============================================================================

library(rpbnb.tmb)

n <- 5000L

# ---- 1. Covariates: 3 continuous + 5 dummy ---------------------------------
set.seed(606)
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

# ---- 3. RANDOM coefficients (one per equation, CONTINUOUS -> identified) -----
random_1 <- list(x_age    = list(dist = "normal", sd = 0.30))
random_2 <- list(x_income = list(dist = "normal", sd = 0.25))

# ---- 4. NB2 dispersions & copula dependence ---------------------------------
dispersion <- c(m1 = 0.50, m2 = 0.60)
cop <- copula("normal", par = 0.5)   # Gaussian copula, rho = 0.5

# ---- 5. Simulate ------------------------------------------------------------
sim <- simulate_rpbnb_tmb(
  n          = n,
  beta1      = beta1,
  beta2      = beta2,
  random_1   = random_1,
  random_2   = random_2,
  dispersion = dispersion,
  dependence = cop,
  covariates = covariates,
  seed       = 707
)

data <- sim$data

# ---- 6. Report --------------------------------------------------------------
cat("=== Complex copula RP-BNB dataset (TMB) ===\n")
cat("Observations :", nrow(data), "\n")
cat("Variables    :", paste(setdiff(names(data), c("y1", "y2")), collapse = ", "), "\n")
cat("Copula       : Gaussian (rho = 0.5)\n")
cat("Random coefs : eq1 {x_age}, eq2 {x_income} (continuous -> identified)\n\n")

cat("First 6 rows:\n"); print(head(data))
cat("\nOutcome summaries:\n"); print(summary(data[, c("y1", "y2")]))
cat(sprintf("\ny1: mean=%.3f var=%.3f\n", mean(data$y1), var(data$y1)))
cat(sprintf("y2: mean=%.3f var=%.3f\n", mean(data$y2), var(data$y2)))
cat(sprintf("Pearson(y1, y2)  = %.4f\n", cor(data$y1, data$y2)))
cat(sprintf("Spearman(y1, y2) = %.4f (positive copula dependence built in)\n",
            cor(data$y1, data$y2, method = "spearman")))

# ---- 7. Persist data + ground-truth parameters ------------------------------
dir.create("data", showWarnings = FALSE)
out_csv <- file.path("data", "simulated_rpbnb_copula_complex.csv")
write.csv(data, out_csv, row.names = FALSE)

truth <- list(
  beta1 = beta1, beta2 = beta2,
  random_1 = random_1, random_2 = random_2,
  dispersion = dispersion,
  copula = "normal", theta = 0.5,
  random_names_1 = names(random_1),
  random_names_2 = names(random_2)
)
out_rds <- file.path("data", "simulated_rpbnb_copula_complex_truth.rds")
saveRDS(truth, out_rds)

cat("\nSaved data  ->", out_csv, "\n")
cat("Saved truth ->", out_rds, "\n")
