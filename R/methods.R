#' @export
coef.rpbnb_tmb_fit <- function(object, ...) object$coef

#' @export
vcov.rpbnb_tmb_fit <- function(object, ...) object$vcov

#' @export
logLik.rpbnb_tmb_fit <- function(object, ...) {
  structure(object$logLik, df = object$npar, nobs = object$nobs, class = "logLik")
}

#' @export
AIC.rpbnb_tmb_fit <- function(object, ...) {
  -2 * object$logLik + 2 * object$npar
}

#' @export
BIC.rpbnb_tmb_fit <- function(object, ...) {
  -2 * object$logLik + log(object$nobs) * object$npar
}

#' @export
print.rpbnb_tmb_fit <- function(x, ...) {
  cat("rpbnb_tmb fit\n")
  cat("  Log-likelihood:", format(x$logLik, digits = 6), "\n")
  cat("  Nobs:", x$nobs, "  Npar:", x$npar, "\n")
  cat("  Dependence:", deparse(x$dependence), "\n")
  cat("\nCoefficients:\n")
  print(round(x$coef, 4))
  invisible(x)
}

#' @export
summary.rpbnb_tmb_fit <- function(object, ...) {
  cat("Summary: rpbnb_tmb fit\n")
  cat("  Log-likelihood:", format(object$logLik, digits = 6),
      "  AIC:", format(AIC(object), digits = 6),
      "  BIC:", format(BIC(object), digits = 6), "\n")
  cat("  Nobs:", object$nobs, "  Npar:", object$npar, "\n\n")

  # ---- Equation 1 coefficients (b1:*) ----
  nm <- names(object$coef)
  eq1 <- grep("^b1:", nm)
  if (length(eq1)) {
    cat("--- Equation 1 (y1) ---\n")
    tbl1 <- cbind(Estimate = object$coef[eq1],
                  `Std. Error` = object$se[eq1],
                  `z value` = object$coef[eq1] / object$se[eq1],
                  `Pr(>|z|)` = 2 * pnorm(-abs(object$coef[eq1] / object$se[eq1])))
    print(tbl1)
    cat("\n")
  }

  # ---- Equation 2 coefficients (b2:*) ----
  eq2 <- grep("^b2:", nm)
  if (length(eq2)) {
    cat("--- Equation 2 (y2) ---\n")
    tbl2 <- cbind(Estimate = object$coef[eq2],
                  `Std. Error` = object$se[eq2],
                  `z value` = object$coef[eq2] / object$se[eq2],
                  `Pr(>|z|)` = 2 * pnorm(-abs(object$coef[eq2] / object$se[eq2])))
    print(tbl2)
    cat("\n")
  }

  # ---- Random-coefficient scales (displayed as sd1/sd2) ----
  s1 <- grep("^(log_sd1|log_s1|log_w1):", nm)
  if (length(s1)) {
    cat("--- Random-coefficient SDs (equation 1) ---\n")
    tbl_s1 <- cbind(Estimate_log = object$coef[s1],
                    Estimate = exp(object$coef[s1]),
                    `Std. Error` = object$se[s1])
    rownames(tbl_s1) <- gsub("^log_sd1:|^log_s1:|^log_w1:", "sd1:", nm[s1])
    print(tbl_s1)
    cat("\n")
  }

  # ---- Random-coefficient scales (displayed as sd1/sd2) ----
  s2 <- grep("^(log_sd2|log_s2|log_w2):", nm)
  if (length(s2)) {
    cat("--- Random-coefficient SDs (equation 2) ---\n")
    tbl_s2 <- cbind(Estimate_log = object$coef[s2],
                    Estimate = exp(object$coef[s2]),
                    `Std. Error` = object$se[s2])
    rownames(tbl_s2) <- gsub("^log_sd2:|^log_s2:|^log_w2:", "sd2:", nm[s2])
    print(tbl_s2)
    cat("\n")
  }

  # ---- Dispersion parameters (natural scale from fit object) ----
  cat("--- Dispersion (m1, m2) ---\n")
  disp <- data.frame(
    Parameter = c("m1", "m2"),
    Estimate  = c(object$m1, object$m2),
    row.names = NULL
  )
  print(disp, digits = 4)
  cat("\n")

  # ---- Dependence parameter (from sdreport) ----
  cat("--- Dependence ---\n")
  sdr <- object$sdreport
  dep <- object$dependence
  dep_name <- "z_dep"
  if (inherits(dep, "rpbnb_copula")) {
    dep_name <- switch(dep$family, frank = "theta", normal = "rho", kimeldorf = "theta")
  } else if (identical(dep, "famoye")) {
    dep_name <- "lam"
  }
  if (!is.null(sdr)) {
    sdr_sum <- suppressWarnings(try(summary(sdr, "report"), silent = TRUE))
    if (!inherits(sdr_sum, "try-error") && dep_name %in% rownames(sdr_sum)) {
      dep_val <- sdr_sum[dep_name, "Estimate"]
      dep_se  <- sdr_sum[dep_name, "Std. Error"]
      cat("  ", dep_name, " =", format(dep_val, digits = 4),
          "  Std. Error =", format(dep_se, digits = 4), "\n")
    } else {
      cat("  (dependence parameter not in sdreport)\n")
    }
  }
  cat("\n")

  invisible(object)
}

#' @export
predict.rpbnb_tmb_fit <- function(object, ...) {
  coefs <- object$coef
  k1 <- sum(grepl("^b1:", names(coefs)))
  k2 <- sum(grepl("^b2:", names(coefs)))
  beta1 <- coefs[seq_len(k1)]
  beta2 <- coefs[seq_len(k2) + k1]
  # Fitted linear predictors — requires model frame; without newdata,
  # return the coefficients as fitted parameter vector.
  structure(coefs, class = "rpbnb_tmb_pred")
}
