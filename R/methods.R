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
  coef_table <- cbind(Estimate = object$coef,
                      `Std. Error` = object$se,
                      `z value` = object$coef / object$se,
                      `Pr(>|z|)` = 2 * pnorm(-abs(object$coef / object$se)))
  cat("Summary: rpbnb_tmb fit\n")
  cat("  Log-likelihood:", format(object$logLik, digits = 6), "\n")
  cat("  AIC:", format(AIC(object), digits = 6), " BIC:", format(BIC(object), digits = 6), "\n")
  cat("  Dependence:", deparse(object$dependence), "\n")
  cat("\nCoefficients:\n")
  print(coef_table)
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
