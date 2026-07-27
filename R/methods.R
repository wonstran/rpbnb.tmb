#' @export
coef.rpbnb_tmb_fit <- function(object, ...) object$coef

#' @export
vcov.rpbnb_tmb_fit <- function(object, ...) {
  if (!is.null(object$vcov)) return(object$vcov)
  nm <- names(object$coef)
  out <- matrix(NA_real_, length(nm), length(nm), dimnames = list(nm, nm))
  if (!is.null(object$vcov_diag)) diag(out) <- object$vcov_diag
  out
}

#' @export
logLik.rpbnb_tmb_fit <- function(object, ...) {
  structure(object$logLik, df = object$npar, nobs = object$nobs, class = "logLik")
}

#' @export
print.rpbnb_tmb_fit <- function(x, ...) {
  cat("rpbnb_tmb fit\n")
  cat("  Log-likelihood:", format(x$logLik, digits = 6), "\n")
  cat("  Nobs:", x$nobs, "  Npar:", x$npar, "\n")
  if (!is.null(x$parallel$realized)) {
    cat(sprintf("  TMB threads: %d\n", as.integer(x$parallel$realized)))
  }
  cat("  Estimator:", if (is.null(x$method)) "sml" else x$method, "\n")
  cat("  Dependence:", deparse(x$dependence), "\n")
  cat("\nCoefficients:\n")
  print(round(x$coef, 4))
  invisible(x)
}

# Significance stars helper
.signif_stars <- function(p) {
  symnum(p, corr = FALSE, na = FALSE,
         cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
         symbols = c("***", "**", "*", ".", " "))
}

# Print a coefficient table: controlled decimal places, no quotes on column names.
# digits = 4 (default): 4 decimal places. digits < 0: full precision.
.print_tbl <- function(x, digits = 4) {
  x <- as.data.frame(x)
  if (!is.null(digits) && digits >= 0) {
    fmt <- sprintf("%%.%df", digits)
    for (j in seq_len(ncol(x))) {
      if (is.numeric(x[[j]])) x[[j]] <- sprintf(fmt, x[[j]])
    }
  }
  print(x, print.gap = 3)
}

#' Summarize a fitted rpbnb.tmb model
#'
#' @param object A fitted model object of class rpbnb_tmb_fit.
#' @param digits Number of decimal places for numeric columns in coefficient
#'   tables. Default 4. Use a negative value (e.g. -1) for full precision.
#' @param ... Not used.
#' @return The fitted object, invisibly.
#' @export
summary.rpbnb_tmb_fit <- function(object, digits = 4L, ...) {
  cat("Summary: rpbnb_tmb fit\n")
  cat("  Log-likelihood:", format(object$logLik, digits = 6),
      "  AIC:", format(AIC(object), digits = 6),
      "  BIC:", format(BIC(object), digits = 6), "\n")
  iterations <- object$optimizer$iterations
  cat("  Nobs:", object$nobs, "  Npar:", object$npar,
      if (!is.null(iterations)) c("  Iterations:", iterations), "\n")
  cat("  Estimator:", if (is.null(object$method)) "sml" else object$method,
      "\n\n")

  # ---- Equation 1 coefficients (b1:*) ----
  nm <- names(object$coef)
  eq1 <- grep("^b1:", nm)
  if (length(eq1)) {
    cat("--- Equation 1 (y1) ---\n")
    p1 <- 2 * pnorm(-abs(object$coef[eq1] / object$se[eq1]))
    tbl1 <- data.frame(Estimate = object$coef[eq1],
                       `Std. Error` = object$se[eq1],
                       `z value` = object$coef[eq1] / object$se[eq1],
                       `Pr(>|z|)` = p1,
                       Signif = .signif_stars(p1),
                       row.names = nm[eq1], check.names = FALSE)
    .print_tbl(tbl1, digits)
    cat("\n")
  }

  # ---- Equation 2 coefficients (b2:*) ----
  eq2 <- grep("^b2:", nm)
  if (length(eq2)) {
    cat("--- Equation 2 (y2) ---\n")
    p2 <- 2 * pnorm(-abs(object$coef[eq2] / object$se[eq2]))
    tbl2 <- data.frame(Estimate = object$coef[eq2],
                       `Std. Error` = object$se[eq2],
                       `z value` = object$coef[eq2] / object$se[eq2],
                       `Pr(>|z|)` = p2,
                       Signif = .signif_stars(p2),
                       row.names = nm[eq2], check.names = FALSE)
    .print_tbl(tbl2, digits)
    cat("\n")
  }

  # ---- Random-coefficient scales (displayed as sd1/sd2) ----
  s1 <- grep("^(log_sd1|log_s1|log_w1):", nm)
  if (length(s1)) {
    cat("--- Random-coefficient SDs (equation 1) ---\n")
    p_s1 <- 2 * pnorm(-abs(object$coef[s1] / object$se[s1]))
    tbl_s1 <- data.frame(Estimate_log = object$coef[s1],
                         Estimate = exp(object$coef[s1]),
                         `Std. Error` = object$se[s1],
                         `Pr(>|z|)` = p_s1,
                         Signif = .signif_stars(p_s1),
                         row.names = gsub("^log_sd1:|^log_s1:|^log_w1:", "sd1:", nm[s1]),
                         check.names = FALSE)
    .print_tbl(tbl_s1, digits)
    cat("\n")
  }

  # ---- Random-coefficient scales (displayed as sd1/sd2) ----
  s2 <- grep("^(log_sd2|log_s2|log_w2):", nm)
  if (length(s2)) {
    cat("--- Random-coefficient SDs (equation 2) ---\n")
    p_s2 <- 2 * pnorm(-abs(object$coef[s2] / object$se[s2]))
    tbl_s2 <- data.frame(Estimate_log = object$coef[s2],
                         Estimate = exp(object$coef[s2]),
                         `Std. Error` = object$se[s2],
                         `Pr(>|z|)` = p_s2,
                         Signif = .signif_stars(p_s2),
                         row.names = gsub("^log_sd2:|^log_s2:|^log_w2:", "sd2:", nm[s2]),
                         check.names = FALSE)
    .print_tbl(tbl_s2, digits)
    cat("\n")
  }

  # ---- Dispersion parameters (natural scale from fit object) ----
  cat("--- Dispersion (m1, m2) ---\n")
  disp <- data.frame(
    Parameter = c("m1", "m2"),
    Estimate  = c(object$m1, object$m2),
    row.names = NULL
  )
  .print_tbl(disp, digits)
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
      cat("  ", dep_name, " =", format(dep_val, digits = if (digits >= 0) digits else 6),
          "  Std. Error =", format(dep_se, digits = if (digits >= 0) digits else 6), "\n")
    } else {
      cat("  (dependence parameter not in sdreport)\n")
    }
  }
  cat("\n")

  invisible(object)
}

#' Predict from a fitted bivariate count model
#'
#' Predictions integrate over the retained simulation draws when random
#' coefficients are present. Link predictions are the log of the integrated
#' response mean, so \code{exp(predict(fit, type = "link"))} equals response
#' predictions.
#'
#' @param object A fitted \code{rpbnb_tmb_fit} object.
#' @param newdata Optional data frame. If omitted, the retained fitting design
#'   is used; compact fits require \code{newdata}.
#' @param type Either \code{"response"} or \code{"link"}.
#' @param which Return both margins or only \code{"y1"} or \code{"y2"}.
#' @param ... Reserved for future use.
#' @return A two-column numeric matrix for \code{which = "both"}, otherwise a
#'   numeric vector.
#' @export
predict.rpbnb_tmb_fit <- function(
  object, newdata = NULL,
  type = c("response", "link"),
  which = c("both", "y1", "y2"), ...
) {
  type <- match.arg(type)
  which <- match.arg(which)
  equations <- if (which == "both") 1:2 else if (which == "y1") 1L else 2L

  predictions <- lapply(equations, function(eq) {
    X <- .prediction_design(object, eq, newdata)
    response <- .prediction_margin(object, eq, X)
    if (type == "link") log(response) else response
  })

  if (which != "both") return(predictions[[1L]])
  out <- do.call(cbind, predictions)
  colnames(out) <- c("y1", "y2")
  out
}

.prediction_design <- function(object, eq, newdata) {
  if (is.null(newdata)) {
    X <- object[[paste0("X", eq)]]
    if (is.null(X)) {
      stop(
        'This compact fit does not retain fitting rows; supply `newdata`.',
        call. = FALSE
      )
    }
    return(X)
  }
  if (!is.data.frame(newdata)) {
    stop("`newdata` must be a data frame.", call. = FALSE)
  }
  meta <- object$model_meta[[paste0("eq", eq)]]
  frame <- stats::model.frame(
    meta$terms, data = newdata, na.action = stats::na.pass,
    xlev = meta$xlevels
  )
  X <- stats::model.matrix(
    meta$terms, data = frame, contrasts.arg = meta$contrasts
  )
  expected <- sub(paste0("^b", eq, ":"), "", grep(
    paste0("^b", eq, ":"), names(object$coef), value = TRUE
  ))
  if (!identical(colnames(X), expected)) {
    stop("`newdata` produced an incompatible design matrix.", call. = FALSE)
  }
  X
}

.prediction_margin <- function(object, eq, X) {
  cn <- colnames(X)
  beta_names <- paste0("b", eq, ":", cn)
  beta <- object$coef[beta_names]
  rand_idx <- object[[paste0("rand_idx", eq)]]

  if (!length(rand_idx)) {
    return(exp(pmin(pmax(as.vector(X %*% beta), -700), 700)))
  }
  if (is.null(object$rp_meta)) {
    stop(
      "Random-coefficient prediction requires retained simulation draws.",
      call. = FALSE
    )
  }

  Z <- object$rp_meta[[paste0("Z", eq)]]
  dist <- object$rp_meta[[paste0("dist", eq)]]
  sign <- object$rp_meta[[paste0("sign", eq)]]
  scale_labels <- vapply(
    dist, function(name) rand_dist_registry[[name]]$scale_label,
    character(1)
  )
  scale_names <- paste0(scale_labels, eq, ":", cn[rand_idx])
  scales <- exp(object$coef[scale_names])
  realized <- rand_realize(
    Z, dist, sign, beta[rand_idx], scales
  )
  xb <- as.vector(X %*% beta)
  xr <- X[, rand_idx, drop = FALSE]
  .draw_mean_exp(xb, xr, realized$dev)
}
