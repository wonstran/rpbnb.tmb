#' Marginal effects for a rpbnb_tmb model
#'
#' Computes average marginal effects (AME) or marginal effects at the mean (MEM)
#' for a fitted random-parameter bivariate negative binomial model.  For
#' continuous covariates the effect is \eqn{\partial E[Y]/\partial x_j}; for
#' binary covariates it is the discrete difference
#' \eqn{E[Y|x_j = 1] - E[Y|x_j = 0]}.  When a coefficient is random the
#' per-draw realized coefficients are used to form the Monte-Carlo integrated
#' mean.
#'
#' @param fit An object of class \code{rpbnb_tmb_fit}.
#' @param which Which margin: \code{"y1"}, \code{"y2"}, or \code{"both"}.
#' @param type \code{"AME"} (average over the sample) or \code{"MEM"} (at the
#'   sample mean of covariates).
#' @param vars Optional character vector of variable names to restrict output.
#' @param include_intercept Logical; include the intercept term in output.
#' @param digits Number of decimal places for printed table entries.
#' @param ... Not used.
#' @return A data frame (single margin) or named list of two data frames
#'   (\code{"both"}).
#' @export
rpbnb_tmb_marginal_effects <- function(fit,
                                       which = c("y1", "y2", "both"),
                                       type  = c("AME", "MEM"),
                                       vars  = NULL,
                                       include_intercept = FALSE,
                                       digits = 4L, ...) {
  which <- match.arg(which)
  type  <- match.arg(type)
  res <- if (which == "both") {
    list(y1 = .me_one_eq(fit, 1L, type, vars, include_intercept, digits),
         y2 = .me_one_eq(fit, 2L, type, vars, include_intercept, digits))
  } else {
    .me_one_eq(fit, if (which == "y1") 1L else 2L,
               type, vars, include_intercept, digits)
  }
  res
}

#' Elasticities for a rpbnb_tmb model
#'
#' Continuous elasticities \eqn{x_j / E[Y] \cdot \partial E[Y]/\partial x_j}
#' and binary semi-elasticities
#' \eqn{E[Y|x_j = 1] / E[Y|x_j = 0] - 1}.  Built on the same Monte-Carlo
#' integrated mean used by \code{\link{rpbnb_tmb_marginal_effects}}.
#'
#' @inheritParams rpbnb_tmb_marginal_effects
#' @return A data frame or named list of data frames.
#' @export
rpbnb_tmb_elasticities <- function(fit,
                                   which = c("y1", "y2", "both"),
                                   type  = c("AME", "MEM"),
                                   vars  = NULL,
                                   include_intercept = FALSE,
                                   digits = 4L, ...) {
  which <- match.arg(which)
  type  <- match.arg(type)
  res <- if (which == "both") {
    list(y1 = .el_one_eq(fit, 1L, type, vars, include_intercept, digits),
         y2 = .el_one_eq(fit, 2L, type, vars, include_intercept, digits))
  } else {
    .el_one_eq(fit, if (which == "y1") 1L else 2L,
               type, vars, include_intercept, digits)
  }
  res
}

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

#' Compute marginal effects for one equation
#' @keywords internal
#' @noRd
.me_one_eq <- function(fit, eq, type, vars, include_intercept, digits) {
  .diag_one_eq(fit, eq, type, vars, include_intercept, digits, quantity = "me")
}

#' Compute elasticities for one equation
#' @keywords internal
#' @noRd
.el_one_eq <- function(fit, eq, type, vars, include_intercept, digits) {
  .diag_one_eq(fit, eq, type, vars, include_intercept, digits, quantity = "elas")
}

#' Core diagnostic engine for one equation
#' @keywords internal
#' @noRd
.diag_one_eq <- function(fit, eq, type, vars, include_intercept, digits,
                         quantity = c("me", "elas")) {
  quantity <- match.arg(quantity)

  # ---- Extract design & metadata ----
  X     <- if (eq == 1L) fit$X1 else fit$X2
  cn    <- colnames(X)
  n     <- nrow(X)
  coef_prefix <- paste0("b", eq, ":")
  b     <- fit$coef[grep(paste0("^", coef_prefix), names(fit$coef))]
  names(b) <- cn

  rand_idx <- if (eq == 1L) fit$rand_idx1 else fit$rand_idx2
  nr     <- length(rand_idx)

  # Variable selection
  sel <- if (is.null(vars)) seq_along(cn) else {
    m <- match(vars, cn)
    if (anyNA(m)) stop("Unknown variable(s) in equation ", eq, ": ",
                       paste(vars[is.na(m)], collapse = ", "))
    m
  }
  if (!include_intercept) sel <- sel[cn[sel] != "(Intercept)"]
  if (!length(sel)) stop("No variables selected (intercept excluded).")

  # Binary detection
  is_bin <- vapply(sel, function(j) all(X[, j] %in% c(0, 1)), logical(1))

  # ---- Design for AME / MEM ----
  X_use <- if (type == "MEM") matrix(colMeans(X), 1,
                                     dimnames = list(NULL, cn)) else X
  xu    <- nrow(X_use)

  # ---- Random coefficient draws ----
  if (nr > 0) {
    Z     <- fit$rp_meta[[paste0("Z", eq)]]
    dist  <- fit$rp_meta[[paste0("dist", eq)]]
    sign  <- fit$rp_meta[[paste0("sign", eq)]]
    scale_key <- c("log_sd", "log_s", "log_w", "log_w")
    scale_names <- paste0(scale_key[dist + 1], eq, ":", cn[rand_idx])
    scales <- exp(fit$coef[scale_names])
    rr    <- rand_realize(Z, dist, sign, b[rand_idx], scales)
    dev   <- rr$dev
    cmat  <- rr$coef        # n_draws x n_rand
    R     <- nrow(Z)
  } else {
    R     <- 1L
    dev   <- matrix(0, 1, 0)
    cmat  <- NULL
    scale_names <- character(0)
  }

  # ---- Compute point estimates ----
  xb   <- as.vector(X_use %*% b)
  xr   <- if (nr > 0) X_use[, rand_idx, drop = FALSE] else NULL

  # Pre-compute per-draw mu matrices for the chosen design
  mu_mat <- matrix(0, xu, R)
  for (r in seq_len(R)) {
    eta <- xb + if (nr > 0) as.vector(xr %*% dev[r, ]) else 0
    mu_mat[, r] <- exp(pmin(pmax(eta, -700), 700))
  }
  mu_bar <- rowMeans(mu_mat)   # x u x 1

  res <- vapply(seq_along(sel), function(m) {
    j <- sel[m]; rj <- match(j, rand_idx)
    if (is_bin[m]) {
      # Discrete difference
      X0 <- X_use; X0[, j] <- 0
      X1 <- X_use; X1[, j] <- 1
      xb0 <- as.vector(X0 %*% b)
      xb1 <- as.vector(X1 %*% b)
      xr0 <- if (nr > 0) X0[, rand_idx, drop = FALSE] else NULL
      xr1 <- if (nr > 0) X1[, rand_idx, drop = FALSE] else NULL
      mu0 <- rowMeans(vapply(seq_len(R), function(r) {
        eta <- if (nr > 0) xb0 + as.vector(xr0 %*% dev[r, ]) else xb0
        exp(pmin(pmax(eta, -700), 700))
      }, numeric(xu)))
      mu1 <- rowMeans(vapply(seq_len(R), function(r) {
        eta <- if (nr > 0) xb1 + as.vector(xr1 %*% dev[r, ]) else xb1
        exp(pmin(pmax(eta, -700), 700))
      }, numeric(xu)))
      if (quantity == "me") mean(mu1 - mu0) else mean(mu1 / mu0 - 1)
    } else {
      if (!is.na(rj)) {
        # Random continuous: average of coef_rj * mu_r over draws
        dmu <- rowMeans(vapply(seq_len(R), function(r) {
          eta <- xb + as.vector(xr %*% dev[r, ])
          cmat[r, rj] * exp(pmin(pmax(eta, -700), 700))
        }, numeric(xu)))
      } else {
        # Fixed continuous: b_j * mu
        dmu <- b[j] * mu_bar
      }
      if (quantity == "me") mean(dmu) else mean(X_use[, j] * dmu / mu_bar)
    }
  }, numeric(1))

  # ---- Delta-method standard errors ----
  theta_names <- c(paste0(coef_prefix, cn), scale_names)
  theta_hat   <- fit$coef[theta_names]
  V <- fit$vcov[theta_names, theta_names, drop = FALSE]

  se <- rep(NA_real_, length(res))
  if (all(is.finite(V))) {
    G <- numDeriv::jacobian(.estimand_t, theta_hat,
                            fit = fit, eq = eq, sel = sel, cn = cn,
                            is_bin = is_bin, quantity = quantity,
                            X_use = X_use, rand_idx = rand_idx,
                            R = R, dev0 = dev, cmat0 = cmat,
                            xu = xu, nr = nr, coef_prefix = coef_prefix)
    for (m in seq_along(res)) {
      g <- G[m, ]
      se[m] <- sqrt(as.numeric(t(g) %*% V %*% g))
    }
  }

  # ---- Assemble output ----
  tab <- data.frame(
    Name = cn[sel],
    Estimate = res,
    `Std. Error` = se,
    `z value` = res / se,
    `Pr(>|z|)` = 2 * pnorm(-abs(res / se)),
    Type = ifelse(is_bin, if (quantity == "me") "binary" else "semi-elas",
                  ifelse(is_bin, "binary", "continuous")),
    row.names = NULL, check.names = FALSE
  )

  label <- if (quantity == "me") "Marginal effects" else "Elasticities"
  cat(sprintf("\n--- %s (equation %d, %s) ---\n", label, eq, type))
  rnd_cols <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  tab_print <- tab
  for (cc in rnd_cols) if (is.numeric(tab_print[[cc]])) tab_print[[cc]] <- round(tab_print[[cc]], digits)
  print(tab_print, row.names = FALSE, print.gap = 3)
  cat("\n")

  invisible(tab)
}

#' Estimand vector for one equation (used in delta-method Jacobian)
#'
#' Returns the same numeric vector as the point-estimate computation above,
#' parameterised on \code{theta} (mean + log-scale parameters).  Called with
#' zero-argument count for \code{numDeriv::jacobian}.
#' @keywords internal
#' @noRd
.estimand_t <- function(theta, fit, eq, sel, cn, is_bin, quantity,
                        X_use, rand_idx, R, dev0, cmat0, xu, nr,
                        coef_prefix) {
  n_par <- length(cn)
  b_t <- theta[seq_len(n_par)]; names(b_t) <- cn

  nr_t <- length(rand_idx)
  if (nr_t > 0) {
    sc_names <- names(theta)[-seq_len(n_par)]
    scales_t <- exp(theta[sc_names])
    dist <- fit$rp_meta[[paste0("dist", eq)]]
    sign <- fit$rp_meta[[paste0("sign", eq)]]
    Z    <- fit$rp_meta[[paste0("Z", eq)]]
    rr   <- rand_realize(Z, dist, sign, b_t[rand_idx], scales_t)
    dev_t  <- rr$dev
    cmat_t <- rr$coef
    R_t    <- nrow(Z)
  } else {
    dev_t  <- matrix(0, 1, 0)
    cmat_t <- NULL
    R_t    <- 1L
  }

  xb_t <- as.vector(X_use %*% b_t)
  xr_t <- if (nr_t > 0) X_use[, rand_idx, drop = FALSE] else NULL

  mu_mat_t <- matrix(0, xu, R_t)
  for (r in seq_len(R_t)) {
    eta <- xb_t + if (nr_t > 0) as.vector(xr_t %*% dev_t[r, ]) else 0
    mu_mat_t[, r] <- exp(pmin(pmax(eta, -700), 700))
  }
  mu_bar_t <- rowMeans(mu_mat_t)

  vapply(seq_along(sel), function(m) {
    j <- sel[m]; rj <- match(j, rand_idx)
    if (is_bin[m]) {
      X0 <- X_use; X0[, j] <- 0
      X1 <- X_use; X1[, j] <- 1
      xb0 <- as.vector(X0 %*% b_t)
      xb1 <- as.vector(X1 %*% b_t)
      xr0 <- if (nr_t > 0) X0[, rand_idx, drop = FALSE] else NULL
      xr1 <- if (nr_t > 0) X1[, rand_idx, drop = FALSE] else NULL
      mu0 <- rowMeans(vapply(seq_len(R_t), function(r) {
        eta <- if (nr_t > 0) xb0 + as.vector(xr0 %*% dev_t[r, ]) else xb0
        exp(pmin(pmax(eta, -700), 700))
      }, numeric(xu)))
      mu1 <- rowMeans(vapply(seq_len(R_t), function(r) {
        eta <- if (nr_t > 0) xb1 + as.vector(xr1 %*% dev_t[r, ]) else xb1
        exp(pmin(pmax(eta, -700), 700))
      }, numeric(xu)))
      if (quantity == "me") mean(mu1 - mu0) else mean(mu1 / mu0 - 1)
    } else {
      if (!is.na(rj)) {
        dmu <- rowMeans(vapply(seq_len(R_t), function(r) {
          eta <- xb_t + as.vector(xr_t %*% dev_t[r, ])
          cmat_t[r, rj] * exp(pmin(pmax(eta, -700), 700))
        }, numeric(xu)))
      } else {
        dmu <- b_t[j] * mu_bar_t
      }
      if (quantity == "me") mean(dmu) else mean(X_use[, j] * dmu / mu_bar_t)
    }
  }, numeric(1))
}
