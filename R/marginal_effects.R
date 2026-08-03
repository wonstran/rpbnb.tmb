#' Centre/scale lookup for a design that was standardized before fitting
#'
#' Returns centre and scale vectors aligned to \code{cn}, defaulting to 0 and 1
#' for any column not named in \code{scaling} (binary indicators, the
#' intercept, anything left in raw units).
#' @keywords internal
#' @noRd
.scaling_vec <- function(scaling, cn) {
  ctr <- setNames(numeric(length(cn)), cn)
  scl <- setNames(rep(1, length(cn)), cn)
  if (is.null(scaling)) return(list(center = ctr, scale = scl))
  if (!is.list(scaling)) {
    stop("`scaling` must be a named list of c(center =, scale =) pairs.",
         call. = FALSE)
  }
  for (v in names(scaling)) {
    if (!v %in% cn) next
    p <- scaling[[v]]
    if (is.null(names(p)) || !all(c("center", "scale") %in% names(p))) {
      stop("scaling[['", v, "']] needs named `center` and `scale` entries.",
           call. = FALSE)
    }
    if (!is.finite(p[["scale"]]) || p[["scale"]] <= 0) {
      stop("scaling[['", v, "']]$scale must be a positive finite number.",
           call. = FALSE)
    }
    ctr[v] <- p[["center"]]
    scl[v] <- p[["scale"]]
  }
  list(center = ctr, scale = scl)
}

#' Flag which selected columns are already on a log scale
#'
#' Validates \code{log_vars} against the design and returns a logical vector
#' aligned to \code{sel}.
#' @keywords internal
#' @noRd
.log_vars_flag <- function(log_vars, cn, sel, is_bin) {
  flag <- rep(FALSE, length(sel))
  if (is.null(log_vars)) return(flag)
  if (!is.character(log_vars)) {
    stop("`log_vars` must be a character vector of column names.", call. = FALSE)
  }
  unknown <- setdiff(log_vars, cn)
  if (length(unknown)) {
    stop("log_vars not found in the design: ",
         paste(unknown, collapse = ", "), call. = FALSE)
  }
  flag <- cn[sel] %in% log_vars
  # A 0/1 column is not a log of anything, and treating it as one would divide
  # by exp(0) = 1 for the zeros and silently report nonsense.
  if (any(flag & is_bin)) {
    stop("log_vars names binary column(s): ",
         paste(cn[sel][flag & is_bin], collapse = ", "),
         ". A 0/1 indicator is not a log-transformed variable.", call. = FALSE)
  }
  flag
}

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
#' @param scaling Optional named list of \code{c(center =, scale =)} pairs, one
#'   per covariate that was centred and/or scaled before fitting, used to report
#'   the results in the covariate's original units.  Columns not named are left
#'   alone, so binary indicators and untransformed variables need no entry.
#'
#'   Only the reported quantity is rescaled; the fitted design is not rebuilt.
#'   That distinction is not cosmetic when a standardized covariate also carries
#'   a random coefficient.  The random term enters as
#'   \code{x_std * dev}, so substituting \code{x_raw = x_std * s + c} would add a
#'   \code{(c/s) * dev} random intercept the model never estimated -- which is
#'   how a centred random slope turns back into the disguised random intercept
#'   that centring was meant to remove.  The chain rule avoids that entirely:
#'   \eqn{\partial\mu/\partial x_{raw} = (\partial\mu/\partial x_{std})/s}, and
#'   the elasticity's leading factor becomes
#'   \eqn{x_{raw}/s = x_{std} + c/s}.  Binary effects are untouched.
#' @param log_vars Optional character vector naming covariates that are ALREADY
#'   a log, so that results are reported per unit of the underlying variable
#'   \eqn{v_j = \exp(x_j)} rather than per unit of \eqn{x_j} itself.
#'
#'   Without this the elasticity formula treats \eqn{\log v} as an ordinary
#'   regressor and returns \eqn{\bar{x} b}, the elasticity with respect to the
#'   LOG -- for the truck model's \code{LNAADT_3} that is 8.59, where the
#'   elasticity with respect to AADT is the coefficient itself, 0.871.  A
#'   ten-fold overstatement that looks perfectly plausible in a table.
#'
#'   Because \eqn{\partial\log\mu/\partial\log v = (\partial\mu/\partial x)/\mu},
#'   the elasticity is \eqn{E[d\mu/\mu]/s} and the marginal effect is
#'   \eqn{E[d\mu/(s v)]}, with \eqn{v = \exp(x_{std} s + c)} picking up any
#'   \code{scaling} applied to the logged column.  Named columns are reported
#'   with type \code{"log-continuous"}.  Binary columns are rejected.
#' @param ... Not used.
#' @return A data frame (single margin) or named list of two data frames
#'   (\code{"both"}).
#' @export
rpbnb_tmb_marginal_effects <- function(fit,
                                       which = c("y1", "y2", "both"),
                                       type  = c("AME", "MEM"),
                                       vars  = NULL,
                                       include_intercept = FALSE,
                                       digits = 4L,
                                       scaling = NULL,
                                       log_vars = NULL, ...) {
  which <- match.arg(which)
  type  <- match.arg(type)
  if (!identical(fit$inference, "full") && !is.null(fit$rp_meta)) {
    warning(
      'Marginal-effect standard errors require the full covariance; ',
      'refit with inference = "full". Point estimates are still returned.',
      call. = FALSE
    )
  }
  res <- if (which == "both") {
    list(y1 = .me_one_eq(fit, 1L, type, vars, include_intercept, digits,
                         scaling, log_vars),
         y2 = .me_one_eq(fit, 2L, type, vars, include_intercept, digits,
                         scaling, log_vars))
  } else {
    .me_one_eq(fit, if (which == "y1") 1L else 2L,
               type, vars, include_intercept, digits, scaling, log_vars)
  }
  invisible(res)
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
                                   digits = 4L,
                                   scaling = NULL,
                                   log_vars = NULL, ...) {
  which <- match.arg(which)
  type  <- match.arg(type)
  if (!identical(fit$inference, "full") && !is.null(fit$rp_meta)) {
    warning(
      'Elasticity standard errors require the full covariance; ',
      'refit with inference = "full". Point estimates are still returned.',
      call. = FALSE
    )
  }
  res <- if (which == "both") {
    list(y1 = .el_one_eq(fit, 1L, type, vars, include_intercept, digits,
                         scaling, log_vars),
         y2 = .el_one_eq(fit, 2L, type, vars, include_intercept, digits,
                         scaling, log_vars))
  } else {
    .el_one_eq(fit, if (which == "y1") 1L else 2L,
               type, vars, include_intercept, digits, scaling, log_vars)
  }
  invisible(res)
}

# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

#' Mean of exponentiated predictors accumulated over random draws
#' @keywords internal
#' @noRd
.draw_mean_exp <- function(xb, xr, dev) {
  total <- numeric(length(xb))
  has_random <- !is.null(xr) && ncol(xr) > 0L
  for (r in seq_len(nrow(dev))) {
    eta <- xb + if (has_random) as.vector(xr %*% dev[r, ]) else 0
    total <- total + exp(pmin(pmax(eta, -700), 700))
  }
  total / nrow(dev)
}

#' Weighted mean of exponentiated predictors over random draws
#' @keywords internal
#' @noRd
.draw_mean_weighted_exp <- function(xb, xr, dev, weights) {
  total <- numeric(length(xb))
  for (r in seq_len(nrow(dev))) {
    eta <- xb + as.vector(xr %*% dev[r, ])
    total <- total +
      weights[r] * exp(pmin(pmax(eta, -700), 700))
  }
  total / nrow(dev)
}

#' Compute marginal effects for one equation
#' @keywords internal
#' @noRd
.me_one_eq <- function(fit, eq, type, vars, include_intercept, digits,
                       scaling = NULL, log_vars = NULL) {
  .diag_one_eq(fit, eq, type, vars, include_intercept, digits,
               quantity = "me", scaling = scaling, log_vars = log_vars)
}

#' Compute elasticities for one equation
#' @keywords internal
#' @noRd
.el_one_eq <- function(fit, eq, type, vars, include_intercept, digits,
                       scaling = NULL, log_vars = NULL) {
  .diag_one_eq(fit, eq, type, vars, include_intercept, digits,
               quantity = "elas", scaling = scaling, log_vars = log_vars)
}

#' Core diagnostic engine for one equation
#' @keywords internal
#' @noRd
.diag_one_eq <- function(fit, eq, type, vars, include_intercept, digits,
                         quantity = c("me", "elas"), scaling = NULL,
                         log_vars = NULL) {
  quantity <- match.arg(quantity)

  # ---- Extract design & metadata ----
  X     <- if (eq == 1L) fit$X1 else fit$X2
  if (is.null(X) || is.null(fit$rp_meta)) {
    stop(
      'Marginal effects require a fit created with keep = "postfit" or "full".',
      call. = FALSE
    )
  }
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

  # Centre/scale of the pre-fit standardization, or 0/1 throughout when the
  # design was never transformed.  X_use itself is deliberately left on the
  # scale the model was fitted on; see the `scaling` documentation.
  sv  <- .scaling_vec(scaling, cn)
  ctr <- sv$center
  scl <- sv$scale
  is_log <- .log_vars_flag(log_vars, cn, sel, is_bin)

  # ---- Random coefficient draws ----
  if (nr > 0) {
    Z     <- fit$rp_meta[[paste0("Z", eq)]]
    dist  <- fit$rp_meta[[paste0("dist", eq)]]
    sign  <- fit$rp_meta[[paste0("sign", eq)]]
    scale_lab <- vapply(dist, function(d) {
      rand_dist_registry[[d]]$scale_label
    }, character(1))
    scale_names <- paste0(scale_lab, eq, ":", cn[rand_idx])
    scales <- exp(fit$coef[scale_names])
    rr    <- rand_realize(Z, dist, sign, b[rand_idx], scales)
    dev   <- rr$dev
    cmat  <- rr$coef        # n_draws x n_rand
  } else {
    dev   <- matrix(0, 1, 0)
    cmat  <- NULL
    scale_names <- character(0)
  }

  # ---- Compute point estimates ----
  xb   <- as.vector(X_use %*% b)
  xr   <- if (nr > 0) X_use[, rand_idx, drop = FALSE] else NULL

  mu_bar <- .draw_mean_exp(xb, xr, dev)

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
      mu0 <- .draw_mean_exp(xb0, xr0, dev)
      mu1 <- .draw_mean_exp(xb1, xr1, dev)
      if (quantity == "me") mean(mu1 - mu0) else mean(mu1 / mu0 - 1)
    } else {
      if (!is.na(rj)) {
        # Random continuous: average of coef_rj * mu_r over draws
        dmu <- .draw_mean_weighted_exp(xb, xr, dev, cmat[, rj])
      } else {
        # Fixed continuous: b_j * mu
        dmu <- b[j] * mu_bar
      }
      if (is_log[m]) {
        # x_j is already log(v_j): report per unit of v_j, not of log v_j.
        # d/dlog v = (1/s) d/dx_std, so the elasticity is E[dmu/mu]/s and the
        # marginal effect divides by v = exp(x_std * s + c) as well.
        if (quantity == "me") {
          mean(dmu / (scl[j] * exp(X_use[, j] * scl[j] + ctr[j])))
        } else {
          mean(dmu / mu_bar) / scl[j]
        }
      } else if (quantity == "me") {
        mean(dmu) / scl[j]
      } else {
        mean((X_use[, j] + ctr[j] / scl[j]) * dmu / mu_bar)
      }
    }
  }, numeric(1))

  # ---- Delta-method standard errors ----
  theta_names <- c(paste0(coef_prefix, cn), scale_names)
  theta_hat   <- fit$coef[theta_names]
  # Extract vcov for these parameters by position (the full vcov matrix has
  # dimnames matching par_names, but mapped/fixed parameters may have NAs).
  idx <- match(theta_names, names(fit$coef))
  V <- vcov(fit)[idx, idx, drop = FALSE]

  se <- rep(NA_real_, length(res))
  if (all(is.finite(V)) && nrow(V) == length(theta_hat)) {
    G <- numDeriv::jacobian(.estimand_t, theta_hat,
                            fit = fit, eq = eq, sel = sel, cn = cn,
                            is_bin = is_bin, quantity = quantity,
                            X_use = X_use, rand_idx = rand_idx,
                            ctr = ctr, scl = scl, is_log = is_log)
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
    Type = ifelse(
      is_bin,
      if (quantity == "me") "binary" else "semi-elas",
      ifelse(is_log, "log-continuous", "continuous")
    ),
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
                        X_use, rand_idx,
                        ctr = setNames(numeric(length(cn)), cn),
                        scl = setNames(rep(1, length(cn)), cn),
                        is_log = rep(FALSE, length(sel))) {
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
  } else {
    dev_t  <- matrix(0, 1, 0)
    cmat_t <- NULL
  }

  xb_t <- as.vector(X_use %*% b_t)
  xr_t <- if (nr_t > 0) X_use[, rand_idx, drop = FALSE] else NULL

  mu_bar_t <- .draw_mean_exp(xb_t, xr_t, dev_t)

  vapply(seq_along(sel), function(m) {
    j <- sel[m]; rj <- match(j, rand_idx)
    if (is_bin[m]) {
      X0 <- X_use; X0[, j] <- 0
      X1 <- X_use; X1[, j] <- 1
      xb0 <- as.vector(X0 %*% b_t)
      xb1 <- as.vector(X1 %*% b_t)
      xr0 <- if (nr_t > 0) X0[, rand_idx, drop = FALSE] else NULL
      xr1 <- if (nr_t > 0) X1[, rand_idx, drop = FALSE] else NULL
      mu0 <- .draw_mean_exp(xb0, xr0, dev_t)
      mu1 <- .draw_mean_exp(xb1, xr1, dev_t)
      if (quantity == "me") mean(mu1 - mu0) else mean(mu1 / mu0 - 1)
    } else {
      if (!is.na(rj)) {
        dmu <- .draw_mean_weighted_exp(
          xb_t, xr_t, dev_t, cmat_t[, rj]
        )
      } else {
        dmu <- b_t[j] * mu_bar_t
      }
      if (is_log[m]) {
        if (quantity == "me") {
          mean(dmu / (scl[j] * exp(X_use[, j] * scl[j] + ctr[j])))
        } else {
          mean(dmu / mu_bar_t) / scl[j]
        }
      } else if (quantity == "me") {
        mean(dmu) / scl[j]
      } else {
        mean((X_use[, j] + ctr[j] / scl[j]) * dmu / mu_bar_t)
      }
    }
  }, numeric(1))
}
