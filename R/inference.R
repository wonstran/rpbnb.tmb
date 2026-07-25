#' Natural-scale report values and one-parameter derivatives
#' @keywords internal
#' @noRd
.rpbnb_natural_report <- function(coef, family_code, lamLo, lamHi) {
  clamp_exp <- function(x) exp(pmin(pmax(x, -20), 20))
  clamp_deriv <- function(x, value) if (x > -20 && x < 20) value else 0
  add <- function(value, source, derivative, boundary = FALSE) {
    list(
      value = value, source = source, derivative = derivative,
      boundary = boundary
    )
  }
  # Each dependence link maps the unbounded working parameter onto a bounded
  # natural range.  Near either end the link's derivative collapses, so the
  # delta-method standard error tends to zero and reports false precision for
  # an estimate that is really pinned against the bound.  Flag those cases so
  # the caller sees NA rather than a confident-looking number.  Callers compare
  # on whichever scale the link is bounded on: the natural scale for the linear
  # links, the working scale for the exponential Clayton link, whose natural
  # bounds span nine orders of magnitude and make a relative margin useless.
  at_boundary <- function(value, lo, hi, margin = 0.02) {
    if (!is.finite(value) || !is.finite(lo) || !is.finite(hi)) return(FALSE)
    if (!(hi > lo)) return(TRUE)
    tol <- margin * (hi - lo)
    (value - lo) <= tol || (hi - value) <= tol
  }

  z_m1 <- unname(coef["log_m1"])
  z_m2 <- unname(coef["log_m2"])
  m1 <- clamp_exp(z_m1)
  m2 <- clamp_exp(z_m2)
  out <- list(
    m1 = add(m1, "log_m1", clamp_deriv(z_m1, m1)),
    m2 = add(m2, "log_m2", clamp_deriv(z_m2, m2))
  )

  if (family_code < 0L) return(out)
  z <- unname(coef["z_dep"])
  if (family_code == 0L) {
    eps <- 1e-6
    sig <- stats::plogis(z)
    span <- lamHi - lamLo
    value <- lamLo + span * (eps + (1 - 2 * eps) * sig)
    derivative <- span * (1 - 2 * eps) * sig * (1 - sig)
    out$lam <- add(
      value, "z_dep", derivative,
      boundary = at_boundary(value, lamLo, lamHi)
    )
  } else if (family_code == 1L) {
    theta <- FRANK_THETA_MAX * tanh(z / FRANK_THETA_MAX)
    derivative <- 1 - tanh(z / FRANK_THETA_MAX)^2
    out$theta <- add(
      theta, "z_dep", derivative,
      boundary = at_boundary(theta, -FRANK_THETA_MAX, FRANK_THETA_MAX)
    )
    tau <- .frank_tau(theta)
    tau_derivative <- if (derivative == 0) {
      0
    } else {
      as.numeric(numDeriv::grad(.frank_tau, theta)) * derivative
    }
    # Kendall's tau is a monotone reparameterisation of the dependence
    # parameter, so it is pinned exactly when that parameter is -- even though
    # tau itself is nowhere near +/-1 (Frank's ceiling is only tau = 0.891).
    out$tau <- add(
      tau, "z_dep", tau_derivative,
      boundary = out$theta$boundary || at_boundary(tau, -1, 1)
    )
  } else if (family_code == 2L) {
    rho <- tanh(z)
    drho <- 1 - rho^2
    out$rho <- add(rho, "z_dep", drho, boundary = at_boundary(rho, -1, 1))
    tau <- 2 / pi * asin(rho)
    out$tau <- add(
      tau, "z_dep", 2 / pi * drho / sqrt(1 - rho^2),
      boundary = out$rho$boundary || at_boundary(tau, -1, 1)
    )
  } else if (family_code == 3L) {
    theta <- clamp_exp(z)
    dtheta <- clamp_deriv(z, theta)
    out$theta <- add(
      theta, "z_dep", dtheta,
      boundary = at_boundary(z, -20, 20)
    )
    tau <- theta / (theta + 2)
    out$tau <- add(
      tau, "z_dep", 2 * dtheta / (theta + 2)^2,
      boundary = out$theta$boundary || at_boundary(tau, -1, 1)
    )
  }
  out
}

#' Warn when a dependence estimate is pinned against its link's bound
#'
#' The estimate is then an artefact of the link, not of the data, and its
#' delta-method standard error has already been set to NA.
#' @keywords internal
#' @noRd
.warn_boundary_report <- function(boundary_report, family_code) {
  flagged <- setdiff(boundary_report, c("m1", "m2"))
  if (!length(flagged)) return(invisible(FALSE))
  cause <- switch(
    as.character(family_code),
    "0" = "the frozen admissible Famoye bound",
    "1" = sprintf("the bounded Frank link (|theta| < %g)", FRANK_THETA_MAX),
    "2" = "the Gaussian link (|rho| < 1)",
    "3" = "the bounded Clayton link",
    "its link bound"
  )
  warning(
    "Dependence estimate(s) ", paste(flagged, collapse = ", "),
    " are pinned at a boundary of ", cause,
    ". The estimates are constrained by the link rather than identified by ",
    "the data, so their standard errors are reported as NA. Refit from ",
    "different starting values, or use a dependence family whose range ",
    "covers the association in these data.",
    call. = FALSE
  )
  invisible(TRUE)
}

#' Frank Kendall tau using the template's fixed quadrature rule
#' @keywords internal
#' @noRd
.frank_tau <- function(theta) {
  if (abs(theta) < 1e-4) {
    return(theta / 9 - theta^3 / 900)
  }
  x20 <- c(
    0.9931285991850949, 0.9639719272779138, 0.9122344282513259,
    0.8391169718222188, 0.7463319064601508, 0.6360536807265150,
    0.5108670019508271, 0.3737060887154196, 0.2277858511416451,
    0.07652652113349733
  )
  w20 <- c(
    0.01761400713915212, 0.04060142980038694, 0.06267204833410906,
    0.08327674157670475, 0.1019301198172404, 0.1181945319615184,
    0.1316886384491766, 0.1420961093183821, 0.1491729864726037,
    0.1527533871307259
  )
  f <- function(t) {
    near_zero <- abs(t) < 1e-6
    value <- numeric(length(t))
    value[near_zero] <- 1 - t[near_zero] / 2 + t[near_zero]^2 / 12
    value[!near_zero] <- t[!near_zero] / expm1(t[!near_zero])
    value
  }
  upper <- theta * 0.5 * (1 + x20)
  lower <- theta * 0.5 * (1 - x20)
  D1 <- 0.5 * sum(w20 * (f(upper) + f(lower)))
  1 - 4 / theta * (1 - D1)
}

#' Inference from the existing TMB objective
#' @keywords internal
#' @noRd
.rpbnb_inference <- function(obj, par, coef, par_names, free, mode,
                              family_code, lamLo, lamHi) {
  p <- length(par_names)
  free_names <- par_names[free]
  gradient <- obj$gr(par)
  names(gradient) <- free_names
  pdHess <- NA
  covariance <- NULL
  covariance_diag <- rep(NA_real_, p)
  names(covariance_diag) <- par_names

  if (mode != "none") {
    hessian <- stats::optimHess(par, obj$fn, obj$gr)
    hessian <- (hessian + t(hessian)) / 2
    factor <- try(chol(hessian), silent = TRUE)
    pdHess <- !inherits(factor, "try-error")

    if (mode == "full") {
      free_covariance <- if (pdHess) {
        chol2inv(factor)
      } else {
        solved <- try(solve(hessian), silent = TRUE)
        if (inherits(solved, "try-error")) hessian * NaN else solved
      }
      covariance <- matrix(
        NA_real_, p, p, dimnames = list(par_names, par_names)
      )
      covariance[free, free] <- free_covariance
      covariance_diag[free] <- diag(free_covariance)
    } else {
      inverse_diagonal <- rep(NA_real_, length(par))
      for (j in seq_along(par)) {
        unit <- rep(0, length(par))
        unit[j] <- 1
        column <- if (pdHess) {
          backsolve(factor, forwardsolve(t(factor), unit))
        } else {
          try(solve(hessian, unit), silent = TRUE)
        }
        if (!inherits(column, "try-error")) inverse_diagonal[j] <- column[j]
      }
      covariance_diag[free] <- inverse_diagonal
    }
  }

  se <- rep(NA_real_, length(covariance_diag))
  names(se) <- names(covariance_diag)
  nonnegative <- is.finite(covariance_diag) & covariance_diag >= 0
  se[nonnegative] <- sqrt(covariance_diag[nonnegative])
  report_items <- .rpbnb_natural_report(
    coef, family_code = family_code, lamLo = lamLo, lamHi = lamHi
  )
  report_value <- vapply(report_items, `[[`, numeric(1), "value")
  report_sd <- vapply(report_items, function(item) {
    index <- match(item$source, par_names)
    abs(item$derivative) * se[index]
  }, numeric(1))
  boundary_report <- vapply(
    report_items, `[[`, logical(1), "boundary"
  )
  report_sd[boundary_report] <- NA_real_

  report_covariance <- matrix(
    NA_real_, length(report_items), length(report_items),
    dimnames = list(names(report_items), names(report_items))
  )
  if (mode == "full") {
    for (i in seq_along(report_items)) {
      for (j in seq_along(report_items)) {
        ii <- match(report_items[[i]]$source, par_names)
        jj <- match(report_items[[j]]$source, par_names)
        report_covariance[i, j] <-
          report_items[[i]]$derivative *
          covariance[ii, jj] *
          report_items[[j]]$derivative
      }
    }
  } else if (mode == "diag") {
    diag(report_covariance) <- report_sd^2
  }

  compact_sdreport <- structure(
    list(
      value = report_value,
      sd = report_sd,
      cov = report_covariance,
      par.fixed = coef,
      cov.fixed = covariance,
      pdHess = pdHess,
      gradient.fixed = gradient,
      fixed.value = coef,
      fixed.sd = se
    ),
    class = "rpbnb_sdreport"
  )

  list(
    vcov = covariance,
    vcov_diag = if (mode == "full") NULL else covariance_diag,
    se = se,
    sdreport = compact_sdreport,
    report = report_value,
    boundary_report = names(report_items)[boundary_report]
  )
}

#' @export
summary.rpbnb_sdreport <- function(object,
                                    select = c("all", "fixed", "report"),
                                    ...) {
  select <- match.arg(select)
  fixed <- cbind(
    Estimate = object$fixed.value,
    `Std. Error` = object$fixed.sd
  )
  report <- cbind(
    Estimate = object$value,
    `Std. Error` = object$sd
  )
  if (select == "fixed") return(fixed)
  if (select == "report") return(report)
  rbind(fixed, report)
}
