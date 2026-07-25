#' Monotone links from the working dependence parameter to natural scale
#'
#' Every family's link is monotone increasing in \code{z}, which is what makes
#' it legitimate to map the endpoints of a working-scale interval through it.
#' \code{.rpbnb_natural_report()} and \code{rpbnb_tmb_dependence_profile()} both
#' read their values from here so the two cannot disagree.
#'
#' \code{map} takes a scalar: \code{.frank_tau()} is not vectorized.
#' @keywords internal
#' @noRd
.rpbnb_dependence_link <- function(family_code,
                                   lamLo = NA_real_, lamHi = NA_real_) {
  if (family_code < 0L) return(NULL)
  if (family_code == 0L) {
    eps <- 1e-6
    span <- lamHi - lamLo
    return(list(
      names = "lam",
      map = function(z) {
        stopifnot(length(z) == 1L)
        c(lam = lamLo + span * (eps + (1 - 2 * eps) * stats::plogis(z)))
      }
    ))
  }
  if (family_code == 1L) {
    return(list(
      names = c("theta", "tau"),
      map = function(z) {
        stopifnot(length(z) == 1L)
        theta <- FRANK_THETA_MAX * tanh(z / FRANK_THETA_MAX)
        c(theta = theta, tau = .frank_tau(theta))
      }
    ))
  }
  if (family_code == 2L) {
    return(list(
      names = c("rho", "tau"),
      map = function(z) {
        stopifnot(length(z) == 1L)
        rho <- tanh(z)
        c(rho = rho, tau = 2 / pi * asin(rho))
      }
    ))
  }
  list(
    names = c("theta", "tau"),
    map = function(z) {
      stopifnot(length(z) == 1L)
      theta <- exp(pmin(pmax(z, -20), 20))
      c(theta = theta, tau = theta / (theta + 2))
    }
  )
}

#' Profile-likelihood interval for the working dependence parameter
#'
#' Returns working-scale endpoints, or \code{NULL} when no profile can be
#' attempted so the caller can fall back to a Wald interval.
#' @keywords internal
#' @noRd
.dependence_profile_ci <- function(fit, level, ...) {
  obj <- fit$obj
  if (is.null(obj)) return(NULL)

  # Presence is not liveness: a fit round-tripped through saveRDS() keeps the
  # obj field but its external pointer is dead.
  probe <- try(obj$fn(fit$optimizer$par), silent = TRUE)
  if (inherits(probe, "try-error") ||
      length(probe) != 1L || !is.finite(probe)) {
    return(NULL)
  }

  # `obj` is fit$obj's environment, not a copy of it, so the last.par.best
  # write below is visible to the caller's fit too. Capture the pre-call
  # value now (before it is ever touched) and restore it on every exit path,
  # mirroring how TMB::tmbprofile() protects its own obj$env writes via
  # on.exit(restore.oldvars()). Nothing re-reads fit$obj after this call
  # today, but leaving the mutation in place would be an unenforced
  # invariant, not a guarantee.
  old_last_par_best <- obj$env$last.par.best
  on.exit(obj$env$last.par.best <- old_last_par_best, add = TRUE)

  # tmbprofile() profiles from last.par.best, not obj$par. Restore it rather
  # than trusting whatever last touched the objective.
  obj$env$last.par.best <- fit$optimizer$par

  dots <- list(...)

  # `trace` is a real, documented tmbprofile() argument and this function's
  # own docs promise `...` is forwarded to it. Building `args` as
  # c(list(obj, "z_dep", trace = FALSE), dots, extra) let a caller-supplied
  # `trace` collide with the hardcoded one: do.call() then threw "formal
  # argument matched by multiple actual arguments", which the try() below
  # swallowed into a silent NULL -- downgrading a legitimate request to a
  # Wald interval with no diagnostic. Resolve `trace` once and drop it from
  # `dots` so no duplicate can ever reach do.call().
  if ("trace" %in% names(dots)) {
    trace <- dots[["trace"]]
    dots[["trace"]] <- NULL
  } else {
    trace <- FALSE
  }

  # `obj` and `"z_dep"` are likewise hardcoded (positionally, into
  # tmbprofile()'s `obj` and `name` formals) and not meant to be overridden --
  # this function's entire contract is profiling *this* objective at
  # *z_dep*. A caller-supplied `name` or `obj` in `...` would hit the same
  # duplicate-argument crash as `trace`; drop them so our fixed values always
  # win instead of erroring.
  dots[["obj"]] <- NULL
  dots[["name"]] <- NULL

  run <- function(extra) {
    args <- c(list(obj, "z_dep", trace = trace), dots, extra)
    result <- try(do.call(TMB::tmbprofile, args), silent = TRUE)
    if (inherits(result, "try-error")) NULL else result
  }

  profile <- run(NULL)
  if (is.null(profile)) return(NULL)
  endpoints <- as.numeric(stats::confint(profile, level = level))

  # confint.tmbprofile() locates endpoints with approx(), so an uncrossed
  # profile yields NA rather than an infinite bound. Widening ytol searches
  # further before giving up.
  if (anyNA(endpoints) && !("ytol" %in% names(dots))) {
    wider <- run(list(ytol = 10))
    if (!is.null(wider)) {
      wider_endpoints <- as.numeric(stats::confint(wider, level = level))
      if (sum(is.na(wider_endpoints)) < sum(is.na(endpoints))) {
        profile <- wider
        endpoints <- wider_endpoints
      }
    }
  }

  structure(endpoints, profile = profile)
}

#' Natural-scale report values and one-parameter derivatives
#' @keywords internal
#' @noRd
.rpbnb_natural_report <- function(coef, family_code, lamLo, lamHi) {
  clamp_exp <- function(x) exp(pmin(pmax(x, -20), 20))
  # Matches the template's CondExp clamp, which passes the input through at
  # exactly +/-20 rather than treating equality as clamped.
  clamp_deriv <- function(x, value) if (x >= -20 && x <= 20) value else 0
  add <- function(value, source, derivative, boundary = FALSE,
                  side = NA_character_) {
    list(
      value = value, source = source, derivative = derivative,
      boundary = isTRUE(boundary),
      side = if (isTRUE(boundary)) side else NA_character_
    )
  }

  # A reported delta-method standard error is untrustworthy for three distinct
  # reasons, and only these three.  Conflating any of them with "the estimate
  # is close to the edge of its range" would suppress good inference: a
  # Gaussian correlation of 0.964 sits near |rho| = 1, but |rho| < 1 is the
  # model's genuine parameter domain and drho/dz is still 0.07.
  #
  # (a) SMOOTH artificial cap.  The Famoye bounds frozen at the starting values
  #     and the Frank overflow guard squash the working parameter through a
  #     logistic or tanh, so the map is already compressed well before the cap
  #     is reached and the estimate is demonstrably being shaped by it.
  #     Proximity is the right test here.  The 2% margin is a judgement call:
  #     close enough that the cap dominates, loose enough not to fire on
  #     ordinary fits.
  near_smooth_cap <- function(value, lo, hi, margin = 0.02) {
    if (!is.finite(value) || !is.finite(lo) || !is.finite(hi)) {
      return(NA_character_)
    }
    if (!(hi > lo)) return("degenerate")
    tol <- margin * (hi - lo)
    if ((value - lo) <= tol) return("lower")
    if ((hi - value) <= tol) return("upper")
    NA_character_
  }
  # (b) HARD clamp reached.  pmin/pmax clamps are the *identity* strictly
  #     inside their range, so a log_m of 19.3 is untouched by a clamp at 20
  #     and flagging it would be a false positive.
  #
  #     The comparison is non-strict, because exact equality is a kink, not an
  #     interior point: at z = -20 the left derivative is 0 and the right is
  #     exp(-20), so the delta method has no unique derivative there.  The
  #     template's CondExp happens to select the pass-through branch, and
  #     clamp_deriv() matches it so the R/C++ contract stays exact -- but which
  #     branch CppAD differentiates is a mechanical detail, not a licence to
  #     report the result as identified inference.
  clamp_reached <- function(z, lo = -20, hi = 20) {
    if (!is.finite(z)) return("degenerate")
    if (z <= lo) return("lower")
    if (z >= hi) return("upper")
    NA_character_
  }
  # (c) The link's derivative has collapsed, so the delta-method standard error
  #     is exactly zero or NaN whatever the working-scale uncertainty is.  This
  #     is what catches a genuinely bounded domain: tanh reaches 1 in double
  #     precision near z = 19.1, and past that rho carries no information.
  flag <- function(derivative, cap_side = NA_character_) {
    if (!is.na(cap_side)) return(cap_side)
    if (!is.finite(derivative) || derivative == 0) return("degenerate")
    NA_character_
  }

  z_m1 <- unname(coef["log_m1"])
  z_m2 <- unname(coef["log_m2"])
  m1 <- clamp_exp(z_m1)
  m2 <- clamp_exp(z_m2)
  d_m1 <- clamp_deriv(z_m1, m1)
  d_m2 <- clamp_deriv(z_m2, m2)
  side_m1 <- flag(d_m1, clamp_reached(z_m1))
  side_m2 <- flag(d_m2, clamp_reached(z_m2))
  out <- list(
    m1 = add(m1, "log_m1", d_m1, !is.na(side_m1), side_m1),
    m2 = add(m2, "log_m2", d_m2, !is.na(side_m2), side_m2)
  )

  if (family_code < 0L) return(out)
  z <- unname(coef["z_dep"])
  link_values <- .rpbnb_dependence_link(
    family_code, lamLo = lamLo, lamHi = lamHi
  )$map(z)
  if (family_code == 0L) {
    eps <- 1e-6
    sig <- stats::plogis(z)
    span <- lamHi - lamLo
    value <- link_values[["lam"]]
    derivative <- span * (1 - 2 * eps) * sig * (1 - sig)
    # Frozen at the starting values, so the interval itself is the artefact,
    # and the logistic compresses lam well before it reaches either end.
    side <- flag(derivative, near_smooth_cap(value, lamLo, lamHi))
    out$lam <- add(value, "z_dep", derivative, !is.na(side), side)
  } else if (family_code == 1L) {
    # |theta| < FRANK_THETA_MAX is an overflow guard, not a property of the
    # Frank family, whose theta is unbounded.
    theta <- link_values[["theta"]]
    derivative <- 1 - tanh(z / FRANK_THETA_MAX)^2
    side <- flag(
      derivative,
      near_smooth_cap(theta, -FRANK_THETA_MAX, FRANK_THETA_MAX)
    )
    out$theta <- add(theta, "z_dep", derivative, !is.na(side), side)
    tau <- link_values[["tau"]]
    tau_derivative <- if (derivative == 0) {
      0
    } else {
      as.numeric(numDeriv::grad(.frank_tau, theta)) * derivative
    }
    # Kendall's tau is a monotone reparameterisation of the dependence
    # parameter, so it is pinned exactly when that parameter is -- even though
    # tau itself is nowhere near +/-1 (Frank's ceiling is only tau = 0.891).
    tau_side <- if (out$theta$boundary) out$theta$side else flag(tau_derivative)
    out$tau <- add(
      tau, "z_dep", tau_derivative, !is.na(tau_side), tau_side
    )
  } else if (family_code == 2L) {
    # |rho| < 1 is the model's own domain, so there is no artificial cap here:
    # only genuine saturation of tanh disqualifies the standard error.
    rho <- link_values[["rho"]]
    drho <- 1 - rho^2
    rho_side <- flag(drho)
    out$rho <- add(rho, "z_dep", drho, !is.na(rho_side), rho_side)
    tau <- link_values[["tau"]]
    # 2/pi * drho / sqrt(1 - rho^2) is 0/0 once tanh saturates; the limit of
    # dtau/dz is 0 there, and leaving NaN would poison the report covariance.
    tau_derivative <- if (drho <= 0) 0 else 2 / pi * drho / sqrt(1 - rho^2)
    tau_side <- if (out$rho$boundary) out$rho$side else flag(tau_derivative)
    out$tau <- add(
      tau, "z_dep", tau_derivative, !is.na(tau_side), tau_side
    )
  } else if (family_code == 3L) {
    # exp() clamped at +/-20 in both R and the template.  The clamp is the
    # identity inside (-20, 20), so proximity to it means nothing.
    theta <- link_values[["theta"]]
    dtheta <- clamp_deriv(z, theta)
    theta_side <- flag(dtheta, clamp_reached(z))
    out$theta <- add(
      theta, "z_dep", dtheta, !is.na(theta_side), theta_side
    )
    tau <- link_values[["tau"]]
    tau_derivative <- 2 * dtheta / (theta + 2)^2
    tau_side <- if (out$theta$boundary) out$theta$side else flag(tau_derivative)
    out$tau <- add(
      tau, "z_dep", tau_derivative, !is.na(tau_side), tau_side
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
.warn_boundary_report <- function(boundary_report, family_code,
                                   sides = NULL) {
  dispersion <- intersect(boundary_report, c("m1", "m2"))
  dependence <- setdiff(boundary_report, c("m1", "m2"))
  if (!length(dispersion) && !length(dependence)) return(invisible(FALSE))
  side_of <- function(item) {
    if (is.null(sides) || !item %in% names(sides)) return(NA_character_)
    unname(sides[[item]])
  }

  # Dispersion and dependence fail for unrelated reasons, so they get separate
  # diagnoses.  Describing a clamped m1 as "pinned against the Frank link" and
  # advising a different dependence family would both be false.  The two sides
  # of the dispersion clamp also need opposite advice: at the lower clamp m -> 0
  # and the margin is Poisson, at the upper clamp m is enormous and the margin
  # is degenerately over-dispersed.  Recommending a Poisson margin there would
  # be the opposite model.
  for (side in unique(vapply(dispersion, side_of, character(1)))) {
    items <- dispersion[vapply(dispersion, side_of, character(1)) %in% side]
    remedy <- if (identical(side, "lower")) {
      paste0(
        "The margin is effectively Poisson; use poisson_1/poisson_2 to fit ",
        "that exact limit instead."
      )
    } else if (identical(side, "upper")) {
      paste0(
        "The margin is degenerately over-dispersed, which usually means the ",
        "mean model is badly misspecified or the counts are near-degenerate; ",
        "a Poisson margin is not the remedy."
      )
    } else {
      paste0(
        "The delta-method derivative has collapsed, so the dispersion is not ",
        "identified at this optimum."
      )
    }
    warning(
      "Dispersion estimate(s) ", paste(items, collapse = ", "),
      " are pinned at the log-dispersion clamp (|log m| <= 20), so the ",
      "estimate is set by the clamp rather than by the data and its standard ",
      "error is reported as NA. ", remedy,
      call. = FALSE
    )
  }
  if (length(dependence)) {
    cause <- switch(
      as.character(family_code),
      "0" = "the Famoye bounds, which are frozen at the starting values",
      "1" = sprintf(
        "the Frank overflow guard (|theta| < %g)", FRANK_THETA_MAX
      ),
      "2" = "numerical saturation of the Gaussian link (|rho| = 1)",
      "3" = "the Clayton exp() clamp",
      "an implementation bound"
    )
    warning(
      "Dependence estimate(s) ", paste(dependence, collapse = ", "),
      " are pinned at a boundary of ", cause,
      ". The estimates are constrained by the implementation rather than ",
      "identified by the data, so their standard errors are reported as NA. ",
      "Refit from different starting values, or use a dependence family whose ",
      "range covers the association in these data. ",
      "rpbnb_tmb_dependence_profile() gives an interval that stays valid here.",
      call. = FALSE
    )
  }
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
  # If a standard error has been withdrawn as unidentified, the corresponding
  # covariances are unidentified too.  Leaving them numeric would let callers
  # reconstruct exactly the quantity we just declared meaningless -- and on a
  # non-positive-definite Hessian those entries can even be negative.
  if (any(boundary_report)) {
    report_covariance[boundary_report, ] <- NA_real_
    report_covariance[, boundary_report] <- NA_real_
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
    boundary_report = names(report_items)[boundary_report],
    boundary_sides = vapply(
      report_items[boundary_report], `[[`, character(1), "side"
    )
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

#' Confidence interval for a fitted dependence parameter
#'
#' Returns a profile-likelihood interval for the dependence parameter of a
#' fitted model, computed on the unconstrained working scale and mapped through
#' the family's monotone link. This is the tool to reach for when
#' \code{summary()} reports \code{NA} for the dependence standard error: at a
#' boundary the delta-method derivative collapses, so
#' \eqn{SE = |d\theta/dz| \cdot SE(z)} is a \eqn{0 \times \infty} product and no
#' symmetric standard error exists. A profile interval needs no derivative.
#'
#' For Famoye this is the usual situation, because the admissible lambda
#' interval is frozen at the starting values (see \code{lambda_bounds} in
#' \code{\link{fit_rpbnb_tmb}}). An interval around a pinned estimate is still
#' an interval around an artefact -- widen the box first by refitting from
#' better starting values.
#'
#' @param fit An object of class \code{rpbnb_tmb_fit}.
#' @param level Coverage level; one number strictly between 0 and 1.
#' @param method \code{"profile"} (default) for a profile-likelihood interval,
#'   or \code{"wald"} for a Wald interval on the working scale mapped through
#'   the same link. \code{"profile"} requires the TMB objective, retained only
#'   under \code{keep = "full"}; when it is unavailable this degrades to
#'   \code{"wald"} with a warning rather than failing.
#' @param ... Passed to \code{\link[TMB]{tmbprofile}}, e.g. \code{ytol} or
#'   \code{parm.range}.
#' @return A data frame with one row per reported dependence quantity and
#'   columns \code{parameter}, \code{estimate}, \code{lower}, \code{upper},
#'   \code{level}, and \code{method} -- the method actually used, which may
#'   differ from the one requested. The raw profile, when one was computed, is
#'   attached as \code{attr(, "profile")} so it can be plotted.
#' @seealso \code{\link{fit_rpbnb_tmb}}
#' @export
#' @examples
#' \dontrun{
#' fit <- fit_rpbnb_tmb(y1 ~ x, y2 ~ x, data = d,
#'                      dependence = "famoye", keep = "full")
#' rpbnb_tmb_dependence_profile(fit)
#' plot(attr(rpbnb_tmb_dependence_profile(fit), "profile"))
#' }
rpbnb_tmb_dependence_profile <- function(fit, level = 0.95,
                                        method = c("profile", "wald"), ...) {
  if (!inherits(fit, "rpbnb_tmb_fit")) {
    stop("`fit` must be an object of class \"rpbnb_tmb_fit\".", call. = FALSE)
  }
  if (length(level) != 1L || !is.numeric(level) || is.na(level) ||
      level <= 0 || level >= 1) {
    stop("`level` must be one number strictly between 0 and 1.", call. = FALSE)
  }
  method <- match.arg(method)

  family_code <- .resolve_family_code(fit$dependence)
  if (family_code < 0L) {
    stop("An independence fit estimates no dependence parameter, so there is ",
         "nothing to profile.", call. = FALSE)
  }

  bounds <- fit$lambda_bounds
  link <- .rpbnb_dependence_link(
    family_code,
    lamLo = if (is.null(bounds)) NA_real_ else bounds[["lower"]],
    lamHi = if (is.null(bounds)) NA_real_ else bounds[["upper"]]
  )
  z_hat <- unname(fit$coef[["z_dep"]])

  profile <- NULL
  endpoints <- NULL
  if (identical(method, "profile")) {
    endpoints <- .dependence_profile_ci(fit, level = level, ...)
    if (is.null(endpoints)) {
      warning(
        "The TMB objective is unavailable, so no profile can be computed; ",
        "falling back to a Wald interval on the working scale. Refit with ",
        "keep = \"full\" for a profile interval.",
        call. = FALSE
      )
      method <- "wald"
    } else {
      profile <- attr(endpoints, "profile")
      endpoints <- as.numeric(endpoints)
    }
  }
  if (identical(method, "wald")) {
    se_z <- unname(fit$se[["z_dep"]])
    if (!is.finite(se_z)) {
      warning(
        "No standard error is available for z_dep, so the interval endpoints ",
        "are NA. Refit with inference = \"full\" or \"diag\".",
        call. = FALSE
      )
      endpoints <- c(NA_real_, NA_real_)
    } else {
      half <- stats::qnorm(1 - (1 - level) / 2) * se_z
      endpoints <- c(z_hat - half, z_hat + half)
    }
  }

  # Every link is monotone increasing, so the lower working endpoint maps to the
  # lower natural endpoint.
  map_or_na <- function(z) {
    if (!is.finite(z)) {
      stats::setNames(rep(NA_real_, length(link$names)), link$names)
    } else {
      link$map(z)
    }
  }
  estimate <- map_or_na(z_hat)
  lower <- map_or_na(endpoints[1])
  upper <- map_or_na(endpoints[2])

  out <- data.frame(
    parameter = link$names,
    estimate  = unname(estimate[link$names]),
    lower     = unname(lower[link$names]),
    upper     = unname(upper[link$names]),
    level     = level,
    method    = method,
    stringsAsFactors = FALSE
  )

  # Only for the profile path: the Wald path already warned about its own NAs.
  if (identical(method, "profile")) {
    if (anyNA(out$lower)) {
      warning("The profile did not cross the cutoff on the lower side, so ",
              "the lower endpoint is NA; the likelihood is flat in that ",
              "direction.", call. = FALSE)
    }
    if (anyNA(out$upper)) {
      warning("The profile did not cross the cutoff on the upper side, so ",
              "the upper endpoint is NA; the likelihood is flat in that ",
              "direction.", call. = FALSE)
    }
  }

  attr(out, "profile") <- profile
  out
}
