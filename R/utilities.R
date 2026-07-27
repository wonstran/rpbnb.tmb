#' @keywords internal
#' @noRd
POISSON_M <- 1e-6

#' Ceiling of the bounded Frank link, shared with the C++ template
#'
#' The template maps the working parameter through
#' \code{theta = FRANK_THETA_MAX * tanh(z_dep / FRANK_THETA_MAX)} so that
#' \code{exp(-theta * u)} cannot overflow.  The value must match
#' \code{FRANK_THETA_MAX} in \code{src/rpbnb.tmb.cpp}.  It caps attainable
#' Frank dependence at Kendall's tau of about 0.891.
#' @keywords internal
#' @noRd
FRANK_THETA_MAX <- 35

#' Ceiling of the Clayton link, shared with the C++ template
#'
#' The template maps the working parameter through
#' \code{theta = exp(clamp(z_dep, -20, 20))} (\code{src/rpbnb.tmb.cpp}, the
#' \code{FAM_CLAYTON} branch), so a fitted natural-scale Clayton theta can be
#' as large as \code{exp(20)}, about 4.9e8. Simulation has to remain accurate
#' that far up, which is why its conditional inverse works in log space --
#' the direct \code{u^(-theta)} form overflows around theta = 700.
#' @keywords internal
#' @noRd
CLAYTON_THETA_MAX <- exp(20)

#' Measured tape-memory calibration
#'
#' The single source for every number in the `max_workload` safety contract:
#' the guard reads the family weights from here, `rpbnb_tmb_control()` derives
#' its default from `budget_gib` and `bytes_per_unit`, and the roxygen
#' documentation is generated from this object by `.calibration_doc()`.
#'
#' Two places restate these figures instead of deriving them, and each has its
#' own check that compares the restatement to this object and fails on drift:
#'
#'   * `README.md` -- checked by `tests/testthat/test-parallel.R`
#'     ("the README's restated calibration figures match the calibration").
#'   * `docs/reference/rpbnb.tmb-reference-manual.html` -- checked by
#'     `docs/reference/verify_reference_manual.R`.
#'
#' Update a restatement together with the constant, or add the new figure to
#' the corresponding check. Do not add a third restatement without one.
#'
#' Regenerate with `Rscript inst/benchmark_memory.R`, which writes the raw
#' measurements to `inst/extdata/memory_calibration.csv` and prints the
#' regression these constants come from. Re-run it after any template, TMB,
#' compiler or allocator change.
#'
#' @keywords internal
#' @noRd
TAPE_CALIBRATION <- list(
  # Retained tape, Famoye: tape_MiB = 13.374 + 0.0011332 * units, R^2 = 0.9994.
  tape_bytes_per_unit = 1221,
  tape_intercept_mib = 13.374,
  tape_r_squared = 0.9994,
  # PEAK working set is what actually causes std::bad_alloc -- the failure that
  # started this guard -- and it runs about 6.1x the retained tape PLUS the
  # growth across one fn()/gr() pair (see peak_over_tape_and_eval below for
  # why the denominator has to be named precisely), because TMB records the
  # full likelihood before pruning to each region. The budget is therefore set
  # on peak, not on tape, and this slope is regressed on peak_mib directly
  # rather than scaled up from the tape slope.
  peak_bytes_per_unit = 12083,
  # Median of peak_mib / (tape_mib + eval_mib) over the grid. The denominator
  # is retained tape PLUS the further growth across one fn()/gr() pair, not
  # tape alone -- against tape alone the same data gives 8.78. Named for what
  # it actually divides by, because the two are far enough apart to mislead.
  peak_over_tape_and_eval = 6.11,
  budget_gib = 8,
  # The largest PEAK ratio to Famoye observed at matched workload, rounded up
  # to the next tenth. Peak, not retained tape: peak is the quantity this
  # guard budgets, and the two disagree materially. Frank retains 2.87x Famoye
  # but peaks at 3.53x, so tape-derived weights under-budgeted Frank's peak by
  # about 22%. Famoye is the numeraire and is fixed at 1.
  family_weight = c(
    independence = 0.7, famoye = 1.0, frank = 3.6,
    gaussian = 0.9, clayton = 1.1
  ),
  measured_families = c(
    "independence", "famoye", "frank", "gaussian", "clayton"
  ),
  source = "inst/benchmark_memory.R"
)

#' Documentation text generated from TAPE_CALIBRATION
#' @keywords internal
#' @noRd
.calibration_doc <- function(calibration = TAPE_CALIBRATION) {
  weights <- calibration$family_weight
  paste0(
    "@param max_workload Maximum weighted observation-draw evaluations ",
    "permitted before TMB tape construction; \\code{Inf} disables the guard. ",
    "The default and every figure here are derived from ",
    "\\code{TAPE_CALIBRATION}, so this text cannot drift from the shipped ",
    "behaviour.\n\n",
    "With the default \\code{parallel_tape = FALSE} the budget is per fit; ",
    "with \\code{parallel_tape = TRUE} the tapes are built concurrently and ",
    "the guard multiplies the workload by the realized thread count.\n\n",
    "One unit is one weighted observation-draw. All figures are measured by ",
    "\\code{", calibration$source, "}, whose raw results are stored in ",
    "\\code{inst/extdata/memory_calibration.csv}.\n\n",
    "Retained tape size depends on \\code{n * draws} alone: ",
    "tape (MiB) = ", format(calibration$tape_intercept_mib, digits = 5),
    " + ", format(calibration$tape_bytes_per_unit / 1024^2, digits = 4),
    " * units, with R^2 = ", format(calibration$tape_r_squared, digits = 5),
    " (", calibration$tape_bytes_per_unit, " bytes per unit over the largest ",
    "workloads; per-unit cost is higher at small workloads because of the ",
    "fixed intercept).\n\n",
    "The budget is set on \\emph{peak} working set, not on retained tape, ",
    "because peak is what exhausts memory: TMB records the whole likelihood ",
    "before pruning it to each parallel region, so peak runs about ",
    calibration$peak_over_tape_and_eval,
    " times the retained tape plus first-evaluation growth, or ",
    calibration$peak_bytes_per_unit,
    " bytes per unit measured directly against peak. The default of ",
    format(.calibration_default_workload(calibration),
           scientific = FALSE, big.mark = ","),
    " units therefore targets a peak of about ", calibration$budget_gib,
    " GiB for one fit.\n\n",
    "Families cost different amounts per unit, so each carries a weight -- ",
    "the largest \\emph{peak} ratio to Famoye observed at matched workload, ",
    "rounded up to the next tenth: ",
    paste(sprintf("%s %g", names(weights), weights), collapse = ", "),
    ". Frank peaks at over three and a half times Famoye, so treating it as ",
    "unweighted would under-budget it more than threefold.\n\n",
    "Raise \\code{max_workload} deliberately against the memory you actually ",
    "have, not to make one particular dataset fit.\n\n",
    "As of \\code{rpbnb_tmb_max_workload()}, the value \\code{rpbnb_tmb_control()} ",
    "actually uses by default is no longer the fixed figure above: it will ",
    "auto-detect available memory and budget 80% of it, falling back to the ",
    "fixed 8 GiB figure only when detection is unavailable on the current ",
    "platform. Call \\code{\\link{rpbnb_tmb_max_workload}} directly to set ",
    "your own budget or detection fraction."
  )
}

#' Default workload budget implied by the calibration
#'
#' Derived, not restated: `rpbnb_tmb_control()` uses this for its default so
#' the published constants and the shipped behaviour cannot disagree.
#' @keywords internal
#' @noRd
.calibration_default_workload <- function(calibration = TAPE_CALIBRATION) {
  signif(calibration$budget_gib * 1024^3 / calibration$peak_bytes_per_unit, 1)
}

#' Available system memory, in GiB
#'
#' Best-effort and platform-specific. Every code path is wrapped so failure
#' -- a missing binary, a locale-mangled number, a sandbox that blocks
#' subprocess execution, an unrecognized OS -- degrades to `NA_real_` rather
#' than propagating an error. This return value feeds a default argument
#' (`rpbnb_tmb_control()`'s `max_workload`, by way of
#' `rpbnb_tmb_max_workload()`); a default argument that can error breaks
#' every caller who doesn't override it.
#' @keywords internal
#' @noRd
.detect_available_memory_gib <- function() {
  sysname <- tryCatch(Sys.info()[["sysname"]], error = function(e) NA_character_)
  if (is.na(sysname)) return(NA_real_)

  gib <- tryCatch({
    if (identical(sysname, "Linux")) {
      .detect_available_memory_gib_linux()
    } else if (identical(sysname, "Darwin")) {
      .detect_available_memory_gib_darwin()
    } else if (identical(sysname, "Windows")) {
      .detect_available_memory_gib_windows()
    } else {
      NA_real_
    }
  }, error = function(e) NA_real_)

  if (length(gib) != 1L || !is.numeric(gib) || is.na(gib) ||
      !is.finite(gib) || gib < 0) {
    return(NA_real_)
  }
  gib
}

#' @keywords internal
#' @noRd
.detect_available_memory_gib_linux <- function() {
  path <- "/proc/meminfo"
  if (!file.exists(path)) return(NA_real_)
  lines <- readLines(path, warn = FALSE)
  # MemAvailable is the kernel's own reclaimable-cache-aware estimate of what
  # can actually be handed to a new allocation; MemFree alone systematically
  # undercounts memory the kernel would reclaim on request. Fall back to
  # MemFree only on very old kernels that lack MemAvailable.
  value_kib <- .parse_meminfo_field(lines, "MemAvailable")
  if (is.na(value_kib)) value_kib <- .parse_meminfo_field(lines, "MemFree")
  if (is.na(value_kib)) return(NA_real_)
  value_kib / 1024^2
}

#' @keywords internal
#' @noRd
.parse_meminfo_field <- function(lines, field) {
  pattern <- paste0("^", field, ":\\s*([0-9]+)\\s*kB")
  hit <- grep(pattern, lines, value = TRUE)
  if (!length(hit)) return(NA_real_)
  as.numeric(sub(pattern, "\\1", hit[1L]))
}

#' @keywords internal
#' @noRd
.detect_available_memory_gib_darwin <- function() {
  vm <- tryCatch(system2("vm_stat", stdout = TRUE, stderr = FALSE),
                 error = function(e) character(0))
  if (!length(vm)) return(NA_real_)

  page_size_line <- grep("page size of", vm, value = TRUE)
  page_size <- if (length(page_size_line)) {
    as.numeric(sub(".*page size of ([0-9]+) bytes.*", "\\1", page_size_line[1L]))
  } else {
    4096  # Documented default when the header line's wording ever changes.
  }
  if (is.na(page_size) || page_size <= 0) return(NA_real_)

  free_pages <- .parse_vm_stat_field(vm, "Pages free")
  inactive_pages <- .parse_vm_stat_field(vm, "Pages inactive")
  if (is.na(free_pages) || is.na(inactive_pages)) return(NA_real_)

  (free_pages + inactive_pages) * page_size / 1024^3
}

#' @keywords internal
#' @noRd
.parse_vm_stat_field <- function(lines, field) {
  pattern <- paste0("^", field, ":\\s*([0-9]+)\\.?")
  hit <- grep(pattern, lines, value = TRUE)
  if (!length(hit)) return(NA_real_)
  as.numeric(sub(pattern, "\\1", hit[1L]))
}

#' @keywords internal
#' @noRd
.detect_available_memory_gib_windows <- function() {
  gib <- .detect_available_memory_gib_windows_powershell()
  if (!is.na(gib)) return(gib)
  .detect_available_memory_gib_windows_wmic()
}

# PowerShell's CIM cmdlets have been available since Windows 8 / Server 2012
# and, unlike wmic.exe, are not deprecated -- Windows 11 24H2 removed wmic.exe
# from the default install entirely, which left this detector permanently
# blind on current systems. Tried first for that reason; wmic below is kept
# only as a fallback for older or stripped-down environments.
#' @keywords internal
#' @noRd
.detect_available_memory_gib_windows_powershell <- function() {
  out <- tryCatch(
    system2(
      "powershell",
      c("-NoProfile", "-Command",
        "(Get-CimInstance Win32_OperatingSystem).FreePhysicalMemory"),
      stdout = TRUE, stderr = FALSE
    ),
    error = function(e) character(0)
  )
  if (!length(out)) return(NA_real_)
  value_kib <- suppressWarnings(as.numeric(trimws(out[1L])))
  if (is.na(value_kib)) return(NA_real_)
  value_kib / 1024^2
}

#' @keywords internal
#' @noRd
.detect_available_memory_gib_windows_wmic <- function() {
  out <- tryCatch(
    system2("wmic", c("OS", "get", "FreePhysicalMemory", "/value"),
            stdout = TRUE, stderr = FALSE),
    error = function(e) character(0)
  )
  if (!length(out)) return(NA_real_)
  hit <- grep("^FreePhysicalMemory=", out, value = TRUE)
  if (!length(hit)) return(NA_real_)
  value_kib <- suppressWarnings(
    as.numeric(sub("^FreePhysicalMemory=", "", hit[1L]))
  )
  if (is.na(value_kib)) return(NA_real_)
  value_kib / 1024^2
}

#' Compute a TMB workload budget from a memory figure
#'
#' Companion to \code{rpbnb_tmb_control()}'s \code{max_workload}: rather than
#' picking a weighted-observation-draw count directly, state a memory budget
#' and let this function do the arithmetic \code{TAPE_CALIBRATION} implies.
#'
#' With \code{budget_gib} omitted, this is also \code{rpbnb_tmb_control()}'s
#' own default: every fit that doesn't set \code{max_workload} explicitly
#' already goes through this function.
#' @param budget_gib Optional memory budget in GiB, stated explicitly. When
#'   supplied, used as-is -- \code{fraction} does not apply, because a number
#'   you state yourself is not second-guessed with a discount. When omitted
#'   (the default), available memory is auto-detected and \code{fraction} of
#'   it is budgeted instead.
#' @param fraction Of auto-detected available memory, the fraction to
#'   actually budget; one number in \code{(0, 1]}. Available memory
#'   fluctuates and competes with other processes, so budgeting all of it
#'   risks the guard passing a fit that then exhausts memory anyway. Ignored
#'   when \code{budget_gib} is supplied.
#' @return One positive numeric workload value, on the same scale as
#'   \code{rpbnb_tmb_control()}'s \code{max_workload}.
#' @export
#' @examples
#' rpbnb_tmb_max_workload(budget_gib = 16)
#' \dontrun{
#' ctrl <- rpbnb_tmb_control(max_workload = rpbnb_tmb_max_workload())
#' }
rpbnb_tmb_max_workload <- function(budget_gib = NULL, fraction = 0.8) {
  if (length(fraction) != 1L || !is.numeric(fraction) || is.na(fraction) ||
      !is.finite(fraction) || fraction <= 0 || fraction > 1) {
    stop("fraction must be one number in (0, 1].", call. = FALSE)
  }

  if (!is.null(budget_gib)) {
    if (length(budget_gib) != 1L || !is.numeric(budget_gib) ||
        is.na(budget_gib) || !is.finite(budget_gib) || budget_gib <= 0) {
      stop("budget_gib must be one positive finite number, or NULL.",
           call. = FALSE)
    }
    return(signif(
      budget_gib * 1024^3 / TAPE_CALIBRATION$peak_bytes_per_unit, 1
    ))
  }

  detected_gib <- .detect_available_memory_gib()
  if (is.na(detected_gib)) {
    warning(
      "Could not detect available memory on this platform; using the ",
      "default 8 GiB calibration budget. Pass `budget_gib` explicitly to ",
      "set your own.",
      call. = FALSE
    )
    return(.calibration_default_workload())
  }

  signif(
    fraction * detected_gib * 1024^3 / TAPE_CALIBRATION$peak_bytes_per_unit, 1
  )
}

#' @keywords internal
#' @noRd
.chk_poisson_flag <- function(x, arg) {
  if (!(is.logical(x) && length(x) == 1L && !is.na(x))) {
    stop("`", arg, "` must be a single non-missing logical.")
  }
}

#' @keywords internal
#' @noRd
.check_counts <- function(y, label) {
  y <- as.numeric(y)
  if (any(!is.finite(y))) stop("Response ", label, " has non-finite values.")
  if (any(y < 0)) stop("Response ", label, " has negative values.")
  if (any(abs(y - round(y)) > 1e-8)) stop("Response ", label, " has non-integer values.")
  as.integer(round(y))
}

#' @keywords internal
#' @noRd
.prepare_bnb_data <- function(formula_1, formula_2, data) {
  vars <- unique(c(all.vars(stats::terms(formula_1, data = data)),
                   all.vars(stats::terms(formula_2, data = data))))
  missing_vars <- vars[!vars %in% names(data)]
  if (length(missing_vars)) stop("Variables not found: ", paste(missing_vars, collapse = ", "))
  ok <- stats::complete.cases(data[, vars, drop = FALSE])
  if (!any(ok)) stop("No complete cases.")
  data_cc <- data[ok, , drop = FALSE]
  # Both equations must be built from ONE row mask. Completeness of the raw
  # variables is not enough: a transformation in the formula can manufacture
  # NA/NaN in one equation only (`log(x)` with a non-positive x), and letting
  # each model.frame() drop its own rows silently pairs equation 1's
  # observation i with equation 2's observation j. The template takes n from
  # Y1 alone and indexes Y2/X2 with it, so nothing downstream would catch it.
  # Evaluate both frames with na.pass, intersect the surviving rows, then
  # rebuild on that shared subset so the terms attributes stay intact.
  mf1 <- stats::model.frame(formula_1, data = data_cc,
                            na.action = stats::na.pass)
  mf2 <- stats::model.frame(formula_2, data = data_cc,
                            na.action = stats::na.pass)
  shared <- stats::complete.cases(mf1) & stats::complete.cases(mf2)
  if (!any(shared)) stop("No complete cases.")
  if (!all(shared)) {
    data_cc <- data_cc[shared, , drop = FALSE]
    mf1 <- stats::model.frame(formula_1, data = data_cc,
                              na.action = stats::na.pass)
    mf2 <- stats::model.frame(formula_2, data = data_cc,
                              na.action = stats::na.pass)
  }
  Y1 <- .check_counts(stats::model.response(mf1), "1")
  Y2 <- .check_counts(stats::model.response(mf2), "2")
  X1 <- stats::model.matrix(formula_1, mf1)
  X2 <- stats::model.matrix(formula_2, mf2)
  # Cheap, and the failure it guards against is silent misalignment rather
  # than an error, so assert rather than trust the construction above.
  if (length(Y1) != length(Y2) || nrow(X1) != nrow(X2) ||
      length(Y1) != nrow(X1) || length(Y1) != nrow(data_cc)) {
    stop("The two equations resolved to different numbers of rows; this is ",
         "an internal error in row alignment.", call. = FALSE)
  }
  model_meta <- list(
    eq1 = list(
      terms = stats::delete.response(stats::terms(mf1)),
      xlevels = lapply(mf1[vapply(mf1, is.factor, logical(1))], levels),
      contrasts = attr(X1, "contrasts")
    ),
    eq2 = list(
      terms = stats::delete.response(stats::terms(mf2)),
      xlevels = lapply(mf2[vapply(mf2, is.factor, logical(1))], levels),
      contrasts = attr(X2, "contrasts")
    )
  )
  list(Y1 = Y1, Y2 = Y2, X1 = X1, X2 = X2,
       cn1 = colnames(X1), cn2 = colnames(X2),
       n = length(Y1), data = data_cc, model_meta = model_meta)
}

#' @keywords internal
#' @noRd
rand_dist_registry <- list(
  normal = list(
    base = "normal",
    u_to_base = function(u) stats::qnorm(u),
    coef = function(b, s, base, sign) b + s * base,
    dev = function(b, s, base, sign) s * base,
    dloc_factor = function(b, s, base, coef) rep(1, length(base)),
    dscale = function(b, s, base, coef) s * base,
    scale_label = "log_sd"
  ),
  lognormal = list(
    base = "normal",
    u_to_base = function(u) stats::qnorm(u),
    coef = function(b, s, base, sign) sign * exp(b + s * base),
    dev = function(b, s, base, sign) sign * exp(b + s * base) - b,
    dloc_factor = function(b, s, base, coef) coef,
    dscale = function(b, s, base, coef) coef * base * s,
    scale_label = "log_s"
  ),
  uniform = list(
    base = "uniform",
    u_to_base = function(u) u,
    coef = function(b, s, base, sign) b + s * (2 * base - 1),
    dev = function(b, s, base, sign) s * (2 * base - 1),
    dloc_factor = function(b, s, base, coef) rep(1, length(base)),
    dscale = function(b, s, base, coef) s * (2 * base - 1),
    scale_label = "log_w"
  ),
  triangular = list(
    base = "uniform",
    u_to_base = function(u) tri_icdf(u),
    coef = function(b, s, base, sign) b + s * base,
    dev = function(b, s, base, sign) s * base,
    dloc_factor = function(b, s, base, coef) rep(1, length(base)),
    dscale = function(b, s, base, coef) s * base,
    scale_label = "log_w"
  )
)

tri_icdf <- function(u) {
  ifelse(u < 0.5, -1 + sqrt(2 * u), 1 - sqrt(2 * (1 - u)))
}

parse_rand_spec <- function(spec) {
  if (is.null(spec) || length(spec) == 0) {
    return(list(names = character(0), dist = character(0),
                sign = numeric(0), scale = numeric(0)))
  }
  valid <- names(rand_dist_registry)
  if (is.character(spec) && is.null(names(spec))) {
    return(list(names = spec, dist = rep("normal", length(spec)),
                sign = rep(1, length(spec)), scale = rep(NA_real_, length(spec))))
  }
  if (!is.list(spec) || is.null(names(spec)) || any(!nzchar(names(spec)))) {
    stop("random spec must be a character vector or named list.")
  }
  nm <- names(spec)
  dist <- character(length(nm))
  sgn <- numeric(length(nm))
  scale <- numeric(length(nm))
  for (i in seq_along(nm)) {
    v <- spec[[i]]
    if (is.character(v) && length(v) == 1L) {
      d <- v; this_sign <- 1; this_scale <- NA_real_
    } else if (is.list(v)) {
      d <- if (is.null(v$dist)) "normal" else v$dist
      this_sign <- if (is.null(v$sign)) 1 else v$sign
      this_scale <- if (!is.null(v$scale)) v$scale else if (!is.null(v$sd)) v$sd else NA_real_
      if (!is.null(v$sign) && d != "lognormal") {
        stop("`sign` is only for lognormal.")
      }
    } else {
      stop("Invalid random spec for '", nm[i], "'.")
    }
    if (!d %in% valid) stop("Unknown distribution '", d, "'.")
    if (!this_sign %in% c(-1, 1)) stop("`sign` must be -1 or 1.")
    dist[i] <- d; sgn[i] <- this_sign; scale[i] <- this_scale
  }
  list(names = nm, dist = dist, sign = sgn, scale = scale)
}

chk_rand_spec <- function(spec, bv, lbl) {
  miss <- spec$names[!spec$names %in% names(bv)]
  if (length(miss)) stop("random name(s) ", paste(miss, collapse = ", "), " not in ", lbl, ".")
  if (length(spec$names)) {
    if (any(is.na(spec$scale))) stop("Each random coefficient needs a `scale`/`sd` in ", lbl, ".")
    if (any(!is.finite(spec$scale)) || any(spec$scale <= 0)) {
      stop("Random scales must be finite and positive in ", lbl, ".")
    }
  }
}

chk_dispersion <- function(dispersion) {
  if (!all(c("m1", "m2") %in% names(dispersion))) {
    stop("`dispersion` must be c(m1 = ., m2 = .).")
  }
  m <- dispersion[c("m1", "m2")]
  if (any(!is.finite(m)) || any(m <= 0)) stop("Dispersions must be finite and positive.")
}

rand_realize <- function(U, dist, sign, b, s) {
  R <- nrow(U); q <- ncol(U)
  base <- coef <- dev <- dloc <- dscale <- matrix(0, nrow = R, ncol = q)
  for (j in seq_len(q)) {
    reg <- rand_dist_registry[[dist[j]]]
    bj <- reg$u_to_base(U[, j])
    cj <- reg$coef(b[j], s[j], bj, sign[j])
    base[, j] <- bj
    coef[, j] <- cj
    dev[, j] <- reg$dev(b[j], s[j], bj, sign[j])
    dloc[, j] <- reg$dloc_factor(b[j], s[j], bj, cj)
    dscale[, j] <- reg$dscale(b[j], s[j], bj, cj)
  }
  list(base = base, coef = coef, dev = dev, dloc = dloc, dscale = dscale)
}

d_const <- function() 1 - exp(-1)

c_val <- function(mu, m) {
  if (m == 0) return(exp(-d_const() * mu))
  (1 + d_const() * m * mu)^(-1 / m)
}

lambda_bounds_vec <- function(c1, c2) {
  lam_min <- -1 / pmax((1 - c1) * (1 - c2), c1 * c2)
  lam_max <-  1 / pmax(c1 * (1 - c2), c2 * (1 - c1))
  c(max(lam_min), min(lam_max))
}

.resolve_start <- function(start, default, par_names, label = "start") {
  if (is.null(start)) { names(default) <- par_names; return(default) }
  if (!is.numeric(start)) stop("`", label, "` must be numeric.")
  nm <- names(start)
  if (is.null(nm) || all(!nzchar(nm))) {
    if (length(start) != length(par_names)) {
      stop("`", label, "` must have length ", length(par_names), ".")
    }
    names(start) <- par_names
    return(start)
  }
  if (any(!nzchar(nm))) stop("Mix of named/unnamed in `", label, "`.")
  if (anyDuplicated(nm)) stop("Duplicate names in `", label, "`.")
  unknown <- setdiff(nm, par_names)
  if (length(unknown)) stop("Unknown names in `", label, "`: ", paste(unknown, collapse = ", "))
  out <- default; names(out) <- par_names
  out[nm] <- start
  out
}

#' Map a dependence specification to its family code
#'
#' The single definition shared by the fitter and by
#' \code{rpbnb_tmb_dependence_profile()}. Keeping one copy is the point: a
#' second inline \code{switch} would silently disagree the moment a family is
#' added.
#' @keywords internal
#' @noRd
.resolve_family_code <- function(dependence) {
  if (inherits(dependence, "rpbnb_copula")) {
    return(switch(dependence$family,
                  frank = 1L, normal = 2L, kimeldorf = 3L))
  }
  if (identical(dependence, "famoye")) return(0L)
  if (identical(dependence, "independence")) return(-1L)
  stop("`dependence` must be \"famoye\", \"independence\", or copula().",
       call. = FALSE)
}

#' Specify a copula dependence model
#'
#' Creates a copula specification for \code{fit_rpbnb_tmb()} or
#' \code{simulate_rpbnb_tmb()}.
#'
#' \code{par} is a \strong{simulation-only} argument. \code{fit_rpbnb_tmb()}
#' always estimates the dependence parameter and reads only \code{family}, so
#' it rejects a non-\code{NULL} \code{par} rather than accepting and ignoring
#' it; use \code{start} there to set a working-scale starting value.
#' \code{simulate_rpbnb_tmb()} uses independence when \code{par} is omitted.
#'
#' @param family Copula family: \code{"frank"}, \code{"normal"} (Gaussian), or
#'   \code{"kimeldorf"} (Clayton).
#' @param par Optional natural-scale dependence parameter for
#'   \emph{simulation}: Frank theta, Gaussian correlation rho, or positive
#'   Clayton theta. Not accepted by \code{fit_rpbnb_tmb()}.
#' @return An object of class \code{rpbnb_copula}.
#' @export
copula <- function(family, par = NULL) {
  family <- match.arg(family, c("frank", "normal", "kimeldorf"))
  if (!is.null(par)) {
    if (!is.numeric(par) || length(par) != 1L || !is.finite(par)) {
      stop("`par` must be one finite numeric value.", call. = FALSE)
    }
    if (family == "normal" && abs(par) >= 1) stop("|rho| must be < 1.")
    if (family == "kimeldorf" && par <= 0) stop("Clayton theta must be > 0.")
  }
  structure(list(family = family, par = par), class = "rpbnb_copula")
}

#' Control parameters for rpbnb_tmb estimators
#' @param iterlim Maximum optimizer iterations.
#' @param reltol Relative convergence tolerance.
#' @param print_level Verbosity for nlminb.
#' @param n_cores Number of OpenMP threads for TMB.
#' @param max_threads Maximum OpenMP threads permitted for one fit. Defaults
#'   to \code{n_cores}, so threads are not capped unless set explicitly below
#'   \code{n_cores}.
#' @eval .calibration_doc()
#' @param parallel_tape Construct per-thread TMB tapes concurrently. The
#'   default \code{FALSE} constructs them sequentially to reduce peak memory;
#'   objective and gradient evaluation remains parallel.
#' @param halton_burn Number of leading Halton points to discard.
#' @export
rpbnb_tmb_control <- function(iterlim = 500L,
                              reltol = 1e-8,
                              print_level = 0L,
                              n_cores = 1L,
                              max_threads = n_cores,
                              max_workload =
                                rpbnb_tmb_max_workload(),
                              parallel_tape = FALSE,
                              halton_burn = 300L) {
  if (length(n_cores) != 1L || !is.numeric(n_cores) ||
      is.na(n_cores) || !is.finite(n_cores) ||
      n_cores < 1 || n_cores != floor(n_cores) ||
      n_cores > .Machine$integer.max) {
    stop("n_cores must be one whole number greater than or equal to 1.",
         call. = FALSE)
  }
  if (length(max_threads) != 1L || !is.numeric(max_threads) ||
      is.na(max_threads) || !is.finite(max_threads) ||
      max_threads < 1 || max_threads != floor(max_threads) ||
      max_threads > .Machine$integer.max) {
    stop("max_threads must be one whole number greater than or equal to 1.",
         call. = FALSE)
  }
  if (length(max_workload) != 1L || !is.numeric(max_workload) ||
      is.na(max_workload) || max_workload <= 0) {
    stop("max_workload must be one positive number or Inf.", call. = FALSE)
  }
  if (!is.logical(parallel_tape) || length(parallel_tape) != 1L ||
      is.na(parallel_tape)) {
    stop("parallel_tape must be one non-missing logical value.",
         call. = FALSE)
  }
  structure(list(iterlim = as.integer(iterlim), reltol = reltol,
                  print_level = as.integer(print_level),
                  n_cores = as.integer(n_cores),
                  max_threads = as.integer(max_threads),
                  max_workload = as.numeric(max_workload),
                  parallel_tape = parallel_tape,
                  halton_burn = as.integer(halton_burn)),
            class = "rpbnb_tmb_control")
}
