#' Configure the TMB model DLL's OpenMP thread count
#' @keywords internal
#' @noRd
.configure_tmb_threads <- function(n_cores, max_threads = 4L,
                                   parallel_tape = FALSE,
                                   DLL = "rpbnb.tmb") {
  supported <- suppressWarnings(TMB::openmp(max = TRUE, DLL = DLL))
  supported <- as.integer(supported[[1L]])
  if (length(supported) != 1L || is.na(supported) || supported < 1L) {
    supported <- 1L
  }

  requested <- as.integer(n_cores)
  policy_capped <- min(requested, as.integer(max_threads))
  realized <- min(policy_capped, supported)
  TMB::openmp(n = realized, DLL = DLL)
  TMB::config(
    tape.parallel = as.integer(isTRUE(parallel_tape)),
    DLL = DLL
  )

  if (policy_capped < requested) {
    warning(
      sprintf(
        "Requested %d TMB threads; using %d because max_threads = %d.",
        requested, realized, as.integer(max_threads)
      ),
      call. = FALSE
    )
  } else if (realized < requested) {
    warning(
      sprintf("Requested %d TMB threads; using %d supported thread%s.",
              requested, realized, if (realized == 1L) "" else "s"),
      call. = FALSE
    )
  }

  realized
}

#' Construct the RP-BNB TMB objective with an explicit thread setting
#' @keywords internal
#' @noRd
.make_rpbnb_tmb_object <- function(data, parameters, map = NULL,
                                   random = NULL, silent = TRUE,
                                   n_cores = 1L, max_threads = 4L,
                                   parallel_tape = FALSE,
                                   DLL = "rpbnb.tmb") {
  realized <- .configure_tmb_threads(
    n_cores, max_threads = max_threads,
    parallel_tape = parallel_tape, DLL = DLL
  )
  obj <- TMB::MakeADFun(
    data = data,
    parameters = parameters,
    map = map,
    random = random,
    DLL = DLL,
    silent = silent
  )
  list(obj = obj, n_cores = realized)
}

#' Reject automatic-differentiation workloads above an explicit budget
#' @keywords internal
#' @noRd
.check_tmb_workload <- function(n, draws, family_code, max_workload,
                                n_threads = 1L, parallel_tape = FALSE) {
  if (is.infinite(max_workload)) return(invisible(0))
  # Weights are measured, not assumed: see TAPE_CALIBRATION and
  # inst/benchmark_memory.R.  Frank costs ~2.9x Famoye per unit; Gaussian and
  # Clayton are close to 1; independence is cheaper.  A previous revision
  # generalised "weight 1" to every family from Famoye and Gaussian
  # measurements alone, which under-budgeted Frank threefold.
  family_weight <- unname(TAPE_CALIBRATION$family_weight[[
    match(family_code, c(-1L, 0L, 1L, 2L, 3L))
  ]])
  tape_multiplier <- if (isTRUE(parallel_tape)) as.double(n_threads) else 1
  workload <- as.double(n) * as.double(draws) *
    family_weight * tape_multiplier
  if (!is.finite(workload) || workload > max_workload) {
    stop(
      sprintf(
        paste0(
          "Weighted TMB workload is %s, above max_workload = %s. ",
          "Reduce observations or draws, or explicitly increase ",
          "control$max_workload (Inf disables this guard)."
        ),
        format(workload, scientific = FALSE, trim = TRUE),
        format(max_workload, scientific = FALSE, trim = TRUE)
      ),
      call. = FALSE
    )
  }
  invisible(workload)
}

#' Build the TMB data list for the RP-BNB model
#' @keywords internal
#' @noRd
.build_tmb_data <- function(Y1, Y2, X1, X2, rand_idx1, rand_idx2,
                            Z1, Z2, dist1, dist2, sign1, sign2,
                            family_code, pois1, pois2,
                            lamLo, lamHi) {
  list(
    Y1 = Y1, Y2 = Y2,
    X1 = unname(as.matrix(X1)), X2 = unname(as.matrix(X2)),
    rand_idx1 = as.integer(rand_idx1) - 1L,  # 0-based for C++
    rand_idx2 = as.integer(rand_idx2) - 1L,
    Z1 = unname(as.matrix(Z1)), Z2 = unname(as.matrix(Z2)),
    dist1 = as.integer(dist1), dist2 = as.integer(dist2),
    sign1 = as.integer(sign1), sign2 = as.integer(sign2),
    family = as.integer(family_code),
    pois1 = as.integer(pois1), pois2 = as.integer(pois2),
    lamLo = as.numeric(lamLo), lamHi = as.numeric(lamHi)
  )
}
