#' Load and compile the TMB template if needed
#' @keywords internal
#' @noRd
.load_tmb_dll <- function() {
  pkg_path <- find.package("rpbnb.tmb")
  TMB::dynlib(file.path(pkg_path, "src", "rpbnb.tmb"))
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
