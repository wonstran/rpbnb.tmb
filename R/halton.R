#' Generate uniform Halton draws with Cranley-Patterson randomization
#'
#' @param n_draws Number of draws.
#' @param d Dimension (columns). If 0, returns 0-column matrix.
#' @param burn Number of leading points to discard.
#' @return n_draws x d matrix of uniforms in (1e-12, 1-1e-12).
#' @keywords internal
#' @noRd
halton_uniform <- function(n_draws, d, burn = 300L) {
  if (d <= 0L) return(matrix(0, nrow = n_draws, ncol = 0L))
  # Prime bases: first d primes
  primes <- c(2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47,
              53, 59, 61, 67, 71, 73, 79, 83, 89, 97)
  if (d > length(primes)) stop("d > ", length(primes), " not supported.")
  n <- burn + n_draws
  U <- matrix(0, nrow = n, ncol = d)
  for (j in seq_len(d)) {
    b <- primes[j]
    U[, j] <- .halton1d(n, b)
  }
  # Cranley-Patterson rotation
  shift <- runif(d)
  U <- sweep(U, 2L, shift, `+`)
  U <- U - floor(U)
  # Clamp away from open endpoints (TMB qnorm needs strict interior)
  U <- U[(burn + 1L):n, , drop = FALSE]
  pmax(pmin(U, 1 - 1e-12), 1e-12)
}

#' 1D Halton sequence
#' @keywords internal
#' @noRd
.halton1d <- function(n, b) {
  # Returns n Halton points in base b in (0, 1)
  v <- numeric(n)
  # Use the radical-inverse construction
  # This is a simple loop; for large n a vectorized approach is used
  # b^k expansion: x = sum a_k * b^{-k-1}
  inv_b <- 1.0 / b
  for (i in seq_len(n)) {
    ii <- i
    f <- inv_b
    x <- 0.0
    while (ii > 0) {
      digit <- ii %% b
      x <- x + digit * f
      ii <- ii %/% b
      f <- f * inv_b
    }
    v[i] <- x
  }
  v
}
