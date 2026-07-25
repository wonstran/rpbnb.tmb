# rpbnb.tmb

Maximum-simulated-likelihood estimation of bivariate random-parameter negative
binomial models with Famoye/Sarmanov or discrete-copula dependence, using TMB.

## Memory-aware fitting

The default fit keeps coefficient inference and the state required for
marginal effects, but releases the TMB objective and response vectors after
fitting:

```r
fit <- fit_rpbnb_tmb(
  y1 ~ x1 + x2,
  y2 ~ x1 + x2,
  data = dat,
  random_1 = "x1",
  draws = 400,
  inference = "full",
  keep = "postfit",
  control = rpbnb_tmb_control(
    n_cores = 4,
    max_threads = 4,
    max_workload = 2e6
  )
)
```

Use `inference = "diag"` when only coefficient standard errors are required.
This avoids retaining a full covariance matrix; `vcov(fit)` then has the
estimated diagonal and `NA` off-diagonal entries. Use `inference = "none"` to
skip Hessian construction entirely. Marginal-effect and elasticity
standard errors require `inference = "full"`; the lighter modes return their
point estimates with `NA` standard errors and an explicit warning.

Use `keep = "compact"` for the smallest returned object when marginal effects
are not needed. Use `keep = "full"` only for low-level diagnostics that require
`fit$obj` or the response vectors.

TMB tapes are constructed sequentially by default while objective and gradient
evaluation remain parallel, which lowers peak memory without changing fitted
results. Advanced users can opt into concurrent tape construction with
`parallel_tape = TRUE`; peak memory then scales with the thread count, so
`max_threads` caps per-fit concurrency.

`max_workload` rejects oversized observation-by-draw workloads before the
automatic-differentiation tape is built. One unit is one observation-draw.
Measured on a 12-point grid (n from 500 to 4000, draws from 50 to 200, single
threaded), tape size depends on `n * draws` alone to within 2.6%:

```
tape (MB) = 17.7 + 0.001117 * units          R^2 = 0.999
```

The slope is 1171 bytes per unit overall, rising to about 1215 over the largest
workloads — the conservative figure the default is set from. All dependence
families carry weight 1; the Gaussian copula tape is 7–9% *smaller* than
Famoye's at matched workload, because the atomic bivariate-normal kernel
compresses its quadrature. The default of `2e6` units admits fits whose tape
reaches roughly 2.3 GB, with objective and gradient evaluation adding about
another 55% for a working footprint near 3.5 GB. Raise it deliberately against
the memory you actually have; `max_workload = Inf` disables the guard.

## Boundary-constrained dependence estimates

Every dependence link is bounded — the Famoye admissible interval, Frank's
`|theta| < 35` (Kendall's tau up to about 0.891), Gaussian's `|rho| < 1`, and
Clayton's clamped exponential. When an estimate is pinned against one of these
bounds it is determined by the link rather than by the data, and its
delta-method standard error collapses toward zero. Such estimates are reported
with `NA` standard errors, listed in `fit$boundary_report`, and accompanied by
a warning. For Famoye fits, `fit$lambda_bounds` records the frozen admissible
interval so the constraint is inspectable. Treat a flagged fit as evidence that
the dependence family does not span the association in the data.

`predict(fit)` returns a two-column matrix of integrated fitted means.
`predict(fit, newdata = ..., type = "link")` supports new data and returns the
log of the integrated means.

The Gaussian copula's 20-point bivariate-normal quadrature is represented by a
thread-safe atomic TMB operation. Its value, gradient, and Hessian remain
consistent with the original quadrature while avoiding expansion of every
quadrature operation on the outer likelihood tape.
