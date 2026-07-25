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
    max_workload = 2.5e6
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
automatic-differentiation tape is built. One weighted unit is one
observation-draw of a Famoye or independence tape, measured at roughly 850
bytes; copula families carry a weight of 1.25. The default of `2.5e6` units
therefore encodes a budget of about 2.1 GB of tape per fit. Raise it
deliberately against the memory you actually have; `max_workload = Inf`
disables the guard.

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
