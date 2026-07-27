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
    max_workload = 7e5
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
automatic-differentiation tape is built. Every constant behind it is measured
by `inst/benchmark_memory.R`, whose raw results are committed to
`inst/extdata/memory_calibration.csv`; `TAPE_CALIBRATION` in `R/utilities.R` is
the single source that the guard, the default, and the generated `?` help all
derive from. The figures restated on this page and in
`docs/reference/rpbnb.tmb-reference-manual.html` are *copies*, not derivations
— both are checked against `TAPE_CALIBRATION` (by `test-parallel.R` and by
`docs/reference/verify_reference_manual.R` respectively) so a copy that drifts
fails rather than misleads.

The budget is set on **peak** working set rather than retained tape, because
peak is what exhausts memory: TMB records the whole likelihood before pruning
it to each parallel region, so peak runs about 6.1x the retained tape *plus
first-evaluation growth* — roughly 12 kB per weighted observation-draw,
measured directly against peak. The default of `7e5` units therefore targets
about 8 GiB of peak memory for one fit.

All of the above bounds the *simulated* likelihood, whose tape scales with
`nrow(data) * draws`. `method = "laplace"` integrates the random coefficients
with TMB's Laplace approximation instead, taping one conditional evaluation per
observation and integrating the latents through a sparse Hessian. Tape size then
scales with `nrow(data)` alone and the draw budget stops binding. That is the
option to reach for when a fit exhausts memory during `MakeADFun()`, rather than
raising `max_workload` against RAM you do not have. It supports normal and
lognormal random coefficients, and it is a different approximation to the same
integral — see `?fit_rpbnb_tmb`.

Families differ, so each carries a measured weight: the largest **peak** ratio
to Famoye at matched workload, rounded up to the next tenth. Peak rather than
retained tape, because peak is the quantity being budgeted — the two disagree
enough to matter.

| family | weight | largest peak ratio | largest tape ratio |
|---|---|---|---|
| independence | 0.7 | 0.602 | 0.647 |
| famoye | 1.0 | — | — |
| frank | **3.6** | 3.530 | 2.867 |
| gaussian | 0.9 | 0.844 | 0.947 |
| clayton | 1.1 | 1.079 | 0.921 |

Frank peaks at over three and a half times Famoye per unit, so a Frank fit buys
proportionally fewer draws for the same memory. Raise `max_workload`
deliberately against the memory you actually have —
`inst/fit_rpbnb_diff_copula.R` shows that opt-in — and `max_workload = Inf`
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
