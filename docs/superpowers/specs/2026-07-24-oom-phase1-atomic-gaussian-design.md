# OOM Phase 1 and Atomic Gaussian Kernel Design

## Goal

Reduce peak and retained memory during RP-BNB fitting without changing the
simulated-likelihood estimator, and compress the Gaussian copula's repeated
bivariate-normal quadrature on the automatic-differentiation tape.

## Public Interface

`fit_rpbnb_tmb()` gains two matched arguments:

- `inference = c("full", "diag", "none")`
- `keep = c("postfit", "compact", "full")`

`full` computes and stores the full free-parameter covariance. `diag` computes
standard errors without retaining off-diagonal covariance. `none` skips
Hessian-based inference. Natural-scale dispersion and dependence estimates are
available in every mode.

`postfit` keeps the design matrices and random-draw metadata needed by marginal
effects, but releases the TMB objective and response vectors. `compact` also
releases post-fit design and draw data. `full` retains all current diagnostic
objects. Defaults are `inference = "full"` and `keep = "postfit"` so existing
coefficient and marginal-effect behavior remains available while the large TMB
objective is released.

The `sdreport` field becomes a compact, package-owned compatibility object with
`pdHess`, fixed-parameter and natural-scale report summaries. Calling
`summary(fit$sdreport, "fixed")` or `summary(fit$sdreport, "report")` remains
supported without asking TMB to retape the likelihood.

## Workload and Thread Guardrails

`rpbnb_tmb_control()` gains:

- `max_threads = 4L`
- `max_workload = 2e6`

The configured TMB threads are capped at `min(n_cores, max_threads)` with a
warning when capped. Users can opt into more threads by increasing
`max_threads`.

Before constructing the TMB object, the fit computes weighted work as
`n * effective_draws * family_weight`, where weights are 1 for independence
and Famoye and 4 for each copula family after atomic compression. Work above
`max_workload` stops before taping and names all relevant controls.
`max_workload = Inf` explicitly disables the guard.

## Inference Architecture

Inference uses the already-recorded objective:

1. Evaluate `optimHess(opt$par, obj$fn, obj$gr)` once.
2. Symmetrize and Cholesky-check the Hessian.
3. For `full`, solve for the full covariance.
4. For `diag`, solve one unit-vector system at a time and retain only inverse
   diagonal entries.
5. For `none`, return NA standard errors without evaluating a Hessian.

Natural-scale reports and their one-parameter delta-method derivatives are
calculated in R. Frank Kendall's tau and its derivative are evaluated using the
same fixed Gauss-Legendre rule as the C++ implementation. No call to
`TMB::sdreport()` is made.

`vcov(fit)` materializes a named matrix on request. In diagonal mode the
off-diagonal entries are `NA`; in no-inference mode all entries are `NA`.

## Streaming Marginal Effects

Marginal-effect means are accumulated one draw at a time. The implementation
does not allocate `n_observations * n_draws` matrices in either the point
estimate or numerical-Jacobian path. Compact fits fail early with a message
that requests `keep = "postfit"` or `"full"`.

## Atomic Gaussian Kernel

The existing 20-point bivariate-normal CDF quadrature is moved to a standalone
templated function. A known-derivative `TMB_ATOMIC_VECTOR_FUNCTION` accepts
`(h, k, rho)` and returns the same quadrature result. Its reverse rule
differentiates the implemented quadrature, including its conditional
integration limit, so its gradient is consistent with its value. Each outer
likelihood call occupies one atomic tape node instead of expanding all
quadrature operations.

`REGISTER_ATOMIC` was rejected during implementation because its static
checkpoint storage is sized from the thread count on first use. Constructing a
one-thread fit before a two-thread fit in the same R process could then access
missing per-thread storage.

Serial and OpenMP objectives, gradients, and Hessians must agree. The atomic
result and derivatives must also agree with the pre-existing numerical
references.

## Compatibility and Failure Behavior

- Existing coefficient names, likelihood values, summaries, and full-mode
  covariance semantics remain unchanged.
- Code that directly requires `fit$obj` must request `keep = "full"`.
- Marginal effects in `diag` mode return point estimates and NA standard errors
  because covariance cross-terms are unavailable.
- Allocation and non-positive-definite Hessian failures produce explicit
  diagnostics rather than silently constructing repaired covariance matrices.

## Verification

Run focused red-green tests for controls, inference modes, retention modes,
streaming marginal effects, and atomic serial/parallel derivatives. Reinstall
the package and run the complete installed-package suite. The two
`test-against-rpbnb.R` failures caused by the malformed external
`rpbnb::rwm1984_clean` fixture are recorded separately and are not treated as
regressions.
