# RP-BNB TMB Implementation Design

## Overview

Implement the random-parameter bivariate negative binomial (RP-BNB) model using the **TMB** (Template Model Builder) R package as the sole computational backend. This is a new package `rpbnb_tmb` that fully replaces the existing `rpbnb` package's R+Rcpp/OpenMP C++ backend with TMB's CppAD automatic differentiation.

**Key advantage:** TMB's AD eliminates the need for hand-coded analytic gradients (Famoye path) and per-draw score scalars (copula path), dramatically reducing the maintenance burden while providing numerically exact derivatives.

## Dependencies

- **R (>= 4.1)**
- **TMB** (LinkingTo + Imports) — the template compiler and AD engine
- **stats** — `nlminb` optimizer, distribution functions
- **parallel** — optional, for cluster-based Halton generation

No Rcpp, no maxLik, no numDeriv, no pbivnorm, no MASS.

## Architecture

### Package Structure

```
rpbnb_tmb/
├── DESCRIPTION
├── NAMESPACE
├── src/
│   ├── rpbnb_tmb.cpp        # Single TMB template (Famoye + Copula paths)
│   └── Makevars              # TMB linking flags
├── R/
│   ├── rpbnb_tmb-package.R   # Package-level docs
│   ├── fit_rpbnb_tmb.R       # Main estimation function
│   ├── simulate_rpbnb_tmb.R  # Data simulator (pure R, no TMB needed)
│   ├── methods.R             # S3 print/summary/coef/vcov/logLik/predict
│   ├── halton.R              # Halton sequence generation
│   ├── utilities.R           # Formula parsing, start resolution, bounds
│   └── tmb_helpers.R         # MakeADFun wrapper, sdreport wrapper
├── man/
└── tests/
    └── testthat/
        ├── test-fit-famoye.R
        ├── test-fit-copula.R
        └── test-against-rpbnb.R
```

### TMB Template (src/rpbnb_tmb.cpp)

Single template file, branching on `family_code` data variable:

| Code | Family |
|------|--------|
| 0 | Famoye/Sarmanov |
| 1 | Frank copula |
| 2 | Gaussian copula |
| 3 | Clayton (Kimeldorf) copula |
| -1 | Independence (no dependence term) |

**Data passed from R:**

| Data object | Type | Description |
|---|---|---|
| `Y1`, `Y2` | `vector<int>` | Response counts |
| `X1`, `X2` | `matrix<Type>` | Design matrices |
| `rand_idx1`, `rand_idx2` | `vector<int>` | 0-based column indices for random coefs |
| `Z1`, `Z2` | `matrix<Type>` | Halton uniform draws (R_draws × q) |
| `dist1`, `dist2` | `vector<int>` | Distribution codes: 0=normal, 1=lognormal, 2=uniform, 3=triangular |
| `sign1`, `sign2` | `vector<int>` | ±1 lognormal sign constraints |
| `family` | `int` | Family code (-1, 0, 1, 2, 3) |
| `pois1`, `pois2` | `int` | 0/1 flags for exact Poisson (m=0) margin |
| `lamLo`, `lamHi` | `Type` | Frozen Famoye bounds (0/0 sentinel for copula) |

When `q1 + q2 = 0` (no random coefficients), `Z1`/`Z2` are 1×0 matrices and `R_draws = 1`.

**Parameters (optimized):**

| Parameter | Type | Count | Notes |
|---|---|---|---|
| `beta1` | `vector<Type>` | k1 | Mean coefficients eq1 |
| `beta2` | `vector<Type>` | k2 | Mean coefficients eq2 |
| `log_sd1` | `vector<Type>` | q1 | Log-scale for random coefs eq1 |
| `log_sd2` | `vector<Type>` | q2 | Log-scale for random coefs eq2 |
| `log_m1` | `vector<Type>` | 1 | Log dispersion eq1 |
| `log_m2` | `vector<Type>` | 1 | Log dispersion eq2 |
| `z_dep` | `vector<Type>` | 1 | Transformed dependence (lambda or theta) |

**Objective:** `f = -sum_i log(1/R * sum_r exp(LL_ir))`

The simulated log-likelihood averaged over R draws.

**REPORTED quantities (for sdreport):**
- `m1 = exp(log_m1)`, `m2 = exp(log_m2)` — natural-scale dispersions
- For Famoye: `lambda` — dependence parameter (with bounds applied)
- For copula: `theta` or `rho` — native copula parameter
- For all copulas: `tau` — Kendall's tau

## Famoye/Sarmanov Path

### Joint pmf per draw

```
P(y1, y2) = P1(y1) × P2(y2) × [1 + λ × (exp(-y1) - c1) × (exp(-y2) - c2)]
```

### Lambda parameterization

Bounds computed in R at starting values (frozen during optimization):

```
lamLo = max over i of -1/max((1-c1_i)(1-c2_i), c1_i*c2_i)
lamHi = min over i of 1/max(c1_i*(1-c2_i), c2_i*(1-c1_i))
```

Inner logistic map (in TMB):

```
λ = lamLo + (lamHi - lamLo) × (ε + (1 - 2ε) × logistic(z_dep))
```

where ε = 1e-6 keeps λ strictly inside the interval.

### Poisson margin (m=0 branch)

When `pois1 = TRUE`:
- `log_m1` is pinned (fixed parameter in R, not optimized)
- TMB uses `dpois` instead of `dnbinom` for the log-pmf
- `c1 = exp(-d × μ1)` (the exact m→0 limit)

## Copula Path

### Discrete-copula joint pmf per draw

```
P(y1, y2) = C(F1(y1), F2(y2)) - C(F1(y1-1), F2(y2)) - C(F1(y1), F2(y2-1)) + C(F1(y1-1), F2(y2-1))
```

### Copula CDFs implemented in TMB template

**Frank:** `C(u,v;θ) = -1/θ × log(1 + (e^{-θu} - 1)(e^{-θv} - 1)/(e^{-θ} - 1))`
- Native parameter: `θ = z_dep` (unbounded real)

**Gaussian:** `C(u,v;ρ) = Φ_2(Φ^{-1}(u), Φ^{-1}(v); ρ)`
- Native parameter: `ρ = tanh(z_dep)` (ensures |ρ| < 1)
- Bivariate normal CDF via Genz's bvnu algorithm (ported from rpbnb/src/copula_parallel.cpp)

**Clayton:** `C(u,v;θ) = max(u^{-θ} + v^{-θ} - 1, 0)^{-1/θ}`
- Native parameter: `θ = exp(z_dep)` (positive constraint)

### Kendall's tau

```
Frank:      τ = 1 - 4/θ × (1 - D_1(θ))     (Debye function order 1)
Gaussian:   τ = 2/π × asin(ρ)
Clayton:    τ = θ/(θ + 2)
```

All `ADREPORT`ed for standard errors via delta method.

### NB2 CDF in TMB

For copula path, the NB2 CDF `F(y) = pnbinom(y, r, mu)` is computed inside the TMB template using the incomplete beta function relationship (or via a series approximation for small y). TMB includes Rmath via the `pbeta` or direct implementation.

## Random Coefficient Transforms

Per-draw Halton uniform draws `U ~ U(0,1)` are transformed to coefficient deviations:

| Distribution | Base | Coefficient | Dev |
|---|---|---|---|
| Normal (0) | `Z = qnorm(U)` | `β + σZ` | `σZ` |
| Lognormal (1) | `Z = qnorm(U)` | `sign × exp(β + σZ)` | `sign × exp(β + σZ) - β` |
| Uniform (2) | `U` | `β + σ(2U - 1)` | `σ(2U - 1)` |
| Triangular (3) | `T = tri_icdf(U)` | `β + σT` | `σT` |

These are implemented as C++ template functions in the TMB template.

## R-Side API

### fit_rpbnb_tmb()

```r
fit_rpbnb_tmb(formula_1, formula_2, data,
              random_1 = NULL, random_2 = NULL,
              draws = 400, seed = 1234, start = NULL,
              dependence = "famoye",
              control = rpbnb_tmb_control(),
              poisson_1 = FALSE, poisson_2 = FALSE)
```

**`dependence`** accepts:
- `"famoye"` — Famoye/Sarmanov dependence
- `"independence"` — two univariate NB2 models
- A `copula()` object: `copula("frank")`, `copula("normal")`, `copula("kimeldorf")`

**Re-exported `copula()` helper** (same interface as rpbnb) for constructing copula specs.

**`random_1`/`random_2`** — same interface as rpbnb:
- `NULL`: all fixed
- Character vector: all normal random coefficients
- Named list: per-variable distribution specs (normal, lognormal, uniform, triangular)

### rpbnb_tmb_control()

```r
rpbnb_tmb_control(iterlim = 500, reltol = 1e-8,
                  print_level = 0L, n_cores = 1L,
                  optimizer = "nlminb")
```

### Workflow

1. Parse formulas → design matrices
2. Parse random specs → indices + distributions
3. Generate Halton draws (Cranley-Patterson rotation)
4. If Famoye: compute frozen lambda bounds at start
5. If copula: set lamLo=lamHi=0 (sentinel)
6. Build TMB data list + parameter list
7. `MakeADFun(data, parameters, DLL = "rpbnb_tmb")`
8. `nlminb(start, obj$fn, obj$gr)`
9. `sdreport(obj)` → SEs
10. Construct `rpbnb_tmb_fit` S3 object

### S3 Methods

- `print()` — parameter estimates, SEs, logLik
- `summary()` — coefficient table, dispersion, dependence, model summary
- `coef()` — full parameter vector
- `vcov()` — variance-covariance matrix (from sdreport)
- `logLik()`, `AIC()`, `BIC()`
- `predict()` — population mean E[exp(x'beta)]
- `residuals()` — randomized quantile residuals
- `simulate()` — Monte Carlo simulation from fitted model

## Simulation

### simulate_rpbnb_tmb()

Pure R function (no TMB dependency). Mirrors `rpbnb::simulate_rpbnb()`:

```r
simulate_rpbnb_tmb(n, beta1, beta2,
                   random_1 = NULL, random_2 = NULL,
                   dispersion = c(m1 = 0.5, m2 = 0.5),
                   dependence = "famoye",
                   copula = NULL,
                   seed = 1234)
```

For the Famoye case, the joint pmf is the product-form Sarmanov:
```
P(y1, y2) = P1(y1) × P2(y2) × [1 + λ × (exp(-y1) - c1) × (exp(-y2) - c2)]
```

Simulation uses the conditional method: draw y1 from the marginal, then draw y2 from the conditional (using the copula CDF for the copula case, or the Sarmanov conditional for Famoye). Falls back to the bivariate acceptance-rejection or the rpbnb approach of joint pmf enumeration for small counts.

## Verification

| Test | Method |
|---|---|
| **Log-likelihood agreement** | Fit same model in rpbnb and rpbnb_tmb; compare LL within 1e-6 |
| **Parameter agreement** | Compare coef estimates within 1e-4 relative difference |
| **SE agreement** | Compare SEs within 5% (method differences expected) |
| **AD gradient check** | TMB AD gradient vs numDeriv::grad on TMB objective |
| **Truth recovery** | Simulate with known params, fit, check truth in 95% CI |
| **Independence boundary** | ρ=0 / θ=0 / λ=0 → same LL as two marginal fits |
| **Poisson limit** | poisson_1=TRUE → log_m1 = log(POISSON_M) pinned, LL matches dpois |
| **Random coef recovery** | Known random SD, check 95% CI covers truth |

## Standard Errors

Primary method: **TMB `sdreport()`** — uses the generalized delta method on the AD Hessian of the negative log-likelihood. Reports SEs for:
- All free parameters (beta, log_sd, log_m, z_dep)
- Derived parameters: m1, m2, lambda/theta/rho, tau

## Comparison with rpbnb

| Feature | rpbnb | rpbnb_tmb |
|---|---|---|
| Gradient | Hand-coded analytic | TMB AD (automatic) |
| Hessian | Analytic, numeric, or OPG | AD Hessian (via sdreport) |
| Optimizer | maxLik::maxLik (BFGS) | stats::nlminb |
| C++ backend | Rcpp + OpenMP | TMB (CppAD + template) |
| Parallelism | OpenMP (per-draw) | OpenMP (per-draw, via TMB) |
| Poisson margins | Exact m=0 branch | Exact m=0 branch |
| Copula CDFs | R + Rcpp | TMB template |
| SE method | Configurable (analytic/num/OPG) | TMB sdreport (default) |
| Standard errors | Optional | Always computed |

## Files to Create

### Core (required for estimation)
1. `rpbnb_tmb/DESCRIPTION`
2. `rpbnb_tmb/NAMESPACE`
3. `rpbnb_tmb/src/rpbnb_tmb.cpp` — TMB template
4. `rpbnb_tmb/src/Makevars`
5. `rpbnb_tmb/R/fit_rpbnb_tmb.R` — main estimator
6. `rpbnb_tmb/R/utilities.R` — formula parsing, Halton, bounds
7. `rpbnb_tmb/R/tmb_helpers.R` — MakeADFun wrapper, sdreport
8. `rpbnb_tmb/R/rpbnb_tmb-package.R` — package-level docs

### Supporting (needed for completeness)
9. `rpbnb_tmb/R/simulate_rpbnb_tmb.R` — data simulator
10. `rpbnb_tmb/R/methods.R` — S3 methods
11. `rpbnb_tmb/R/halton.R` — Halton sequence (standalone)

### Tests
12. `rpbnb_tmb/tests/testthat/test-fit-famoye.R`
13. `rpbnb_tmb/tests/testthat/test-fit-copula.R`
14. `rpbnb_tmb/tests/testthat/test-against-rpbnb.R`

## Implementation Order

1. **Package skeleton** — DESCRIPTION, NAMESPACE, Makevars, package doc
2. **R utilities** — formula parsing, Halton generation, start resolution
3. **TMB template** — NB2 math, random coef transforms, Famoye likelihood, copula CDFs
4. **Main estimator** — fit_rpbnb_tmb() with nlminb optimization
5. **Standard errors** — sdreport integration
6. **S3 methods** — print, summary, coef, vcov, logLik, predict
7. **Simulator** — simulate_rpbnb_tmb()
8. **Tests** — cross-validation against rpbnb
9. **Poisson margins** — poisson_1/poisson_2 support
10. **Independence path** — dependence = "independence"
