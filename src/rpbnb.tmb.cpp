// src/rpbnb.tmb.cpp — main TMB template for rpbnb.tmb
// Override TMB_LIB_INIT to avoid '.' in C++ identifier (package name has a dot)
// Command-line -DTMB_LIB_INIT=R_init_rpbnb.tmb is problematic; we override it here
#undef TMB_LIB_INIT
#define TMB_LIB_INIT R_init_rpbnb_tmb
#include <TMB.hpp>
#include <cmath>

// Family codes
#define FAM_INDEP     -1
#define FAM_FAMOYE     0
#define FAM_FRANK      1
#define FAM_GAUSSIAN   2
#define FAM_CLAYTON    3

// Distribution codes
#define DIST_NORMAL     0
#define DIST_LOGNORMAL  1
#define DIST_UNIFORM    2
#define DIST_TRIANGULAR 3

template<class Type>
Type objective_function<Type>::operator() () {

  // ---- Data ----
  DATA_VECTOR(Y1);
  DATA_VECTOR(Y2);
  DATA_MATRIX(X1);
  DATA_MATRIX(X2);
  DATA_IVECTOR(rand_idx1);
  DATA_IVECTOR(rand_idx2);
  DATA_MATRIX(Z1);
  DATA_MATRIX(Z2);
  DATA_IVECTOR(dist1);
  DATA_IVECTOR(dist2);
  DATA_IVECTOR(sign1);
  DATA_IVECTOR(sign2);
  DATA_INTEGER(family);
  DATA_INTEGER(pois1);
  DATA_INTEGER(pois2);
  DATA_SCALAR(lamLo);
  DATA_SCALAR(lamHi);

  // ---- Parameters ----
  PARAMETER_VECTOR(beta1);
  PARAMETER_VECTOR(beta2);
  PARAMETER_VECTOR(log_sd1);
  PARAMETER_VECTOR(log_sd2);
  PARAMETER(log_m1);
  PARAMETER(log_m2);
  PARAMETER(z_dep);

  // ---- Dimensions ----
  int n = Y1.size();
  int k1 = X1.cols();
  int k2 = X2.cols();
  int q1 = rand_idx1.size();
  int q2 = rand_idx2.size();
  int R = (q1 + q2 > 0) ? Z1.rows() : 1;

  // ---- Natural-scale parameters ----
  Type m1 = exp(log_m1);
  Type m2 = exp(log_m2);
  Type r1 = 1.0 / m1;
  Type r2 = 1.0 / m2;
  vector<Type> sd1 = exp(log_sd1);
  vector<Type> sd2 = exp(log_sd2);

  // Linear predictors (fixed part)
  vector<Type> xb1 = X1 * beta1;
  vector<Type> xb2 = X2 * beta2;

  // Dependence transform (Famoye: logistic map; Copula: identity/tanh/exp)
  Type eps = 1e-6;
  Type lam, theta, rho;
  if (family == FAM_FAMOYE) {
    Type sig = invlogit(z_dep);  // logistic(0,1) = 1/(1+exp(-x))
    lam = lamLo + (lamHi - lamLo) * (eps + (1.0 - 2.0 * eps) * sig);
  } else if (family == FAM_FRANK) {
    theta = z_dep;
  } else if (family == FAM_GAUSSIAN) {
    rho = tanh(z_dep);
  } else if (family == FAM_CLAYTON) {
    theta = exp(z_dep);
  }

  // ---- Random coefficient transforms (per draw) ----
  // Returns the deviation (eta + dev = eta with random perturbation)
  // Using the TMB Type so AD flows through everything

  // Helper functions inside operator()():
  auto u_to_base = [](Type u, int dist_code) -> Type {
    if (dist_code == DIST_NORMAL || dist_code == DIST_LOGNORMAL)
      return qnorm(u);  // TMB's qnorm for Type
    // DIST_UNIFORM: base = u (identity)
    // DIST_TRIANGULAR: symmetric triangular on [-1, 1]
    if (dist_code == DIST_TRIANGULAR) {
      // tri_icdf: if u < 0.5 then -1+sqrt(2u) else 1-sqrt(2(1-u))
      Type two_u = 2.0 * u;
      return CppAD::CondExpLt(u, Type(0.5), Type(-1.0) + sqrt(two_u),
                              Type(1.0) - sqrt(2.0 * (1.0 - u)));
    }
    return u;  // uniform
  };

  auto compute_dev = [](Type b, Type s, Type base, int dist_code, int sign_code) -> Type {
    if (dist_code == DIST_NORMAL)
      return s * base;
    if (dist_code == DIST_LOGNORMAL)
      return Type(sign_code) * exp(b + s * base) - b;
    if (dist_code == DIST_UNIFORM)
      return s * (2.0 * base - 1.0);
    // triangular
    return s * base;
  };

  // LL matrix: n x R, per-draw per-obs log-likelihood
  // Use parallel_accumulator for OpenMP
  parallel_accumulator<Type> nll_accum(this);
  nll_accum = 0;

  // Pre-compute exp(-y) for Famoye path
  vector<Type> ey1(n), ey2(n);
  if (family == FAM_FAMOYE || family == FAM_INDEP) {
    for (int i = 0; i < n; i++) {
      ey1(i) = exp(-Y1(i));
      ey2(i) = exp(-Y2(i));
    }
  }

  // ---- Per-draw loop (MSL: average likelihood over draws) ----
  const Type logR = log(Type(R));

  // For the independence case: single draw, no dependence term
  // For random coefficients: average over draws

  // Storage for log-lik per draw per obs
  vector<Type> ll_draw(n);
  matrix<Type> ll_all(n, R);

  for (int r = 0; r < R; r++) {
    // Compute per-draw eta = xb + XR * dev
    vector<Type> eta1 = xb1;
    vector<Type> eta2 = xb2;
    if (q1 > 0) {
      for (int j = 0; j < q1; j++) {
        int col = rand_idx1(j);  // 0-based column in X1
        Type u = Z1(r, j);
        Type base = u_to_base(u, dist1(j));
        Type dev = compute_dev(beta1(col), sd1(j), base, dist1(j), sign1(j));
        for (int i = 0; i < n; i++) eta1(i) += X1(i, col) * dev;
      }
    }
    if (q2 > 0) {
      for (int j = 0; j < q2; j++) {
        int col = rand_idx2(j);
        Type u = Z2(r, j);
        Type base = u_to_base(u, dist2(j));
        Type dev = compute_dev(beta2(col), sd2(j), base, dist2(j), sign2(j));
        for (int i = 0; i < n; i++) eta2(i) += X2(i, col) * dev;
      }
    }

    // Per-draw means
    vector<Type> mu1(n), mu2(n);
    for (int i = 0; i < n; i++) {
      mu1(i) = exp(eta1(i));
      if (mu1(i) > Type(1e15)) mu1(i) = Type(1e15);
      if (mu1(i) < Type(1e-300)) mu1(i) = Type(1e-300);
      mu2(i) = exp(eta2(i));
      if (mu2(i) > Type(1e15)) mu2(i) = Type(1e15);
      if (mu2(i) < Type(1e-300)) mu2(i) = Type(1e-300);
    }

    // Per-observation log-likelihood for this draw
    if (family == FAM_FAMOYE) {
      // Famoye: P1 * P2 * [1 + lam * (e^-y - c1) * (e^-y - c2)]
      // Use the exact Poisson branch when pois1/pois2 is set
      for (int i = 0; i < n; i++) {
        // NB2 log-pmf (or Poisson if pois flag set)
        Type lnb1, lnb2;
        if (pois1) {
          lnb1 = dpois(Y1(i), mu1(i), true);
        } else {
          lnb1 = dnbinom2(Y1(i), mu1(i), mu1(i) + m1 * mu1(i) * mu1(i), true);
        }
        if (pois2) {
          lnb2 = dpois(Y2(i), mu2(i), true);
        } else {
          lnb2 = dnbinom2(Y2(i), mu2(i), mu2(i) + m2 * mu2(i) * mu2(i), true);
        }

        // c_val = (1 + d * m * mu)^(-1/m) or exp(-d * mu) for Poisson
        Type d = Type(1.0) - exp(Type(-1.0));
        Type c1, c2;
        if (pois1) {
          c1 = exp(-d * mu1(i));
        } else {
          c1 = pow(Type(1.0) + d * m1 * mu1(i), Type(-1.0) / m1);
        }
        if (pois2) {
          c2 = exp(-d * mu2(i));
        } else {
          c2 = pow(Type(1.0) + d * m2 * mu2(i), Type(-2.0) / m2);
        }

        Type dep = Type(1.0) + lam * (ey1(i) - c1) * (ey2(i) - c2);
        if (dep < Type(1e-300)) dep = Type(1e-300);

        ll_draw(i) = lnb1 + lnb2 + log(dep);
      }
    } else if (family == FAM_INDEP) {
      // Independence: just product of two NB2 (no dependence term)
      for (int i = 0; i < n; i++) {
        Type lnb1 = pois1 ? dpois(Y1(i), mu1(i), true)
                          : dnbinom2(Y1(i), mu1(i), mu1(i) + m1 * mu1(i) * mu1(i), true);
        Type lnb2 = pois2 ? dpois(Y2(i), mu2(i), true)
                          : dnbinom2(Y2(i), mu2(i), mu2(i) + m2 * mu2(i) * mu2(i), true);
        ll_draw(i) = lnb1 + lnb2;
      }
    }

    // Store per-draw LL
    for (int i = 0; i < n; i++) ll_all(i, r) = ll_draw(i);
  }

  // ---- Row-wise log-sum-exp over draws ----
  vector<Type> log_contrib(n);
  for (int i = 0; i < n; i++) {
    Type mx = ll_all(i, 0);
    for (int r = 1; r < R; r++) if (ll_all(i, r) > mx) mx = ll_all(i, r);
    Type s = 0;
    for (int r = 0; r < R; r++) s += exp(ll_all(i, r) - mx);
    log_contrib(i) = mx + log(s) - logR;
  }

  Type nll = -sum(log_contrib);

  // ---- REPORT derived parameters for sdreport ----
  ADREPORT(m1);
  ADREPORT(m2);
  if (family == FAM_FAMOYE) {
    ADREPORT(lam);
  } else if (family == FAM_FRANK) {
    ADREPORT(theta);
  } else if (family == FAM_GAUSSIAN) {
    ADREPORT(rho);
  } else if (family == FAM_CLAYTON) {
    ADREPORT(theta);
  }

  return nll;
}
