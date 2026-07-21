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
  int k1 = X1.cols(); (void)k1;
  int k2 = X2.cols(); (void)k2;
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

  // ---- Bivariate normal CDF (for Gaussian copula) ----
  auto bvncdf = [](Type h, Type k, Type r) -> Type {
    // P(Z1 <= h, Z2 <= k) for standard bivariate normal with correlation r
    // Uses Gauss-Legendre quadrature on the integral:
    // Phi2(h,k;r) = int_0^{Phi(h)} Phi((k - r*Phi^{-1}(p))/sqrt(1-r^2)) dp
    if (r == Type(0)) return pnorm(h) * pnorm(k);
    if (!CppAD::isfinite(h)) return (h > Type(0)) ? pnorm(k) : Type(0);
    if (!CppAD::isfinite(k)) return (k > Type(0)) ? pnorm(h) : Type(0);
    if (r >= Type(1.0)) {
      if (h < k) return pnorm(h); else return pnorm(k);
    }
    if (r <= Type(-1.0)) {
      Type s = pnorm(h) + pnorm(k) - Type(1.0);
      if (s < Type(0.0)) s = Type(0.0);
      return s;
    }
    Type sig = sqrt(Type(1.0) - r * r);
    if (sig < Type(1e-15)) sig = Type(1e-15);
    // Choose shorter integration path: use min(Phi(h), Phi(k))
    Type ph = pnorm(h);
    Type pk = pnorm(k);
    Type a, oth;
    if (ph <= pk) { a = ph; oth = k; } else { a = pk; oth = h; }
    if (a >= Type(1.0)) a = Type(1.0 - 1e-15);
    if (a <= Type(0.0)) return Type(0);
    // 20-point Gauss-Legendre quadrature on [0, a]
    static const double x20[10] = {0.9931285991850949, 0.9639719272779138,
      0.9122344282513259, 0.8391169718222188, 0.7463319064601508,
      0.6360536807265150, 0.5108670019508271, 0.3737060887154196,
      0.2277858511416451, 0.07652652113349733};
    static const double w20[10] = {0.01761400713915212, 0.04060142980038694,
      0.06267204833410906, 0.08327674157670475, 0.1019301198172404,
      0.1181945319615184, 0.1316886384491766, 0.1420961093183821,
      0.1491729864726037, 0.1527533871307259};
    Type result = 0;
    Type hw = a / Type(2.0);
    Type mid = a / Type(2.0);
    for (int i = 0; i < 10; i++) {
      Type z1 = qnorm(mid + hw * x20[i]);
      Type z2 = qnorm(mid - hw * x20[i]);
      Type arg1 = (oth - r * z1) / sig;
      Type arg2 = (oth - r * z2) / sig;
      result += w20[i] * (pnorm(arg1) + pnorm(arg2));
    }
    result = result * hw;
    return result;
  };

  // Famoye constant d = 1 - exp(-1) (loop-invariant)
  Type famoye_d = Type(1.0) - exp(Type(-1.0));

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
        Type c1, c2;
        if (pois1) {
          c1 = exp(-famoye_d * mu1(i));
        } else {
          c1 = pow(Type(1.0) + famoye_d * m1 * mu1(i), Type(-1.0) / m1);
        }
        if (pois2) {
          c2 = exp(-famoye_d * mu2(i));
        } else {
          c2 = pow(Type(1.0) + famoye_d * m2 * mu2(i), Type(-1.0) / m2);
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
    } else if (family == FAM_FRANK) {
      // Frank copula: C(u,v;theta) = -1/theta * log(1 + (e^{-theta*u}-1)(e^{-theta*v}-1)/(e^{-theta}-1))
      for (int i = 0; i < n; i++) {
        // NB2 CDF corners
        Type a1, a1m, b1, b1m;
        Type r1i = r1, r2i = r2;
        // Eq1 CDF
        if (pois1) {
          a1 = ppois(Y1(i), mu1(i));
          a1m = (Y1(i) > Type(0)) ? ppois(Y1(i) - Type(1), mu1(i)) : Type(0);
        } else {
          Type p1 = r1i / (r1i + mu1(i));
          a1 = pbeta(p1, r1i, Y1(i) + Type(1));
          a1m = (Y1(i) > Type(0)) ? pbeta(p1, r1i, Y1(i)) : Type(0);
        }
        // Eq2 CDF
        if (pois2) {
          b1 = ppois(Y2(i), mu2(i));
          b1m = (Y2(i) > Type(0)) ? ppois(Y2(i) - Type(1), mu2(i)) : Type(0);
        } else {
          Type p2 = r2i / (r2i + mu2(i));
          b1 = pbeta(p2, r2i, Y2(i) + Type(1));
          b1m = (Y2(i) > Type(0)) ? pbeta(p2, r2i, Y2(i)) : Type(0);
        }
        // Frank CDF lambda
        auto frank_cdf = [](Type u, Type v, Type th) -> Type {
          if (th > Type(-1e-10) && th < Type(1e-10)) return u * v;
          Type et = exp(-th);
          Type num = (exp(-th * u) - Type(1)) * (exp(-th * v) - Type(1));
          Type denom = et - Type(1);
          Type arg = Type(1) + num / denom;
          if (arg < Type(1e-300)) arg = Type(1e-300);
          return -log(arg) / th;
        };
        Type C_ab   = frank_cdf(a1,   b1,  theta);
        Type C_amb  = frank_cdf(a1m,  b1,  theta);
        Type C_abm  = frank_cdf(a1,   b1m, theta);
        Type C_ambm = frank_cdf(a1m,  b1m, theta);
        Type p_obs = C_ab - C_amb - C_abm + C_ambm;
        if (p_obs < Type(1e-300)) p_obs = Type(1e-300);
        ll_draw(i) = log(p_obs);
      }
    } else if (family == FAM_GAUSSIAN) {
      // Gaussian copula: C(u,v;rho) = Phi_2(Phi^{-1}(u), Phi^{-1}(v); rho)
      for (int i = 0; i < n; i++) {
        // NB2 CDF corners
        Type a1, a1m, b1, b1m;
        Type r1i = r1, r2i = r2;
        if (pois1) {
          a1 = ppois(Y1(i), mu1(i));
          a1m = (Y1(i) > Type(0)) ? ppois(Y1(i) - Type(1), mu1(i)) : Type(0);
        } else {
          Type p1 = r1i / (r1i + mu1(i));
          a1 = pbeta(p1, r1i, Y1(i) + Type(1));
          a1m = (Y1(i) > Type(0)) ? pbeta(p1, r1i, Y1(i)) : Type(0);
        }
        if (pois2) {
          b1 = ppois(Y2(i), mu2(i));
          b1m = (Y2(i) > Type(0)) ? ppois(Y2(i) - Type(1), mu2(i)) : Type(0);
        } else {
          Type p2 = r2i / (r2i + mu2(i));
          b1 = pbeta(p2, r2i, Y2(i) + Type(1));
          b1m = (Y2(i) > Type(0)) ? pbeta(p2, r2i, Y2(i)) : Type(0);
        }
        // qnorm of each corner, clamped for numerical stability
        auto safe_qnorm = [](Type p) -> Type {
          if (p < Type(1e-15)) p = Type(1e-15);
          if (p > Type(1.0 - 1e-15)) p = Type(1.0 - 1e-15);
          return qnorm(p);
        };
        Type qa  = safe_qnorm(a1);
        Type qam = safe_qnorm(a1m);
        Type qb  = safe_qnorm(b1);
        Type qbm = safe_qnorm(b1m);
        Type C_ab   = bvncdf(qa,   qb,  rho);
        Type C_amb  = bvncdf(qam,  qb,  rho);
        Type C_abm  = bvncdf(qa,   qbm, rho);
        Type C_ambm = bvncdf(qam,  qbm, rho);
        Type p_obs = C_ab - C_amb - C_abm + C_ambm;
        if (p_obs < Type(1e-300)) p_obs = Type(1e-300);
        ll_draw(i) = log(p_obs);
      }
    } else if (family == FAM_CLAYTON) {
      // Clayton copula: C(u,v;th) = max(u^{-th} + v^{-th} - 1, 0)^{-1/th}
      for (int i = 0; i < n; i++) {
        // NB2 CDF corners
        Type a1, a1m, b1, b1m;
        Type r1i = r1, r2i = r2;
        if (pois1) {
          a1 = ppois(Y1(i), mu1(i));
          a1m = (Y1(i) > Type(0)) ? ppois(Y1(i) - Type(1), mu1(i)) : Type(0);
        } else {
          Type p1 = r1i / (r1i + mu1(i));
          a1 = pbeta(p1, r1i, Y1(i) + Type(1));
          a1m = (Y1(i) > Type(0)) ? pbeta(p1, r1i, Y1(i)) : Type(0);
        }
        if (pois2) {
          b1 = ppois(Y2(i), mu2(i));
          b1m = (Y2(i) > Type(0)) ? ppois(Y2(i) - Type(1), mu2(i)) : Type(0);
        } else {
          Type p2 = r2i / (r2i + mu2(i));
          b1 = pbeta(p2, r2i, Y2(i) + Type(1));
          b1m = (Y2(i) > Type(0)) ? pbeta(p2, r2i, Y2(i)) : Type(0);
        }
        // Clayton CDF lambda
        auto clayton_cdf = [](Type u, Type v, Type th) -> Type {
          if (th < Type(1e-10)) return u * v;
          Type inner = pow(u, -th) + pow(v, -th) - Type(1);
          if (inner < Type(0)) inner = Type(0);
          return pow(inner, Type(-1.0) / th);
        };
        Type C_ab   = clayton_cdf(a1,   b1,  theta);
        Type C_amb  = clayton_cdf(a1m,  b1,  theta);
        Type C_abm  = clayton_cdf(a1,   b1m, theta);
        Type C_ambm = clayton_cdf(a1m,  b1m, theta);
        Type p_obs = C_ab - C_amb - C_abm + C_ambm;
        if (p_obs < Type(1e-300)) p_obs = Type(1e-300);
        ll_draw(i) = log(p_obs);
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

  // ---- Kendall's tau for copula families ----
  if (family == FAM_FRANK || family == FAM_GAUSSIAN || family == FAM_CLAYTON) {
    Type tau;
    if (family == FAM_GAUSSIAN) {
      tau = Type(2.0) / M_PI * asin(rho);
    } else if (family == FAM_CLAYTON) {
      tau = theta / (theta + Type(2.0));
    } else {  // Frank
      // Frank tau: 1 - 4/th * (1 - D1(th)) where D1 is Debye function order 1
      if (fabs(theta) < Type(1e-10)) {
        tau = Type(0);
      } else {
        // 20-point Gauss-Legendre quadrature on [0, theta]
        static const double x20[10] = {0.9931285991850949, 0.9639719272779138,
          0.9122344282513259, 0.8391169718222188, 0.7463319064601508,
          0.6360536807265150, 0.5108670019508271, 0.3737060887154196,
          0.2277858511416451, 0.07652652113349733};
        static const double w20[10] = {0.01761400713915212, 0.04060142980038694,
          0.06267204833410906, 0.08327674157670475, 0.1019301198172404,
          0.1181945319615184, 0.1316886384491766, 0.1420961093183821,
          0.1491729864726037, 0.1527533871307259};
        Type D1 = 0;
        for (int qq = 0; qq < 10; qq++) {
          Type t = theta * Type(0.5) * (Type(1.0) + Type(x20[qq]));
          Type f = t / (exp(t) - Type(1.0));
          if (!CppAD::isfinite(f)) f = Type(0);
          D1 += Type(w20[qq]) * f;
          t = theta * Type(0.5) * (Type(1.0) - Type(x20[qq]));
          f = t / (exp(t) - Type(1.0));
          if (!CppAD::isfinite(f)) f = Type(0);
          D1 += Type(w20[qq]) * f;
        }
        D1 = D1 * Type(0.5);
        tau = Type(1.0) - Type(4.0) / theta * (Type(1.0) - D1);
      }
    }
    ADREPORT(tau);
  }

  return nll;
}
