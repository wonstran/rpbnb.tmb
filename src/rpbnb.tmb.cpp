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

// Ceiling of the bounded Frank link.  exp(-theta * u) must stay finite, so
// theta is squashed into (-FRANK_THETA_MAX, FRANK_THETA_MAX).  This caps
// attainable Frank dependence at Kendall's tau of about 0.891; R reports a
// boundary warning when an estimate is pinned there.  Must match
// FRANK_THETA_MAX in R/utilities.R.
#define FRANK_THETA_MAX 35.0

static const double gaussian_x20[10] = {
  0.9931285991850949, 0.9639719272779138, 0.9122344282513259,
  0.8391169718222188, 0.7463319064601508, 0.6360536807265150,
  0.5108670019508271, 0.3737060887154196, 0.2277858511416451,
  0.07652652113349733
};
static const double gaussian_w20[10] = {
  0.01761400713915212, 0.04060142980038694, 0.06267204833410906,
  0.08327674157670475, 0.1019301198172404, 0.1181945319615184,
  0.1316886384491766, 0.1420961093183821, 0.1491729864726037,
  0.1527533871307259
};

// CppAD in the supported TMB toolchain does not overload std::expm1/log1p.
// These series-backed equivalents retain precision and AD derivatives near 0.
template<class Type>
Type stable_expm1(Type x) {
  Type x2 = x * x;
  Type series = x + x2 / Type(2) + x2 * x / Type(6) +
    x2 * x2 / Type(24);
  Type regular = exp(x) - Type(1);
  return CppAD::CondExpLt(fabs(x), Type(1e-4), series, regular);
}

template<class Type>
Type stable_log1p(Type x) {
  Type x2 = x * x;
  Type series = x - x2 / Type(2) + x2 * x / Type(3) -
    x2 * x2 / Type(4);
  Type regular = log(Type(1) + x);
  return CppAD::CondExpLt(fabs(x), Type(1e-4), series, regular);
}

// Tape-compressed 20-point bivariate-normal CDF quadrature.
template<class Type>
Type gaussian_bvn_quadrature(Type h, Type k, Type r) {
  Type sig2 = Type(1.0) - r * r;
  sig2 = CppAD::CondExpLt(sig2, Type(1e-12), Type(1e-12), sig2);
  Type sig = sqrt(sig2);
  Type ph = pnorm(h);
  Type pk = pnorm(k);
  Type a = CppAD::CondExpLe(ph, pk, ph, pk);
  Type oth = CppAD::CondExpLe(ph, pk, k, h);
  oth = CppAD::CondExpGt(oth, Type(10.0), Type(10.0), oth);
  oth = CppAD::CondExpLt(oth, Type(-10.0), Type(-10.0), oth);
  a = CppAD::CondExpGe(a, Type(1.0), Type(1.0 - 1e-10), a);
  a = CppAD::CondExpLe(a, Type(1e-10), Type(1e-10), a);

  Type result = Type(0);
  Type hw = a / Type(2.0);
  Type mid = a / Type(2.0);
  for (int i = 0; i < 10; i++) {
    Type z1 = qnorm(mid + hw * Type(gaussian_x20[i]));
    Type z2 = qnorm(mid - hw * Type(gaussian_x20[i]));
    Type arg1 = (oth - r * z1) / sig;
    Type arg2 = (oth - r * z2) / sig;
    result += Type(gaussian_w20[i]) * (pnorm(arg1) + pnorm(arg2));
  }

  return result * hw;
}

// Differentiate the quadrature itself so the atomic gradient remains exactly
// consistent with its value, including its conditional integration limits.
template<class Type>
vector<Type> gaussian_bvn_quadrature_gradient(Type h, Type k, Type r) {
  Type raw_sig2 = Type(1.0) - r * r;
  Type sig2 = CppAD::CondExpLt(
    raw_sig2, Type(1e-12), Type(1e-12), raw_sig2
  );
  Type dsig2_dr = CppAD::CondExpLt(
    raw_sig2, Type(1e-12), Type(0), Type(-2) * r
  );
  Type sig = sqrt(sig2);
  Type dsig_dr = dsig2_dr / (Type(2) * sig);

  Type ph = pnorm(h);
  Type pk = pnorm(k);
  Type raw_a = CppAD::CondExpLe(ph, pk, ph, pk);
  Type raw_oth = CppAD::CondExpLe(ph, pk, k, h);
  Type a = CppAD::CondExpGe(
    raw_a, Type(1.0), Type(1.0 - 1e-10), raw_a
  );
  a = CppAD::CondExpLe(a, Type(1e-10), Type(1e-10), a);
  Type oth = CppAD::CondExpGt(
    raw_oth, Type(10.0), Type(10.0), raw_oth
  );
  oth = CppAD::CondExpLt(oth, Type(-10.0), Type(-10.0), oth);

  Type a_active = CppAD::CondExpGe(
    raw_a, Type(1.0), Type(0),
    CppAD::CondExpLe(raw_a, Type(1e-10), Type(0), Type(1))
  );
  Type oth_active = CppAD::CondExpGt(
    raw_oth, Type(10.0), Type(0),
    CppAD::CondExpLt(raw_oth, Type(-10.0), Type(0), Type(1))
  );
  Type da_dh = a_active * CppAD::CondExpLe(
    ph, pk, dnorm(h, Type(0), Type(1), false), Type(0)
  );
  Type da_dk = a_active * CppAD::CondExpLe(
    ph, pk, Type(0), dnorm(k, Type(0), Type(1), false)
  );
  Type doth_dh = oth_active * CppAD::CondExpLe(
    ph, pk, Type(0), Type(1)
  );
  Type doth_dk = oth_active * CppAD::CondExpLe(
    ph, pk, Type(1), Type(0)
  );

  Type sum_f = Type(0);
  Type sum_da = Type(0);
  Type sum_doth = Type(0);
  Type sum_dr = Type(0);
  for (int i = 0; i < 10; i++) {
    for (int side = -1; side <= 1; side += 2) {
      Type c = (Type(1) + Type(side) * Type(gaussian_x20[i])) /
        Type(2);
      Type z = qnorm(a * c);
      Type numerator = oth - r * z;
      Type arg = numerator / sig;
      Type phi_arg = dnorm(arg, Type(0), Type(1), false);
      Type phi_z = dnorm(z, Type(0), Type(1), false);
      Type weight = Type(gaussian_w20[i]);

      sum_f += weight * pnorm(arg);
      sum_da += weight * phi_arg * (-r / sig) * c / phi_z;
      sum_doth += weight * phi_arg / sig;
      sum_dr += weight * phi_arg * (
        -z / sig - numerator * dsig_dr / sig2
      );
    }
  }

  Type dvalue_da = sum_f / Type(2) + a * sum_da / Type(2);
  Type dvalue_doth = a * sum_doth / Type(2);
  Type dvalue_dr = a * sum_dr / Type(2);
  vector<Type> gradient(3);
  gradient(0) = dvalue_da * da_dh + dvalue_doth * doth_dh;
  gradient(1) = dvalue_da * da_dk + dvalue_doth * doth_dk;
  gradient(2) = dvalue_dr;
  return gradient;
}

// A known-derivative primitive is thread safe and compresses the quadrature
// to one tape operation. TMB differentiates its reverse rule for Hessians.
TMB_ATOMIC_VECTOR_FUNCTION(gaussian_bvn_atomic,
  1,
  ty[0] = gaussian_bvn_quadrature(tx[0], tx[1], tx[2]);
  ,
  vector<Type> gradient = gaussian_bvn_quadrature_gradient(
    tx[0], tx[1], tx[2]
  );
  px[0] = py[0] * gradient(0);
  px[1] = py[0] * gradient(1);
  px[2] = py[0] * gradient(2);
)

template<class Type>
Type gaussian_bvn_cdf(Type h, Type k, Type r) {
  CppAD::vector<Type> input(3);
  input[0] = h;
  input[1] = k;
  input[2] = r;
  return gaussian_bvn_atomic(input)[0];
}

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
  int openmp_compiled = 0;
#ifdef _OPENMP
  openmp_compiled = 1;
#endif
  REPORT(openmp_compiled);

  // ---- Natural-scale parameters ----
  auto clamp_ad = [](Type x, Type lo, Type hi) -> Type {
    x = CppAD::CondExpLt(x, lo, lo, x);
    return CppAD::CondExpGt(x, hi, hi, x);
  };
  Type log_m1_c = clamp_ad(log_m1, Type(-20.0), Type(20.0));
  Type log_m2_c = clamp_ad(log_m2, Type(-20.0), Type(20.0));
  Type m1 = exp(log_m1_c);
  Type m2 = exp(log_m2_c);
  Type r1 = 1.0 / m1;
  Type r2 = 1.0 / m2;

  // dnbinom2(y, mu, mu + m*mu*mu) evaluates log(var - mu).  The variance
  // increment m*mu*mu is lost to rounding once it falls below ulp(mu), i.e.
  // once log(m) + log(mu) drops under about -36.04 = log(2^-52).  Clamping
  // the linear predictor alone cannot enforce that, because m is estimated:
  // the floor has to move with m.  Keep -35 as the ceiling on the floor so
  // that over-dispersed fits are unaffected.
  auto nb2_eta_floor = [&](Type log_m_clamped) -> Type {
    Type negative_log_m = CppAD::CondExpLt(
      log_m_clamped, Type(0), log_m_clamped, Type(0)
    );
    return Type(-35.0) - negative_log_m;
  };
  Type eta_floor1 = pois1 ? Type(-35.0) : nb2_eta_floor(log_m1_c);
  Type eta_floor2 = pois2 ? Type(-35.0) : nb2_eta_floor(log_m2_c);
  vector<Type> sd1(log_sd1.size()), sd2(log_sd2.size());
  for (int j = 0; j < log_sd1.size(); j++)
    sd1(j) = exp(clamp_ad(log_sd1(j), Type(-20.0), Type(20.0)));
  for (int j = 0; j < log_sd2.size(); j++)
    sd2(j) = exp(clamp_ad(log_sd2(j), Type(-20.0), Type(20.0)));

  // Linear predictors (fixed part)
  vector<Type> xb1 = X1 * beta1;
  vector<Type> xb2 = X2 * beta2;

  // Dependence transform (Famoye: logistic map; Copula: identity/tanh/exp)
  Type eps = 1e-6;
  Type lam = Type(0), theta = Type(0), rho = Type(0);
  if (family == FAM_FAMOYE) {
    Type sig = invlogit(z_dep);  // logistic(0,1) = 1/(1+exp(-x))
    lam = lamLo + (lamHi - lamLo) * (eps + (1.0 - 2.0 * eps) * sig);
  } else if (family == FAM_FRANK) {
    // A smooth bounded link prevents exponential overflow while retaining
    // theta = 0 and unit derivative at independence.
    theta = Type(FRANK_THETA_MAX) * tanh(z_dep / Type(FRANK_THETA_MAX));
  } else if (family == FAM_GAUSSIAN) {
    rho = tanh(z_dep);
  } else if (family == FAM_CLAYTON) {
    theta = exp(clamp_ad(z_dep, Type(-20.0), Type(20.0)));
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

  // Precompute parameter-dependent deviations once per simulation draw.  The
  // matrices are read-only while observation contributions are accumulated.
  matrix<Type> dev1(R, q1), dev2(R, q2);
  for (int r = 0; r < R; r++) {
    for (int j = 0; j < q1; j++) {
      Type base = u_to_base(Z1(r, j), dist1(j));
      int col = rand_idx1(j);
      dev1(r, j) = compute_dev(beta1(col), sd1(j), base,
                               dist1(j), sign1(j));
    }
    for (int j = 0; j < q2; j++) {
      Type base = u_to_base(Z2(r, j), dist2(j));
      int col = rand_idx2(j);
      dev2(r, j) = compute_dev(beta2(col), sd2(j), base,
                               dist2(j), sign2(j));
    }
  }

  // TMB partitions these independent observation contributions across its
  // configured OpenMP regions while keeping each draw reduction local.
  parallel_accumulator<Type> nll(this);
  const Type logR = log(Type(R));

  for (int i = 0; i < n; i++) {
    vector<Type> log_draw(R);

    for (int r = 0; r < R; r++) {
      Type eta1 = xb1(i);
      Type eta2 = xb2(i);
      for (int j = 0; j < q1; j++) {
        eta1 += X1(i, rand_idx1(j)) * dev1(r, j);
      }
      for (int j = 0; j < q2; j++) {
        eta2 += X2(i, rand_idx2(j)) * dev2(r, j);
      }

      Type mu1 = exp(clamp_ad(eta1, eta_floor1,
                              Type(34.538776394910684)));
      Type mu2 = exp(clamp_ad(eta2, eta_floor2,
                              Type(34.538776394910684)));

      if (family == FAM_FAMOYE) {
        Type lnb1 = pois1 ? dpois(Y1(i), mu1, true)
                          : dnbinom2(Y1(i), mu1,
                                     mu1 + m1 * mu1 * mu1, true);
        Type lnb2 = pois2 ? dpois(Y2(i), mu2, true)
                          : dnbinom2(Y2(i), mu2,
                                     mu2 + m2 * mu2 * mu2, true);

        Type c1 = pois1
          ? exp(-famoye_d * mu1)
          : exp(-stable_log1p(famoye_d * m1 * mu1) / m1);
        Type c2 = pois2
          ? exp(-famoye_d * mu2)
          : exp(-stable_log1p(famoye_d * m2 * mu2) / m2);
        Type dep = Type(1.0) + lam * (ey1(i) - c1) * (ey2(i) - c2);
        // The Sarmanov factor can go non-positive because lamLo/lamHi are
        // frozen at the starting values rather than recomputed at the current
        // mu.  Penalising makes that visible in the objective instead of
        // hiding it behind a probability clamp.
        //
        // This is a value-only barrier, NOT a constraint: CondExpLe is a step,
        // so the penalty term contributes exactly zero gradient on both sides.
        // A gradient-driven optimizer is repelled only by the function value,
        // and cannot be steered out of the invalid region by the score.  The
        // real fix is parameter-dependent bounds evaluated here on the tape;
        // until then, treat a fit that lands on the penalty as invalid rather
        // than as a converged optimum.
        Type invalid_dep = CppAD::CondExpLe(
          dep, Type(0), Type(1), Type(0)
        );
        Type safe_dep = CppAD::CondExpLe(
          dep, Type(0), Type(1e-300), dep
        );
        log_draw(r) = lnb1 + lnb2 + log(safe_dep) -
          invalid_dep * Type(1e10);
      } else if (family == FAM_INDEP) {
        Type lnb1 = pois1 ? dpois(Y1(i), mu1, true)
                          : dnbinom2(Y1(i), mu1,
                                     mu1 + m1 * mu1 * mu1, true);
        Type lnb2 = pois2 ? dpois(Y2(i), mu2, true)
                          : dnbinom2(Y2(i), mu2,
                                     mu2 + m2 * mu2 * mu2, true);
        log_draw(r) = lnb1 + lnb2;
      } else {
        Type a1, a1m, b1, b1m;
        if (pois1) {
          a1 = ppois(Y1(i), mu1);
          a1m = (Y1(i) > Type(0))
            ? ppois(Y1(i) - Type(1), mu1) : Type(0);
        } else {
          Type p1 = r1 / (r1 + mu1);
          a1 = pbeta(p1, r1, Y1(i) + Type(1));
          a1m = (Y1(i) > Type(0)) ? pbeta(p1, r1, Y1(i)) : Type(0);
        }
        if (pois2) {
          b1 = ppois(Y2(i), mu2);
          b1m = (Y2(i) > Type(0))
            ? ppois(Y2(i) - Type(1), mu2) : Type(0);
        } else {
          Type p2 = r2 / (r2 + mu2);
          b1 = pbeta(p2, r2, Y2(i) + Type(1));
          b1m = (Y2(i) > Type(0)) ? pbeta(p2, r2, Y2(i)) : Type(0);
        }

        Type C_ab, C_amb, C_abm, C_ambm;
        if (family == FAM_FRANK) {
          auto frank_cdf = [](Type u, Type v, Type th) -> Type {
            Type signed_eps = CppAD::CondExpGe(
              th, Type(0), Type(1e-5), Type(-1e-5)
            );
            Type safe_th = CppAD::CondExpLt(
              fabs(th), Type(1e-5), signed_eps, th
            );
            Type ratio = stable_expm1(-safe_th * u) *
                         stable_expm1(-safe_th * v) /
                         stable_expm1(-safe_th);
            ratio = CppAD::CondExpLt(ratio, Type(-1.0 + 1e-15),
                                     Type(-1.0 + 1e-15), ratio);
            Type regular = -stable_log1p(ratio) / safe_th;
            Type near_independence = u * v +
              th * u * v * (Type(1) - u) * (Type(1) - v) / Type(2);
            return CppAD::CondExpLt(
              fabs(th), Type(1e-5), near_independence, regular
            );
          };
          C_ab   = frank_cdf(a1,   b1,  theta);
          C_amb  = frank_cdf(a1m,  b1,  theta);
          C_abm  = frank_cdf(a1,   b1m, theta);
          C_ambm = frank_cdf(a1m,  b1m, theta);
        } else if (family == FAM_GAUSSIAN) {
          auto safe_qnorm = [](Type p) -> Type {
            p = CppAD::CondExpLt(p, Type(1e-15), Type(1e-15), p);
            p = CppAD::CondExpGt(p, Type(1.0 - 1e-15),
                                 Type(1.0 - 1e-15), p);
            return qnorm(p);
          };
          Type qa = safe_qnorm(a1);
          Type qam = safe_qnorm(a1m);
          Type qb = safe_qnorm(b1);
          Type qbm = safe_qnorm(b1m);
          C_ab   = gaussian_bvn_cdf(qa,   qb,  rho);
          C_amb  = gaussian_bvn_cdf(qam,  qb,  rho);
          C_abm  = gaussian_bvn_cdf(qa,   qbm, rho);
          C_ambm = gaussian_bvn_cdf(qam,  qbm, rho);
        } else {
          auto clayton_cdf = [](Type u, Type v, Type th) -> Type {
            // These exact-zero branches depend only on observed counts, so
            // they remain valid on a non-retaped objective.
            if (u == Type(0.0) || v == Type(0.0)) return Type(0.0);
            u = CppAD::CondExpLt(u, Type(1e-15), Type(1e-15), u);
            v = CppAD::CondExpLt(v, Type(1e-15), Type(1e-15), v);
            Type inner = pow(u, -th) + pow(v, -th) - Type(1);
            inner = CppAD::CondExpLt(inner, Type(1e-300),
                                     Type(1e-300), inner);
            return pow(inner, Type(-1.0) / th);
          };
          C_ab   = clayton_cdf(a1,   b1,  theta);
          C_amb  = clayton_cdf(a1m,  b1,  theta);
          C_abm  = clayton_cdf(a1,   b1m, theta);
          C_ambm = clayton_cdf(a1m,  b1m, theta);
        }

        Type p_obs = C_ab - C_amb - C_abm + C_ambm;
        p_obs = CppAD::CondExpLt(p_obs, Type(1e-300),
                                 Type(1e-300), p_obs);
        log_draw(r) = log(p_obs);
      }
    }

    Type max_log = log_draw(0);
    for (int r = 1; r < R; r++) {
      max_log = CppAD::CondExpGt(log_draw(r), max_log,
                                 log_draw(r), max_log);
    }
    Type scaled_sum = Type(0);
    for (int r = 0; r < R; r++) {
      scaled_sum += exp(log_draw(r) - max_log);
    }
    Type log_contribution = max_log + log(scaled_sum) - logR;
    nll -= log_contribution;
  }

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
      tau = Type(2.0) / Type(3.14159265358979323846) * asin(rho);
    } else if (family == FAM_CLAYTON) {
      tau = theta / (theta + Type(2.0));
    } else {  // Frank
      // Frank tau: 1 - 4/th * (1 - D1(th)) where D1 is Debye function order 1
      {
        // 20-point Gauss-Legendre quadrature on [0, theta], with a
        // small-theta series that keeps the value and derivatives smooth.
        static const double x20[10] = {0.9931285991850949, 0.9639719272779138,
          0.9122344282513259, 0.8391169718222188, 0.7463319064601508,
          0.6360536807265150, 0.5108670019508271, 0.3737060887154196,
          0.2277858511416451, 0.07652652113349733};
        static const double w20[10] = {0.01761400713915212, 0.04060142980038694,
          0.06267204833410906, 0.08327674157670475, 0.1019301198172404,
          0.1181945319615184, 0.1316886384491766, 0.1420961093183821,
          0.1491729864726037, 0.1527533871307259};
        Type signed_eps = CppAD::CondExpGe(
          theta, Type(0), Type(1e-4), Type(-1e-4)
        );
        Type quadrature_theta = CppAD::CondExpLt(
          fabs(theta), Type(1e-4), signed_eps, theta
        );
        Type D1 = 0;
        for (int qq = 0; qq < 10; qq++) {
          Type t = quadrature_theta * Type(0.5) *
            (Type(1.0) + Type(x20[qq]));
          Type f = t / stable_expm1(t);
          D1 += Type(w20[qq]) * f;
          t = quadrature_theta * Type(0.5) *
            (Type(1.0) - Type(x20[qq]));
          f = t / stable_expm1(t);
          D1 += Type(w20[qq]) * f;
        }
        D1 = D1 * Type(0.5);
        Type quadrature_tau = Type(1.0) -
          Type(4.0) / quadrature_theta * (Type(1.0) - D1);
        Type series_tau = theta / Type(9.0) -
          theta * theta * theta / Type(900.0);
        tau = CppAD::CondExpLt(
          fabs(theta), Type(1e-4), series_tau, quadrature_tau
        );
      }
    }
    ADREPORT(tau);
  }

  return nll;
}
