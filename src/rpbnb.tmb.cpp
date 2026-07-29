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

// 16-point Gauss-Legendre, positive half of the symmetric node set. Applied
// per panel by gaussian_cell_prob(), which splits its interval into five.
static const double gauss_x16[8] = {
  0.09501250983763769, 0.28160355077925908, 0.45801677765722731,
  0.61787624440264388, 0.75540440835500322, 0.86563120238783220,
  0.94457502307323249, 0.98940093499165027
};
static const double gauss_w16[8] = {
  0.18945061045506834, 0.18260341504492425, 0.16915651939500245,
  0.14959598881657588, 0.12462897125553447, 0.09515851168249304,
  0.06225352393864833, 0.02715245941175411
};
// Panel geometry, in units of the conditional SD divided by |rho| -- the width
// of the transition the integrand makes at each edge. GAUSS_EDGE brackets each
// edge; GAUSS_WINDOW bounds the region where the integrand is not negligible.
// Both are measured, not guessed: at GAUSS_WINDOW = 6 the window truncates
// real mass and the worst-case relative error over a 3,600-cell grid is 1.0,
// while 10 and 14 agree to every digit; GAUSS_EDGE = 4 with 16 nodes holds the
// worst case to 7e-7 out to |rho| = 0.99999.
#define GAUSS_EDGE 4.0
#define GAUSS_WINDOW 10.0

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

// NB2 distribution function by direct summation of the mass function, writing
// F(y), F(y - 1) and P(Y = y) from one recursion because the discrete-copula
// likelihood needs all three.
//
// P(Y = y) is returned in its own right rather than left to the caller as
// F(y) - F(y - 1). In the far tail both CDFs have saturated at 1 and their
// difference is exactly zero in double precision, while the recursion still
// holds the mass to full relative precision.
//
// The textbook route is F(y) = pbeta(r / (r + mu), r, y + 1), and that is what
// this template used until the Laplace estimator exposed it.  TMB's pbeta()
// wraps TOMS 708, whose branches truncate a shape parameter to its integer part
// and reassign the remainder as a constant.  The VALUE stays accurate; the
// derivatives taken through those branches do not.  Two failures follow:
//
//   * the first derivative with respect to shape1 collapses to exactly zero at
//     some ordinary arguments -- including r = 2, which is where this package's
//     own default dispersion start (m = 0.5) puts it -- so the score for
//     log_m1/log_m2 was silently wrong there;
//   * third derivatives are NaN over wide regions of ordinary parameter values,
//     and the Laplace outer gradient differentiates the joint negative
//     log-likelihood three times.  Every copula family therefore failed under
//     method = "laplace" with "inner newton optimization failed during gradient
//     calculation", one or two nlminb steps in.
//
// Summing the mass function is exact and differentiable to every order.  Cost
// is O(y) per margin per evaluation, against one atomic call; on the truck
// workload that is about 12 extra terms per observation.
//
// Seeding the recursion with the linear-space P(Y = 0) is unsound: that term
// underflows to exact 0 whenever r * log(p) drops below about -745 (e.g.
// mu = r = 2000, an ordinary low-dispersion/high-mean region -- p^r =
// 0.5^2000), and once the seed is zero every later term stays zero too,
// because each step only ever multiplies the previous one. The recursion
// then reports P(Y = y) = 0 even where the true mass near the mode is
// perfectly representable (0.0063 at mu = r = y = 2000), handing the
// optimizer a flat, wrong objective instead of the real likelihood.
//
// Tracking the same recursion in log space avoids this: log P(Y = 0) =
// r * log(p) is an ordinary finite double even when P(Y = 0) itself
// underflows, and each step adds a log increment that is well-conditioned
// wherever the true mass is well-conditioned, regardless of how small the
// k = 0 term was. The running cdf is accumulated with a log-sum-exp so it
// never has to pass through the linear-space representation of a term
// until the final result -- which only underflows to 0 when the true
// probability truly is negligible.
template<class Type>
Type log_add_exp(Type a, Type b) {
  Type hi = CppAD::CondExpGt(a, b, a, b);
  Type lo = CppAD::CondExpGt(a, b, b, a);
  return hi + stable_log1p(exp(lo - hi));
}

// The log-space accumulators are returned alongside the linear ones because
// Clayton needs them. Its cell probability is built from a^-theta, and with
// theta reaching 4.85e8 that quantity only stays finite in logs; exp()ing the
// CDF first and taking the log again would also lose every cell whose CDF has
// underflowed, which is exactly the regime this recursion exists to keep.
// log_cdf_ym1 is the log of P(Y <= y - 1) for y > 0 and is left at log P(Y = 0)
// for y = 0, where the caller must not use it -- there is no representable log
// of zero to return, and every caller selects that case on the observed count.
template<class Type>
void nb2_cdf_pair(int y, Type mu, Type r,
                  Type &cdf_y, Type &cdf_ym1, Type &pmf_y,
                  Type &log_cdf_y, Type &log_cdf_ym1, Type &log_pmf_y) {
  Type log_p = log(r) - log(r + mu);
  Type log_q = log(mu) - log(r + mu);
  Type log_term = r * log_p;  // log P(Y = 0)
  Type log_cum = log_term;
  Type lcm1 = log_term;
  cdf_ym1 = Type(0);
  for (int k = 1; k <= y; k++) {
    lcm1 = log_cum;
    // log P(Y = k) = log P(Y = k - 1) + log(r + k - 1) - log(k) + log(q)
    log_term += log(r + Type(k - 1)) - log(Type(k)) + log_q;
    log_cum = log_add_exp(log_cum, log_term);
  }
  if (y > 0) cdf_ym1 = exp(lcm1);
  cdf_y = exp(log_cum);
  pmf_y = exp(log_term);
  log_cdf_y = log_cum;
  log_cdf_ym1 = lcm1;
  log_pmf_y = log_term;
}

// Frank's joint probability of the cell (a', a] x (b', b], evaluated as ONE
// log1p rather than as the second difference
// C(a,b) - C(a',b) - C(a,b') + C(a',b').
//
// With C(u,v) = -log1p(A(u) B(v) / D) / th, A(u) = expm1(-th u) and
// D = expm1(-th), the four logarithms telescope exactly:
//
//   p     = -log1p( dA * dB / (D * M) ) / th
//   dA    = A(a) - A(a') = exp(-th a') * expm1(-th * pmf_a)
//   M     = (1 + A(a') B(b) / D) * (1 + A(a) B(b') / D)
//
// so every cancelling difference becomes an expm1/log1p of a small argument,
// and dA is built from the marginal mass pmf_a directly instead of from two
// CDFs that have both saturated at 1.
//
// The naive second difference is unusable in the tail. At the truck fit's own
// starting values -- all slopes zero, so mu = 1 against counts running to 266
// -- it returns pure rounding noise for counts from about 26 and exactly zero
// above about 40, where the true probabilities are 1e-13 and 1e-20. Under SML
// that noise only corrupts the objective; under Laplace it puts negative
// curvature into the inner Hessian (161 of 27,896 latent rows on the truck
// data), and TMB's inner Newton cannot take even its first step.
//
// This is Frank-specific: it relies on C being a log of a bilinear form in
// A(u) and B(v). The Gaussian and Clayton branches below still take the naive
// second difference and remain subject to the same cancellation.
template<class Type>
Type frank_cell_prob(Type a, Type am, Type pmf_a,
                     Type b, Type bm, Type pmf_b, Type th) {
  Type signed_eps = CppAD::CondExpGe(th, Type(0), Type(1e-5), Type(-1e-5));
  Type safe_th = CppAD::CondExpLt(fabs(th), Type(1e-5), signed_eps, th);

  Type D = stable_expm1(-safe_th);
  Type A_a = stable_expm1(-safe_th * a);
  Type A_am = stable_expm1(-safe_th * am);
  Type B_b = stable_expm1(-safe_th * b);
  Type B_bm = stable_expm1(-safe_th * bm);
  Type dA = exp(-safe_th * am) * stable_expm1(-safe_th * pmf_a);
  Type dB = exp(-safe_th * bm) * stable_expm1(-safe_th * pmf_b);
  // Each factor is exp(-th * C(.,.)) and so is positive for any admissible
  // Frank argument; the guard below mirrors the one the naive form carried.
  Type M = (Type(1) + A_am * B_b / D) * (Type(1) + A_a * B_bm / D);
  Type ratio = dA * dB / (D * M);
  ratio = CppAD::CondExpLt(ratio, Type(-1.0 + 1e-15),
                           Type(-1.0 + 1e-15), ratio);
  Type regular = -stable_log1p(ratio) / safe_th;

  // Second difference of the near-independence expansion the naive form used,
  // C(u,v) ~ u v + th u v (1-u) (1-v) / 2, which telescopes to this in closed
  // form and so carries no cancellation either.
  Type near_independence = pmf_a * pmf_b *
    (Type(1) + th * (Type(1) - a - am) * (Type(1) - b - bm) / Type(2));

  return CppAD::CondExpLt(fabs(th), Type(1e-5), near_independence, regular);
}

// Clayton's joint probability of the cell (a', a] x (b', b], rearranged so
// that no step subtracts two nearly-equal numbers.
//
// Clayton has the same tail-cancellation defect the naive second difference
// gave Frank, and it is not hypothetical.  On the truck data at the fit's own
// starting values -- all slopes zero, so mu = 1 against counts running to 266
// -- C(a,b) - C(a',b) - C(a,b') + C(a',b') returns a NON-POSITIVE probability
// for 243 of 3,487 observations and a strictly NEGATIVE one for 11, where the
// true cell probabilities run down to 1e-136.  Under SML those only corrupt
// the objective (each clamped cell contributes the 1e-300 floor, about 691
// nats).  Under Laplace the negative ones put negative curvature into the
// inner Hessian, and TMB's inner Newton cannot take even its first step:
// method = "laplace" with copula("kimeldorf") failed outright with "inner
// newton optimization failed during gradient calculation".
//
// Write C(u,v) = s^k with s = 1 + A(u) + A(v), A(u) = u^-th - 1 and
// k = -1/th.  Factoring out the corner s00 = 1 + A(a) + A(b) and setting
// x = dA / s00, y = dB / s00 with dA = A(a') - A(a) >= 0 leaves
//
//   p  = s00^k * [ exp(u2) * expm1(u1 - u2) + ex * expm1(u1) ]
//   ex = expm1(k * log1p(x))
//   u1 = k * log1p(y / (1 + x))
//   u2 = k * log1p(y)
//   u1 - u2 = k * log1p(-x*y / ((1 + x)(1 + y)))
//
// The last identity is what removes the cancellation.  The naive difference
// has to recover an O(xy) second-order term by subtracting four O(1)
// quantities; here that term is written in closed form.  Because k < 0 and
// x, y >= 0, every factor above is sign-determined and both bracket terms are
// positive, so p is positive BY CONSTRUCTION rather than by clamping -- which
// is exactly what the inner Newton needs.
//
// x and y are carried as LOGARITHMS throughout.  th = exp(z_dep) with z_dep
// clamped to [-20, 20], so th reaches 4.85e8 and a^-th = exp(-th log a)
// overflows a double for any a bounded away from 1 -- log x, by contrast, is
// an ordinary number near 1e8.  Every quantity the result needs is a function
// of log1p(x), log1p(y) and log(1 + x + y), each of which follows from log x
// and log y by log-sum-exp, so x and y themselves are never formed.
//
// An earlier version instead capped x and y at 1e15, on the reasoning that
// (1 + x)^k is indistinguishable from x^k beyond that point.  That reasoning
// was wrong, and the error was not small.  The cap is applied BEFORE the
// result raises the ratio to the power k = -1/th, and k * log x is materially
// different from k * log(1e15) whenever th is large: capping loses
// |k| * (log x - 34.5) in the exponent.  For the symmetric cell y1 = y2 = 1 at
// mu = 1, m = 0.5, the exact probability at z_dep = 5 is 0.29077 and the
// capped form returned 0.15036.  Worse, Clayton tends to the comonotonic
// copula as th grows, so this cell must approach P(Y = 1) = 0.29630; the cap
// instead drove it to 2.4e-10 at z_dep = 20, reversing the likelihood's
// behaviour across the whole strong-dependence region that R/inference.R
// reports as interior.
//
// The one place the cap did real work was keeping 1 - xy/((1+x)(1+y)) off the
// rounding floor.  That is no longer needed either, because the quantity has
// an exact closed form with no subtraction at all:
//
//   1 - xy/((1+x)(1+y)) = (1 + x + y) / ((1+x)(1+y))
//
// so log1p(-xy/((1+x)(1+y))) = log(1+x+y) - log1p(x) - log1p(y), which is what
// the code computes.
//
// The y = 0 branches are separate because A(0) is infinite: there the cell is
// bounded by the axis and the second difference degenerates to a first
// difference.  They are selected from the OBSERVED COUNTS, passed in as
// y1_zero/y2_zero.  They previously tested asDouble(am) == 0.0, which is a
// different thing in two ways: asDouble() resolves when the tape is built, so
// the branch was frozen at whatever the starting values implied; and am is a
// parameter-dependent CDF that underflows to exactly zero for positive counts
// in ordinary regions (P(Y <= 0) at mu = r = 2000 is 0.5^2000).  A tape built
// at such a start kept the axis formula permanently, and the same parameter
// vector then scored differently depending on where its tape had been made --
// 2.256 against 1.592 nats on a single observation.
template<class Type>
Type clayton_cell_prob(Type log_a, Type log_am, Type log_pmf_a,
                       Type log_b, Type log_bm, Type log_pmf_b,
                       Type th, bool y1_zero, bool y2_zero) {
  auto vmax = [](Type u, Type v) { return CppAD::CondExpGt(u, v, u, v); };
  Type k = Type(-1.0) / th;

  // log s00, s00 = a^-th + b^-th - 1 >= 1 (a, b <= 1 keeps both terms >= 1).
  Type La = -th * log_a;
  Type Lb = -th * log_b;
  Type M = vmax(La, Lb);
  Type log_s00 = M + log(exp(La - M) + exp(Lb - M) - exp(-M));
  Type C00 = exp(k * log_s00);  // C(a, b)

  // log x for x = (u'^-th - u^-th) / s00 > 0.
  //
  // The ratio log(u/u') is built from the MARGINAL MASS, not from
  // log_u - log_um.  That difference is the whole reason nb2_cdf_pair()
  // returns a mass at all: once the CDF saturates, log F(y) and log F(y - 1)
  // are equal to the last bit and their difference is exactly zero, which
  // sends log() to -infinity and every derivative through it to NaN.  On the
  // truck margins at mu = 1 that happens for all 197 counts from y = 70 up,
  // covering 50 of the 3,487 observations, and the difference is already 32%
  // wrong at y = 69 before it collapses.  The mass route stays finite and
  // accurate there: at y = 266 it gives 1.45e-125 where the difference gives
  // 0.  This is the same cancellation the CDF-difference form of this
  // function was written to remove, and computing the ratio from two log CDFs
  // reintroduced it in log space.
  //
  // dA = u'^-th - u^-th = u'^-th * (1 - exp(-th * log(u/u'))), so with
  // LR = log(u/u') = log1p(pmf_u / u') the bracket is -expm1(-th * LR).
  auto log_ratio = [&](Type log_um, Type log_pmf_u) -> Type {
    Type t = log_pmf_u - log_um;
    Type LR = log_add_exp(Type(0), t);        // log(u / u') > 0
    // log LR without forming LR when the mass is far below the CDF: there
    // LR -> exp(t), so log LR -> t.  Keeps the small-S branch below finite
    // even if LR itself underflows.
    Type log_LR = CppAD::CondExpLt(t, Type(-30), t, log(vmax(LR, Type(1e-300))));
    Type S = th * LR;
    Type log_S = log(th) + log_LR;
    // Both CondExp branches are evaluated, so the exact branch gets a floored
    // argument; for S below the switch, -expm1(-S)/S differs from 1 by less
    // than half an ulp, so log_S is the accurate form anyway.
    Type S_safe = vmax(S, Type(1e-300));
    Type log_bracket = CppAD::CondExpLt(
      S, Type(1e-8), log_S, log(-stable_expm1(-S_safe))
    );
    return -th * log_um + log_bracket - log_s00;
  };

  if (y1_zero && y2_zero) return C00;

  if (y2_zero) {                      // cell bounded by the y2 axis
    Type L1x = log_add_exp(Type(0), log_ratio(log_am, log_pmf_a));
    return -C00 * stable_expm1(k * L1x);
  }
  Type log_y = log_ratio(log_bm, log_pmf_b);
  Type L1y = log_add_exp(Type(0), log_y);
  if (y1_zero) return -C00 * stable_expm1(k * L1y);

  Type log_x = log_ratio(log_am, log_pmf_a);
  Type L1x = log_add_exp(Type(0), log_x);

  // u1 = k * log1p(y / (1 + x)) and du = k * log1p(-w),
  // w = xy / ((1+x)(1+y)).  Both are written so that neither end of the range
  // cancels, which needs more care than it looks.
  //
  // Writing u1 as k * (log(1+x+y) - log1p(x)) is exact on paper and useless in
  // arithmetic: for x >> y both logs are O(x) and their difference is O(y), so
  // y is lost entirely.  The log-sum-exp form below never forms that
  // difference.
  //
  // du is the term that carries the cell.  Expanding the bracket for small
  // x, y gives (1+x+y)^k - (1+x)^k - (1+y)^k + 1 = k(k-1)xy + O(3), and the
  // two terms of the return are -k*xy and k^2*xy -- the SAME order.  Dropping
  // either one is not a rounding error but a factor of k/(k-1).  Computing
  // du as log(1+x+y) - log1p(x) - log1p(y) does exactly that: all three are
  // O(x+y) and their difference is O(xy), so it underflows to zero and the
  // cell comes back as C00*k^2*xy instead of C00*k(k-1)*xy -- a factor of 3.7
  // at theta = e, and it gets worse as theta falls.
  //
  // So the small-w side is computed from log1p(-w) directly, with log w built
  // additively (no cancellation), and only the large-w side -- where x and y
  // are big and the three logs are well separated -- uses the difference.
  Type u1 = k * log_add_exp(Type(0), log_y - L1x);
  Type u2 = k * L1y;
  Type log_w = log_x + log_y - L1x - L1y;      // log of xy/((1+x)(1+y)) < 0
  Type Lsum = log_add_exp(L1x, log_y);         // log(1 + x + y)
  // Both CondExp branches are evaluated, so the log1p branch gets an argument
  // floored away from -1 even when it is not the one selected.
  Type w_safe = exp(CppAD::CondExpGt(log_w, Type(-0.7), Type(-0.7), log_w));
  Type du = k * CppAD::CondExpLt(log_w, Type(-0.7),
                                 stable_log1p(-w_safe),
                                 Lsum - L1x - L1y);
  Type ex = stable_expm1(k * L1x);
  return C00 * (exp(u2) * stable_expm1(du) + ex * stable_expm1(u1));
}

// Gaussian's joint probability of the cell (a', a] x (b', b], evaluated as ONE
// strip integral instead of the second difference of four corner CDFs.
//
// Gaussian was the last family still taking the naive second difference, and
// it cancels there exactly as Frank and Clayton did.  On the truck data at the
// fit's own starting values (mu = 1 against counts running to 266) the
// four-corner form returns a non-positive cell probability for 457 to 600 of
// the 3,487 observations depending on rho, of which 212 to 355 are strictly
// NEGATIVE -- the negative-curvature source that stops TMB's inner Newton
// under method = "laplace".  The strip integral below leaves 245 floored
// cells and no negative ones, which lowers the objective at those starting
// values from about 334,000-432,000 to about 196,000-203,000.
//
// Conditioning the bivariate normal on the first margin gives
//
//   P = int_{q(a')}^{q(a)} phi(z) [ Phi((q(b) - rho z)/s)
//                                 - Phi((q(b') - rho z)/s) ] dz
//
// with s = sqrt(1 - rho^2).  The integrand is a product of non-negative
// factors (q(b) >= q(b') makes the bracket non-negative) and the limits are
// ordered, so P is non-negative BY CONSTRUCTION rather than by clamping.
//
// Integrating in z rather than in the probability variable matters: the same
// rule applied to int_{a'}^{a} ... dt has to evaluate qnorm() near the ends of
// the interval, where it is singular, and Gauss-Legendre then converges only
// as O(1/n) -- 1.5e-4 relative error at 20 points, against 1.5e-10 for the
// form below.  The probability-space version buys an exactly-known interval
// width (the marginal mass) at the cost of that singularity; it is not worth
// the trade here.
//
// The inner difference is taken on whichever tail is not saturated, since
// Phi(A) - Phi(B) cancels when both arguments are large and positive -- and
// they are: the truck data's second margin reaches counts of 47 against
// mu = 1, putting b at 1 - 1e-22.  The branch is applied to the ARGUMENTS
// rather than to the two pnorm() results: CppAD::CondExp evaluates both of its
// value branches, so selecting after the fact would put four pnorm() calls on
// the tape per node instead of two.
//
// A single fixed rule across the whole quantile interval is not enough, and
// this is the one part of the function that is not a matter of taste.  The
// integrand is a smoothed indicator of [q(b')/rho, q(b)/rho] whose edges have
// width s/|rho|.  As |rho| -> 1 that width collapses while the interval does
// not, and the nodes simply step over the plateau.  At rho = 0.9999 on an
// ordinary observation -- y = (0, 17), mu = (0.193, 17.6), a cell whose true
// probability is 0.0401, nowhere near any rounding floor -- one 20-point rule
// over [-7.94, 1.06] returns 1.7e-25 and loses 53.8 log-likelihood units on
// that observation alone.  Making the integrand non-negative had not made it
// accurate.  This is inside the supported domain: R/inference.R reports
// rho = 0.9999 as an interior estimate, not a boundary.
//
// So the interval is first cut down to the window where the integrand is not
// negligible, then split at +/- GAUSS_EDGE around each edge centre, giving
// five panels each smooth on its own scale.  The cut points are forced
// monotone by a running maximum, which both keeps every panel width
// non-negative (a negative width would subtract mass) and lets panels
// collapse to nothing when the plateau is narrower than the brackets or when
// the window clips an edge away.  The union of the panels is exactly the
// retained interval regardless of how many collapse.
//
// Dividing by rho needs |rho| bounded away from zero; the floor is set where
// the resulting window is already far wider than any quantile interval, so
// the clipping makes the choice of floor immaterial.
//
// Remaining limitation, deliberately not papered over: when a and a' are close
// enough that the safe_qnorm() clamp maps both to the same point (245 of the
// 3,487 truck cells at the starting values), q(a) == q(a') and this returns 0,
// which the caller floors.  Gaussian must pass through qnorm(), which is
// singular at 1, so once the NB2 CDF reaches the clamp the cell cannot be
// recovered from the CDF at all -- that would take a separately accumulated
// survival function.  Frank and Clayton are not affected because their
// generators stay analytic at u = 1.  Unlike the second difference, this at
// least fails to zero rather than to a negative number.
template<class Type>
Type gaussian_cell_prob(Type qa, Type qam, Type qb, Type qbm, Type rho) {
  Type sig2 = Type(1.0) - rho * rho;
  sig2 = CppAD::CondExpLt(sig2, Type(1e-12), Type(1e-12), sig2);
  Type sig = sqrt(sig2);

  auto vmax = [](Type u, Type v) { return CppAD::CondExpGt(u, v, u, v); };
  auto vmin = [](Type u, Type v) { return CppAD::CondExpLt(u, v, u, v); };

  // Signed rho with |rho| floored, so the edge centres below stay finite.
  Type rho_abs = vmax(rho, -rho);
  Type rho_faf = vmax(rho_abs, Type(1e-6));
  Type rho_sgn = CppAD::CondExpLt(rho, Type(0), -rho_faf, rho_faf);

  Type edge = sig / rho_faf;               // transition width, in z
  Type c1 = qb / rho_sgn;
  Type c2 = qbm / rho_sgn;
  Type clo = vmin(c1, c2);
  Type chi = vmax(c1, c2);

  Type lo = vmax(qam, clo - Type(GAUSS_WINDOW) * edge);
  Type hi = vmin(qa, chi + Type(GAUSS_WINDOW) * edge);
  hi = vmax(hi, lo);                       // empty window => every panel empty

  // Cut points clamped into [lo, hi] and then forced non-decreasing.
  Type cuts[6];
  cuts[0] = lo;
  cuts[1] = clo - Type(GAUSS_EDGE) * edge;
  cuts[2] = clo + Type(GAUSS_EDGE) * edge;
  cuts[3] = chi - Type(GAUSS_EDGE) * edge;
  cuts[4] = chi + Type(GAUSS_EDGE) * edge;
  cuts[5] = hi;
  for (int j = 1; j < 6; j++) {
    cuts[j] = vmin(vmax(cuts[j], lo), hi);
    cuts[j] = vmax(cuts[j], cuts[j - 1]);
  }

  Type total = Type(0);
  for (int p = 0; p < 5; p++) {
    Type half = (cuts[p + 1] - cuts[p]) / Type(2);
    Type mid = (cuts[p + 1] + cuts[p]) / Type(2);
    Type acc = Type(0);
    for (int i = 0; i < 8; i++) {
      for (int side = -1; side <= 1; side += 2) {
        Type z = mid + half * Type(side) * Type(gauss_x16[i]);
        Type A = (qb - rho * z) / sig;
        Type B = (qbm - rho * z) / sig;
        Type flip = A + B;
        Type Au = CppAD::CondExpGt(flip, Type(0), -B, A);
        Type Bu = CppAD::CondExpGt(flip, Type(0), -A, B);
        acc += Type(gauss_w16[i]) *
          dnorm(z, Type(0), Type(1), false) * (pnorm(Au) - pnorm(Bu));
      }
    }
    total += half * acc;
  }
  return total;
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
  // 0 = simulated maximum likelihood (Halton draws)
  // 1 = Laplace approximation (latent u1/u2 integrated by TMB)
  DATA_INTEGER(est_method);

  // ---- Parameters ----
  PARAMETER_VECTOR(beta1);
  PARAMETER_VECTOR(beta2);
  PARAMETER_VECTOR(log_sd1);
  PARAMETER_VECTOR(log_sd2);
  PARAMETER(log_m1);
  PARAMETER(log_m2);
  PARAMETER(z_dep);
  // Latent standard normals, one row per observation. Under est_method == 0
  // these are map-fixed at zero in R and never read; under est_method == 1
  // they are TMB random effects.
  PARAMETER_MATRIX(u1);
  PARAMETER_MATRIX(u2);

  // ---- Dimensions ----
  int n = Y1.size();
  int k1 = X1.cols(); (void)k1;
  int k2 = X2.cols(); (void)k2;
  int q1 = rand_idx1.size();
  int q2 = rand_idx2.size();
  int R = (q1 + q2 > 0) ? Z1.rows() : 1;
  // Laplace evaluates the conditional density once at the current latent
  // values; there is no draw dimension to average over.
  if (est_method == 1) R = 1;
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

  // Integer responses for the copula margins: nb2_cdf_pair() sums the mass
  // function up to y, so the loop bound must be an int rather than a Type.
  // Y1/Y2 are data and R has already checked they are whole and non-negative,
  // so this conversion costs nothing on the tape.
  vector<int> Y1_int(n), Y2_int(n);
  if (family >= FAM_FRANK) {
    for (int i = 0; i < n; i++) {
      Y1_int(i) = (int)asDouble(Y1(i));
      Y2_int(i) = (int)asDouble(Y2(i));
    }
  }

  // Precompute parameter-dependent deviations once per simulation draw.  The
  // matrices are read-only while observation contributions are accumulated.
  matrix<Type> dev1(R, q1), dev2(R, q2);
  if (est_method != 1) {
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
        int col = rand_idx1(j);
        // u_to_base() is deliberately NOT applied to the latent: it is a
        // per-distribution inverse CDF, not a general uniform-to-normal map,
        // and skipping it is only valid because u_to_base == qnorm for the
        // normal/lognormal distributions the Laplace path is restricted to.
        // fit_rpbnb_tmb() enforces that restriction on the R side by
        // rejecting method = "laplace" with uniform/triangular coefficients.
        Type d = (est_method == 1)
          ? compute_dev(beta1(col), sd1(j), u1(i, j), dist1(j), sign1(j))
          : dev1(r, j);
        eta1 += X1(i, col) * d;
      }
      for (int j = 0; j < q2; j++) {
        int col = rand_idx2(j);
        // Same restriction as the u1 loop above: u_to_base() is skipped
        // because the Laplace path only allows normal/lognormal
        // coefficients, for which u_to_base == qnorm.
        Type d = (est_method == 1)
          ? compute_dev(beta2(col), sd2(j), u2(i, j), dist2(j), sign2(j))
          : dev2(r, j);
        eta2 += X2(i, col) * d;
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
        // Y1_int/Y2_int are cast from data, so these two flags are ordinary
        // compile-time branches that carry no parameter dependence onto the
        // tape.  Clayton selects its axis cases from them; see
        // clayton_cell_prob() for why reading them off a CDF instead made the
        // taped objective depend on the starting values.
        const bool y1_zero = (Y1_int(i) == 0);
        const bool y2_zero = (Y2_int(i) == 0);

        Type a1, a1m, b1, b1m, pmf1, pmf2;
        Type la1, la1m, lpmf1, lb1, lb1m, lpmf2;
        if (pois1) {
          a1 = ppois(Y1(i), mu1);
          a1m = y1_zero ? Type(0) : ppois(Y1(i) - Type(1), mu1);
          pmf1 = dpois(Y1(i), mu1, false);
          la1 = log(a1);
          // a1m is exactly 0 when y = 0; log() of it is never read, because
          // every consumer selects that case on y1_zero.
          la1m = y1_zero ? Type(0) : log(a1m);
          lpmf1 = log(pmf1);
        } else {
          nb2_cdf_pair(Y1_int(i), mu1, r1, a1, a1m, pmf1, la1, la1m, lpmf1);
        }
        if (pois2) {
          b1 = ppois(Y2(i), mu2);
          b1m = y2_zero ? Type(0) : ppois(Y2(i) - Type(1), mu2);
          pmf2 = dpois(Y2(i), mu2, false);
          lb1 = log(b1);
          lb1m = y2_zero ? Type(0) : log(b1m);
          lpmf2 = log(pmf2);
        } else {
          nb2_cdf_pair(Y2_int(i), mu2, r2, b1, b1m, pmf2, lb1, lb1m, lpmf2);
        }

        // Every copula family now builds the cell probability directly rather
        // than as a second difference of corner CDFs, so none of them can
        // return a negative probability.  See frank_cell_prob(),
        // clayton_cell_prob() and gaussian_cell_prob().
        Type p_obs = Type(0);
        if (family == FAM_FRANK) {
          p_obs = frank_cell_prob(a1, a1m, pmf1, b1, b1m, pmf2, theta);
        } else if (family == FAM_GAUSSIAN) {
          auto safe_qnorm = [](Type p) -> Type {
            p = CppAD::CondExpLt(p, Type(1e-15), Type(1e-15), p);
            p = CppAD::CondExpGt(p, Type(1.0 - 1e-15),
                                 Type(1.0 - 1e-15), p);
            return qnorm(p);
          };
          p_obs = gaussian_cell_prob(safe_qnorm(a1), safe_qnorm(a1m),
                                     safe_qnorm(b1), safe_qnorm(b1m), rho);
        } else {
          p_obs = clayton_cell_prob(la1, la1m, lpmf1, lb1, lb1m, lpmf2,
                                    theta, y1_zero, y2_zero);
        }
        p_obs = CppAD::CondExpLt(p_obs, Type(1e-300),
                                 Type(1e-300), p_obs);
        log_draw(r) = log(p_obs);
      }
    }

    if (est_method == 1) {
      Type obs_ll = log_draw(0);
      for (int j = 0; j < q1; j++)
        obs_ll += dnorm(u1(i, j), Type(0), Type(1), true);
      for (int j = 0; j < q2; j++)
        obs_ll += dnorm(u2(i, j), Type(0), Type(1), true);
      nll -= obs_ll;
    } else {
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
