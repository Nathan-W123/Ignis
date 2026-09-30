// SPDX-License-Identifier: MIT
#include "ignis/nozzle/NozzleFlow.hpp"

#include "ignis/nozzle/Atmosphere.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <iomanip>
#include <sstream>

namespace ignis {
namespace {

/// Bracketed root find using the Illinois variant of regula falsi: it keeps the
/// bracket (so it cannot diverge) but converges super-linearly, which matters
/// because every function evaluation here is a Gibbs minimisation.
template <typename F>
double illinois(F&& f, double lo, double hi, double f_lo, double f_hi, double x_tol,
                double f_tol, int max_iter = 100) {
  double a = lo, b = hi, fa = f_lo, fb = f_hi;
  double x = 0.5 * (a + b);
  for (int i = 0; i < max_iter; ++i) {
    if (fb == fa) break;
    x = b - fb * (b - a) / (fb - fa);
    // Keep the secant step inside the bracket with a small safety margin.
    const double lo_guard = a + 0.01 * (b - a);
    const double hi_guard = b - 0.01 * (b - a);
    if (!(x > lo_guard && x < hi_guard)) x = 0.5 * (a + b);
    const double fx = f(x);
    if (std::abs(fx) <= f_tol || (b - a) <= x_tol) return x;
    if (fx * fb < 0.0) {
      a = b;
      fa = fb;
    } else {
      fa *= 0.5;  // Illinois down-weighting of the stagnant endpoint
    }
    b = x;
    fb = fx;
    if (a > b) { std::swap(a, b); std::swap(fa, fb); }
  }
  return x;
}

}  // namespace

std::string toString(ExpansionRegime r) {
  switch (r) {
    case ExpansionRegime::kUnderExpanded: return "under-expanded";
    case ExpansionRegime::kIdeallyExpanded: return "ideally-expanded";
    case ExpansionRegime::kOverExpanded: return "over-expanded";
    case ExpansionRegime::kSubsonicExit: return "subsonic-exit";
  }
  return "?";
}

std::string toString(SeparationCriterion c) {
  switch (c) {
    case SeparationCriterion::kNone: return "none";
    case SeparationCriterion::kSummerfield: return "summerfield";
    case SeparationCriterion::kSchmucker: return "schmucker";
  }
  return "?";
}

SeparationCriterion separationCriterionFromString(const std::string& s) {
  if (s == "none") return SeparationCriterion::kNone;
  if (s == "summerfield") return SeparationCriterion::kSummerfield;
  if (s == "schmucker") return SeparationCriterion::kSchmucker;
  throw ConfigError("separation criterion must be 'none', 'summerfield' or 'schmucker', got '" +
                    s + "'");
}

// ---------------------------------------------------------------------------
// NozzleFlow
// ---------------------------------------------------------------------------

NozzleFlow::NozzleFlow(const EquilibriumSolver& solver, CompositionModel model,
                       ChamberReference ref)
    : solver_(&solver), model_(model), ref_(std::move(ref)) {
  h0_ = ref_.stagnation.h;
  s0_ = ref_.stagnation.s;
  p0_ = ref_.stagnation.p;
  if (!(p0_ > 0.0)) throw ConfigError("nozzle flow: chamber pressure must be positive");

  // Locate the sonic point by bisection on M(p) - 1.  M decreases monotonically
  // with pressure, so the bracket [p_lo, p_hi] is easy to establish.
  double p_hi = 0.95 * p0_;
  double p_lo = 0.2 * p0_;
  auto mach = [&](double p) { return atPressure(p).mach; };
  double m_hi = mach(p_hi);
  if (m_hi > 1.0) {
    // very low gamma: sonic point sits above 0.95 p0
    for (int i = 0; i < 40 && m_hi > 1.0; ++i) {
      p_lo = p_hi;
      p_hi = p0_ - 0.5 * (p0_ - p_hi);
      m_hi = mach(p_hi);
    }
  }
  double m_lo = mach(p_lo);
  for (int i = 0; i < 60 && m_lo < 1.0; ++i) {
    p_hi = p_lo;
    p_lo *= 0.7;
    m_lo = mach(p_lo);
  }
  if (!(m_lo >= 1.0 && m_hi <= 1.0))
    throw ConvergenceError("nozzle flow: could not bracket the sonic point");

  const double p_sonic = illinois([&](double p) { return mach(p) - 1.0; }, p_lo, p_hi,
                                  m_lo - 1.0, m_hi - 1.0, 1.0e-13 * p0_, 1.0e-12);
  throat_ = atPressure(p_sonic);
  throat_.area_ratio = 1.0;
  if (std::abs(throat_.mach - 1.0) > 1.0e-7) {
    std::ostringstream os;
    os << "nozzle flow: throat search converged to M = " << throat_.mach << ", not 1";
    throw ConvergenceError(os.str());
  }

}


void NozzleFlow::computeFloor() const {
  if (floor_known_) return;
  // Walk down in halvings until the expansion leaves the validated property
  // range, then refine geometrically.  This is only needed when a caller asks
  // how far the nozzle *could* expand; ordinary station queries bracket
  // themselves and never reach here.
  double ok = throat_.gas.p, bad = 0.0, probe = throat_.gas.p;
  for (int i = 0; i < 200; ++i) {
    probe *= 0.5;
    bool good = true;
    try {
      (void)atPressure(probe);
    } catch (const IgnisError&) {
      good = false;
    }
    if (good) {
      ok = probe;
      if (probe < 1.0e-12 * p0_) break;
    } else {
      bad = probe;
      break;
    }
  }
  if (bad > 0.0) {
    for (int i = 0; i < 40 && ok / bad > 1.02; ++i) {
      const double mid = std::sqrt(ok * bad);
      bool good = true;
      try {
        (void)atPressure(mid);
      } catch (const IgnisError&) {
        good = false;
      }
      if (good) ok = mid; else bad = mid;
    }
  }
  p_floor_ = ok;
  max_area_ratio_ = throat_.mass_flux / atPressure(p_floor_).mass_flux;
  floor_known_ = true;
}

ExpansionState NozzleFlow::buildState(const GasState& gas) const {
  ExpansionState st;
  st.gas = gas;
  const double dh = h0_ - gas.h;
  if (dh < 0.0) {
    // Numerically this only happens within round-off of the chamber itself.
    if (dh > -1.0e-6 * std::max(1.0, std::abs(h0_))) {
      st.u = 0.0;
    } else {
      std::ostringstream os;
      os << "nozzle flow: static enthalpy exceeds the stagnation value by " << -dh
         << " J/kg at p = " << gas.p << " Pa -- the expansion is not physical";
      throw InfeasibleError(os.str());
    }
  } else {
    st.u = std::sqrt(2.0 * dh);
  }
  st.mach = st.u / gas.a;
  st.mass_flux = gas.rho * st.u;
  st.supersonic = st.mach > 1.0;
  st.area_ratio = (throat_.mass_flux > 0.0) ? throat_.mass_flux / st.mass_flux : 0.0;
  return st;
}

ExpansionState NozzleFlow::atPressure(double p) const {
  if (!(p > 0.0)) throw RangeError("nozzle flow: static pressure must be positive");
  if (p > p0_ * (1.0 + 1e-12))
    throw RangeError("nozzle flow: static pressure exceeds the chamber pressure");

  if (model_ == CompositionModel::kFrozen) {
    const auto& mix = solver_->mixture();
    const double T = mix.temperatureFromEntropy(ref_.stagnation.n, s0_, p, ref_.stagnation.T);
    return buildState(mix.frozenState(ref_.stagnation.n, T, p));
  }
  // Warm start from the previous solve: successive pressures in a station march
  // or a root find are close, and the Gibbs solver converges in a handful of
  // Newton steps instead of twenty from the chamber composition.
  const Eigen::VectorXd* guess = has_last_ ? &last_n_ : &ref_.stagnation.n;
  const double T_guess = has_last_ ? last_T_ : ref_.stagnation.T;
  const auto r = solver_->sp(ref_.b, s0_, p, T_guess, guess);
  last_n_ = r.state.n;
  last_T_ = r.state.T;
  has_last_ = true;
  return buildState(r.state);
}

ExpansionState NozzleFlow::atAreaRatio(double area_ratio, bool supersonic) const {
  if (!(area_ratio >= 1.0))
    throw RangeError("nozzle flow: area ratio must be at least 1, got " +
                     std::to_string(area_ratio));
  if (std::abs(area_ratio - 1.0) < 1e-12) return throat_;

  const double p_star = throat_.gas.p;
  const double target = throat_.mass_flux / area_ratio;  // required rho u

  auto flux = [&](double p) { return atPressure(p).mass_flux; };

  double lo, hi, f_lo, f_hi;
  if (supersonic) {
    // Mass flux falls monotonically as the pressure drops on the supersonic
    // branch, so the bracket is found by halving from the throat.  Descending
    // only as far as this particular area ratio needs keeps the search away
    // from the low-pressure end of the property data, which is both expensive
    // to probe and irrelevant here.
    hi = p_star;
    f_hi = flux(hi) - target;
    lo = p_star;
    f_lo = f_hi;
    bool bracketed = false;
    for (int i = 0; i < 200; ++i) {
      lo *= 0.5;
      try {
        f_lo = flux(lo) - target;
      } catch (const IgnisError& e) {
        computeFloor();
        std::ostringstream os;
        os << "nozzle flow: area ratio " << area_ratio << " exceeds the largest expansion the "
           << "species data supports (" << max_area_ratio_ << ", reached at p = " << p_floor_
           << " Pa where the static temperature hits the lower limit of the NASA polynomial "
           << "fits). Reduce the expansion ratio or extend the species database. Underlying "
           << "limit: " << e.what();
        throw InfeasibleError(os.str());
      }
      if (f_lo < 0.0) { bracketed = true; break; }
      hi = lo;
      f_hi = f_lo;
    }
    if (!bracketed)
      throw ConvergenceError("nozzle flow: could not bracket the supersonic solution for A/At = " +
                             std::to_string(area_ratio));
  } else {
    lo = p_star;
    hi = p0_ * (1.0 - 1.0e-12);
    f_lo = flux(lo) - target;
    f_hi = flux(hi) - target;
    if (f_lo * f_hi > 0.0) {
      std::ostringstream os;
      os << "nozzle flow: could not bracket the subsonic solution for A/At = " << area_ratio
         << " (residuals " << f_lo << ", " << f_hi << ")";
      throw ConvergenceError(os.str());
    }
  }
  const double p_sol = illinois([&](double p) { return flux(p) - target; }, lo, hi, f_lo, f_hi,
                                1.0e-13 * p0_, 1.0e-11 * target);
  ExpansionState st = atPressure(p_sol);
  st.area_ratio = area_ratio;  // exact by construction of the search target
  return st;
}

double NozzleFlow::massFlowResidual(const ExpansionState& st) const {
  const double mdot_ref = throat_.mass_flux;                 // per unit throat area
  const double mdot_here = st.mass_flux * st.area_ratio;     // per unit throat area
  return std::abs(mdot_here - mdot_ref) / mdot_ref;
}

double NozzleFlow::energyResidual(const ExpansionState& st) const {
  return std::abs(st.gas.h + 0.5 * st.u * st.u - h0_) / std::max(1.0, std::abs(h0_));
}

// ---------------------------------------------------------------------------
// Normal shock
// ---------------------------------------------------------------------------
namespace {

/// Stagnation pressure of a stream: the pressure reached by isentropically
/// bringing (s, h0) to rest.  Newton on dh/dp|_s = v.
double stagnationPressure(const GasMixture& mix, const Eigen::VectorXd& n, double s, double h0,
                          double T_guess, double p_guess) {
  double p = p_guess;
  for (int i = 0; i < 300; ++i) {
    const double T = mix.temperatureFromEntropy(n, s, p, T_guess);
    const double h = mix.enthalpy(n, T);
    const double err = h - h0;
    if (std::abs(err) < 1.0e-11 * std::max(1.0, std::abs(h0))) return p;
    const double M = mix.molarMass(n);
    const double v = (constants::R_universal / M) * T / p;   // specific volume
    double p_new = p - err / v;
    if (!(p_new > 0.0)) p_new = 0.5 * p;
    if (p_new > 1.0e4 * p_guess) p_new = 1.0e4 * p_guess;
    if (std::abs(p_new - p) < 1e-14 * p) return p_new;
    p = p_new;
  }
  throw ConvergenceError("stagnation pressure iteration did not converge");
}

}  // namespace

NormalShockResult solveNormalShock(const GasMixture& mix, const ExpansionState& up) {
  if (!(up.mach > 1.0))
    throw InfeasibleError("normal shock: the upstream state is not supersonic (M = " +
                          std::to_string(up.mach) + ")");
  const double rho1 = up.gas.rho, u1 = up.u, p1 = up.gas.p, h1 = up.gas.h;
  const Eigen::VectorXd n = up.gas.n;
  const double R_specific = up.gas.R;

  // Parameterise by v = u2/u1 = rho1/rho2 in (0, 1):
  //   p2 = p1 + rho1 u1^2 (1 - v)      (momentum)
  //   h2 = h1 + u1^2 (1 - v^2)/2       (energy)
  //   T2 = p2 v / (rho1 R)             (mass + ideal-gas EOS)
  // and require the caloric equation h(T2) = h2.  v = 1 is the trivial
  // (shock-free) root; the physical root is the other sign change.
  auto residual = [&](double v) { 
    const double p2 = p1 + rho1 * u1 * u1 * (1.0 - v);
    const double h2 = h1 + 0.5 * u1 * u1 * (1.0 - v * v);
    const double T2 = p2 * v / (rho1 * R_specific);
    return mix.enthalpy(n, T2) - h2;
  };

  const double v_top = 1.0 - 1.0e-7;
  double f_top;
  try {
    f_top = residual(v_top);
  } catch (const IgnisError& e) {
    throw ConvergenceError(std::string("normal shock: upstream branch unusable -- ") + e.what());
  }
  double lo = 0.0, hi = v_top, f_lo = 0.0, f_hi = f_top;
  bool bracketed = false;
  const int nscan = 600;
  for (int i = nscan - 1; i >= 1; --i) {
    const double v = v_top * i / nscan;
    double f;
    try {
      f = residual(v);
    } catch (const IgnisError&) {
      continue;  // T2 left the polynomial range -- keep scanning downward
    }
    if (f * f_hi < 0.0) { lo = v; f_lo = f; bracketed = true; break; }
    hi = v;
    f_hi = f;
  }
  if (!bracketed) {
    std::ostringstream os;
    os << "normal shock: no downstream branch found for M = " << up.mach << " at p = "
       << p1 << " Pa";
    throw ConvergenceError(os.str());
  }
  for (int i = 0; i < 300; ++i) {
    const double vm = 0.5 * (lo + hi);
    const double fm = residual(vm);
    if (fm * f_lo > 0.0) { lo = vm; f_lo = fm; } else { hi = vm; f_hi = fm; }
    if (hi - lo < 1.0e-15) break;
  }
  const double v = 0.5 * (lo + hi);
  const double p2 = p1 + rho1 * u1 * u1 * (1.0 - v);
  const double T2 = p2 * v / (rho1 * R_specific);

  NormalShockResult out;
  out.upstream = up;
  out.downstream = mix.frozenState(n, T2, p2);
  out.u_downstream = v * u1;
  out.mach_downstream = out.u_downstream / out.downstream.a;
  if (out.mach_downstream >= 1.0)
    throw ConvergenceError("normal shock: the solver returned a supersonic downstream state");
  out.entropy_rise = out.downstream.s - up.gas.s;
  if (out.entropy_rise < -1.0e-9 * std::abs(up.gas.s))
    throw InfeasibleError("normal shock: the computed jump decreases entropy, which means the "
                          "upstream state or the property data is inconsistent");

  const double h01 = up.gas.h + 0.5 * u1 * u1;
  const double h02 = out.downstream.h + 0.5 * out.u_downstream * out.u_downstream;
  const double p01 = stagnationPressure(mix, n, up.gas.s, h01, up.gas.T, p1 * 2.0);
  const double p02 = stagnationPressure(mix, n, out.downstream.s, h02, out.downstream.T, p2 * 1.5);
  out.stagnation_pressure_ratio = p02 / p01;
  return out;
}

// ---------------------------------------------------------------------------
// Steady performance
// ---------------------------------------------------------------------------
namespace {

/// Static pressure at the nozzle exit for a normal shock standing at
/// `area_ratio_shock`, with the post-shock composition frozen.
double exitPressureAfterShock(const GasMixture& mix, const NormalShockResult& sh,
                              double area_ratio_shock, double area_ratio_exit, double& u_exit) {
  const Eigen::VectorXd n = sh.downstream.n;
  const double h0 = sh.downstream.h + 0.5 * sh.u_downstream * sh.u_downstream;
  const double s0 = sh.downstream.s;
  const double G_shock = sh.downstream.rho * sh.u_downstream;
  const double G_exit = G_shock * area_ratio_shock / area_ratio_exit;

  auto flux = [&](double p, double* u_out) {
    const double T = mix.temperatureFromEntropy(n, s0, p, sh.downstream.T);
    const double h = mix.enthalpy(n, T);
    const double dh = h0 - h;
    if (dh <= 0.0) { if (u_out) *u_out = 0.0; return 0.0; }
    const double u = std::sqrt(2.0 * dh);
    if (u_out) *u_out = u;
    return (p / (sh.downstream.R * T)) * u;
  };

  // Subsonic branch: raise p until the mass flux drops below the target.
  double lo = sh.downstream.p, hi = lo;
  bool found = false;
  for (int i = 0; i < 400; ++i) {
    hi *= 1.02;
    if (flux(hi, nullptr) < G_exit) { found = true; break; }
  }
  if (!found)
    throw ConvergenceError("post-shock expansion: could not bracket the subsonic exit state");
  for (int i = 0; i < 300; ++i) {
    const double pm = 0.5 * (lo + hi);
    if (flux(pm, nullptr) > G_exit) lo = pm; else hi = pm;
    if (hi - lo < 1.0e-12 * hi) break;
  }
  const double p = 0.5 * (lo + hi);
  flux(p, &u_exit);
  return p;
}

/// Wall pressure below which the empirical criterion predicts separation.
double separationPressure(SeparationCriterion c, double p_ambient, double mach) {
  switch (c) {
    case SeparationCriterion::kNone:
      return 0.0;
    case SeparationCriterion::kSummerfield:
      return 0.4 * p_ambient;
    case SeparationCriterion::kSchmucker: {
      const double arg = 1.88 * mach - 1.0;
      if (arg <= 0.0) return 0.0;
      return p_ambient * std::pow(arg, -0.64);
    }
  }
  return 0.0;
}

}  // namespace

NozzlePerformance evaluateNozzle(const NozzleFlow& flow, const NozzleGeometry& geom,
                                 double p_ambient, const NozzlePerformanceOptions& opts) {
  if (p_ambient < 0.0) throw ConfigError("nozzle: ambient pressure must be non-negative");
  if (!(opts.eta_c_star > 0.0 && opts.eta_c_star <= 1.0))
    throw ConfigError("nozzle: eta_c_star must lie in (0, 1]");

  NozzlePerformance perf;
  const GasState& cham = flow.chamber().stagnation;
  const GasMixture& mix = flow.mixture();
  perf.composition_model = toString(flow.model());
  perf.p_chamber = cham.p;
  perf.p_ambient = p_ambient;
  perf.throat_area = geom.throatArea();
  perf.exit_area = geom.exitArea();
  perf.area_ratio = geom.expansionRatio();
  perf.eta_c_star = opts.eta_c_star;
  perf.separation_criterion = opts.separation;

  perf.c_star_ideal = flow.cStarIdeal();
  perf.c_star = perf.eta_c_star * perf.c_star_ideal;
  perf.mdot = perf.p_chamber * perf.throat_area / perf.c_star;

  // Shock-free supersonic exit state -- always computed, since the vacuum
  // performance and the regime classification are defined from it.
  const ExpansionState exit_sup = flow.atAreaRatio(perf.area_ratio, true);
  perf.mass_flow_residual = flow.massFlowResidual(exit_sup);
  perf.energy_residual = flow.energyResidual(exit_sup);

  double p_exit = exit_sup.gas.p;
  double u_exit = exit_sup.u;
  double t_exit = exit_sup.gas.T;
  double mach_exit = exit_sup.mach;

  // --- internal normal shock -------------------------------------------
  if (opts.resolve_internal_shocks && p_ambient > exit_sup.gas.p * (1.0 + 1e-12)) {
    NormalShockResult sh_exit;
    bool have_exit_shock = false;
    try {
      sh_exit = solveNormalShock(mix, exit_sup);
      have_exit_shock = true;
    } catch (const IgnisError& e) {
      perf.shock_note = std::string("shock-at-exit solve failed: ") + e.what();
    }
    if (have_exit_shock && p_ambient >= sh_exit.downstream.p) {
      // A normal shock somewhere inside the nozzle can raise the exit pressure
      // to p_ambient.  Its post-shock exit pressure increases monotonically as
      // the shock moves upstream (weaker shock, more subsonic diffusion).
      auto exitPressureForShockAt = [&](double eps_s, double* u_out) {
        const ExpansionState up = flow.atAreaRatio(eps_s, true);
        const NormalShockResult sh = solveNormalShock(mix, up);
        double u = 0.0;
        const double pe = exitPressureAfterShock(mix, sh, eps_s, perf.area_ratio, u);
        if (u_out) *u_out = u;
        return pe;
      };
      // Moving the shock upstream raises the exit pressure, so the interior
      // position that matches the ambient pressure is bracketed by a shock at
      // the throat and a shock at the exit.
      double lo = 1.0 + 1e-6, hi = perf.area_ratio;
      double p_lo = 0.0;
      bool ok = true;
      try {
        p_lo = exitPressureForShockAt(lo, nullptr);
      } catch (const IgnisError&) {
        ok = false;
      }
      if (ok && p_ambient <= p_lo) {
        for (int i = 0; i < 200; ++i) {
          const double mid = 0.5 * (lo + hi);
          double pm;
          try {
            pm = exitPressureForShockAt(mid, nullptr);
          } catch (const IgnisError&) {
            lo = mid;
            continue;
          }
          if (pm > p_ambient) { lo = mid; p_lo = pm; } else { hi = mid; }
          if (hi - lo < 1e-10 * perf.area_ratio) break;
        }
        const double eps_s = 0.5 * (lo + hi);
        double u = 0.0;
        const double pe = exitPressureForShockAt(eps_s, &u);
        perf.shock_in_nozzle = true;
        perf.shock_area_ratio = eps_s;
        // Axial position of that area ratio in the divergent.
        perf.shock_x = geom.throatPosition();
        for (const auto& st : geom.stations())
          if (st.x >= geom.throatPosition() && st.area_ratio <= eps_s) perf.shock_x = st.x;
        p_exit = pe;
        u_exit = u;
        const ExpansionState up = flow.atAreaRatio(eps_s, true);
        const NormalShockResult sh = solveNormalShock(mix, up);
        const double T = mix.temperatureFromEntropy(sh.downstream.n, sh.downstream.s, pe,
                                                    sh.downstream.T);
        t_exit = T;
        mach_exit = u / mix.frozenState(sh.downstream.n, T, pe).a;
        perf.shock_note = "normal shock resolved inside the divergent section";
      } else {
        perf.shock_unresolved = true;
        perf.shock_note =
            "ambient pressure exceeds the highest exit pressure a single interior normal "
            "shock can produce; the real flow is fully subsonic in the divergent or the "
            "nozzle is unchoked. No supersonic solution is reported.";
      }
    } else if (have_exit_shock) {
      perf.shock_note =
          "over-expanded: the internal flow stays supersonic and the compression to ambient "
          "occurs outside the nozzle through oblique shocks (not modelled in quasi-1D)";
    }
  }

  perf.p_exit = p_exit;
  perf.t_exit = t_exit;
  perf.u_exit = u_exit;
  perf.mach_exit = mach_exit;

  // --- loss factors -----------------------------------------------------
  if (opts.auto_divergence) {
    const double theta = (geom.spec().divergent == DivergentType::kConical)
                             ? geom.spec().cone_half_angle
                             : (geom.spec().bell_exit_angle > 0.0 ? geom.spec().bell_exit_angle
                                                                  : 8.0);
    perf.lambda_divergence = 0.5 * (1.0 + std::cos(theta * constants::pi / 180.0));
  } else {
    perf.lambda_divergence = opts.lambda_divergence;
  }
  perf.eta_nozzle = opts.eta_nozzle;

  // --- thrust -----------------------------------------------------------
  perf.thrust_momentum = perf.mdot * perf.u_exit;
  perf.thrust_pressure = (perf.p_exit - p_ambient) * perf.exit_area;
  perf.thrust_ideal = perf.thrust_momentum + perf.thrust_pressure;
  perf.thrust = perf.lambda_divergence * perf.eta_nozzle * perf.thrust_momentum +
                perf.thrust_pressure;
  perf.c_effective = perf.thrust / perf.mdot;
  perf.isp_ideal = perf.thrust_ideal / (perf.mdot * constants::g0);
  perf.isp = perf.thrust / (perf.mdot * constants::g0);
  perf.cf_ideal = perf.thrust_ideal / (perf.p_chamber * perf.throat_area);
  perf.cf = perf.thrust / (perf.p_chamber * perf.throat_area);
  perf.isp_vacuum = (perf.lambda_divergence * perf.eta_nozzle * perf.mdot * exit_sup.u +
                     exit_sup.gas.p * perf.exit_area) / (perf.mdot * constants::g0);

  // --- ambient-pressure family -------------------------------------------
  // F(p_a) = F_vac - p_a A_e exactly, for a full-flowing nozzle, so every
  // ambient pressure of interest is one subtraction away from the vacuum
  // thrust.  Nothing below re-solves the flow.
  const double thrust_vacuum = perf.isp_vacuum * perf.mdot * constants::g0;
  const auto thrustAt = [&](double p_a) { return thrust_vacuum - p_a * perf.exit_area; };
  perf.thrust_sea_level = thrustAt(constants::atm);
  perf.isp_sea_level = perf.thrust_sea_level / (perf.mdot * constants::g0);

  if (!opts.ascent_altitudes.empty()) {
    if (opts.ascent_weights.size() != opts.ascent_altitudes.size())
      throw ConfigError("ascent profile: " + std::to_string(opts.ascent_altitudes.size()) +
                        " altitudes but " + std::to_string(opts.ascent_weights.size()) +
                        " weights; they must match");
    double w_sum = 0.0, wp_sum = 0.0, wf_sum = 0.0;
    for (std::size_t i = 0; i < opts.ascent_altitudes.size(); ++i) {
      const double w = opts.ascent_weights[i];
      if (w < 0.0) throw ConfigError("ascent profile: weights must not be negative");
      const double p_a = Atmosphere::at(opts.ascent_altitudes[i]).pressure;
      w_sum += w;
      wp_sum += w * p_a;
      wf_sum += w * thrustAt(p_a);
    }
    if (w_sum <= 0.0) throw ConfigError("ascent profile: the weights sum to zero");
    perf.ascent_mean_ambient = wp_sum / w_sum;
    perf.isp_ascent = (wf_sum / w_sum) / (perf.mdot * constants::g0);
  }

  // --- regime -----------------------------------------------------------
  if (perf.shock_in_nozzle) {
    perf.regime = ExpansionRegime::kSubsonicExit;
  } else if (p_ambient <= 0.0) {
    perf.regime = ExpansionRegime::kUnderExpanded;
  } else if (std::abs(perf.p_exit / p_ambient - 1.0) <= opts.ideal_tolerance) {
    perf.regime = ExpansionRegime::kIdeallyExpanded;
  } else if (perf.p_exit > p_ambient) {
    perf.regime = ExpansionRegime::kUnderExpanded;
  } else {
    perf.regime = ExpansionRegime::kOverExpanded;
  }

  // --- empirical separation diagnostic ----------------------------------
  // The exit-plane margin is the optimisable form of the same criterion: the
  // wall pressure falls and the Mach number rises monotonically through the
  // divergent, so p_wall - p_sep(M) is monotone and the exit plane is where it
  // is smallest.  Positive => attached to the lip; negative => separated
  // somewhere inside.  Normalising by p_ambient makes it a percentage of
  // ambient and so comparable across altitudes.
  perf.separation_margin = NozzlePerformance::kNoSeparationRisk;
  if (opts.separation != SeparationCriterion::kNone && p_ambient > 0.0 &&
      !perf.shock_in_nozzle) {
    perf.separation_margin =
        (perf.p_exit - separationPressure(opts.separation, p_ambient, perf.mach_exit)) / p_ambient;
  }

  if (opts.separation != SeparationCriterion::kNone && p_ambient > 0.0 &&
      !perf.shock_in_nozzle && perf.regime == ExpansionRegime::kOverExpanded) {
    // Wall pressure falls and Mach number rises monotonically along the
    // divergent, so the separation margin p_wall - p_sep(M) is monotone and the
    // onset can be bisected instead of marched station by station.
    auto margin = [&](double eps) {
      const auto e = flow.atAreaRatio(eps, true);
      return e.gas.p - separationPressure(opts.separation, p_ambient, e.mach);
    };
    double lo = 1.0 + 1.0e-9, hi = perf.area_ratio;
    bool ok = true;
    double f_lo = 0.0, f_hi = 0.0;
    try {
      f_lo = margin(lo);
      f_hi = margin(hi);
    } catch (const IgnisError&) {
      ok = false;
    }
    if (ok && f_hi < 0.0) {
      if (f_lo < 0.0) {
        // Separated from the throat onwards.
        perf.separation_area_ratio = lo;
      } else {
        for (int i = 0; i < 60; ++i) {
          const double mid = 0.5 * (lo + hi);
          double fm;
          try {
            fm = margin(mid);
          } catch (const IgnisError&) {
            break;
          }
          if (fm > 0.0) lo = mid; else hi = mid;
          if (hi - lo < 1.0e-4 * perf.area_ratio) break;
        }
        perf.separation_area_ratio = hi;
      }
      perf.separation_predicted = true;
      perf.separation_x = geom.exitPosition();
      for (const auto& st : geom.stations())
        if (st.x >= geom.throatPosition() && st.area_ratio <= perf.separation_area_ratio)
          perf.separation_x = st.x;
    }
  }
  return perf;
}

AxialProfile sampleAxialProfile(const NozzleFlow& flow, const NozzleGeometry& geom) {
  AxialProfile prof;
  const double xt = geom.throatPosition();
  for (const auto& st : geom.stations()) {
    ExpansionState e;
    if (std::abs(st.x - xt) < 1e-14 * std::max(1.0, xt)) {
      e = flow.throat();
    } else {
      e = flow.atAreaRatio(std::max(st.area_ratio, 1.0), st.x > xt);
    }
    prof.x.push_back(st.x);
    prof.radius.push_back(st.r);
    prof.area_ratio.push_back(st.area_ratio);
    prof.p.push_back(e.gas.p);
    prof.T.push_back(e.gas.T);
    prof.rho.push_back(e.gas.rho);
    prof.u.push_back(e.u);
    prof.mach.push_back(e.mach);
    prof.a.push_back(e.gas.a);
    prof.gamma_s.push_back(e.gas.gamma_s);
    prof.cp.push_back(e.gas.cp_eff);
    prof.molar_mass.push_back(e.gas.M);
    prof.supersonic.push_back(st.x > xt ? 1 : 0);
  }
  return prof;
}

}  // namespace ignis
