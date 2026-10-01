// SPDX-License-Identifier: MIT
#include "ignis/kinetics/NozzleKinetics.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <iomanip>
#include <sstream>

#include "ignis/core/Constants.hpp"
#include "ignis/core/Exceptions.hpp"

namespace ignis {

FrozenExpansion frozenExpansion(const GasMixture& mix, const Eigen::VectorXd& n, double h0,
                                double t0, double p0, double mdot, double area) {
  const double r_gas = constants::R_universal * n.sum();
  const double s0 = mix.entropy(n, t0, p0);
  // At fixed composition s(T, p) = s(T, p0) - R ln(p/p0), so the isentrope
  // gives p(T) directly; energy gives u(T); mass gives the area.
  auto at = [&](double T, FrozenExpansion& f) {
    f.t = T;
    f.p = p0 * std::exp((mix.entropy(n, T, p0) - s0) / r_gas);
    f.u = std::sqrt(std::max(2.0 * (h0 - mix.enthalpy(n, T)), 0.0));
    f.rho = f.p / (r_gas * T);
    return (f.u > 0.0) ? mdot / (f.rho * f.u) : 1.0e300;
  };
  FrozenExpansion f;
  double hi = t0, lo = mix.database().tMinCommon() + 1.0;
  if (at(hi, f) > area * (1.0 + 1e-12))
    throw ConfigError("frozen expansion: the start state is already past the target area");
  if (at(lo, f) < area)
    throw RangeError("frozen expansion: the gas reaches the species data's lower temperature "
                     "limit before the target area");
  // Supersonic branch: the area grows as the gas cools.
  for (int it = 0; it < 200 && hi - lo > 1.0e-10 * hi; ++it) {
    const double mid = 0.5 * (lo + hi);
    if (at(mid, f) < area) hi = mid; else lo = mid;
  }
  at(0.5 * (lo + hi), f);
  f.isp_vacuum = (f.u + f.p / (f.rho * f.u)) / constants::g0;
  return f;
}

namespace {

/// The kinetic nozzle as an ODE in x, y = [n_k of the mechanism's species,
/// u].  The other species ride along at their start values.
class NozzleOde {
 public:
  NozzleOde(const GasMixture& mix, const Mechanism& mech,
            const std::function<double(double)>& area, double x0, double x1,
            Eigen::VectorXd n_all, double h0, double mdot, double rate_multiplier,
            double t_start)
      : mix_(mix),
        mech_(mech),
        area_(area),
        x0_(x0),
        x1_(x1),
        n_all_(std::move(n_all)),
        h0_(h0),
        mdot_(mdot),
        mult_(rate_multiplier),
        t_hint_(t_start),
        active_(mech.activeSpecies()) {
    dx_ = 1.0e-7 * std::max(x1 - x0, 1.0e-3);
  }

  int size() const { return static_cast<int>(active_.size()) + 1; }
  const std::vector<int>& active() const { return active_; }
  int evaluations() const { return evaluations_; }

  Eigen::VectorXd pack(const Eigen::VectorXd& n, double u) const {
    Eigen::VectorXd y(size());
    for (std::size_t i = 0; i < active_.size(); ++i) y(static_cast<Eigen::Index>(i)) = n(active_[i]);
    y(size() - 1) = u;
    return y;
  }
  void composition(const Eigen::VectorXd& y, Eigen::VectorXd& n) const {
    n = n_all_;
    for (std::size_t i = 0; i < active_.size(); ++i) n(active_[i]) = y(static_cast<Eigen::Index>(i));
  }

  /// T from the energy equation, by Newton from the last temperature found.
  double temperature(const Eigen::VectorXd& n, double u) const {
    const double h = h0_ - 0.5 * u * u;
    double T = t_hint_;
    for (int it = 0; it < 40; ++it) {
      const double dT = (mix_.enthalpy(n, T) - h) / mix_.cpFrozen(n, T);
      T -= dT;
      if (!(T > 0.0)) break;
      if (std::abs(dT) < 1.0e-13 * T) {
        t_hint_ = T;
        return T;
      }
    }
    T = mix_.temperatureFromEnthalpy(n, h, t_hint_);
    t_hint_ = T;
    return T;
  }

  double area(double x) const { return area_(x); }
  double areaSlope(double x) const {
    const double xa = std::max(x - dx_, x0_);
    const double xb = std::min(x + dx_, x1_);
    return (area(xb) - area(xa)) / (xb - xa);
  }

  struct Local {
    double T = 0.0, p = 0.0, rho = 0.0, u = 0.0, mach_frozen = 0.0, area = 0.0;
  };

  void rhs(double x, const Eigen::VectorXd& y, Eigen::VectorXd& f, Local* local = nullptr) const {
    ++evaluations_;
    composition(y, n_);
    const double u = y(size() - 1);
    const double T = temperature(n_, u);
    const double A = area(x);
    const double rho = mdot_ / (u * A);
    const double N = n_.sum();
    const double r_gas = constants::R_universal * N;
    c_ = rho * n_;
    mech_.productionRates(T, c_, wdot_, mult_);
    f.resize(size());
    double sum_dn = 0.0, sum_hdn = 0.0;
    const auto& db = mix_.database();
    for (std::size_t i = 0; i < active_.size(); ++i) {
      const int k = active_[i];
      const double dn = wdot_(k) / (rho * u);
      f(static_cast<Eigen::Index>(i)) = dn;
      sum_dn += dn;
      sum_hdn += db[static_cast<std::size_t>(k)].h(T) * dn;
    }
    const double cp = mix_.cpFrozen(n_, T);
    const double psi = sum_dn / N - sum_hdn / (cp * T);
    const double gamma_f = cp / (cp - r_gas);
    const double m2 = u * u / (gamma_f * r_gas * T);
    f(size() - 1) = u * (-areaSlope(x) / A + psi) / (1.0 - m2);
    if (local != nullptr) {
      local->T = T;
      local->rho = rho;
      local->p = rho * r_gas * T;
      local->u = u;
      local->mach_frozen = std::sqrt(m2);
      local->area = A;
    }
  }

  double dx() const { return dx_; }

 private:
  const GasMixture& mix_;
  const Mechanism& mech_;
  const std::function<double(double)>& area_;
  double x0_, x1_;
  Eigen::VectorXd n_all_;
  double h0_, mdot_, mult_;
  mutable double t_hint_;
  std::vector<int> active_;
  double dx_ = 0.0;
  mutable Eigen::VectorXd n_, c_, wdot_;
  mutable int evaluations_ = 0;
};

double frozenMach(const ExpansionState& s) {
  return s.u / std::sqrt(s.gas.gamma_frozen * s.gas.R * s.gas.T);
}

}  // namespace

KineticNozzleResult marchKinetic(const GasMixture& mix, const Mechanism& mech,
                                 const std::function<double(double)>& area, double x0,
                                 double x_end, const Eigen::VectorXd& n0, double t0, double p0,
                                 double u0, const KineticNozzleOptions& opts) {
  if (&mix.database() != &mech.database())
    throw ConfigError("finite-rate nozzle: the mechanism was read against a different species "
                      "database from the gas's");
  if (!(opts.rate_multiplier >= 0.0))
    throw ConfigError("finite-rate nozzle: the rate multiplier must be >= 0");
  if (!(x_end > x0)) throw ConfigError("finite-rate nozzle: the march must run downstream");
  if (!(t0 > 0.0 && p0 > 0.0 && u0 > 0.0))
    throw ConfigError("finite-rate nozzle: the start state needs T, p and u > 0");
  const auto& db = mix.database();
  KineticNozzleResult res;
  res.start_x = x0;
  const double h0 = mix.enthalpy(n0, t0) + 0.5 * u0 * u0;
  const double rho0 = p0 / (constants::R_universal * n0.sum() * t0);
  res.mass_flow = rho0 * u0 * area(x0);
  NozzleOde ode(mix, mech, area, x0, x_end, n0, h0, res.mass_flow, opts.rate_multiplier, t0);
  const int m = ode.size();
  Eigen::VectorXd y = ode.pack(n0, u0);
  const Eigen::VectorXd b0 = db.elementMatrix() * n0;

  // Error scales: species absolutely and relatively, velocity relatively.
  Eigen::VectorXd atol = Eigen::VectorXd::Constant(m, opts.absolute_tolerance);
  atol(m - 1) = 1.0e-6;
  const double rtol = opts.relative_tolerance;

  std::vector<double> rx, rT, rp, ru, rm, ra;
  std::vector<std::vector<double>> rX(ode.active().size());
  Eigen::VectorXd n_now;
  auto record = [&](double x, const Eigen::VectorXd& yy, const NozzleOde::Local& loc) {
    rx.push_back(x);
    rT.push_back(loc.T);
    rp.push_back(loc.p);
    ru.push_back(loc.u);
    rm.push_back(loc.mach_frozen);
    ra.push_back(loc.area);
    ode.composition(yy, n_now);
    const double N = n_now.sum();
    for (std::size_t i = 0; i < ode.active().size(); ++i)
      rX[i].push_back(std::max(n_now(ode.active()[i]), 0.0) / N);
    const double h = mix.enthalpy(n_now, loc.T) + 0.5 * loc.u * loc.u;
    res.energy_residual = std::max(res.energy_residual, std::abs(h - h0) / std::abs(h0));
  };

  Eigen::VectorXd f0, f1, fa, fx;
  NozzleOde::Local loc;
  double x = x0;
  ode.rhs(x, y, f0, &loc);
  if (!(loc.mach_frozen > 1.0))
    throw ConfigError("finite-rate nozzle: the march must start supersonic in the frozen sense "
                      "(its momentum equation is singular at frozen Mach 1)");
  res.start_frozen_mach = loc.mach_frozen;
  record(x, y, loc);
  Eigen::MatrixXd J(m, m);
  const double length = x_end - x0;
  double h = 1.0e-6 * length;
  const Eigen::MatrixXd I = Eigen::MatrixXd::Identity(m, m);
  while (x < x_end - 1.0e-12 * length) {
    if (++res.steps > opts.max_steps)
      throw ConvergenceError("finite-rate nozzle: more than the allowed number of steps");
    h = std::min(h, x_end - x);
    // Jacobian (finite differences) and the explicit x-dependence.
    for (int j = 0; j < m; ++j) {
      Eigen::VectorXd yp = y;
      const double scale = (j == m - 1) ? std::max(std::abs(y(j)), 1.0)
                                        : std::max(std::abs(y(j)), 1.0e-6);
      const double d = 1.0e-7 * scale;
      yp(j) += d;
      ode.rhs(x, yp, f1);
      J.col(j) = (f1 - f0) / d;
    }
    {
      const double dxx = std::min(ode.dx(), 0.5 * (x_end - x));
      if (dxx > 0.0) {
        ode.rhs(x + dxx, y, f1);
        fx = (f1 - f0) / dxx;
      } else {
        fx.setZero(m);
      }
    }
    // One full step and two half steps of linearly implicit Euler.
    const Eigen::PartialPivLU<Eigen::MatrixXd> full(I - h * J);
    const Eigen::VectorXd y1 = y + full.solve(h * f0 + h * h * fx);
    const double hh = 0.5 * h;
    const Eigen::PartialPivLU<Eigen::MatrixXd> half(I - hh * J);
    const Eigen::VectorXd ya = y + half.solve(hh * f0 + hh * hh * fx);
    bool ok = ya.allFinite() && ya(m - 1) > 0.0;
    Eigen::VectorXd y2;
    if (ok) {
      ode.rhs(x + hh, ya, fa);
      y2 = ya + half.solve(hh * fa + hh * hh * fx);
      ok = y2.allFinite() && y1.allFinite() && y2(m - 1) > 0.0;
    }
    double err = 1.0e10;
    if (ok) {
      double sum = 0.0;
      for (int i = 0; i < m; ++i) {
        const double sc = atol(i) + rtol * std::max(std::abs(y(i)), std::abs(y2(i)));
        const double e = (y2(i) - y1(i)) / sc;
        sum += e * e;
      }
      err = std::sqrt(sum / m);
    }
    if (ok && err <= 1.0) {
      y = 2.0 * y2 - y1;  // Richardson: second order, still L-stable
      x += h;
      if (x > x_end - 1.0e-12 * length) x = x_end;
      ode.rhs(x, y, f0, &loc);
      record(x, y, loc);
    } else {
      ++res.rejected_steps;
    }
    const double grow = ok ? 0.9 / std::sqrt(std::max(err, 1.0e-10)) : 0.2;
    h *= std::min(4.0, std::max(0.2, grow));
    if (h < 1.0e-14 * length)
      throw ConvergenceError("finite-rate nozzle: the step size underflowed");
  }
  res.rhs_evaluations = ode.evaluations();

  ode.composition(y, n_now);
  res.n_exit = n_now;
  res.t_exit = loc.T;
  res.p_exit = loc.p;
  res.rho_exit = loc.rho;
  res.u_exit = loc.u;
  res.mach_frozen_exit = loc.mach_frozen;
  res.isp_vacuum = (loc.u + loc.p / (loc.rho * loc.u)) / constants::g0;
  {
    const Eigen::VectorXd b1 = db.elementMatrix() * n_now;
    res.element_residual = (b1 - b0).cwiseAbs().maxCoeff() / b0.cwiseAbs().maxCoeff();
  }

  // The march, resampled evenly in x.
  const int np = std::max(2, opts.profile_points);
  for (const int k : ode.active()) res.species.push_back(db[static_cast<std::size_t>(k)].name());
  res.mole_fraction.assign(ode.active().size(), {});
  const double a_start = area(x0);
  std::size_t j = 0;
  for (int i = 0; i < np; ++i) {
    const double xi = x0 + (x_end - x0) * i / (np - 1);
    while (j + 2 < rx.size() && rx[j + 1] < xi) ++j;
    const std::size_t j1 = std::min(j + 1, rx.size() - 1);
    const double w =
        (rx[j1] > rx[j]) ? std::clamp((xi - rx[j]) / (rx[j1] - rx[j]), 0.0, 1.0) : 0.0;
    auto lerp = [&](const std::vector<double>& v) { return v[j] + w * (v[j1] - v[j]); };
    res.x.push_back(xi);
    res.area_ratio.push_back(lerp(ra) / a_start);
    res.temperature.push_back(lerp(rT));
    res.pressure.push_back(lerp(rp));
    res.velocity.push_back(lerp(ru));
    res.mach_frozen.push_back(lerp(rm));
    for (std::size_t s2 = 0; s2 < rX.size(); ++s2) res.mole_fraction[s2].push_back(lerp(rX[s2]));
  }
  if (!mech.dropped().empty()) {
    std::ostringstream os;
    os << mech.dropped().size() << " of the mechanism's reactions involve species this gas "
       << "does not contain and were left out";
    res.warnings.push_back(os.str());
  }
  return res;
}

KineticNozzleResult integrateKineticNozzle(const NozzleFlow& flow, const NozzleGeometry& geom,
                                           const Mechanism& mech,
                                           const KineticNozzleOptions& opts) {
  if (flow.model() != CompositionModel::kEquilibrium)
    throw ConfigError("finite-rate nozzle: it starts from shifting equilibrium, so the chamber's "
                      "composition model must be 'equilibrium'");
  if (!(opts.start_frozen_mach > 1.0))
    throw ConfigError("finite-rate nozzle: the start frozen Mach number must exceed 1");
  const auto& mix = flow.mixture();
  const double eps_exit = geom.expansionRatio();
  const double at = geom.throatArea();

  // --- the start station: frozen Mach number reaches the target ----------
  double lo = 1.0001, hi = eps_exit;
  if (frozenMach(flow.atAreaRatio(hi, true)) <= opts.start_frozen_mach)
    throw ConfigError("finite-rate nozzle: the frozen Mach number never reaches the start value "
                      "inside this nozzle; it is too short to integrate");
  if (frozenMach(flow.atAreaRatio(lo, true)) >= opts.start_frozen_mach) {
    hi = lo;
  } else {
    for (int it = 0; it < 100 && hi - lo > 1e-12 * hi; ++it) {
      const double mid = 0.5 * (lo + hi);
      if (frozenMach(flow.atAreaRatio(mid, true)) < opts.start_frozen_mach) lo = mid; else hi = mid;
    }
  }
  const double eps_start = hi;
  const ExpansionState start = flow.atAreaRatio(eps_start, true);
  double xa = geom.throatPosition(), xb = geom.exitPosition();
  for (int it = 0; it < 200 && xb - xa > 1e-13 * geom.exitPosition(); ++it) {
    const double mid = 0.5 * (xa + xb);
    const double r = geom.radius(mid);
    if (constants::pi * r * r / at < eps_start) xa = mid; else xb = mid;
  }
  const std::function<double(double)> area = [&geom](double x) {
    const double r = geom.radius(x);
    return constants::pi * r * r;
  };
  KineticNozzleResult res = marchKinetic(mix, mech, area, xb, geom.exitPosition(), start.gas.n,
                                         start.gas.T, start.gas.p, start.u, opts);
  // The march's area ratios are relative to its start; the nozzle's to the throat.
  for (auto& e : res.area_ratio) e *= eps_start;
  res.start_area_ratio = eps_start;
  res.exit_area_ratio = eps_exit;

  // --- the two limits from the same geometry -----------------------------
  const ExpansionState shift = flow.atAreaRatio(eps_exit, true);
  res.isp_vacuum_shifting =
      (shift.u + shift.gas.p / (shift.gas.rho * shift.u)) / constants::g0;
  const double h0 = start.gas.h + 0.5 * start.u * start.u;
  const FrozenExpansion frozen = frozenExpansion(mix, start.gas.n, h0, start.gas.T, start.gas.p,
                                                 res.mass_flow, geom.exitArea());
  res.isp_vacuum_frozen = frozen.isp_vacuum;
  res.kinetic_efficiency = res.isp_vacuum / res.isp_vacuum_shifting;
  const double span = res.isp_vacuum_shifting - res.isp_vacuum_frozen;
  res.recovered_fraction =
      (std::abs(span) > 1e-9) ? (res.isp_vacuum - res.isp_vacuum_frozen) / span : 1.0;
  return res;
}

std::string KineticNozzleResult::summary() const {
  std::ostringstream os;
  os << std::fixed;
  os << "finite-rate nozzle (one-dimensional kinetics)\n"
     << "  start                 A/At " << std::setprecision(4) << start_area_ratio
     << " (frozen Mach " << std::setprecision(3) << start_frozen_mach << "), x = "
     << std::setprecision(1) << start_x * 1e3 << " mm, from shifting equilibrium\n"
     << "  exit                  " << std::setprecision(1) << t_exit << " K, "
     << std::setprecision(1) << p_exit << " Pa, " << u_exit << " m/s\n"
     << "  vacuum Isp (inviscid) " << std::setprecision(2) << isp_vacuum << " s finite rate; "
     << isp_vacuum_shifting << " s shifting; " << isp_vacuum_frozen << " s frozen from the start\n"
     << "  kinetic efficiency    " << std::setprecision(5) << kinetic_efficiency
     << " (recovers " << std::setprecision(1) << 100.0 * recovered_fraction
     << " % of the shifting-over-frozen impulse)\n"
     << "  integration           " << steps << " steps (" << rejected_steps << " rejected), "
     << rhs_evaluations << " evaluations; residuals: energy " << std::scientific
     << std::setprecision(1) << energy_residual << ", elements " << element_residual << "\n";
  for (const auto& w : warnings) os << "  note: " << w << "\n";
  return os.str();
}

}  // namespace ignis
