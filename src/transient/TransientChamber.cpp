// SPDX-License-Identifier: MIT
#include "ignis/transient/TransientChamber.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <sstream>

#include "ignis/core/Constants.hpp"
#include "ignis/core/Exceptions.hpp"

namespace ignis {

Schedule::Schedule(std::vector<double> times, std::vector<double> values)
    : t_(std::move(times)), v_(std::move(values)) {
  if (t_.size() != v_.size() || t_.empty())
    throw ConfigError("schedule: times and values must be non-empty and the same length");
  for (std::size_t i = 1; i < t_.size(); ++i)
    if (!(t_[i] > t_[i - 1]))
      throw ConfigError("schedule: times must be strictly increasing");
}

Schedule Schedule::ramp(double t_open, double ramp_up, double steady_value, double t_close,
                        double ramp_down) {
  if (!(ramp_up > 0.0) || !(ramp_down > 0.0))
    throw ConfigError("schedule: ramp durations must be positive");
  if (!(t_close >= t_open + ramp_up))
    throw ConfigError("schedule: the close command must come after the opening ramp finishes");
  return Schedule({t_open, t_open + ramp_up, t_close, t_close + ramp_down},
                  {0.0, steady_value, steady_value, 0.0});
}

Schedule Schedule::constant(double value) { return Schedule({0.0, 1.0e9}, {value, value}); }

double Schedule::at(double t) const {
  if (t_.empty()) return 0.0;
  if (t <= t_.front()) return v_.front();
  if (t >= t_.back()) return v_.back();
  const auto it = std::upper_bound(t_.begin(), t_.end(), t);
  const std::size_t i = static_cast<std::size_t>(std::distance(t_.begin(), it)) - 1;
  const double f = (t - t_[i]) / (t_[i + 1] - t_[i]);
  return (1.0 - f) * v_[i] + f * v_[i + 1];
}

std::string TransientResult::summary() const {
  std::ostringstream os;
  os << std::fixed << std::setprecision(4);
  os << "transient (" << integrator << ")\n"
     << "  samples                " << samples.size() << "\n"
     << "  accepted / rejected    " << steps << " / " << rejected_steps << "\n"
     << "  right-hand sides       " << rhs_evaluations << "\n"
     << "  peak chamber pressure  " << max_pressure * 1e-6 << " MPa\n";
  if (time_to_90_percent >= 0.0)
    os << "  time to 90 % of peak   " << time_to_90_percent * 1e3 << " ms\n";
  os << std::fixed << std::setprecision(3)
     << "  chamber O/F range      " << min_mixture_ratio << " to " << max_mixture_ratio << "\n"
     << "  chamber T range        " << std::setprecision(1) << min_temperature << " to "
     << max_temperature << " K\n"
     << std::scientific << std::setprecision(3)
     << "  mass conservation err  " << mass_conservation_error << "\n"
     << "  energy conservation err " << energy_conservation_error;
  return os.str();
}

namespace {

/// The integration state: oxidiser mass, fuel mass, internal energy.
/// Integration state.  `Im` and `Ie` accumulate the net mass and energy flux
/// with exactly the same Runge-Kutta weights as the physical states, so
/// comparing them against the change in (m, U) is a direct check that the
/// conservative form of the equations survived the implementation.
struct State {
  double m_ox = 0.0, m_f = 0.0, U = 0.0, Im = 0.0, Ie = 0.0;
  State operator+(const State& o) const {
    return {m_ox + o.m_ox, m_f + o.m_f, U + o.U, Im + o.Im, Ie + o.Ie};
  }
  State operator*(double s) const { return {m_ox * s, m_f * s, U * s, Im * s, Ie * s}; }
};

/// Thermodynamic closure of a state.
/// Chemical energy per second that the incoming streams *could* release.
///
/// Only the part of the stream that has a reaction partner counts: pure
/// oxidiser or pure fuel carries no heat release, and a stream outside the
/// tabulated mixture-ratio window is split into a reactive core at the nearest
/// tabulated ratio plus an inert remainder.  Withholding (1 - eta) of this is
/// what the ignition ramp does; with eta = 0 the injected stream is
/// energetically equivalent to injecting combustion products at 298.15 K,
/// which is a well-defined and bounded no-heat-release limit.
double withheldHeatRate(const EquilibriumTable& tab, double mdot_ox, double mdot_f) {
  if (!(mdot_ox > 0.0) || !(mdot_f > 0.0)) return 0.0;
  const double mr_in = mdot_ox / mdot_f;
  if (mr_in > tab.mrMax()) {
    const double mr = tab.mrMax();
    return mdot_f * (1.0 + mr) * tab.heatOfCombustion(mr);
  }
  if (mr_in < tab.mrMin()) {
    const double mr = tab.mrMin();
    return mdot_ox * (1.0 + 1.0 / mr) * tab.heatOfCombustion(mr);
  }
  return (mdot_ox + mdot_f) * tab.heatOfCombustion(mr_in);
}

struct Closure {
  double m = 0.0, rho = 0.0, u = 0.0, T = 0.0, p = 0.0, mr = 0.0, mr_raw = 0.0;
  bool mr_clamped = false;
  double M = 0.0, gamma = 0.0, c_star = 0.0, h = 0.0;
  double mdot_out = 0.0;
  bool choked = false;
};

Closure close(const EquilibriumTable& tab, const TransientSpec& spec, const State& y) {
  Closure c;
  c.m = y.m_ox + y.m_f;
  if (!(c.m > spec.mass_floor)) {
    std::ostringstream os;
    os << "transient: chamber mass fell to " << c.m << " kg, below the " << spec.mass_floor
       << " kg floor. Give the chamber a finite initial charge or shorten the time step.";
    throw InfeasibleError(os.str());
  }
  if (!(y.m_f > 0.0) || !(y.m_ox > 0.0))
    throw InfeasibleError("transient: the chamber ran out of one propellant component, so the "
                          "mixture ratio is undefined. Keep a small residual charge of both.");
  c.mr_raw = y.m_ox / y.m_f;
  if (c.mr_raw < tab.mrMin() || c.mr_raw > tab.mrMax()) {
    // Clamping here would silently pair the chamber's internal energy with the
    // wrong composition, so the run stops with an actionable message instead.
    std::ostringstream os;
    os << "transient: the chamber mixture ratio reached " << c.mr_raw
       << ", outside the tabulated range [" << tab.mrMin() << ", " << tab.mrMax()
       << "]. Widen the equilibrium table's mixture-ratio axis, or use a valve schedule "
          "that does not drive the chamber that far off ratio (a fuel-lead start and a "
          "fuel-lag shutdown keep it fuel-rich).";
    throw InfeasibleError(os.str());
  }
  c.mr = c.mr_raw;
  c.rho = c.m / spec.chamber_volume;
  c.u = y.U / c.m;
  tab.solveState(c.mr, c.u, c.rho, c.T, c.p);
  c.M = tab.molarMass(c.mr, c.T, c.p);
  c.gamma = tab.gammaS(c.mr, c.T, c.p);
  // Use the *actual* chamber temperature, not the adiabatic flame value: during
  // a start-up the chamber is far from equilibrium temperature and the outflow
  // must reflect that.
  c.c_star = tab.cStarAt(c.mr, c.T, c.p);
  c.h = c.u + c.p / c.rho;

  // Nozzle outflow.
  const double pa = spec.ambient_pressure;
  const double g = c.gamma;
  const double crit = std::pow(2.0 / (g + 1.0), g / (g - 1.0));
  if (c.p <= pa) {
    c.mdot_out = 0.0;
    c.choked = false;
  } else if (pa / c.p <= crit) {
    c.mdot_out = c.p * spec.throat_area / c.c_star;
    c.choked = true;
  } else {
    const double R = constants::R_universal / c.M;
    const double pr = pa / c.p;
    const double term = std::pow(pr, 2.0 / g) - std::pow(pr, (g + 1.0) / g);
    c.mdot_out = spec.throat_area * c.p *
                 std::sqrt(std::max(0.0, 2.0 * g / ((g - 1.0) * R * c.T) * term));
    c.choked = false;
  }
  return c;
}

State rhs(const EquilibriumTable& tab, const TransientSpec& spec, double t, const State& y,
          Closure* out) {
  Closure c;
  try {
    c = close(tab, spec, y);
  } catch (const IgnisError& e) {
    std::ostringstream os;
    os << "transient failed at t = " << t << " s (chamber mass " << (y.m_ox + y.m_f)
       << " kg, O/F " << (y.m_f > 0.0 ? y.m_ox / y.m_f : -1.0) << "): " << e.what();
    throw InfeasibleError(os.str());
  }
  if (out != nullptr) *out = c;
  const double mdot_ox = std::max(0.0, spec.oxidizer_flow.at(t));
  const double mdot_f = std::max(0.0, spec.fuel_flow.at(t));
  const double eta = std::min(1.0, std::max(0.0, spec.combustion_efficiency.at(t)));
  State d;
  d.m_ox = mdot_ox - c.mdot_out * (y.m_ox / c.m);
  d.m_f = mdot_f - c.mdot_out * (y.m_f / c.m);

  d.U = mdot_ox * spec.oxidizer_inlet_enthalpy + mdot_f * spec.fuel_inlet_enthalpy -
        c.mdot_out * c.h - (1.0 - eta) * withheldHeatRate(tab, mdot_ox, mdot_f) -
        spec.wall_heat_rate;
  d.Im = d.m_ox + d.m_f;
  d.Ie = d.U;
  return d;
}

}  // namespace

TransientResult simulateTransient(const EquilibriumTable& table, const TransientSpec& spec) {
  if (!(spec.chamber_volume > 0.0)) throw ConfigError("transient: chamber_volume must be positive");
  if (!(spec.throat_area > 0.0)) throw ConfigError("transient: throat_area must be positive");
  if (!(spec.t_end > 0.0)) throw ConfigError("transient: t_end must be positive");
  if (!(spec.dt > 0.0)) throw ConfigError("transient: dt must be positive");
  if (!(spec.dt_output > 0.0)) throw ConfigError("transient: dt_output must be positive");
  if (spec.integrator != "rk4" && spec.integrator != "rk45")
    throw ConfigError("transient: integrator must be 'rk4' or 'rk45', got '" + spec.integrator + "'");

  // --- initial charge ---------------------------------------------------
  double mr0 = spec.initial_mixture_ratio;
  if (!(mr0 > 0.0)) mr0 = 0.5 * (table.mrMin() + table.mrMax());
  mr0 = std::min(std::max(mr0, table.mrMin()), table.mrMax());
  const double M0 = table.molarMass(mr0, spec.initial_temperature, spec.initial_pressure);
  const double rho0 = spec.initial_pressure * M0 /
                      (constants::R_universal * spec.initial_temperature);
  const double m0 = rho0 * spec.chamber_volume;
  if (!(m0 > spec.mass_floor))
    throw ConfigError("transient: the initial charge is below the mass floor; raise "
                      "initial_pressure or the chamber volume");
  const double u0 = table.internalEnergy(mr0, spec.initial_temperature, spec.initial_pressure);

  State y;
  y.m_f = m0 / (1.0 + mr0);
  y.m_ox = m0 - y.m_f;
  y.U = m0 * u0;

  TransientResult res;
  res.integrator = spec.integrator;

  // Cash-Karp embedded RK4(5) coefficients.
  static const double a[6] = {0.0, 1.0 / 5, 3.0 / 10, 3.0 / 5, 1.0, 7.0 / 8};
  static const double b[6][5] = {
      {0, 0, 0, 0, 0},
      {1.0 / 5, 0, 0, 0, 0},
      {3.0 / 40, 9.0 / 40, 0, 0, 0},
      {3.0 / 10, -9.0 / 10, 6.0 / 5, 0, 0},
      {-11.0 / 54, 5.0 / 2, -70.0 / 27, 35.0 / 27, 0},
      {1631.0 / 55296, 175.0 / 512, 575.0 / 13824, 44275.0 / 110592, 253.0 / 4096}};
  static const double c5[6] = {37.0 / 378, 0, 250.0 / 621, 125.0 / 594, 0, 512.0 / 1771};
  static const double c4[6] = {2825.0 / 27648, 0, 18575.0 / 48384, 13525.0 / 55296,
                               277.0 / 14336, 1.0 / 4};

  double t = 0.0, h = spec.dt;
  double next_output = 0.0;
  const double m_start = y.m_ox + y.m_f;
  const double U_start = y.U;

  auto record = [&](double time, const State& s) {
    Closure c;
    rhs(table, spec, time, s, &c);
    TransientSample smp;
    smp.t = time;
    smp.mass = c.m;
    smp.mass_oxidizer = s.m_ox;
    smp.mass_fuel = s.m_f;
    smp.internal_energy = s.U;
    smp.pressure = c.p;
    smp.temperature = c.T;
    smp.density = c.rho;
    smp.mixture_ratio = c.mr;
    smp.molar_mass = c.M;
    smp.gamma_s = c.gamma;
    smp.c_star = c.c_star;
    smp.mdot_out = c.mdot_out;
    smp.choked = c.choked;
    const double mr_in = (smp.mdot_fuel_in > 0.0 && smp.mdot_ox_in > 0.0)
                             ? smp.mdot_ox_in / smp.mdot_fuel_in
                             : -1.0;
    smp.inlet_mixture_ratio_clamped =
        (mr_in > 0.0) && (mr_in < table.mrMin() || mr_in > table.mrMax());
    if (smp.inlet_mixture_ratio_clamped) ++res.inlet_clamped_samples;
    smp.mdot_ox_in = std::max(0.0, spec.oxidizer_flow.at(time));
    smp.mdot_fuel_in = std::max(0.0, spec.fuel_flow.at(time));
    smp.eta_heat = std::min(1.0, std::max(0.0, spec.combustion_efficiency.at(time)));
    // Quasi-steady thrust estimate from the instantaneous chamber state.
    const double g = c.gamma;
    const double eps = (spec.exit_area > 0.0) ? spec.exit_area / spec.throat_area : 0.0;
    if (c.choked && eps > 1.0) {
      // Solve the constant-gamma area relation for the supersonic exit Mach.
      double Me = 2.0;
      for (int i = 0; i < 200; ++i) {
        const double f = (1.0 / Me) * std::pow((2.0 / (g + 1.0)) *
                                               (1.0 + 0.5 * (g - 1.0) * Me * Me),
                                               (g + 1.0) / (2.0 * (g - 1.0))) - eps;
        const double d = 1.0e-6 * Me;
        const double f2 = (1.0 / (Me + d)) * std::pow((2.0 / (g + 1.0)) *
                              (1.0 + 0.5 * (g - 1.0) * (Me + d) * (Me + d)),
                              (g + 1.0) / (2.0 * (g - 1.0))) - eps;
        const double dm = -f * d / (f2 - f);
        Me += std::max(-0.5 * Me, std::min(0.5 * Me, dm));
        if (std::abs(dm) < 1e-12) break;
      }
      const double pe = c.p * std::pow(1.0 + 0.5 * (g - 1.0) * Me * Me, -g / (g - 1.0));
      const double Te = c.T / (1.0 + 0.5 * (g - 1.0) * Me * Me);
      const double R = constants::R_universal / c.M;
      const double ue = Me * std::sqrt(g * R * Te);
      smp.thrust = c.mdot_out * ue + (pe - spec.ambient_pressure) * spec.exit_area;
    } else {
      smp.thrust = 0.0;
    }
    res.max_pressure = std::max(res.max_pressure, c.p);
    if (res.samples.empty()) {
      res.min_mixture_ratio = res.max_mixture_ratio = c.mr;
      res.min_temperature = res.max_temperature = c.T;
    } else {
      res.min_mixture_ratio = std::min(res.min_mixture_ratio, c.mr);
      res.max_mixture_ratio = std::max(res.max_mixture_ratio, c.mr);
      res.min_temperature = std::min(res.min_temperature, c.T);
      res.max_temperature = std::max(res.max_temperature, c.T);
    }
    res.samples.push_back(smp);
  };

  record(0.0, y);
  next_output = spec.dt_output;

  const int max_steps = 20000000;
  while (t < spec.t_end && res.steps < max_steps) {
    if (spec.integrator == "rk4") h = spec.dt;
    h = std::min(h, spec.t_end - t);
    h = std::min(h, std::max(spec.dt_min, next_output - t));

    std::array<State, 6> k{};
    State accepted;
    bool ok = true;
    double err_norm = 0.0;

    if (spec.integrator == "rk4") {
      const State k1 = rhs(table, spec, t, y, nullptr);
      const State k2 = rhs(table, spec, t + 0.5 * h, y + k1 * (0.5 * h), nullptr);
      const State k3 = rhs(table, spec, t + 0.5 * h, y + k2 * (0.5 * h), nullptr);
      const State k4 = rhs(table, spec, t + h, y + k3 * h, nullptr);
      accepted = y + (k1 + k2 * 2.0 + k3 * 2.0 + k4) * (h / 6.0);
      res.rhs_evaluations += 4;
    } else {
      for (int s = 0; s < 6; ++s) {
        State ys = y;
        for (int j = 0; j < s; ++j) ys = ys + k[static_cast<std::size_t>(j)] * (h * b[s][j]);
        k[static_cast<std::size_t>(s)] = rhs(table, spec, t + a[s] * h, ys, nullptr);
      }
      res.rhs_evaluations += 6;
      State y5 = y, y4 = y;
      for (int s = 0; s < 6; ++s) {
        y5 = y5 + k[static_cast<std::size_t>(s)] * (h * c5[s]);
        y4 = y4 + k[static_cast<std::size_t>(s)] * (h * c4[s]);
      }
      const double sc_ox = spec.atol + spec.rtol * std::max(std::abs(y.m_ox), std::abs(y5.m_ox));
      const double sc_f = spec.atol + spec.rtol * std::max(std::abs(y.m_f), std::abs(y5.m_f));
      const double sc_U = spec.atol + spec.rtol * std::max(std::abs(y.U), std::abs(y5.U));
      const double e1 = (y5.m_ox - y4.m_ox) / sc_ox;
      const double e2 = (y5.m_f - y4.m_f) / sc_f;
      const double e3 = (y5.U - y4.U) / sc_U;
      err_norm = std::sqrt((e1 * e1 + e2 * e2 + e3 * e3) / 3.0);
      ok = (err_norm <= 1.0) || (h <= spec.dt_min * 1.0000001);
      accepted = y5;
    }

    if (!ok) {
      ++res.rejected_steps;
      const double factor = std::max(0.2, 0.9 * std::pow(1.0 / std::max(err_norm, 1e-16), 0.25));
      h = std::max(spec.dt_min, h * factor);
      if (res.rejected_steps > 2000000)
        throw ConvergenceError("transient: the adaptive integrator kept rejecting steps");
      continue;
    }

    y = accepted;
    t += h;
    ++res.steps;

    if (spec.integrator == "rk45") {
      const double factor = std::min(5.0, 0.9 * std::pow(1.0 / std::max(err_norm, 1e-16), 0.2));
      h = std::min(spec.dt_max, std::max(spec.dt_min, h * std::max(0.2, factor)));
    }

    if (t >= next_output - 1e-15) {
      record(t, y);
      next_output += spec.dt_output;
      while (next_output < t) next_output += spec.dt_output;
    }
  }
  if (res.steps >= max_steps)
    throw ConvergenceError("transient: step limit reached before t_end");
  if (res.samples.empty() || res.samples.back().t < t) record(t, y);

  const double m_end = y.m_ox + y.m_f;
  const double scale_m = std::max(m_end, std::max(m_start, 1e-12));
  res.mass_conservation_error = std::abs((m_end - m_start) - y.Im) / scale_m;
  const double scale_U = std::max(std::abs(y.U), std::max(std::abs(U_start), 1.0));
  res.energy_conservation_error = std::abs((y.U - U_start) - y.Ie) / scale_U;

  // Time to reach 90 % of the peak chamber pressure.
  const double target = 0.9 * res.max_pressure;
  for (const auto& s : res.samples)
    if (s.pressure >= target) { res.time_to_90_percent = s.t; break; }
  return res;
}

}  // namespace ignis
