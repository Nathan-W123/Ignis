// SPDX-License-Identifier: MIT
#include "ignis/cycle/FeedSystem.hpp"

#include <cmath>
#include <iomanip>
#include <sstream>

#include "ignis/core/Constants.hpp"
#include "ignis/core/Exceptions.hpp"

namespace ignis {
namespace {

double colebrook(double reynolds, double relative_roughness) {
  if (reynolds < 1.0) return 64.0;
  if (reynolds < 2300.0) return 64.0 / reynolds;
  double inv = -2.0 * std::log10(relative_roughness / 3.7 + 5.74 / std::pow(reynolds, 0.9));
  for (int i = 0; i < 100; ++i) {
    const double next = -2.0 * std::log10(relative_roughness / 3.7 + 2.51 * inv / reynolds);
    if (std::abs(next - inv) < 1e-12) { inv = next; break; }
    inv = next;
  }
  return 1.0 / (inv * inv);
}

void validate(const FeedLeg& leg, const char* name) {
  if (!(leg.density > 0.0)) throw ConfigError(std::string("feed system: ") + name + " density must be positive");
  if (!(leg.line_diameter > 0.0)) throw ConfigError(std::string("feed system: ") + name + " line diameter must be positive");
  if (!(leg.line_length >= 0.0)) throw ConfigError(std::string("feed system: ") + name + " line length must be non-negative");
  if (!(leg.injector_cd > 0.0 && leg.injector_cd <= 1.0))
    throw ConfigError(std::string("feed system: ") + name + " injector Cd must lie in (0, 1]");
}

/// Line losses for a given mass flow.
void lineLosses(const FeedLeg& leg, double mdot, FeedLegResult& out) {
  const double area = constants::pi * leg.line_diameter * leg.line_diameter / 4.0;
  out.line_velocity = mdot / (leg.density * area);
  out.reynolds = leg.density * out.line_velocity * leg.line_diameter / std::max(leg.viscosity, 1e-12);
  out.friction_factor = colebrook(out.reynolds, leg.line_roughness / leg.line_diameter);
  const double head = 0.5 * leg.density * out.line_velocity * out.line_velocity;
  out.dp_lines = (out.friction_factor * leg.line_length / leg.line_diameter + leg.fitting_k) * head;
  out.dp_dynamic = head;
}

}  // namespace

std::string FeedSystemResult::summary() const {
  std::ostringstream os;
  os << std::fixed << std::setprecision(4);
  auto leg = [&](const char* name, const FeedLegResult& r) {
    os << "  " << name << "\n"
       << "    mass flow          " << r.mass_flow << " kg/s\n"
       << "    line velocity      " << std::setprecision(2) << r.line_velocity << " m/s (Re "
       << std::scientific << std::setprecision(3) << r.reynolds << std::fixed
       << std::setprecision(2) << ", f " << std::setprecision(5) << r.friction_factor << ")\n"
       << "    injector velocity  " << std::setprecision(2) << r.injector_velocity << " m/s\n"
       << "    dp lines           " << std::setprecision(4) << r.dp_lines * 1e-6 << " MPa\n"
       << "    dp dynamic         " << r.dp_dynamic * 1e-6 << " MPa\n"
       << "    dp injector        " << r.dp_injector * 1e-6 << " MPa\n"
       << "    injector area      " << std::setprecision(3) << r.injector_area * 1e6 << " mm^2\n"
       << "    tank pressure      " << std::setprecision(4) << r.required_tank_pressure * 1e-6
       << " MPa\n"
       << "    stiffness          " << std::setprecision(4) << r.stiffness << "\n"
       << "    pump power (equiv) " << std::setprecision(2) << r.pump_power * 1e-3 << " kW\n";
  };
  os << "feed system (pressure-fed)\n"
     << "  chamber pressure     " << chamber_pressure * 1e-6 << " MPa\n"
     << "  mixture ratio        " << mixture_ratio << "\n"
     << "  total mass flow      " << total_mass_flow << " kg/s\n";
  leg("oxidizer", oxidizer);
  leg("fuel", fuel);
  os << "  total pump power     " << std::setprecision(2) << total_pump_power * 1e-3 << " kW";
  if (!stiffness_ok) os << "\n  WARNING: injector stiffness below the configured minimum";
  if (!notes.empty()) os << "\n  note: " << notes;
  return os.str();
}

FeedSystemResult sizeFeedSystem(const FeedSystemSpec& spec, double p_chamber,
                                double target_stiffness) {
  if (!(p_chamber > 0.0)) throw ConfigError("feed system: chamber pressure must be positive");
  if (!(target_stiffness > 0.0 && target_stiffness < 1.0))
    throw ConfigError("feed system: target stiffness must lie in (0, 1)");
  validate(spec.oxidizer, "oxidizer");
  validate(spec.fuel, "fuel");

  FeedSystemResult res;
  res.chamber_pressure = p_chamber;
  auto doLeg = [&](const FeedLeg& leg, FeedLegResult& out) {
    if (!(leg.mass_flow > 0.0)) throw ConfigError("feed system: leg mass flow must be positive");
    out.mass_flow = leg.mass_flow;
    lineLosses(leg, leg.mass_flow, out);
    out.dp_injector = target_stiffness * p_chamber;
    out.injector_velocity = std::sqrt(2.0 * out.dp_injector / leg.density);
    out.injector_area = leg.mass_flow / (leg.injector_cd * leg.density * out.injector_velocity);
    out.required_tank_pressure = p_chamber + out.dp_injector + out.dp_lines + out.dp_dynamic;
    out.stiffness = out.dp_injector / p_chamber;
    const double dp_pump = out.required_tank_pressure - spec.pump_inlet_pressure;
    out.pump_power = (dp_pump > 0.0)
                         ? leg.mass_flow * dp_pump / (leg.density * spec.pump_efficiency)
                         : 0.0;
  };
  doLeg(spec.oxidizer, res.oxidizer);
  doLeg(spec.fuel, res.fuel);
  res.mixture_ratio = spec.oxidizer.mass_flow / spec.fuel.mass_flow;
  res.total_mass_flow = spec.oxidizer.mass_flow + spec.fuel.mass_flow;
  res.total_pump_power = res.oxidizer.pump_power + res.fuel.pump_power;
  res.stiffness_ok = (target_stiffness >= spec.minimum_stiffness);
  if (!res.stiffness_ok)
    res.notes = "requested injector stiffness is below the configured minimum";
  return res;
}

FeedSystemResult evaluateFeedSystem(const FeedSystemSpec& spec, double p_chamber) {
  if (!(p_chamber > 0.0)) throw ConfigError("feed system: chamber pressure must be positive");
  validate(spec.oxidizer, "oxidizer");
  validate(spec.fuel, "fuel");

  FeedSystemResult res;
  res.chamber_pressure = p_chamber;
  auto doLeg = [&](const FeedLeg& leg, FeedLegResult& out, const char* name) {
    if (!(leg.injector_area > 0.0))
      throw ConfigError(std::string("feed system: ") + name +
                        " injector area must be positive in capability mode");
    if (!(leg.tank_pressure > p_chamber)) {
      std::ostringstream os;
      os << "feed system: the " << name << " tank at " << leg.tank_pressure * 1e-6
         << " MPa cannot supply a chamber at " << p_chamber * 1e-6 << " MPa";
      throw InfeasibleError(os.str());
    }
    // Solve  p_tank - p_chamber = dp_lines(mdot) + dp_dynamic(mdot) + dp_inj(mdot)
    // by bisection; every term rises monotonically with mdot.
    const double budget = leg.tank_pressure - p_chamber;
    auto deficit = [&](double mdot) {
      FeedLegResult tmp;
      lineLosses(leg, mdot, tmp);
      const double v_inj = mdot / (leg.injector_cd * leg.density * leg.injector_area);
      const double dp_inj = 0.5 * leg.density * v_inj * v_inj;
      return tmp.dp_lines + tmp.dp_dynamic + dp_inj - budget;
    };
    double lo = 1e-9, hi = 1.0;
    for (int i = 0; i < 200 && deficit(hi) < 0.0; ++i) hi *= 2.0;
    if (deficit(hi) < 0.0)
      throw InfeasibleError(std::string("feed system: could not bracket the ") + name +
                            " flow; the line and injector losses are implausibly small");
    for (int i = 0; i < 300; ++i) {
      const double mid = 0.5 * (lo + hi);
      if (deficit(mid) < 0.0) lo = mid; else hi = mid;
      if (hi - lo < 1e-12 * hi) break;
    }
    out.mass_flow = 0.5 * (lo + hi);
    lineLosses(leg, out.mass_flow, out);
    out.injector_velocity = out.mass_flow / (leg.injector_cd * leg.density * leg.injector_area);
    out.dp_injector = 0.5 * leg.density * out.injector_velocity * out.injector_velocity;
    out.injector_area = leg.injector_area;
    out.required_tank_pressure = leg.tank_pressure;
    out.stiffness = out.dp_injector / p_chamber;
    const double dp_pump = leg.tank_pressure - spec.pump_inlet_pressure;
    out.pump_power = (dp_pump > 0.0)
                         ? out.mass_flow * dp_pump / (leg.density * spec.pump_efficiency)
                         : 0.0;
  };
  doLeg(spec.oxidizer, res.oxidizer, "oxidizer");
  doLeg(spec.fuel, res.fuel, "fuel");
  res.mixture_ratio = res.oxidizer.mass_flow / res.fuel.mass_flow;
  res.total_mass_flow = res.oxidizer.mass_flow + res.fuel.mass_flow;
  res.total_pump_power = res.oxidizer.pump_power + res.fuel.pump_power;
  res.stiffness_ok = res.oxidizer.stiffness >= spec.minimum_stiffness &&
                     res.fuel.stiffness >= spec.minimum_stiffness;
  if (!res.stiffness_ok) {
    std::ostringstream os;
    os << "injector stiffness (ox " << res.oxidizer.stiffness << ", fuel " << res.fuel.stiffness
       << ") is below the configured minimum of " << spec.minimum_stiffness;
    res.notes = os.str();
  }
  return res;
}

double pressurantMass(double tank_pressure, double propellant_volume, double pressurant_molar_mass,
                      double pressurant_temperature) {
  if (!(tank_pressure > 0.0 && propellant_volume > 0.0 && pressurant_molar_mass > 0.0 &&
        pressurant_temperature > 0.0))
    throw ConfigError("pressurant mass: every argument must be positive");
  return tank_pressure * propellant_volume * pressurant_molar_mass /
         (constants::R_universal * pressurant_temperature);
}

}  // namespace ignis
