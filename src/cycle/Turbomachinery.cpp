// SPDX-License-Identifier: MIT
#include "ignis/cycle/Turbomachinery.hpp"

#include <algorithm>
#include <cmath>
#include <sstream>

#include "ignis/core/Exceptions.hpp"

namespace ignis {

namespace {

/// March an isentrope along a property table from (t0, p0) to p1:
/// dh = dp / rho(T, p), midpoint in each step, T recovered from (h, p).
/// Returns the isentropic enthalpy change h(p1) - h(p0).
double tableIsentrope(const CoolantFluid& fluid, double t0, double p0, double p1, double* t_end) {
  constexpr int kSteps = 400;
  double h = fluid.at(t0, p0).h;
  double T = t0;
  const double h0 = h;
  for (int k = 0; k < kSteps; ++k) {
    const double pa = p0 + (p1 - p0) * k / kSteps;
    const double pb = p0 + (p1 - p0) * (k + 1) / kSteps;
    const double pm = 0.5 * (pa + pb);
    const double rho_a = fluid.at(T, pa).rho;
    const double t_mid = fluid.temperatureFromEnthalpy(h + 0.5 * (pm - pa) / rho_a, pm, T);
    const double rho_m = fluid.at(t_mid, pm).rho;
    h += (pb - pa) / rho_m;
    T = fluid.temperatureFromEnthalpy(h, pb, t_mid);
  }
  if (t_end != nullptr) *t_end = T;
  return h - h0;
}

}  // namespace

PumpResult pumpLiquid(const CoolantFluid* table, double density, double mass_flow, double p_in,
                      double t_in, double p_out, double efficiency) {
  if (!(mass_flow > 0.0)) throw ConfigError("pump: mass flow must be positive");
  if (!(p_out > p_in)) throw ConfigError("pump: discharge pressure must exceed inlet pressure");
  if (!(efficiency > 0.0 && efficiency <= 1.0))
    throw ConfigError("pump: efficiency must lie in (0, 1]");
  PumpResult r;
  r.mass_flow = mass_flow;
  r.p_in = p_in;
  r.p_out = p_out;
  r.t_in = t_in;
  r.efficiency = efficiency;
  if (table != nullptr) {
    r.real_fluid = true;
    r.isentropic_head = tableIsentrope(*table, t_in, p_in, p_out, nullptr);
    const double h_out = table->at(t_in, p_in).h + r.isentropic_head / efficiency;
    r.t_out = table->temperatureFromEnthalpy(h_out, p_out, t_in);
  } else {
    if (!(density > 0.0)) throw ConfigError("pump: density must be positive");
    r.isentropic_head = (p_out - p_in) / density;
    r.t_out = t_in;   // no caloric model without a table
  }
  r.power = mass_flow * r.isentropic_head / efficiency;
  return r;
}

TurbineResult gasTurbine(const GasMixture& mix, const Eigen::VectorXd& n, double mass_flow,
                         double t_in, double p_in, double pressure_ratio, double efficiency) {
  if (!(mass_flow > 0.0)) throw ConfigError("turbine: mass flow must be positive");
  if (!(pressure_ratio > 1.0)) throw ConfigError("turbine: pressure ratio must exceed 1");
  if (!(efficiency > 0.0 && efficiency <= 1.0))
    throw ConfigError("turbine: efficiency must lie in (0, 1]");
  TurbineResult r;
  r.mass_flow = mass_flow;
  r.p_in = p_in;
  r.p_out = p_in / pressure_ratio;
  r.t_in = t_in;
  r.pressure_ratio = pressure_ratio;
  r.efficiency = efficiency;
  const double s_in = mix.entropy(n, t_in, p_in);
  const double h_in = mix.enthalpy(n, t_in);
  const double t_s = mix.temperatureFromEntropy(n, s_in, r.p_out, 0.8 * t_in);
  r.isentropic_drop = h_in - mix.enthalpy(n, t_s);
  r.power = mass_flow * efficiency * r.isentropic_drop;
  r.t_out = mix.temperatureFromEnthalpy(n, h_in - efficiency * r.isentropic_drop, t_s);
  return r;
}

TurbineResult fluidTurbine(const CoolantFluid& fluid, double mass_flow, double t_in, double p_in,
                           double pressure_ratio, double efficiency) {
  if (!(mass_flow > 0.0)) throw ConfigError("turbine: mass flow must be positive");
  if (!(pressure_ratio > 1.0)) throw ConfigError("turbine: pressure ratio must exceed 1");
  if (!(efficiency > 0.0 && efficiency <= 1.0))
    throw ConfigError("turbine: efficiency must lie in (0, 1]");
  TurbineResult r;
  r.mass_flow = mass_flow;
  r.p_in = p_in;
  r.p_out = p_in / pressure_ratio;
  r.t_in = t_in;
  r.pressure_ratio = pressure_ratio;
  r.efficiency = efficiency;
  r.isentropic_drop = -tableIsentrope(fluid, t_in, p_in, r.p_out, nullptr);
  r.power = mass_flow * efficiency * r.isentropic_drop;
  const double h_out = fluid.at(t_in, p_in).h - efficiency * r.isentropic_drop;
  r.t_out = fluid.temperatureFromEnthalpy(h_out, r.p_out, t_in);
  return r;
}

CombustorResult combustorAtTemperature(const EquilibriumSolver& solver, const Propellant& oxidizer,
                                       const Propellant& fuel, double t_oxidizer, double t_fuel,
                                       double pressure, double target, bool fuel_rich,
                                       double oxidizer_enthalpy_gain, double fuel_enthalpy_gain) {
  if (!(target > 300.0)) throw ConfigError("combustor: turbine inlet temperature must exceed 300 K");
  if (!(pressure > 0.0)) throw ConfigError("combustor: pressure must be positive");
  const auto& db = solver.database();
  // Reactant enthalpy per kg of mixture: each propellant at its storage
  // state, plus whatever the pumps and the jacket have added on the way.
  auto reactants = [&](const PropellantMixture& mix) {
    return mix.enthalpy(db) + mix.oxidizerMassFraction() * oxidizer_enthalpy_gain +
           mix.fuelMassFraction() * fuel_enthalpy_gain;
  };
  // The mixture ratio whose adiabatic flame temperature is the target is the
  // one whose reactant enthalpy equals that of its own equilibrium products
  // AT the target: h_reactants(MR) = h_eq(MR; T_target, p).  Solving that
  // balance with fixed-temperature equilibria, instead of comparing adiabatic
  // flame temperatures with the target, never asks for an adiabatic state
  // colder than the species data reach -- a cryogenic hydrogen-rich mixture
  // far from stoichiometric burns below their 200 K floor.
  auto balance = [&](double mr) {
    const PropellantMixture mix(oxidizer, fuel, mr, t_oxidizer, t_fuel);
    const auto b = mix.elementMoles(db);
    return reactants(mix) - solver.tp(b, target, pressure).state.h;
  };
  // Bracket on the requested side of stoichiometric.  Towards stoichiometric
  // the adiabatic temperature rises above any turbine's limit (the balance is
  // positive); far from it the products at the target carry more enthalpy
  // than the cold reactants can supply (negative).
  const double mr_st = PropellantMixture(oxidizer, fuel, 1.0, t_oxidizer, t_fuel)
                           .stoichiometricMixtureRatio();
  double lo, hi;
  if (fuel_rich) {
    lo = 1.0e-3 * mr_st;
    hi = 0.9 * mr_st;
  } else {
    lo = 1.1 * mr_st;
    hi = 1.0e3 * mr_st;
  }
  double f_lo = balance(lo), f_hi = balance(hi);
  if (f_lo * f_hi > 0.0) {
    std::ostringstream os;
    os << "combustor: no " << (fuel_rich ? "fuel" : "oxidiser") << "-rich mixture ratio of "
       << oxidizer.name << "/" << fuel.name << " burns at " << target << " K at "
       << pressure * 1e-6 << " MPa";
    throw InfeasibleError(os.str());
  }
  // Bisection in log(MR): the temperature is steep near stoichiometric.
  double a = std::log(lo), b = std::log(hi);
  for (int it = 0; it < 200 && b - a > 1.0e-13; ++it) {
    const double m = 0.5 * (a + b);
    const double fm = balance(std::exp(m));
    if ((fm > 0.0) == (f_hi > 0.0)) { b = m; f_hi = fm; } else { a = m; f_lo = fm; }
  }
  CombustorResult r;
  r.mixture_ratio = std::exp(0.5 * (a + b));
  const PropellantMixture mix(oxidizer, fuel, r.mixture_ratio, t_oxidizer, t_fuel);
  const auto st = solver.tp(mix.elementMoles(db), target, pressure).state;
  r.temperature = st.T;
  r.pressure = pressure;
  r.n = st.n;
  r.molar_mass = st.M;
  r.cp_frozen = st.cp_frozen;
  r.enthalpy_residual = reactants(mix) - st.h;
  return r;
}

}  // namespace ignis
