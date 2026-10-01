// SPDX-License-Identifier: MIT
#include "ignis/cycle/Cycle.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <iomanip>
#include <memory>
#include <sstream>

#include "ignis/core/Constants.hpp"
#include "ignis/core/Exceptions.hpp"

namespace ignis {

std::string toString(CycleType type) {
  switch (type) {
    case CycleType::kGasGenerator: return "gas_generator";
    case CycleType::kStagedCombustion: return "staged_combustion";
    case CycleType::kExpander: return "expander";
  }
  return "unknown";
}

CycleType parseCycleType(const std::string& name) {
  if (name == "gas_generator") return CycleType::kGasGenerator;
  if (name == "staged_combustion") return CycleType::kStagedCombustion;
  if (name == "expander") return CycleType::kExpander;
  throw ConfigError("cycle: unknown type '" + name +
                    "' (expected gas_generator, staged_combustion or expander)");
}

namespace {

/// How one propellant is pumped: along its real-fluid table where the table
/// is the propellant's own, otherwise as an incompressible liquid at the
/// library density.  A surrogate table (RP-1's n-dodecane, 7 % light) would
/// put its own density error straight into the pump work, while a kerosene's
/// compressibility over a few tens of MPa is under one percent.
struct PumpFluid {
  std::unique_ptr<CoolantFluid> table;
  double density = 0.0;
  std::string model;
  double pMax() const { return table ? table->pMax() : 1.0e300; }
};

PumpFluid pumpFluid(const Propellant& prop, double p_in) {
  if (!prop.liquid)
    throw ConfigError("cycle: " + prop.name + " is not a liquid, and a pump-fed cycle pumps liquids");
  if (!(prop.density > 0.0)) throw ConfigError("cycle: " + prop.name + " has no density");
  PumpFluid f;
  f.density = prop.density;
  if (prop.coolant_table.empty()) {
    f.model = "incompressible at the library density (no property table)";
    return f;
  }
  auto table = std::make_unique<CoolantFluid>(CoolantFluid::load(prop.coolant_table));
  std::ostringstream os;
  os << std::fixed << std::setprecision(1);
  const double p_check = std::min(std::max(p_in, table->pMin()), table->pMax());
  double rho = 0.0;
  try {
    rho = table->at(prop.reference_temperature, p_check).rho;
  } catch (const RangeError&) {
    os << "incompressible at the library density: the '" << prop.coolant_table
       << "' table has no liquid at the propellant's storage state";
    f.model = os.str();
    return f;
  }
  const double deviation = rho / prop.density - 1.0;
  if (std::abs(deviation) > 0.02) {
    os << "incompressible at the library density: the '" << prop.coolant_table
       << "' table is a surrogate, " << 100.0 * deviation << " % off the propellant's density";
    f.model = os.str();
    return f;
  }
  if (p_in < table->pMin()) {
    os << "incompressible at the library density: the pump inlet is below the '"
       << prop.coolant_table << "' table's lowest pressure";
    f.model = os.str();
    return f;
  }
  f.table = std::move(table);
  f.model = "real fluid, integrated along the '" + prop.coolant_table + "' table";
  return f;
}

PumpResult pump(const PumpFluid& f, double mass_flow, double p_in, double t_in, double p_out,
                double efficiency) {
  return pumpLiquid(f.table.get(), f.density, mass_flow, p_in, t_in, p_out, efficiency);
}

/// Work per kg a pump puts into its liquid: the isentropic head over the
/// efficiency, which is also the liquid's actual enthalpy rise.
double workPerKg(const PumpResult& r) { return r.isentropic_head / r.efficiency; }

void checkFraction(double v, const char* what, bool allow_zero, bool allow_one) {
  const bool ok = (allow_zero ? v >= 0.0 : v > 0.0) && (allow_one ? v <= 1.0 : v < 1.0);
  if (!ok) {
    std::ostringstream os;
    os << "cycle: " << what << " must lie in " << (allow_zero ? "[0, " : "(0, ")
       << (allow_one ? "1]" : "1)") << ", got " << v;
    throw ConfigError(os.str());
  }
}

void validate(const CycleSpec& s, const CycleInputs& in) {
  checkFraction(s.pump_efficiency_oxidizer, "pump_efficiency_oxidizer", false, true);
  checkFraction(s.pump_efficiency_fuel, "pump_efficiency_fuel", false, true);
  checkFraction(s.boost_pump_efficiency, "boost_pump_efficiency", false, true);
  checkFraction(s.turbine_efficiency, "turbine_efficiency", false, true);
  checkFraction(s.mechanical_efficiency, "mechanical_efficiency", false, true);
  checkFraction(s.preburner_flow_fraction, "preburner_flow_fraction", false, true);
  checkFraction(s.turbine_bypass_fraction, "turbine_bypass_fraction", true, false);
  if (!(s.injector_stiffness >= 0.0)) throw ConfigError("cycle: injector_stiffness must be >= 0");
  if (!(s.line_loss_fraction >= 0.0)) throw ConfigError("cycle: line_loss_fraction must be >= 0");
  if (!(s.preburner_injector_stiffness >= 0.0))
    throw ConfigError("cycle: preburner_injector_stiffness must be >= 0");
  if (!(s.hot_gas_injector_stiffness >= 0.0))
    throw ConfigError("cycle: hot_gas_injector_stiffness must be >= 0");
  if (!(s.pump_inlet_pressure_oxidizer > 0.0) || !(s.pump_inlet_pressure_fuel > 0.0))
    throw ConfigError("cycle: pump inlet pressures must be positive");
  if (in.oxidizer == nullptr || in.fuel == nullptr)
    throw ConfigError("cycle: both propellants are required");
  if (!(in.chamber_pressure > 0.0)) throw ConfigError("cycle: chamber pressure must be positive");
  if (!(in.oxidizer_mass_flow > 0.0) || !(in.fuel_mass_flow > 0.0))
    throw ConfigError("cycle: propellant flows must be positive");
  if (!(in.oxidizer_temperature > 0.0) || !(in.fuel_temperature > 0.0))
    throw ConfigError("cycle: pump inlet temperatures must be positive");
  if (in.jacket && !(in.jacket_pressure_drop >= 0.0))
    throw ConfigError("cycle: the jacket pressure drop must be >= 0");
  const double p_main = in.chamber_pressure * (1.0 + s.injector_stiffness + s.line_loss_fraction);
  if (!(s.pump_inlet_pressure_oxidizer < p_main) || !(s.pump_inlet_pressure_fuel < p_main))
    throw ConfigError("cycle: a pump inlet pressure is above the injector supply pressure");
}

/// Smallest root of g on [x_lo, x_hi] (g(x_lo) < 0), found on a grid in
/// log(x) and closed by bisection; and, over the whole range, the largest
/// `ratio`, refined by golden section.
struct RootSearch {
  bool found = false;
  double x = 0.0;
  double best_ratio = 0.0;
  double best_x = 0.0;
};

RootSearch smallestRoot(const std::function<double(double)>& g,
                        const std::function<double(double)>& ratio, double x_lo, double x_hi) {
  constexpr int kGrid = 48;
  RootSearch out;
  const double la = std::log(x_lo), lb = std::log(x_hi);
  const double step = (lb - la) / kGrid;
  double prev_l = la;
  int root_cell = -1;
  for (int k = 1; k <= kGrid; ++k) {
    const double l = la + step * k;
    const double x = std::exp(l);
    const double r = ratio(x);
    if (r > out.best_ratio) { out.best_ratio = r; out.best_x = x; }
    if (root_cell < 0 && g(x) >= 0.0) root_cell = k;
    if (root_cell < 0) prev_l = l;
  }
  if (root_cell > 0) {
    double a = prev_l, b = la + step * root_cell;
    for (int it = 0; it < 200 && b - a > 1.0e-12; ++it) {
      const double m = 0.5 * (a + b);
      if (g(std::exp(m)) >= 0.0) b = m; else a = m;
    }
    out.found = true;
    out.x = std::exp(b);
  }
  // The best power ratio, refined within a grid cell either side of the
  // best sample.
  double a = std::max(la, std::log(out.best_x) - step);
  double b = std::min(lb, std::log(out.best_x) + step);
  const double phi = 0.5 * (std::sqrt(5.0) - 1.0);
  for (int it = 0; it < 40; ++it) {
    const double c = b - phi * (b - a), d = a + phi * (b - a);
    if (ratio(std::exp(c)) > ratio(std::exp(d))) b = d; else a = c;
  }
  const double x = std::exp(0.5 * (a + b));
  const double r = ratio(x);
  if (r > out.best_ratio) { out.best_ratio = r; out.best_x = x; }
  return out;
}

/// Ideal frozen expansion of the generator exhaust from the turbine outlet
/// to `p_target`, stopped where its isentrope reaches the species data's
/// lower temperature limit.
void exhaust(const GasMixture& mix, const SpeciesDatabase& db, const CycleInputs& in,
             CycleResult& r) {
  const auto& n = r.combustor.n;
  const double p_t = r.turbine.p_out, t_t = r.turbine.t_out;
  double p_x = in.exit_pressure;
  r.exhaust_exit_pressure = p_x;
  if (!(p_x > 0.0) || p_t <= p_x) {
    r.warnings.push_back(
        "the gas-generator turbine exhausts at or below the main nozzle's exit pressure, so no "
        "exhaust thrust is counted");
    return;
  }
  const double s_t = mix.entropy(n, t_t, p_t);
  const double r_gas = constants::R_universal * n.sum();
  const double t_floor = db.tMinCommon() + 1.0;
  const double s_floor = mix.entropy(n, t_floor, p_x);
  if (s_floor > s_t) {
    // On an ideal-gas mixture s depends on pressure only through -R ln p.
    p_x *= std::exp((s_floor - s_t) / r_gas);
    std::ostringstream os;
    os << "the generator exhaust's frozen isentrope reaches the species data's " << t_floor
       << " K floor at " << std::setprecision(3) << p_x * 1e-3
       << " kPa; its expansion is stopped there (water would condense first in reality)";
    r.warnings.push_back(os.str());
    r.exhaust_exit_pressure = p_x;
  }
  const double t_x = mix.temperatureFromEntropy(n, s_t, p_x, std::max(t_floor, 0.5 * t_t));
  const double dh = mix.enthalpy(n, t_t) - mix.enthalpy(n, t_x);
  r.exhaust_velocity = std::sqrt(std::max(2.0 * dh, 0.0));
  if (!(r.exhaust_velocity > 0.0)) return;
  const double rho = p_x / (r_gas * t_x);
  const double area = r.gg_mass_flow / (rho * r.exhaust_velocity);
  r.exhaust_thrust_vacuum = r.gg_mass_flow * r.exhaust_velocity + p_x * area;
  r.exhaust_thrust = r.exhaust_thrust_vacuum - in.ambient_pressure * area;
}

void delivered(const CycleInputs& in, CycleResult& r) {
  const double m_main = in.oxidizer_mass_flow + in.fuel_mass_flow;
  r.total_mass_flow = m_main + r.gg_mass_flow;
  r.overall_mixture_ratio = (in.oxidizer_mass_flow + r.gg_oxidizer_flow) /
                            (in.fuel_mass_flow + r.gg_fuel_flow);
  r.thrust_vacuum_delivered = in.thrust_vacuum + r.exhaust_thrust_vacuum;
  r.thrust_delivered = in.thrust + r.exhaust_thrust;
  const double w = r.total_mass_flow * constants::g0;
  r.isp_vacuum_delivered = r.thrust_vacuum_delivered / w;
  r.isp_delivered = r.thrust_delivered / w;
  r.isp_vacuum_dumped = in.thrust_vacuum / w;
}

void finishBalance(const CycleSpec& spec, CycleResult& r) {
  r.pump_power = r.oxidizer_pump.power + r.fuel_pump.power +
                 (r.has_boost_pump ? r.boost_pump.power : 0.0);
  r.turbine_power = r.turbine.power;
  r.power_balance_residual =
      (spec.mechanical_efficiency * r.turbine_power - r.pump_power) / r.pump_power;
  r.turbine_pressure_ratio = r.turbine.pressure_ratio;
  r.turbine_inlet_temperature = r.turbine.t_in;
  r.oxidizer_discharge_pressure = r.oxidizer_pump.p_out;
  r.fuel_discharge_pressure = r.fuel_pump.p_out;
}

// --- gas generator -------------------------------------------------------

void gasGenerator(const CycleSpec& spec, const CycleInputs& in, const EquilibriumSolver& solver,
                  const PumpFluid& oxf, const PumpFluid& fuf, CycleResult& r) {
  const double pc = in.chamber_pressure;
  if (!(spec.turbine_pressure_ratio > 1.0))
    throw ConfigError("cycle: the gas-generator turbine pressure ratio must exceed 1");
  const double p_gg = spec.gas_generator_pressure > 0.0 ? spec.gas_generator_pressure : pc;
  if (p_gg > pc * (1.0 + 1e-12))
    throw ConfigError(
        "cycle: the gas generator is fed from the main pump discharges, which deliver the main "
        "chamber's pressure; it cannot run above it");
  const double p_main = pc * (1.0 + spec.injector_stiffness + spec.line_loss_fraction);
  const double p_d_ox = p_main;
  const double p_need_f = p_main + (in.jacket ? in.jacket_pressure_drop : 0.0);
  const double p_d_f = std::max(p_need_f, in.fuel_discharge_floor);
  r.fuel_throttle_loss = p_d_f - p_need_f;

  // Work per kg at these discharge pressures (flow-independent).
  const auto ox1 = pump(oxf, 1.0, spec.pump_inlet_pressure_oxidizer, in.oxidizer_temperature,
                        p_d_ox, spec.pump_efficiency_oxidizer);
  const auto fu1 = pump(fuf, 1.0, spec.pump_inlet_pressure_fuel, in.fuel_temperature, p_d_f,
                        spec.pump_efficiency_fuel);
  const double w_ox = workPerKg(ox1), w_f = workPerKg(fu1);

  // The generator's propellants leave the pumps before the jacket.
  r.combustor = combustorAtTemperature(solver, *in.oxidizer, *in.fuel, in.oxidizer_temperature,
                                       in.fuel_temperature, p_gg, spec.turbine_inlet_temperature,
                                       spec.fuel_rich, w_ox, w_f);
  r.has_combustor = true;
  r.combustor_pressure = p_gg;
  const double mr = r.combustor.mixture_ratio;
  const double x_ox = mr / (1.0 + mr), x_f = 1.0 / (1.0 + mr);
  const auto& mix = solver.mixture();
  const auto t1 = gasTurbine(mix, r.combustor.n, 1.0, r.combustor.temperature, p_gg,
                             spec.turbine_pressure_ratio, spec.turbine_efficiency);
  const double dh_t = t1.power;  // J/kg

  const double need = in.oxidizer_mass_flow * w_ox + in.fuel_mass_flow * w_f;
  const double net = spec.mechanical_efficiency * dh_t - x_ox * w_ox - x_f * w_f;
  if (!(net > 0.0)) {
    std::ostringstream os;
    os << "each kilogram through the gas generator's turbine yields "
       << spec.mechanical_efficiency * dh_t * 1e-3 << " kJ at the shaft but costs "
       << (x_ox * w_ox + x_f * w_f) * 1e-3 << " kJ to pump: no generator flow closes the cycle";
    r.feasible = false;
    r.infeasibility = os.str();
    r.oxidizer_pump = pump(oxf, in.oxidizer_mass_flow, spec.pump_inlet_pressure_oxidizer,
                           in.oxidizer_temperature, p_d_ox, spec.pump_efficiency_oxidizer);
    r.fuel_pump = pump(fuf, in.fuel_mass_flow, spec.pump_inlet_pressure_fuel, in.fuel_temperature,
                       p_d_f, spec.pump_efficiency_fuel);
    r.turbine = t1;
    finishBalance(spec, r);
    delivered(in, r);
    return;
  }
  r.gg_mass_flow = need / net;
  r.gg_oxidizer_flow = x_ox * r.gg_mass_flow;
  r.gg_fuel_flow = x_f * r.gg_mass_flow;
  r.oxidizer_pump = pump(oxf, in.oxidizer_mass_flow + r.gg_oxidizer_flow,
                         spec.pump_inlet_pressure_oxidizer, in.oxidizer_temperature, p_d_ox,
                         spec.pump_efficiency_oxidizer);
  r.fuel_pump = pump(fuf, in.fuel_mass_flow + r.gg_fuel_flow, spec.pump_inlet_pressure_fuel,
                     in.fuel_temperature, p_d_f, spec.pump_efficiency_fuel);
  r.turbine = gasTurbine(mix, r.combustor.n, r.gg_mass_flow, r.combustor.temperature, p_gg,
                         spec.turbine_pressure_ratio, spec.turbine_efficiency);
  finishBalance(spec, r);
  exhaust(mix, solver.database(), in, r);
  delivered(in, r);
  r.gg_flow_fraction = r.gg_mass_flow / r.total_mass_flow;
}

// --- staged combustion ---------------------------------------------------

void stagedCombustion(const CycleSpec& spec, const CycleInputs& in,
                      const EquilibriumSolver& solver, const PumpFluid& oxf, const PumpFluid& fuf,
                      CycleResult& r) {
  const double pc = in.chamber_pressure;
  const bool fuel_rich = spec.fuel_rich;
  const double p_main = pc * (1.0 + spec.injector_stiffness + spec.line_loss_fraction);
  const double dp_jacket = in.jacket ? in.jacket_pressure_drop : 0.0;
  const double p_t_out = pc * (1.0 + spec.hot_gas_injector_stiffness);
  const double m_major = fuel_rich ? in.fuel_mass_flow : in.oxidizer_mass_flow;
  const double m_major_pb = spec.preburner_flow_fraction * m_major;
  // Heat the jacket put into the fuel, per kg of all of it: on a fuel-rich
  // preburner the jacket comes first and its flow mixes with any bypass.
  const double q_fuel = (in.jacket && fuel_rich) ? in.jacket_heat / in.fuel_mass_flow : 0.0;
  const PumpFluid& minor_fluid = fuel_rich ? oxf : fuf;
  const double eta_minor = spec.boost_pump_efficiency;

  // The preburner supply pressure for a turbine pressure ratio, and the
  // largest ratio the property tables can follow.
  auto p_feed = [&](double pr) {
    return pr * p_t_out * (1.0 + spec.preburner_injector_stiffness) +
           spec.line_loss_fraction * pc;
  };
  // Which pump hits its property table's ceiling first as the preburner
  // pressure rises: the major propellant's main pump (through the jacket, on
  // a fuel-rich cycle) or the minor propellant's boost stage.
  const PumpFluid& major_fluid = fuel_rich ? fuf : oxf;
  const double major_limit = major_fluid.pMax() - (fuel_rich ? dp_jacket : 0.0);
  const double minor_limit = minor_fluid.pMax();
  const bool major_binds = major_limit <= minor_limit;
  const double feed_limit = 0.999 * std::min(major_limit, minor_limit);
  const double pr_table = (feed_limit - spec.line_loss_fraction * pc) /
                          (p_t_out * (1.0 + spec.preburner_injector_stiffness));
  const double pr_max = std::min(10.0, pr_table);
  if (!(pr_max > 1.0 + 1e-6)) {
    std::ostringstream os;
    os << "the preburner would need a supply pressure above the property tables' "
       << feed_limit * 1e-6 << " MPa at any turbine pressure ratio";
    r.feasible = false;
    r.infeasibility = os.str();
    return;
  }

  struct State {
    PumpResult ox, fu, boost;
    bool boosted = false;
    double minor_flow = 0.0;
    TurbineResult turbine;
  };
  // All pumps at one turbine pressure ratio, for a given preburner gas.
  // An oxidiser-rich cycle fixes the fuel pump's discharge by the injector
  // and jacket alone, so a configured jacket pressure above that is throttled.
  const double p_need_f_oxrich = p_main + dp_jacket;
  const double p_d_f_oxrich = std::max(p_need_f_oxrich, in.fuel_discharge_floor);
  if (!fuel_rich) r.fuel_throttle_loss = p_d_f_oxrich - p_need_f_oxrich;
  auto evaluate = [&](double pr, const CombustorResult& c) {
    State s;
    const double feed = p_feed(pr);
    const double p_d_f = fuel_rich ? feed + dp_jacket : p_d_f_oxrich;
    const double p_d_ox = fuel_rich ? p_main : feed;
    s.ox = pump(oxf, in.oxidizer_mass_flow, spec.pump_inlet_pressure_oxidizer,
                in.oxidizer_temperature, p_d_ox, spec.pump_efficiency_oxidizer);
    s.fu = pump(fuf, in.fuel_mass_flow, spec.pump_inlet_pressure_fuel, in.fuel_temperature,
                p_d_f, spec.pump_efficiency_fuel);
    s.minor_flow = fuel_rich ? c.mixture_ratio * m_major_pb : m_major_pb / c.mixture_ratio;
    const PumpResult& from = fuel_rich ? s.ox : s.fu;
    if (feed > from.p_out) {
      s.boost = pump(minor_fluid, s.minor_flow, from.p_out, from.t_out, feed, eta_minor);
      s.boosted = true;
    }
    s.turbine = gasTurbine(solver.mixture(), c.n, m_major_pb + s.minor_flow, c.temperature,
                           pr * p_t_out, pr, spec.turbine_efficiency);
    return s;
  };
  auto pumpPower = [](const State& s) {
    return s.ox.power + s.fu.power + (s.boosted ? s.boost.power : 0.0);
  };
  // The preburner gas for a turbine pressure ratio: its propellants arrive
  // with the pump work (and on a fuel-rich preburner the jacket heat) in them.
  auto preburner = [&](double pr, const State* s) {
    double dh_ox = 0.0, dh_f = 0.0;
    if (s != nullptr) {
      dh_ox = workPerKg(s->ox) + (fuel_rich && s->boosted ? workPerKg(s->boost) : 0.0);
      dh_f = workPerKg(s->fu) + q_fuel + (!fuel_rich && s->boosted ? workPerKg(s->boost) : 0.0);
    }
    return combustorAtTemperature(solver, *in.oxidizer, *in.fuel, in.oxidizer_temperature,
                                  in.fuel_temperature, pr * p_t_out,
                                  spec.turbine_inlet_temperature, fuel_rich, dh_ox, dh_f);
  };

  // Alternate: the preburner gas at the current pressure ratio, then the
  // pressure ratio that balances with that gas, until the two agree.
  double pr = std::min(1.5, 0.5 * (1.0 + pr_max));
  CombustorResult gas = preburner(pr, nullptr);
  {
    const State s0 = evaluate(pr, gas);
    gas = preburner(pr, &s0);
  }
  RootSearch root;
  for (int pass = 0; pass < 8; ++pass) {
    auto g = [&](double x) {
      const State s = evaluate(x, gas);
      return spec.mechanical_efficiency * s.turbine.power - pumpPower(s);
    };
    auto ratio = [&](double x) {
      const State s = evaluate(x, gas);
      return spec.mechanical_efficiency * s.turbine.power / pumpPower(s);
    };
    root = smallestRoot(g, ratio, 1.0 + 1e-9, pr_max);
    const double pr_new = root.found ? root.x : root.best_x;
    const State s = evaluate(pr_new, gas);
    const CombustorResult next = preburner(pr_new, &s);
    const bool settled = std::abs(next.mixture_ratio - gas.mixture_ratio) <=
                         1e-10 * gas.mixture_ratio;
    gas = next;
    pr = pr_new;
    if (settled) break;
  }
  r.max_power_ratio = root.best_ratio;
  if (!root.found) {
    std::ostringstream os;
    os << std::setprecision(4) << "no turbine pressure ratio up to " << pr_max;
    if (pr_table <= 10.0)
      os << " (where the " << (fuel_rich == major_binds ? "fuel" : "oxidizer")
         << (major_binds ? " pump" : " boost stage") << " reaches its property table's "
         << (major_binds ? major_fluid : minor_fluid).pMax() * 1e-6 << " MPa ceiling)";
    os << " balances: at best (PR " << root.best_x << ") the turbine makes "
       << 100.0 * root.best_ratio << " % of the pump power";
    r.feasible = false;
    r.infeasibility = os.str();
  }
  const State s = evaluate(pr, gas);
  r.combustor = gas;
  r.has_combustor = true;
  r.combustor_pressure = pr * p_t_out;
  r.oxidizer_pump = s.ox;
  r.fuel_pump = s.fu;
  r.has_boost_pump = s.boosted;
  if (s.boosted) {
    r.boost_pump = s.boost;
    r.boost_propellant = fuel_rich ? "oxidizer" : "fuel";
  }
  r.turbine = s.turbine;
  r.preburner_major_flow = m_major_pb;
  r.preburner_minor_flow = s.minor_flow;
  finishBalance(spec, r);
  delivered(in, r);
}

// --- expander ------------------------------------------------------------

void expander(const CycleSpec& spec, const CycleInputs& in, const PumpFluid& oxf,
              const PumpFluid& fuf, CycleResult& r) {
  if (!in.jacket || in.jacket_fluid == nullptr)
    throw ConfigError("cycle: an expander is driven by the fuel the cooling jacket heats, so it "
                      "needs a fuel-cooled jacket ('cooling.enabled' with the fuel as coolant)");
  if (!(in.jacket_mass_flow > 0.0) || !(in.jacket_outlet_temperature > 0.0))
    throw ConfigError("cycle: the expander needs the jacket's flow and outlet temperature");
  const double pc = in.chamber_pressure;
  const double p_main = pc * (1.0 + spec.injector_stiffness + spec.line_loss_fraction);
  const double dp_jacket = in.jacket_pressure_drop;
  const double p_t_out = p_main;  // turbine -> valves and lines -> main injector
  const double m_t = (1.0 - spec.turbine_bypass_fraction) * in.jacket_mass_flow;
  const double t_in = in.jacket_outlet_temperature;
  const auto ox = pump(oxf, in.oxidizer_mass_flow, spec.pump_inlet_pressure_oxidizer,
                       in.oxidizer_temperature, p_main, spec.pump_efficiency_oxidizer);

  const double x_max = 0.999 * std::min(fuf.pMax(), in.jacket_fluid->pMax());
  const double pr_table = (x_max - dp_jacket) / p_t_out;
  const double pr_max = std::min(20.0, pr_table);
  if (!(pr_max > 1.0 + 1e-6)) {
    r.feasible = false;
    r.infeasibility = "the jacket's pressure drop alone takes the fuel pump past the property "
                      "table's highest pressure";
    r.oxidizer_pump = ox;
    return;
  }
  struct State {
    PumpResult fu;
    TurbineResult turbine;
  };
  auto evaluate = [&](double pr) {
    State s;
    const double p_t_in = pr * p_t_out;
    s.fu = pump(fuf, in.fuel_mass_flow, spec.pump_inlet_pressure_fuel, in.fuel_temperature,
                p_t_in + dp_jacket, spec.pump_efficiency_fuel);
    s.turbine = fluidTurbine(*in.jacket_fluid, m_t, t_in, p_t_in, pr, spec.turbine_efficiency);
    return s;
  };
  auto g = [&](double pr) {
    const State s = evaluate(pr);
    return spec.mechanical_efficiency * s.turbine.power - (ox.power + s.fu.power);
  };
  auto ratio = [&](double pr) {
    const State s = evaluate(pr);
    return spec.mechanical_efficiency * s.turbine.power / (ox.power + s.fu.power);
  };
  const RootSearch root = smallestRoot(g, ratio, 1.0 + 1e-9, pr_max);
  r.max_power_ratio = root.best_ratio;
  const double pr = root.found ? root.x : root.best_x;
  if (!root.found) {
    std::ostringstream os;
    os << std::setprecision(4) << "the jacket does not heat the fuel enough: at best (turbine PR "
       << root.best_x << ") the turbine makes " << 100.0 * root.best_ratio
       << " % of the pump power";
    if (pr_table <= 20.0)
      os << " (pressure ratios are searched up to " << pr_max
         << ", where the fuel pump reaches its property table's " << x_max * 1e-6
         << " MPa ceiling)";
    r.feasible = false;
    r.infeasibility = os.str();
  }
  const State s = evaluate(pr);
  r.oxidizer_pump = ox;
  r.fuel_pump = s.fu;
  r.turbine = s.turbine;
  finishBalance(spec, r);
  delivered(in, r);
}

}  // namespace

double pumpDischargeTemperature(const Propellant& prop, double t_in, double p_in, double p_out,
                                double efficiency) {
  const PumpFluid f = pumpFluid(prop, p_in);
  return pump(f, 1.0, p_in, t_in, p_out, efficiency).t_out;
}

CycleResult solveCycle(const CycleSpec& spec, const CycleInputs& in,
                       const EquilibriumSolver& solver) {
  validate(spec, in);
  CycleResult r;
  r.type = spec.type;
  r.chamber_pressure = in.chamber_pressure;
  const PumpFluid oxf = pumpFluid(*in.oxidizer, spec.pump_inlet_pressure_oxidizer);
  const PumpFluid fuf = pumpFluid(*in.fuel, spec.pump_inlet_pressure_fuel);
  r.oxidizer_pump_model = oxf.model;
  r.fuel_pump_model = fuf.model;
  switch (spec.type) {
    case CycleType::kGasGenerator: gasGenerator(spec, in, solver, oxf, fuf, r); break;
    case CycleType::kStagedCombustion: stagedCombustion(spec, in, solver, oxf, fuf, r); break;
    case CycleType::kExpander: expander(spec, in, oxf, fuf, r); break;
  }
  if (r.has_combustor) {
    const bool hydrocarbon = in.fuel->composition.count("C") > 0;
    if (hydrocarbon && r.combustor.mixture_ratio <
                           PropellantMixture(*in.oxidizer, *in.fuel, 1.0, in.oxidizer_temperature,
                                             in.fuel_temperature)
                               .stoichiometricMixtureRatio())
      r.warnings.push_back(
          "the fuel-rich hydrocarbon generator gas is taken as chemical equilibrium; a real one "
          "makes soot and cracked fuel instead, so its composition (not its temperature) is "
          "approximate");
  }
  return r;
}

std::string CycleResult::summary() const {
  std::ostringstream os;
  os << std::fixed;
  os << "engine cycle (" << toString(type) << ")";
  if (!feasible) os << "  -- DOES NOT CLOSE";
  os << "\n";
  if (!feasible) os << "  " << infeasibility << "\n";
  auto pumpLine = [&](const char* name, const PumpResult& p) {
    os << "  " << std::left << std::setw(20) << name << std::right << std::setprecision(3)
       << p.mass_flow << " kg/s, " << std::setprecision(3) << p.p_in * 1e-6 << " -> "
       << p.p_out * 1e-6 << " MPa, " << std::setprecision(1) << p.power * 1e-3 << " kW (eta "
       << std::setprecision(3) << p.efficiency << ")\n";
  };
  pumpLine("oxidizer pump", oxidizer_pump);
  os << "                      " << oxidizer_pump_model << "\n";
  pumpLine("fuel pump", fuel_pump);
  os << "                      " << fuel_pump_model << "\n";
  if (has_boost_pump) pumpLine((boost_propellant + " boost").c_str(), boost_pump);
  if (has_combustor) {
    os << "  " << (type == CycleType::kGasGenerator ? "gas generator       " : "preburner           ")
       << "O/F " << std::setprecision(4) << combustor.mixture_ratio << " at "
       << std::setprecision(1) << combustor.temperature << " K, " << std::setprecision(3)
       << combustor_pressure * 1e-6 << " MPa\n";
  }
  os << "  turbine             " << std::setprecision(3) << turbine.mass_flow << " kg/s, "
     << std::setprecision(1) << turbine.t_in << " K in, PR " << std::setprecision(4)
     << turbine.pressure_ratio << ", " << std::setprecision(1) << turbine_power * 1e-3
     << " kW (eta " << std::setprecision(3) << turbine.efficiency << ")\n";
  os << "  power balance       pumps " << std::setprecision(1) << pump_power * 1e-3
     << " kW; residual " << std::scientific << std::setprecision(2) << power_balance_residual
     << std::fixed << "\n";
  if (max_power_ratio > 0.0)
    os << "  best power ratio    " << std::setprecision(3) << max_power_ratio
       << " (turbine at its best pressure ratio / pumps)\n";
  if (fuel_throttle_loss > 0.0)
    os << "  fuel throttled      " << std::setprecision(3) << fuel_throttle_loss * 1e-6
       << " MPa: the pump delivers the jacket's configured inlet pressure, more than the "
          "injector needs\n";
  if (type == CycleType::kGasGenerator && feasible) {
    os << "  generator flow      " << std::setprecision(4) << gg_mass_flow << " kg/s ("
       << std::setprecision(2) << 100.0 * gg_flow_fraction << " % of the engine's)\n"
       << "  exhaust             " << std::setprecision(1) << exhaust_velocity
       << " m/s ideal, to " << std::setprecision(2) << exhaust_exit_pressure * 1e-3 << " kPa\n"
       << "  delivered Isp (vac) " << std::setprecision(2) << isp_vacuum_delivered
       << " s with the exhaust's ideal thrust, " << isp_vacuum_dumped << " s dumped\n"
       << "  overall O/F         " << std::setprecision(4) << overall_mixture_ratio << "\n";
  }
  for (const auto& w : warnings) os << "  WARNING: " << w << "\n";
  return os.str();
}

}  // namespace ignis
