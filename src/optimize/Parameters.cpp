// SPDX-License-Identifier: MIT
#include "ignis/optimize/Parameters.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <sstream>

#include "ignis/nozzle/Atmosphere.hpp"

namespace ignis {
namespace {

struct ParamEntry {
  const char* units;
  void (*set)(EngineConfig&, double);
  double (*get)(const EngineConfig&);
};

const std::map<std::string, ParamEntry>& paramTable() {
  static const std::map<std::string, ParamEntry> table = {
      {"chamber.pressure",
       {"Pa", [](EngineConfig& c, double v) { c.chamber_pressure = v; },
        [](const EngineConfig& c) { return c.chamber_pressure; }}},
      {"chamber.eta_c_star",
       {"-", [](EngineConfig& c, double v) { c.eta_c_star = v; c.performance.eta_c_star = v; },
        [](const EngineConfig& c) { return c.eta_c_star; }}},
      {"propellants.mixture_ratio",
       {"-", [](EngineConfig& c, double v) { c.mixture_ratio = v; },
        [](const EngineConfig& c) { return c.mixture_ratio; }}},
      {"propellants.oxidizer_temperature",
       {"K", [](EngineConfig& c, double v) { c.oxidizer_temperature = v; },
        [](const EngineConfig& c) { return c.oxidizer_temperature; }}},
      {"propellants.fuel_temperature",
       {"K", [](EngineConfig& c, double v) { c.fuel_temperature = v; },
        [](const EngineConfig& c) { return c.fuel_temperature; }}},
      {"nozzle.throat_radius",
       {"m", [](EngineConfig& c, double v) { c.nozzle.throat_radius = v; c.nozzle.throat_area = 0.0; },
        [](const EngineConfig& c) { return c.nozzle.throat_radius; }}},
      {"nozzle.throat_area",
       {"m^2", [](EngineConfig& c, double v) { c.nozzle.throat_area = v; c.nozzle.throat_radius = 0.0; },
        [](const EngineConfig& c) { return c.nozzle.throat_area; }}},
      {"nozzle.expansion_ratio",
       {"-", [](EngineConfig& c, double v) { c.nozzle.expansion_ratio = v; },
        [](const EngineConfig& c) { return c.nozzle.expansion_ratio; }}},
      {"nozzle.contraction_ratio",
       {"-", [](EngineConfig& c, double v) { c.nozzle.contraction_ratio = v; c.nozzle.chamber_radius = 0.0; },
        [](const EngineConfig& c) { return c.nozzle.contraction_ratio; }}},
      {"nozzle.chamber_length",
       {"m", [](EngineConfig& c, double v) { c.nozzle.chamber_length = v; },
        [](const EngineConfig& c) { return c.nozzle.chamber_length; }}},
      {"nozzle.bell_length_fraction",
       {"-", [](EngineConfig& c, double v) { c.nozzle.bell_length_fraction = v; },
        [](const EngineConfig& c) { return c.nozzle.bell_length_fraction; }}},
      {"nozzle.bell_exit_angle",
       {"deg", [](EngineConfig& c, double v) { c.nozzle.bell_exit_angle = v; },
        [](const EngineConfig& c) { return c.nozzle.bell_exit_angle; }}},
      {"nozzle.cone_half_angle",
       {"deg", [](EngineConfig& c, double v) { c.nozzle.cone_half_angle = v; },
        [](const EngineConfig& c) { return c.nozzle.cone_half_angle; }}},
      {"performance.ambient_pressure",
       {"Pa", [](EngineConfig& c, double v) { c.ambient_pressure = v; c.use_altitude = false; },
        [](const EngineConfig& c) { return c.resolveAmbientPressure(); }}},
      {"performance.altitude",
       {"m", [](EngineConfig& c, double v) { c.altitude = v; c.use_altitude = true; },
        [](const EngineConfig& c) { return c.use_altitude ? c.altitude : 0.0; }}},
      {"cooling.channel_height",
       {"m", [](EngineConfig& c, double v) { c.cooling.channel_height = v; },
        [](const EngineConfig& c) { return c.cooling.channel_height; }}},
      {"cooling.width_fraction",
       {"-", [](EngineConfig& c, double v) { c.cooling.width_fraction = v; },
        [](const EngineConfig& c) { return c.cooling.width_fraction; }}},
      {"cooling.num_channels",
       {"-", [](EngineConfig& c, double v) { c.cooling.num_channels = static_cast<int>(std::lround(v)); },
        [](const EngineConfig& c) { return static_cast<double>(c.cooling.num_channels); }}},
      {"cooling.wall_thickness",
       {"m", [](EngineConfig& c, double v) { c.cooling.wall_thickness = v; },
        [](const EngineConfig& c) { return c.cooling.wall_thickness; }}},
      {"cooling.wall_conductivity",
       {"W/(m K)",
        [](EngineConfig& c, double v) { c.cooling.wall.conductivity_reference = v; },
        [](const EngineConfig& c) { return c.cooling.wall.conductivity_reference; }}},
      {"cooling.nusselt_multiplier",
       {"-", [](EngineConfig& c, double v) { c.cooling.nusselt_multiplier = v; },
        [](const EngineConfig& c) { return c.cooling.nusselt_multiplier; }}},
      {"cooling.bartz_multiplier",
       {"-", [](EngineConfig& c, double v) { c.cooling.bartz_multiplier = v; },
        [](const EngineConfig& c) { return c.cooling.bartz_multiplier; }}},
      {"cooling.inlet_temperature",
       {"K", [](EngineConfig& c, double v) { c.cooling.inlet_temperature = v; },
        [](const EngineConfig& c) { return c.cooling.inlet_temperature; }}},
      {"cooling.inlet_pressure",
       {"Pa", [](EngineConfig& c, double v) { c.cooling.inlet_pressure = v; },
        [](const EngineConfig& c) { return c.cooling.inlet_pressure; }}},
      {"cooling.coolant_fuel_fraction",
       {"-", [](EngineConfig& c, double v) { c.coolant_fuel_fraction = v; },
        [](const EngineConfig& c) { return c.coolant_fuel_fraction; }}},
      {"cooling.roughness",
       {"m", [](EngineConfig& c, double v) { c.cooling.roughness = v; },
        [](const EngineConfig& c) { return c.cooling.roughness; }}},
      {"feed.injector_stiffness",
       {"-", [](EngineConfig& c, double v) { c.injector_stiffness = v; },
        [](const EngineConfig& c) { return c.injector_stiffness; }}},
  };
  return table;
}

struct MetricEntry {
  const char* units;
  double (*get)(const SteadyEngineResult&);
};

double coolingValue(const SteadyEngineResult& r, double CoolingResult::*m) {
  if (!r.has_cooling)
    throw ConfigError("metric requires a cooling analysis, but 'cooling.enabled' is false");
  return r.cooling.*m;
}

const std::map<std::string, MetricEntry>& metricTable() {
  static const std::map<std::string, MetricEntry> table = {
      {"chamber.temperature", {"K", [](const SteadyEngineResult& r) { return r.chamber.state.T; }}},
      {"chamber.pressure", {"Pa", [](const SteadyEngineResult& r) { return r.chamber.state.p; }}},
      {"chamber.molar_mass", {"kg/mol", [](const SteadyEngineResult& r) { return r.chamber.state.M; }}},
      {"chamber.density", {"kg/m^3", [](const SteadyEngineResult& r) { return r.chamber.state.rho; }}},
      {"chamber.gamma_s", {"-", [](const SteadyEngineResult& r) { return r.chamber.state.gamma_s; }}},
      {"chamber.gamma_frozen", {"-", [](const SteadyEngineResult& r) { return r.chamber.state.gamma_frozen; }}},
      {"chamber.cp", {"J/(kg K)", [](const SteadyEngineResult& r) { return r.chamber.state.cp_eff; }}},
      {"chamber.sound_speed", {"m/s", [](const SteadyEngineResult& r) { return r.chamber.state.a; }}},
      {"chamber.c_star_ideal", {"m/s", [](const SteadyEngineResult& r) { return r.chamber.c_star_ideal; }}},
      {"chamber.c_star", {"m/s", [](const SteadyEngineResult& r) { return r.chamber.c_star; }}},
      {"chamber.equivalence_ratio", {"-", [](const SteadyEngineResult& r) { return r.chamber.equivalence_ratio; }}},
      {"chamber.throat_temperature", {"K", [](const SteadyEngineResult& r) { return r.chamber.throat_state.T; }}},
      {"chamber.element_residual", {"-", [](const SteadyEngineResult& r) { return r.chamber.diagnostics.element_residual_rel; }}},
      {"chamber.gibbs_residual", {"-", [](const SteadyEngineResult& r) { return r.chamber.diagnostics.gibbs_residual; }}},
      {"chamber.enthalpy_residual", {"-", [](const SteadyEngineResult& r) { return r.chamber.diagnostics.state_residual; }}},
      {"performance.thrust", {"N", [](const SteadyEngineResult& r) { return r.performance.thrust; }}},
      {"performance.thrust_ideal", {"N", [](const SteadyEngineResult& r) { return r.performance.thrust_ideal; }}},
      {"performance.thrust_momentum", {"N", [](const SteadyEngineResult& r) { return r.performance.thrust_momentum; }}},
      {"performance.thrust_pressure", {"N", [](const SteadyEngineResult& r) { return r.performance.thrust_pressure; }}},
      {"performance.isp", {"s", [](const SteadyEngineResult& r) { return r.performance.isp; }}},
      {"performance.isp_ideal", {"s", [](const SteadyEngineResult& r) { return r.performance.isp_ideal; }}},
      {"performance.isp_vacuum", {"s", [](const SteadyEngineResult& r) { return r.performance.isp_vacuum; }}},
      {"performance.cf", {"-", [](const SteadyEngineResult& r) { return r.performance.cf; }}},
      {"performance.cf_ideal", {"-", [](const SteadyEngineResult& r) { return r.performance.cf_ideal; }}},
      {"performance.c_effective", {"m/s", [](const SteadyEngineResult& r) { return r.performance.c_effective; }}},
      {"performance.mdot", {"kg/s", [](const SteadyEngineResult& r) { return r.mdot; }}},
      {"performance.exit_pressure", {"Pa", [](const SteadyEngineResult& r) { return r.performance.p_exit; }}},
      {"performance.exit_temperature", {"K", [](const SteadyEngineResult& r) { return r.performance.t_exit; }}},
      {"performance.exit_mach", {"-", [](const SteadyEngineResult& r) { return r.performance.mach_exit; }}},
      {"performance.exit_velocity", {"m/s", [](const SteadyEngineResult& r) { return r.performance.u_exit; }}},
      {"performance.pressure_ratio", {"-", [](const SteadyEngineResult& r) {
         return r.ambient_pressure > 0.0 ? r.performance.p_exit / r.ambient_pressure : 0.0; }}},
      {"performance.expansion_ratio", {"-", [](const SteadyEngineResult& r) { return r.performance.area_ratio; }}},
      {"performance.mass_flow_residual", {"-", [](const SteadyEngineResult& r) { return r.performance.mass_flow_residual; }}},
      {"performance.thrust_sea_level",
       {"N", [](const SteadyEngineResult& r) { return r.performance.thrust_sea_level; }}},
      {"performance.isp_sea_level",
       {"s", [](const SteadyEngineResult& r) { return r.performance.isp_sea_level; }}},
      {"performance.isp_ascent",
       {"s", [](const SteadyEngineResult& r) { return r.performance.isp_ascent; }}},
      {"performance.separation_margin",
       {"-", [](const SteadyEngineResult& r) { return r.performance.separation_margin; }}},
      {"performance.energy_residual", {"-", [](const SteadyEngineResult& r) { return r.performance.energy_residual; }}},
      {"geometry.l_star", {"m", [](const SteadyEngineResult& r) { return r.l_star; }}},
      {"geometry.residence_time", {"s", [](const SteadyEngineResult& r) { return r.residence_time; }}},
      {"geometry.throat_area", {"m^2", [](const SteadyEngineResult& r) { return r.throat_area; }}},
      {"geometry.exit_radius",
       {"m", [](const SteadyEngineResult& r) { return r.exit_radius; }}},
      {"geometry.exit_diameter",
       {"m", [](const SteadyEngineResult& r) { return r.exit_diameter; }}},
      {"geometry.total_length",
       {"m", [](const SteadyEngineResult& r) { return r.total_length; }}},
      {"geometry.divergent_length",
       {"m", [](const SteadyEngineResult& r) { return r.divergent_length; }}},
      {"geometry.exit_area", {"m^2", [](const SteadyEngineResult& r) { return r.exit_area; }}},
      {"cooling.max_wall_temperature", {"K", [](const SteadyEngineResult& r) { return coolingValue(r, &CoolingResult::max_wall_temperature); }}},
      {"cooling.max_heat_flux", {"W/m^2", [](const SteadyEngineResult& r) { return coolingValue(r, &CoolingResult::max_heat_flux); }}},
      {"cooling.total_heat_load", {"W", [](const SteadyEngineResult& r) { return coolingValue(r, &CoolingResult::total_heat_load); }}},
      {"cooling.pressure_drop", {"Pa", [](const SteadyEngineResult& r) { return coolingValue(r, &CoolingResult::coolant_pressure_drop); }}},
      {"cooling.outlet_temperature", {"K", [](const SteadyEngineResult& r) { return coolingValue(r, &CoolingResult::coolant_outlet_temperature); }}},
      {"cooling.temperature_rise", {"K", [](const SteadyEngineResult& r) { return coolingValue(r, &CoolingResult::coolant_temperature_rise); }}},
      {"cooling.energy_balance_residual", {"-", [](const SteadyEngineResult& r) { return coolingValue(r, &CoolingResult::energy_balance_residual); }}},
      {"cooling.flux_residual", {"-", [](const SteadyEngineResult& r) { return coolingValue(r, &CoolingResult::max_flux_residual); }}},
      {"feed.oxidizer_tank_pressure", {"Pa", [](const SteadyEngineResult& r) {
         if (!r.has_feed) throw ConfigError("metric requires a feed-system analysis");
         return r.feed.oxidizer.required_tank_pressure; }}},
      {"feed.fuel_tank_pressure", {"Pa", [](const SteadyEngineResult& r) {
         if (!r.has_feed) throw ConfigError("metric requires a feed-system analysis");
         return r.feed.fuel.required_tank_pressure; }}},
      {"feed.total_pump_power", {"W", [](const SteadyEngineResult& r) {
         if (!r.has_feed) throw ConfigError("metric requires a feed-system analysis");
         return r.feed.total_pump_power; }}},
  };
  return table;
}

[[noreturn]] void unknown(const std::string& what, const std::string& name,
                          const std::vector<std::string>& valid) {
  std::ostringstream os;
  os << "unknown " << what << " '" << name << "'. Valid names are:";
  for (const auto& v : valid) os << "\n  " << v;
  throw ConfigError(os.str());
}

}  // namespace

std::vector<std::string> parameterNames() {
  std::vector<std::string> out;
  for (const auto& kv : paramTable()) out.push_back(kv.first);
  return out;
}

std::vector<std::string> metricNames() {
  std::vector<std::string> out;
  for (const auto& kv : metricTable()) out.push_back(kv.first);
  return out;
}

std::string parameterUnits(const std::string& name) {
  const auto it = paramTable().find(name);
  if (it == paramTable().end()) unknown("parameter", name, parameterNames());
  return it->second.units;
}

std::string metricUnits(const std::string& name) {
  const auto it = metricTable().find(name);
  if (it == metricTable().end()) unknown("metric", name, metricNames());
  return it->second.units;
}

void applyParameter(EngineConfig& cfg, const std::string& name, double value) {
  const auto it = paramTable().find(name);
  if (it == paramTable().end()) unknown("parameter", name, parameterNames());
  if (!std::isfinite(value))
    throw ConfigError("parameter '" + name + "' was given a non-finite value");
  it->second.set(cfg, value);
}

double readParameter(const EngineConfig& cfg, const std::string& name) {
  const auto it = paramTable().find(name);
  if (it == paramTable().end()) unknown("parameter", name, parameterNames());
  return it->second.get(cfg);
}

double readMetric(const SteadyEngineResult& result, const std::string& name) {
  const auto it = metricTable().find(name);
  if (it == metricTable().end()) unknown("metric", name, metricNames());
  return it->second.get(result);
}

}  // namespace ignis
