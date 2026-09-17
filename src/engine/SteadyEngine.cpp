// SPDX-License-Identifier: MIT
#include "ignis/engine/SteadyEngine.hpp"

#include <cmath>
#include <iomanip>
#include <sstream>

namespace ignis {

SteadyEngine::SteadyEngine(EngineConfig config) : config_(std::move(config)) {
  const std::string db_path = config_.species_database.empty() ? findDefaultSpeciesDatabase()
                                                               : config_.species_database;
  auto full = std::make_shared<SpeciesDatabase>(SpeciesDatabase::loadYaml(db_path));
  if (!config_.species_include.empty()) {
    db_ = std::make_shared<SpeciesDatabase>(full->subset(config_.species_include));
  } else {
    db_ = std::make_shared<SpeciesDatabase>(full->restrictToElements(config_.elements));
  }
  solver_ = std::make_shared<EquilibriumSolver>(*db_);
  transport_ = std::make_shared<TransportModel>(*db_);
  const std::string prop_path = config_.propellant_library.empty()
                                    ? findDefaultPropellantLibrary()
                                    : config_.propellant_library;
  propellants_ = std::make_shared<PropellantLibrary>(PropellantLibrary::loadYaml(prop_path));
}

PropellantMixture SteadyEngine::mixture(const EngineConfig& cfg) const {
  const auto& ox = propellants_->at(cfg.oxidizer);
  const auto& fu = propellants_->at(cfg.fuel);
  const double t_ox = cfg.oxidizer_temperature > 0.0 ? cfg.oxidizer_temperature
                                                     : ox.reference_temperature;
  const double t_fu = cfg.fuel_temperature > 0.0 ? cfg.fuel_temperature
                                                 : fu.reference_temperature;
  return PropellantMixture(ox, fu, cfg.mixture_ratio, t_ox, t_fu);
}

SteadyEngineResult SteadyEngine::run() const {
  return analyse(config_, config_.resolveAmbientPressure());
}

SteadyEngineResult SteadyEngine::runAt(double ambient_pressure) const {
  return analyse(config_, ambient_pressure);
}

SteadyEngineResult SteadyEngine::runWith(const EngineConfig& cfg) const {
  return analyse(cfg, cfg.resolveAmbientPressure());
}

SteadyEngineResult SteadyEngine::analyse(const EngineConfig& cfg, double p_ambient) const {
  SteadyEngineResult res;
  res.name = cfg.name;
  res.ambient_pressure = p_ambient;
  if (cfg.use_altitude) {
    res.altitude = cfg.altitude;
    res.altitude_known = true;
  } else if (p_ambient > 0.0) {
    try {
      res.altitude = Atmosphere::altitudeForPressure(p_ambient);
      res.altitude_known = true;
    } catch (const IgnisError&) {
      res.altitude_known = false;
    }
  }

  const auto mix = mixture(cfg);
  const CombustionChamber chamber(*solver_, cfg.composition);
  res.chamber = chamber.solve(mix, cfg.chamber_pressure, cfg.eta_c_star);

  const auto geom = NozzleGeometry::build(cfg.nozzle);
  res.throat_area = geom.throatArea();
  res.exit_area = geom.exitArea();
  res.chamber_volume = geom.chamberVolume();
  res.l_star = geom.characteristicLength();

  const auto flow = chamber.makeFlow(res.chamber);
  auto perf_opts = cfg.performance;
  perf_opts.eta_c_star = cfg.eta_c_star;
  res.performance = evaluateNozzle(flow, geom, p_ambient, perf_opts);

  res.mdot = res.performance.mdot;
  res.mdot_oxidizer = res.mdot * mix.oxidizerMassFraction();
  res.mdot_fuel = res.mdot * mix.fuelMassFraction();
  res.residence_time = res.chamber.residenceTime(res.chamber_volume, res.mdot);
  if (cfg.sample_profile) res.profile = sampleAxialProfile(flow, geom);

  if (cfg.cooling_enabled) {
    CoolingSpec cs = cfg.cooling;
    if (!(cs.coolant_mass_flow > 0.0))
      cs.coolant_mass_flow = cfg.coolant_fuel_fraction * res.mdot_fuel;
    if (!(cs.inlet_temperature > 0.0))
      cs.inlet_temperature = propellants_->at(cfg.fuel).reference_temperature;
    if (!(cs.inlet_pressure > 0.0)) cs.inlet_pressure = 1.6 * cfg.chamber_pressure;
    res.cooling = solveRegenerativeCooling(flow, geom, res.chamber, cs, *transport_);
    res.has_cooling = true;
  }

  if (cfg.feed_enabled) {
    FeedSystemSpec fs = cfg.feed;
    fs.oxidizer.mass_flow = res.mdot_oxidizer;
    fs.fuel.mass_flow = res.mdot_fuel;
    fs.oxidizer.density = propellants_->at(cfg.oxidizer).density;
    fs.fuel.density = propellants_->at(cfg.fuel).density;
    if (cfg.feed_capability_mode) {
      res.feed = evaluateFeedSystem(fs, cfg.chamber_pressure);
    } else {
      res.feed = sizeFeedSystem(fs, cfg.chamber_pressure, cfg.injector_stiffness);
    }
    res.has_feed = true;
  }
  return res;
}

std::string SteadyEngineResult::summary() const {
  std::ostringstream os;
  os << "=== " << name << " ===\n\n" << chamber.summary() << "\n\n";
  os << std::fixed;
  os << "nozzle performance (" << performance.composition_model << " expansion)\n"
     << "  ambient pressure      " << std::setprecision(1) << ambient_pressure << " Pa";
  if (altitude_known) os << "  (" << std::setprecision(0) << altitude << " m)";
  os << "\n"
     << "  expansion ratio       " << std::setprecision(3) << performance.area_ratio << "\n"
     << "  mass flow             " << std::setprecision(4) << mdot << " kg/s  (ox "
     << mdot_oxidizer << ", fuel " << mdot_fuel << ")\n"
     << "  exit Mach             " << std::setprecision(4) << performance.mach_exit << "\n"
     << "  exit pressure         " << std::setprecision(1) << performance.p_exit << " Pa\n"
     << "  exit temperature      " << std::setprecision(2) << performance.t_exit << " K\n"
     << "  exit velocity         " << performance.u_exit << " m/s\n"
     << "  momentum thrust       " << std::setprecision(1) << performance.thrust_momentum * 1e-3
     << " kN\n"
     << "  pressure thrust       " << performance.thrust_pressure * 1e-3 << " kN\n"
     << "  thrust (ideal)        " << performance.thrust_ideal * 1e-3 << " kN\n"
     << "  thrust (corrected)    " << performance.thrust * 1e-3 << " kN\n"
     << "  Isp (ideal)           " << std::setprecision(2) << performance.isp_ideal << " s\n"
     << "  Isp (corrected)       " << performance.isp << " s\n"
     << "  Isp (vacuum)          " << performance.isp_vacuum << " s\n"
     << "  effective velocity    " << std::setprecision(1) << performance.c_effective << " m/s\n"
     << "  Cf (ideal / corr.)    " << std::setprecision(4) << performance.cf_ideal << " / "
     << performance.cf << "\n"
     << "  loss factors          eta_c* " << performance.eta_c_star << ", lambda_div "
     << performance.lambda_divergence << ", eta_nozzle " << performance.eta_nozzle << "\n"
     << "  regime                " << toString(performance.regime) << "\n";
  if (performance.shock_in_nozzle)
    os << "  internal shock        at A/At = " << std::setprecision(3)
       << performance.shock_area_ratio << " (x = " << std::setprecision(1)
       << performance.shock_x * 1e3 << " mm)\n";
  if (performance.shock_unresolved) os << "  internal shock        UNRESOLVED\n";
  if (!performance.shock_note.empty()) os << "  note                  " << performance.shock_note << "\n";
  if (performance.separation_predicted)
    os << "  separation (" << toString(performance.separation_criterion) << ")  predicted at A/At = "
       << std::setprecision(3) << performance.separation_area_ratio << " (x = "
       << std::setprecision(1) << performance.separation_x * 1e3 << " mm)\n";
  os << "  chamber volume        " << std::setprecision(2) << chamber_volume * 1e6 << " cm^3\n"
     << "  characteristic L*     " << std::setprecision(4) << l_star << " m\n"
     << "  residence time        " << std::setprecision(3) << residence_time * 1e3 << " ms\n"
     << std::scientific << std::setprecision(2)
     << "  mass-flow residual    " << performance.mass_flow_residual << "\n"
     << "  energy residual       " << performance.energy_residual << "\n";
  if (has_cooling) os << "\n" << cooling.summary() << "\n";
  if (has_feed) os << "\n" << feed.summary() << "\n";
  return os.str();
}

Table SteadyEngineResult::profileTable() const {
  Table t("axial_profile");
  t.addColumn("x", profile.x, "m");
  t.addColumn("radius", profile.radius, "m");
  t.addColumn("area_ratio", profile.area_ratio, "-");
  t.addColumn("pressure", profile.p, "Pa");
  t.addColumn("temperature", profile.T, "K");
  t.addColumn("density", profile.rho, "kg/m^3");
  t.addColumn("velocity", profile.u, "m/s");
  t.addColumn("mach", profile.mach, "-");
  t.addColumn("sound_speed", profile.a, "m/s");
  t.addColumn("gamma_s", profile.gamma_s, "-");
  t.addColumn("cp", profile.cp, "J/(kg K)");
  t.addColumn("molar_mass", profile.molar_mass, "kg/mol");
  std::vector<double> sup(profile.supersonic.begin(), profile.supersonic.end());
  t.addColumn("supersonic", sup, "-");
  return t;
}

Table SteadyEngineResult::coolingTable() const {
  Table t("thermal_stations");
  if (!has_cooling) return t;
  auto col = [&](const char* name, double ThermalStation::*m, const char* units) {
    std::vector<double> v;
    v.reserve(cooling.stations.size());
    for (const auto& s : cooling.stations) v.push_back(s.*m);
    t.addColumn(name, v, units);
  };
  col("x", &ThermalStation::x, "m");
  col("radius", &ThermalStation::radius, "m");
  col("area_ratio", &ThermalStation::area_ratio, "-");
  col("gas_T", &ThermalStation::gas_T, "K");
  col("gas_p", &ThermalStation::gas_p, "Pa");
  col("mach", &ThermalStation::gas_mach, "-");
  col("t_adiabatic_wall", &ThermalStation::t_adiabatic_wall, "K");
  col("h_gas", &ThermalStation::h_gas, "W/(m^2 K)");
  col("q_convective", &ThermalStation::q_convective, "W/m^2");
  col("q_radiative", &ThermalStation::q_radiative, "W/m^2");
  col("q_total", &ThermalStation::q_total, "W/m^2");
  col("t_wall_hot", &ThermalStation::t_wall_hot, "K");
  col("t_wall_cold", &ThermalStation::t_wall_cold, "K");
  col("wall_conductivity", &ThermalStation::wall_conductivity, "W/(m K)");
  col("coolant_T", &ThermalStation::coolant_T, "K");
  col("coolant_p", &ThermalStation::coolant_p, "Pa");
  col("coolant_rho", &ThermalStation::coolant_rho, "kg/m^3");
  col("coolant_cp", &ThermalStation::coolant_cp, "J/(kg K)");
  col("coolant_velocity", &ThermalStation::coolant_velocity, "m/s");
  col("reynolds", &ThermalStation::reynolds, "-");
  col("prandtl", &ThermalStation::prandtl, "-");
  col("nusselt", &ThermalStation::nusselt, "-");
  col("h_coolant", &ThermalStation::h_coolant, "W/(m^2 K)");
  col("fin_efficiency", &ThermalStation::fin_efficiency, "-");
  col("friction_factor", &ThermalStation::friction_factor, "-");
  col("channel_width", &ThermalStation::channel_width, "m");
  col("land_width", &ThermalStation::land_width, "m");
  col("hydraulic_diameter", &ThermalStation::hydraulic_diameter, "m");
  return t;
}

Json SteadyEngineResult::toJson() const {
  Json j = Json::object();
  j["name"] = Json(name);
  Json c = Json::object();
  c["mixture_ratio"] = Json(chamber.mixture_ratio);
  c["equivalence_ratio"] = Json(chamber.equivalence_ratio);
  c["stoichiometric_mixture_ratio"] = Json(chamber.stoichiometric_mixture_ratio);
  c["reactant_enthalpy"] = Json(chamber.reactant_enthalpy);
  c["pressure"] = Json(chamber.state.p);
  c["temperature"] = Json(chamber.state.T);
  c["molar_mass"] = Json(chamber.state.M);
  c["gas_constant"] = Json(chamber.state.R);
  c["density"] = Json(chamber.state.rho);
  c["cp_frozen"] = Json(chamber.state.cp_frozen);
  c["cp_equilibrium"] = Json(chamber.state.cp_eff);
  c["cv_frozen"] = Json(chamber.state.cv_frozen);
  c["cv_equilibrium"] = Json(chamber.state.cv_eff);
  c["gamma_frozen"] = Json(chamber.state.gamma_frozen);
  c["gamma_s"] = Json(chamber.state.gamma_s);
  c["sound_speed"] = Json(chamber.state.a);
  c["entropy"] = Json(chamber.state.s);
  c["c_star_ideal"] = Json(chamber.c_star_ideal);
  c["c_star"] = Json(chamber.c_star);
  c["eta_c_star"] = Json(chamber.eta_c_star);
  c["throat_mass_flux"] = Json(chamber.throat_mass_flux);
  c["throat_temperature"] = Json(chamber.throat_state.T);
  c["throat_pressure"] = Json(chamber.throat_state.p);
  Json d = Json::object();
  d["iterations"] = Json(chamber.diagnostics.iterations);
  d["element_residual_rel"] = Json(chamber.diagnostics.element_residual_rel);
  d["gibbs_residual"] = Json(chamber.diagnostics.gibbs_residual);
  d["enthalpy_residual"] = Json(chamber.diagnostics.state_residual);
  d["mass_residual"] = Json(chamber.diagnostics.mass_residual);
  c["diagnostics"] = d;
  j["chamber"] = c;

  Json p = Json::object();
  p["composition_model"] = Json(performance.composition_model);
  p["ambient_pressure"] = Json(ambient_pressure);
  if (altitude_known) p["altitude"] = Json(altitude);
  p["area_ratio"] = Json(performance.area_ratio);
  p["throat_area"] = Json(throat_area);
  p["exit_area"] = Json(exit_area);
  p["mdot"] = Json(mdot);
  p["mdot_oxidizer"] = Json(mdot_oxidizer);
  p["mdot_fuel"] = Json(mdot_fuel);
  p["exit_mach"] = Json(performance.mach_exit);
  p["exit_pressure"] = Json(performance.p_exit);
  p["exit_temperature"] = Json(performance.t_exit);
  p["exit_velocity"] = Json(performance.u_exit);
  p["thrust_momentum"] = Json(performance.thrust_momentum);
  p["thrust_pressure"] = Json(performance.thrust_pressure);
  p["thrust_ideal"] = Json(performance.thrust_ideal);
  p["thrust"] = Json(performance.thrust);
  p["isp_ideal"] = Json(performance.isp_ideal);
  p["isp"] = Json(performance.isp);
  p["isp_vacuum"] = Json(performance.isp_vacuum);
  p["c_effective"] = Json(performance.c_effective);
  p["cf_ideal"] = Json(performance.cf_ideal);
  p["cf"] = Json(performance.cf);
  p["lambda_divergence"] = Json(performance.lambda_divergence);
  p["eta_nozzle"] = Json(performance.eta_nozzle);
  p["regime"] = Json(toString(performance.regime));
  p["shock_in_nozzle"] = Json(performance.shock_in_nozzle);
  if (performance.shock_in_nozzle) {
    p["shock_area_ratio"] = Json(performance.shock_area_ratio);
    p["shock_x"] = Json(performance.shock_x);
  }
  p["shock_unresolved"] = Json(performance.shock_unresolved);
  if (!performance.shock_note.empty()) p["shock_note"] = Json(performance.shock_note);
  p["separation_predicted"] = Json(performance.separation_predicted);
  if (performance.separation_predicted) {
    p["separation_criterion"] = Json(toString(performance.separation_criterion));
    p["separation_area_ratio"] = Json(performance.separation_area_ratio);
    p["separation_x"] = Json(performance.separation_x);
  }
  p["mass_flow_residual"] = Json(performance.mass_flow_residual);
  p["energy_residual"] = Json(performance.energy_residual);
  j["performance"] = p;

  Json g = Json::object();
  g["chamber_volume"] = Json(chamber_volume);
  g["l_star"] = Json(l_star);
  g["residence_time"] = Json(residence_time);
  j["geometry"] = g;

  if (has_cooling) {
    Json t = Json::object();
    t["max_wall_temperature"] = Json(cooling.max_wall_temperature);
    t["max_wall_temperature_x"] = Json(cooling.max_wall_temperature_x);
    t["max_heat_flux"] = Json(cooling.max_heat_flux);
    t["max_heat_flux_x"] = Json(cooling.max_heat_flux_x);
    t["total_heat_load"] = Json(cooling.total_heat_load);
    t["coolant_inlet_temperature"] = Json(cooling.coolant_inlet_temperature);
    t["coolant_outlet_temperature"] = Json(cooling.coolant_outlet_temperature);
    t["coolant_temperature_rise"] = Json(cooling.coolant_temperature_rise);
    t["coolant_pressure_drop"] = Json(cooling.coolant_pressure_drop);
    t["coolant_outlet_pressure"] = Json(cooling.coolant_outlet_pressure);
    t["cooled_area"] = Json(cooling.cooled_area);
    t["energy_balance_residual"] = Json(cooling.energy_balance_residual);
    t["max_flux_residual"] = Json(cooling.max_flux_residual);
    t["boiling_detected"] = Json(cooling.boiling_detected);
    t["wall_limit_exceeded"] = Json(cooling.wall_limit_exceeded);
    t["warnings"] = Json::of(cooling.warnings);
    j["cooling"] = t;
  }
  if (has_feed) {
    Json f = Json::object();
    auto leg = [](const FeedLegResult& r) {
      Json l = Json::object();
      l["mass_flow"] = Json(r.mass_flow);
      l["line_velocity"] = Json(r.line_velocity);
      l["reynolds"] = Json(r.reynolds);
      l["friction_factor"] = Json(r.friction_factor);
      l["dp_lines"] = Json(r.dp_lines);
      l["dp_dynamic"] = Json(r.dp_dynamic);
      l["dp_injector"] = Json(r.dp_injector);
      l["injector_area"] = Json(r.injector_area);
      l["injector_velocity"] = Json(r.injector_velocity);
      l["tank_pressure"] = Json(r.required_tank_pressure);
      l["stiffness"] = Json(r.stiffness);
      l["pump_power"] = Json(r.pump_power);
      return l;
    };
    f["oxidizer"] = leg(feed.oxidizer);
    f["fuel"] = leg(feed.fuel);
    f["total_pump_power"] = Json(feed.total_pump_power);
    f["stiffness_ok"] = Json(feed.stiffness_ok);
    if (!feed.notes.empty()) f["notes"] = Json(feed.notes);
    j["feed"] = f;
  }
  return j;
}

}  // namespace ignis
