// SPDX-License-Identifier: MIT
#include "ignis/engine/SteadyEngine.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <memory>
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
  res.exit_radius = geom.exitRadius();
  res.exit_diameter = 2.0 * geom.exitRadius();
  res.total_length = geom.exitPosition();
  res.divergent_length = geom.divergentLength();
  res.chamber_volume = geom.chamberVolume();
  res.l_star = geom.characteristicLength();

  const auto flow = chamber.makeFlow(res.chamber);
  auto perf_opts = cfg.performance;
  perf_opts.eta_c_star = cfg.eta_c_star;
  perf_opts.viscous = {};
  res.performance = evaluateNozzle(flow, geom, p_ambient, perf_opts);

  // Mass flows.  The cooling solve below takes its coolant from these
  // inviscid flows; the boundary layer then moves the engine's flow by
  // (1 - C_d), a fraction of a percent, which the jacket does not see.
  auto setFlows = [&]() {
    res.mdot = res.performance.mdot;
    res.mdot_oxidizer = res.mdot * mix.oxidizerMassFraction();
    res.mdot_fuel = res.mdot * mix.fuelMassFraction();
  };
  setFlows();
  res.residence_time = res.chamber.residenceTime(res.chamber_volume, res.mdot);
  if (cfg.sample_profile) res.profile = sampleAxialProfile(flow, geom);

  CoolingSpec cs = cfg.cooling;
  // The cycle sees the jacket only when the fuel cools it, straight from the
  // fuel pump.
  const bool fuel_cooled =
      cfg.cooling_enabled && cs.coolant == propellants_->at(cfg.fuel).coolant_table;
  std::unique_ptr<CoolantFluid> jacket_fluid;
  // The jacket's share of the fuel.  The jacket is solved on the inviscid
  // flow; the cycle carries the same share of the final one, so that the
  // turbine never passes more fuel than the pump delivers.
  double jacket_fuel_share = 0.0;
  auto cycleInputs = [&]() {
    CycleInputs ci;
    ci.chamber_pressure = cfg.chamber_pressure;
    ci.oxidizer_mass_flow = res.mdot_oxidizer;
    ci.fuel_mass_flow = res.mdot_fuel;
    ci.oxidizer = &propellants_->at(cfg.oxidizer);
    ci.fuel = &propellants_->at(cfg.fuel);
    ci.oxidizer_temperature = mix.oxidizerTemperature();
    ci.fuel_temperature = mix.fuelTemperature();
    if (res.has_cooling && fuel_cooled) {
      ci.jacket = true;
      // A film drawn from the jacket outlet leaves before the turbine.
      const bool film_from_jacket = res.cooling.has_film && cs.film_temperature <= 0.0 &&
                                    (cs.film_coolant.empty() || cs.film_coolant == cs.coolant);
      ci.jacket_mass_flow = jacket_fuel_share * res.mdot_fuel *
                            (1.0 - (film_from_jacket ? res.cooling.film_mass_flow /
                                                           cs.coolant_mass_flow
                                                     : 0.0));
      ci.jacket_pressure_drop = res.cooling.coolant_pressure_drop;
      ci.jacket_outlet_temperature = res.cooling.coolant_outlet_temperature;
      ci.jacket_heat = res.cooling.total_heat_load;
      ci.jacket_fluid = jacket_fluid.get();
    }
    if (res.has_cooling && fuel_cooled && cfg.cooling.inlet_pressure > 0.0)
      ci.fuel_discharge_floor = cfg.cooling.inlet_pressure;
    ci.exit_pressure = res.performance.p_exit;
    ci.ambient_pressure = p_ambient;
    ci.thrust = res.performance.thrust;
    ci.thrust_vacuum = res.performance.isp_vacuum * constants::g0 * res.mdot;
    return ci;
  };

  if (cfg.cooling_enabled) {
    if (!(cs.coolant_mass_flow > 0.0))
      cs.coolant_mass_flow = cfg.coolant_fuel_fraction * res.mdot_fuel;
    if (!(cs.inlet_temperature > 0.0))
      cs.inlet_temperature = propellants_->at(cfg.fuel).reference_temperature;
    if (!(cs.inlet_pressure > 0.0)) cs.inlet_pressure = 1.6 * cfg.chamber_pressure;
    // Under a cycle an unconfigured inlet temperature is the fuel pump's
    // outlet; start from it, since a cryogen at its storage temperature can
    // be solid at a turbopump's discharge pressure.
    if (cfg.cycle_enabled && fuel_cooled && !(cfg.cooling.inlet_temperature > 0.0))
      cs.inlet_temperature = pumpDischargeTemperature(
          propellants_->at(cfg.fuel), mix.fuelTemperature(), cfg.cycle.pump_inlet_pressure_fuel,
          cs.inlet_pressure, cfg.cycle.pump_efficiency_fuel);
    if (cfg.film_fuel_fraction > 0.0) cs.film_mass_flow = cfg.film_fuel_fraction * res.mdot_fuel;
    jacket_fuel_share = cs.coolant_mass_flow / res.mdot_fuel;
    res.cooling = solveRegenerativeCooling(flow, geom, res.chamber, cs, *transport_);
    res.has_cooling = true;
    if (cfg.cycle_enabled && fuel_cooled)
      jacket_fluid = std::make_unique<CoolantFluid>(CoolantFluid::load(cs.coolant));
    // Unless they are configured, the jacket's inlet pressure and temperature
    // are the fuel pump's discharge, which the cycle decides and which depends
    // on the jacket's own pressure drop (and, in an expander, its outlet
    // temperature): iterate.
    const bool set_p = !(cfg.cooling.inlet_pressure > 0.0);
    const bool set_t = !(cfg.cooling.inlet_temperature > 0.0);
    if (cfg.cycle_enabled && fuel_cooled && (set_p || set_t)) {
      auto& cc = res.cycle_coupling;
      cc.jacket_inlet_from_cycle = set_p;
      for (int pass = 0; pass < 8; ++pass) {
        const auto trial = solveCycle(cfg.cycle, cycleInputs(), *solver_);
        cc.passes = pass + 1;
        if (!trial.feasible || !(trial.fuel_discharge_pressure > 0.0)) break;
        const double p_new = set_p ? trial.fuel_discharge_pressure : cs.inlet_pressure;
        const double t_new = set_t ? trial.fuel_pump.t_out : cs.inlet_temperature;
        if (std::abs(p_new - cs.inlet_pressure) <= 1.0e-3 * cs.inlet_pressure &&
            std::abs(t_new - cs.inlet_temperature) <= 0.05)
          break;
        cs.inlet_pressure = p_new;
        cs.inlet_temperature = t_new;
        res.cooling = solveRegenerativeCooling(flow, geom, res.chamber, cs, *transport_);
      }
    }
  }

  // --- boundary-layer losses ---------------------------------------------
  if (cfg.boundary_layer_losses) {
    // The jacket's hot-wall temperature where there is a jacket.
    std::vector<double> wx, wt;
    if (res.has_cooling)
      for (const auto& st : res.cooling.stations) { wx.push_back(st.x); wt.push_back(st.t_wall_hot); }
    auto wall = [&](double x) {
      if (wx.empty() || x < wx.front() || x > wx.back()) return 0.0;
      const auto it = std::upper_bound(wx.begin(), wx.end(), x);
      if (it == wx.begin()) return wt.front();
      if (it == wx.end()) return wt.back();
      const std::size_t i = static_cast<std::size_t>(it - wx.begin());
      const double t = (x - wx[i - 1]) / (wx[i] - wx[i - 1]);
      return wt[i - 1] + t * (wt[i] - wt[i - 1]);
    };
    BoundaryLayerLossOptions bo;
    bo.uncooled_wall_temperature = cfg.uncooled_wall_temperature;
    res.boundary_layer = computeBoundaryLayerLosses(flow, geom, res.chamber, *transport_, wall, bo);
    res.has_boundary_layer = true;
    perf_opts.viscous.apply = true;
    perf_opts.viscous.throat_displacement_thickness = res.boundary_layer.throat_displacement_thickness;
    perf_opts.viscous.exit_displacement_thickness = res.boundary_layer.exit_displacement_thickness;
    perf_opts.viscous.exit_momentum_thickness = res.boundary_layer.exit_momentum_thickness;
    res.performance = evaluateNozzle(flow, geom, p_ambient, perf_opts);
    setFlows();
    res.residence_time = res.chamber.residenceTime(res.chamber_volume, res.mdot);
  }

  // --- film-cooling performance bracket ----------------------------------
  if (res.has_cooling && res.cooling.has_film) {
    // The unmixed limit of the film's cost (see FilmPerformance).  The core
    // carries the same boundary-layer correction as the engine.
    auto& fp = res.film;
    fp.present = true;
    fp.film_mass_flow = res.cooling.film_mass_flow;
    if (!(fp.film_mass_flow < res.mdot_fuel))
      throw ConfigError("film cooling: the film cannot take the whole fuel flow");
    fp.film_fraction_of_total = fp.film_mass_flow / res.mdot;
    fp.core_mixture_ratio = res.mdot_oxidizer / (res.mdot_fuel - fp.film_mass_flow);
    EngineConfig core_cfg = cfg;
    core_cfg.mixture_ratio = fp.core_mixture_ratio;
    const auto core_chamber = chamber.solve(mixture(core_cfg), cfg.chamber_pressure, cfg.eta_c_star);
    const auto core_flow = chamber.makeFlow(core_chamber);
    const auto core_perf = evaluateNozzle(core_flow, geom, p_ambient, perf_opts);

    const CoolantFluid fluid =
        CoolantFluid::load(cs.film_coolant.empty() ? cs.coolant : cs.film_coolant);
    const double t_film = res.cooling.film_temperature;
    // cp at the table's lowest pressure: the nearest it has to the ideal gas.
    const double cp = fluid.at(t_film, fluid.pMin()).cp;
    const double r_gas = constants::R_universal / fluid.molarMass();
    const double p_exit = core_perf.p_exit;
    const double t_exit = t_film * std::pow(p_exit / cfg.chamber_pressure, r_gas / cp);
    fp.film_exhaust_velocity = std::sqrt(std::max(2.0 * cp * (t_film - t_exit), 0.0));
    const double rho_exit = p_exit / (r_gas * t_exit);
    const double film_isp_vac = (fp.film_exhaust_velocity +
                                 p_exit / (rho_exit * std::max(fp.film_exhaust_velocity, 1e-9))) /
                                constants::g0;
    const double f = fp.film_fraction_of_total;
    fp.isp_vacuum_mixed = res.performance.isp_vacuum;
    fp.isp_vacuum_unmixed = (1.0 - f) * core_perf.isp_vacuum + f * film_isp_vac;
    fp.isp_vacuum_penalty = 1.0 - fp.isp_vacuum_unmixed / fp.isp_vacuum_mixed;
  }

  // --- turbopump cycle ---------------------------------------------------
  if (cfg.cycle_enabled) {
    res.cycle = solveCycle(cfg.cycle, cycleInputs(), *solver_);
    res.has_cycle = true;
    auto& cc = res.cycle_coupling;
    if (res.has_cooling && fuel_cooled) {
      cc.jacket_inlet_pressure = cs.inlet_pressure;
      cc.mismatch = res.cycle.fuel_discharge_pressure / cs.inlet_pressure - 1.0;
      if (!cc.jacket_inlet_from_cycle && std::abs(cc.mismatch) > 0.05) {
        std::ostringstream os;
        os << std::fixed << std::setprecision(2) << "the jacket was solved at its configured "
           << cs.inlet_pressure * 1e-6 << " MPa inlet pressure, but the cycle's fuel pump delivers "
           << res.cycle.fuel_discharge_pressure * 1e-6
           << " MPa; leave cooling.inlet_pressure unset to let the cycle decide it";
        res.cycle.warnings.push_back(os.str());
      }
    } else if (res.has_cooling) {
      res.cycle.warnings.push_back("the jacket is not cooled by the fuel, so the cycle carries "
                                   "neither its pressure drop nor its heat");
    }
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
  if (has_boundary_layer) {
    const auto& bl = boundary_layer;
    os << "  boundary layer        C_d " << std::setprecision(5) << performance.discharge_coefficient
       << ";  delta* throat " << std::setprecision(3) << bl.throat_displacement_thickness * 1e3
       << " mm, exit " << bl.exit_displacement_thickness * 1e3 << " mm;  theta exit "
       << bl.exit_momentum_thickness * 1e3 << " mm\n"
       << "                        core expands to A/A* " << std::setprecision(3)
       << performance.effective_area_ratio << " (geometric " << bl.geometric_area_ratio << ")\n"
       << "  boundary-layer loss   " << std::setprecision(2)
       << (performance.thrust_inviscid - performance.thrust) * 1e-3 << " kN, "
       << performance.isp_vacuum_inviscid - performance.isp_vacuum << " s vacuum Isp ("
       << 100.0 * (1.0 - performance.isp_vacuum / performance.isp_vacuum_inviscid) << " %)\n";
  }
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
  if (film.present) {
    os << "\nfilm-cooling performance bracket (vacuum Isp)\n" << std::fixed
       << "  film                  " << std::setprecision(2) << 100.0 * film.film_fraction_of_total
       << " % of the total flow; core O/F " << std::setprecision(3) << film.core_mixture_ratio << "\n"
       << "  fully mixed and burnt " << std::setprecision(2) << film.isp_vacuum_mixed << " s\n"
       << "  unmixed (two-stream)  " << film.isp_vacuum_unmixed << " s  (at most "
       << 100.0 * film.isp_vacuum_penalty << " % lower)\n";
  }
  if (has_cycle) {
    os << "\n" << cycle.summary();
    if (cycle_coupling.jacket_inlet_from_cycle)
      os << "  jacket inlet          set by the cycle: " << std::fixed << std::setprecision(3)
         << cycle_coupling.jacket_inlet_pressure * 1e-6 << " MPa after "
         << cycle_coupling.passes << " pass(es), " << std::setprecision(3)
         << 100.0 * cycle_coupling.mismatch << " % from the final discharge\n";
  }
  if (has_feed) {
    os << "\n" << feed.summary() << "\n";
    if (has_cycle)
      os << "  (the feed block sizes a pressure-fed alternative; the turbopump cycle above is "
            "separate)\n";
  }
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
  auto col = [&](const char* col_name, double ThermalStation::*m, const char* units) {
    std::vector<double> v;
    v.reserve(cooling.stations.size());
    for (const auto& s : cooling.stations) v.push_back(s.*m);
    t.addColumn(col_name, v, units);
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
  col("channel_height", &ThermalStation::channel_height, "m");
  col("t_drive", &ThermalStation::t_drive, "K");
  col("film_effectiveness", &ThermalStation::film_effectiveness, "-");
  // Zero throughout when the hot-gas model is Bartz, which has no layer.
  col("bl_momentum_thickness", &ThermalStation::bl_momentum_thickness, "m");
  col("bl_enthalpy_thickness", &ThermalStation::bl_enthalpy_thickness, "m");
  col("bl_re_theta", &ThermalStation::bl_re_theta, "-");
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
  p["thrust_sea_level"] = Json(performance.thrust_sea_level);
  p["isp_sea_level"] = Json(performance.isp_sea_level);
  p["isp_ascent"] = Json(performance.isp_ascent);
  p["ascent_mean_ambient"] = Json(performance.ascent_mean_ambient);
  p["separation_margin"] = Json(performance.separation_margin);
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
  if (has_boundary_layer) {
    Json b = Json::object();
    b["discharge_coefficient"] = Json(performance.discharge_coefficient);
    b["throat_displacement_thickness"] = Json(boundary_layer.throat_displacement_thickness);
    b["throat_momentum_thickness"] = Json(boundary_layer.throat_momentum_thickness);
    b["exit_displacement_thickness"] = Json(boundary_layer.exit_displacement_thickness);
    b["exit_momentum_thickness"] = Json(boundary_layer.exit_momentum_thickness);
    b["exit_re_theta"] = Json(boundary_layer.exit_re_theta);
    b["effective_area_ratio"] = Json(performance.effective_area_ratio);
    b["geometric_area_ratio"] = Json(boundary_layer.geometric_area_ratio);
    b["momentum_deficit"] = Json(performance.momentum_deficit);
    b["mdot_inviscid"] = Json(performance.mdot_inviscid);
    b["thrust_inviscid"] = Json(performance.thrust_inviscid);
    b["isp_vacuum_inviscid"] = Json(performance.isp_vacuum_inviscid);
    b["thrust_loss"] = Json(performance.thrust_inviscid - performance.thrust);
    b["isp_vacuum_loss"] = Json(performance.isp_vacuum_inviscid - performance.isp_vacuum);
    b["cooled_fraction"] = Json(boundary_layer.cooled_fraction);
    b["warnings"] = Json::of(boundary_layer.warnings);
    p["boundary_layer"] = b;
  }
  j["performance"] = p;

  Json g = Json::object();
  g["chamber_volume"] = Json(chamber_volume);
  g["l_star"] = Json(l_star);
  g["exit_radius"] = Json(exit_radius);
  g["exit_diameter"] = Json(exit_diameter);
  g["total_length"] = Json(total_length);
  g["divergent_length"] = Json(divergent_length);
  g["residence_time"] = Json(residence_time);
  j["geometry"] = g;

  if (has_cooling) {
    Json t = Json::object();
    t["hot_gas_model"] = Json(cooling.hot_gas_model);
    t["wall_iterations"] = Json(cooling.wall_iterations);
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
    t["wall_limit_temperature"] = Json(cooling.wall_limit_temperature);
    t["wall_material"] = Json(cooling.wall_material);
    t["max_coolant_side_wall_temperature"] = Json(cooling.max_coolant_side_wall_temperature);
    t["max_coolant_side_wall_temperature_x"] = Json(cooling.max_coolant_side_wall_temperature_x);
    t["coolant_wall_limit_exceeded"] = Json(cooling.coolant_wall_limit_exceeded);
    t["warnings"] = Json::of(cooling.warnings);
    if (cooling.has_film) {
      Json fj = Json::object();
      fj["mass_flow"] = Json(cooling.film_mass_flow);
      fj["temperature"] = Json(cooling.film_temperature);
      fj["coolant_velocity"] = Json(cooling.film_coolant_velocity);
      fj["velocity_ratio"] = Json(cooling.film_velocity_ratio);
      fj["length_to_half"] = Json(cooling.film_length_half);
      fj["length_to_fifth"] = Json(cooling.film_length_fifth);
      if (film.present) {
        fj["fraction_of_total_flow"] = Json(film.film_fraction_of_total);
        fj["core_mixture_ratio"] = Json(film.core_mixture_ratio);
        fj["film_exhaust_velocity"] = Json(film.film_exhaust_velocity);
        fj["isp_vacuum_mixed"] = Json(film.isp_vacuum_mixed);
        fj["isp_vacuum_unmixed"] = Json(film.isp_vacuum_unmixed);
        fj["isp_vacuum_penalty_max"] = Json(film.isp_vacuum_penalty);
      }
      t["film"] = fj;
    }
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
  if (has_cycle) {
    const auto& cy = cycle;
    Json y = Json::object();
    y["type"] = Json(toString(cy.type));
    y["feasible"] = Json(cy.feasible);
    if (!cy.feasible) y["infeasibility"] = Json(cy.infeasibility);
    auto pumpJson = [](const PumpResult& pump, const std::string& model) {
      Json q = Json::object();
      q["mass_flow"] = Json(pump.mass_flow);
      q["inlet_pressure"] = Json(pump.p_in);
      q["discharge_pressure"] = Json(pump.p_out);
      q["inlet_temperature"] = Json(pump.t_in);
      q["outlet_temperature"] = Json(pump.t_out);
      q["isentropic_head"] = Json(pump.isentropic_head);
      q["efficiency"] = Json(pump.efficiency);
      q["power"] = Json(pump.power);
      q["real_fluid"] = Json(pump.real_fluid);
      if (!model.empty()) q["model"] = Json(model);
      return q;
    };
    y["oxidizer_pump"] = pumpJson(cy.oxidizer_pump, cy.oxidizer_pump_model);
    y["fuel_pump"] = pumpJson(cy.fuel_pump, cy.fuel_pump_model);
    if (cy.has_boost_pump) {
      y["boost_pump"] = pumpJson(cy.boost_pump, "");
      y["boost_propellant"] = Json(cy.boost_propellant);
    }
    if (cy.has_combustor) {
      Json gas = Json::object();
      gas["mixture_ratio"] = Json(cy.combustor.mixture_ratio);
      gas["temperature"] = Json(cy.combustor.temperature);
      gas["pressure"] = Json(cy.combustor_pressure);
      gas["molar_mass"] = Json(cy.combustor.molar_mass);
      gas["enthalpy_residual"] = Json(cy.combustor.enthalpy_residual);
      y[cy.type == CycleType::kGasGenerator ? "gas_generator" : "preburner"] = gas;
    }
    Json t = Json::object();
    t["mass_flow"] = Json(cy.turbine.mass_flow);
    t["inlet_temperature"] = Json(cy.turbine.t_in);
    t["outlet_temperature"] = Json(cy.turbine.t_out);
    t["inlet_pressure"] = Json(cy.turbine.p_in);
    t["outlet_pressure"] = Json(cy.turbine.p_out);
    t["pressure_ratio"] = Json(cy.turbine.pressure_ratio);
    t["efficiency"] = Json(cy.turbine.efficiency);
    t["power"] = Json(cy.turbine_power);
    y["turbine"] = t;
    y["pump_power"] = Json(cy.pump_power);
    y["power_balance_residual"] = Json(cy.power_balance_residual);
    if (cy.max_power_ratio > 0.0) y["max_power_ratio"] = Json(cy.max_power_ratio);
    y["oxidizer_discharge_pressure"] = Json(cy.oxidizer_discharge_pressure);
    y["fuel_discharge_pressure"] = Json(cy.fuel_discharge_pressure);
    if (cy.type == CycleType::kGasGenerator) {
      y["gas_generator_mass_flow"] = Json(cy.gg_mass_flow);
      y["gas_generator_flow_fraction"] = Json(cy.gg_flow_fraction);
      y["exhaust_velocity"] = Json(cy.exhaust_velocity);
      y["exhaust_exit_pressure"] = Json(cy.exhaust_exit_pressure);
      y["exhaust_thrust_vacuum"] = Json(cy.exhaust_thrust_vacuum);
    }
    if (cy.type == CycleType::kStagedCombustion) {
      y["preburner_major_flow"] = Json(cy.preburner_major_flow);
      y["preburner_minor_flow"] = Json(cy.preburner_minor_flow);
    }
    y["total_mass_flow"] = Json(cy.total_mass_flow);
    y["overall_mixture_ratio"] = Json(cy.overall_mixture_ratio);
    y["thrust_vacuum_delivered"] = Json(cy.thrust_vacuum_delivered);
    y["thrust_delivered"] = Json(cy.thrust_delivered);
    y["isp_vacuum_delivered"] = Json(cy.isp_vacuum_delivered);
    y["isp_delivered"] = Json(cy.isp_delivered);
    y["isp_vacuum_dumped"] = Json(cy.isp_vacuum_dumped);
    if (cycle_coupling.jacket_inlet_from_cycle || cycle_coupling.jacket_inlet_pressure > 0.0) {
      Json k = Json::object();
      k["jacket_inlet_from_cycle"] = Json(cycle_coupling.jacket_inlet_from_cycle);
      k["passes"] = Json(static_cast<double>(cycle_coupling.passes));
      k["jacket_inlet_pressure"] = Json(cycle_coupling.jacket_inlet_pressure);
      k["mismatch"] = Json(cycle_coupling.mismatch);
      y["jacket_coupling"] = k;
    }
    y["warnings"] = Json::of(cy.warnings);
    j["cycle"] = y;
  }
  return j;
}

}  // namespace ignis
