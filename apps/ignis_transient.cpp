// SPDX-License-Identifier: MIT
/// \file ignis_transient.cpp
/// \brief Startup and shutdown transient simulation.

#include <chrono>
#include <iomanip>
#include <iostream>

#include "AppCommon.hpp"
#include "ignis/engine/SteadyEngine.hpp"
#include "ignis/transient/TransientChamber.hpp"

using namespace ignis;

int main(int argc, char** argv) {
  return app::run([&]() -> int {
    CommandLine cli;
    const std::vector<CliOption> options = {
        {"integrator", "NAME", "override the integrator: rk4 or rk45", "from the configuration"},
        {"dt", "S", "override the fixed (rk4) or initial (rk45) step", "from the configuration"},
        {"no-interp-check", "", "skip the table interpolation-error measurement", ""},
        {"export-table", "", "also export the equilibrium property table as CSV", ""},
    };
    if (!parseCommandLine(argc, argv, "ignis_transient",
                          "Integrate the zero-dimensional chamber mass and energy equations\n"
                          "through a startup or shutdown sequence, reporting the pressure,\n"
                          "temperature and mixture-ratio histories and the conservation\n"
                          "residuals.",
                          options, cli))
      return 0;

    EngineConfig cfg = EngineConfig::load(cli.config);
    app::applyCommandLine(cfg, cli);
    auto tc = parseTransient(cfg.root(), cfg);
    if (!tc.enabled)
      throw ConfigError("configuration '" + cli.config +
                        "' has no enabled 'transient:' section");
    if (flagPresent(cli, "integrator")) tc.spec.integrator = flagString(cli, "integrator", "rk45");
    if (flagPresent(cli, "dt")) tc.spec.dt = flagDouble(cli, "dt", tc.spec.dt);

    const SteadyEngine engine(cfg);
    const auto geom = NozzleGeometry::build(cfg.nozzle);
    if (!(tc.spec.chamber_volume > 0.0)) tc.spec.chamber_volume = geom.chamberVolume();
    if (!(tc.spec.throat_area > 0.0)) tc.spec.throat_area = geom.throatArea();
    if (!(tc.spec.exit_area > 0.0)) tc.spec.exit_area = geom.exitArea();

    const auto mix = engine.mixture(cfg);
    const auto& ox = engine.propellants().at(cfg.oxidizer);
    const auto& fu = engine.propellants().at(cfg.fuel);
    tc.spec.oxidizer_inlet_enthalpy =
        ox.molarEnthalpy(mix.oxidizerTemperature(), engine.database()) / ox.molar_mass;
    tc.spec.fuel_inlet_enthalpy =
        fu.molarEnthalpy(mix.fuelTemperature(), engine.database()) / fu.molar_mass;

    if (!cli.quiet) {
      std::cout << versionBanner() << "\n\ncase " << cfg.name << "\n"
                << "building the equilibrium property table ("
                << tc.grid.mr_points << " x " << tc.grid.t_points << " x " << tc.grid.p_points
                << " = " << tc.grid.mr_points * tc.grid.t_points * tc.grid.p_points
                << " equilibrium solves) ...\n"
                << std::flush;
    }
    const auto t0 = std::chrono::steady_clock::now();
    const auto table = EquilibriumTable::build(engine.solver(), mix, tc.grid);
    const auto t1 = std::chrono::steady_clock::now();
    const double build_s = std::chrono::duration<double>(t1 - t0).count();

    InterpolationError interp;
    if (!flagPresent(cli, "no-interp-check"))
      interp = table.measureError(engine.solver(), mix, tc.interpolation_samples,
                                  tc.interpolation_seed);

    const auto t2 = std::chrono::steady_clock::now();
    const auto res = simulateTransient(table, tc.spec);
    const auto t3 = std::chrono::steady_clock::now();
    const double run_s = std::chrono::duration<double>(t3 - t2).count();

    Table t("transient");
    std::vector<double> time, p, T, mr, m, mox, mf, mdot_o, mdot_f, mdot_out, cstar, gam, eta,
        thrust, molar, rho, choked, clamped;
    for (const auto& s : res.samples) {
      time.push_back(s.t);
      p.push_back(s.pressure);
      T.push_back(s.temperature);
      mr.push_back(s.mixture_ratio);
      m.push_back(s.mass);
      mox.push_back(s.mass_oxidizer);
      mf.push_back(s.mass_fuel);
      mdot_o.push_back(s.mdot_ox_in);
      mdot_f.push_back(s.mdot_fuel_in);
      mdot_out.push_back(s.mdot_out);
      cstar.push_back(s.c_star);
      gam.push_back(s.gamma_s);
      eta.push_back(s.eta_heat);
      thrust.push_back(s.thrust);
      molar.push_back(s.molar_mass);
      rho.push_back(s.density);
      choked.push_back(s.choked ? 1.0 : 0.0);
      clamped.push_back(s.inlet_mixture_ratio_clamped ? 1.0 : 0.0);
    }
    t.addColumn("time", time, "s");
    t.addColumn("pressure", p, "Pa");
    t.addColumn("temperature", T, "K");
    t.addColumn("mixture_ratio", mr, "-");
    t.addColumn("mass", m, "kg");
    t.addColumn("mass_oxidizer", mox, "kg");
    t.addColumn("mass_fuel", mf, "kg");
    t.addColumn("mdot_ox_in", mdot_o, "kg/s");
    t.addColumn("mdot_fuel_in", mdot_f, "kg/s");
    t.addColumn("mdot_out", mdot_out, "kg/s");
    t.addColumn("c_star", cstar, "m/s");
    t.addColumn("gamma_s", gam, "-");
    t.addColumn("eta_heat", eta, "-");
    t.addColumn("thrust", thrust, "N");
    t.addColumn("molar_mass", molar, "kg/mol");
    t.addColumn("density", rho, "kg/m^3");
    t.addColumn("choked", choked, "-");
    t.addColumn("inlet_mr_clamped", clamped, "-");
    t.writeCsv(app::outputPath(cfg, "transient.csv"));
    if (flagPresent(cli, "export-table"))
      table.exportCsv(app::outputPath(cfg, "equilibrium_table.csv"));

    Json j = app::provenance(cfg);
    Json tab = Json::object();
    tab["mr_points"] = Json(tc.grid.mr_points);
    tab["t_points"] = Json(tc.grid.t_points);
    tab["p_points"] = Json(tc.grid.p_points);
    tab["mr_range"] = Json::of(std::vector<double>{tc.grid.mr_min, tc.grid.mr_max});
    tab["t_range"] = Json::of(std::vector<double>{tc.grid.t_min, tc.grid.t_max});
    tab["p_range"] = Json::of(std::vector<double>{tc.grid.p_min, tc.grid.p_max});
    tab["build_seconds"] = Json(build_s);
    if (interp.samples > 0) {
      Json ie = Json::object();
      ie["samples"] = Json(interp.samples);
      ie["max_rel_temperature"] = Json(interp.max_rel_temperature);
      ie["rms_rel_temperature"] = Json(interp.rms_rel_temperature);
      ie["max_rel_pressure"] = Json(interp.max_rel_pressure);
      ie["rms_rel_pressure"] = Json(interp.rms_rel_pressure);
      ie["max_rel_molar_mass"] = Json(interp.max_rel_molar_mass);
      ie["max_rel_cstar"] = Json(interp.max_rel_cstar);
      tab["interpolation_error"] = ie;
    }
    j["equilibrium_table"] = tab;
    Json r = Json::object();
    r["integrator"] = Json(res.integrator);
    r["steps"] = Json(res.steps);
    r["rejected_steps"] = Json(res.rejected_steps);
    r["rhs_evaluations"] = Json(res.rhs_evaluations);
    r["wall_seconds"] = Json(run_s);
    r["max_pressure"] = Json(res.max_pressure);
    r["time_to_90_percent"] = Json(res.time_to_90_percent);
    r["mass_conservation_error"] = Json(res.mass_conservation_error);
    r["energy_conservation_error"] = Json(res.energy_conservation_error);
    r["min_mixture_ratio"] = Json(res.min_mixture_ratio);
    r["max_mixture_ratio"] = Json(res.max_mixture_ratio);
    r["min_temperature"] = Json(res.min_temperature);
    r["max_temperature"] = Json(res.max_temperature);
    r["inlet_clamped_samples"] = Json(res.inlet_clamped_samples);
    j["transient"] = r;
    j["history"] = t.toJson();
    j.writeFile(app::outputPath(cfg, "transient.json"));

    if (!cli.quiet) {
      std::cout << "table built in " << std::fixed << std::setprecision(2) << build_s << " s\n";
      if (interp.samples > 0) std::cout << interp.summary() << "\n";
      std::cout << "\n" << res.summary() << "\n  wall time              " << std::fixed
                << std::setprecision(3) << run_s << " s\n\n"
                << std::setw(10) << "t_ms" << std::setw(12) << "p_MPa" << std::setw(10) << "T_K"
                << std::setw(9) << "O/F" << std::setw(12) << "mdot_out" << std::setw(12)
                << "F_kN" << "\n";
      const std::size_t step = std::max<std::size_t>(1, res.samples.size() / 20);
      for (std::size_t i = 0; i < res.samples.size(); i += step) {
        const auto& s = res.samples[i];
        std::cout << std::setw(10) << std::setprecision(2) << s.t * 1e3 << std::setw(12)
                  << std::setprecision(5) << s.pressure * 1e-6 << std::setw(10)
                  << std::setprecision(1) << s.temperature << std::setw(9)
                  << std::setprecision(3) << s.mixture_ratio << std::setw(12)
                  << std::setprecision(3) << s.mdot_out << std::setw(12) << s.thrust * 1e-3
                  << "\n";
      }
      std::cout << "\nwrote " << app::outputPath(cfg, "transient.csv") << "\n";
    }
    return 0;
  });
}
