// SPDX-License-Identifier: MIT
/// \file ignis_nozzle.cpp
/// \brief Nozzle performance across ambient pressures or altitudes.

#include <cmath>
#include <iomanip>
#include <iostream>

#include "AppCommon.hpp"
#include "ignis/engine/SteadyEngine.hpp"
#include "ignis/nozzle/Atmosphere.hpp"

using namespace ignis;

int main(int argc, char** argv) {
  return app::run([&]() -> int {
    CommandLine cli;
    const std::vector<CliOption> options = {
        {"altitude-min", "M", "lowest geometric altitude", "0"},
        {"altitude-max", "M", "highest geometric altitude", "80000"},
        {"points", "N", "number of altitude points", "81"},
        {"vacuum", "", "append a true-vacuum point (ambient pressure zero)", ""},
        {"frozen", "", "use frozen composition instead of shifting equilibrium", ""},
        {"contour", "", "also export the contour and the axial flow profile", ""},
    };
    if (!parseCommandLine(argc, argv, "ignis_nozzle",
                          "Analyse a nozzle across a range of ambient pressures: thrust, specific\n"
                          "impulse, the momentum and pressure contributions, the exit-to-ambient\n"
                          "pressure ratio and the expansion regime at each altitude.",
                          options, cli))
      return 0;

    EngineConfig cfg = EngineConfig::load(cli.config);
    app::applyCommandLine(cfg, cli);
    cfg.cooling_enabled = false;   // the sweep is a pure nozzle study
    cfg.feed_enabled = false;
    if (flagPresent(cli, "frozen")) cfg.composition = CompositionModel::kFrozen;

    const double z0 = flagDouble(cli, "altitude-min", 0.0);
    const double z1 = flagDouble(cli, "altitude-max", 80000.0);
    const int n = flagInt(cli, "points", 81);
    if (n < 2) throw ConfigError("--points must be at least 2");

    const SteadyEngine engine(cfg);
    const auto geom = NozzleGeometry::build(cfg.nozzle);

    std::vector<double> alt, pa, thrust, thrust_i, tm, tp, isp, isp_i, cf, pe, pr, mach, sep;
    std::vector<double> shock, ue;
    std::vector<std::string> regime;
    std::vector<double> points;
    for (int i = 0; i < n; ++i) points.push_back(z0 + (z1 - z0) * i / (n - 1));

    double cstar_ideal = 0.0, mdot = 0.0;
    for (double z : points) {
      const double p_amb = Atmosphere::at(z).pressure;
      const auto res = engine.runAt(p_amb);
      cstar_ideal = res.chamber.c_star_ideal;
      mdot = res.mdot;
      alt.push_back(z);
      pa.push_back(p_amb);
      thrust.push_back(res.performance.thrust);
      thrust_i.push_back(res.performance.thrust_ideal);
      tm.push_back(res.performance.thrust_momentum);
      tp.push_back(res.performance.thrust_pressure);
      isp.push_back(res.performance.isp);
      isp_i.push_back(res.performance.isp_ideal);
      cf.push_back(res.performance.cf);
      pe.push_back(res.performance.p_exit);
      pr.push_back(p_amb > 0.0 ? res.performance.p_exit / p_amb : std::nan(""));
      mach.push_back(res.performance.mach_exit);
      ue.push_back(res.performance.u_exit);
      sep.push_back(res.performance.separation_predicted ? 1.0 : 0.0);
      shock.push_back(res.performance.shock_in_nozzle ? 1.0 : 0.0);
      regime.push_back(toString(res.performance.regime));
    }
    if (flagPresent(cli, "vacuum")) {
      const auto res = engine.runAt(0.0);
      alt.push_back(std::nan(""));
      pa.push_back(0.0);
      thrust.push_back(res.performance.thrust);
      thrust_i.push_back(res.performance.thrust_ideal);
      tm.push_back(res.performance.thrust_momentum);
      tp.push_back(res.performance.thrust_pressure);
      isp.push_back(res.performance.isp);
      isp_i.push_back(res.performance.isp_ideal);
      cf.push_back(res.performance.cf);
      pe.push_back(res.performance.p_exit);
      pr.push_back(std::nan(""));
      mach.push_back(res.performance.mach_exit);
      ue.push_back(res.performance.u_exit);
      sep.push_back(0.0);
      shock.push_back(0.0);
      regime.push_back(toString(res.performance.regime));
    }

    Table t("altitude_performance");
    t.addColumn("altitude", alt, "m");
    t.addColumn("ambient_pressure", pa, "Pa");
    t.addColumn("thrust", thrust, "N");
    t.addColumn("thrust_ideal", thrust_i, "N");
    t.addColumn("thrust_momentum", tm, "N");
    t.addColumn("thrust_pressure", tp, "N");
    t.addColumn("isp", isp, "s");
    t.addColumn("isp_ideal", isp_i, "s");
    t.addColumn("cf", cf, "-");
    t.addColumn("exit_pressure", pe, "Pa");
    t.addColumn("exit_ambient_ratio", pr, "-");
    t.addColumn("exit_mach", mach, "-");
    t.addColumn("exit_velocity", ue, "m/s");
    t.addColumn("separation_predicted", sep, "-");
    t.addColumn("shock_in_nozzle", shock, "-");
    t.addColumn("regime", regime);
    t.writeCsv(app::outputPath(cfg, "altitude.csv"));

    if (flagPresent(cli, "contour")) {
      const auto res = engine.run();
      res.profileTable().writeCsv(app::outputPath(cfg, "profile.csv"));
    }

    Json j = app::provenance(cfg);
    j["composition_model"] = Json(toString(cfg.composition));
    j["expansion_ratio"] = Json(geom.expansionRatio());
    j["throat_area"] = Json(geom.throatArea());
    j["c_star_ideal"] = Json(cstar_ideal);
    j["mdot"] = Json(mdot);
    j["altitude_sweep"] = t.toJson();
    j.writeFile(app::outputPath(cfg, "altitude.json"));

    if (!cli.quiet) {
      std::cout << versionBanner() << "\n\ncase " << cfg.name << ": " << toString(cfg.composition)
                << " expansion, Ae/At = " << std::fixed << std::setprecision(2)
                << geom.expansionRatio() << ", mdot = " << std::setprecision(3) << mdot
                << " kg/s\n\n"
                << std::setw(10) << "alt_km" << std::setw(12) << "p_amb_kPa" << std::setw(12)
                << "F_kN" << std::setw(10) << "Isp_s" << std::setw(10) << "Cf"
                << std::setw(12) << "pe/pa" << "  regime\n";
      const std::size_t step = std::max<std::size_t>(1, alt.size() / 12);
      for (std::size_t i = 0; i < alt.size(); i += step)
        std::cout << std::setw(10) << std::setprecision(2) << alt[i] * 1e-3 << std::setw(12)
                  << std::setprecision(4) << pa[i] * 1e-3 << std::setw(12) << std::setprecision(3)
                  << thrust[i] * 1e-3 << std::setw(10) << std::setprecision(2) << isp[i]
                  << std::setw(10) << std::setprecision(4) << cf[i] << std::setw(12) << pr[i]
                  << "  " << regime[i] << (sep[i] > 0.5 ? "  [separation predicted]" : "")
                  << (shock[i] > 0.5 ? "  [internal shock]" : "") << "\n";
      std::cout << "\nwrote " << app::outputPath(cfg, "altitude.csv") << "\n";
    }
    return 0;
  });
}
