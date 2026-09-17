// SPDX-License-Identifier: MIT
/// \file ignis_engine.cpp
/// \brief One complete steady engine analysis.

#include <iostream>

#include "AppCommon.hpp"
#include "ignis/engine/SteadyEngine.hpp"

using namespace ignis;

int main(int argc, char** argv) {
  return app::run([&]() -> int {
    CommandLine cli;
    const std::vector<CliOption> options = {
        {"altitude", "M", "evaluate at this geometric altitude instead of the configured ambient", ""},
        {"ambient", "PA", "evaluate at this ambient pressure instead of the configured one", ""},
        {"no-cooling", "", "skip the regenerative-cooling analysis", ""},
        {"no-feed", "", "skip the feed-system analysis", ""},
        {"contour", "", "also export the bare nozzle contour", ""},
    };
    if (!parseCommandLine(argc, argv, "ignis_engine",
                          "Run one complete steady-state engine analysis: chamber equilibrium,\n"
                          "quasi-1D nozzle expansion, thrust and specific impulse, regenerative\n"
                          "cooling and feed-system sizing.",
                          options, cli))
      return 0;

    EngineConfig cfg = EngineConfig::load(cli.config);
    app::applyCommandLine(cfg, cli);
    if (flagPresent(cli, "altitude")) {
      cfg.use_altitude = true;
      cfg.altitude = flagDouble(cli, "altitude", 0.0);
    }
    if (flagPresent(cli, "ambient")) {
      cfg.use_altitude = false;
      cfg.ambient_pressure = flagDouble(cli, "ambient", 0.0);
    }
    if (flagPresent(cli, "no-cooling")) cfg.cooling_enabled = false;
    if (flagPresent(cli, "no-feed")) cfg.feed_enabled = false;

    const SteadyEngine engine(cfg);
    const auto geom = NozzleGeometry::build(cfg.nozzle);
    const auto res = engine.run();

    if (!cli.quiet) {
      std::cout << versionBanner() << "\n\n" << geom.summary() << "\n\n" << res.summary();
    }

    res.profileTable().writeCsv(app::outputPath(cfg, "profile.csv"));
    if (res.has_cooling) res.coolingTable().writeCsv(app::outputPath(cfg, "thermal.csv"));
    if (flagPresent(cli, "contour")) {
      Table c("contour");
      std::vector<double> x, r, a;
      for (const auto& s : geom.stations()) { x.push_back(s.x); r.push_back(s.r); a.push_back(s.area_ratio); }
      c.addColumn("x", x, "m");
      c.addColumn("radius", r, "m");
      c.addColumn("area_ratio", a, "-");
      c.writeCsv(app::outputPath(cfg, "contour.csv"));
    }

    Json j = app::provenance(cfg);
    j["result"] = res.toJson();
    Json g = Json::object();
    g["throat_radius"] = Json(geom.throatRadius());
    g["chamber_radius"] = Json(geom.chamberRadius());
    g["exit_radius"] = Json(geom.exitRadius());
    g["throat_position"] = Json(geom.throatPosition());
    g["exit_position"] = Json(geom.exitPosition());
    g["divergent_length"] = Json(geom.divergentLength());
    g["reference_conical_length"] = Json(geom.referenceConicalLength());
    g["bartz_curvature_radius"] = Json(geom.bartzCurvatureRadius());
    double dr = 0.0, ds = 0.0;
    geom.checkContinuity(dr, ds);
    g["max_radius_jump"] = Json(dr);
    g["max_slope_jump"] = Json(ds);
    j["contour"] = g;
    j.writeFile(app::outputPath(cfg, "engine.json"));

    if (!cli.quiet)
      std::cout << "\nwrote " << app::outputPath(cfg, "engine.json") << "\n      "
                << app::outputPath(cfg, "profile.csv") << "\n";
    return 0;
  });
}
