// SPDX-License-Identifier: MIT
/// \file ignis_sweep.cpp
/// \brief Parameter sweeps and constrained design optimisation.

#include <chrono>
#include <iomanip>
#include <iostream>

#include "AppCommon.hpp"
#include "ignis/engine/SteadyEngine.hpp"
#include "ignis/optimize/Optimizer.hpp"
#include "ignis/optimize/Sweep.hpp"

using namespace ignis;

int main(int argc, char** argv) {
  return app::run([&]() -> int {
    CommandLine cli;
    const std::vector<CliOption> options = {
        {"mode", "NAME", "sweep, optimize or both", "both when the configuration has both"},
        {"threads", "N", "worker threads for the sweep", "hardware concurrency"},
        {"max-evaluations", "N", "cap the optimiser's objective evaluations", "from the configuration"},
        {"starts", "N", "Latin-hypercube multi-start count", "from the configuration"},
        {"list-parameters", "", "print every settable parameter and readable metric, then exit", ""},
    };
    if (!parseCommandLine(argc, argv, "ignis_sweep",
                          "Run a full-factorial parameter sweep and/or a constrained design\n"
                          "optimisation over the engine model defined by a configuration file.",
                          options, cli, false))
      return 0;

    if (flagPresent(cli, "list-parameters")) {
      std::cout << "settable parameters:\n";
      for (const auto& p : parameterNames())
        std::cout << "  " << std::left << std::setw(38) << p << parameterUnits(p) << "\n";
      std::cout << "\nreadable metrics:\n";
      for (const auto& m : metricNames())
        std::cout << "  " << std::left << std::setw(38) << m << metricUnits(m) << "\n";
      return 0;
    }
    if (cli.config.empty())
      throw ConfigError("no configuration file given; use --config PATH (or --help)");

    EngineConfig cfg = EngineConfig::load(cli.config);
    app::applyCommandLine(cfg, cli);
    const auto root = cfg.root();
    const bool has_sweep = root.has("sweep");
    const bool has_opt = root.has("optimization");
    std::string mode = flagString(cli, "mode", "");
    if (mode.empty()) mode = (has_sweep && has_opt) ? "both" : (has_sweep ? "sweep" : "optimize");
    if (mode != "sweep" && mode != "optimize" && mode != "both")
      throw ConfigError("--mode must be 'sweep', 'optimize' or 'both'");
    if ((mode == "sweep" || mode == "both") && !has_sweep && mode == "sweep")
      throw ConfigError("configuration '" + cli.config + "' has no 'sweep:' section");
    if ((mode == "optimize" || mode == "both") && !has_opt && mode == "optimize")
      throw ConfigError("configuration '" + cli.config + "' has no 'optimization:' section");

    const SteadyEngine engine(cfg);
    Json j = app::provenance(cfg);
    if (!cli.quiet) std::cout << versionBanner() << "\n\ncase " << cfg.name << "\n\n";

    if (has_sweep && (mode == "sweep" || mode == "both")) {
      auto spec = parseSweep(root);
      if (flagPresent(cli, "threads")) spec.threads = flagInt(cli, "threads", 0);
      const auto t0 = std::chrono::steady_clock::now();
      const auto res = runSweep(engine, spec);
      const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
      res.table.writeCsv(app::outputPath(cfg, "sweep.csv"));
      Json s = res.table.toJson();
      s["points"] = Json(res.points);
      s["failures"] = Json(res.failures);
      s["wall_seconds"] = Json(secs);
      s["failure_messages"] = Json::of(res.failure_messages);
      j["sweep"] = s;
      if (!cli.quiet)
        std::cout << res.summary() << "\n  wall time " << std::fixed << std::setprecision(2)
                  << secs << " s\n  wrote " << app::outputPath(cfg, "sweep.csv") << "\n\n";
    }

    if (has_opt && (mode == "optimize" || mode == "both")) {
      auto spec = parseOptimization(root);
      // Both exist so that a short run can be asked for without editing the
      // configuration: run_all.sh --quick and CI use them.
      if (flagPresent(cli, "max-evaluations")) {
        spec.max_evaluations = flagInt(cli, "max-evaluations", spec.max_evaluations);
        if (spec.max_evaluations < 1)
          throw ConfigError("--max-evaluations must be at least 1");
      }
      if (flagPresent(cli, "starts")) {
        spec.starts = flagInt(cli, "starts", spec.starts);
        if (spec.starts < 1) throw ConfigError("--starts must be at least 1");
      }
      const auto t0 = std::chrono::steady_clock::now();
      const auto res = optimize(engine, spec);
      const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
      res.historyTable(spec).writeCsv(app::outputPath(cfg, "optimization_history.csv"));
      Json o = Json::object();
      o["objective_metric"] = Json(spec.objective);
      o["sense"] = Json(spec.maximize ? std::string("maximize") : std::string("minimize"));
      o["objective"] = Json(res.objective);
      o["feasible"] = Json(res.feasible);
      o["evaluations"] = Json(res.evaluations);
      o["failed_evaluations"] = Json(res.failures);
      o["wall_seconds"] = Json(secs);
      Json vars = Json::object();
      for (std::size_t i = 0; i < res.x.size(); ++i) vars[res.variable_names[i]] = Json(res.x[i]);
      o["variables"] = vars;
      Json cons = Json::array();
      for (std::size_t c = 0; c < spec.constraints.size(); ++c) {
        Json e = Json::object();
        e["metric"] = Json(spec.constraints[c].metric);
        e["op"] = Json(spec.constraints[c].op);
        e["bound"] = Json(spec.constraints[c].bound);
        e["value"] = Json(res.constraint_values[c]);
        e["normalised_violation"] = Json(res.constraint_violations[c]);
        cons.push(e);
      }
      o["constraints"] = cons;
      o["best_design"] = res.best.toJson();
      j["optimization"] = o;

      // The optimum as a ready-to-run override list.
      Table best("optimum");
      std::vector<std::string> names;
      std::vector<double> values;
      for (std::size_t i = 0; i < res.x.size(); ++i) {
        names.push_back(res.variable_names[i]);
        values.push_back(res.x[i]);
      }
      best.addColumn("parameter", names);
      best.addColumn("value", values, "-");
      best.writeCsv(app::outputPath(cfg, "optimum.csv"));

      if (!cli.quiet)
        std::cout << res.summary(spec) << "  wall time            " << std::fixed
                  << std::setprecision(2) << secs << " s\n\n"
                  << res.best.summary() << "\n  wrote "
                  << app::outputPath(cfg, "optimization_history.csv") << "\n";
    }
    j.writeFile(app::outputPath(cfg, "sweep.json"));
    return 0;
  });
}
