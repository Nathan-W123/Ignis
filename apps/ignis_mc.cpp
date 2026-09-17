// SPDX-License-Identifier: MIT
/// \file ignis_mc.cpp
/// \brief Deterministic Monte Carlo uncertainty and sensitivity analysis.

#include <chrono>
#include <iomanip>
#include <iostream>

#include "AppCommon.hpp"
#include "ignis/engine/SteadyEngine.hpp"
#include "ignis/uncertainty/MonteCarlo.hpp"

using namespace ignis;

int main(int argc, char** argv) {
  return app::run([&]() -> int {
    CommandLine cli;
    const std::vector<CliOption> options = {
        {"samples", "N", "override the sample count", "from the configuration"},
        {"seed", "N", "override the campaign seed", "from the configuration"},
        {"threads", "N", "worker threads", "hardware concurrency"},
        {"no-samples", "", "do not write the per-sample table", ""},
    };
    if (!parseCommandLine(argc, argv, "ignis_mc",
                          "Run a deterministic, multithreaded Monte Carlo campaign over the\n"
                          "uncertain inputs declared in a configuration file, and report output\n"
                          "distributions, failure modes and sensitivity rankings.\n\n"
                          "Results depend only on the seed and the sample index, never on the\n"
                          "thread count.",
                          options, cli))
      return 0;

    EngineConfig cfg = EngineConfig::load(cli.config);
    app::applyCommandLine(cfg, cli);
    const auto root = cfg.root();
    if (!root.has("monte_carlo"))
      throw ConfigError("configuration '" + cli.config + "' has no 'monte_carlo:' section");
    auto spec = parseMonteCarlo(root);
    if (flagPresent(cli, "samples")) spec.samples = flagInt(cli, "samples", spec.samples);
    if (flagPresent(cli, "seed"))
      spec.seed = static_cast<std::uint64_t>(flagInt(cli, "seed", 0));
    if (flagPresent(cli, "threads")) spec.threads = flagInt(cli, "threads", 0);
    if (flagPresent(cli, "no-samples")) spec.record_samples = false;

    const SteadyEngine engine(cfg);
    if (!cli.quiet)
      std::cout << versionBanner() << "\n\ncase " << cfg.name << "\nrunning " << spec.samples
                << " samples, seed " << spec.seed << " ...\n" << std::flush;

    const auto t0 = std::chrono::steady_clock::now();
    const auto res = runMonteCarlo(engine, spec);
    const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();

    if (spec.record_samples) res.samples.writeCsv(app::outputPath(cfg, "mc_samples.csv"));
    res.statisticsTable().writeCsv(app::outputPath(cfg, "mc_statistics.csv"));
    res.sensitivityTable().writeCsv(app::outputPath(cfg, "mc_sensitivity.csv"));

    Json j = app::provenance(cfg);
    Json m = Json::object();
    m["samples_requested"] = Json(res.requested);
    m["samples_succeeded"] = Json(res.succeeded);
    m["samples_failed"] = Json(res.failed);
    m["seed"] = Json(static_cast<double>(spec.seed));
    m["wall_seconds"] = Json(secs);
    m["throughput_per_second"] = Json(secs > 0.0 ? res.requested / secs : 0.0);
    Json inputs = Json::array();
    for (const auto& in : spec.inputs) {
      Json e = Json::object();
      e["parameter"] = Json(in.parameter);
      e["distribution"] = Json(toString(in.distribution));
      e["a"] = Json(in.a);
      e["b"] = Json(in.b);
      if (in.distribution == DistributionType::kTriangular) e["c"] = Json(in.c);
      inputs.push(e);
    }
    m["inputs"] = inputs;
    Json fm = Json::object();
    for (const auto& kv : res.failure_modes) fm[kv.first] = Json(kv.second);
    m["failure_modes"] = fm;
    Json fl = Json::object();
    for (const auto& kv : res.flags) fl[kv.first] = Json(kv.second);
    m["constraint_flags"] = fl;
    m["statistics"] = res.statisticsTable().toJson();
    m["sensitivity"] = res.sensitivityTable().toJson();
    j["monte_carlo"] = m;
    j.writeFile(app::outputPath(cfg, "mc.json"));

    if (!cli.quiet)
      std::cout << "\n" << res.summary() << "\nwall time " << std::fixed << std::setprecision(2)
                << secs << " s (" << std::setprecision(1) << res.requested / std::max(secs, 1e-9)
                << " samples/s)\nwrote " << app::outputPath(cfg, "mc_statistics.csv") << "\n";
    if (res.succeeded == 0) {
      std::cerr << "ignis: every Monte Carlo sample failed\n";
      return 4;
    }
    return 0;
  });
}
