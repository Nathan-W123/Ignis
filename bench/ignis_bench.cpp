// SPDX-License-Identifier: MIT
/// \file ignis_bench.cpp
/// \brief Measured runtime, scaling and residual benchmarks.
///
/// Every number this program prints is measured on the machine it runs on.  It
/// records the build and host configuration alongside the timings so a result
/// can be compared honestly with another machine's.

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <sstream>
#include <thread>
#include <vector>

#include "ignis/engine/SteadyEngine.hpp"
#include "ignis/core/Version.hpp"
#include "ignis/io/Cli.hpp"
#include "ignis/io/Json.hpp"
#include "ignis/io/Table.hpp"
#include "ignis/optimize/Parameters.hpp"
#include "ignis/transient/TransientChamber.hpp"
#include "ignis/uncertainty/MonteCarlo.hpp"

using namespace ignis;
using Clock = std::chrono::steady_clock;

namespace {

struct Timing {
  std::string name;
  std::string unit_of_work;
  int repeats = 0;
  double mean_ms = 0.0;
  double min_ms = 0.0;
  double max_ms = 0.0;
  double extra = 0.0;          ///< optional secondary figure (iterations, points, ...)
  std::string extra_label;
};

std::vector<Timing> g_timings;

/// Time `fn` `repeats` times and record the distribution.
template <typename Fn>
Timing measure(const std::string& name, const std::string& work, int repeats, Fn&& fn) {
  std::vector<double> samples;
  samples.reserve(static_cast<std::size_t>(repeats));
  for (int i = 0; i < repeats; ++i) {
    const auto t0 = Clock::now();
    fn(i);
    const auto t1 = Clock::now();
    samples.push_back(1e3 * std::chrono::duration<double>(t1 - t0).count());
  }
  Timing t;
  t.name = name;
  t.unit_of_work = work;
  t.repeats = repeats;
  t.mean_ms = std::accumulate(samples.begin(), samples.end(), 0.0) / samples.size();
  t.min_ms = *std::min_element(samples.begin(), samples.end());
  t.max_ms = *std::max_element(samples.begin(), samples.end());
  return t;
}

void record(Timing t) {
  std::cout << "  " << std::left << std::setw(46) << t.name << std::right << std::fixed
            << std::setprecision(4) << std::setw(12) << t.mean_ms << " ms  (min "
            << std::setw(10) << t.min_ms << ", max " << std::setw(10) << t.max_ms << ", n "
            << t.repeats << ")";
  if (!t.extra_label.empty())
    std::cout << "  [" << t.extra_label << " " << std::setprecision(3) << t.extra << "]";
  std::cout << "\n";
  g_timings.push_back(std::move(t));
}

std::string readFirstLine(const char* path, const char* key) {
  std::ifstream in(path);
  std::string line;
  while (std::getline(in, line))
    if (line.rfind(key, 0) == 0) {
      const auto colon = line.find(':');
      if (colon != std::string::npos) {
        std::string v = line.substr(colon + 1);
        while (!v.empty() && (v.front() == ' ' || v.front() == '\t')) v.erase(v.begin());
        return v;
      }
    }
  return "unknown";
}

}  // namespace

int main(int argc, char** argv) {
  try {
    CommandLine cli;
    const std::vector<CliOption> options = {
        {"repeats", "N", "repeat count for the fast benchmarks", "50"},
        {"mc-samples", "N", "Monte Carlo sample count", "400"},
        {"skip-transient", "", "skip the transient benchmark (it builds a table)", ""},
    };
    if (!parseCommandLine(argc, argv, "ignis_bench",
                          "Measure Ignis runtime, scaling and residuals, and record the build\n"
                          "and host configuration alongside them.",
                          options, cli, false))
      return 0;

    const std::string out_dir = cli.output_dir.empty() ? std::string("results/benchmarks")
                                                  : cli.output_dir;
    const int repeats = flagInt(cli, "repeats", 50);
    const int mc_samples = flagInt(cli, "mc-samples", 400);
    const auto t_start = Clock::now();

    // ---------------------------------------------------------------- config
    const unsigned hw = std::max(1u, std::thread::hardware_concurrency());
    const std::string cpu = readFirstLine("/proc/cpuinfo", "model name");
    std::cout << versionBanner() << "\n"
              << "  hardware threads : " << hw << "\n"
              << "  cpu              : " << cpu << "\n"
              << "  eigen            : " << EIGEN_WORLD_VERSION << "." << EIGEN_MAJOR_VERSION
              << "." << EIGEN_MINOR_VERSION << "\n\n";

    auto cfg = EngineConfig::load("configs/methane_nominal.yaml");
    const auto full_db = SpeciesDatabase::loadYaml(findDefaultSpeciesDatabase());

    // ------------------------------------------------------------ equilibrium
    std::cout << "equilibrium\n";
    {
      const auto db = full_db.restrictToElements({"C", "H", "O"});
      const EquilibriumSolver solver(db);
      const auto lib = PropellantLibrary::loadYaml(findDefaultPropellantLibrary());
      const PropellantMixture mix(lib.at("LOX"), lib.at("LCH4"), 3.4, 90.18, 111.66);
      const Eigen::VectorXd b = mix.elementMoles(db);
      const double h0 = mix.enthalpy(db);

      double iters = 0.0, elem = 0.0, gibbs = 0.0, energy = 0.0;
      auto t = measure("adiabatic (HP) solve, 26 species", "solve", repeats * 10, [&](int i) {
        const auto r = solver.hp(b, h0, 5.5e6 + i);
        iters += r.diagnostics.iterations;
        elem = std::max(elem, r.diagnostics.element_residual_rel);
        gibbs = std::max(gibbs, r.diagnostics.gibbs_residual);
        energy = std::max(energy, r.diagnostics.state_residual);
      });
      t.extra = iters / (repeats * 10);
      t.extra_label = "mean Newton iterations";
      record(t);
      std::cout << "    worst residuals over " << repeats * 10 << " solves: element "
                << std::scientific << std::setprecision(2) << elem << ", Gibbs " << gibbs
                << ", enthalpy " << energy << std::fixed << "\n";

      auto t2 = measure("isentropic (SP) solve, 26 species", "solve", repeats * 10, [&](int i) {
        const auto ref = solver.hp(b, h0, 5.5e6);
        (void)solver.sp(b, ref.state.s, 1.0e6 + 1.0e3 * i, ref.state.T);
      });
      record(t2);

      // Scaling with species count.
      std::cout << "  scaling with species count\n";
      std::vector<double> ns, ms;
      for (const auto& names : std::vector<std::vector<std::string>>{
               {"H2", "H", "O", "O2", "OH", "H2O", "CO", "CO2"},
               {"H2", "H", "O", "O2", "OH", "H2O", "HO2", "H2O2", "CO", "CO2", "CH4", "CH3"},
               {"H2", "H", "O", "O2", "OH", "H2O", "HO2", "H2O2", "CO", "CO2", "C", "CH", "CH2",
                "CH3", "CH4", "HCO", "CH2O", "CH2OH"},
               db.names()}) {
        const auto sub = full_db.subset(names);
        const EquilibriumSolver s2(sub);
        const auto b2 = mix.elementMoles(sub);
        const double h2 = mix.enthalpy(sub);
        auto ts = measure("  " + std::to_string(sub.size()) + " species", "solve", repeats * 5,
                          [&](int i) { (void)s2.hp(b2, h2, 5.5e6 + i); });
        ts.extra = static_cast<double>(sub.size());
        ts.extra_label = "species";
        record(ts);
        ns.push_back(static_cast<double>(sub.size()));
        ms.push_back(ts.mean_ms);
      }
      Table t_species("equilibrium_species_scaling");
      t_species.addColumn("species", ns, "-");
      t_species.addColumn("mean_time", ms, "ms");
      t_species.writeCsv(out_dir + "/species_scaling.csv");
    }

    // ---------------------------------------------------------------- nozzle
    std::cout << "\nnozzle\n";
    {
      auto c = cfg;
      c.cooling_enabled = false;
      c.feed_enabled = false;
      c.sample_profile = false;
      const SteadyEngine engine(c);
      record(measure("chamber + throat + exit (no profile, no cooling)", "analysis", repeats,
                     [&](int i) {
                       auto k = c;
                       k.chamber_pressure = 5.5e6 + i;
                       (void)engine.runWith(k);
                     }));

      std::cout << "  scaling with station count\n";
      std::vector<double> stations, times;
      for (int n : {50, 100, 200, 400, 800, 1600}) {
        auto k = c;
        k.nozzle.num_stations = n;
        k.sample_profile = true;
        const SteadyEngine e2(k);
        auto ts = measure("  " + std::to_string(n) + " stations", "profile",
                          std::max(3, repeats / 8), [&](int i) {
                            auto q = k;
                            q.chamber_pressure = 5.5e6 + i;
                            (void)e2.runWith(q);
                          });
        ts.extra = n;
        ts.extra_label = "stations";
        record(ts);
        stations.push_back(n);
        times.push_back(ts.mean_ms);
      }
      Table t_st("nozzle_station_scaling");
      t_st.addColumn("stations", stations, "-");
      t_st.addColumn("mean_time", times, "ms");
      t_st.writeCsv(out_dir + "/station_scaling.csv");
    }

    // --------------------------------------------------------- whole analysis
    std::cout << "\ncomplete steady analysis\n";
    {
      const SteadyEngine engine(cfg);
      double res_energy = 0.0, res_flux = 0.0, res_mass = 0.0;
      auto t = measure("chamber + nozzle + profile + cooling + feed", "analysis",
                       std::max(5, repeats / 4), [&](int i) {
                         auto k = cfg;
                         k.chamber_pressure = 5.5e6 + i;
                         const auto r = engine.runWith(k);
                         res_energy = std::max(res_energy, r.cooling.energy_balance_residual);
                         res_flux = std::max(res_flux, r.cooling.max_flux_residual);
                         res_mass = std::max(res_mass, r.performance.mass_flow_residual);
                       });
      record(t);
      std::cout << "    worst residuals: cooling energy balance " << std::scientific
                << std::setprecision(2) << res_energy << ", local flux " << res_flux
                << ", nozzle mass flow " << res_mass << std::fixed << "\n";

      auto c2 = cfg;
      c2.cooling_enabled = false;
      c2.feed_enabled = false;
      const SteadyEngine e2(c2);
      record(measure("chamber + nozzle + profile (no cooling)", "analysis",
                     std::max(5, repeats / 4), [&](int i) {
                       auto k = c2;
                       k.chamber_pressure = 5.5e6 + i;
                       (void)e2.runWith(k);
                     }));
    }

    // ------------------------------------------------------------- transient
    if (!flagPresent(cli, "skip-transient")) {
      std::cout << "\ntransient\n";
      auto tcfg = EngineConfig::load("configs/startup_transient.yaml");
      const SteadyEngine engine(tcfg);
      auto tc = parseTransient(tcfg.root(), tcfg);
      const auto geom = NozzleGeometry::build(tcfg.nozzle);
      tc.spec.chamber_volume = geom.chamberVolume();
      tc.spec.throat_area = geom.throatArea();
      tc.spec.exit_area = geom.exitArea();
      const auto mix = engine.mixture(tcfg);
      const auto& ox = engine.propellants().at(tcfg.oxidizer);
      const auto& fu = engine.propellants().at(tcfg.fuel);
      tc.spec.oxidizer_inlet_enthalpy =
          ox.molarEnthalpy(mix.oxidizerTemperature(), engine.database()) / ox.molar_mass;
      tc.spec.fuel_inlet_enthalpy =
          fu.molarEnthalpy(mix.fuelTemperature(), engine.database()) / fu.molar_mass;

      EquilibriumTable table = [&] {
        const auto t0 = Clock::now();
        auto tab = EquilibriumTable::build(engine.solver(), mix, tc.grid);
        const auto t1 = Clock::now();
        Timing tt;
        tt.name = "equilibrium table build";
        tt.unit_of_work = "table";
        tt.repeats = 1;
        tt.mean_ms = tt.min_ms = tt.max_ms = 1e3 * std::chrono::duration<double>(t1 - t0).count();
        tt.extra = tc.grid.mr_points * tc.grid.t_points * tc.grid.p_points;
        tt.extra_label = "equilibrium solves";
        record(tt);
        return tab;
      }();

      const auto err = table.measureError(engine.solver(), mix, 100, 20260917u);
      std::cout << "    interpolation error: max |dT|/T " << std::scientific
                << std::setprecision(2) << err.max_rel_temperature << ", rms "
                << err.rms_rel_temperature << " over " << err.samples << " points\n"
                << std::fixed;

      double steps = 0.0, mass_err = 0.0, energy_err = 0.0;
      auto t45 = measure("adaptive RK4(5), 0.8 s of physical time", "run", 3, [&](int) {
        auto s = tc.spec;
        s.integrator = "rk45";
        const auto r = simulateTransient(table, s);
        steps = r.steps;
        mass_err = r.mass_conservation_error;
        energy_err = r.energy_conservation_error;
      });
      t45.extra = steps;
      t45.extra_label = "accepted steps";
      record(t45);
      std::cout << "    conservation residuals: mass " << std::scientific << std::setprecision(2)
                << mass_err << ", energy " << energy_err << std::fixed << "\n";

      auto t4 = measure("fixed-step RK4 at 2 us, 0.8 s of physical time", "run", 1, [&](int) {
        auto s = tc.spec;
        s.integrator = "rk4";
        s.dt = 2.0e-6;
        const auto r = simulateTransient(table, s);
        steps = r.steps;
      });
      t4.extra = steps;
      t4.extra_label = "steps";
      record(t4);
    }

    // ----------------------------------------------------------- Monte Carlo
    std::cout << "\nMonte Carlo\n";
    {
      auto mcfg = EngineConfig::load("configs/monte_carlo.yaml");
      const SteadyEngine engine(mcfg);
      auto spec = parseMonteCarlo(mcfg.root());
      spec.samples = mc_samples;
      spec.record_samples = false;
      std::vector<double> threads_col, rate_col;
      for (int threads : {1, 2, 4, static_cast<int>(hw)}) {
        if (threads > static_cast<int>(hw)) continue;
        if (!threads_col.empty() && threads_col.back() == threads) continue;
        auto s = spec;
        s.threads = threads;
        int ok = 0;
        auto t = measure("  " + std::to_string(threads) + " thread(s)", "campaign", 1,
                         [&](int) { ok = runMonteCarlo(engine, s).succeeded; });
        t.extra = 1e3 * spec.samples / t.mean_ms;
        t.extra_label = "samples/s";
        record(t);
        std::cout << "      " << ok << " of " << spec.samples << " samples succeeded\n";
        threads_col.push_back(threads);
        rate_col.push_back(t.extra);
      }
      Table t_mc("monte_carlo_scaling");
      t_mc.addColumn("threads", threads_col, "-");
      t_mc.addColumn("throughput", rate_col, "samples/s");
      t_mc.writeCsv(out_dir + "/monte_carlo_scaling.csv");
      if (rate_col.size() > 1)
        std::cout << "    parallel speed-up at " << threads_col.back() << " threads: "
                  << std::setprecision(2) << rate_col.back() / rate_col.front() << "x\n";
    }

    // --------------------------------------------------------------- export
    const double total_s = std::chrono::duration<double>(Clock::now() - t_start).count();
    std::cout << "\nbenchmark suite wall time: " << std::setprecision(2) << total_s << " s\n";

    Table summary("benchmarks");
    std::vector<std::string> names, units;
    std::vector<double> mean, lo, hi, reps, extra;
    std::vector<std::string> extra_labels;
    for (const auto& t : g_timings) {
      names.push_back(t.name);
      units.push_back(t.unit_of_work);
      mean.push_back(t.mean_ms);
      lo.push_back(t.min_ms);
      hi.push_back(t.max_ms);
      reps.push_back(t.repeats);
      extra.push_back(t.extra);
      extra_labels.push_back(t.extra_label);
    }
    summary.addColumn("benchmark", names);
    summary.addColumn("unit_of_work", units);
    summary.addColumn("mean_time", mean, "ms");
    summary.addColumn("min_time", lo, "ms");
    summary.addColumn("max_time", hi, "ms");
    summary.addColumn("repeats", reps, "-");
    summary.addColumn("secondary_value", extra, "-");
    summary.addColumn("secondary_label", extra_labels);
    summary.writeCsv(out_dir + "/benchmarks.csv");

    Json j = Json::object();
    j["ignis_version"] = Json(std::string(kVersion));
    j["git"] = Json(std::string(kGitDescribe));
    j["build_type"] = Json(std::string(kBuildType));
    j["compiler"] = Json(std::string(kCompilerId) + " " + kCompilerVersion);
    j["cpu"] = Json(cpu);
    j["hardware_threads"] = Json(static_cast<double>(hw));
    {
      std::ostringstream os;
      os << EIGEN_WORLD_VERSION << "." << EIGEN_MAJOR_VERSION << "." << EIGEN_MINOR_VERSION;
      j["eigen"] = Json(os.str());
    }
    j["suite_wall_seconds"] = Json(total_s);
    j["benchmarks"] = summary.toJson();
    j.writeFile(out_dir + "/benchmarks.json");
    std::cout << "wrote " << out_dir << "/benchmarks.csv and benchmarks.json\n";
    return 0;
  } catch (const std::exception& e) {
    std::cerr << "ignis_bench: " << e.what() << "\n";
    return 2;
  }
}
