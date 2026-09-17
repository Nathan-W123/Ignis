// SPDX-License-Identifier: MIT
/// \file test_engine.cpp
/// \brief End-to-end engine, sweep, optimisation and feed-system behaviour.

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <cmath>

#include "TestHelpers.hpp"
#include "ignis/engine/SteadyEngine.hpp"
#include "ignis/optimize/Optimizer.hpp"
#include "ignis/optimize/Parameters.hpp"
#include "ignis/optimize/Sweep.hpp"
#include "ignis/nozzle/Atmosphere.hpp"

using namespace ignis;
using namespace ignis_test;
using Catch::Approx;

TEST_CASE("the nominal methane engine runs end to end and is self-consistent",
          "[integration][engine]") {
  const auto cfg = EngineConfig::load(configDir() + "/methane_nominal.yaml");
  const SteadyEngine engine(cfg);
  const auto res = engine.run();

  SECTION("identities between the reported performance quantities") {
    REQUIRE(res.performance.thrust_ideal ==
            Approx(res.performance.thrust_momentum + res.performance.thrust_pressure)
                .epsilon(1e-12));
    REQUIRE(res.performance.isp_ideal ==
            Approx(res.performance.thrust_ideal / (res.mdot * constants::g0)).epsilon(1e-12));
    REQUIRE(res.performance.c_effective ==
            Approx(res.performance.thrust / res.mdot).epsilon(1e-12));
    REQUIRE(res.performance.cf ==
            Approx(res.performance.thrust / (cfg.chamber_pressure * res.throat_area))
                .epsilon(1e-12));
    REQUIRE(res.performance.thrust_momentum ==
            Approx(res.mdot * res.performance.u_exit).epsilon(1e-12));
    REQUIRE(res.mdot == Approx(res.mdot_oxidizer + res.mdot_fuel).epsilon(1e-12));
    REQUIRE(res.mdot_oxidizer / res.mdot_fuel == Approx(cfg.mixture_ratio).epsilon(1e-12));
    REQUIRE(res.chamber.c_star == Approx(cfg.eta_c_star * res.chamber.c_star_ideal).epsilon(1e-12));
    REQUIRE(res.mdot == Approx(cfg.chamber_pressure * res.throat_area / res.chamber.c_star)
                            .epsilon(1e-12));
    REQUIRE(res.l_star == Approx(res.chamber_volume / res.throat_area).epsilon(1e-12));
  }
  SECTION("solver residuals are small") {
    REQUIRE(res.chamber.diagnostics.element_residual_rel < 1.0e-10);
    REQUIRE(res.chamber.diagnostics.gibbs_residual < 1.0e-6);
    REQUIRE(res.chamber.diagnostics.state_residual < 1.0e-10);
    REQUIRE(res.performance.mass_flow_residual < 1.0e-9);
    REQUIRE(res.performance.energy_residual < 1.0e-12);
    REQUIRE(res.cooling.energy_balance_residual < 1.0e-8);
    REQUIRE(res.cooling.max_flux_residual < 1.0e-6);
  }
  SECTION("the design point is in the intended regime") {
    REQUIRE(res.has_cooling);
    REQUIRE(res.has_feed);
    REQUIRE_FALSE(res.cooling.wall_limit_exceeded);
    REQUIRE_FALSE(res.cooling.boiling_detected);
    REQUIRE(res.feed.stiffness_ok);
    REQUIRE(res.performance.regime == ExpansionRegime::kOverExpanded);
  }
  SECTION("exports are complete") {
    const auto profile = res.profileTable();
    REQUIRE(profile.rows() == res.profile.x.size());
    REQUIRE(profile.rows() > 300);
    const auto thermal = res.coolingTable();
    REQUIRE(thermal.rows() == res.cooling.stations.size());
    const auto json = res.toJson().dump();
    REQUIRE(json.find("\"c_star_ideal\"") != std::string::npos);
    REQUIRE(json.find("\"max_wall_temperature\"") != std::string::npos);
    REQUIRE(json.find("\"oxidizer\"") != std::string::npos);
  }
  SECTION("the summary mentions every major section") {
    const auto text = res.summary();
    for (const char* needle : {"chamber", "nozzle performance", "regenerative cooling",
                               "feed system", "c* ideal", "Isp (vacuum)"})
      REQUIRE(text.find(needle) != std::string::npos);
  }
}

TEST_CASE("the nominal hydrogen engine runs end to end", "[integration][engine][hydrogen]") {
  const auto cfg = EngineConfig::load(configDir() + "/hydrogen_nominal.yaml");
  const SteadyEngine engine(cfg);
  const auto res = engine.run();
  REQUIRE(res.chamber.state.T > 3000.0);
  REQUIRE(res.chamber.state.T < 3800.0);
  // A hydrogen-rich exhaust has a much lower molar mass than a methane one.
  REQUIRE(res.chamber.state.M * 1e3 < 16.0);
  REQUIRE(res.chamber.c_star_ideal > 2200.0);
  REQUIRE(res.performance.isp_vacuum > 420.0);
  REQUIRE(res.performance.regime == ExpansionRegime::kUnderExpanded);  // vacuum
  REQUIRE(res.has_cooling);
  REQUIRE(res.cooling.energy_balance_residual < 1.0e-8);
  REQUIRE(res.cooling.coolant_temperature_rise > 0.0);
}

TEST_CASE("altitude changes performance in the expected direction",
          "[integration][engine]") {
  auto cfg = EngineConfig::load(configDir() + "/methane_vacuum.yaml");
  cfg.cooling_enabled = false;
  cfg.feed_enabled = false;
  cfg.sample_profile = false;
  const SteadyEngine engine(cfg);
  double previous_thrust = 0.0, previous_isp = 0.0;
  for (double z : {0.0, 5000.0, 15000.0, 30000.0, 60000.0}) {
    const auto res = engine.runAt(Atmosphere::at(z).pressure);
    INFO("altitude " << z << " m");
    REQUIRE(res.performance.thrust > previous_thrust);
    REQUIRE(res.performance.isp > previous_isp);
    previous_thrust = res.performance.thrust;
    previous_isp = res.performance.isp;
    // Mass flow and c* do not depend on the ambient pressure.
    REQUIRE(res.mdot == Approx(engine.runAt(101325.0).mdot).epsilon(1e-12));
  }
  const auto vac = engine.runAt(0.0);
  REQUIRE(vac.performance.isp == Approx(vac.performance.isp_vacuum).epsilon(1e-9));
}

TEST_CASE("a sweep reproduces the single-point analysis exactly",
          "[integration][sweep]") {
  auto cfg = EngineConfig::load(configDir() + "/methane_nominal.yaml");
  cfg.sample_profile = false;
  const SteadyEngine engine(cfg);

  SweepSpec spec;
  spec.axes.push_back(SweepAxis::linear("propellants.mixture_ratio", 3.0, 3.8, 5));
  spec.metrics = {"chamber.temperature", "performance.isp", "cooling.max_wall_temperature"};
  spec.threads = 3;
  const auto res = runSweep(engine, spec);
  REQUIRE(res.points == 5);
  REQUIRE(res.failures == 0);
  REQUIRE(res.table.rows() == 5);

  // Row 2 is O/F = 3.4, which is the configuration's own design point.
  const auto direct = engine.run();
  const auto& mr = res.table.numeric(0);
  const auto& T = res.table.numeric(1);
  REQUIRE(mr[2] == Approx(3.4).epsilon(1e-12));
  REQUIRE(T[2] == Approx(direct.chamber.state.T).epsilon(1e-12));

  SECTION("a sweep is independent of the thread count") {
    auto s2 = spec;
    s2.threads = 1;
    const auto single = runSweep(engine, s2);
    for (std::size_t c = 0; c + 1 < res.table.columns(); ++c)
      for (std::size_t r = 0; r < res.table.rows(); ++r)
        REQUIRE(res.table.numeric(c)[r] == single.table.numeric(c)[r]);
  }
  SECTION("a two-dimensional sweep visits every combination") {
    auto nozzle_only = cfg;
    nozzle_only.cooling_enabled = false;   // a pure performance grid
    nozzle_only.feed_enabled = false;
    const SteadyEngine nozzle_engine(nozzle_only);
    SweepSpec two;
    two.axes.push_back(SweepAxis::linear("propellants.mixture_ratio", 3.0, 3.8, 3));
    two.axes.push_back(SweepAxis::logarithmic("chamber.pressure", 3.0e6, 9.0e6, 4));
    two.metrics = {"performance.isp", "chamber.temperature"};
    const auto r2 = runSweep(nozzle_engine, two);
    REQUIRE(r2.points == 12);
    REQUIRE(r2.table.rows() == 12);
    REQUIRE(r2.failures == 0);
    // The last axis varies fastest, so the first four rows share a mixture ratio.
    const auto& mr_col = r2.table.numeric(0);
    const auto& pc_col = r2.table.numeric(1);
    for (int i = 0; i < 4; ++i) REQUIRE(mr_col[static_cast<std::size_t>(i)] == Approx(3.0));
    REQUIRE(pc_col[0] == Approx(3.0e6));
    REQUIRE(pc_col[3] == Approx(9.0e6));
    REQUIRE(mr_col[4] == Approx(3.4));
  }
  SECTION("failed points are recorded rather than dropped") {
    SweepSpec bad;
    bad.axes.push_back(SweepAxis::linear("nozzle.contraction_ratio", 0.5, 3.0, 6));
    bad.metrics = {"performance.isp"};
    const auto rb = runSweep(engine, bad);
    REQUIRE(rb.points == 6);
    REQUIRE(rb.failures > 0);
    REQUIRE(rb.table.rows() == 6);
    REQUIRE_FALSE(rb.failure_messages.empty());
  }
  SECTION("an unknown parameter or metric is rejected before any work") {
    SweepSpec bad = spec;
    bad.axes[0].parameter = "chamber.presure";
    REQUIRE_THROWS_AS(runSweep(engine, bad), ConfigError);
    bad = spec;
    bad.metrics = {"performance.thurst"};
    REQUIRE_THROWS_AS(runSweep(engine, bad), ConfigError);
  }
}

TEST_CASE("constrained optimisation respects its constraints", "[integration][optimize]") {
  auto cfg = EngineConfig::load(configDir() + "/optimization.yaml");
  cfg.sample_profile = false;
  const SteadyEngine engine(cfg);
  auto spec = parseOptimization(cfg.root());
  // Keep the test affordable; the shipped campaign is larger.
  spec.starts = 2;
  spec.max_evaluations = 400;
  spec.outer_iterations = 4;

  const auto res = optimize(engine, spec);
  INFO(res.summary(spec));
  REQUIRE(res.feasible);
  REQUIRE(res.evaluations > 20);
  REQUIRE(res.evaluations <= spec.max_evaluations + 50);
  REQUIRE(res.x.size() == spec.variables.size());

  SECTION("the reported optimum satisfies every constraint") {
    for (std::size_t c = 0; c < spec.constraints.size(); ++c) {
      INFO("constraint on " << spec.constraints[c].metric);
      REQUIRE(res.constraint_violations[c] <= 1.0e-9);
      if (spec.constraints[c].op == "<=")
        REQUIRE(res.constraint_values[c] <= spec.constraints[c].bound * (1.0 + 1e-9));
      else
        REQUIRE(res.constraint_values[c] >= spec.constraints[c].bound * (1.0 - 1e-9));
    }
  }
  SECTION("every variable is inside its bounds") {
    for (std::size_t i = 0; i < res.x.size(); ++i) {
      REQUIRE(res.x[i] >= spec.variables[i].min - 1e-9);
      REQUIRE(res.x[i] <= spec.variables[i].max + 1e-9);
    }
  }
  SECTION("the optimum is at least as good as the starting design") {
    auto start = engine.config();
    for (const auto& v : spec.variables)
      if (v.start > 0.0) applyParameter(start, v.parameter, v.start);
    const auto base = engine.runWith(start);
    REQUIRE(res.objective >= readMetric(base, spec.objective) - 1e-9);
  }
  SECTION("rerunning the reported design reproduces the reported objective") {
    auto best = engine.config();
    for (std::size_t i = 0; i < res.x.size(); ++i)
      applyParameter(best, res.variable_names[i], res.x[i]);
    const auto check = engine.runWith(best);
    REQUIRE(readMetric(check, spec.objective) == Approx(res.objective).epsilon(1e-9));
  }
  SECTION("the convergence history is exported") {
    const auto t = res.historyTable(spec);
    REQUIRE(t.rows() == res.history.size());
    REQUIRE(t.rows() == static_cast<std::size_t>(res.evaluations));
  }
  SECTION("an impossible constraint set is reported, not silently relaxed") {
    auto impossible = spec;
    impossible.constraints.push_back({"cooling.max_wall_temperature", "<=", 100.0, 0.0});
    impossible.starts = 1;
    impossible.max_evaluations = 120;
    REQUIRE_THROWS_AS(optimize(engine, impossible), InfeasibleError);
  }
}

TEST_CASE("the feed system sizes and evaluates consistently", "[integration][feed]") {
  const auto cfg = EngineConfig::load(configDir() + "/methane_nominal.yaml");
  const SteadyEngine engine(cfg);
  const auto res = engine.run();
  REQUIRE(res.has_feed);

  SECTION("sizing hits the requested stiffness") {
    REQUIRE(res.feed.oxidizer.stiffness == Approx(cfg.injector_stiffness).epsilon(1e-12));
    REQUIRE(res.feed.fuel.stiffness == Approx(cfg.injector_stiffness).epsilon(1e-12));
    REQUIRE(res.feed.oxidizer.dp_injector ==
            Approx(cfg.injector_stiffness * cfg.chamber_pressure).epsilon(1e-12));
  }
  SECTION("the tank pressure is the sum of every loss") {
    const auto& l = res.feed.oxidizer;
    REQUIRE(l.required_tank_pressure ==
            Approx(cfg.chamber_pressure + l.dp_injector + l.dp_lines + l.dp_dynamic)
                .epsilon(1e-12));
  }
  SECTION("capability mode recovers the sized flows") {
    FeedSystemSpec fs = cfg.feed;
    fs.oxidizer.density = engine.propellants().at(cfg.oxidizer).density;
    fs.fuel.density = engine.propellants().at(cfg.fuel).density;
    fs.oxidizer.mass_flow = res.mdot_oxidizer;
    fs.fuel.mass_flow = res.mdot_fuel;
    fs.oxidizer.injector_area = res.feed.oxidizer.injector_area;
    fs.fuel.injector_area = res.feed.fuel.injector_area;
    fs.oxidizer.tank_pressure = res.feed.oxidizer.required_tank_pressure;
    fs.fuel.tank_pressure = res.feed.fuel.required_tank_pressure;
    const auto cap = evaluateFeedSystem(fs, cfg.chamber_pressure);
    REQUIRE(cap.oxidizer.mass_flow == Approx(res.mdot_oxidizer).epsilon(1e-6));
    REQUIRE(cap.fuel.mass_flow == Approx(res.mdot_fuel).epsilon(1e-6));
    REQUIRE(cap.mixture_ratio == Approx(cfg.mixture_ratio).epsilon(1e-5));
  }
  SECTION("a tank that cannot supply the chamber is reported") {
    FeedSystemSpec fs = cfg.feed;
    fs.oxidizer.density = 1141.0;
    fs.fuel.density = 422.0;
    fs.oxidizer.injector_area = 1e-3;
    fs.fuel.injector_area = 1e-3;
    fs.oxidizer.tank_pressure = 1.0e6;   // below the chamber pressure
    fs.fuel.tank_pressure = 8.0e6;
    REQUIRE_THROWS_WITH(evaluateFeedSystem(fs, cfg.chamber_pressure),
                        Catch::Matchers::ContainsSubstring("cannot supply"));
  }
  SECTION("pressurant mass is the isothermal ideal-gas bound") {
    const double m = pressurantMass(6.0e6, 2.0, 0.028, 300.0);
    REQUIRE(m == Approx(6.0e6 * 2.0 * 0.028 / (constants::R_universal * 300.0)).epsilon(1e-12));
    REQUIRE_THROWS_AS(pressurantMass(-1.0, 2.0, 0.028, 300.0), ConfigError);
  }
}
