// SPDX-License-Identifier: MIT
/// \file test_uncertainty.cpp
/// \brief Monte Carlo determinism, statistics and sensitivity.

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <cmath>

#include "TestHelpers.hpp"
#include "ignis/engine/SteadyEngine.hpp"
#include "ignis/uncertainty/MonteCarlo.hpp"

using namespace ignis;
using namespace ignis_test;
using Catch::Approx;

namespace {

/// A deliberately cheap campaign: cooling off, few samples, so the determinism
/// checks stay fast.
MonteCarloSpec cheapSpec(int samples = 64) {
  MonteCarloSpec spec;
  spec.samples = samples;
  spec.seed = 987654321ull;
  spec.record_samples = true;
  spec.inputs = {
      {"chamber.pressure", DistributionType::kNormal, 5.5e6, 1.1e5, 0.0,
       1.0e6, 2.0e7, false},
      {"propellants.mixture_ratio", DistributionType::kNormal, 3.4, 0.08, 0.0, 2.5, 4.5, false},
      {"nozzle.throat_radius", DistributionType::kNormal, 0.07, 2.0e-4, 0.0, 0.05, 0.09, false},
      {"chamber.eta_c_star", DistributionType::kTriangular, 0.93, 0.96, 0.98, 0.9, 1.0, false},
  };
  spec.outputs = {"performance.thrust", "performance.isp", "chamber.temperature",
                  "performance.mdot"};
  return spec;
}

SteadyEngine cheapEngine() {
  auto cfg = EngineConfig::load(configDir() + "/methane_nominal.yaml");
  cfg.cooling_enabled = false;
  cfg.feed_enabled = false;
  cfg.sample_profile = false;
  return SteadyEngine(cfg);
}

}  // namespace

TEST_CASE("Monte Carlo results are independent of the thread count",
          "[uncertainty][determinism]") {
  const auto engine = cheapEngine();
  auto spec = cheapSpec();

  spec.threads = 1;
  const auto single = runMonteCarlo(engine, spec);
  spec.threads = 4;
  const auto four = runMonteCarlo(engine, spec);
  spec.threads = 7;
  const auto seven = runMonteCarlo(engine, spec);

  REQUIRE(single.succeeded == four.succeeded);
  REQUIRE(single.succeeded == seven.succeeded);
  REQUIRE(single.samples.rows() == four.samples.rows());

  // Every drawn input and every computed output must be bit-for-bit identical.
  for (std::size_t c = 0; c < single.samples.columns(); ++c) {
    if (single.samples.headers()[c] == "status") continue;
    const auto& a = single.samples.numeric(c);
    const auto& b = four.samples.numeric(c);
    const auto& d = seven.samples.numeric(c);
    INFO("column " << single.samples.headers()[c]);
    for (std::size_t r = 0; r < a.size(); ++r) {
      if (std::isnan(a[r])) {
        REQUIRE(std::isnan(b[r]));
        REQUIRE(std::isnan(d[r]));
      } else {
        REQUIRE(a[r] == b[r]);
        REQUIRE(a[r] == d[r]);
      }
    }
  }
  for (const auto& name : spec.outputs) {
    REQUIRE(single.statistics.at(name).mean == four.statistics.at(name).mean);
    REQUIRE(single.statistics.at(name).stddev == seven.statistics.at(name).stddev);
  }
}

TEST_CASE("a different seed gives a different sample, the same seed the same one",
          "[uncertainty][determinism]") {
  const auto engine = cheapEngine();
  auto spec = cheapSpec(32);
  const auto a = runMonteCarlo(engine, spec);
  const auto again = runMonteCarlo(engine, spec);
  spec.seed += 1;
  const auto b = runMonteCarlo(engine, spec);

  REQUIRE(a.statistics.at("performance.thrust").mean ==
          again.statistics.at("performance.thrust").mean);
  REQUIRE(a.statistics.at("performance.thrust").mean !=
          b.statistics.at("performance.thrust").mean);
}

TEST_CASE("drawn inputs follow the requested distributions", "[uncertainty][slow]") {
  const auto engine = cheapEngine();
  auto spec = cheapSpec(4000);
  spec.threads = 0;
  // One campaign feeds every check below.  As SECTIONs they re-ran the
  // whole 4000-sample campaign once each.
  const auto res = runMonteCarlo(engine, spec);
  REQUIRE(res.succeeded > 3900);

  {  // the normal inputs reproduce their mean and standard deviation
    INFO("the normal inputs reproduce their mean and standard deviation");
    const auto& pc = res.samples.numeric(0);
    const auto st = computeStatistics(pc);
    REQUIRE(st.mean == Approx(5.5e6).epsilon(0.01));
    REQUIRE(st.stddev == Approx(1.1e5).epsilon(0.06));
    REQUIRE(std::abs(st.skewness) < 0.15);
  }
  {  // bounds are respected
    INFO("bounds are respected");
    const auto& mr = res.samples.numeric(1);
    for (double v : mr) {
      REQUIRE(v >= 2.5);
      REQUIRE(v <= 4.5);
    }
  }
  {  // the triangular input stays inside its support and is skewed
    INFO("the triangular input stays inside its support and is skewed");
    const auto& eta = res.samples.numeric(3);
    const auto st = computeStatistics(eta);
    REQUIRE(st.minimum >= 0.93);
    REQUIRE(st.maximum <= 0.98);
    // Mean of a triangular distribution is (low + mode + high)/3.
    REQUIRE(st.mean == Approx((0.93 + 0.96 + 0.98) / 3.0).epsilon(0.01));
  }
  {  // percentiles are ordered
    INFO("percentiles are ordered");
    for (const auto& name : spec.outputs) {
      const auto& s = res.statistics.at(name);
      INFO("output " << name);
      REQUIRE(s.minimum <= s.p5);
      REQUIRE(s.p5 <= s.p50);
      REQUIRE(s.p50 <= s.p95);
      REQUIRE(s.p95 <= s.p99);
      REQUIRE(s.p99 <= s.maximum);
      REQUIRE(s.stddev > 0.0);
    }
  }
}

TEST_CASE("sensitivity rankings recover known dependencies", "[uncertainty][sensitivity][slow]") {
  const auto engine = cheapEngine();
  auto spec = cheapSpec(1500);
  // One campaign feeds every ranking check below; SECTIONs would re-run it.
  const auto res = runMonteCarlo(engine, spec);
  REQUIRE(res.sensitivity.size() == spec.outputs.size());

  auto find = [&](const std::string& output) {
    for (const auto& s : res.sensitivity)
      if (s.output == output) return s;
    throw std::runtime_error("missing sensitivity for " + output);
  };
  auto indexOf = [&](const std::string& in) {
    for (std::size_t i = 0; i < res.input_names.size(); ++i)
      if (res.input_names[i] == in) return i;
    throw std::runtime_error("missing input " + in);
  };

  {  // mass flow is driven by chamber pressure and throat area
    INFO("mass flow is driven by chamber pressure and throat area");
    const auto s = find("performance.mdot");
    REQUIRE(s.r_squared > 0.95);
    const double src_pc = std::abs(s.src[indexOf("chamber.pressure")]);
    const double src_rt = std::abs(s.src[indexOf("nozzle.throat_radius")]);
    const double src_mr = std::abs(s.src[indexOf("propellants.mixture_ratio")]);
    REQUIRE(src_pc > src_mr);
    REQUIRE(src_rt > src_mr);
    // mdot scales as p_c and as r_t^2, so the elasticities are known exactly.
    REQUIRE(s.elasticity[indexOf("chamber.pressure")] == Approx(1.0).epsilon(0.05));
    REQUIRE(s.elasticity[indexOf("nozzle.throat_radius")] == Approx(2.0).epsilon(0.05));
  }
  {  // flame temperature is driven by mixture ratio
    INFO("flame temperature is driven by mixture ratio");
    const auto s = find("chamber.temperature");
    const double src_mr = std::abs(s.src[indexOf("propellants.mixture_ratio")]);
    for (std::size_t i = 0; i < res.input_names.size(); ++i)
      if (res.input_names[i] != "propellants.mixture_ratio")
        REQUIRE(src_mr > std::abs(s.src[i]));
  }
  {  // specific impulse is insensitive to throat size at fixed area ratio
    INFO("specific impulse is insensitive to throat size at fixed area ratio");
    const auto s = find("performance.isp");
    REQUIRE(std::abs(s.elasticity[indexOf("nozzle.throat_radius")]) < 0.05);
  }
  {  // Spearman and SRC agree in sign for the dominant inputs
    INFO("Spearman and SRC agree in sign for the dominant inputs");
    for (const auto& s : res.sensitivity)
      for (std::size_t i = 0; i < res.input_names.size(); ++i)
        if (std::abs(s.src[i]) > 0.3) REQUIRE(s.src[i] * s.spearman[i] > 0.0);
  }
}

TEST_CASE("failed samples are counted and classified, not silently dropped",
          "[uncertainty][errors]") {
  const auto engine = cheapEngine();
  MonteCarloSpec spec;
  spec.samples = 60;
  spec.seed = 13ull;
  // A contraction ratio below one is geometrically impossible, so part of this
  // uniform range must fail and be reported rather than quietly disappearing.
  spec.inputs = {{"nozzle.contraction_ratio", DistributionType::kUniform, 0.8, 3.0, 0.0,
                  0.8, 3.0, false}};
  spec.outputs = {"performance.isp"};
  const auto res = runMonteCarlo(engine, spec);
  INFO(res.summary());
  REQUIRE(res.requested == 60);
  REQUIRE(res.succeeded + res.failed == 60);
  REQUIRE(res.failed > 0);
  REQUIRE_FALSE(res.failure_modes.empty());
  int counted = 0;
  for (const auto& kv : res.failure_modes) counted += kv.second;
  REQUIRE(counted == res.failed);
  // The failed rows must still be present, flagged, with NaN outputs.
  int nan_rows = 0;
  const auto& isp = res.samples.numeric(1);
  for (double v : isp) if (std::isnan(v)) ++nan_rows;
  REQUIRE(nan_rows == res.failed);
}

TEST_CASE("Monte Carlo input validation", "[uncertainty][errors]") {
  const auto engine = cheapEngine();
  MonteCarloSpec spec = cheapSpec(8);
  SECTION("no inputs") {
    spec.inputs.clear();
    REQUIRE_THROWS_AS(runMonteCarlo(engine, spec), ConfigError);
  }
  SECTION("no outputs") {
    spec.outputs.clear();
    REQUIRE_THROWS_AS(runMonteCarlo(engine, spec), ConfigError);
  }
  SECTION("too few samples") {
    spec.samples = 1;
    REQUIRE_THROWS_AS(runMonteCarlo(engine, spec), ConfigError);
  }
  SECTION("an unknown parameter") {
    spec.inputs[0].parameter = "chamber.presure";
    REQUIRE_THROWS_AS(runMonteCarlo(engine, spec), ConfigError);
  }
  SECTION("bounds that exclude the distribution") {
    spec.inputs[0].lower_bound = 1.0e9;
    spec.inputs[0].upper_bound = 2.0e9;
    REQUIRE_THROWS_WITH(runMonteCarlo(engine, spec),
                        Catch::Matchers::ContainsSubstring("inconsistent"));
  }
}

TEST_CASE("statistics helper handles edge cases", "[uncertainty]") {
  REQUIRE(computeStatistics({}).count == 0);
  const auto one = computeStatistics({5.0});
  REQUIRE(one.count == 1);
  REQUIRE(one.mean == Approx(5.0));
  REQUIRE(one.stddev == Approx(0.0));
  const auto with_nan = computeStatistics({1.0, std::nan(""), 3.0});
  REQUIRE(with_nan.count == 2);
  REQUIRE(with_nan.mean == Approx(2.0));
  const auto uniform = computeStatistics({1.0, 2.0, 3.0, 4.0, 5.0});
  REQUIRE(uniform.p50 == Approx(3.0));
  REQUIRE(uniform.minimum == Approx(1.0));
  REQUIRE(uniform.maximum == Approx(5.0));
}
