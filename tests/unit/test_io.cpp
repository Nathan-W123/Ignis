// SPDX-License-Identifier: MIT
/// \file test_io.cpp
/// \brief Configuration parsing, validation and result export.

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <cstdio>
#include <filesystem>
#include <fstream>

#include "TestHelpers.hpp"
#include "ignis/io/Config.hpp"
#include "ignis/io/Json.hpp"
#include "ignis/io/Table.hpp"
#include "ignis/engine/SteadyEngine.hpp"
#include "ignis/optimize/Parameters.hpp"

using namespace ignis;
using namespace ignis_test;
using Catch::Approx;

namespace {

/// Write a configuration to a scratch file and load it.
class TempConfig {
 public:
  explicit TempConfig(const std::string& yaml) {
    path_ = std::filesystem::temp_directory_path() /
            ("ignis_test_" + std::to_string(reinterpret_cast<std::uintptr_t>(this)) + ".yaml");
    std::ofstream out(path_);
    out << yaml;
  }
  ~TempConfig() { std::error_code ec; std::filesystem::remove(path_, ec); }
  std::string path() const { return path_.string(); }

 private:
  std::filesystem::path path_;
};

const char* kMinimal = R"(
name: test
propellants:
  oxidizer: LOX
  fuel: LCH4
  mixture_ratio: 3.4
chamber:
  pressure: 5.0e6
nozzle:
  throat_radius: 0.05
  contraction_ratio: 3.0
  chamber_length: 0.2
  expansion_ratio: 20.0
)";

}  // namespace

TEST_CASE("every shipped configuration loads and is self-consistent", "[io][config]") {
  int seen = 0;
  for (const auto& entry : std::filesystem::directory_iterator(configDir())) {
    if (entry.path().extension() != ".yaml") continue;
    INFO("configuration " << entry.path().string());
    const auto cfg = EngineConfig::load(entry.path().string());
    REQUIRE_FALSE(cfg.name.empty());
    REQUIRE_FALSE(cfg.description.empty());
    REQUIRE(cfg.chamber_pressure > 0.0);
    REQUIRE(cfg.mixture_ratio > 0.0);
    REQUIRE(cfg.nozzle.expansion_ratio > 1.0);
    REQUIRE(cfg.resolveAmbientPressure() >= 0.0);
    REQUIRE_FALSE(cfg.output_directory.empty());
    REQUIRE_FALSE(cfg.output_prefix.empty());
    // The geometry must build.
    REQUIRE_NOTHROW(NozzleGeometry::build(cfg.nozzle));
    ++seen;
  }
  REQUIRE(seen >= 12);
}

TEST_CASE("a minimal configuration loads with documented defaults", "[io][config]") {
  const TempConfig tc(kMinimal);
  const auto cfg = EngineConfig::load(tc.path());
  REQUIRE(cfg.eta_c_star == Approx(1.0));
  REQUIRE(cfg.composition == CompositionModel::kEquilibrium);
  REQUIRE(cfg.ambient_pressure == Approx(101325.0));
  REQUIRE_FALSE(cfg.cooling_enabled);
  REQUIRE_FALSE(cfg.feed_enabled);
  REQUIRE(cfg.nozzle.converging_half_angle == Approx(30.0));
  REQUIRE(cfg.nozzle.bell_length_fraction == Approx(0.8));
  REQUIRE(cfg.elements == std::vector<std::string>{"C", "H", "O"});
}

TEST_CASE("configuration errors are specific and actionable", "[io][config][errors]") {
  SECTION("a missing file") {
    REQUIRE_THROWS_WITH(EngineConfig::load("/nonexistent/path.yaml"),
                        Catch::Matchers::ContainsSubstring("cannot open"));
  }
  SECTION("malformed YAML") {
    const TempConfig tc("name: [unclosed\n");
    REQUIRE_THROWS_AS(EngineConfig::load(tc.path()), ConfigError);
  }
  SECTION("a missing required key names the key") {
    const TempConfig tc("name: t\npropellants:\n  oxidizer: LOX\n  fuel: LCH4\n");
    REQUIRE_THROWS_WITH(EngineConfig::load(tc.path()),
                        Catch::Matchers::ContainsSubstring("mixture_ratio"));
  }
  SECTION("a misspelt key is rejected, not ignored") {
    std::string yaml = kMinimal;
    yaml += "performance:\n  ambeint_pressure: 1000.0\n";
    const TempConfig tc(yaml);
    REQUIRE_THROWS_WITH(EngineConfig::load(tc.path()),
                        Catch::Matchers::ContainsSubstring("ambeint_pressure"));
    REQUIRE_THROWS_WITH(EngineConfig::load(tc.path()),
                        Catch::Matchers::ContainsSubstring("Recognised keys"));
  }
  SECTION("an out-of-range value names the range") {
    std::string yaml = kMinimal;
    yaml += "chamber:\n";  // duplicate key: YAML keeps the last
    const TempConfig tc(std::string(kMinimal) + "\nperformance:\n  eta_nozzle: 3.0\n");
    REQUIRE_THROWS_WITH(EngineConfig::load(tc.path()),
                        Catch::Matchers::ContainsSubstring("outside the allowed range"));
  }
  SECTION("altitude and ambient pressure cannot both be given") {
    const TempConfig tc(std::string(kMinimal) +
                        "\nperformance:\n  altitude: 0.0\n  ambient_pressure: 101325.0\n");
    REQUIRE_THROWS_WITH(EngineConfig::load(tc.path()),
                        Catch::Matchers::ContainsSubstring("either"));
  }
  SECTION("an unknown propellant lists the available ones") {
    const TempConfig tc(R"(
name: t
propellants:
  oxidizer: LOX
  fuel: kerosene
  mixture_ratio: 2.5
chamber:
  pressure: 5.0e6
nozzle:
  throat_radius: 0.05
  contraction_ratio: 3.0
  chamber_length: 0.2
  expansion_ratio: 20.0
)");
    const auto cfg = EngineConfig::load(tc.path());
    REQUIRE_THROWS_WITH(SteadyEngine(cfg).run(),
                        Catch::Matchers::ContainsSubstring("Available:"));
  }
  SECTION("an unknown composition model") {
    const TempConfig tc(R"(
name: t
propellants: {oxidizer: LOX, fuel: LCH4, mixture_ratio: 3.4}
chamber: {pressure: 5.0e6, composition: magic}
nozzle: {throat_radius: 0.05, contraction_ratio: 3.0, chamber_length: 0.2, expansion_ratio: 20.0}
)");
    REQUIRE_THROWS_WITH(EngineConfig::load(tc.path()),
                        Catch::Matchers::ContainsSubstring("frozen"));
  }
}

TEST_CASE("propellant inlet temperatures are validated", "[io][config][errors]") {
  const auto lib =
      PropellantLibrary::loadYaml(sourceDir() + "/data/propellants/ignis_propellants.yaml");
  const auto& db = fullDatabase();
  const auto& ch4 = lib.at("LCH4");
  REQUIRE_NOTHROW(ch4.molarEnthalpy(111.66, db));
  REQUIRE_NOTHROW(ch4.molarEnthalpy(120.0, db));
  // Far from the tabulated reference the linear sensible correction is refused.
  REQUIRE_THROWS_WITH(ch4.molarEnthalpy(300.0, db),
                      Catch::Matchers::ContainsSubstring("linear sensible correction"));
  // Gaseous entries use the polynomial directly and accept any valid T.
  REQUIRE_NOTHROW(lib.at("GCH4").molarEnthalpy(400.0, db));
}

TEST_CASE("named parameters and metrics round trip", "[io][parameters]") {
  auto cfg = EngineConfig::load(configDir() + "/methane_nominal.yaml");
  REQUIRE(parameterNames().size() > 20);
  REQUIRE(metricNames().size() > 30);
  for (const auto& name : parameterNames()) {
    INFO("parameter " << name);
    REQUIRE_FALSE(parameterUnits(name).empty());
    const double original = readParameter(cfg, name);
    const double probe = (original != 0.0) ? 1.05 * original : 1.0;
    applyParameter(cfg, name, probe);
    REQUIRE(readParameter(cfg, name) == Approx(probe).epsilon(1e-9));
    applyParameter(cfg, name, original == 0.0 ? 0.0 : original);
  }
  for (const auto& name : metricNames()) REQUIRE_FALSE(metricUnits(name).empty());

  SECTION("unknown names list the valid options") {
    REQUIRE_THROWS_WITH(applyParameter(cfg, "chamber.presure", 1.0),
                        Catch::Matchers::ContainsSubstring("chamber.pressure"));
    REQUIRE_THROWS_WITH(metricUnits("performance.thurst"),
                        Catch::Matchers::ContainsSubstring("performance.thrust"));
  }
  SECTION("a non-finite value is rejected") {
    REQUIRE_THROWS_AS(applyParameter(cfg, "chamber.pressure",
                                     std::numeric_limits<double>::quiet_NaN()),
                      ConfigError);
  }
}

TEST_CASE("JSON output is well formed", "[io][json]") {
  Json j = Json::object();
  j["name"] = Json("case");
  j["value"] = Json(1.5);
  j["flag"] = Json(true);
  j["list"] = Json::of(std::vector<double>{1.0, 2.0, 3.0});
  j["nested"]["deep"] = Json(42);
  const std::string text = j.dump();
  REQUIRE(text.find("\"name\": \"case\"") != std::string::npos);
  REQUIRE(text.find("\"flag\": true") != std::string::npos);
  REQUIRE(text.find("[1, 2, 3]") != std::string::npos);
  REQUIRE(text.find("\"deep\": 42") != std::string::npos);
  // Braces and brackets must balance.
  int depth = 0;
  for (char c : text) {
    if (c == '{' || c == '[') ++depth;
    if (c == '}' || c == ']') --depth;
    REQUIRE(depth >= 0);
  }
  REQUIRE(depth == 0);

  SECTION("non-finite numbers become null rather than invalid JSON") {
    Json k = Json::object();
    k["nan"] = Json(std::nan(""));
    REQUIRE(k.dump().find("null") != std::string::npos);
  }
  SECTION("strings are escaped") {
    Json k = Json::object();
    k["s"] = Json(std::string("a\"b\\c\nd"));
    const auto d = k.dump();
    REQUIRE(d.find("\\\"") != std::string::npos);
    REQUIRE(d.find("\\\\") != std::string::npos);
    REQUIRE(d.find("\\n") != std::string::npos);
  }
}

TEST_CASE("tables export consistent CSV and JSON", "[io][table]") {
  Table t("demo");
  t.addColumn("x", std::vector<double>{1.0, 2.0, 3.0}, "m");
  t.addColumn("y", std::vector<double>{10.0, 20.0, 30.0}, "K");
  t.addColumn("label", std::vector<std::string>{"a", "b", "c"});
  REQUIRE(t.rows() == 3);
  REQUIRE(t.columns() == 3);
  REQUIRE(t.numeric(0)[2] == Approx(3.0));
  REQUIRE_THROWS_AS(t.numeric(2), ConfigError);
  REQUIRE_THROWS_AS(t.addColumn("bad", std::vector<double>{1.0}), ConfigError);

  const auto path = (std::filesystem::temp_directory_path() / "ignis_table_test.csv").string();
  t.writeCsv(path);
  std::ifstream in(path);
  std::string line;
  std::getline(in, line);
  REQUIRE(line.rfind("# units:", 0) == 0);
  std::getline(in, line);
  REQUIRE(line == "x,y,label");
  std::getline(in, line);
  REQUIRE(line == "1,10,a");
  std::error_code ec;
  std::filesystem::remove(path, ec);

  const auto j = t.toJson().dump();
  REQUIRE(j.find("\"rows\": 3") != std::string::npos);
  REQUIRE(j.find("\"x\": [1, 2, 3]") != std::string::npos);
}
