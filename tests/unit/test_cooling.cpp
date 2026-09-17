// SPDX-License-Identifier: MIT
/// \file test_cooling.cpp
/// \brief Regenerative-cooling verification.

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <cmath>

#include "TestHelpers.hpp"
#include "ignis/combustion/Chamber.hpp"
#include "ignis/thermal/RegenerativeCooling.hpp"

using namespace ignis;
using namespace ignis_test;
using Catch::Approx;

namespace {

struct Fixture {
  SpeciesDatabase db = rocketDatabase();
  EquilibriumSolver solver{db};
  TransportModel transport{db};
  PropellantLibrary lib =
      PropellantLibrary::loadYaml(sourceDir() + "/data/propellants/ignis_propellants.yaml");
  PropellantMixture mix{lib.at("LOX"), lib.at("LCH4"), 3.4, 90.18, 111.66};
  CombustionChamber chamber{solver, CompositionModel::kEquilibrium};
  ChamberResult ch = chamber.solve(mix, 5.5e6, 0.96);
  NozzleGeometry geom = [] {
    NozzleGeometrySpec s;
    s.throat_radius = 0.070;
    s.contraction_ratio = 2.8;
    s.chamber_length = 0.22;
    s.expansion_ratio = 20.0;
    s.num_stations = 400;
    return NozzleGeometry::build(s);
  }();
  NozzleFlow flow = chamber.makeFlow(ch);

  CoolingSpec spec() const {
    CoolingSpec s;
    s.coolant = "methane";
    s.num_channels = 300;
    s.width_fraction = 0.50;
    s.channel_height = 5.0e-3;
    s.wall_thickness = 0.6e-3;
    s.coolant_mass_flow = ch.massFlowFor(geom.throatArea()) / 4.4;
    s.inlet_temperature = 111.66;
    s.inlet_pressure = 15.0e6;
    s.counterflow = true;
    s.x_end = -1.0;
    s.num_segments = 200;
    s.wall = MaterialLibrary::loadDefault().at("CuCrZr");
    return s;
  }
};

}  // namespace

TEST_CASE("regenerative cooling closes its energy balance", "[cooling]") {
  Fixture f;
  const auto res = solveRegenerativeCooling(f.flow, f.geom, f.ch, f.spec(), f.transport);

  SECTION("heat leaving the gas equals heat entering the coolant") {
    INFO("energy balance residual " << res.energy_balance_residual);
    REQUIRE(res.energy_balance_residual < 1.0e-8);
  }
  SECTION("every station balances gas-side and coolant-side flux") {
    INFO("largest local flux residual " << res.max_flux_residual);
    REQUIRE(res.max_flux_residual < 1.0e-6);
    for (const auto& s : res.stations) REQUIRE(s.flux_residual < 1.0e-5);
  }
  SECTION("the summed segment heat equals the reported total") {
    double sum = 0.0;
    for (const auto& s : res.stations) sum += s.segment_heat;
    REQUIRE(sum == Approx(res.total_heat_load).epsilon(1e-12));
  }
  SECTION("the temperature ordering is physical at every station") {
    for (const auto& s : res.stations) {
      INFO("station at x = " << s.x);
      REQUIRE(s.t_adiabatic_wall > s.t_wall_hot);
      REQUIRE(s.t_wall_hot > s.t_wall_cold);
      REQUIRE(s.t_wall_cold > s.coolant_T);
      REQUIRE(s.q_total > 0.0);
      REQUIRE(s.h_gas > 0.0);
      REQUIRE(s.h_coolant > 0.0);
      REQUIRE(s.fin_efficiency > 0.0);
      REQUIRE(s.fin_efficiency <= 1.0);
      REQUIRE(s.reynolds > 4000.0);
      REQUIRE(s.land_width > 0.0);
    }
  }
  SECTION("the peak flux sits near the throat") {
    REQUIRE(std::abs(res.max_heat_flux_x - f.geom.throatPosition()) < 0.05);
  }
}

TEST_CASE("coolant temperature rises monotonically along its own path", "[cooling]") {
  Fixture f;
  auto spec = f.spec();

  SECTION("counterflow: the coolant enters at the exit and heats towards the injector") {
    spec.counterflow = true;
    const auto res = solveRegenerativeCooling(f.flow, f.geom, f.ch, spec, f.transport);
    // Stations are stored by increasing x, so a counterflow coolant is hottest
    // at the smallest x.
    for (std::size_t i = 1; i < res.stations.size(); ++i)
      REQUIRE(res.stations[i - 1].coolant_T > res.stations[i].coolant_T);
    REQUIRE(res.coolant_temperature_rise > 0.0);
    REQUIRE(res.stations.front().coolant_T ==
            Approx(res.coolant_outlet_temperature).epsilon(0.02));
  }
  SECTION("co-flow: the coolant enters at the injector") {
    spec.counterflow = false;
    const auto res = solveRegenerativeCooling(f.flow, f.geom, f.ch, spec, f.transport);
    for (std::size_t i = 1; i < res.stations.size(); ++i)
      REQUIRE(res.stations[i].coolant_T > res.stations[i - 1].coolant_T);
    REQUIRE(res.coolant_temperature_rise > 0.0);
  }
  SECTION("pressure falls monotonically along the flow direction") {
    const auto res = solveRegenerativeCooling(f.flow, f.geom, f.ch, spec, f.transport);
    REQUIRE(res.coolant_pressure_drop > 0.0);
    REQUIRE(res.coolant_outlet_pressure < spec.inlet_pressure);
    for (std::size_t i = 1; i < res.stations.size(); ++i)
      REQUIRE(res.stations[i - 1].coolant_p < res.stations[i].coolant_p);
  }
}

TEST_CASE("cooling results converge with station count", "[cooling][convergence]") {
  Fixture f;
  double previous_q = 0.0, previous_T = 0.0, previous_dp = 0.0;
  double last_dq = 1e30, last_dT = 1e30;
  for (int n : {50, 100, 200, 400, 800}) {
    auto spec = f.spec();
    spec.num_segments = n;
    const auto res = solveRegenerativeCooling(f.flow, f.geom, f.ch, spec, f.transport);
    if (previous_q > 0.0) {
      const double dq = relativeError(res.total_heat_load, previous_q);
      const double dT = relativeError(res.max_wall_temperature, previous_T);
      const double ddp = relativeError(res.coolant_pressure_drop, previous_dp);
      INFO("n = " << n << ": dQ " << dq << ", dT_wall " << dT << ", dp " << ddp);
      // Successive refinements must keep shrinking the change.
      REQUIRE(dq < last_dq);
      REQUIRE(dq < 0.02);
      REQUIRE(dT < 0.02);
      last_dq = dq;
      last_dT = dT;
    }
    previous_q = res.total_heat_load;
    previous_T = res.max_wall_temperature;
    previous_dp = res.coolant_pressure_drop;
  }
  INFO("final refinement changes: heat load " << last_dq << ", wall temperature " << last_dT);
  REQUIRE(last_dq < 2.0e-3);
  REQUIRE(last_dT < 2.0e-3);
}

TEST_CASE("cooling responds to design changes in the expected direction", "[cooling]") {
  Fixture f;
  const auto base = solveRegenerativeCooling(f.flow, f.geom, f.ch, f.spec(), f.transport);

  SECTION("a thicker wall runs hotter") {
    auto s = f.spec();
    s.wall_thickness *= 2.0;
    const auto r = solveRegenerativeCooling(f.flow, f.geom, f.ch, s, f.transport);
    REQUIRE(r.max_wall_temperature > base.max_wall_temperature);
  }
  SECTION("a less conductive wall runs hotter") {
    auto s = f.spec();
    s.wall = MaterialLibrary::loadDefault().at("Inconel718");
    s.wall.max_temperature = 5000.0;   // so the run is not flagged, only measured
    const auto r = solveRegenerativeCooling(f.flow, f.geom, f.ch, s, f.transport);
    REQUIRE(r.max_wall_temperature > base.max_wall_temperature + 200.0);
  }
  SECTION("larger channels drop less pressure but cool less well") {
    auto s = f.spec();
    s.channel_height *= 1.6;
    const auto r = solveRegenerativeCooling(f.flow, f.geom, f.ch, s, f.transport);
    REQUIRE(r.coolant_pressure_drop < base.coolant_pressure_drop);
    REQUIRE(r.max_wall_temperature > base.max_wall_temperature);
  }
  SECTION("the Bartz multiplier scales the heat load nearly linearly") {
    auto s = f.spec();
    s.bartz_multiplier = 1.25;
    const auto r = solveRegenerativeCooling(f.flow, f.geom, f.ch, s, f.transport);
    const double ratio = r.total_heat_load / base.total_heat_load;
    REQUIRE(ratio > 1.10);
    REQUIRE(ratio < 1.25);   // sub-linear because the wall heats up
  }
  SECTION("adding radiation increases the flux and is reported separately") {
    auto s = f.spec();
    s.gas_emissivity = 0.15;
    const auto r = solveRegenerativeCooling(f.flow, f.geom, f.ch, s, f.transport);
    REQUIRE(r.total_heat_load > base.total_heat_load);
    double rad = 0.0;
    for (const auto& st : r.stations) rad = std::max(rad, st.q_radiative);
    REQUIRE(rad > 0.0);
    for (const auto& st : base.stations) REQUIRE(st.q_radiative == 0.0);
  }
  SECTION("Gnielinski and Dittus-Boelter bracket the coolant-side uncertainty") {
    auto s = f.spec();
    s.nusselt_correlation = "gnielinski";
    const auto r = solveRegenerativeCooling(f.flow, f.geom, f.ch, s, f.transport);
    const double spread = relativeError(r.max_wall_temperature, base.max_wall_temperature);
    INFO("peak wall temperature: Dittus-Boelter " << base.max_wall_temperature
         << " K, Gnielinski " << r.max_wall_temperature << " K (" << 100.0 * spread << " %)");
    REQUIRE(r.max_wall_temperature != base.max_wall_temperature);
    // The two published correlations differ by around ten percent in peak wall
    // temperature for this design.  That spread is the real coolant-side
    // correlation uncertainty and is why nusselt_multiplier is a dispersed
    // Monte Carlo input; it is recorded in docs/validation.md.
    REQUIRE(spread < 0.20);
  }
}

TEST_CASE("cooling reports its failure modes", "[cooling][errors]") {
  Fixture f;

  SECTION("a wall above the material limit is flagged, not hidden") {
    auto s = f.spec();
    s.wall.max_temperature = 400.0;
    const auto r = solveRegenerativeCooling(f.flow, f.geom, f.ch, s, f.transport);
    REQUIRE(r.wall_limit_exceeded);
    REQUIRE_FALSE(r.warnings.empty());
  }
  SECTION("channels that leave no land are rejected") {
    auto s = f.spec();
    s.num_channels = 4000;
    REQUIRE_THROWS_WITH(solveRegenerativeCooling(f.flow, f.geom, f.ch, s, f.transport),
                        Catch::Matchers::ContainsSubstring("land"));
  }
  SECTION("a flow the channels cannot pass is rejected") {
    auto s = f.spec();
    s.channel_height = 2.0e-4;
    s.width_fraction = 0.12;
    REQUIRE_THROWS_AS(solveRegenerativeCooling(f.flow, f.geom, f.ch, s, f.transport),
                      IgnisError);
  }
  SECTION("a boiling inlet is rejected") {
    auto s = f.spec();
    s.inlet_pressure = 2.0e6;      // below the 4.6 MPa critical pressure
    s.inlet_temperature = 180.0;   // above T_sat at that pressure
    REQUIRE_THROWS_WITH(solveRegenerativeCooling(f.flow, f.geom, f.ch, s, f.transport),
                        Catch::Matchers::ContainsSubstring("vapour"));
  }
  SECTION("invalid configuration is rejected before any work is done") {
    auto s = f.spec();
    s.num_segments = 2;
    REQUIRE_THROWS_AS(solveRegenerativeCooling(f.flow, f.geom, f.ch, s, f.transport),
                      ConfigError);
    s = f.spec();
    s.coolant_mass_flow = 0.0;
    REQUIRE_THROWS_AS(solveRegenerativeCooling(f.flow, f.geom, f.ch, s, f.transport),
                      ConfigError);
  }
}

TEST_CASE("the hot-gas survey isolates the gas side", "[cooling]") {
  Fixture f;
  auto s = f.spec();
  const auto r = surveyHotGasSide(f.flow, f.geom, f.ch, s, f.transport, 800.0);
  for (const auto& st : r.stations) {
    REQUIRE(st.t_wall_hot == Approx(800.0));
    REQUIRE(st.q_total == Approx(st.h_gas * (st.t_adiabatic_wall - 800.0)).epsilon(1e-10));
    REQUIRE(st.t_wall_cold < st.t_wall_hot);
  }
  REQUIRE(r.total_heat_load > 0.0);
  REQUIRE_THROWS_AS(surveyHotGasSide(f.flow, f.geom, f.ch, s, f.transport, 0.0), ConfigError);
}
