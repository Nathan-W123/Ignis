// SPDX-License-Identifier: MIT
/// \file test_film_cooling.cpp
/// \brief The Hatch & Papell film correlation, and films and tapered channels
///        inside the coupled cooling solve.

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cmath>

#include "TestHelpers.hpp"
#include "ignis/combustion/Chamber.hpp"
#include "ignis/core/Constants.hpp"
#include "ignis/core/Exceptions.hpp"
#include "ignis/thermal/FilmCooling.hpp"
#include "ignis/thermal/RegenerativeCooling.hpp"

using namespace ignis;
using namespace ignis_test;
using Catch::Approx;

TEST_CASE("the Hatch and Papell velocity function is the report's", "[thermal][film]") {
  SECTION("both branches meet at a velocity ratio of one") {
    CHECK(hatchPapellVelocityFunction(1.0) == Approx(1.0));
    CHECK(hatchPapellVelocityFunction(1.0 + 1e-9) == Approx(1.0).epsilon(1e-8));
    CHECK(hatchPapellVelocityFunction(1.0 - 1e-9) == Approx(1.0).epsilon(1e-8));
  }
  SECTION("eq. (10) above it, with the angle in radians") {
    CHECK(hatchPapellVelocityFunction(2.0) == Approx(1.0 + 0.4 * constants::pi / 4.0));
  }
  SECTION("eq. (11) below it") {
    // (V_c/V_g)^(1.5 (V_c/V_g - 1)) at V_g/V_c = 0.5: 2^1.5.
    CHECK(hatchPapellVelocityFunction(0.5) == Approx(std::pow(2.0, 1.5)));
  }
  SECTION("mismatched velocities always cost effectiveness") {
    for (double r : {0.5, 0.8, 1.5, 4.0, 20.0}) CHECK(hatchPapellVelocityFunction(r) > 1.0);
  }
}

TEST_CASE("the effectiveness is eq. (12)", "[thermal][film]") {
  SECTION("the film is perfect until the heat-capacity group reaches 0.04") {
    CHECK(hatchPapellEffectiveness(0.0, 3.0) == 1.0);
    CHECK(hatchPapellEffectiveness(0.04, 3.0) == 1.0);
  }
  SECTION("beyond it ln(eta) falls linearly with slope -factor") {
    const double a = std::log(hatchPapellEffectiveness(0.5, 3.0));
    const double b = std::log(hatchPapellEffectiveness(1.0, 3.0));
    CHECK(a == Approx(-(0.5 - 0.04) * 3.0));
    CHECK((b - a) / 0.5 == Approx(-3.0));
  }
  SECTION("bad input is refused") {
    CHECK_THROWS_AS(hatchPapellEffectiveness(-0.1, 3.0), ConfigError);
    CHECK_THROWS_AS(hatchPapellEffectiveness(0.5, 0.0), ConfigError);
  }
}

TEST_CASE("the slot factor holds the velocity ratio inside the fitted range", "[thermal][film]") {
  FilmSlot s;
  s.mass_flow = 0.01;
  s.slot_height = 3.0e-3;
  s.perimeter = 0.2;
  s.temperature = 400.0;
  s.density = 0.1;
  s.cp = 5193.0;
  s.conductivity = 0.2;
  s.gas_velocity = 100.0;
  // V_c = 0.01 / (0.1 * 0.2 * 0.003) = 166.7 m/s, so V_g/V_c = 0.6: inside.
  CHECK(s.coolantVelocity() == Approx(0.01 / (0.1 * 0.2 * 3.0e-3)));
  CHECK(s.diffusivity() == Approx(0.2 / (0.1 * 5193.0)));
  bool clamped = true;
  const double inside = hatchPapellFactor(s, &clamped);
  CHECK_FALSE(clamped);
  CHECK(inside == Approx(std::pow(3.0e-3 * 100.0 / s.diffusivity(), 0.125) *
                         hatchPapellVelocityFunction(0.6)));

  SECTION("a coolant jet much faster than the gas is clamped and flagged") {
    s.mass_flow = 0.1;   // V_g/V_c = 0.06
    const double f = hatchPapellFactor(s, &clamped);
    CHECK(clamped);
    CHECK(f == Approx(std::pow(3.0e-3 * 100.0 / s.diffusivity(), 0.125) *
                      hatchPapellVelocityFunction(kHatchPapellMinVelocityRatio)));
  }
  SECTION("and a very slow one likewise") {
    s.mass_flow = 1.0e-4;  // V_g/V_c = 60
    hatchPapellFactor(s, &clamped);
    CHECK(clamped);
  }
}

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
  double mdot_fuel = ch.massFlowFor(geom.throatArea()) / 4.4;

  CoolingSpec spec() const {
    CoolingSpec s;
    s.coolant = "methane";
    s.num_channels = 300;
    s.width_fraction = 0.50;
    s.channel_height = 5.0e-3;
    s.wall_thickness = 0.6e-3;
    s.coolant_mass_flow = mdot_fuel;
    s.inlet_temperature = 111.66;
    s.inlet_pressure = 15.0e6;
    s.counterflow = true;
    s.num_segments = 160;
    s.wall = MaterialLibrary::loadDefault().at("CuCrZr");
    return s;
  }
};

}  // namespace

TEST_CASE("a wall film inside the coupled solve", "[cooling][film]") {
  Fixture f;
  const auto bare = solveRegenerativeCooling(f.flow, f.geom, f.ch, f.spec(), f.transport);
  auto spec = f.spec();
  spec.film_mass_flow = 0.03 * f.mdot_fuel;
  const auto filmed = solveRegenerativeCooling(f.flow, f.geom, f.ch, spec, f.transport);

  SECTION("the energy balance still closes") {
    CHECK(filmed.energy_balance_residual < 1.0e-8);
    CHECK(filmed.max_flux_residual < 1.0e-6);
  }
  SECTION("the film is perfect at the slot and fades downstream") {
    REQUIRE(filmed.has_film);
    CHECK(filmed.stations.front().film_effectiveness == Approx(1.0));
    CHECK(filmed.stations.back().film_effectiveness < 0.01);
    for (std::size_t i = 1; i < filmed.stations.size(); ++i)
      CHECK(filmed.stations[i].film_effectiveness <= filmed.stations[i - 1].film_effectiveness);
    CHECK(filmed.film_length_half > 0.0);
    CHECK(filmed.film_length_fifth > filmed.film_length_half);
  }
  SECTION("it lowers the driving temperature, never raises it") {
    for (const auto& s : filmed.stations) {
      CHECK(s.t_drive <= s.t_adiabatic_wall + 1e-9);
      CHECK(s.t_drive >= filmed.film_temperature - 1e-9);
    }
  }
  SECTION("it cools the wall it covers and the total heat load falls") {
    CHECK(filmed.stations.front().t_wall_hot < bare.stations.front().t_wall_hot - 100.0);
    CHECK(filmed.total_heat_load < bare.total_heat_load);
  }
  SECTION("by default the film is the fuel leaving the jacket") {
    // Counterflow: the jacket's outlet is the injector end the film leaves
    // from.  The passes iterate the two to 0.5 K.
    CHECK(filmed.film_temperature == Approx(filmed.coolant_outlet_temperature).margin(0.5));
    CHECK(filmed.wall_iterations > 1);
  }
  SECTION("a survey has no outlet to draw a film from, so it must be told") {
    CHECK_THROWS_AS(surveyHotGasSide(f.flow, f.geom, f.ch, spec, f.transport, 600.0),
                    ConfigError);
    auto told = spec;
    told.film_temperature = 300.0;
    const auto s = surveyHotGasSide(f.flow, f.geom, f.ch, told, f.transport, 600.0);
    CHECK(s.film_temperature == 300.0);
  }
  SECTION("a film colder than the coolant reverses the local heat flow without failing") {
    // At the injector end of a counterflow jacket the coolant is near its
    // outlet temperature; a film far colder than that draws heat out of it.
    auto cold = spec;
    cold.film_temperature = 120.0;
    const auto r = solveRegenerativeCooling(f.flow, f.geom, f.ch, cold, f.transport);
    CHECK(r.stations.front().q_total < 0.0);
    CHECK(r.energy_balance_residual < 1.0e-8);
  }
  SECTION("a film at the coolant's own temperature leaves nothing to bracket, and still solves") {
    // A film drawn from the jacket outlet sits, at its slot, at the very
    // temperature of the coolant flowing past behind the wall.  The local
    // driving difference is then ~zero and the wall balance degenerates; the
    // root bracket used to pinch shut there.
    auto same = spec;
    same.film_temperature = filmed.coolant_outlet_temperature;
    const auto r = solveRegenerativeCooling(f.flow, f.geom, f.ch, same, f.transport);
    CHECK(std::abs(r.stations.front().q_total) < 0.05 * std::abs(bare.stations.front().q_total));
    CHECK(r.max_flux_residual < 1.0e-6);
    CHECK(r.energy_balance_residual < 1.0e-8);
  }
  SECTION("a film injected downstream does nothing upstream of its slot") {
    auto late = spec;
    late.film_x = 0.10;
    late.film_temperature = 400.0;
    const auto r = solveRegenerativeCooling(f.flow, f.geom, f.ch, late, f.transport);
    for (const auto& s : r.stations) {
      if (s.x < 0.10) CHECK(s.film_effectiveness == 0.0);
    }
  }
}

TEST_CASE("tapered channels", "[cooling][taper]") {
  Fixture f;
  const auto base = solveRegenerativeCooling(f.flow, f.geom, f.ch, f.spec(), f.transport);

  SECTION("a flat profile is the constant channel, exactly") {
    auto flat = f.spec();
    flat.channel_height_profile = {{0.0, 5.0e-3}, {0.5, 5.0e-3}};
    const auto r = solveRegenerativeCooling(f.flow, f.geom, f.ch, flat, f.transport);
    CHECK(r.total_heat_load == Approx(base.total_heat_load).epsilon(1e-12));
    CHECK(r.coolant_pressure_drop == Approx(base.coolant_pressure_drop).epsilon(1e-12));
  }

  SECTION("narrowing the channel through the throat speeds the coolant and cools the wall there") {
    auto taper = f.spec();
    const double xt = f.geom.throatPosition();
    taper.channel_height_profile = {{xt - 0.10, 5.0e-3}, {xt, 2.5e-3}, {xt + 0.10, 5.0e-3}};
    const auto r = solveRegenerativeCooling(f.flow, f.geom, f.ch, taper, f.transport);
    std::size_t k = 0;
    for (std::size_t i = 0; i < r.stations.size(); ++i)
      if (std::abs(r.stations[i].x - xt) < std::abs(r.stations[k].x - xt)) k = i;
    CHECK(r.stations[k].channel_height < 2.7e-3);
    CHECK(r.stations[k].coolant_velocity > 1.6 * base.stations[k].coolant_velocity);
    CHECK(r.stations[k].h_coolant > base.stations[k].h_coolant);
    CHECK(r.stations[k].t_wall_hot < base.stations[k].t_wall_hot);
    CHECK(r.coolant_pressure_drop > base.coolant_pressure_drop);
    CHECK(r.energy_balance_residual < 1.0e-8);
  }

  SECTION("a width taper in fraction-of-pitch mode") {
    auto taper = f.spec();
    taper.channel_width_profile = {{0.0, 0.5}, {f.geom.throatPosition(), 0.35}};
    const auto r = solveRegenerativeCooling(f.flow, f.geom, f.ch, taper, f.transport);
    const auto& last = r.stations.back();
    const double pitch = 2.0 * constants::pi * last.radius / taper.num_channels;
    CHECK(last.channel_width == Approx(0.35 * pitch).epsilon(1e-12));
  }

  SECTION("malformed profiles are refused") {
    auto bad = f.spec();
    bad.channel_height_profile = {{0.2, 3e-3}, {0.1, 3e-3}};
    CHECK_THROWS_AS(solveRegenerativeCooling(f.flow, f.geom, f.ch, bad, f.transport), ConfigError);
    bad.channel_height_profile = {{0.1, 3e-3}, {0.2, -1e-3}};
    CHECK_THROWS_AS(solveRegenerativeCooling(f.flow, f.geom, f.ch, bad, f.transport), ConfigError);
    bad = f.spec();
    bad.channel_width_profile = {{0.1, 1.2}};
    CHECK_THROWS_AS(solveRegenerativeCooling(f.flow, f.geom, f.ch, bad, f.transport), ConfigError);
  }
}
