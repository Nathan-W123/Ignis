// SPDX-License-Identifier: MIT
/// \file test_cycle.cpp
/// \brief The turbopump cycle closures: their identities and their limits.
///
/// The pumps and turbines themselves are checked against the RS-25's
/// published turbopump data in tests/validation/test_cycle_validation.cpp;
/// these tests check that each closure balances, carries its pressures and
/// flows through consistently, and moves the way the physics says it must.

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <cmath>

#include "TestHelpers.hpp"
#include "ignis/core/Constants.hpp"
#include "ignis/core/Exceptions.hpp"
#include "ignis/cycle/Cycle.hpp"

using namespace ignis;
using namespace ignis_test;
using Catch::Approx;

namespace {

const PropellantLibrary& library() {
  static const PropellantLibrary lib =
      PropellantLibrary::loadYaml(sourceDir() + "/data/propellants/ignis_propellants.yaml");
  return lib;
}

/// A methane engine's main-chamber flows at a given chamber pressure, scaled
/// with it as a fixed throat would.
CycleInputs methaneEngine(double pc, double dp_jacket) {
  CycleInputs in;
  const double scale = pc / 5.5e6;
  in.chamber_pressure = pc;
  in.oxidizer_mass_flow = 37.07 * scale;
  in.fuel_mass_flow = 10.90 * scale;
  in.oxidizer = &library().at("LOX");
  in.fuel = &library().at("LCH4");
  in.oxidizer_temperature = 90.18;
  in.fuel_temperature = 111.66;
  in.jacket = dp_jacket > 0.0;
  in.jacket_mass_flow = in.fuel_mass_flow;
  in.jacket_pressure_drop = dp_jacket;
  in.jacket_outlet_temperature = 400.0;
  in.jacket_heat = 8.0e6 * scale;
  in.exit_pressure = 35.0e3;
  in.ambient_pressure = 0.0;
  in.thrust_vacuum = 346.0 * constants::g0 * (in.oxidizer_mass_flow + in.fuel_mass_flow);
  in.thrust = in.thrust_vacuum;
  return in;
}

CycleInputs hydrogenEngine(double pc) {
  CycleInputs in;
  const double scale = pc / 5.5e6;
  in.chamber_pressure = pc;
  in.oxidizer_mass_flow = 23.18 * scale;
  in.fuel_mass_flow = 4.214 * scale;
  in.oxidizer = &library().at("LOX");
  in.fuel = &library().at("LH2");
  in.oxidizer_temperature = 90.18;
  in.fuel_temperature = 20.27;
  in.jacket = true;
  in.jacket_mass_flow = in.fuel_mass_flow;
  in.jacket_pressure_drop = 1.37e6 * scale;
  in.jacket_outlet_temperature = 215.0;
  in.jacket_heat = 3.4e6 * in.fuel_mass_flow;  // ~ J/kg from 32 K to 215 K
  in.exit_pressure = 2.0e3;
  in.thrust_vacuum = 453.0 * constants::g0 * (in.oxidizer_mass_flow + in.fuel_mass_flow);
  in.thrust = in.thrust_vacuum;
  return in;
}

/// The reactant enthalpy the combustor balanced, rebuilt from the pumps.
double reactantEnthalpy(const CycleResult& r, const CycleInputs& in, double dh_ox, double dh_f) {
  const PropellantMixture mix(*in.oxidizer, *in.fuel, r.combustor.mixture_ratio,
                              in.oxidizer_temperature, in.fuel_temperature);
  return mix.enthalpy(rocketDatabase()) + mix.oxidizerMassFraction() * dh_ox +
         mix.fuelMassFraction() * dh_f;
}

double workPerKg(const PumpResult& p) { return p.isentropic_head / p.efficiency; }

}  // namespace

TEST_CASE("a gas-generator cycle balances its turbine against its pumps", "[cycle][gg]") {
  const auto& db = rocketDatabase();
  const EquilibriumSolver solver(db);
  CycleSpec spec;
  spec.type = CycleType::kGasGenerator;
  spec.turbine_inlet_temperature = 900.0;
  const auto in = methaneEngine(5.5e6, 2.0e6);
  const auto r = solveCycle(spec, in, solver);
  INFO(r.summary());
  REQUIRE(r.feasible);

  SECTION("power, flows and pressures add up") {
    CHECK(std::abs(r.power_balance_residual) < 1e-10);
    CHECK(r.turbine.mass_flow == Approx(r.gg_mass_flow).epsilon(1e-14));
    CHECK(r.gg_oxidizer_flow + r.gg_fuel_flow == Approx(r.gg_mass_flow).epsilon(1e-14));
    CHECK(r.gg_oxidizer_flow / r.gg_fuel_flow == Approx(r.combustor.mixture_ratio).epsilon(1e-12));
    // The pumps lift the generator's propellant as well as the chamber's.
    CHECK(r.oxidizer_pump.mass_flow == Approx(in.oxidizer_mass_flow + r.gg_oxidizer_flow));
    CHECK(r.fuel_pump.mass_flow == Approx(in.fuel_mass_flow + r.gg_fuel_flow));
    CHECK(r.oxidizer_discharge_pressure == Approx(5.5e6 * 1.25).epsilon(1e-14));
    CHECK(r.fuel_discharge_pressure == Approx(5.5e6 * 1.25 + 2.0e6).epsilon(1e-14));
    CHECK(r.turbine.p_in == Approx(5.5e6));
    CHECK(r.turbine.p_out == Approx(5.5e6 / 20.0));
    CHECK(r.total_mass_flow ==
          Approx(in.oxidizer_mass_flow + in.fuel_mass_flow + r.gg_mass_flow).epsilon(1e-14));
  }
  SECTION("the generator burns at the turbine inlet temperature") {
    // An adiabatic solve at the mixture ratio found, with the reactants
    // carrying their pump work, lands on the target: the balance the solver
    // closed with fixed-temperature equilibria, checked another way.
    const PropellantMixture mix(*in.oxidizer, *in.fuel, r.combustor.mixture_ratio,
                                in.oxidizer_temperature, in.fuel_temperature);
    const double h = reactantEnthalpy(r, in, workPerKg(r.oxidizer_pump), workPerKg(r.fuel_pump));
    const auto st = solver.hp(mix.elementMoles(db), h, r.combustor_pressure, 900.0).state;
    CHECK(st.T == Approx(900.0).epsilon(1e-6));
    CHECK(r.combustor.mixture_ratio < 1.0);  // far fuel-rich
  }
  SECTION("the dumped exhaust costs specific impulse, and its thrust gives some back") {
    const double main_isp = in.thrust_vacuum /
                            ((in.oxidizer_mass_flow + in.fuel_mass_flow) * constants::g0);
    CHECK(r.isp_vacuum_dumped < r.isp_vacuum_delivered);
    CHECK(r.isp_vacuum_delivered < main_isp);
    CHECK(r.exhaust_velocity > 0.0);
    CHECK(r.gg_flow_fraction > 0.01);
    CHECK(r.gg_flow_fraction < 0.06);
    CHECK(r.overall_mixture_ratio <
          in.oxidizer_mass_flow / in.fuel_mass_flow);  // the generator burns fuel-rich
  }
  SECTION("a hotter generator needs less flow; a higher chamber pressure needs more") {
    auto hot = spec;
    hot.turbine_inlet_temperature = 1100.0;
    CHECK(solveCycle(hot, in, solver).gg_flow_fraction < r.gg_flow_fraction);
    CHECK(solveCycle(spec, methaneEngine(10.0e6, 3.6e6), solver).gg_flow_fraction >
          r.gg_flow_fraction);
  }
  SECTION("a configured jacket pressure above the need is delivered, and throttled") {
    auto floored = in;
    floored.fuel_discharge_floor = 15.0e6;
    const auto f = solveCycle(spec, floored, solver);
    CHECK(f.fuel_discharge_pressure == Approx(15.0e6).epsilon(1e-14));
    CHECK(f.fuel_throttle_loss == Approx(15.0e6 - 5.5e6 * 1.25 - 2.0e6).epsilon(1e-12));
    CHECK(f.gg_mass_flow > r.gg_mass_flow);
    CHECK(r.fuel_throttle_loss == 0.0);
  }
}

TEST_CASE("fuel-rich staged combustion solves its turbine pressure ratio",
          "[cycle][staged]") {
  const auto& db = hydrogenDatabase();
  const EquilibriumSolver solver(db);
  CycleSpec spec;
  spec.type = CycleType::kStagedCombustion;
  spec.fuel_rich = true;
  spec.turbine_inlet_temperature = 900.0;
  spec.turbine_efficiency = 0.75;
  const auto in = hydrogenEngine(10.0e6);
  const auto r = solveCycle(spec, in, solver);
  INFO(r.summary());
  REQUIRE(r.feasible);

  SECTION("the pressure chain and the flows are consistent") {
    CHECK(std::abs(r.power_balance_residual) < 1e-9);
    const double pr = r.turbine_pressure_ratio;
    CHECK(pr > 1.0);
    CHECK(r.combustor_pressure == Approx(pr * 10.0e6 * 1.10).epsilon(1e-12));
    CHECK(r.turbine.p_out == Approx(10.0e6 * 1.10).epsilon(1e-12));
    const double feed = r.combustor_pressure * 1.15 + 0.05 * 10.0e6;
    CHECK(r.fuel_discharge_pressure == Approx(feed + in.jacket_pressure_drop).epsilon(1e-12));
    CHECK(r.oxidizer_discharge_pressure == Approx(10.0e6 * 1.25).epsilon(1e-12));
    REQUIRE(r.has_boost_pump);
    CHECK(r.boost_propellant == "oxidizer");
    CHECK(r.boost_pump.p_out == Approx(feed).epsilon(1e-12));
    // All the fuel goes through the preburner with its oxygen.
    CHECK(r.preburner_major_flow == Approx(in.fuel_mass_flow).epsilon(1e-14));
    CHECK(r.turbine.mass_flow ==
          Approx(in.fuel_mass_flow * (1.0 + r.combustor.mixture_ratio)).epsilon(1e-12));
    CHECK(r.boost_pump.mass_flow == Approx(r.preburner_minor_flow).epsilon(1e-14));
    // A closed cycle delivers what its chamber does.
    CHECK(r.isp_vacuum_delivered ==
          Approx(in.thrust_vacuum / ((in.oxidizer_mass_flow + in.fuel_mass_flow) * constants::g0))
              .epsilon(1e-12));
    CHECK(r.max_power_ratio > 1.0);
  }
  SECTION("jacket heat in the fuel means less oxygen for the same temperature") {
    auto cold = in;
    cold.jacket_heat = 0.0;
    CHECK(solveCycle(spec, cold, solver).combustor.mixture_ratio > r.combustor.mixture_ratio);
  }
  SECTION("headroom shrinks with chamber pressure until the cycle cannot close") {
    const auto r20 = solveCycle(spec, hydrogenEngine(20.0e6), solver);
    CHECK(r20.max_power_ratio < r.max_power_ratio);
    const auto r30 = solveCycle(spec, hydrogenEngine(30.0e6), solver);
    INFO(r30.summary());
    CHECK_FALSE(r30.feasible);
    CHECK(r30.max_power_ratio < 1.0);
    CHECK_FALSE(r30.infeasibility.empty());
  }
}

TEST_CASE("oxidiser-rich staged combustion boosts its fuel and not its oxygen",
          "[cycle][staged]") {
  const auto& db = rocketDatabase();
  const EquilibriumSolver solver(db);
  CycleSpec spec;
  spec.type = CycleType::kStagedCombustion;
  spec.fuel_rich = false;
  spec.turbine_inlet_temperature = 750.0;
  spec.turbine_efficiency = 0.75;
  auto in = methaneEngine(15.0e6, 3.0e6);
  const auto r = solveCycle(spec, in, solver);
  INFO(r.summary());
  REQUIRE(r.feasible);
  CHECK(std::abs(r.power_balance_residual) < 1e-9);
  CHECK(r.combustor.mixture_ratio > 20.0);  // far oxidiser-rich
  CHECK(r.preburner_major_flow == Approx(in.oxidizer_mass_flow).epsilon(1e-14));
  CHECK(r.turbine.mass_flow ==
        Approx(in.oxidizer_mass_flow * (1.0 + 1.0 / r.combustor.mixture_ratio)).epsilon(1e-12));
  CHECK(r.oxidizer_discharge_pressure ==
        Approx(r.combustor_pressure * 1.15 + 0.05 * 15.0e6).epsilon(1e-12));
  CHECK(r.fuel_discharge_pressure == Approx(15.0e6 * 1.25 + 3.0e6).epsilon(1e-12));
  if (r.has_boost_pump) CHECK(r.boost_propellant == "fuel");
}

TEST_CASE("an expander closes on the jacket's heat alone", "[cycle][expander]") {
  const auto& db = hydrogenDatabase();
  const EquilibriumSolver solver(db);
  const auto hydrogen = CoolantFluid::load("hydrogen");
  CycleSpec spec;
  spec.type = CycleType::kExpander;
  spec.turbine_efficiency = 0.70;
  auto in = hydrogenEngine(5.5e6);
  in.jacket_fluid = &hydrogen;
  const auto r = solveCycle(spec, in, solver);
  INFO(r.summary());
  REQUIRE(r.feasible);

  SECTION("the turbine runs from the jacket outlet to the injector supply") {
    CHECK(std::abs(r.power_balance_residual) < 1e-9);
    CHECK_FALSE(r.has_combustor);
    CHECK(r.turbine.t_in == Approx(215.0));
    CHECK(r.turbine.p_out == Approx(5.5e6 * 1.25).epsilon(1e-12));
    CHECK(r.turbine.p_in == Approx(r.fuel_discharge_pressure - in.jacket_pressure_drop).epsilon(1e-12));
    CHECK(r.turbine.mass_flow == Approx(in.jacket_mass_flow).epsilon(1e-14));
    CHECK(r.fuel_pump.real_fluid);
  }
  SECTION("a hotter jacket needs less pressure; a cold one cannot drive the pumps") {
    auto hot = in;
    hot.jacket_outlet_temperature = 300.0;
    CHECK(solveCycle(spec, hot, solver).fuel_discharge_pressure < r.fuel_discharge_pressure);
    auto cold = in;
    cold.jacket_outlet_temperature = 60.0;
    const auto c = solveCycle(spec, cold, solver);
    CHECK_FALSE(c.feasible);
    CHECK(c.max_power_ratio < 1.0);
  }
  SECTION("bypassing the turbine costs pressure") {
    auto bypass = spec;
    bypass.turbine_bypass_fraction = 0.05;
    CHECK(solveCycle(bypass, in, solver).fuel_discharge_pressure > r.fuel_discharge_pressure);
  }
  SECTION("without a fuel-cooled jacket there is nothing to drive it") {
    auto none = in;
    none.jacket = false;
    CHECK_THROWS_AS(solveCycle(spec, none, solver), ConfigError);
  }
}

TEST_CASE("pumps follow a propellant's own table and never a surrogate's", "[cycle][pump]") {
  const EquilibriumSolver solver(rocketDatabase());
  CycleSpec spec;
  auto in = methaneEngine(6.0e6, 8.0e6);
  in.fuel = &library().at("RP1");
  in.fuel_temperature = 298.15;
  const auto r = solveCycle(spec, in, solver);
  CHECK(r.oxidizer_pump.real_fluid);
  CHECK_FALSE(r.fuel_pump.real_fluid);
  CHECK_THAT(r.fuel_pump_model, Catch::Matchers::ContainsSubstring("surrogate"));
  // Incompressible: the head is exactly dp / rho at RP-1's measured density.
  CHECK(r.fuel_pump.isentropic_head ==
        Approx((r.fuel_pump.p_out - r.fuel_pump.p_in) / 800.96).epsilon(1e-12));
  // ... and the fuel-rich kerosene generator gas is flagged as approximate.
  bool flagged = false;
  for (const auto& w : r.warnings)
    if (w.find("soot") != std::string::npos) flagged = true;
  CHECK(flagged);
}

TEST_CASE("a generator or preburner reaches its temperature from either side",
          "[cycle][combustor]") {
  const auto& db = hydrogenDatabase();
  const EquilibriumSolver solver(db);
  const auto& ox = library().at("LOX");
  const auto& fu = library().at("LH2");
  for (const bool fuel_rich : {true, false}) {
    INFO((fuel_rich ? "fuel-rich" : "oxidiser-rich"));
    const auto c = combustorAtTemperature(solver, ox, fu, 90.18, 20.27, 10.0e6, 1000.0, fuel_rich);
    const PropellantMixture mix(ox, fu, c.mixture_ratio, 90.18, 20.27);
    const auto st = solver.hp(mix.elementMoles(db), mix.enthalpy(db), 10.0e6, 1000.0).state;
    CHECK(st.T == Approx(1000.0).epsilon(1e-6));
    CHECK(c.temperature == Approx(1000.0).epsilon(1e-12));
    const double mr_st = mix.stoichiometricMixtureRatio();
    CHECK((fuel_rich ? c.mixture_ratio < mr_st : c.mixture_ratio > mr_st));
  }
  // Warmer hydrogen needs less oxygen to reach the same temperature.
  const auto cold = combustorAtTemperature(solver, ox, fu, 90.18, 20.27, 10.0e6, 1000.0, true);
  const auto warm =
      combustorAtTemperature(solver, ox, fu, 90.18, 20.27, 10.0e6, 1000.0, true, 0.0, 3.0e6);
  CHECK(warm.mixture_ratio < cold.mixture_ratio);
  // No mixture of these propellants burns at 6000 K.
  CHECK_THROWS_AS(combustorAtTemperature(solver, ox, fu, 90.18, 20.27, 10.0e6, 6000.0, true),
                  InfeasibleError);
}

TEST_CASE("cycle inputs are checked", "[cycle][errors]") {
  const EquilibriumSolver solver(rocketDatabase());
  const auto in = methaneEngine(5.5e6, 2.0e6);
  CycleSpec bad;
  bad.turbine_efficiency = 0.0;
  CHECK_THROWS_AS(solveCycle(bad, in, solver), ConfigError);
  CycleSpec above;
  above.gas_generator_pressure = 6.0e6;
  CHECK_THROWS_AS(solveCycle(above, in, solver), ConfigError);
  CHECK(parseCycleType("expander") == CycleType::kExpander);
  CHECK_THROWS_AS(parseCycleType("tap_off"), ConfigError);
}
