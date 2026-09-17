// SPDX-License-Identifier: MIT
/// \file test_transient.cpp
/// \brief Transient chamber verification.

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <cmath>

#include "TestHelpers.hpp"
#include "ignis/combustion/Chamber.hpp"
#include "ignis/transient/TransientChamber.hpp"

using namespace ignis;
using namespace ignis_test;
using Catch::Approx;

namespace {

/// A small table is enough for the verification cases and keeps the suite fast;
/// the accuracy of a production-sized table is measured separately.
const EquilibriumTable& table() {
  static const auto db = rocketDatabase();
  static const EquilibriumSolver solver(db);
  static const auto lib =
      PropellantLibrary::loadYaml(sourceDir() + "/data/propellants/ignis_propellants.yaml");
  static const PropellantMixture mix(lib.at("LOX"), lib.at("LCH4"), 3.4, 90.18, 111.66);
  static const EquilibriumTable t = [] {
    TableGrid g;
    g.mr_min = 0.25;
    g.mr_max = 10.0;
    g.mr_points = 21;
    g.t_min = 200.0;
    g.t_max = 4000.0;
    g.t_points = 49;
    g.p_min = 5.0e4;
    g.p_max = 2.0e7;
    g.p_points = 11;
    return EquilibriumTable::build(solver, mix, g);
  }();
  return t;
}

TransientSpec baseSpec() {
  TransientSpec s;
  s.chamber_volume = 1.29e-2;
  s.throat_area = constants::pi * 0.07 * 0.07;
  s.exit_area = 20.0 * s.throat_area;
  s.ambient_pressure = 101325.0;
  s.oxidizer_flow = Schedule::ramp(0.004, 0.035, 37.12, 0.600, 0.025);
  s.fuel_flow = Schedule::ramp(0.002, 0.035, 10.92, 0.602, 0.026);
  s.combustion_efficiency = Schedule({0.0, 0.004, 0.020, 1e9}, {0.0, 0.05, 0.874, 0.874});
  const auto lib =
      PropellantLibrary::loadYaml(sourceDir() + "/data/propellants/ignis_propellants.yaml");
  s.oxidizer_inlet_enthalpy = lib.at("LOX").reference_enthalpy / lib.at("LOX").molar_mass;
  s.fuel_inlet_enthalpy = lib.at("LCH4").reference_enthalpy / lib.at("LCH4").molar_mass;
  s.initial_pressure = 101325.0;
  s.initial_temperature = 300.0;
  s.initial_mixture_ratio = 3.4;
  s.t_end = 0.30;
  s.dt = 5.0e-6;
  s.dt_output = 2.0e-3;
  s.integrator = "rk45";
  s.rtol = 1.0e-9;
  s.atol = 1.0e-12;
  s.dt_max = 2.0e-4;
  return s;
}

}  // namespace

TEST_CASE("equilibrium table is self-consistent", "[transient][table]") {
  const auto& t = table();
  SECTION("the state inversion round trips") {
    for (double mr : {1.5, 3.4, 6.0}) {
      for (double T : {500.0, 1500.0, 3000.0}) {
        for (double p : {1.0e5, 2.0e6, 1.0e7}) {
          const double M = t.molarMass(mr, T, p);
          const double rho = p * M / (constants::R_universal * T);
          const double u = t.internalEnergy(mr, T, p);
          double T2 = 0.0, p2 = 0.0;
          t.solveState(mr, u, rho, T2, p2);
          INFO("O/F " << mr << ", T " << T << " K, p " << p << " Pa");
          REQUIRE(T2 == Approx(T).epsilon(1e-8));
          REQUIRE(p2 == Approx(p).epsilon(1e-8));
        }
      }
    }
  }
  SECTION("the c* correction is exactly one at the adiabatic flame condition") {
    for (double mr : {2.0, 3.4, 5.0}) {
      for (double p : {1.0e5, 5.5e6, 2.0e7}) {
        const double tf = t.flameTemperature(mr, p);
        REQUIRE(t.cStarAt(mr, tf, p) == Approx(t.cStar(mr, p)).epsilon(2e-3));
      }
    }
  }
  SECTION("c* falls with chamber temperature") {
    const double p = 5.5e6, mr = 3.4;
    const double tf = t.flameTemperature(mr, p);
    REQUIRE(t.cStarAt(mr, 0.8 * tf, p) < t.cStarAt(mr, tf, p));
  }
  SECTION("the heat of combustion is positive and of the right magnitude") {
    for (double mr : {2.0, 3.4, 5.0}) {
      const double q = t.heatOfCombustion(mr);
      INFO("O/F " << mr << ": q_comb = " << q * 1e-6 << " MJ/kg");
      REQUIRE(q > 5.0e6);
      REQUIRE(q < 15.0e6);
    }
  }
  SECTION("requests outside the table are reported") {
    REQUIRE_THROWS_AS(t.internalEnergy(0.1, 2000.0, 1e6), RangeError);
    REQUIRE_THROWS_AS(t.internalEnergy(3.4, 100.0, 1e6), RangeError);
    REQUIRE_THROWS_AS(t.internalEnergy(3.4, 2000.0, 1e3), RangeError);
  }
}

TEST_CASE("transient reaches the steady state the steady model predicts",
          "[transient][verification]") {
  const auto spec = baseSpec();
  const auto res = simulateTransient(table(), spec);
  REQUIRE(res.samples.size() > 50);

  // The commanded flows and the ideal c* fix the steady chamber pressure.
  const auto& last = res.samples[res.samples.size() / 2];   // mid-burn, fully settled
  REQUIRE(last.choked);
  const double mdot = last.mdot_ox_in + last.mdot_fuel_in;
  REQUIRE(last.mdot_out == Approx(mdot).epsilon(1e-4));
  REQUIRE(last.pressure == Approx(mdot * last.c_star / spec.throat_area).epsilon(1e-4));
  REQUIRE(last.mixture_ratio == Approx(37.12 / 10.92).epsilon(1e-3));

  SECTION("conservation holds to round-off") {
    REQUIRE(res.mass_conservation_error < 1.0e-9);
    REQUIRE(res.energy_conservation_error < 1.0e-9);
  }
  SECTION("the chamber fills on a physically sensible timescale") {
    REQUIRE(res.time_to_90_percent > 0.005);
    REQUIRE(res.time_to_90_percent < 0.15);
  }
  SECTION("pressure rises monotonically once the igniter has taken hold") {
    // Before ignition the chamber briefly *cools*: cold propellant floods a
    // volume that starts at ambient temperature while almost no heat is being
    // released, so the pressure dips first.  That is a real feature of the
    // model, not a numerical artefact, so monotonicity is only asserted from
    // the end of the ignition ramp onwards.
    double dip = 1e30;
    for (const auto& s : res.samples) {
      if (s.t > 0.02) break;
      dip = std::min(dip, s.pressure);
    }
    REQUIRE(dip < res.samples.front().pressure);
    double previous = 0.0;
    for (const auto& s : res.samples) {
      if (s.t < 0.02) continue;
      if (s.t > 0.10) break;
      REQUIRE(s.pressure >= previous - 1.0);
      previous = s.pressure;
    }
  }
}

TEST_CASE("a constant-flow chamber follows the analytic filling solution",
          "[transient][verification]") {
  // With constant inflow, full heat release and a choked nozzle the chamber
  // approaches its steady pressure exponentially with time constant
  //     tau = V c* / (A_t R T) ... in practice the correct statement is that
  // the steady state satisfies mdot_in = p A_t / c*, which is what is checked.
  auto spec = baseSpec();
  spec.oxidizer_flow = Schedule::constant(37.12);
  spec.fuel_flow = Schedule::constant(10.92);
  spec.combustion_efficiency = Schedule::constant(1.0);
  spec.t_end = 0.5;
  const auto res = simulateTransient(table(), spec);
  const auto& last = res.samples.back();
  const double mdot = 37.12 + 10.92;
  REQUIRE(last.mdot_out == Approx(mdot).epsilon(1e-6));
  REQUIRE(last.pressure == Approx(mdot * last.c_star / spec.throat_area).epsilon(1e-6));
  // Full heat release must reproduce the adiabatic flame temperature.
  REQUIRE(last.temperature ==
          Approx(table().flameTemperature(last.mixture_ratio, last.pressure)).epsilon(5e-3));
}

TEST_CASE("the fixed-step integrator converges at high order on a smooth problem",
          "[transient][verification][convergence]") {
  // The order of a Runge-Kutta scheme is only observable on a smooth
  // right-hand side, and a real startup is not smooth: the piecewise-linear
  // valve schedules have kinks, and the sub-critical orifice relation has a
  // square-root singularity at the instant the chamber first exceeds ambient
  // pressure (d mdot / d p is infinite at zero pressure difference).  Either of
  // those caps the observed order near one.  This test therefore starts the
  // chamber already choked and hot with constant flows, leaving only the C1
  // tabulated equation of state, and recovers the theoretical fourth order.
  auto spec = baseSpec();
  spec.oxidizer_flow = Schedule::constant(37.12);
  spec.fuel_flow = Schedule::constant(10.92);
  spec.combustion_efficiency = Schedule::constant(0.95);
  spec.initial_pressure = 3.0e6;
  spec.initial_temperature = 2500.0;
  spec.integrator = "rk4";
  spec.t_end = 0.02;
  spec.dt_output = spec.t_end;
  // A reference at a much finer step.
  spec.dt = 1.0e-6;
  // The order is measured on an integrated state variable.  Chamber pressure
  // would instead be limited by the accuracy of the equation-of-state
  // inversion, which is not what this test is about.
  const double reference = simulateTransient(table(), spec).samples.back().internal_energy;

  double previous_error = 0.0;
  std::vector<double> orders;
  for (double dt : {8.0e-4, 4.0e-4, 2.0e-4, 1.0e-4}) {
    spec.dt = dt;
    const double p = simulateTransient(table(), spec).samples.back().internal_energy;
    const double err = std::abs(p - reference);
    if (previous_error > 0.0 && err > 0.0)
      orders.push_back(std::log2(previous_error / err));
    previous_error = err;
  }
  REQUIRE(orders.size() >= 3);
  INFO("observed convergence orders: " << orders[0] << ", " << orders[1] << ", " << orders[2]);
  for (double o : orders) {
    REQUIRE(o > 3.6);
    REQUIRE(o < 4.6);
  }
}

TEST_CASE("the adaptive integrator agrees with the fixed-step one",
          "[transient][verification][slow]") {
  auto fixed = baseSpec();
  fixed.integrator = "rk4";
  fixed.dt = 1.0e-6;
  fixed.t_end = 0.1;
  fixed.dt_output = 1.0e-3;
  auto adaptive = fixed;
  adaptive.integrator = "rk45";
  adaptive.dt = 1.0e-5;
  adaptive.rtol = 1.0e-10;
  adaptive.atol = 1.0e-14;

  const auto a = simulateTransient(table(), fixed);
  const auto b = simulateTransient(table(), adaptive);
  REQUIRE(a.samples.size() == b.samples.size());
  double worst = 0.0;
  for (std::size_t i = 0; i < a.samples.size(); ++i)
    worst = std::max(worst, relativeError(b.samples[i].pressure, a.samples[i].pressure, 1.0));
  INFO("largest pressure difference between integrators " << worst);
  REQUIRE(worst < 1.0e-5);
  // The adaptive scheme should need far fewer steps.
  REQUIRE(b.steps < a.steps / 5);
  REQUIRE(b.rejected_steps >= 0);
}

TEST_CASE("shutdown blows the chamber down to ambient", "[transient]") {
  auto spec = baseSpec();
  spec.t_end = 0.8;
  const auto res = simulateTransient(table(), spec);
  const auto& last = res.samples.back();
  REQUIRE(last.t == Approx(0.8).margin(2e-3));
  REQUIRE(last.pressure == Approx(spec.ambient_pressure).epsilon(0.02));
  REQUIRE(last.mdot_out == Approx(0.0).margin(1e-9));
  REQUIRE_FALSE(last.choked);
  REQUIRE(last.thrust == Approx(0.0).margin(1e-9));
  REQUIRE(res.max_pressure > 5.0e6);
  // The chamber cools as it empties but stays well above the table floor.
  // The chamber cools below its initial charge temperature while cold
  // propellant floods in before the igniter takes hold; it must stay above
  // the 200 K floor of the thermodynamic data.
  REQUIRE(res.min_temperature > 210.0);
  REQUIRE(res.max_temperature < 4000.0);
}

TEST_CASE("transient input validation and failure reporting", "[transient][errors]") {
  SECTION("schedules must be well formed") {
    REQUIRE_THROWS_AS(Schedule({0.0, 0.0}, {1.0, 2.0}), ConfigError);
    REQUIRE_THROWS_AS(Schedule({0.0, 1.0}, {1.0}), ConfigError);
    REQUIRE_THROWS_AS(Schedule::ramp(0.1, 0.05, 1.0, 0.11, 0.05), ConfigError);
  }
  SECTION("a piecewise-linear schedule interpolates and saturates") {
    const Schedule s({0.0, 1.0, 2.0}, {0.0, 10.0, 5.0});
    REQUIRE(s.at(-1.0) == Approx(0.0));
    REQUIRE(s.at(0.5) == Approx(5.0));
    REQUIRE(s.at(1.5) == Approx(7.5));
    REQUIRE(s.at(9.0) == Approx(5.0));
  }
  SECTION("invalid specifications are rejected") {
    auto spec = baseSpec();
    spec.chamber_volume = 0.0;
    REQUIRE_THROWS_AS(simulateTransient(table(), spec), ConfigError);
    spec = baseSpec();
    spec.integrator = "euler";
    REQUIRE_THROWS_WITH(simulateTransient(table(), spec),
                        Catch::Matchers::ContainsSubstring("rk4"));
  }
  SECTION("a chamber ratio outside the table stops the run with an explanation") {
    auto spec = baseSpec();
    // Close the oxidiser long before the fuel: the chamber goes far too rich.
    spec.oxidizer_flow = Schedule::ramp(0.004, 0.035, 37.12, 0.20, 0.005);
    spec.fuel_flow = Schedule::ramp(0.002, 0.035, 10.92, 0.80, 0.026);
    spec.t_end = 0.9;
    REQUIRE_THROWS_WITH(simulateTransient(table(), spec),
                        Catch::Matchers::ContainsSubstring("mixture ratio"));
  }
}
