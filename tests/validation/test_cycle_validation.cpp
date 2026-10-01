// SPDX-License-Identifier: MIT
/// \file test_cycle_validation.cpp
/// \brief Pumps and turbines against the RS-25's published turbopump data.
///
/// SOURCE.  "Space Shuttle Main Engine Orientation", Boeing/Rocketdyne
/// Propulsion & Power (BC98-04, 1998), at 104.5 % of rated power level -- the
/// operating point Ignis's RS-25 preset reproduces: the Block IIA "Key
/// Performance Parameters" page of each turbopump (printed pp. 53, 57, 59,
/// 67), the propellant flow schematic (p. 19) and the propellant flow
/// analysis (pp. 18-21).  Values are the manual's, in its units, converted
/// below.
///
/// Each turbopump is checked twice and independently.  From the PUMP side:
/// mass flow, inlet and discharge pressure and pump efficiency give the power
/// the pumps absorb.  From the TURBINE side: turbine flow, inlet temperature,
/// pressure ratio and efficiency, with the equilibrium preburner gas at the
/// mixture ratio the manual's flows imply, give the power the turbine
/// delivers.  In steady running both equal the published turbine horsepower,
/// less the bearing and seal losses the manual does not tabulate; nothing is
/// fitted.  The tolerances were set before the first run.
#include <cmath>
#include <utility>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "TestHelpers.hpp"
#include "ignis/cycle/Turbomachinery.hpp"

using namespace ignis;
using namespace ignis_test;

namespace {

constexpr double kPsi = 6894.757293168361;      // Pa
constexpr double kPound = 0.45359237;           // kg
constexpr double kHorsepower = 745.69987158227; // W, mechanical hp
double fahrenheit(double f) { return (f - 32.0) * 5.0 / 9.0 + 273.15; }

/// Propellant inlet states.  The manual gives pressures, not temperatures:
/// the external tank holds its liquids near their normal boiling points, and
/// the low-pressure pumps add a fraction of a kelvin, so each high-pressure
/// pump is fed at the low-pressure pump's computed outlet temperature.
constexpr double kLh2Tank = 20.3;   // K
constexpr double kLoxTank = 90.2;   // K

}  // namespace

TEST_CASE("pumps reproduce the RS-25's turbopump horsepower", "[validation][cycle]") {
  const auto hydrogen = CoolantFluid::load("hydrogen");
  const auto oxygen = CoolantFluid::load("oxygen");

  // Low-pressure fuel pump, then the high-pressure one it feeds.
  const auto lpfp = pumpLiquid(&hydrogen, 0.0, 155.0 * kPound, 30.0 * kPsi, kLh2Tank,
                               290.0 * kPsi, 0.713);
  const auto hpfp = pumpLiquid(&hydrogen, 0.0, 155.0 * kPound, 250.0 * kPsi, lpfp.t_out,
                               5950.0 * kPsi, 0.750);
  // Low-pressure oxidiser pump; the high-pressure turbopump's main pump and
  // its preburner boost pump share one turbine.
  const auto lpop = pumpLiquid(&oxygen, 0.0, 935.0 * kPound, 100.0 * kPsi, kLoxTank,
                               417.0 * kPsi, 0.677);
  const auto hpop = pumpLiquid(&oxygen, 0.0, 1122.0 * kPound, 380.0 * kPsi, lpop.t_out,
                               4045.0 * kPsi, 0.718);
  const auto pbp = pumpLiquid(&oxygen, 0.0, 111.0 * kPound, 3910.0 * kPsi, hpop.t_out,
                              6970.0 * kPsi, 0.758);

  const double hpftp = hpfp.power / kHorsepower;
  const double hpotp = (hpop.power + pbp.power) / kHorsepower;
  const double lpftp = lpfp.power / kHorsepower;
  const double lpotp = lpop.power / kHorsepower;
  INFO("pump-side power, hp: HPFTP " << hpftp << " (turbine 63080), HPOTP " << hpotp
       << " (turbine 22880), LPFTP " << lpftp << " (turbine 3330), LPOTP " << lpotp
       << " (turbine 1614); LH2 into the HPFTP at " << lpfp.t_out << " K, out at " << hpfp.t_out
       << " K");
  CHECK(hpfp.real_fluid);

  SECTION("the high-pressure pumps") {
    // The two pumps that set the cycle, at 41 and 48 MPa.  A turbine also pays
    // for bearings and seals, so the pump side should sit at or a little
    // below the turbine figure.
    CHECK(hpftp == Catch::Approx(63080.0).epsilon(0.06));
    CHECK(hpotp == Catch::Approx(22880.0).epsilon(0.06));
  }
  SECTION("the low-pressure pumps") {
    // Smaller machines, larger relative losses, and inlet temperatures the
    // manual does not give: a looser check.
    CHECK(lpftp == Catch::Approx(3330.0).epsilon(0.15));
    CHECK(lpotp == Catch::Approx(1614.0).epsilon(0.15));
  }
  SECTION("liquid hydrogen is compressible enough to matter at 41 MPa") {
    // dp/rho at the inlet density overstates the HPFTP's work; integrating
    // along the table is what lands it on the turbine figure.
    const double rho_in = hydrogen.at(lpfp.t_out, 250.0 * kPsi).rho;
    const auto naive = pumpLiquid(nullptr, rho_in, 155.0 * kPound, 250.0 * kPsi, lpfp.t_out,
                                  5950.0 * kPsi, 0.750);
    INFO("incompressible at the inlet density: " << naive.power / kHorsepower << " hp");
    CHECK(naive.power > 1.05 * hpfp.power);
  }
}

TEST_CASE("turbines on preburner gas reproduce the RS-25's turbine horsepower",
          "[validation][cycle]") {
  // The preburner gas is fixed by the manual's own flows, not by a combustion
  // solve.  "76 percent of the hydrogen flow and 11 percent of the oxygen
  // flow are injected into two preburners ... the fuel preburner receiving,
  // by far, the larger share (50 percent versus 26 percent)" of the 155 lb/s
  // of hydrogen, and each turbine's flow is its preburner's whole output, so
  //     fuel preburner      O/F = (145 - 0.50 x 155) / (0.50 x 155) = 0.871
  //     oxidiser preburner  O/F = ( 62 - 0.26 x 155) / (0.26 x 155) = 0.538
  // -- "less than one pound of LOX to one pound of hydrogen", as the manual
  // says.  At those mixture ratios and the published turbine inlet
  // temperatures the equilibrium gas is hydrogen and steam.
  //
  // Solving the mixture ratio from the temperature instead, with both
  // propellants at their tank temperatures, gives 1.02 and 0.75 and turbine
  // powers 11 % and 14 % low: the real preburner hydrogen arrives warm,
  // having partly passed through the 1,080 nozzle coolant tubes, so it needs
  // less oxygen to reach the same temperature (docs/validation.md, Sec. 5d, and the
  // preburner energy-balance test below).
  const auto& db = hydrogenDatabase();
  const EquilibriumSolver solver(db);
  const auto lib = PropellantLibrary::loadYaml(sourceDir() + "/data/propellants/ignis_propellants.yaml");
  const GasMixture mix(db);
  // The manual gives the high-pressure oxidiser turbine's inlet pressure,
  // 4785 psia, and no other; at ~1000 K the hydrogen-steam equilibrium does
  // not depend on it, and an ideal-gas turbine sees only the ratio.
  const double p_pb = 4785.0 * kPsi;
  const double h2 = 155.0 * kPound;
  auto gas = [&](double turbine_flow, double fuel_share, double t_in) {
    const double fuel = fuel_share * h2;
    const double mr = (turbine_flow - fuel) / fuel;
    const PropellantMixture m(lib.at("LOX"), lib.at("LH2"), mr, 90.18, 20.27);
    return std::make_pair(mr, solver.tp(m.elementMoles(db), t_in, p_pb).state);
  };
  const auto [mr_f, fpb] = gas(145.0 * kPound, 0.50, fahrenheit(1330.0));
  const auto [mr_o, opb] = gas(62.0 * kPound, 0.26, fahrenheit(870.0));
  const auto hpft = gasTurbine(mix, fpb.n, 145.0 * kPound, fpb.T, p_pb, 1.50, 0.811);
  const auto hpot = gasTurbine(mix, opb.n, 62.0 * kPound, opb.T, p_pb, 1.53, 0.746);
  INFO("preburner O/F from the published flows: fuel " << mr_f << ", oxidiser " << mr_o
       << ";  turbine power, hp: HPFTP " << hpft.power / kHorsepower << " (published 63080), HPOTP "
       << hpot.power / kHorsepower << " (published 22880)");
  CHECK(mr_f == Catch::Approx(0.871).epsilon(1e-3));
  CHECK(mr_o == Catch::Approx(0.538).epsilon(2e-3));
  CHECK(hpft.power / kHorsepower == Catch::Approx(63080.0).epsilon(0.08));
  CHECK(hpot.power / kHorsepower == Catch::Approx(22880.0).epsilon(0.08));
}

TEST_CASE("a hydrogen expander turbine reproduces the RS-25's low-pressure fuel turbine",
          "[validation][cycle]") {
  // The low-pressure fuel turbopump is driven by the 29 lb/s of hydrogen that
  // has cooled the main combustion chamber -- an expander turbine inside a
  // staged-combustion engine.  Its Key Performance Parameters page (p. 57):
  // 29 lb/s, pressure ratio 1.30, efficiency 58.0 %, 3,330 hp, and an inlet
  // temperature of -17 F with no inlet pressure.  The flow schematic (p. 19)
  // labels the same 29 lb/s at the turbine inlet "+17 F, 4,217 psia", and
  // the chamber's coolant leaves at +17 F (MCC operating parameters, p. 45).
  // The page's -17 F is taken to be a sign slip: between the chamber and the
  // turbine the hydrogen is throttled from 4,441 to 4,217 psia, and hydrogen
  // above its ~200 K inversion temperature warms when throttled, so it cannot
  // arrive 19 K colder than it left.  At -17 F the power would be ~6 % lower.
  // At this pressure hydrogen is 22 % off the ideal gas (Z = 1.22), which the
  // real-fluid table carries.
  const auto hydrogen = CoolantFluid::load("hydrogen");
  const auto lpft = fluidTurbine(hydrogen, 29.0 * kPound, fahrenheit(17.0), 4217.0 * kPsi, 1.30,
                                 0.580);
  INFO("LPFT power " << lpft.power / kHorsepower << " hp (published 3330); outlet "
       << lpft.t_out << " K");
  CHECK(lpft.power / kHorsepower == Catch::Approx(3330.0).epsilon(0.08));
  // The turbine cools the gas it expands.
  CHECK(lpft.t_out < lpft.t_in);
}

TEST_CASE("preburner energy balances reproduce the RS-25's preburner mixture ratios",
          "[validation][cycle]") {
  // The flow schematic (p. 19) gives each preburner's inlet streams and its
  // outlet temperature:
  //   fuel preburner      77 lb/s H2 at -193 F, 5,310 psia; 67 lb/s O2 at
  //                       -252 F, 5,636 psia; out at 1,310 F, 4,793 psia
  //   oxidiser preburner  42 lb/s H2 at -193 F; 25 lb/s O2 at -252 F,
  //                       5,734 psia; out at 871 F, 4,812 psia
  // so O/F 0.870 and 0.595.  Ignis solves the mixture ratio that burns at the
  // published outlet temperature, with each propellant carrying the enthalpy
  // it has gained since the tank (from the real-fluid tables), and the
  // published flows are the check.  Taken from 2 bar rather than 1 atm, the
  // tank liquids' enthalpy is off by v dp, ~1 kJ/kg of ~3 MJ/kg gained.  The
  // tolerance, 5 %, is what a 20 K uncertainty in the outlet temperature
  // means for the mixture ratio here (~500 K per unit O/F).
  const auto hydrogen = CoolantFluid::load("hydrogen");
  const auto oxygen = CoolantFluid::load("oxygen");
  const auto& db = hydrogenDatabase();
  const EquilibriumSolver solver(db);
  const auto lib = PropellantLibrary::loadYaml(sourceDir() + "/data/propellants/ignis_propellants.yaml");
  const double tank = 0.2e6;
  const double dh_fuel = hydrogen.at(fahrenheit(-193.0), 5310.0 * kPsi).h -
                         hydrogen.at(lib.at("LH2").reference_temperature, tank).h;
  auto dh_ox = [&](double psia) {
    return oxygen.at(fahrenheit(-252.0), psia * kPsi).h -
           oxygen.at(lib.at("LOX").reference_temperature, tank).h;
  };
  const auto fpb = combustorAtTemperature(solver, lib.at("LOX"), lib.at("LH2"), 90.18, 20.27,
                                          4793.0 * kPsi, fahrenheit(1310.0), true, dh_ox(5636.0),
                                          dh_fuel);
  const auto opb = combustorAtTemperature(solver, lib.at("LOX"), lib.at("LH2"), 90.18, 20.27,
                                          4812.0 * kPsi, fahrenheit(871.0), true, dh_ox(5734.0),
                                          dh_fuel);
  INFO("hydrogen gains " << dh_fuel * 1e-6 << " MJ/kg before the preburners;  O/F: fuel preburner "
       << fpb.mixture_ratio << " (published flows 0.870), oxidiser preburner "
       << opb.mixture_ratio << " (0.595)");
  CHECK(fpb.mixture_ratio == Catch::Approx(67.0 / 77.0).epsilon(0.05));
  CHECK(opb.mixture_ratio == Catch::Approx(25.0 / 42.0).epsilon(0.05));
}
