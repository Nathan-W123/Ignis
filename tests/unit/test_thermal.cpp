// SPDX-License-Identifier: MIT
/// \file test_thermal.cpp
/// \brief Heat transfer, wall conduction, coolant properties and cooling.

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <cmath>

#include "TestHelpers.hpp"
#include "ignis/combustion/Chamber.hpp"
#include "ignis/thermal/CoolantFluid.hpp"
#include "ignis/thermal/HeatTransfer.hpp"
#include "ignis/thermal/RegenerativeCooling.hpp"

using namespace ignis;
using namespace ignis_test;
using Catch::Approx;

TEST_CASE("cylindrical wall conduction matches the analytical solution", "[thermal]") {
  const double k = 320.0;
  SECTION("resistance equals the closed-form expression") {
    const double ri = 0.05, t = 1.0e-3;
    const double expected = ri * std::log((ri + t) / ri) / k;
    REQUIRE(cylindricalWallResistance(ri, t, k) == Approx(expected).epsilon(1e-14));
  }
  SECTION("the thin-wall limit recovers the planar result") {
    const double ri = 1.0;
    for (double t : {1.0e-3, 1.0e-4, 1.0e-5}) {
      const double cyl = cylindricalWallResistance(ri, t, k);
      const double plane = planarWallResistance(t, k);
      INFO("t/ri = " << t / ri);
      REQUIRE(relativeError(cyl, plane) < t / ri);
    }
  }
  SECTION("the temperature drop reproduces a hand calculation") {
    // A 1 mm copper wall at 50 MW/m^2 on a 50 mm radius.
    const double q = 50.0e6, ri = 0.05, t = 1.0e-3;
    const double dT = q * cylindricalWallResistance(ri, t, k);
    // Analytical: dT = q ri ln(ro/ri) / k
    REQUIRE(dT == Approx(q * ri * std::log(0.051 / 0.05) / k).epsilon(1e-12));
    REQUIRE(dT > 150.0);
    REQUIRE(dT < 160.0);
  }
  SECTION("invalid inputs are rejected") {
    REQUIRE_THROWS_AS(cylindricalWallResistance(0.0, 1e-3, k), ConfigError);
    REQUIRE_THROWS_AS(cylindricalWallResistance(0.05, 0.0, k), ConfigError);
    REQUIRE_THROWS_AS(cylindricalWallResistance(0.05, 1e-3, 0.0), ConfigError);
  }
}

TEST_CASE("the recovery temperature has the right limits", "[thermal]") {
  const double gamma = 1.2, Pr = 0.6;
  REQUIRE(recoveryTemperature(3000.0, 0.0, gamma, Pr) == Approx(3000.0));
  // At M = 1 the recovery temperature sits between static and stagnation.
  const double Tstat = 3000.0;
  const double Tstag = Tstat * (1.0 + 0.5 * (gamma - 1.0));
  const double Taw = recoveryTemperature(Tstat, 1.0, gamma, Pr);
  REQUIRE(Taw > Tstat);
  REQUIRE(Taw < Tstag);
  // The recovery factor is Pr^(1/3).
  REQUIRE((Taw - Tstat) / (Tstag - Tstat) == Approx(std::cbrt(Pr)).epsilon(1e-12));
  // With Pr = 1 the recovery temperature is the stagnation temperature.
  REQUIRE(recoveryTemperature(Tstat, 3.0, gamma, 1.0) ==
          Approx(Tstat * (1.0 + 0.5 * (gamma - 1.0) * 9.0)).epsilon(1e-12));
}

TEST_CASE("the Bartz correlation scales as documented", "[thermal]") {
  BartzReference ref;
  ref.throat_diameter = 0.14;
  ref.curvature_radius = 0.10;
  ref.chamber_pressure = 5.5e6;
  ref.c_star = 1770.0;
  ref.chamber_temperature = 3500.0;
  ref.viscosity = 1.0e-4;
  ref.cp = 2300.0;
  ref.prandtl = 0.6;

  const double h_throat = bartzFilmCoefficient(ref, 1.0, 1.0, 1.2, 800.0);
  REQUIRE(h_throat > 0.0);

  SECTION("the film coefficient falls as (At/A)^0.9 away from the throat") {
    const double h4 = bartzFilmCoefficient(ref, 4.0, 1.0, 1.2, 800.0);
    REQUIRE(h4 / h_throat == Approx(std::pow(0.25, 0.9)).epsilon(1e-12));
  }
  SECTION("the film coefficient scales as p_c^0.8") {
    auto hot = ref;
    hot.chamber_pressure = 2.0 * ref.chamber_pressure;
    REQUIRE(bartzFilmCoefficient(hot, 1.0, 1.0, 1.2, 800.0) / h_throat ==
            Approx(std::pow(2.0, 0.8)).epsilon(1e-12));
  }
  SECTION("the multiplier is applied linearly") {
    auto scaled = ref;
    scaled.multiplier = 1.37;
    REQUIRE(bartzFilmCoefficient(scaled, 1.0, 1.0, 1.2, 800.0) ==
            Approx(1.37 * h_throat).epsilon(1e-12));
  }
  SECTION("a hotter wall reduces the film coefficient through sigma") {
    REQUIRE(bartzFilmCoefficient(ref, 1.0, 1.0, 1.2, 1200.0) < h_throat);
  }
  SECTION("dimensional consistency: SI in, W/(m^2 K) out, magnitude plausible") {
    REQUIRE(h_throat > 3.0e3);
    REQUIRE(h_throat < 1.0e5);
  }
  SECTION("invalid inputs are rejected") {
    auto bad = ref;
    bad.c_star = 0.0;
    REQUIRE_THROWS_AS(bartzFilmCoefficient(bad, 1.0, 1.0, 1.2, 800.0), ConfigError);
    REQUIRE_THROWS_AS(bartzFilmCoefficient(ref, 0.5, 1.0, 1.2, 800.0), ConfigError);
  }
}

TEST_CASE("gray-gas radiation is separate from convection and has the right form",
          "[thermal]") {
  REQUIRE(grayGasRadiation(3000.0, 800.0, 0.0, 0.3) == 0.0);
  const double q = grayGasRadiation(3000.0, 800.0, 0.2, 0.3);
  REQUIRE(q == Approx(constants::sigma_SB * 0.2 * 0.3 *
                      (std::pow(3000.0, 4) - std::pow(800.0, 4))).epsilon(1e-12));
  REQUIRE(grayGasRadiation(1000.0, 1000.0, 0.2, 0.3) == Approx(0.0).margin(1e-9));
  REQUIRE_THROWS_AS(grayGasRadiation(3000.0, 800.0, 1.5, 0.3), ConfigError);
}

TEST_CASE("material library carries cited, physically sensible data", "[thermal]") {
  const auto lib = MaterialLibrary::loadDefault();
  REQUIRE(lib.names().size() >= 4);
  for (const auto& name : lib.names()) {
    const auto& m = lib.at(name);
    INFO("material " << name);
    REQUIRE_FALSE(m.source.empty());
    REQUIRE(m.conductivity(300.0) > 5.0);
    REQUIRE(m.conductivity(300.0) < 500.0);
    REQUIRE(m.max_temperature > 500.0);
    REQUIRE(m.emissivity > 0.0);
    REQUIRE(m.emissivity <= 1.0);
    REQUIRE(m.density > 1000.0);
    REQUIRE(m.conductivity(m.valid_min) > 0.0);
    REQUIRE(m.conductivity(m.valid_max) > 0.0);
  }
  REQUIRE(lib.at("CuCrZr").conductivity(300.0) > lib.at("Inconel718").conductivity(300.0));
  REQUIRE_THROWS_WITH(lib.at("Unobtainium"),
                      Catch::Matchers::ContainsSubstring("unknown wall material"));
}

TEST_CASE("coolant tables reproduce their reference equations of state",
          "[thermal][coolant][verification]") {
  const ReferenceTable ref(referenceDir() + "/coolant_reference.csv");
  double worst_rho = 0.0, worst_cp = 0.0, worst_mu = 0.0, worst_k = 0.0;
  std::string last;
  CoolantFluid fluid = CoolantFluid::load("methane");
  for (std::size_t r = 0; r < ref.rows(); ++r) {
    const auto& name = ref.text("fluid", r);
    if (name != last) { fluid = CoolantFluid::load(name); last = name; }
    const auto st = fluid.at(ref.num("T", r), ref.num("p", r));
    worst_rho = std::max(worst_rho, relativeError(st.rho, ref.num("rho", r)));
    worst_cp = std::max(worst_cp, relativeError(st.cp, ref.num("cp", r)));
    worst_mu = std::max(worst_mu, relativeError(st.mu, ref.num("mu", r)));
    worst_k = std::max(worst_k, relativeError(st.k, ref.num("k", r)));
  }
  INFO("interpolation error: rho " << worst_rho << ", cp " << worst_cp << ", mu " << worst_mu
       << ", k " << worst_k);
  // Bilinear interpolation on a 10 K by geometric-pressure grid.
  REQUIRE(worst_rho < 0.02);
  REQUIRE(worst_cp < 0.06);
  REQUIRE(worst_mu < 0.03);
  REQUIRE(worst_k < 0.04);
}

TEST_CASE("Peng-Robinson agrees with the reference tables where a cubic can",
          "[thermal][coolant]") {
  const ReferenceTable ref(referenceDir() + "/coolant_reference.csv");
  double worst_dilute = 0.0, worst_dense = 0.0;
  for (std::size_t r = 0; r < ref.rows(); ++r) {
    if (ref.text("fluid", r) != "methane") continue;
    const PengRobinsonFluid pr(ref.num("Tcrit", r), ref.num("pcrit", r),
                               ref.num("acentric", r), ref.num("molar_mass", r));
    const double T = ref.num("T", r), p = ref.num("p", r);
    const double rho_ref = ref.num("rho", r);
    const double e = relativeError(pr.density(T, p), rho_ref);
    if (T > 250.0) worst_dilute = std::max(worst_dilute, e);
    else worst_dense = std::max(worst_dense, e);
  }
  INFO("Peng-Robinson density error: " << worst_dilute << " above 250 K, " << worst_dense
       << " below");
  // Documented in CoolantFluid.hpp: a cubic equation of state is good to a
  // couple of percent in the gas-like region and degrades badly in the dense
  // liquid-like region, which is why Ignis interpolates reference tables.
  REQUIRE(worst_dilute < 0.05);
  REQUIRE(worst_dense < 0.20);
  REQUIRE(worst_dense > worst_dilute);
}

TEST_CASE("coolant property lookups report their own limits", "[thermal][coolant][errors]") {
  const auto methane = CoolantFluid::load("methane");
  REQUIRE(methane.criticalTemperature() == Approx(190.564).epsilon(1e-4));
  REQUIRE(methane.criticalPressure() == Approx(4.5992e6).epsilon(1e-4));

  SECTION("out-of-range temperature and pressure") {
    REQUIRE_THROWS_AS(methane.at(50.0, 1.0e7), RangeError);
    REQUIRE_THROWS_AS(methane.at(500.0, 1.0e9), RangeError);
  }
  SECTION("the melting line raises the lower temperature limit with pressure") {
    double lo_low = 0.0, hi_low = 0.0, lo_high = 0.0, hi_high = 0.0;
    methane.validTemperatureRange(1.0e6, lo_low, hi_low);
    methane.validTemperatureRange(3.5e7, lo_high, hi_high);
    INFO("valid T from " << lo_low << " K at 1 MPa to " << lo_high << " K at 35 MPa");
    REQUIRE(lo_high > lo_low);
    REQUIRE_THROWS_WITH(methane.at(lo_high - 15.0, 3.5e7),
                        Catch::Matchers::ContainsSubstring("melting line"));
  }
  SECTION("phase classification") {
    // Well above the critical pressure there is no phase boundary.
    REQUIRE(methane.at(200.0, 1.5e7).phase == CoolantPhase::kSupercriticalPressure);
    // Below the critical pressure and temperature the fluid is a liquid or a
    // vapour depending on which side of the saturation line it is on.
    const double tsat = methane.saturationTemperature(2.0e6);
    REQUIRE(tsat > 100.0);
    REQUIRE(tsat < methane.criticalTemperature());
    REQUIRE(methane.at(tsat - 10.0, 2.0e6).phase == CoolantPhase::kLiquid);
    REQUIRE(methane.at(tsat + 10.0, 2.0e6).phase == CoolantPhase::kVapor);
    // Above the critical temperature there is no liquid branch to boil out of.
    REQUIRE(methane.at(400.0, 2.0e6).phase == CoolantPhase::kSupercriticalTemperature);
    REQUIRE_THROWS_AS(methane.saturationTemperature(1.0e7), RangeError);
  }
  SECTION("enthalpy inversion round trips") {
    const auto st = methane.at(300.0, 1.5e7);
    REQUIRE(methane.temperatureFromEnthalpy(st.h, 1.5e7, 200.0) == Approx(300.0).epsilon(1e-6));
  }
}
