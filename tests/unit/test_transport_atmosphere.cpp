// SPDX-License-Identifier: MIT
/// \file test_transport_atmosphere.cpp
/// \brief Chapman-Enskog transport and U.S. Standard Atmosphere verification.

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <sstream>

#include "TestHelpers.hpp"
#include "ignis/nozzle/Atmosphere.hpp"
#include "ignis/thermo/Transport.hpp"

using namespace ignis;
using namespace ignis_test;
using Catch::Approx;

TEST_CASE("collision integrals behave correctly", "[transport]") {
  // Omega(2,2)* falls monotonically with reduced temperature and tends towards
  // the rigid-sphere value of one at high T*.
  double previous = 1e30;
  for (double ts : {0.3, 0.5, 1.0, 2.0, 5.0, 10.0, 50.0, 100.0}) {
    const double o = omega22LJ(ts);
    REQUIRE(o > 0.0);
    REQUIRE(o < previous);
    previous = o;
  }
  REQUIRE(omega22LJ(100.0) < 1.0);
  REQUIRE(omega22LJ(100.0) > 0.5);
  // Omega(1,1)* is smaller than Omega(2,2)* over the useful range.
  for (double ts : {0.5, 1.0, 5.0, 20.0}) REQUIRE(omega11LJ(ts) < omega22LJ(ts));
}

TEST_CASE("pure-species viscosity matches an independent kinetic-theory implementation",
          "[transport][verification]") {
  const ReferenceTable ref(referenceDir() + "/transport_reference.csv");
  const auto& db = fullDatabase();
  const TransportModel model(db);
  int checked = 0;
  double worst_nonpolar = 0.0, worst_polar = 0.0;
  std::string where_nonpolar, where_polar;
  for (std::size_t r = 0; r < ref.rows(); ++r) {
    const auto& name = ref.text("mixture", r);
    if (!db.has(name)) continue;   // the rows named after a single species
    const double T = ref.num("T", r);
    const std::size_t j = static_cast<std::size_t>(db.index(name));
    const double mu = model.speciesViscosity(j, T);
    const double e = relativeError(mu, ref.num("viscosity", r));
    const std::string tag = name + " at " + std::to_string(T) + " K";
    if (db[j].transport().dipole > 0.0) {
      if (e > worst_polar) { worst_polar = e; where_polar = tag; }
    } else {
      if (e > worst_nonpolar) { worst_nonpolar = e; where_nonpolar = tag; }
    }
    ++checked;
  }
  INFO("non-polar worst " << worst_nonpolar << " at " << where_nonpolar
       << "; polar worst " << worst_polar << " at " << where_polar << " over " << checked
       << " points");
  REQUIRE(checked >= 30);
  // For non-polar species both codes evaluate the same Chapman-Enskog
  // expression with the same Lennard-Jones parameters, so they agree to well
  // under a percent.  For polar species Ignis uses Brokaw's closed-form
  // correction to the collision integral whereas Cantera interpolates the
  // Monchick-Mason Stockmayer table; the resulting few-percent difference for
  // water is a documented modelling choice, not an error.
  REQUIRE(worst_nonpolar < 0.01);
  REQUIRE(worst_polar < 0.05);
}

TEST_CASE("mixture viscosity matches an independent Wilke implementation",
          "[transport][verification]") {
  const ReferenceTable ref(referenceDir() + "/transport_reference.csv");
  const auto& db = fullDatabase();
  const TransportModel model(db);
  int checked = 0;
  double worst_mu = 0.0, worst_k = 0.0;
  for (std::size_t r = 0; r < ref.rows(); ++r) {
    const auto& name = ref.text("mixture", r);
    if (name.rfind("products", 0) != 0) continue;
    Eigen::VectorXd X = Eigen::VectorXd::Zero(static_cast<Eigen::Index>(db.size()));
    std::istringstream ss(ref.text("composition", r));
    std::string token;
    while (ss >> token) {
      const auto colon = token.find(':');
      const std::string sp = token.substr(0, colon);
      REQUIRE(db.has(sp));
      X(db.index(sp)) = std::atof(token.c_str() + colon + 1);
    }
    const double T = ref.num("T", r);
    const auto props = model.mixture(X, T, ref.num("cp_mass", r));
    worst_mu = std::max(worst_mu, relativeError(props.viscosity, ref.num("viscosity", r)));
    worst_k = std::max(worst_k, relativeError(props.conductivity, ref.num("conductivity", r)));
    REQUIRE(props.covered_mole_fraction == Approx(1.0).epsilon(1e-12));
    ++checked;
  }
  INFO("mixture viscosity worst " << worst_mu << ", conductivity worst " << worst_k
       << " over " << checked << " mixtures");
  REQUIRE(checked >= 8);
  // Same mixing rule (Wilke) and same Lennard-Jones data: viscosity agrees
  // closely.  Conductivity uses the modified Eucken correlation here against
  // Cantera's full kinetic-theory expression, so a larger, documented spread is
  // expected -- see docs/validation.md.
  REQUIRE(worst_mu < 0.05);
  REQUIRE(worst_k < 0.25);
}

TEST_CASE("transport reports missing data instead of guessing", "[transport][errors]") {
  const auto& db = fullDatabase();
  const TransportModel model(db);
  const int c2 = db.index("C2");
  REQUIRE(c2 >= 0);
  REQUIRE_FALSE(db[static_cast<std::size_t>(c2)].transport().valid);
  REQUIRE_THROWS_AS(model.speciesViscosity(static_cast<std::size_t>(c2), 1000.0), RangeError);

  SECTION("a mixture that is mostly C2 reports low coverage") {
    Eigen::VectorXd X = Eigen::VectorXd::Zero(static_cast<Eigen::Index>(db.size()));
    X(c2) = 0.9;
    X(db.index("H2O")) = 0.1;
    const auto props = model.mixture(X, 2000.0, 2000.0);
    REQUIRE(props.covered_mole_fraction == Approx(0.1).epsilon(1e-12));
  }
}

TEST_CASE("U.S. Standard Atmosphere 1976 reproduces its published values",
          "[atmosphere][verification]") {
  // Values from the standard itself.  Sea level is definitional; 86 km
  // geometric (84.852 km geopotential) is the top of the layer table.
  const auto sl = Atmosphere::at(0.0);
  REQUIRE(sl.pressure == Approx(101325.0).epsilon(1e-12));
  REQUIRE(sl.temperature == Approx(288.15).epsilon(1e-12));
  REQUIRE(sl.density == Approx(1.225).epsilon(2e-4));
  REQUIRE(sl.speed_of_sound == Approx(340.294).epsilon(1e-4));

  const auto top = Atmosphere::at(86000.0);
  REQUIRE(top.pressure == Approx(0.37338).epsilon(2e-4));
  REQUIRE(top.temperature == Approx(186.946).epsilon(1e-4));

  SECTION("pressure and temperature are continuous across the layer joins") {
    for (double z : {11019.0, 20063.0, 32162.0, 47350.0, 51413.0, 71802.0}) {
      const double eps = 1.0e-3;
      const auto a = Atmosphere::at(z - eps);
      const auto b = Atmosphere::at(z + eps);
      REQUIRE(relativeError(a.pressure, b.pressure) < 1.0e-6);
      REQUIRE(relativeError(a.temperature, b.temperature) < 1.0e-6);
    }
  }
  SECTION("pressure falls monotonically with altitude") {
    double previous = 1e30;
    for (double z = -4000.0; z <= 120000.0; z += 500.0) {
      const auto s = Atmosphere::at(z);
      REQUIRE(s.pressure < previous);
      REQUIRE(s.pressure > 0.0);
      previous = s.pressure;
    }
  }
  SECTION("the altitude inverse round trips") {
    for (double z : {0.0, 1000.0, 11000.0, 30000.0, 60000.0}) {
      const double p = Atmosphere::at(z).pressure;
      REQUIRE(Atmosphere::altitudeForPressure(p) == Approx(z).margin(0.5));
    }
  }
  SECTION("out-of-range requests are rejected") {
    REQUIRE_THROWS_AS(Atmosphere::at(-6000.0), RangeError);
    REQUIRE_THROWS_AS(Atmosphere::altitudeForPressure(-1.0), RangeError);
    REQUIRE_THROWS_AS(Atmosphere::altitudeForPressure(2.0e5), RangeError);
  }
  SECTION("extrapolation above the layer table is flagged") {
    REQUIRE_FALSE(Atmosphere::at(80000.0).extrapolated);
    REQUIRE(Atmosphere::at(100000.0).extrapolated);
  }
}
