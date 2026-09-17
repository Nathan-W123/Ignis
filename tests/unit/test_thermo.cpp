// SPDX-License-Identifier: MIT
/// \file test_thermo.cpp
/// \brief NASA-7 polynomial verification.

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <algorithm>
#include <cmath>

#include "TestHelpers.hpp"
#include "ignis/thermo/GasMixture.hpp"

using namespace ignis;
using namespace ignis_test;
using Catch::Approx;

TEST_CASE("species database loads with the expected content", "[thermo]") {
  const auto& db = fullDatabase();
  REQUIRE(db.size() == 40);
  REQUIRE(db.numElements() == 5);
  for (const char* required : {"CH4", "O2", "H2", "H2O", "CO2", "CO", "OH", "O", "H", "NO", "N2"})
    REQUIRE(db.has(required));
  // Every species must carry provenance and a positive molar mass.
  for (const auto& s : db.species()) {
    REQUIRE_FALSE(s.source().empty());
    REQUIRE(s.molarMass() > 0.0);
    REQUIRE(s.tMin() <= 200.0);
    REQUIRE(s.tMax() >= 5000.0);
  }
  REQUIRE_FALSE(db.provenance().sources.empty());
}

TEST_CASE("molar masses follow from the elemental composition", "[thermo]") {
  const auto& db = fullDatabase();
  for (const auto& s : db.species()) {
    double m = 0.0;
    for (const auto& kv : s.composition()) m += kv.second * db.atomicWeight(kv.first);
    // The stored value is rounded to six decimals in g/mol.
    REQUIRE(relativeError(s.molarMass(), m) < 1.0e-6);
  }
}

TEST_CASE("cp equals dh/dT to the accuracy of a central difference", "[thermo]") {
  const auto& db = fullDatabase();
  for (const auto& s : db.species()) {
    for (double T : {250.0, 500.0, 900.0, 1500.0, 3000.0, 5000.0}) {
      const double dT = 1.0e-4 * T;
      // Stay inside one polynomial branch: the derivative identity holds on
      // each piece, and the join is tested separately.
      if (std::abs(T - s.tMid()) < 2.0 * dT) continue;
      const double numeric = (s.h(T + dT) - s.h(T - dT)) / (2.0 * dT);
      INFO(s.name() << " at " << T << " K");
      REQUIRE(relativeError(numeric, s.cp(T)) < 1.0e-8);
    }
  }
}

TEST_CASE("cp/T equals ds/dT to the accuracy of a central difference", "[thermo]") {
  const auto& db = fullDatabase();
  for (const auto& s : db.species()) {
    for (double T : {300.0, 1200.0, 3500.0}) {
      const double dT = 1.0e-4 * T;
      if (std::abs(T - s.tMid()) < 2.0 * dT) continue;
      const double numeric = (s.s0(T + dT) - s.s0(T - dT)) / (2.0 * dT);
      INFO(s.name() << " at " << T << " K");
      REQUIRE(relativeError(numeric, s.cp(T) / T) < 1.0e-8);
    }
  }
}

TEST_CASE("Gibbs energy is h - T s by construction", "[thermo]") {
  const auto& db = fullDatabase();
  for (const auto& s : db.species())
    for (double T : {200.0, 1000.0, 4000.0, 6000.0})
      REQUIRE(s.g0(T) == Approx(s.h(T) - T * s.s0(T)).epsilon(1e-14));
}

TEST_CASE("the two polynomial branches join continuously", "[thermo]") {
  const auto& db = fullDatabase();
  double worst_cp = 0.0, worst_h = 0.0, worst_s = 0.0;
  std::string worst_name;
  for (const auto& s : db.species()) {
    const double Tm = s.tMid();
    const double eps = 1.0e-7 * Tm;
    const double dcp = relativeError(s.cp(Tm - eps), s.cp(Tm + eps));
    const double dh = relativeError(s.h(Tm - eps), s.h(Tm + eps), 1.0e3);
    const double ds = relativeError(s.s0(Tm - eps), s.s0(Tm + eps));
    if (dcp > worst_cp) { worst_cp = dcp; worst_name = s.name(); }
    worst_h = std::max(worst_h, dh);
    worst_s = std::max(worst_s, ds);
  }
  INFO("largest cp jump on species " << worst_name);
  // NASA fits are constructed to be continuous at the join, but only to the
  // precision of the fitting procedure -- a few parts in 10^4 is normal.
  REQUIRE(worst_cp < 2.0e-4);
  REQUIRE(worst_h < 2.0e-5);
  REQUIRE(worst_s < 2.0e-5);
}

TEST_CASE("properties match an independent implementation of the same coefficients",
          "[thermo][verification]") {
  const ReferenceTable ref(referenceDir() + "/thermo_reference.csv");
  REQUIRE(ref.rows() > 300);
  const auto& db = fullDatabase();
  double worst_cp = 0.0, worst_h = 0.0, worst_s = 0.0;
  int checked = 0;
  for (std::size_t r = 0; r < ref.rows(); ++r) {
    const auto& name = ref.text("species", r);
    if (!db.has(name)) continue;
    const auto& s = db.at(name);
    const double T = ref.num("T", r);
    // Exactly at the branch join the two implementations may pick different
    // sides of the piecewise fit.  The size of that step is measured by the
    // continuity test above; here it would only mask the real comparison.
    if (std::abs(T - s.tMid()) < 1.0e-9) continue;
    worst_cp = std::max(worst_cp, relativeError(s.cp(T), ref.num("cp_nasa", r)));
    worst_h = std::max(worst_h, relativeError(s.h(T), ref.num("h_nasa", r), 1.0e3));
    worst_s = std::max(worst_s, relativeError(s.s0(T), ref.num("s_nasa", r)));
    ++checked;
  }
  INFO("checked " << checked << " species/temperature pairs against " << ref.path());
  // Same coefficients, different code: this must agree to round-off.
  // The reference file carries twelve significant digits.
  REQUIRE(checked > 300);
  REQUIRE(worst_cp < 1.0e-11);
  REQUIRE(worst_h < 1.0e-11);
  REQUIRE(worst_s < 1.0e-11);
}

TEST_CASE("properties are close to an independent fit of the same species",
          "[thermo][validation]") {
  // GRI-Mech 3.0 is a separate evaluation of the same thermochemistry.  The
  // spread between published data sets is a real modelling uncertainty, so this
  // test records its size rather than demanding agreement.
  const ReferenceTable ref(referenceDir() + "/thermo_reference.csv");
  const auto& db = fullDatabase();
  const std::vector<std::string> major = {"H2O", "CO2", "CO",  "H2", "O2", "OH",
                                          "O",   "H",   "N2",  "CH4", "NO"};
  double worst_major = 0.0, worst_any = 0.0;
  std::string where_major, where_any;
  int checked = 0;
  for (std::size_t r = 0; r < ref.rows(); ++r) {
    const auto& name = ref.text("species", r);
    const double cp_gri = ref.num("cp_gri", r);
    if (!db.has(name) || !std::isfinite(cp_gri)) continue;
    const double T = ref.num("T", r);
    if (T > 3500.0) continue;   // the GRI-Mech fits stop at 3500 K
    if (std::abs(T - db.at(name).tMid()) < 1.0e-9) continue;
    const double e = relativeError(db.at(name).cp(T), cp_gri);
    const std::string tag = name + " at " + std::to_string(T) + " K";
    if (e > worst_any) { worst_any = e; where_any = tag; }
    if (std::find(major.begin(), major.end(), name) != major.end() && e > worst_major) {
      worst_major = e;
      where_major = tag;
    }
    ++checked;
  }
  INFO("major-species disagreement with GRI-Mech 3.0: " << worst_major << " at " << where_major
       << "; worst over all species " << worst_any << " at " << where_any << " ("
       << checked << " points)");
  REQUIRE(checked > 200);
  // The species that dominate rocket combustion are well characterised and the
  // two data sets agree closely.  Minor radicals such as CH3O and CH2OH differ
  // by up to ten percent between published evaluations; that spread is a real
  // modelling uncertainty and is recorded in docs/validation.md.
  REQUIRE(worst_major < 0.01);
  REQUIRE(worst_any < 0.15);
}

TEST_CASE("mixture properties are consistent and dimensionally correct", "[thermo]") {
  const auto& db = rocketDatabase();
  const GasMixture mix(db);
  Eigen::VectorXd X = Eigen::VectorXd::Zero(static_cast<Eigen::Index>(db.size()));
  X(db.index("H2O")) = 0.5;
  X(db.index("CO2")) = 0.2;
  X(db.index("CO")) = 0.2;
  X(db.index("H2")) = 0.1;
  const Eigen::VectorXd n = mix.fromMoleFractions(X);

  SECTION("the per-kilogram basis is exact") {
    REQUIRE(std::abs(mix.massResidual(n)) < 1.0e-15);
  }
  SECTION("mole and mass fractions round trip") {
    const auto Y = n.cwiseProduct(db.molarMasses());
    const Eigen::VectorXd n2 = mix.fromMassFractions(Y);
    REQUIRE((n2 - n).cwiseAbs().maxCoeff() < 1.0e-12);
  }
  SECTION("state assembly satisfies the ideal-gas relations") {
    const auto st = mix.frozenState(n, 2500.0, 3.0e6);
    REQUIRE(st.rho == Approx(st.p / (st.R * st.T)).epsilon(1e-14));
    REQUIRE(st.R == Approx(constants::R_universal / st.M).epsilon(1e-14));
    REQUIRE(st.u == Approx(st.h - st.p / st.rho).epsilon(1e-12));
    REQUIRE(st.g == Approx(st.h - st.T * st.s).epsilon(1e-12));
    REQUIRE(st.cp_frozen - st.cv_frozen == Approx(st.R).epsilon(1e-12));
    REQUIRE(st.a == Approx(std::sqrt(st.gamma_frozen * st.R * st.T)).epsilon(1e-12));
  }
  SECTION("temperature inversions are exact") {
    const double T = 2137.0;
    const double h = mix.enthalpy(n, T);
    REQUIRE(mix.temperatureFromEnthalpy(n, h, 1000.0) == Approx(T).epsilon(1e-10));
    const double u = mix.internalEnergy(n, T);
    REQUIRE(mix.temperatureFromInternalEnergy(n, u, 4000.0) == Approx(T).epsilon(1e-10));
    const double s = mix.entropy(n, T, 2.0e6);
    REQUIRE(mix.temperatureFromEntropy(n, s, 2.0e6, 600.0) == Approx(T).epsilon(1e-10));
  }
}

TEST_CASE("out-of-range property evaluation is reported, not extrapolated", "[thermo]") {
  const auto& db = fullDatabase();
  const auto& h2o = db.at("H2O");
  REQUIRE_THROWS_AS(h2o.requireInRange(100.0), RangeError);
  REQUIRE_THROWS_AS(h2o.requireInRange(7000.0), RangeError);
  const GasMixture mix(db);
  REQUIRE_THROWS_AS(mix.requireInRange(50.0), RangeError);
  REQUIRE_THROWS_WITH(mix.requireInRange(50.0),
                      Catch::Matchers::ContainsSubstring("outside the fit range"));
}
