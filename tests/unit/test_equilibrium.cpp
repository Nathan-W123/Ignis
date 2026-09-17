// SPDX-License-Identifier: MIT
/// \file test_equilibrium.cpp
/// \brief Gibbs-minimisation verification.

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <cmath>
#include <random>

#include "TestHelpers.hpp"
#include "ignis/combustion/Propellant.hpp"
#include "ignis/equilibrium/Equilibrium.hpp"

using namespace ignis;
using namespace ignis_test;
using Catch::Approx;

namespace {

PropellantLibrary library() {
  return PropellantLibrary::loadYaml(sourceDir() + "/data/propellants/ignis_propellants.yaml");
}

}  // namespace

TEST_CASE("equilibrium conserves elements and mass and satisfies optimality",
          "[equilibrium]") {
  const auto& db = rocketDatabase();
  const EquilibriumSolver solver(db);
  const auto lib = library();
  const PropellantMixture base(lib.at("LOX"), lib.at("LCH4"), 3.4, 90.18, 111.66);

  for (double mr : {2.0, 3.0, 3.4, 4.0, 5.0}) {
    for (double p : {1.0e5, 1.0e6, 1.0e7, 2.0e7}) {
      const auto mix = base.withMixtureRatio(mr);
      const Eigen::VectorXd b = mix.elementMoles(db);
      const auto r = solver.hp(b, mix.enthalpy(db), p);
      INFO("O/F = " << mr << ", p = " << p << " Pa");

      REQUIRE(r.diagnostics.converged);
      REQUIRE(r.diagnostics.element_residual_rel < 1.0e-10);
      REQUIRE(std::abs(r.diagnostics.mass_residual) < 1.0e-12);
      REQUIRE(r.diagnostics.gibbs_residual < 1.0e-6);
      REQUIRE(r.diagnostics.state_residual < 1.0e-10);
      // Nonnegativity holds by construction of the log-space formulation.
      REQUIRE((r.state.n.array() >= 0.0).all());
      // Element balance measured directly rather than trusting the diagnostic.
      const Eigen::VectorXd achieved = db.elementMatrix() * r.state.n;
      REQUIRE((achieved - b).cwiseAbs().maxCoeff() / b.maxCoeff() < 1.0e-10);
      // Enthalpy closure.
      REQUIRE(r.state.h == Approx(mix.enthalpy(db)).epsilon(1e-10));
    }
  }
}

TEST_CASE("equilibrium satisfies the stationarity condition species by species",
          "[equilibrium]") {
  const auto& db = rocketDatabase();
  const EquilibriumSolver solver(db);
  const auto lib = library();
  const PropellantMixture mix(lib.at("LOX"), lib.at("LCH4"), 3.4, 90.18, 111.66);
  const auto r = solver.hp(mix.elementMoles(db), mix.enthalpy(db), 5.5e6);

  const double RT = constants::R_universal * r.state.T;
  const auto X = r.state.moleFractions();
  double worst = 0.0;
  for (std::size_t j = 0; j < db.size(); ++j) {
    const auto jj = static_cast<Eigen::Index>(j);
    if (X(jj) < 1.0e-12) continue;
    // mu_j / RT must equal sum_i a_ij pi_i.
    const double mu = db[j].g0(r.state.T) +
                      RT * std::log(X(jj) * r.state.p / constants::p_reference);
    double lhs = mu / RT;
    double rhs = 0.0;
    for (std::size_t i = 0; i < db.numElements(); ++i)
      rhs += db.elementMatrix()(static_cast<Eigen::Index>(i), jj) *
             r.pi(static_cast<Eigen::Index>(i));
    worst = std::max(worst, std::abs(lhs - rhs));
  }
  INFO("largest stationarity residual " << worst);
  REQUIRE(worst < 1.0e-8);
}

TEST_CASE("equilibrium is independent of the initial guess", "[equilibrium]") {
  const auto& db = rocketDatabase();
  const auto lib = library();
  const PropellantMixture mix(lib.at("LOX"), lib.at("LCH4"), 3.4, 90.18, 111.66);
  const Eigen::VectorXd b = mix.elementMoles(db);
  const double h0 = mix.enthalpy(db);

  EquilibriumOptions opts;
  const EquilibriumSolver solver(db, opts);
  const auto base = solver.hp(b, h0, 5.5e6, 3500.0);

  std::mt19937_64 rng(4242);
  std::uniform_real_distribution<double> uni(0.0, 1.0);
  for (int trial = 0; trial < 12; ++trial) {
    // Wildly different compositions and temperatures must all land on the same
    // answer: the Gibbs minimum is unique for a convex ideal-gas mixture.
    Eigen::VectorXd guess(static_cast<Eigen::Index>(db.size()));
    for (Eigen::Index j = 0; j < guess.size(); ++j) guess(j) = 1.0e-6 + 40.0 * uni(rng);
    const double T0 = 800.0 + 4000.0 * uni(rng);
    const auto r = solver.hp(b, h0, 5.5e6, T0, &guess);
    INFO("trial " << trial << " starting from T = " << T0 << " K");
    REQUIRE(r.state.T == Approx(base.state.T).epsilon(1e-8));
    REQUIRE(r.state.M == Approx(base.state.M).epsilon(1e-8));
    REQUIRE((r.state.moleFractions() - base.state.moleFractions()).cwiseAbs().maxCoeff() < 1e-8);
  }
}

TEST_CASE("equilibrium trends follow the physics", "[equilibrium]") {
  const auto& db = rocketDatabase();
  const EquilibriumSolver solver(db);
  const auto lib = library();
  const PropellantMixture base(lib.at("LOX"), lib.at("LCH4"), 3.4, 90.18, 111.66);

  SECTION("flame temperature peaks slightly rich of stoichiometric") {
    double best_T = 0.0, best_mr = 0.0;
    for (double mr = 2.0; mr <= 5.0; mr += 0.1) {
      const auto mix = base.withMixtureRatio(mr);
      const auto r = solver.hp(mix.elementMoles(db), mix.enthalpy(db), 5.5e6);
      if (r.state.T > best_T) { best_T = r.state.T; best_mr = mr; }
    }
    const double stoich = base.stoichiometricMixtureRatio();
    INFO("peak at O/F = " << best_mr << ", stoichiometric = " << stoich);
    REQUIRE(best_mr < stoich);           // dissociation moves the peak fuel-rich
    REQUIRE(best_mr > 0.75 * stoich);
  }
  SECTION("raising the pressure suppresses dissociation and raises the temperature") {
    double previous_T = 0.0, previous_X_radicals = 1.0;
    for (double p : {1.0e5, 1.0e6, 5.5e6, 2.0e7}) {
      const auto r = solver.hp(base.elementMoles(db), base.enthalpy(db), p);
      const auto X = r.state.moleFractions();
      const double radicals = X(db.index("OH")) + X(db.index("H")) + X(db.index("O"));
      INFO("p = " << p << " Pa, T = " << r.state.T << " K, radicals = " << radicals);
      REQUIRE(r.state.T > previous_T);
      REQUIRE(radicals < previous_X_radicals);
      previous_T = r.state.T;
      previous_X_radicals = radicals;
    }
  }
  SECTION("fuel-rich mixtures make CO and H2, lean mixtures make CO2, H2O and O2") {
    const auto rich = solver.hp(base.withMixtureRatio(2.5).elementMoles(db),
                                base.withMixtureRatio(2.5).enthalpy(db), 5.5e6);
    const auto lean = solver.hp(base.withMixtureRatio(5.0).elementMoles(db),
                                base.withMixtureRatio(5.0).enthalpy(db), 5.5e6);
    const auto Xr = rich.state.moleFractions();
    const auto Xl = lean.state.moleFractions();
    REQUIRE(Xr(db.index("CO")) > Xl(db.index("CO")));
    REQUIRE(Xr(db.index("H2")) > Xl(db.index("H2")));
    REQUIRE(Xl(db.index("CO2")) > Xr(db.index("CO2")));
    REQUIRE(Xl(db.index("O2")) > Xr(db.index("O2")));
  }
  SECTION("equilibrium cp exceeds frozen cp wherever the mixture dissociates") {
    const auto r = solver.hp(base.elementMoles(db), base.enthalpy(db), 5.5e6);
    REQUIRE(r.state.cp_eff > r.state.cp_frozen);
    REQUIRE(r.state.gamma_s < r.state.gamma_frozen);
    REQUIRE(r.state.dlnV_dlnT_p > 1.0);
    REQUIRE(r.state.dlnV_dlnp_T < -1.0);
  }
}

TEST_CASE("equilibrium derivative terms match finite differences", "[equilibrium]") {
  const auto& db = rocketDatabase();
  const EquilibriumSolver solver(db);
  const auto lib = library();
  const PropellantMixture mix(lib.at("LOX"), lib.at("LCH4"), 3.4, 90.18, 111.66);
  const Eigen::VectorXd b = mix.elementMoles(db);
  const double T = 3400.0, p = 5.5e6;
  const auto st = solver.tp(b, T, p).state;

  SECTION("(d ln V / d ln T)_p") {
    const double dT = 1.0e-4 * T;
    const auto hi = solver.tp(b, T + dT, p).state;
    const auto lo = solver.tp(b, T - dT, p).state;
    const double numeric = (std::log(1.0 / hi.rho) - std::log(1.0 / lo.rho)) /
                           (std::log(T + dT) - std::log(T - dT));
    REQUIRE(st.dlnV_dlnT_p == Approx(numeric).epsilon(2e-5));
  }
  SECTION("(d ln V / d ln p)_T") {
    const double dp = 1.0e-4 * p;
    const auto hi = solver.tp(b, T, p + dp).state;
    const auto lo = solver.tp(b, T, p - dp).state;
    const double numeric = (std::log(1.0 / hi.rho) - std::log(1.0 / lo.rho)) /
                           (std::log(p + dp) - std::log(p - dp));
    REQUIRE(st.dlnV_dlnp_T == Approx(numeric).epsilon(2e-5));
  }
  SECTION("equilibrium cp is dh/dT at constant pressure") {
    const double dT = 1.0e-3 * T;
    const auto hi = solver.tp(b, T + dT, p).state;
    const auto lo = solver.tp(b, T - dT, p).state;
    const double numeric = (hi.h - lo.h) / (2.0 * dT);
    REQUIRE(st.cp_eff == Approx(numeric).epsilon(5e-5));
  }
  SECTION("the isentropic exponent reproduces an isentropic pressure perturbation") {
    const double dp = 1.0e-4 * p;
    const auto hi = solver.sp(b, st.s, p + dp, st.T).state;
    const auto lo = solver.sp(b, st.s, p - dp, st.T).state;
    const double numeric = -(std::log(p + dp) - std::log(p - dp)) /
                           (std::log(1.0 / hi.rho) - std::log(1.0 / lo.rho));
    REQUIRE(st.gamma_s == Approx(numeric).epsilon(1e-4));
  }
}

TEST_CASE("the constant-entropy and constant-volume problems agree with the base solve",
          "[equilibrium]") {
  const auto& db = rocketDatabase();
  const EquilibriumSolver solver(db);
  const auto lib = library();
  const PropellantMixture mix(lib.at("LOX"), lib.at("LCH4"), 3.4, 90.18, 111.66);
  const Eigen::VectorXd b = mix.elementMoles(db);
  const auto base = solver.hp(b, mix.enthalpy(db), 5.5e6);

  SECTION("SP returns to the same point") {
    const auto sp = solver.sp(b, base.state.s, 5.5e6, 2000.0);
    REQUIRE(sp.state.T == Approx(base.state.T).epsilon(1e-9));
    REQUIRE(sp.state.s == Approx(base.state.s).epsilon(1e-10));
  }
  SECTION("UV returns to the same point") {
    const auto uv = solver.uv(b, base.state.u, base.state.v, 2500.0);
    REQUIRE(uv.state.T == Approx(base.state.T).epsilon(1e-7));
    REQUIRE(uv.state.p == Approx(base.state.p).epsilon(1e-7));
  }
  SECTION("TV returns to the same point") {
    const auto tv = solver.tv(b, base.state.T, base.state.v, 1.0e6);
    REQUIRE(tv.state.p == Approx(base.state.p).epsilon(1e-9));
  }
}

TEST_CASE("infeasible and non-convergent equilibrium problems are reported",
          "[equilibrium][errors]") {
  const auto& db = rocketDatabase();
  const EquilibriumSolver solver(db);
  Eigen::VectorXd b = Eigen::VectorXd::Zero(static_cast<Eigen::Index>(db.numElements()));

  SECTION("an all-zero element vector") {
    REQUIRE_THROWS_AS(solver.tp(b, 3000.0, 1.0e6), ConfigError);
  }
  SECTION("a negative element amount") {
    b(db.elementIndex("H")) = -1.0;
    REQUIRE_THROWS_AS(solver.tp(b, 3000.0, 1.0e6), ConfigError);
  }
  SECTION("a non-positive pressure") {
    b(db.elementIndex("H")) = 100.0;
    b(db.elementIndex("O")) = 50.0;
    REQUIRE_THROWS_AS(solver.tp(b, 3000.0, 0.0), ConfigError);
  }
  SECTION("a temperature outside the polynomial range") {
    b(db.elementIndex("H")) = 100.0;
    b(db.elementIndex("O")) = 50.0;
    REQUIRE_THROWS_AS(solver.tp(b, 50.0, 1.0e6), RangeError);
    REQUIRE_THROWS_AS(solver.tp(b, 9000.0, 1.0e6), RangeError);
  }
  SECTION("an element that no available species can carry") {
    const auto ho = hydrogenDatabase();
    const EquilibriumSolver ho_solver(ho);
    Eigen::VectorXd bb = Eigen::VectorXd::Zero(static_cast<Eigen::Index>(ho.numElements()));
    bb(ho.elementIndex("H")) = 100.0;
    bb(ho.elementIndex("O")) = 50.0;
    REQUIRE_NOTHROW(ho_solver.tp(bb, 3000.0, 1.0e6));
  }
  SECTION("an unreachable enthalpy") {
    b(db.elementIndex("H")) = 100.0;
    b(db.elementIndex("O")) = 50.0;
    b(db.elementIndex("C")) = 10.0;
    REQUIRE_THROWS_AS(solver.hp(b, 1.0e12, 1.0e6), ConvergenceError);
  }
}

TEST_CASE("solver diagnostics record real work", "[equilibrium]") {
  const auto& db = rocketDatabase();
  const EquilibriumSolver solver(db);
  const auto lib = library();
  const PropellantMixture mix(lib.at("LOX"), lib.at("LCH4"), 3.4, 90.18, 111.66);
  const auto r = solver.hp(mix.elementMoles(db), mix.enthalpy(db), 5.5e6);
  REQUIRE(r.diagnostics.iterations > 1);
  REQUIRE(r.diagnostics.iterations < 100);
  REQUIRE(r.diagnostics.history.size() == static_cast<std::size_t>(r.diagnostics.iterations));
  // The correction history must be decreasing overall.
  REQUIRE(r.diagnostics.history.back() < r.diagnostics.history.front());
  REQUIRE_FALSE(r.diagnostics.message.empty());
}
