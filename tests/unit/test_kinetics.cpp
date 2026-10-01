// SPDX-License-Identifier: MIT
/// \file test_kinetics.cpp
/// \brief The reaction mechanism and the finite-rate nozzle: identities and
///        limits.

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <cmath>

#include "TestHelpers.hpp"
#include "ignis/combustion/Chamber.hpp"
#include "ignis/core/Constants.hpp"
#include "ignis/core/Exceptions.hpp"
#include "ignis/kinetics/NozzleKinetics.hpp"

using namespace ignis;
using namespace ignis_test;
using Catch::Approx;

namespace {

std::string mechanismPath() { return sourceDir() + "/data/kinetics/gri30_nozzle.yaml"; }

const PropellantLibrary& library() {
  static const PropellantLibrary lib =
      PropellantLibrary::loadYaml(sourceDir() + "/data/propellants/ignis_propellants.yaml");
  return lib;
}

/// A LOX/methane nozzle: chamber, equilibrium flow and a bell contour.
struct MethaneNozzle {
  const SpeciesDatabase& db = rocketDatabase();
  EquilibriumSolver solver{db};
  PropellantMixture mix{library().at("LOX"), library().at("LCH4"), 3.4, 90.18, 111.66};
  CombustionChamber chamber{solver, CompositionModel::kEquilibrium};
  ChamberResult ch;
  NozzleGeometry geom;
  NozzleFlow flow;
  Mechanism mech;
  explicit MethaneNozzle(double pc = 5.5e6, double eps = 20.0)
      : ch(chamber.solve(mix, pc, 1.0)),
        geom([eps] {
          NozzleGeometrySpec s;
          s.throat_radius = 0.070;
          s.contraction_ratio = 2.8;
          s.chamber_length = 0.22;
          s.expansion_ratio = eps;
          s.num_stations = 300;
          return NozzleGeometry::build(s);
        }()),
        flow(chamber.makeFlow(ch)),
        mech(Mechanism::loadYaml(mechanismPath(), db)) {}
  KineticNozzleResult run(double multiplier = 1.0, double start_mach = 1.10) const {
    KineticNozzleOptions o;
    o.rate_multiplier = multiplier;
    o.start_frozen_mach = start_mach;
    return integrateKineticNozzle(flow, geom, mech, o);
  }
};

}  // namespace

TEST_CASE("the nozzle mechanism maps onto each propellant's species", "[kinetics][mechanism]") {
  const auto cho = Mechanism::loadYaml(mechanismPath(), rocketDatabase());
  // GRI-Mech's H + O2 + N2 and H + O2 + AR have no collider in a C/H/O gas.
  CHECK(cho.reactions().size() == 39);
  CHECK(cho.dropped().size() == 2);
  const auto ho = Mechanism::loadYaml(mechanismPath(), hydrogenDatabase());
  // No carbon: CO, CO2 and HCO reactions cannot occur either.
  CHECK(ho.reactions().size() < cho.reactions().size());
  for (const auto& r : ho.reactions())
    CHECK(r.equation.find("CO") == std::string::npos);
  CHECK(ho.reactions().size() > 20);
}

TEST_CASE("rate constants are GRI-Mech 3.0's in mol-based SI units", "[kinetics][mechanism]") {
  const auto m = Mechanism::loadYaml(mechanismPath(), rocketDatabase());
  auto find = [&](const std::string& eq) -> std::size_t {
    for (std::size_t i = 0; i < m.reactions().size(); ++i)
      if (m.reactions()[i].equation == eq) return i;
    FAIL("no reaction " << eq);
    return 0;
  };
  // Straight from the GRI-Mech 3.0 listing (cm, mol, s, cal/mol):
  //   H+O2<=>O+OH       2.650E+16  -0.6707  17041.00
  //   2O+M<=>O2+M       1.200E+17  -1.000       0.00
  const double cal = 4.184;
  const double T = 2000.0;
  double kf = 0.0, kr = 0.0;
  m.rateConstants(find("H + O2 <=> O + OH"), T, 1.0, kf, kr);
  const double gri = 2.650e16 * std::pow(T, -0.6707) *
                     std::exp(-17041.0 * cal / (constants::R_universal * T)) * 1e-6;
  CHECK(kf == Approx(gri).epsilon(1e-9));
  m.rateConstants(find("2 O + M <=> O2 + M"), T, 1.0, kf, kr);
  CHECK(kf == Approx(1.200e17 / T * 1e-12).epsilon(1e-9));
}

TEST_CASE("at chemical equilibrium every reaction is balanced", "[kinetics][mechanism]") {
  // Reverse rates come from the same species data as the equilibrium solver,
  // so at the solver's equilibrium each reaction's forward and reverse rates
  // cancel.  This is detailed balance, checked reaction by reaction.
  const auto& db = rocketDatabase();
  const EquilibriumSolver solver(db);
  const auto m = Mechanism::loadYaml(mechanismPath(), db);
  const PropellantMixture mix(library().at("LOX"), library().at("LCH4"), 3.4, 90.18, 111.66);
  for (const double T : {1500.0, 2500.0, 3500.0}) {
    const auto st = solver.tp(mix.elementMoles(db), T, 2.0e6).state;
    const Eigen::VectorXd C = st.rho * st.n;
    Eigen::VectorXd wdot;
    m.productionRates(T, C, wdot);
    double total = 0.0;
    for (const auto& r : m.reactions()) {
      double kf = 0.0, kr = 0.0;
      double M = 1.0;
      if (r.type != ReactionType::kElementary) {
        M = r.default_efficiency * C.sum();
        for (const auto& [j, e] : r.efficiencies) M += (e - r.default_efficiency) * C(j);
      }
      m.rateConstants(static_cast<std::size_t>(&r - &m.reactions()[0]), T, M, kf, kr);
      double fwd = kf, rev = kr;
      for (std::size_t i = 0; i < r.reactants.size(); ++i)
        fwd *= std::pow(C(r.reactants[i]), r.reactant_nu[i]);
      for (std::size_t i = 0; i < r.products.size(); ++i)
        rev *= std::pow(C(r.products[i]), r.product_nu[i]);
      INFO(r.equation << " at " << T << " K");
      if (fwd > 0.0) CHECK(std::abs(fwd - rev) / fwd < 1e-7);
      total += fwd;
    }
    // ... and so the net production of every species vanishes.
    INFO("T = " << T);
    CHECK(wdot.cwiseAbs().maxCoeff() < 1e-7 * total);
  }
}

TEST_CASE("the finite-rate nozzle lies between frozen and shifting", "[kinetics][nozzle]") {
  const MethaneNozzle n;
  const auto k = n.run();
  INFO(k.summary());
  CHECK(k.start_frozen_mach == Approx(1.10).epsilon(1e-6));
  CHECK(k.isp_vacuum_frozen < k.isp_vacuum);
  CHECK(k.isp_vacuum < k.isp_vacuum_shifting);
  CHECK(k.recovered_fraction > 0.0);
  CHECK(k.recovered_fraction < 1.0);
  CHECK(k.kinetic_efficiency == Approx(k.isp_vacuum / k.isp_vacuum_shifting).epsilon(1e-14));

  SECTION("mass, energy and every element are conserved") {
    CHECK(k.energy_residual < 1e-10);
    CHECK(k.element_residual < 1e-10);
    const double r = n.geom.exitRadius();
    CHECK(k.rho_exit * k.u_exit * constants::pi * r * r == Approx(k.mass_flow).epsilon(1e-12));
  }
}

TEST_CASE("the finite-rate nozzle tends to its two limits", "[kinetics][nozzle]") {
  // Verification tolerances follow from the integration tolerance (1e-6
  // relative): the march carries errors of that order, not smaller.
  const MethaneNozzle n;
  SECTION("with no reactions it is a frozen expansion from its start") {
    const auto k = n.run(0.0);
    INFO("frozen march " << k.isp_vacuum << " s, frozen isentrope " << k.isp_vacuum_frozen << " s");
    CHECK(k.isp_vacuum == Approx(k.isp_vacuum_frozen).epsilon(1e-5));
    CHECK(k.element_residual < 1e-12);
  }
  SECTION("faster reactions close the gap to shifting equilibrium") {
    double previous = 1e30;
    for (const double m : {1.0, 10.0, 100.0, 1000.0}) {
      const auto k = n.run(m);
      const double gap = k.isp_vacuum_shifting - k.isp_vacuum;
      INFO("rates x" << m << ": " << gap << " s short of shifting");
      CHECK(gap < previous / 3.0);
      previous = gap;
    }
    // At a thousand times GRI's rates the gas tracks equilibrium to 1e-4.
    CHECK(std::abs(previous) < 1e-4 * n.run(1.0).isp_vacuum_shifting);
  }
  SECTION("where the march starts barely matters") {
    const auto early = n.run(1.0, 1.05);
    const auto late = n.run(1.0, 1.20);
    INFO("start at frozen Mach 1.05: " << early.isp_vacuum << " s; at 1.20: " << late.isp_vacuum
         << " s");
    CHECK(early.isp_vacuum == Approx(late.isp_vacuum).epsilon(1e-5));
  }
  SECTION("the step control converges") {
    KineticNozzleOptions tight;
    tight.relative_tolerance = 1e-8;
    const auto a = n.run(1.0);
    const auto b = integrateKineticNozzle(n.flow, n.geom, n.mech, tight);
    CHECK(a.isp_vacuum == Approx(b.isp_vacuum).epsilon(1e-6));
    CHECK(b.steps > a.steps);
  }
}

TEST_CASE("a hydrogen nozzle recombines too without the carbon reactions", "[kinetics][nozzle]") {
  const auto& db = hydrogenDatabase();
  const EquilibriumSolver solver(db);
  const PropellantMixture mix(library().at("LOX"), library().at("LH2"), 5.5, 90.18, 20.27);
  const CombustionChamber chamber(solver, CompositionModel::kEquilibrium);
  const auto ch = chamber.solve(mix, 5.5e6, 1.0);
  NozzleGeometrySpec s;
  s.throat_radius = 0.060;
  s.contraction_ratio = 3.2;
  s.chamber_length = 0.25;
  s.expansion_ratio = 60.0;
  s.num_stations = 300;
  const auto geom = NozzleGeometry::build(s);
  const auto flow = chamber.makeFlow(ch);
  const auto mech = Mechanism::loadYaml(mechanismPath(), db);
  const auto k = integrateKineticNozzle(flow, geom, mech);
  INFO(k.summary());
  CHECK(k.isp_vacuum_frozen < k.isp_vacuum);
  CHECK(k.isp_vacuum < k.isp_vacuum_shifting);
  CHECK(k.energy_residual < 1e-10);
  CHECK(k.element_residual < 1e-10);
  for (const auto& sp : k.species) CHECK(sp.find('C') == std::string::npos);
}

TEST_CASE("finite-rate inputs are checked", "[kinetics][errors]") {
  const MethaneNozzle n;
  KineticNozzleOptions o;
  o.start_frozen_mach = 0.9;
  CHECK_THROWS_AS(integrateKineticNozzle(n.flow, n.geom, n.mech, o), ConfigError);
  o = {};
  o.rate_multiplier = -1.0;
  CHECK_THROWS_AS(integrateKineticNozzle(n.flow, n.geom, n.mech, o), ConfigError);
  const auto other = Mechanism::loadYaml(mechanismPath(), hydrogenDatabase());
  CHECK_THROWS_AS(integrateKineticNozzle(n.flow, n.geom, other), ConfigError);
}
