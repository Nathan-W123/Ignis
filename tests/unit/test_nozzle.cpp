// SPDX-License-Identifier: MIT
/// \file test_nozzle.cpp
/// \brief Nozzle geometry and quasi-1D compressible-flow verification.
///
/// The compressible-flow tests build a *calorically perfect* gas out of a
/// NASA-7 polynomial whose higher coefficients are zero, so cp/R is exactly
/// a1 and gamma is exactly a1/(a1-1).  The general variable-property solver is
/// then run on it and compared against the closed-form isentropic and
/// normal-shock relations.  Nothing in the solver knows it is being handed a
/// constant-gamma gas, so this is a genuine check of the general path.

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <cmath>

#include "TestHelpers.hpp"
#include "ignis/combustion/Chamber.hpp"
#include "ignis/nozzle/NozzleFlow.hpp"
#include "ignis/nozzle/NozzleGeometry.hpp"

using namespace ignis;
using namespace ignis_test;
using Catch::Approx;

namespace {

/// gamma = 1.4, M = 28 g/mol calorically perfect gas as a one-species database.
SpeciesDatabase perfectGasDatabase(double gamma = 1.4, double molar_mass_g = 28.0) {
  const double cp_R = gamma / (gamma - 1.0);
  std::array<double, 7> a{};
  a[0] = cp_R;
  a[5] = 0.0;   // zero enthalpy of formation
  a[6] = 0.0;
  TransportData tr;
  tr.geometry = TransportData::Geometry::kLinear;
  tr.well_depth = 97.53;
  tr.diameter = 3.621;
  tr.valid = true;
  // A single fictitious element keeps the element matrix non-singular.
  const Species s("PG", {{"N", 2}}, molar_mass_g * 1e-3, 50.0, 1000.0, 20000.0, a, a, tr,
                  "synthetic calorically perfect gas for verification");
  // Build a database by writing a temporary YAML file would be circular, so the
  // subset machinery is used on a one-species vector instead.
  return SpeciesDatabase::fromSpecies({s}, {{"N", 14.007e-3}});
}

struct Perfect {
  double gamma, R;
  double areaRatio(double M) const {
    const double t = 1.0 + 0.5 * (gamma - 1.0) * M * M;
    return (1.0 / M) * std::pow(2.0 * t / (gamma + 1.0),
                                0.5 * (gamma + 1.0) / (gamma - 1.0));
  }
  double pressureRatio(double M) const {
    return std::pow(1.0 + 0.5 * (gamma - 1.0) * M * M, gamma / (gamma - 1.0));
  }
  double temperatureRatio(double M) const { return 1.0 + 0.5 * (gamma - 1.0) * M * M; }
  double cStar(double T0) const {
    return std::sqrt(R * T0 / gamma) * std::pow(0.5 * (gamma + 1.0),
                                                0.5 * (gamma + 1.0) / (gamma - 1.0));
  }
  double machBehindShock(double M1) const {
    return std::sqrt((1.0 + 0.5 * (gamma - 1.0) * M1 * M1) /
                     (gamma * M1 * M1 - 0.5 * (gamma - 1.0)));
  }
  double shockPressureRatio(double M1) const {
    return (2.0 * gamma * M1 * M1 - (gamma - 1.0)) / (gamma + 1.0);
  }
  double shockDensityRatio(double M1) const {
    return (gamma + 1.0) * M1 * M1 / ((gamma - 1.0) * M1 * M1 + 2.0);
  }
  double shockTemperatureRatio(double M1) const {
    return shockPressureRatio(M1) / shockDensityRatio(M1);
  }
};

}  // namespace

TEST_CASE("nozzle contour is valid, smooth and monotone", "[nozzle][geometry]") {
  NozzleGeometrySpec spec;
  spec.throat_radius = 0.05;
  spec.contraction_ratio = 3.0;
  spec.chamber_length = 0.18;
  spec.expansion_ratio = 25.0;
  spec.num_stations = 600;
  const auto g = NozzleGeometry::build(spec);

  SECTION("key dimensions follow from the inputs") {
    REQUIRE(g.throatArea() == Approx(constants::pi * 0.05 * 0.05).epsilon(1e-14));
    REQUIRE(g.expansionRatio() == Approx(25.0).epsilon(1e-10));
    REQUIRE(g.contractionRatio() == Approx(3.0).epsilon(1e-10));
    REQUIRE(g.throatPosition() > spec.chamber_length);
    REQUIRE(g.exitPosition() > g.throatPosition());
  }
  SECTION("segments join with C0 and C1 continuity") {
    double dr = 0.0, ds = 0.0;
    g.checkContinuity(dr, ds);
    INFO("radius jump " << dr << ", slope jump " << ds);
    REQUIRE(dr < 1.0e-12);
    REQUIRE(ds < 1.0e-5);
  }
  SECTION("the throat is the unique minimum and the area is monotone either side") {
    double r_min = 1e30;
    for (const auto& s : g.stations()) r_min = std::min(r_min, s.r);
    REQUIRE(r_min == Approx(g.throatRadius()).epsilon(1e-12));
    for (std::size_t i = 1; i < g.stations().size(); ++i) {
      const auto& a = g.stations()[i - 1];
      const auto& b = g.stations()[i];
      if (b.x <= g.throatPosition()) REQUIRE(b.r <= a.r + 1e-12);
      if (a.x >= g.throatPosition()) REQUIRE(b.r >= a.r - 1e-12);
    }
  }
  SECTION("analytic radius agrees with the sampled contour") {
    for (const auto& s : g.stations())
      REQUIRE(g.radius(s.x) == Approx(s.r).epsilon(1e-12));
  }
  SECTION("chamber volume matches an independent quadrature") {
    const int n = 20000;
    double v = 0.0;
    const double dx = g.throatPosition() / n;
    for (int i = 0; i < n; ++i) {
      const double r = g.radius((i + 0.5) * dx);
      v += constants::pi * r * r * dx;
    }
    REQUIRE(g.chamberVolume() == Approx(v).epsilon(2e-4));
    REQUIRE(g.characteristicLength() == Approx(g.chamberVolume() / g.throatArea()).epsilon(1e-14));
  }
}

TEST_CASE("invalid nozzle geometry is rejected with an explanatory message",
          "[nozzle][geometry][errors]") {
  NozzleGeometrySpec spec;
  spec.throat_radius = 0.05;
  spec.contraction_ratio = 3.0;
  spec.chamber_length = 0.15;
  spec.expansion_ratio = 25.0;

  SECTION("expansion ratio below one") {
    auto s = spec;
    s.expansion_ratio = 0.9;
    REQUIRE_THROWS_AS(NozzleGeometry::build(s), ConfigError);
  }
  SECTION("chamber narrower than the throat") {
    auto s = spec;
    s.contraction_ratio = 0.0;
    s.chamber_radius = 0.02;
    REQUIRE_THROWS_AS(NozzleGeometry::build(s), ConfigError);
  }
  SECTION("both throat radius and throat area") {
    auto s = spec;
    s.throat_area = 1.0e-3;
    REQUIRE_THROWS_WITH(NozzleGeometry::build(s),
                        Catch::Matchers::ContainsSubstring("not both"));
  }
  SECTION("converging section that cannot fit") {
    auto s = spec;
    s.throat_upstream_ratio = 9.0;
    s.chamber_fillet_ratio = 3.0;
    REQUIRE_THROWS_WITH(NozzleGeometry::build(s),
                        Catch::Matchers::ContainsSubstring("does not fit"));
  }
  SECTION("bell angles in the wrong order") {
    auto s = spec;
    s.bell_initial_angle = 8.0;
    s.bell_exit_angle = 30.0;
    REQUIRE_THROWS_AS(NozzleGeometry::build(s), ConfigError);
  }
  SECTION("too few stations") {
    auto s = spec;
    s.num_stations = 5;
    REQUIRE_THROWS_AS(NozzleGeometry::build(s), ConfigError);
  }
}

TEST_CASE("quasi-1D solver reproduces the analytic constant-gamma nozzle",
          "[nozzle][compressible][verification]") {
  const double gamma = 1.4, mw = 28.0e-3;
  const auto db = perfectGasDatabase(gamma, 28.0);
  const EquilibriumSolver solver(db);
  const GasMixture mix(db);
  const Perfect exact{gamma, constants::R_universal / mw};

  const double T0 = 3000.0, p0 = 5.0e6;
  Eigen::VectorXd n(1);
  n(0) = 1.0 / mw;
  const auto chamber = mix.frozenState(n, T0, p0);
  ChamberReference ref;
  ref.b = db.elementMatrix() * n;
  ref.stagnation = chamber;
  const NozzleFlow flow(solver, CompositionModel::kFrozen, ref);

  SECTION("the sonic point sits where the analytic relations put it") {
    const auto th = flow.throat();
    REQUIRE(th.mach == Approx(1.0).epsilon(1e-7));
    REQUIRE(th.gas.p == Approx(p0 / exact.pressureRatio(1.0)).epsilon(1e-8));
    REQUIRE(th.gas.T == Approx(T0 / exact.temperatureRatio(1.0)).epsilon(1e-8));
    REQUIRE(flow.cStarIdeal() == Approx(exact.cStar(T0)).epsilon(1e-8));
  }
  SECTION("area-Mach relation on both branches") {
    for (double M : {0.1, 0.3, 0.6, 0.9, 1.5, 2.0, 3.0, 4.0, 5.0}) {
      const double eps = exact.areaRatio(M);
      const auto st = flow.atAreaRatio(eps, M > 1.0);
      INFO("target Mach " << M << " at A/At = " << eps);
      REQUIRE(st.mach == Approx(M).epsilon(2e-6));
      REQUIRE(st.gas.p == Approx(p0 / exact.pressureRatio(M)).epsilon(2e-6));
      REQUIRE(st.gas.T == Approx(T0 / exact.temperatureRatio(M)).epsilon(2e-6));
    }
  }
  SECTION("mass and stagnation enthalpy are conserved") {
    for (double eps : {1.2, 2.0, 5.0, 15.0, 40.0}) {
      for (bool sup : {true, false}) {
        const auto st = flow.atAreaRatio(eps, sup);
        REQUIRE(flow.massFlowResidual(st) < 1.0e-9);
        REQUIRE(flow.energyResidual(st) < 1.0e-12);
      }
    }
  }
  SECTION("gamma really is constant, so the isentropic exponent equals cp/cv") {
    for (double eps : {1.0, 3.0, 20.0}) {
      const auto st = flow.atAreaRatio(eps, true);
      REQUIRE(st.gas.gamma_frozen == Approx(gamma).epsilon(1e-12));
      REQUIRE(st.gas.gamma_s == Approx(gamma).epsilon(1e-12));
    }
  }
  SECTION("normal-shock jump matches the Rankine-Hugoniot relations") {
    for (double M1 : {1.5, 2.0, 3.0, 4.0}) {
      const auto up = flow.atAreaRatio(exact.areaRatio(M1), true);
      const auto sh = solveNormalShock(mix, up);
      INFO("upstream Mach " << M1);
      REQUIRE(sh.mach_downstream == Approx(exact.machBehindShock(M1)).epsilon(5e-6));
      REQUIRE(sh.downstream.p / up.gas.p == Approx(exact.shockPressureRatio(M1)).epsilon(5e-6));
      REQUIRE(sh.downstream.T / up.gas.T ==
              Approx(exact.shockTemperatureRatio(M1)).epsilon(5e-6));
      REQUIRE(sh.downstream.rho / up.gas.rho ==
              Approx(exact.shockDensityRatio(M1)).epsilon(5e-6));
      // Mass, momentum and energy across the jump.
      REQUIRE(sh.downstream.rho * sh.u_downstream ==
              Approx(up.gas.rho * up.u).epsilon(1e-9));
      REQUIRE(sh.downstream.p + sh.downstream.rho * sh.u_downstream * sh.u_downstream ==
              Approx(up.gas.p + up.gas.rho * up.u * up.u).epsilon(1e-9));
      REQUIRE(sh.downstream.h + 0.5 * sh.u_downstream * sh.u_downstream ==
              Approx(up.gas.h + 0.5 * up.u * up.u).epsilon(1e-9));
      REQUIRE(sh.entropy_rise > 0.0);
      REQUIRE(sh.stagnation_pressure_ratio < 1.0);
    }
  }
  SECTION("a subsonic upstream state cannot support a shock") {
    const auto sub = flow.atAreaRatio(2.0, false);
    REQUIRE_THROWS_AS(solveNormalShock(mix, sub), InfeasibleError);
  }
}

TEST_CASE("thrust reduces to the momentum term when the exit is matched",
          "[nozzle][compressible]") {
  const double gamma = 1.25, mw = 22.0e-3;
  const auto db = perfectGasDatabase(gamma, 22.0);
  const EquilibriumSolver solver(db);
  const GasMixture mix(db);
  Eigen::VectorXd n(1);
  n(0) = 1.0 / mw;
  ChamberReference ref;
  ref.b = db.elementMatrix() * n;
  ref.stagnation = mix.frozenState(n, 3200.0, 6.0e6);
  const NozzleFlow flow(solver, CompositionModel::kFrozen, ref);

  NozzleGeometrySpec gs;
  gs.throat_radius = 0.06;
  gs.contraction_ratio = 3.0;
  gs.chamber_length = 0.2;
  gs.expansion_ratio = 12.0;
  const auto geom = NozzleGeometry::build(gs);

  NozzlePerformanceOptions opts;
  opts.auto_divergence = false;
  const auto exit_state = flow.atAreaRatio(geom.expansionRatio(), true);
  const auto perf = evaluateNozzle(flow, geom, exit_state.gas.p, opts);

  REQUIRE(perf.regime == ExpansionRegime::kIdeallyExpanded);
  REQUIRE(std::abs(perf.thrust_pressure) < 1.0e-6 * perf.thrust_momentum);
  REQUIRE(perf.thrust == Approx(perf.mdot * perf.u_exit).epsilon(1e-10));
  REQUIRE(perf.isp == Approx(perf.u_exit / constants::g0).epsilon(1e-10));
  REQUIRE(perf.c_effective == Approx(perf.u_exit).epsilon(1e-10));
  // Thrust coefficient definition.
  REQUIRE(perf.cf == Approx(perf.thrust / (perf.p_chamber * perf.throat_area)).epsilon(1e-12));
  // c* Cf = Isp g0.
  REQUIRE(perf.c_star * perf.cf == Approx(perf.isp * constants::g0).epsilon(1e-10));
}

TEST_CASE("expansion regime classification follows the exit pressure", "[nozzle][compressible]") {
  const auto& db = rocketDatabase();
  const EquilibriumSolver solver(db);
  const auto lib = PropellantLibrary::loadYaml(sourceDir() +
                                               "/data/propellants/ignis_propellants.yaml");
  const PropellantMixture mix(lib.at("LOX"), lib.at("LCH4"), 3.4, 90.18, 111.66);
  const CombustionChamber chamber(solver, CompositionModel::kEquilibrium);
  const auto ch = chamber.solve(mix, 5.5e6, 1.0);
  const auto flow = chamber.makeFlow(ch);

  NozzleGeometrySpec gs;
  gs.throat_radius = 0.07;
  gs.contraction_ratio = 2.8;
  gs.chamber_length = 0.22;
  gs.expansion_ratio = 20.0;
  const auto geom = NozzleGeometry::build(gs);
  NozzlePerformanceOptions opts;
  opts.auto_divergence = false;

  const double pe = flow.atAreaRatio(20.0, true).gas.p;
  REQUIRE(evaluateNozzle(flow, geom, 0.5 * pe, opts).regime == ExpansionRegime::kUnderExpanded);
  REQUIRE(evaluateNozzle(flow, geom, pe, opts).regime == ExpansionRegime::kIdeallyExpanded);
  REQUIRE(evaluateNozzle(flow, geom, 2.0 * pe, opts).regime == ExpansionRegime::kOverExpanded);

  SECTION("vacuum thrust exceeds sea-level thrust by exactly the exit-area term") {
    const auto vac = evaluateNozzle(flow, geom, 0.0, opts);
    const auto sl = evaluateNozzle(flow, geom, 101325.0, opts);
    REQUIRE(vac.thrust - sl.thrust == Approx(101325.0 * geom.exitArea()).epsilon(1e-9));
  }
  SECTION("separation is only reported when over-expanded") {
    REQUIRE_FALSE(evaluateNozzle(flow, geom, 0.5 * pe, opts).separation_predicted);
    const auto deep = evaluateNozzle(flow, geom, 6.0 * pe, opts);
    REQUIRE(deep.separation_predicted);
    REQUIRE(deep.separation_area_ratio > 1.0);
    REQUIRE(deep.separation_area_ratio <= geom.expansionRatio());
  }
}

TEST_CASE("an expansion beyond the property range is reported, not extrapolated",
          "[nozzle][compressible][errors]") {
  const auto& db = rocketDatabase();
  const EquilibriumSolver solver(db);
  const auto lib = PropellantLibrary::loadYaml(sourceDir() +
                                               "/data/propellants/ignis_propellants.yaml");
  const PropellantMixture mix(lib.at("LOX"), lib.at("LCH4"), 3.4, 90.18, 111.66);
  const CombustionChamber chamber(solver, CompositionModel::kEquilibrium);
  const auto ch = chamber.solve(mix, 5.5e6, 1.0);
  const auto flow = chamber.makeFlow(ch);
  // The reachable expansion is bounded by the 200 K floor of the NASA fits.
  const double eps_max = flow.maxAreaRatio();
  INFO("largest reachable area ratio " << eps_max << " at p = " << flow.pressureFloor() << " Pa");
  REQUIRE(eps_max > 50.0);
  REQUIRE_NOTHROW(flow.atAreaRatio(0.5 * eps_max, true));
  REQUIRE_THROWS_AS(flow.atAreaRatio(2.0 * eps_max, true), InfeasibleError);
  REQUIRE_THROWS_WITH(flow.atAreaRatio(2.0 * eps_max, true),
                      Catch::Matchers::ContainsSubstring("exceeds the largest expansion"));
}
