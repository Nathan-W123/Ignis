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
#include "ignis/nozzle/Atmosphere.hpp"
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

TEST_CASE("the ambient-pressure family follows the exit state exactly",
          "[nozzle][compressible][verification]") {
  // Thrust is linear in ambient pressure for a full-flowing nozzle:
  //     F(p_a) = F_vac - p_a A_e
  // so the sea-level and ascent-averaged quantities must be reproducible from
  // the vacuum thrust and the exit area alone.  If they ever stop being, the
  // shortcut in evaluateNozzle is wrong.
  const double gamma = 1.2, mw = 22.0e-3;
  const auto db = perfectGasDatabase(gamma, 22.0);
  const EquilibriumSolver solver(db);
  const GasMixture mix(db);
  Eigen::VectorXd n(1);
  n(0) = 1.0 / mw;
  const auto chamber = mix.frozenState(n, 3000.0, 5.0e6);
  ChamberReference ref;
  ref.b = db.elementMatrix() * n;
  ref.stagnation = chamber;
  const NozzleFlow flow(solver, CompositionModel::kFrozen, ref);

  NozzleGeometrySpec gs;
  gs.throat_radius = 0.05;
  gs.contraction_ratio = 3.0;
  gs.chamber_length = 0.2;
  gs.expansion_ratio = 20.0;
  const auto geom = NozzleGeometry::build(gs);

  NozzlePerformanceOptions opts;
  opts.auto_divergence = false;
  opts.separation = SeparationCriterion::kNone;
  opts.ascent_altitudes = {0.0, 10000.0, 30000.0};
  opts.ascent_weights = {2.0, 1.0, 1.0};   // deliberately not normalised

  const auto perf = evaluateNozzle(flow, geom, constants::atm, opts);
  const double f_vac = perf.isp_vacuum * perf.mdot * constants::g0;

  SECTION("sea-level thrust is the vacuum thrust less the exit-plane term") {
    REQUIRE(perf.thrust_sea_level ==
            Approx(f_vac - constants::atm * perf.exit_area).epsilon(1e-12));
    REQUIRE(perf.isp_sea_level ==
            Approx(perf.thrust_sea_level / (perf.mdot * constants::g0)).epsilon(1e-12));
  }
  SECTION("the ascent average equals the average over the declared altitudes") {
    double w = 0.0, wf = 0.0, wp = 0.0;
    for (std::size_t i = 0; i < opts.ascent_altitudes.size(); ++i) {
      const double p_a = Atmosphere::at(opts.ascent_altitudes[i]).pressure;
      w += opts.ascent_weights[i];
      wf += opts.ascent_weights[i] * (f_vac - p_a * perf.exit_area);
      wp += opts.ascent_weights[i] * p_a;
    }
    REQUIRE(perf.ascent_mean_ambient == Approx(wp / w).epsilon(1e-12));
    REQUIRE(perf.isp_ascent ==
            Approx((wf / w) / (perf.mdot * constants::g0)).epsilon(1e-12));
    // Linearity means the weighted mean is the value at the mean pressure.
    REQUIRE(perf.isp_ascent ==
            Approx((f_vac - (wp / w) * perf.exit_area) /
                   (perf.mdot * constants::g0)).epsilon(1e-12));
  }
  SECTION("the ascent Isp is bracketed by the sea-level and vacuum values") {
    REQUIRE(perf.isp_ascent > perf.isp_sea_level);
    REQUIRE(perf.isp_ascent < perf.isp_vacuum);
  }
  SECTION("a malformed ascent profile is rejected") {
    auto bad = opts;
    bad.ascent_weights = {1.0, 1.0};   // one short
    REQUIRE_THROWS_AS(evaluateNozzle(flow, geom, constants::atm, bad), ConfigError);
    bad = opts;
    bad.ascent_weights = {0.0, 0.0, 0.0};
    REQUIRE_THROWS_AS(evaluateNozzle(flow, geom, constants::atm, bad), ConfigError);
  }
  SECTION("the separation margin is off when the criterion is off") {
    REQUIRE(perf.separation_margin == Approx(NozzlePerformance::kNoSeparationRisk));
  }
  SECTION("the separation margin changes sign where the criterion does") {
    auto on = opts;
    on.separation = SeparationCriterion::kSummerfield;
    // Summerfield separates when the wall pressure drops below 0.4 p_ambient,
    // so raising the ambient pressure past p_e / 0.4 must flip the margin.
    const auto attached = evaluateNozzle(flow, geom, 1.0e3, on);
    REQUIRE(attached.separation_margin > 0.0);
    const double p_flip = attached.p_exit / 0.4;
    const auto separated = evaluateNozzle(flow, geom, 2.0 * p_flip, on);
    REQUIRE(separated.separation_margin < 0.0);
    const auto marginal = evaluateNozzle(flow, geom, p_flip, on);
    REQUIRE(marginal.separation_margin == Approx(0.0).margin(1e-9));
  }
}


// ===========================================================================
// Tabulated contours
// ===========================================================================

namespace {

/// Sample an analytic geometry into the (x, r) table a file would carry.
std::vector<ignis::ContourSample> tabulate(const ignis::NozzleGeometry& g) {
  std::vector<ignis::ContourSample> out;
  out.reserve(g.stations().size());
  for (const auto& p : g.stations()) out.push_back({p.x, p.r});
  return out;
}

ignis::NozzleGeometrySpec m1Spec() {
  ignis::NozzleGeometrySpec spec;
  spec.throat_radius = 0.070;
  spec.contraction_ratio = 2.8;
  spec.chamber_length = 0.22;
  spec.converging_half_angle = 30.0;
  spec.throat_upstream_ratio = 1.5;
  spec.throat_downstream_ratio = 0.382;
  spec.expansion_ratio = 20.0;
  spec.divergent = ignis::DivergentType::kBell;
  spec.bell_length_fraction = 0.8;
  spec.bell_initial_angle = 33.0;
  spec.bell_exit_angle = 8.0;
  spec.num_stations = 400;
  return spec;
}

}  // namespace

TEST_CASE("a tabulated contour reproduces the analytic one it was sampled from",
          "[nozzle][geometry][contour]") {
  // The strongest check available without external hardware: take a contour
  // Ignis built analytically, throw the parameterisation away, hand back only
  // the (x, r) table a CAD revolve would give, and see whether the geometry
  // that comes back is the same engine.
  const ignis::NozzleGeometry analytic = ignis::NozzleGeometry::build(m1Spec());
  const ignis::NozzleGeometry table =
      ignis::NozzleGeometry::fromContour(tabulate(analytic), m1Spec());

  SECTION("the measured throat and exit match the analytic ones") {
    CHECK(table.isTabulated());
    CHECK_FALSE(analytic.isTabulated());
    // The throat is the minimum of the table, so it can only be as good as the
    // sampling.  400 stations over this contour put it within a tenth of a
    // percent, and the error is one-sided: a sampled minimum is never below
    // the true one.
    CHECK(table.throatRadius() >= analytic.throatRadius());
    CHECK(table.throatRadius() ==
          Approx(analytic.throatRadius()).epsilon(1e-3));
    CHECK(table.exitRadius() == Approx(analytic.exitRadius()).epsilon(1e-12));
    CHECK(table.exitPosition() == Approx(analytic.exitPosition()).epsilon(1e-12));
    CHECK(table.expansionRatio() == Approx(analytic.expansionRatio()).epsilon(2e-3));
    CHECK(table.contractionRatio() == Approx(analytic.contractionRatio()).epsilon(2e-3));
  }

  SECTION("the chamber volume and L* survive the round trip") {
    // Both are integrals over the same stations, so they should agree far
    // better than the throat does.
    CHECK(table.chamberVolume() ==
          Approx(analytic.chamberVolume()).epsilon(1e-6));
    CHECK(table.characteristicLength() ==
          Approx(analytic.characteristicLength()).epsilon(3e-3));
  }

  SECTION("the throat curvature is recovered by the circle fit, not assumed") {
    // The analytic throat is two arcs: Ru = 1.5 rt upstream, Rd = 0.382 rt
    // downstream, and Bartz is handed their mean.  A circle fitted across both
    // should land between them rather than on either.
    const double ru = 1.5 * analytic.throatRadius();
    const double rd = 0.382 * analytic.throatRadius();
    INFO("fit gave " << table.bartzCurvatureRadius() << " m; Rd = " << rd
                     << ", Ru = " << ru << "; " << table.curvatureProvenance());
    CHECK(table.curvatureProvenance().find("circle fitted") != std::string::npos);
    CHECK(table.bartzCurvatureRadius() > 0.5 * rd);
    CHECK(table.bartzCurvatureRadius() < 1.5 * ru);
  }

  SECTION("the exit wall angle is measured, not defaulted") {
    // This is the one that would bite silently: NozzleFlow turns the exit
    // angle into the divergence loss, so a contour that arrived without one
    // would quietly be charged an 8-degree bell's loss whatever its shape.
    CHECK(table.spec().bell_exit_angle > 0.0);
    CHECK(table.spec().bell_exit_angle ==
          Approx(analytic.spec().bell_exit_angle).margin(1.5));
  }

  SECTION("the wall departs from the analytic one only by chord error") {
    // A polyline through samples of a smooth curve misses it between samples.
    // The gap is the chord error of an arc, h^2/(8R), so it is worst where the
    // curvature is tightest -- the Rao downstream throat arc, Rd = 0.382 rt.
    // Asserting that, rather than a tolerance picked to pass, is what
    // distinguishes "the representation has a known accuracy limit" from "the
    // conversion has a bug".
    double worst = 0.0;
    double worst_x = 0.0;
    double spacing_at_worst = 0.0;
    const double x0 = analytic.stations().front().x;
    for (int i = 0; i <= 2000; ++i) {
      const double x = x0 + (analytic.exitPosition() - x0) * i / 2000.0;
      const double err = std::abs(table.radius(x) - analytic.radius(x));
      if (err > worst) { worst = err; worst_x = x; }
    }
    for (std::size_t i = 1; i < table.stations().size(); ++i) {
      if (table.stations()[i].x >= worst_x) {
        spacing_at_worst = table.stations()[i].x - table.stations()[i - 1].x;
        break;
      }
    }
    const double rd = 0.382 * analytic.throatRadius();
    const double predicted = spacing_at_worst * spacing_at_worst / (8.0 * rd);
    INFO("worst " << worst * 1e6 << " um at x = " << worst_x << " m; spacing "
                  << spacing_at_worst * 1e3 << " mm; chord error predicted "
                  << predicted * 1e6 << " um");

    // It is small in absolute terms.
    CHECK(worst < 1e-3 * analytic.throatRadius());
    // It sits in the throat region, where the wall curves hardest.
    CHECK(std::abs(worst_x - analytic.throatPosition()) <
          0.5 * analytic.throatRadius());
    // And it is the size chord error says it should be, not something else
    // wearing chord error's clothes.
    CHECK(worst == Approx(predicted).epsilon(0.6));
  }
}

TEST_CASE("a tabulated contour is rejected when it is not a nozzle",
          "[nozzle][geometry][contour]") {
  const auto good = tabulate(ignis::NozzleGeometry::build(m1Spec()));

  SECTION("too few samples") {
    std::vector<ignis::ContourSample> few(good.begin(), good.begin() + 5);
    REQUIRE_THROWS_WITH(ignis::NozzleGeometry::fromContour(few, m1Spec()),
                        Catch::Matchers::ContainsSubstring("at least 9 samples"));
  }

  SECTION("a repeated axial station") {
    auto bad = good;
    bad[40].x = bad[39].x;
    REQUIRE_THROWS_WITH(ignis::NozzleGeometry::fromContour(bad, m1Spec()),
                        Catch::Matchers::ContainsSubstring("strictly increasing"));
  }

  SECTION("a non-positive radius") {
    auto bad = good;
    bad[40].r = 0.0;
    REQUIRE_THROWS_WITH(ignis::NozzleGeometry::fromContour(bad, m1Spec()),
                        Catch::Matchers::ContainsSubstring("non-positive radius"));
  }

  SECTION("a bump in the divergent section") {
    // The kind of defect a noisy STL silhouette produces.  It must be
    // reported, with the station, rather than quietly smoothed away.
    auto bad = good;
    const std::size_t n = bad.size();
    bad[n - 10].r = bad[n - 11].r * 0.98;
    REQUIRE_THROWS_WITH(
        ignis::NozzleGeometry::fromContour(bad, m1Spec()),
        Catch::Matchers::ContainsSubstring("not monotonically expanding"));
  }

  SECTION("a throat at the very end") {
    // A pure converging contour is not a nozzle Ignis can expand through.
    std::vector<ignis::ContourSample> converging;
    for (const auto& s : good) {
      if (s.x > 0.3455) break;
      converging.push_back(s);
    }
    REQUIRE(converging.size() > 9);
    REQUIRE_THROWS_WITH(ignis::NozzleGeometry::fromContour(converging, m1Spec()),
                        Catch::Matchers::ContainsSubstring("must lie inside it"));
  }
}
