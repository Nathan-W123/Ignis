// SPDX-License-Identifier: MIT
/// \file test_boundary_layer_losses.cpp
/// \brief The boundary layer's discharge coefficient and thrust deficit.

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cmath>

#include "TestHelpers.hpp"
#include "ignis/combustion/Chamber.hpp"
#include "ignis/core/Constants.hpp"
#include "ignis/core/Exceptions.hpp"
#include "ignis/thermal/BoundaryLayerLosses.hpp"

using namespace ignis;
using namespace ignis_test;
using Catch::Approx;

namespace {

struct Fixture {
  SpeciesDatabase db = rocketDatabase();
  EquilibriumSolver solver{db};
  TransportModel transport{db};
  PropellantLibrary lib =
      PropellantLibrary::loadYaml(sourceDir() + "/data/propellants/ignis_propellants.yaml");
  PropellantMixture mix{lib.at("LOX"), lib.at("LCH4"), 3.4, 90.18, 111.66};
  CombustionChamber chamber{solver, CompositionModel::kEquilibrium};
  ChamberResult ch = chamber.solve(mix, 5.5e6, 0.96);
  NozzleGeometry geom = [] {
    NozzleGeometrySpec s;
    s.throat_radius = 0.070;
    s.contraction_ratio = 2.8;
    s.chamber_length = 0.22;
    s.expansion_ratio = 20.0;
    s.num_stations = 400;
    return NozzleGeometry::build(s);
  }();
  NozzleFlow flow = chamber.makeFlow(ch);

  BoundaryLayerLossResult losses(double t_wall) const {
    BoundaryLayerLossOptions o;
    o.uncooled_wall_temperature = t_wall;
    return computeBoundaryLayerLosses(flow, geom, ch, transport, nullptr, o);
  }
  NozzlePerformance perform(const BoundaryLayerLossResult* bl, double p_ambient) const {
    NozzlePerformanceOptions o;
    o.eta_c_star = 0.96;
    if (bl != nullptr) {
      o.viscous.apply = true;
      o.viscous.throat_displacement_thickness = bl->throat_displacement_thickness;
      o.viscous.exit_displacement_thickness = bl->exit_displacement_thickness;
      o.viscous.exit_momentum_thickness = bl->exit_momentum_thickness;
    }
    return evaluateNozzle(flow, geom, p_ambient, o);
  }
};

}  // namespace

TEST_CASE("the boundary layer narrows the throat and costs thrust", "[nozzle][boundary_layer]") {
  Fixture f;
  const auto bl = f.losses(1000.0);

  SECTION("a hot-walled layer displaces flow and carries a momentum deficit") {
    CHECK(bl.throat_displacement_thickness > 0.0);
    CHECK(bl.exit_displacement_thickness > bl.throat_displacement_thickness);
    CHECK(bl.exit_momentum_thickness > 0.0);
    CHECK(bl.discharge_coefficient < 1.0);
    CHECK(bl.discharge_coefficient > 0.99);
    CHECK(bl.effective_area_ratio < bl.geometric_area_ratio);
    CHECK(bl.discharge_coefficient ==
          Approx(std::pow(1.0 - bl.throat_displacement_thickness / f.geom.throatRadius(), 2)));
  }
  SECTION("with no jacket temperature every station assumes the uncooled wall, and says so") {
    CHECK(bl.cooled_fraction == 0.0);
    bool said = false;
    for (const auto& w : bl.warnings)
      if (w.find("assume 1000 K") != std::string::npos) said = true;
    CHECK(said);
  }

  const auto inviscid = f.perform(nullptr, 101325.0);
  const auto viscous = f.perform(&bl, 101325.0);

  SECTION("the thrust is the sum of its stated parts, exactly") {
    REQUIRE(viscous.viscous);
    const double lam = viscous.lambda_divergence * viscous.eta_nozzle;
    CHECK(viscous.mdot == Approx(bl.discharge_coefficient * inviscid.mdot).epsilon(1e-14));
    CHECK(viscous.momentum_deficit ==
          Approx(2.0 * constants::pi * f.geom.exitRadius() *
                 f.flow.atAreaRatio(viscous.effective_area_ratio, true).gas.rho *
                 viscous.u_exit * viscous.u_exit * bl.exit_momentum_thickness)
              .epsilon(1e-12));
    CHECK(viscous.thrust ==
          Approx(lam * (viscous.mdot * viscous.u_exit - viscous.momentum_deficit) +
                 (viscous.p_exit - 101325.0) * viscous.exit_area)
              .epsilon(1e-12));
  }
  SECTION("the inviscid and ideal fields are those of the inviscid nozzle") {
    CHECK(viscous.mdot_inviscid == Approx(inviscid.mdot).epsilon(1e-14));
    CHECK(viscous.thrust_inviscid == Approx(inviscid.thrust).epsilon(1e-12));
    CHECK(viscous.isp_vacuum_inviscid == Approx(inviscid.isp_vacuum).epsilon(1e-12));
    CHECK(viscous.thrust_ideal == Approx(inviscid.thrust_ideal).epsilon(1e-12));
    CHECK(viscous.isp_ideal == Approx(inviscid.isp_ideal).epsilon(1e-12));
  }
  SECTION("it costs specific impulse, by the order of a percent") {
    const double loss = 1.0 - viscous.isp_vacuum / inviscid.isp_vacuum;
    INFO("vacuum Isp loss " << 100.0 * loss << " %");
    CHECK(loss > 0.002);
    CHECK(loss < 0.03);
    // A displaced exit expands less, so the core leaves at a higher pressure.
    CHECK(viscous.p_exit > inviscid.p_exit);
  }
  SECTION("the vacuum-to-ambient family stays exact") {
    CHECK(viscous.thrust_sea_level ==
          Approx(viscous.isp_vacuum * viscous.mdot * constants::g0 -
                 constants::atm * viscous.exit_area).epsilon(1e-12));
  }
}

TEST_CASE("a colder wall makes a denser, thinner displacement layer", "[nozzle][boundary_layer]") {
  // The gas next to a cooled wall is denser than the edge gas, so it carries
  // more of the mass the velocity deficit removed and delta* shrinks.  In the
  // limit it can go negative (Elliott, Bartz & Silver, JPL TR 32-387, Fig. 11,
  // for their sample nozzle); on this engine a 150 K wall halves it and more,
  // but it stays positive.  The code does not clip either way.
  Fixture f;
  const auto hot = f.losses(1500.0);
  const auto warm = f.losses(600.0);
  const auto cold = f.losses(150.0);
  INFO("throat delta*: 1500 K " << hot.throat_displacement_thickness * 1e3 << " mm, 600 K "
       << warm.throat_displacement_thickness * 1e3 << " mm, 150 K "
       << cold.throat_displacement_thickness * 1e3 << " mm");
  CHECK(warm.throat_displacement_thickness < hot.throat_displacement_thickness);
  CHECK(cold.throat_displacement_thickness < warm.throat_displacement_thickness);
  CHECK(cold.discharge_coefficient > hot.discharge_coefficient);
}

TEST_CASE("boundary-layer inputs are checked", "[nozzle][boundary_layer]") {
  Fixture f;
  BoundaryLayerLossOptions o;
  o.uncooled_wall_temperature = 0.0;
  CHECK_THROWS_AS(computeBoundaryLayerLosses(f.flow, f.geom, f.ch, f.transport, nullptr, o),
                  ConfigError);
  NozzlePerformanceOptions po;
  po.viscous.apply = true;
  po.viscous.throat_displacement_thickness = 2.0 * f.geom.throatRadius();
  CHECK_THROWS_AS(evaluateNozzle(f.flow, f.geom, 0.0, po), ConfigError);
}
