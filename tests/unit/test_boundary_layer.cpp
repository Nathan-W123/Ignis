// SPDX-License-Identifier: MIT
/// \file test_boundary_layer.cpp
/// \brief Verification of the integral turbulent boundary layer.
///
/// Every case here has a closed-form answer.  With properties held fixed the
/// momentum and energy integrals are Bernoulli equations in theta and Delta,
/// and Y = theta^(5/4) turns each into a linear first-order ODE whose
/// integrating factor is a power of u_e, r or the driving temperature
/// difference.  So the march is checked against exact solutions -- for a flat
/// plate, for an accelerating stream, for a widening duct and for a wall whose
/// temperature varies -- rather than against itself.  Agreement with
/// measurements is a separate question, answered in
/// tests/validation/test_heat_transfer_validation.cpp.

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <vector>

#include "ignis/core/Exceptions.hpp"
#include "ignis/thermal/BoundaryLayer.hpp"

using namespace ignis;
using Catch::Approx;

namespace {

constexpr double kRho = 1.0;       // kg/m^3
constexpr double kMu = 4.0e-5;     // Pa s
constexpr double kCp = 1100.0;     // J/(kg K)
constexpr double kT = 1000.0;      // K, edge static temperature

/// Properties that ignore temperature, so the closure's coefficients are
/// exactly the constants of the analytic solutions.
ReferenceProperties constantProperties(double prandtl) {
  return [prandtl](std::size_t, double) {
    ReferenceState s;
    s.rho = kRho;
    s.viscosity = kMu;
    s.cp = kCp;
    s.prandtl = prandtl;
    return s;
  };
}

/// A low-speed edge with everything uniform; callers then vary one thing.
std::vector<BoundaryLayerEdge> uniformEdge(double length, int n, double u = 100.0) {
  std::vector<BoundaryLayerEdge> e(static_cast<std::size_t>(n));
  for (int i = 0; i < n; ++i) {
    auto& s = e[static_cast<std::size_t>(i)];
    s.x = length * i / (n - 1);
    s.radius = 1.0;
    s.u = u;
    s.rho = kRho;
    s.T = kT;
    s.p = 1.0e5;
    s.mach = 0.0;
    s.viscosity = kMu;
    s.t_adiabatic_wall = kT;
    s.t_wall = 0.5 * kT;
  }
  return e;
}

/// 0.0128 (rho u / mu)^(-1/4): the closure's coefficient with rho* = rho_e.
double closureCoefficient(double u) { return 0.0128 * std::pow(kRho * u / kMu, -0.25); }

}  // namespace

TEST_CASE("the shape factor has the textbook limits", "[thermal][boundary_layer]") {
  SECTION("incompressible 1/n profile: H = (n + 2) / n") {
    // delta*/delta = 1/(n+1) and theta/delta = n/((n+1)(n+2)).
    CHECK(compressibleShapeFactor(1.0, 1.0, 7.0) == Approx(9.0 / 7.0).epsilon(5e-5));
    CHECK(compressibleShapeFactor(1.0, 1.0, 9.0) == Approx(11.0 / 9.0).epsilon(5e-5));
  }
  SECTION("a cold wall packs mass near the wall and lowers H") {
    // Denser gas where the velocity deficit is: the mass deficit shrinks
    // faster than the momentum deficit.
    const double h_cold = compressibleShapeFactor(0.3, 1.0);
    const double h_adiabatic = compressibleShapeFactor(1.0, 1.0);
    CHECK(h_cold < h_adiabatic);
    CHECK(h_cold > 0.0);
  }
  SECTION("aerodynamic heating thickens the layer and raises H") {
    // An adiabatic wall at Mach 3 recovers to about 2.6 T_e.
    CHECK(compressibleShapeFactor(2.6, 2.6) > compressibleShapeFactor(1.0, 1.0) + 0.5);
  }
  SECTION("Eckert's reference temperature") {
    CHECK(eckertReferenceTemperature(1000.0, 1000.0, 1000.0) == Approx(1000.0));
    CHECK(eckertReferenceTemperature(2000.0, 600.0, 2500.0) ==
          Approx(2000.0 + 0.5 * (600.0 - 2000.0) + 0.22 * 500.0));
  }
  SECTION("bad input is refused") {
    CHECK_THROWS_AS(compressibleShapeFactor(0.0, 1.0), ConfigError);
    CHECK_THROWS_AS(compressibleShapeFactor(1.0, 1.0, 0.0), ConfigError);
  }
}

TEST_CASE("on a flat plate the march is the 1/7-law solution", "[thermal][boundary_layer]") {
  // d theta/dx = A theta^-1/4  =>  theta^(5/4) = theta_0^(5/4) + (5/4) A x.
  const double u = 100.0;
  const auto edge = uniformEdge(0.5, 201, u);
  const double A = closureCoefficient(u);

  SECTION("momentum and energy thickness, station by station") {
    const double prandtl = 0.7;
    BoundaryLayerOptions o;
    o.initial_momentum_thickness = 2.0e-5;
    o.initial_enthalpy_thickness = 1.0e-5;
    const auto sol = marchBoundaryLayer(edge, constantProperties(prandtl), o);
    const double B = A * std::pow(prandtl, -2.0 / 3.0);
    for (const auto& p : sol.points) {
      const double theta = std::pow(std::pow(2.0e-5, 1.25) + 1.25 * A * p.x, 0.8);
      const double delta = std::pow(std::pow(1.0e-5, 1.25) + 1.25 * B * p.x, 0.8);
      INFO("x = " << p.x);
      // RK4 with four substeps per station: ~1e-6 where the start is steep.
      CHECK(p.momentum_thickness == Approx(theta).epsilon(1e-5));
      CHECK(p.enthalpy_thickness == Approx(delta).epsilon(1e-5));
      // Reported closure quantities are the closure at the marched thickness.
      CHECK(p.skin_friction == Approx(2.0 * A * std::pow(p.momentum_thickness, -0.25)).epsilon(1e-12));
      CHECK(p.h_gas == Approx(B * std::pow(p.enthalpy_thickness, -0.25) * kRho * u * kCp)
                           .epsilon(1e-12));
      CHECK(p.displacement_thickness == Approx(p.shape_factor * p.momentum_thickness));
    }
  }

  SECTION("grown from a negligible start it is theta/x = 0.0366 Re_x^-0.2") {
    // The textbook turbulent flat-plate result (the constant is 0.016^0.8).
    BoundaryLayerOptions o;
    o.initial_momentum_thickness = 1.0e-8;
    o.initial_enthalpy_thickness = 1.0e-8;
    const auto sol = marchBoundaryLayer(edge, constantProperties(1.0), o);
    const auto& end = sol.points.back();
    const double re_x = kRho * u * end.x / kMu;
    CHECK(end.momentum_thickness / end.x ==
          Approx(std::pow(0.016, 0.8) * std::pow(re_x, -0.2)).epsilon(1e-3));
    // Pr = 1 and equal starts: the energy integral is the momentum one.
    CHECK(end.enthalpy_thickness == Approx(end.momentum_thickness).epsilon(1e-9));
    CHECK(end.stanton == Approx(0.5 * end.skin_friction).epsilon(1e-9));
  }
}

TEST_CASE("the march carries the pressure-gradient, area and wall-temperature terms exactly",
          "[thermal][boundary_layer]") {
  // With Y = theta^(5/4) the momentum integral is
  //   Y' + (5/4)[(2 + H - M^2) u'/u + r'/r] Y = (5/4) a(x),
  // and the energy integral is the same with Delta, St and ln(G cp dT).
  const double L = 0.2;
  const int n = 401;

  SECTION("an accelerating stream: theta u^c is the integral of a u^c") {
    auto edge = uniformEdge(L, n);
    const double u0 = 100.0;
    for (auto& e : edge) {
      e.u = u0 * (1.0 + e.x / L);
      // Constant-property incompressible edge; T_aw a hair above T_e so the
      // energy integral is defined and H is the incompressible 9/7.
      e.t_adiabatic_wall = kT * (1.0 + 1e-9);
      e.t_wall = kT;
    }
    BoundaryLayerOptions o;
    o.initial_momentum_thickness = 3.0e-4;
    o.initial_enthalpy_thickness = 3.0e-4;
    const auto sol = marchBoundaryLayer(edge, constantProperties(1.0), o);
    const double H = compressibleShapeFactor(1.0, 1.0 + 1e-9);
    const double c = 1.25 * (2.0 + H);
    const double k = 0.0128 * std::pow(kRho / kMu, -0.25);   // a = k u^-1/4
    for (std::size_t i = 0; i < sol.points.size(); i += 40) {
      const double x = sol.points[i].x;
      const double s = 1.0 + x / L;
      const double integral = k * std::pow(u0, c - 0.25) * L * (std::pow(s, c + 0.75) - 1.0) / (c + 0.75);
      const double Y = (std::pow(3.0e-4, 1.25) * std::pow(u0, c) + 1.25 * integral) /
                       std::pow(u0 * s, c);
      INFO("x = " << x);
      CHECK(sol.points[i].momentum_thickness == Approx(std::pow(Y, 0.8)).epsilon(2e-4));
    }
    // Acceleration thins the layer: the same plate at the inlet velocity is thicker.
    const auto flat = marchBoundaryLayer(uniformEdge(L, n, u0), constantProperties(1.0), o);
    CHECK(sol.points.back().momentum_thickness < flat.points.back().momentum_thickness);
  }

  SECTION("a widening duct: Y r^(5/4) is the integral of a r^(5/4)") {
    auto edge = uniformEdge(L, n);
    const double r0 = 0.05;
    for (auto& e : edge) e.radius = r0 * (1.0 + e.x / L);
    BoundaryLayerOptions o;
    o.initial_momentum_thickness = 3.0e-4;
    o.initial_enthalpy_thickness = 3.0e-4;
    const auto sol = marchBoundaryLayer(edge, constantProperties(1.0), o);
    const double a = closureCoefficient(100.0);
    for (std::size_t i = 0; i < sol.points.size(); i += 40) {
      const double s = 1.0 + sol.points[i].x / L;
      // (Y r^(5/4))' = (5/4) a r^(5/4); the same holds for the energy integral.
      const double Y = (std::pow(3.0e-4, 1.25) + 1.25 * a * L * (std::pow(s, 2.25) - 1.0) / 2.25) /
                       std::pow(s, 1.25);
      INFO("x = " << sol.points[i].x);
      CHECK(sol.points[i].momentum_thickness == Approx(std::pow(Y, 0.8)).epsilon(2e-4));
      CHECK(sol.points[i].enthalpy_thickness == Approx(std::pow(Y, 0.8)).epsilon(2e-4));
    }
  }

  SECTION("a wall whose temperature falls: the driving difference grows downstream") {
    // dT = T_aw - T_w rising linearly; (Z dT^(5/4))' = (5/4) b dT^(5/4).
    auto edge = uniformEdge(L, n);
    const double dT0 = 400.0;
    for (auto& e : edge) e.t_wall = kT - dT0 * (1.0 + e.x / L);
    BoundaryLayerOptions o;
    o.initial_momentum_thickness = 3.0e-4;
    o.initial_enthalpy_thickness = 3.0e-4;
    const auto sol = marchBoundaryLayer(edge, constantProperties(1.0), o);
    const double b = closureCoefficient(100.0);
    for (std::size_t i = 0; i < sol.points.size(); i += 40) {
      const double s = 1.0 + sol.points[i].x / L;
      const double Z = (std::pow(3.0e-4, 1.25) + 1.25 * b * L * (std::pow(s, 2.25) - 1.0) / 2.25) /
                       std::pow(s, 1.25);
      INFO("x = " << sol.points[i].x);
      CHECK(sol.points[i].enthalpy_thickness == Approx(std::pow(Z, 0.8)).epsilon(2e-4));
    }
  }
}

TEST_CASE("the layer starts as the thinnest turbulent layer and then forgets it",
          "[thermal][boundary_layer]") {
  const auto edge = uniformEdge(0.5, 201);

  SECTION("by default the first station sits at Preston's Re_theta = 320") {
    const auto sol = marchBoundaryLayer(edge, constantProperties(0.7));
    CHECK(sol.points.front().re_theta == Approx(320.0).epsilon(1e-12));
    CHECK(sol.points.front().momentum_thickness == Approx(320.0 * kMu / (kRho * 100.0)));
    CHECK(sol.points.front().enthalpy_thickness == sol.points.front().momentum_thickness);
    // Starting on the minimum must not trip the "below the minimum" flag.
    CHECK(sol.min_re_theta == Approx(320.0).epsilon(1e-12));
    for (const auto& w : sol.warnings) CHECK(w.find("Reynolds number falls") == std::string::npos);
  }

  SECTION("an explicit start below the minimum is flagged") {
    BoundaryLayerOptions o;
    o.initial_momentum_thickness = 1.0e-7;
    o.initial_enthalpy_thickness = 1.0e-7;
    const auto sol = marchBoundaryLayer(edge, constantProperties(0.7), o);
    bool flagged = false;
    for (const auto& w : sol.warnings)
      if (w.find("Reynolds number falls") != std::string::npos) flagged = true;
    CHECK(flagged);
  }

  SECTION("downstream, a thinner start changes nothing that matters") {
    // Elliott, Bartz & Silver found the initial momentum thickness has little
    // effect on downstream heat flux.  At a combustion chamber's mass flux
    // (G/mu = 2.5e7 per metre here, against 1.6e7 in the M1's chamber) the
    // minimum turbulent layer is ~0.01 mm with a virtual origin ~3 mm
    // upstream, so a start ten times thinner moves h at 0.5 m by ~0.1 %.  (On
    // a low-speed plate the same start is a sizeable fraction of the length;
    // the claim is about rocket Reynolds numbers, not every flow.)
    const auto fast = uniformEdge(0.5, 201, 1000.0);
    const auto base = marchBoundaryLayer(fast, constantProperties(0.7));
    BoundaryLayerOptions o;
    o.initial_momentum_thickness = 0.1 * base.points.front().momentum_thickness;
    o.initial_enthalpy_thickness = o.initial_momentum_thickness;
    const auto thin = marchBoundaryLayer(fast, constantProperties(0.7), o);
    CHECK(thin.points.back().h_gas == Approx(base.points.back().h_gas).epsilon(5e-3));
    CHECK(thin.points.back().h_gas > base.points.back().h_gas);
  }
}

TEST_CASE("the boundary layer reports the limits of its own closure", "[thermal][boundary_layer]") {
  SECTION("strong acceleration raises the laminarisation flag") {
    // K = (nu / u^2) du/dx = 4e-5 * 1000 / 100^2 = 4e-6 at the inlet.
    auto edge = uniformEdge(0.1, 101);
    for (auto& e : edge) e.u = 100.0 * (1.0 + e.x / 0.1);
    const auto sol = marchBoundaryLayer(edge, constantProperties(0.7));
    // (One-sided difference at the inlet: ln(1.01)/0.001 against 10.)
    CHECK(sol.max_acceleration == Approx(4.0e-6).epsilon(1e-2));
    CHECK(sol.max_acceleration_x == Approx(0.0).margin(1e-12));
    bool flagged = false;
    for (const auto& w : sol.warnings)
      if (w.find("laminarise") != std::string::npos) flagged = true;
    CHECK(flagged);
  }
  SECTION("a gentle one does not") {
    auto edge = uniformEdge(0.1, 101);
    for (auto& e : edge) e.u = 100.0 * (1.0 + 0.1 * e.x / 0.1);
    const auto sol = marchBoundaryLayer(edge, constantProperties(0.7));
    CHECK(sol.max_acceleration < 3.0e-6);
    CHECK(sol.warnings.empty());
  }
}

TEST_CASE("the film coefficient a wall balance calls is the march's own closure",
          "[thermal][boundary_layer]") {
  // The cooling solver evaluates h at its trial wall temperature from the
  // marched enthalpy thickness; at the march's wall temperature that must be
  // exactly the h the march reported.
  auto edge = uniformEdge(0.3, 61);
  for (auto& e : edge) e.u = 100.0 * (1.0 + e.x);
  const auto props = constantProperties(0.7);
  const auto sol = marchBoundaryLayer(edge, props);
  for (std::size_t i = 0; i < edge.size(); ++i) {
    const double ts = eckertReferenceTemperature(edge[i].T, edge[i].t_wall, edge[i].t_adiabatic_wall);
    const double h = boundaryLayerFilmCoefficient(edge[i], sol.points[i].enthalpy_thickness,
                                                  edge[i].t_wall, props(i, ts));
    CHECK(h == Approx(sol.points[i].h_gas).epsilon(1e-12));
  }
  SECTION("interpolation between stations is linear and clamped") {
    const auto mid = sol.at(0.5 * (edge[3].x + edge[4].x));
    CHECK(mid.momentum_thickness ==
          Approx(0.5 * (sol.points[3].momentum_thickness + sol.points[4].momentum_thickness)));
    CHECK(sol.at(-1.0).x == sol.points.front().x);
    CHECK(sol.at(10.0).x == sol.points.back().x);
  }
}

TEST_CASE("the march refuses what it cannot integrate", "[thermal][boundary_layer]") {
  const auto props = constantProperties(0.7);
  auto edge = uniformEdge(0.3, 11);
  SECTION("too few stations") {
    edge.resize(2);
    CHECK_THROWS_AS(marchBoundaryLayer(edge, props), ConfigError);
  }
  SECTION("x not increasing") {
    edge[5].x = edge[4].x;
    CHECK_THROWS_AS(marchBoundaryLayer(edge, props), ConfigError);
  }
  SECTION("a wall at the adiabatic wall temperature has no driving difference") {
    edge[3].t_wall = edge[3].t_adiabatic_wall;
    CHECK_THROWS_AS(marchBoundaryLayer(edge, props), ConfigError);
  }
  SECTION("negative starting thickness") {
    BoundaryLayerOptions o;
    o.initial_momentum_thickness = -1.0;
    CHECK_THROWS_AS(marchBoundaryLayer(edge, props, o), ConfigError);
  }
}
