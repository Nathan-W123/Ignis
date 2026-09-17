// SPDX-License-Identifier: MIT
#pragma once
/// \file NozzleGeometry.hpp
/// \brief Axisymmetric converging-diverging nozzle contour.
///
/// CONTOUR CONSTRUCTION
/// --------------------
/// The contour follows the classical layout of Huzel & Huang, "Modern
/// Engineering for Design of Liquid-Propellant Rocket Engines", AIAA Progress
/// in Astronautics and Aeronautics Vol. 147 (1992), Ch. 4, and Rao's
/// parabolic bell approximation (G. V. R. Rao, "Exhaust Nozzle Contour for
/// Optimum Thrust", Jet Propulsion 28, 377-382, 1958).
///
/// The axial coordinate x increases from the injector face (x = 0) towards the
/// exit.  Six segments are laid down in order:
///
///   1. cylindrical chamber                r = Rc,          0 <= x <= Lc
///   2. chamber-to-cone fillet arc         radius R1
///   3. straight converging cone           half-angle beta
///   4. throat upstream arc                radius Ru  (default 1.5 Rt)
///   5. throat downstream arc              radius Rd  (default 0.382 Rt)
///   6. divergent: straight cone, or a quadratic Bezier bell that leaves the
///      throat arc at theta_n and arrives at the exit at theta_e
///
/// Every joint is C1 continuous by construction (each arc is tangent to its
/// neighbours), which `checkContinuity()` verifies numerically.  The bell is a
/// Rao-style *approximation*: a true optimum contour requires the method of
/// characteristics, so theta_n and theta_e are user inputs rather than being
/// read off Rao's charts.  This is documented in docs/limitations.md.
///
/// The contour is smooth and monotone in area on each side of the throat,
/// which is all the quasi-1D solver requires.  It is not a manufacturing
/// drawing: no injector face detail, fillet radii at the chamber head, wall
/// thickness taper or throat insert geometry are represented.

#include <cmath>
#include <string>
#include <vector>

#include "ignis/core/Constants.hpp"

namespace ignis {

enum class DivergentType { kConical, kBell };
std::string toString(DivergentType t);
DivergentType divergentTypeFromString(const std::string& s);

/// Which analytic piece a station belongs to.
enum class NozzleSegment {
  kChamber, kChamberFillet, kConvergingCone, kThroatUpstream, kThroatDownstream, kDivergent
};
std::string toString(NozzleSegment s);

/// User-facing geometry inputs.  Exactly one of `throat_radius` and
/// `throat_area` must be positive; likewise for `chamber_radius` and
/// `contraction_ratio`.
struct NozzleGeometrySpec {
  double throat_radius = 0.0;          ///< m
  double throat_area = 0.0;            ///< m^2 (alternative to throat_radius)
  double chamber_radius = 0.0;         ///< m
  double contraction_ratio = 0.0;      ///< Ac/At (alternative to chamber_radius)
  double chamber_length = 0.0;         ///< m, cylindrical portion
  double converging_half_angle = 30.0; ///< deg
  double chamber_fillet_ratio = 0.5;   ///< R1 / Rc
  double throat_upstream_ratio = 1.5;  ///< Ru / Rt
  double throat_downstream_ratio = 0.382;  ///< Rd / Rt  (Rao)
  double expansion_ratio = 40.0;       ///< Ae / At
  DivergentType divergent = DivergentType::kBell;
  double cone_half_angle = 15.0;       ///< deg, conical divergent only
  double bell_length_fraction = 0.8;   ///< L / L(15 deg conical)
  double bell_initial_angle = 0.0;     ///< theta_n, deg; 0 => default schedule
  double bell_exit_angle = 0.0;        ///< theta_e, deg; 0 => default schedule
  int num_stations = 400;              ///< stations written to the contour
  /// Which curvature radius Bartz's correlation should use at the throat:
  /// "upstream", "downstream" or "mean".
  std::string bartz_curvature = "mean";
};

/// One point on the contour.
struct ContourPoint {
  double x = 0.0;           ///< m from the injector face
  double r = 0.0;           ///< m, wall radius
  double area = 0.0;        ///< m^2
  double area_ratio = 0.0;  ///< A / At
  double drdx = 0.0;        ///< wall slope
  NozzleSegment segment = NozzleSegment::kChamber;
};

/// A validated, sampled nozzle contour.
class NozzleGeometry {
 public:
  // The contour segments are an implementation detail held by value, so the
  // special members are declared here and defined where Segment is complete.
  NozzleGeometry();
  ~NozzleGeometry();
  NozzleGeometry(const NozzleGeometry&);
  NozzleGeometry& operator=(const NozzleGeometry&);
  NozzleGeometry(NozzleGeometry&&) noexcept;
  NozzleGeometry& operator=(NozzleGeometry&&) noexcept;

  /// Build and validate a contour.  Throws ConfigError on any inconsistency.
  static NozzleGeometry build(const NozzleGeometrySpec& spec);

  const NozzleGeometrySpec& spec() const { return spec_; }
  const std::vector<ContourPoint>& stations() const { return stations_; }

  double throatRadius() const { return rt_; }
  double throatArea() const { return at_; }
  double chamberRadius() const { return rc_; }
  double chamberArea() const { return constants::pi * rc_ * rc_; }
  double exitRadius() const { return re_; }
  double exitArea() const { return constants::pi * re_ * re_; }
  double expansionRatio() const { return exitArea() / at_; }
  double contractionRatio() const { return chamberArea() / at_; }
  double throatPosition() const { return xt_; }
  double exitPosition() const { return xe_; }
  double divergentLength() const { return xe_ - xt_; }
  /// Length of a 15-degree conical divergent with the same expansion ratio, m.
  double referenceConicalLength() const { return l15_; }
  /// Throat curvature radius handed to the Bartz correlation, m.
  double bartzCurvatureRadius() const { return bartz_rc_; }
  /// Enclosed volume from the injector face to the throat, m^3.
  double chamberVolume() const { return v_chamber_; }
  /// Characteristic chamber length L* = V_chamber / At, m.
  double characteristicLength() const { return v_chamber_ / at_; }
  /// Wetted (lateral) surface area from x_from to x_to, m^2.
  double wettedArea(double x_from, double x_to) const;

  /// Wall radius at an arbitrary axial position, m.  Evaluated analytically
  /// from the underlying segment, not interpolated.
  double radius(double x) const;
  double area(double x) const { const double r = radius(x); return constants::pi * r * r; }
  /// Which segment contains x.
  NozzleSegment segmentAt(double x) const;

  /// Largest C0 and C1 mismatch across the segment joints.  Used by the test
  /// suite to confirm the contour really is smooth.
  void checkContinuity(double& max_radius_jump, double& max_slope_jump) const;

  /// Human-readable summary.
  std::string summary() const;

 private:
  struct Segment;
  const Segment& segmentFor(double x) const;

  NozzleGeometrySpec spec_;
  std::vector<Segment> segments_;
  std::vector<ContourPoint> stations_;
  double rt_ = 0.0, rc_ = 0.0, re_ = 0.0, at_ = 0.0;
  double xt_ = 0.0, xe_ = 0.0, l15_ = 0.0;
  double bartz_rc_ = 0.0, v_chamber_ = 0.0;
  double theta_n_ = 0.0, theta_e_ = 0.0;
};

}  // namespace ignis
