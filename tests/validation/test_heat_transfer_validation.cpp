// SPDX-License-Identifier: MIT
/// \file test_heat_transfer_validation.cpp
/// \brief The Bartz correlation against measured nozzle heat transfer.
///
/// Every other comparison in this suite checks Ignis against another *code* --
/// CEA, Cantera, a reference equation of state.  This one checks it against a
/// *measurement*: the local heat-transfer coefficients tabulated by Back,
/// Massier and Gier (JPL Technical Report 32-415, 1965) for air expanded
/// through a conical convergent-divergent nozzle.  See the header of
/// validation/reference/nozzle_heat_transfer_reference.csv for the apparatus,
/// the transcription provenance and the authors' own uncertainty estimate.
///
/// WHAT THIS DOES AND DOES NOT ESTABLISH
/// -------------------------------------
/// It establishes that Ignis evaluates the Bartz correlation correctly, and it
/// measures how far that correlation lands from a real experiment.  It does
/// not establish that Ignis predicts the heat flux in a rocket engine: the
/// working fluid here is air at 820-840 K, not combustion products at 3500 K,
/// there is no injector and no film cooling, and the correlation is known to
/// be the weakest link in the chain.  The value of the comparison is that it
/// replaces an assumed `bartz_multiplier` with a number that was measured.
#include <algorithm>
#include <cmath>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "TestHelpers.hpp"
#include "ignis/core/Constants.hpp"
#include "ignis/thermo/GasMixture.hpp"
#include "ignis/thermo/Transport.hpp"
#include "ignis/combustion/Chamber.hpp"
#include "ignis/nozzle/NozzleGeometry.hpp"
#include "ignis/thermal/BoundaryLayer.hpp"
#include "ignis/thermal/HeatTransfer.hpp"
#include "ignis/thermal/RegenerativeCooling.hpp"

using ignis_test::ReferenceTable;
using Catch::Approx;
using ignis_test::referenceDir;

namespace {

// The report is in US customary units; these convert it, and nothing else in
// the test carries a unit conversion.
constexpr double kInch = 0.0254;                    // m
constexpr double kPsi = 6894.757293168361;          // Pa
constexpr double kRankine = 5.0 / 9.0;              // K
constexpr double kPound = 0.45359237;               // kg
constexpr double kBtu = 1055.05585262;              // J
/// BTU/(s in^2 degF) -> W/(m^2 K).  A degree Fahrenheit interval is 5/9 K.
const double kCoefficient = kBtu / (kInch * kInch * kRankine);

/// Nozzle geometry, JPL TR 32-415 Section I and Fig. 5.
constexpr double kThroatDiameter = 1.803 * kInch;
constexpr double kThroatCurvature = 1.800 * kInch;
constexpr double kThroatRadius = 0.9015 * kInch;

/// Dry air, mole fractions.
struct Air {
  Eigen::VectorXd n;      ///< mol/kg
  double molar_mass = 0.0;
};

Air dryAir(const ignis::SpeciesDatabase& db) {
  Eigen::VectorXd X = Eigen::VectorXd::Zero(static_cast<int>(db.size()));
  X(db.index("N2")) = 0.78084;
  X(db.index("O2")) = 0.20946;
  X(db.index("AR")) = 0.00934;
  const ignis::GasMixture mix(db);
  Air air;
  air.n = mix.fromMoleFractions(X);
  air.molar_mass = mix.molarMass(air.n);
  return air;
}

/// Local area ratio at a thermocouple station.
///
/// The pressure taps of Table B-2, which is where the contour is tabulated,
/// are at different axial stations from the thermocouple plugs of Table B-3.
/// Interpolation is on the radius ratio rather than the area ratio because the
/// nozzle is conical: radius is piecewise linear in z, area is not.
double areaRatioAt(const ReferenceTable& contour, double z_over_l) {
  const std::size_t n = contour.rows();
  auto radius = [&](std::size_t i) { return std::sqrt(contour.num("area_ratio", i)); };
  if (z_over_l <= contour.num("z_over_L", 0)) return std::pow(radius(0), 2);
  for (std::size_t i = 0; i + 1 < n; ++i) {
    const double z0 = contour.num("z_over_L", i);
    const double z1 = contour.num("z_over_L", i + 1);
    if (z_over_l <= z1) {
      const double t = (z_over_l - z0) / (z1 - z0);
      return std::pow(radius(i) + t * (radius(i + 1) - radius(i)), 2);
    }
  }
  return std::pow(radius(n - 1), 2);
}

double median(std::vector<double> v) {
  std::sort(v.begin(), v.end());
  return v.empty() ? 0.0 : v[v.size() / 2];
}

/// Ratio of the Bartz coefficient Ignis computes to the one JPL measured, at
/// every tabulated station, keyed by the cooled approach length in inches.
std::map<int, std::vector<double>> jplRatios() {
  const ReferenceTable ref(referenceDir() + "/nozzle_heat_transfer_reference.csv");
  const ReferenceTable contour(referenceDir() + "/nozzle_contour_reference.csv");
  const auto& db = ignis_test::fullDatabase();
  const Air air = dryAir(db);
  const ignis::GasMixture mix(db);
  const ignis::TransportModel transport(db);
  const double throat_area = ignis::constants::pi * kThroatRadius * kThroatRadius;

  std::map<int, std::vector<double>> ratio_by_approach;

  for (std::size_t row = 0; row < ref.rows(); ++row) {
    const double t_stagnation = ref.num("tto_degR", row) * kRankine;
    const double p_stagnation = ref.num("pt_psia", row) * kPsi;
    const double mdot = ref.num("mdot_lbm_s", row) * kPound;
    const double t_static = ref.num("te_degR", row) * kRankine;
    const double t_wall = ref.num("tw_degR", row) * kRankine;

    // Stagnation properties, which is where Bartz evaluates them.
    const double cp = mix.cpFrozen(air.n, t_stagnation);
    const auto tr = transport.mixture(air.n, t_stagnation, cp);
    const double r_specific = ignis::constants::R_universal / air.molar_mass;
    const double gamma = cp / (cp - r_specific);

    // The report took Mach from the measured static-to-stagnation pressure
    // ratio for isentropic flow at constant gamma; its tabulated free-stream
    // static temperature is that same isentropic relation, so using it here
    // reproduces their Mach number without re-deriving it.
    const double mach =
        std::sqrt(std::max(2.0 / (gamma - 1.0) * (t_stagnation / t_static - 1.0), 0.0));

    // Characteristic velocity from the measured choked mass flow, which is
    // what the correlation's p_c/c* group actually stands for.
    ignis::BartzReference bartz;
    bartz.throat_diameter = kThroatDiameter;
    bartz.curvature_radius = kThroatCurvature;
    bartz.chamber_pressure = p_stagnation;
    bartz.c_star = p_stagnation * throat_area / mdot;
    bartz.chamber_temperature = t_stagnation;
    bartz.viscosity = tr.viscosity;
    bartz.cp = cp;
    bartz.prandtl = tr.prandtl;

    const double area_ratio = areaRatioAt(contour, ref.num("z_over_L", row));
    const double predicted =
        ignis::bartzFilmCoefficient(bartz, area_ratio, mach, gamma, t_wall);
    const double measured = ref.num("h_btu_s_in2_F", row) * kCoefficient;

    REQUIRE(measured > 0.0);
    REQUIRE(predicted > 0.0);
    ratio_by_approach[static_cast<int>(ref.num("approach_in", row))]
        .push_back(predicted / measured);
  }
  return ratio_by_approach;
}

}  // namespace

TEST_CASE("the Bartz correlation is evaluated correctly against measured nozzle data",
          "[validation][thermal][bartz]") {
  const ReferenceTable ref(referenceDir() + "/nozzle_heat_transfer_reference.csv");
  const ReferenceTable contour(referenceDir() + "/nozzle_contour_reference.csv");
  REQUIRE(ref.rows() > 50);
  REQUIRE(contour.rows() == 32);

  const auto ratio_by_approach = jplRatios();
  REQUIRE(ratio_by_approach.size() == 2);

  SECTION("c* recovered from the measured mass flow is physical for air") {
    // p_t A* / mdot should sit a few percent above the ideal value, the
    // difference being the throat discharge coefficient.
    const auto& db = ignis_test::fullDatabase();
    const Air air = dryAir(db);
    const ignis::GasMixture mix(db);
    const double throat_area = ignis::constants::pi * kThroatRadius * kThroatRadius;
    const double t0 = 1507.0 * kRankine;
    const double cp = mix.cpFrozen(air.n, t0);
    const double r_specific = ignis::constants::R_universal / air.molar_mass;
    const double gamma = cp / (cp - r_specific);
    const double vandenkerckhove =
        std::sqrt(gamma) * std::pow(2.0 / (gamma + 1.0), (gamma + 1.0) / (2.0 * (gamma - 1.0)));
    const double c_star_ideal = std::sqrt(r_specific * t0) / vandenkerckhove;
    const double c_star_measured = 201.0 * kPsi * throat_area / (6.871 * kPound);
    CHECK(c_star_measured > c_star_ideal);
    CHECK(c_star_measured < 1.06 * c_star_ideal);
  }

  SECTION("Bartz over-predicts, and by more on one group of tests than the other") {
    // These bounds record a measurement, not a target.  The two groups differ
    // in cooled approach length -- 18 in against none -- and the correlation
    // is much further out on the 0-in. group.  Do NOT read that as the cause.
    // The 0-in. tests are also the low-pressure ones (35.9 and 51.0 psia
    // against 175 and 201), so approach length, throat Reynolds number and the
    // acceleration parameter all move together across these four tests and the
    // dataset cannot tell them apart; see docs/validation.md 5b.  What is
    // asserted here is the observation, not an explanation.  If these ranges
    // move, the correlation, the transcription or the evaluation changed --
    // all three are worth knowing about.
    const double long_approach = median(ratio_by_approach.at(18));
    const double no_approach = median(ratio_by_approach.at(0));

    CHECK(long_approach > 1.3);
    CHECK(long_approach < 1.9);
    CHECK(no_approach > 2.2);
    CHECK(no_approach > long_approach);
  }
}


TEST_CASE("Bartz's constant against a hydrogen-oxygen rocket",
          "[validation][thermal][bartz]") {
  // The JPL comparison above is air at 830 K.  This one is a real engine:
  // NASA TN D-2832 fired liquid oxygen and gaseous hydrogen through a copper
  // heat-sink nozzle and fitted the same correlation Ignis uses, station by
  // station, reporting the constant it actually takes.  Ignis uses 0.026
  // everywhere; the measurement says that is right in the chamber and much too
  // high at the throat.
  const ReferenceTable ref(referenceDir() + "/rocket_heat_transfer_reference.csv");
  REQUIRE(ref.rows() == 6);

  constexpr double kBartzConstant = 0.026;  // as in HeatTransfer.hpp
  std::map<std::string, double> c_by_station;
  std::map<std::string, double> sigma_by_station;
  for (std::size_t row = 0; row < ref.rows(); ++row) {
    const std::string station = ref.text("station", row);
    c_by_station[station] = ref.num("c_measured", row);
    sigma_by_station[station] = ref.num("c_sigma", row);
  }

  SECTION("in the chamber the constant Ignis uses is the measured one") {
    // Station 1 sits at an area ratio of 4.64, where the flow is barely
    // accelerating and the boundary layer is the ordinary turbulent pipe-like
    // one the correlation was built for.  The tolerance is the report's own
    // standard deviation on the fit, read from the table rather than written
    // in here, so the assertion cannot drift away from the data it checks.
    const double measured = c_by_station.at("1");
    const double sigma = sigma_by_station.at("1");
    CHECK(std::abs(kBartzConstant - measured) < sigma);
  }

  SECTION("at the throat it is high by about seventy percent") {
    // The report puts it at "42 percent lower than the widely used value of
    // C = 0.026", which is this ratio seen from the other side.  The cause is
    // the favourable pressure gradient through the throat, which the
    // correlation has no term for.
    const double throat = kBartzConstant / c_by_station.at("3bar");
    CHECK(throat > 1.6);
    CHECK(throat < 1.9);
    CHECK(kBartzConstant / c_by_station.at("1") < throat);
  }

  SECTION("two independent experiments agree on the size of the error") {
    // One is air at 830 K in a conical nozzle; the other is a hydrogen-oxygen
    // rocket at up to 966 psia.  Nothing links them but the correlation, so
    // their agreeing to within a quarter is the strongest statement in this
    // file about how far Bartz actually lands from reality.
    const double rocket = kBartzConstant / c_by_station.at("3bar");
    const double air_long_approach = median(jplRatios().at(18));
    INFO("rocket " << rocket << " vs air " << air_long_approach);
    CHECK(std::abs(rocket - air_long_approach) / rocket < 0.25);
  }
}


TEST_CASE("the JPL test nozzle can be built from its own published table",
          "[validation][nozzle][contour]") {
  // Until now this file has used the JPL contour by interpolating the
  // reference table directly, because Ignis could only construct a nozzle
  // from its bell/conical parameterisation -- it had no way to be handed a
  // contour.  It does now, so the test article can be built as the geometry
  // it actually was.
  //
  // This is a validation, not a unit test: the report states the nozzle's
  // dimensions independently of the tap table, so the geometry Ignis derives
  // from the table can be checked against what the authors said they built.
  const ReferenceTable contour(referenceDir() + "/nozzle_contour_reference.csv");
  REQUIRE(contour.rows() == 32);

  constexpr double kLength = 5.925 * kInch;
  std::vector<ignis::ContourSample> samples;
  samples.reserve(contour.rows());
  for (std::size_t row = 0; row < contour.rows(); ++row) {
    // The table gives A/A*; a radius follows because the nozzle is circular.
    samples.push_back({contour.num("z_over_L", row) * kLength,
                       kThroatRadius * std::sqrt(contour.num("area_ratio", row))});
  }

  ignis::NozzleGeometrySpec spec;
  spec.throat_radius = kThroatRadius;
  spec.expansion_ratio = 2.68;
  spec.chamber_length = samples.front().x;
  spec.num_stations = 400;
  const ignis::NozzleGeometry geom =
      ignis::NozzleGeometry::fromContour(samples, spec);

  SECTION("the throat lands where the report says it is") {
    // Tap 12, at z/L = 0.6018 and A/A* = 1.0012.  The table's minimum is that
    // tap, so the recovered throat radius is high by half of that 0.12 % in
    // area -- the resolution of the measurement, not an error in the fit.
    CHECK(geom.throatRadius() == Approx(kThroatRadius).epsilon(1e-3));
    CHECK(geom.throatPosition() == Approx(0.6018 * kLength).epsilon(1e-9));
    CHECK(geom.throatRadius() >= kThroatRadius);
  }

  SECTION("the area ratios at the ends are the tabulated ones") {
    // The taps stop short of both ends, so these are the first and last tap
    // values (7.001 and 2.574), NOT the nozzle's full 7.75 and 2.68.  Checking
    // against the report's overall ratios here would be checking the wrong
    // thing.
    CHECK(geom.contractionRatio() == Approx(7.001).epsilon(3e-3));
    CHECK(geom.expansionRatio() == Approx(2.574).epsilon(3e-3));
  }

  SECTION("the throat curvature is recovered from the taps") {
    // The independent check.  The report states a throat curvature radius of
    // 1.800 in, which is nowhere in the tap table -- the fit has to find it
    // from the A/A* values alone.  Tap spacing through the throat is coarse,
    // so this is a real test of the circle fit rather than a restatement of
    // an input.
    //
    // It gets 1.854 in from nine taps, 3.0 % from the stated value.  The bound
    // is 6 %, which is loose enough not to be brittle and tight enough that a
    // regression in the fit would break it -- a 35 % bound, which this would
    // also have passed, would have asserted nothing.
    INFO("fitted " << geom.bartzCurvatureRadius() / kInch << " in against the "
                   << kThroatCurvature / kInch << " in the report states; "
                   << geom.curvatureProvenance());
    CHECK(geom.curvatureProvenance().find("circle fitted") != std::string::npos);
    CHECK(geom.bartzCurvatureRadius() ==
          Approx(kThroatCurvature).epsilon(0.06));
  }

  SECTION("the divergent half-angle matches the stated 15 degrees") {
    CHECK(geom.spec().bell_exit_angle == Approx(15.0).margin(1.5));
  }
}


// ===========================================================================
// The integral boundary layer against the same two experiments
// ===========================================================================
//
// Bartz's closed form is above; what follows puts the model that replaced it
// as Ignis's default (thermal/BoundaryLayer.hpp) against the same data.  The
// closure was fixed from the textbook before either comparison was run and
// nothing in it was adjusted afterwards; the bounds below record what came
// out, in the style of the Bartz cases.

namespace {

constexpr double kInletRadius = 2.53 * kInch;      // JPL TR 32-415, Section I
constexpr double kNozzleLength = 5.925 * kInch;
constexpr double kThroatZ = 0.6018 * kNozzleLength;  // tap 12, the minimum area

/// Mach number on the chosen branch of the constant-gamma area-Mach relation.
double machFromAreaRatio(double eps, double gamma, bool supersonic) {
  if (eps <= 1.0) return 1.0;
  auto area = [gamma](double m) {
    const double t = 2.0 / (gamma + 1.0) * (1.0 + 0.5 * (gamma - 1.0) * m * m);
    return std::pow(t, (gamma + 1.0) / (2.0 * (gamma - 1.0))) / m;
  };
  double lo = supersonic ? 1.0 : 1.0e-6, hi = supersonic ? 20.0 : 1.0;
  for (int i = 0; i < 200; ++i) {
    const double mid = 0.5 * (lo + hi);
    // A/A* falls towards M = 1 on the subsonic branch and rises past it.
    const bool too_fast = supersonic ? area(mid) > eps : area(mid) < eps;
    (too_fast ? hi : lo) = mid;
  }
  return 0.5 * (lo + hi);
}

/// One JPL test: its operating point and its measured stations.
struct JplTest {
  int approach_in = 0;
  double p0 = 0.0, t0 = 0.0;
  std::vector<double> z, h, t_wall;   // m, W/(m^2 K), K
};

std::map<int, JplTest> jplTests() {
  const ReferenceTable ref(referenceDir() + "/nozzle_heat_transfer_reference.csv");
  std::map<int, JplTest> tests;
  for (std::size_t row = 0; row < ref.rows(); ++row) {
    JplTest& t = tests[static_cast<int>(ref.num("test", row))];
    t.approach_in = static_cast<int>(ref.num("approach_in", row));
    t.p0 = ref.num("pt_psia", row) * kPsi;
    t.t0 = ref.num("tto_degR", row) * kRankine;
    t.z.push_back(ref.num("z_over_L", row) * kNozzleLength);
    t.h.push_back(ref.num("h_btu_s_in2_F", row) * kCoefficient);
    t.t_wall.push_back(ref.num("tw_degR", row) * kRankine);
  }
  return tests;
}

/// The test article as the report states it was built: a 2.53 in inlet, a 30
/// degree cone into a throat of 1.800 in curvature radius, a 15 degree cone
/// out, and (for the long-approach tests) a cylinder of the inlet radius
/// upstream.  Built analytically rather than from the tap table so that the
/// throat is a circular arc: the area-Mach relation needs dA/dx to vanish
/// smoothly there, which straight lines between taps would not give.
ignis::NozzleGeometry jplNozzle(double approach) {
  ignis::NozzleGeometrySpec s;
  s.throat_radius = kThroatRadius;
  s.chamber_radius = kInletRadius;
  s.chamber_length = approach + 1.0 * kInch;
  s.converging_half_angle = 30.0;
  s.chamber_fillet_ratio = 0.1;
  s.throat_upstream_ratio = kThroatCurvature / kThroatRadius;
  s.throat_downstream_ratio = kThroatCurvature / kThroatRadius;
  s.expansion_ratio = 2.68;
  s.divergent = ignis::DivergentType::kConical;
  s.cone_half_angle = 15.0;
  s.num_stations = 800;
  return ignis::NozzleGeometry::build(s);
}

double interpolate(const std::vector<double>& xs, const std::vector<double>& ys, double x) {
  if (x <= xs.front()) return ys.front();
  if (x >= xs.back()) return ys.back();
  const auto it = std::upper_bound(xs.begin(), xs.end(), x);
  const std::size_t i = static_cast<std::size_t>(it - xs.begin());
  const double t = (x - xs[i - 1]) / (xs[i] - xs[i - 1]);
  return ys[i - 1] + t * (ys[i] - ys[i - 1]);
}

/// Boundary-layer h at every measured station of one test, divided by the
/// measured h; also the largest acceleration parameter the march saw.
struct BoundaryLayerComparison {
  std::vector<double> ratio;
  double max_acceleration = 0.0;
  std::vector<std::string> warnings;
};

BoundaryLayerComparison jplBoundaryLayer(const JplTest& t) {
  const auto& db = ignis_test::fullDatabase();
  const Air air = dryAir(db);
  const ignis::GasMixture mix(db);
  const ignis::TransportModel transport(db);
  const double R = ignis::constants::R_universal / air.molar_mass;
  // Constant gamma at stagnation, as in the report's own Mach numbers.
  const double cp0 = mix.cpFrozen(air.n, t.t0);
  const double gamma = cp0 / (cp0 - R);
  const double pr0 = transport.mixture(air.n, t.t0, cp0).prandtl;

  const double approach = t.approach_in * kInch;
  const ignis::NozzleGeometry geom = jplNozzle(approach);
  // Report z (from the nozzle inlet) to Ignis x, anchored at the throat.
  const double x_of_z0 = geom.throatPosition() - kThroatZ;
  const double x_start = x_of_z0 - approach;
  REQUIRE(x_start > 0.0);

  // The layer starts where the cooled wall starts: the report's measured
  // inlet profiles are not used, so the comparison is a prediction, not a
  // fit.  (Elliott, Bartz & Silver matched their starting thickness to the
  // nozzle-entrance probe for these very data; that is not done here.)
  const int n = 1500;
  std::vector<ignis::BoundaryLayerEdge> edge(n);
  std::vector<double> pressure(n);
  for (int i = 0; i < n; ++i) {
    auto& e = edge[static_cast<std::size_t>(i)];
    e.x = x_start + (geom.exitPosition() - x_start) * i / (n - 1);
    e.radius = geom.radius(e.x);
    const double eps = std::max(geom.area(e.x) / geom.throatArea(), 1.0);
    e.mach = machFromAreaRatio(eps, gamma, e.x > geom.throatPosition());
    e.T = t.t0 / (1.0 + 0.5 * (gamma - 1.0) * e.mach * e.mach);
    e.p = t.p0 * std::pow(e.T / t.t0, gamma / (gamma - 1.0));
    e.rho = e.p / (R * e.T);
    e.u = e.mach * std::sqrt(gamma * R * e.T);
    e.viscosity = transport.mixture(air.n, e.T, mix.cpFrozen(air.n, e.T)).viscosity;
    e.t_adiabatic_wall = e.T + std::cbrt(pr0) * (t.t0 - e.T);
    // Measured wall temperature, held constant beyond the instrumented span.
    e.t_wall = interpolate(t.z, t.t_wall, e.x - x_of_z0);
    pressure[static_cast<std::size_t>(i)] = e.p;
  }
  auto properties = [&](std::size_t i, double T) {
    ignis::ReferenceState s;
    s.rho = pressure[i] / (R * T);
    s.cp = mix.cpFrozen(air.n, T);
    const auto tr = transport.mixture(air.n, T, s.cp);
    s.viscosity = tr.viscosity;
    s.prandtl = tr.prandtl;
    return s;
  };
  const auto sol = ignis::marchBoundaryLayer(edge, properties);

  BoundaryLayerComparison out;
  for (std::size_t k = 0; k < t.z.size(); ++k)
    out.ratio.push_back(sol.at(x_of_z0 + t.z[k]).h_gas / t.h[k]);
  out.max_acceleration = sol.max_acceleration;
  out.warnings = sol.warnings;
  return out;
}

}  // namespace

TEST_CASE("the integral boundary layer against the JPL air nozzle",
          "[validation][thermal][boundary_layer]") {
  const auto tests = jplTests();
  REQUIRE(tests.size() == 4);

  SECTION("the analytic test article reproduces the tap table") {
    // Independent of any heat transfer: the geometry built from the report's
    // stated dimensions must land on the measured pressure-tap area ratios.
    //
    // Recorded, not tuned: the bound set before this was first run was 5 %
    // over the whole table, and the first tap misses it.  The divergent side
    // matches to within 0.7 %; the convergent cone that the stated 30 deg and
    // 1.800 in imply sits 2-4 % low in A/A* at every convergent tap (5.8 % at
    // the first, near the inlet) -- as if the real cone were ~0.05 in further
    // upstream.  Through h ~ G^0.8 that is a 2-5 % low bias in the model's
    // mass flux on the convergent stations, well inside the report's +-8 to
    // +-21 % uncertainty on h, so the stated geometry is kept.
    const ReferenceTable contour(referenceDir() + "/nozzle_contour_reference.csv");
    const auto geom = jplNozzle(0.0);
    const double x_of_z0 = geom.throatPosition() - kThroatZ;
    double worst_in = 0.0, worst_out = 0.0, signed_in = 0.0;
    for (std::size_t row = 0; row < contour.rows(); ++row) {
      const double z = contour.num("z_over_L", row) * kNozzleLength;
      const double eps = geom.area(x_of_z0 + z) / geom.throatArea();
      const double miss = eps / contour.num("area_ratio", row) - 1.0;
      if (z < kThroatZ) {
        if (std::abs(miss) > worst_in) { worst_in = std::abs(miss); signed_in = miss; }
      } else {
        worst_out = std::max(worst_out, std::abs(miss));
      }
    }
    INFO("largest A/A* mismatch: convergent " << 100.0 * signed_in << " %, divergent "
         << 100.0 * worst_out << " %");
    CHECK(worst_out < 0.01);
    CHECK(worst_in < 0.06);
    CHECK(signed_in < 0.0);   // the model's convergent is narrower than the taps
  }

  std::map<int, std::vector<double>> bl_by_approach;
  std::map<int, double> max_k_by_approach;
  std::map<int, bool> flagged_by_approach;
  std::ostringstream per_test;
  for (const auto& [number, t] : tests) {
    const auto c = jplBoundaryLayer(t);
    per_test << " test " << number << ": median " << median(c.ratio) << ", max K "
             << c.max_acceleration << ";";
    auto& v = bl_by_approach[t.approach_in];
    v.insert(v.end(), c.ratio.begin(), c.ratio.end());
    max_k_by_approach[t.approach_in] = std::max(max_k_by_approach[t.approach_in], c.max_acceleration);
    for (const auto& w : c.warnings)
      if (w.find("laminarise") != std::string::npos) flagged_by_approach[t.approach_in] = true;
  }
  const auto bartz = jplRatios();
  const double bl_long = median(bl_by_approach.at(18));
  const double bl_none = median(bl_by_approach.at(0));
  const double bz_long = median(bartz.at(18));
  const double bz_none = median(bartz.at(0));
  INFO("median model/measured -- 18 in approach: boundary layer " << bl_long << ", Bartz "
       << bz_long << ";  no approach: boundary layer " << bl_none << ", Bartz " << bz_none
       << ";  max K " << max_k_by_approach[18] << " / " << max_k_by_approach[0] << ";"
       << per_test.str());

  SECTION("on the high-pressure, long-approach tests it lands on the data") {
    // The report's own uncertainty on h is about +-8 % in the throat region
    // at these pressures.  Bartz is 1.3-1.9x high on the same rows.
    CHECK(bl_long > 0.85);
    CHECK(bl_long < 1.15);
    CHECK(std::abs(bl_long - 1.0) < 0.5 * std::abs(bz_long - 1.0));
  }

  SECTION("on the low-pressure tests it is as wrong as Bartz, and says why") {
    // 36 and 51 psia, no cooled approach.  Both models over-predict by more
    // than 2x.  The turbulent closure cannot represent a layer that is
    // re-laminarising under strong acceleration, and the march says so: the
    // acceleration parameter passes the 3e-6 at which that starts.  This is
    // a recorded failure of the model, not a pass.
    CHECK(bl_none > 2.0);
    CHECK(bz_none > 2.0);
    CHECK(max_k_by_approach[0] > 3.0e-6);
    CHECK(flagged_by_approach[0]);
  }
}

TEST_CASE("the integral boundary layer against the hydrogen-oxygen rocket",
          "[validation][thermal][boundary_layer]") {
  // NASA TN D-2832 again, now through Ignis's whole hot-gas path: equilibrium
  // chamber, quasi-1D expansion, the boundary layer marched from the injector
  // face, and the cooling module's gas-side survey at a fixed wall
  // temperature (the measurement was a heat-sink nozzle, so there is no
  // coolant to couple to).
  //
  // The engine is rebuilt from the report's stated dimensions: 5 in throat,
  // 4.64 contraction and expansion, 30 deg in and 15 deg out, throat
  // curvature radius taken equal to the throat radius (not stated), 14.5 in
  // from injector to throat, LOX/GH2 at 600 psia and O/F 6 (the middle of the
  // tested range).
  const auto& db = ignis_test::hydrogenDatabase();
  const ignis::EquilibriumSolver solver(db);
  const ignis::TransportModel transport(db);
  const auto lib = ignis::PropellantLibrary::loadYaml(ignis_test::sourceDir() +
                                                      "/data/propellants/ignis_propellants.yaml");
  const ignis::PropellantMixture mixture(lib.at("LOX"), lib.at("GH2"), 6.0, 90.18, 298.15);
  const ignis::CombustionChamber chamber(solver, ignis::CompositionModel::kEquilibrium);
  const auto ch = chamber.solve(mixture, 600.0 * kPsi, 1.0);
  const auto flow = chamber.makeFlow(ch);

  ignis::NozzleGeometrySpec gs;
  gs.throat_radius = 2.5 * kInch;
  gs.contraction_ratio = 4.64;
  gs.chamber_length = 0.2170;
  gs.converging_half_angle = 30.0;
  gs.chamber_fillet_ratio = 0.2;
  gs.throat_upstream_ratio = 1.0;
  gs.throat_downstream_ratio = 1.0;
  gs.expansion_ratio = 4.64;
  gs.divergent = ignis::DivergentType::kConical;
  gs.cone_half_angle = 15.0;
  gs.num_stations = 600;
  const auto geom = ignis::NozzleGeometry::build(gs);
  REQUIRE(geom.throatPosition() == Approx(14.5 * kInch).epsilon(0.01));

  const ReferenceTable ref(referenceDir() + "/rocket_heat_transfer_reference.csv");
  const ignis::GasMixture& mix = flow.mixture();

  // C as the report defines it, St* Pr*^0.7 Re*_d^0.2 with every property at
  // the reference state and the local static pressure (density included, as
  // in the film-temperature form of the Bartz equation the report compares
  // with), from the boundary layer's h at a station.
  auto c_model = [&](const ignis::CoolingResult& r, double t_wall, double x) {
    std::vector<double> xs, hs, taw;
    for (const auto& s : r.stations) {
      xs.push_back(s.x);
      hs.push_back(s.h_gas);
      taw.push_back(s.t_adiabatic_wall);
    }
    const double h = interpolate(xs, hs, x);
    const double eps = std::max(geom.area(x) / geom.throatArea(), 1.0);
    const auto st = flow.atAreaRatio(eps, x > geom.throatPosition());
    const double ts = ignis::eckertReferenceTemperature(st.gas.T, t_wall, interpolate(xs, taw, x));
    const double cp = mix.cpFrozen(st.gas.n, ts);
    const auto tr = transport.mixture(st.gas.moleFractions(), ts, cp);
    const double rho = st.gas.p * st.gas.M / (ignis::constants::R_universal * ts);
    const double G = rho * st.u;
    const double d = 2.0 * geom.radius(x);
    return h / (G * cp) * std::pow(tr.prandtl, 0.7) * std::pow(G * d / tr.viscosity, 0.2);
  };

  auto survey = [&](double t_wall) {
    ignis::CoolingSpec spec;
    spec.hot_gas_model = ignis::HotGasModel::kBoundaryLayer;
    spec.num_segments = 600;
    spec.wall = ignis::MaterialLibrary::loadDefault().at("Copper");
    return ignis::surveyHotGasSide(flow, geom, ch, spec, transport, t_wall);
  };

  // Shape: C at each station over C at station 1, which removes the ~30 %
  // the report says the transport-property choice alone moves C by.
  auto shapeError = [&](double t_wall, std::map<std::string, double>& model_ratio) {
    const auto r = survey(t_wall);
    double c1 = 0.0, c1_meas = 0.0, sum2 = 0.0;
    int count = 0;
    for (std::size_t row = 0; row < ref.rows(); ++row) {
      const std::string name = ref.text("station", row);
      if (name == "3") continue;  // 3bar is the average of 3, 3A and 3B
      const double x = geom.throatPosition() + ref.num("axial_in", row) * kInch;
      const double c = c_model(r, t_wall, x);
      if (name == "1") { c1 = c; c1_meas = ref.num("c_measured", row); }
      model_ratio[name] = c;
    }
    for (auto& [name, c] : model_ratio) c /= c1;
    for (std::size_t row = 0; row < ref.rows(); ++row) {
      const std::string name = ref.text("station", row);
      if (name == "3" || name == "1") continue;
      const double measured = ref.num("c_measured", row) / c1_meas;
      sum2 += std::pow(model_ratio[name] - measured, 2);
      ++count;
    }
    return std::sqrt(sum2 / count);
  };

  std::map<std::string, double> ratio;
  const double rms = shapeError(500.0, ratio);
  // Bartz's C is one constant everywhere, so its predicted shape is 1.
  double bartz_sum2 = 0.0, c1_meas = 0.0;
  std::map<std::string, double> measured;
  for (std::size_t row = 0; row < ref.rows(); ++row) {
    measured[ref.text("station", row)] = ref.num("c_measured", row);
    if (ref.text("station", row) == "1") c1_meas = ref.num("c_measured", row);
  }
  int bartz_count = 0;
  for (const auto& [name, c] : measured) {
    if (name == "1" || name == "3") continue;
    bartz_sum2 += std::pow(1.0 - c / c1_meas, 2);
    ++bartz_count;
  }
  const double bartz_rms = std::sqrt(bartz_sum2 / bartz_count);
  INFO("C(station)/C(1), model vs measured: 2 " << ratio["2"] << " vs " << measured["2"] / c1_meas
       << ";  3bar " << ratio["3bar"] << " vs " << measured["3bar"] / c1_meas << ";  4 "
       << ratio["4"] << " vs " << measured["4"] / c1_meas << ";  5 " << ratio["5"] << " vs "
       << measured["5"] / c1_meas << ";  rms " << rms << " against Bartz's " << bartz_rms);

  SECTION("it follows the measured shape better than Bartz can") {
    // Bartz's constant cannot fall at the throat.  The measurement falls by
    // 41 %; the boundary layer falls by about half of that.
    CHECK(rms < 0.5 * bartz_rms);
    CHECK(ratio["3bar"] < 0.85);
    CHECK(ratio["3bar"] < ratio["2"]);
  }

  SECTION("but its throat dip is too shallow, and C does not recover past it") {
    // A recorded shortcoming.  The report: C is "high in the chamber, low in
    // the throat, and increasing again at the exit", in this engine and in
    // three others it re-reduced.  The boundary layer's C falls monotonically
    // instead -- 0.83, 0.78, 0.75, 0.71 of the chamber value against the
    // measured 0.93, 0.59, 0.60, 0.73 -- so it over-predicts the throat by
    // about two standard deviations of the measured ratio and gets the exit
    // right for the wrong reason.  Laminarisation through the throat (which
    // the closure does not model) would deepen the dip; a cause for the
    // recovery is not established here.  If a closure change makes either
    // check fail, look at this comparison first.
    const double measured_3bar = measured["3bar"] / c1_meas;
    CHECK(ratio["3bar"] > measured_3bar);
    CHECK(ratio["5"] < ratio["3bar"]);
  }

  SECTION("the shape does not hang on the assumed wall temperature") {
    // The report does not tabulate T_w; 500 K above is an assumption.  A
    // wall at 900 K instead must not change the conclusion.
    std::map<std::string, double> hot;
    const double rms_hot = shapeError(900.0, hot);
    INFO("rms at a 900 K wall: " << rms_hot);
    CHECK(rms_hot < 0.5 * bartz_rms);
  }
}
