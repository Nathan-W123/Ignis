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
#include <string>
#include <vector>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "TestHelpers.hpp"
#include "ignis/core/Constants.hpp"
#include "ignis/thermo/GasMixture.hpp"
#include "ignis/thermo/Transport.hpp"
#include "ignis/thermal/HeatTransfer.hpp"

using ignis_test::ReferenceTable;
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

  SECTION("Bartz over-predicts, and by how much depends on the inlet boundary layer") {
    // These bounds record a measurement, not a target.  Bartz assumes a
    // turbulent boundary layer that has grown from the injector; the report's
    // 18-in. cooled approach supplies one, and its 0-in. cases deliberately do
    // not.  The correlation is correspondingly further out on the thin-layer
    // cases, which is the expected failure mode and not an implementation
    // fault.  If these ranges move, the correlation, the transcription or the
    // evaluation changed -- all three are worth knowing about.
    const double developed = median(ratio_by_approach.at(18));
    const double thin = median(ratio_by_approach.at(0));

    CHECK(developed > 1.3);
    CHECK(developed < 1.9);
    CHECK(thin > 2.2);
    CHECK(thin > developed);
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
    const double air_developed = median(jplRatios().at(18));
    INFO("rocket " << rocket << " vs air " << air_developed);
    CHECK(std::abs(rocket - air_developed) / rocket < 0.25);
  }
}
