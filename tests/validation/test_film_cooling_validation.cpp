// SPDX-License-Identifier: MIT
/// \file test_film_cooling_validation.cpp
/// \brief The film-cooling correlation against the measurements it came from.
///
/// Hatch & Papell published their correlation with the helium data it was
/// fitted to (NASA TN D-130, Table I).  Re-deriving their effectiveness at
/// every station, from their stated gas-side coefficient, Ignis's air
/// properties and a published helium conductivity, and comparing it with the
/// wall temperatures they measured checks two things: that eq. (12) and its
/// velocity function were transcribed correctly, and that the inputs are read
/// the way the authors meant them.  It does NOT validate the correlation for a
/// rocket: these films are a few hundred kelvin of air and helium.  See
/// FilmCooling.hpp for that caveat.
#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "TestHelpers.hpp"
#include "ignis/core/Constants.hpp"
#include "ignis/thermal/FilmCooling.hpp"
#include "ignis/thermo/GasMixture.hpp"
#include "ignis/thermo/Transport.hpp"

using ignis_test::ReferenceTable;
using ignis_test::referenceDir;

namespace {

constexpr double kFoot = 0.3048;
constexpr double kPound = 0.45359237;
constexpr double kRankine = 5.0 / 9.0;
constexpr double kInch = 0.0254;
constexpr double kSlotLength = 0.667 * kFoot;   // L, Table I
constexpr double kDuct = 8.0 * kInch;           // 8 in square: D_h = 8 in, A = 64 in^2

/// Helium, an ideal monatomic gas: c_p = (5/2) R / M.
constexpr double kHeliumMolarMass = 4.002602e-3;
const double kHeliumCp = 2.5 * ignis::constants::R_universal / kHeliumMolarMass;

/// Helium thermal conductivity, W/(m K): H. Petersen, "The properties of
/// helium", Risoe Report 224 (1970), eq. (7-1), p in bar.  Standard deviation
/// about 1 % at 273 K.
double heliumConductivity(double T, double p_bar) {
  return 2.682e-3 * (1.0 + 1.123e-3 * p_bar) * std::pow(T, 0.71 * (1.0 - 2.0e-4 * p_bar));
}

const char* kStations[] = {"0.132", "0.209", "0.376", "0.543", "0.710", "0.876",
                           "1.042", "1.293", "1.501", "1.751", "2.000", "2.822"};

struct Point {
  double error;           ///< t_w predicted / measured - 1, absolute temperatures
  double eta_measured;
  double velocity_ratio;  ///< V_g / V_c
  double slot;            ///< in
};

std::vector<Point> compare() {
  const ReferenceTable ref(referenceDir() + "/film_cooling_reference.csv");
  const auto& db = ignis_test::fullDatabase();
  const ignis::GasMixture mix(db);
  const ignis::TransportModel transport(db);
  Eigen::VectorXd X = Eigen::VectorXd::Zero(static_cast<int>(db.size()));
  X(db.index("N2")) = 0.78084;
  X(db.index("O2")) = 0.20946;
  X(db.index("AR")) = 0.00934;
  const Eigen::VectorXd n = mix.fromMoleFractions(X);
  const double r_air = ignis::constants::R_universal / mix.molarMass(n);

  std::vector<Point> out;
  for (std::size_t row = 0; row < ref.rows(); ++row) {
    const double S = ref.num("slot_in", row) * kInch;
    const double wc = ref.num("wc_lb_s", row) * kPound;
    const double wg = ref.num("wg_lb_s", row) * kPound;
    const double vc = ref.num("vc_ft_s", row) * kFoot;
    const double vg = ref.num("vg_ft_s", row) * kFoot;
    const double t_ad = ref.num("tad_R", row) * kRankine;
    const double t_c = ref.num("tc_R", row) * kRankine;

    // The report's T_ad = t_g + Pr^(1/3) V^2 / (2 c_p): recover t_g.
    double t_g = t_ad;
    for (int it = 0; it < 6; ++it) {
      const double cp = mix.cpFrozen(n, t_g);
      const double pr = transport.mixture(X, t_g, cp).prandtl;
      t_g = t_ad - std::cbrt(pr) * vg * vg / (2.0 * cp);
    }
    // Static pressure from the duct's continuity, which the table implies.
    const double p = wg / (kDuct * kDuct * vg) * r_air * t_g;

    // h of the report's assumption 6, every property at (t_g + t_c)/2.
    const double t_f = 0.5 * (t_g + t_c);
    const double cp_f = mix.cpFrozen(n, t_f);
    const auto tr = transport.mixture(X, t_f, cp_f);
    const double h = ignis::hatchPapellGasCoefficient(p / (r_air * t_f), vg, kDuct, tr.viscosity,
                                                      tr.conductivity, tr.prandtl);

    // The helium density the authors' own V_c implies.
    ignis::FilmSlot slot;
    slot.mass_flow = wc;
    slot.slot_height = S;
    slot.perimeter = kSlotLength;
    slot.temperature = t_c;
    slot.density = wc / (vc * kSlotLength * S);
    slot.cp = kHeliumCp;
    slot.conductivity = heliumConductivity(t_c, p * 1e-5);
    slot.gas_velocity = vg;
    // Unclamped, on purpose: this is the correlation exactly as published,
    // including where its velocity function is extrapolated.
    const double factor = std::pow(S * vg / slot.diffusivity(), 0.125) *
                          ignis::hatchPapellVelocityFunction(vg / vc);

    for (const char* station : kStations) {
      const double x = std::stod(station) * kFoot;
      const double eta = ignis::hatchPapellEffectiveness(h * kSlotLength * x / (wc * kHeliumCp),
                                                         factor);
      const double t_w_model = t_ad - eta * (t_ad - t_c);
      const double t_w = ref.num(std::string("tw_") + station, row) * kRankine;
      out.push_back({t_w_model / t_w - 1.0, (t_ad - t_w) / (t_ad - t_c), vg / vc,
                     ref.num("slot_in", row)});
    }
  }
  return out;
}

double meanAbs(const std::vector<double>& v) {
  double s = 0.0;
  for (double e : v) s += std::abs(e);
  return v.empty() ? 0.0 : s / static_cast<double>(v.size());
}

}  // namespace

TEST_CASE("Hatch and Papell reproduce their own helium data",
          "[validation][thermal][film]") {
  const auto points = compare();
  REQUIRE(points.size() == 27 * 12);

  // The report claims about +-5 % in wall temperature for effectiveness 0.2
  // to 1, within its velocity-ratio range.  Sort the points accordingly.
  std::vector<double> in_range, below_range;
  for (const auto& p : points) {
    if (p.eta_measured < 0.2) continue;
    if (p.velocity_ratio >= ignis::kHatchPapellMinVelocityRatio &&
        p.velocity_ratio <= ignis::kHatchPapellMaxVelocityRatio)
      in_range.push_back(p.error);
    else
      below_range.push_back(p.error);
  }
  std::vector<double> sorted = in_range;
  std::sort(sorted.begin(), sorted.end());
  const double median = sorted[sorted.size() / 2];
  int within_five = 0;
  for (double e : in_range) within_five += std::abs(e) <= 0.05 ? 1 : 0;
  INFO("in range: " << in_range.size() << " points, mean |error| " << meanAbs(in_range)
       << ", median " << median << ", " << within_five << " within 5 %;  below the range: "
       << below_range.size() << " points, mean |error| " << meanAbs(below_range));

  SECTION("inside the stated range the transcription matches the paper's claim") {
    // Recorded, not tuned: 257 points, mean |error| 4.0 %, median +1.4 %,
    // 65 % of them inside +-5 %.  "Roughly +-5 %" is the report's word.
    CHECK(in_range.size() > 250);
    CHECK(meanAbs(in_range) < 0.05);
    CHECK(std::abs(median) < 0.02);
    CHECK(within_five > static_cast<int>(0.6 * static_cast<double>(in_range.size())));
  }

  SECTION("below it, where eq. (11) is extrapolated, the correlation fails") {
    // Three 1/8 in runs with a helium jet 2.3-3.8 times faster than the air.
    // Their velocity function reaches f = 300, and the predicted wall runs
    // 40 % hot on average (36 points).  This is why Ignis clamps the ratio and
    // warns.
    CHECK(below_range.size() >= 30);
    CHECK(meanAbs(below_range) > 0.2);
  }
}
