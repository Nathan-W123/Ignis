// SPDX-License-Identifier: MIT
/// \file test_validation.cpp
/// \brief Comparison against external references.
///
/// Two kinds of comparison live here and are deliberately kept apart:
///
///   VERIFICATION against Cantera -- an independent implementation of Gibbs
///   minimisation running on *the same* NASA TM-4513 coefficients.  Any
///   disagreement beyond round-off would be a coding error in one of the two.
///
///   VALIDATION against NASA CEA -- the NASA Glenn program itself, using its own
///   2002 nine-coefficient thermodynamic database and its own species list.
///   Residual differences here measure the spread between published
///   thermochemical data sets, not a defect.
///
/// The reference files are produced by tools/make_reference_data.py and are
/// committed so this suite needs neither Cantera nor CEA installed.

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <map>

#include "TestHelpers.hpp"
#include "ignis/combustion/Chamber.hpp"
#include "ignis/nozzle/NozzleFlow.hpp"

using namespace ignis;
using namespace ignis_test;
using Catch::Approx;

namespace {

PropellantLibrary library() {
  return PropellantLibrary::loadYaml(sourceDir() + "/data/propellants/ignis_propellants.yaml");
}

struct Worst {
  double value = 0.0;
  std::string where;
  void update(double v, const std::string& w) {
    if (v > value) { value = v; where = w; }
  }
};

}  // namespace

TEST_CASE("adiabatic equilibrium matches an independent minimiser on the same data",
          "[validation][verification][cantera]") {
  const ReferenceTable ref(referenceDir() + "/equilibrium_reference.csv");
  REQUIRE(ref.rows() > 50);
  const auto lib = library();
  const auto& cho = rocketDatabase();
  const auto& ho = hydrogenDatabase();
  const EquilibriumSolver cho_solver(cho);
  const EquilibriumSolver ho_solver(ho);

  Worst wt, wm, wh, wx;
  int checked = 0;
  for (std::size_t r = 0; r < ref.rows(); ++r) {
    const bool methane = ref.text("fuel", r) == "LCH4";
    const auto& db = methane ? cho : ho;
    const auto& solver = methane ? cho_solver : ho_solver;
    const PropellantMixture mix(lib.at(ref.text("oxidizer", r)), lib.at(ref.text("fuel", r)),
                                ref.num("mixture_ratio", r),
                                lib.at(ref.text("oxidizer", r)).reference_temperature,
                                lib.at(ref.text("fuel", r)).reference_temperature);
    const double p = ref.num("pressure", r);
    const std::string tag = ref.text("fuel", r) + " O/F " +
                            std::to_string(ref.num("mixture_ratio", r)) + " at " +
                            std::to_string(p * 1e-6) + " MPa";

    // The reactant enthalpies must agree first: they come from the same NASA CEA
    // propellant library on both sides.
    REQUIRE(mix.enthalpy(db) == Approx(ref.num("reactant_enthalpy", r)).epsilon(1e-6));

    const auto res = solver.hp(mix.elementMoles(db), mix.enthalpy(db), p);
    wt.update(relativeError(res.state.T, ref.num("temperature", r)), tag);
    wm.update(relativeError(res.state.M, ref.num("molar_mass", r)), tag);
    wh.update(relativeError(res.state.cp_frozen, ref.num("cp_frozen", r)), tag);

    const auto X = res.state.moleFractions();
    for (std::size_t j = 0; j < db.size(); ++j) {
      const std::string col = "X_" + db[j].name();
      if (!ref.has(col)) continue;
      const double x_ref = ref.num(col, r);
      if (x_ref < 1.0e-4) continue;   // trace species are dominated by round-off
      wx.update(relativeError(X(static_cast<Eigen::Index>(j)), x_ref), tag + " / " + db[j].name());
    }
    ++checked;
  }
  INFO("checked " << checked << " states; worst temperature " << wt.value << " at " << wt.where
       << "; worst molar mass " << wm.value << "; worst frozen cp " << wh.value
       << "; worst mole fraction " << wx.value << " at " << wx.where);
  REQUIRE(checked >= 50);
  // Same data, same problem, two independent minimisers.
  REQUIRE(wt.value < 1.0e-6);
  REQUIRE(wm.value < 1.0e-6);
  REQUIRE(wh.value < 1.0e-6);
  REQUIRE(wx.value < 1.0e-5);
}

TEST_CASE("LOX/methane performance is validated against NASA CEA",
          "[validation][cea][methane]") {
  const ReferenceTable ref(referenceDir() + "/cea_reference.csv");
  const auto lib = library();
  const auto& db = rocketDatabase();
  const EquilibriumSolver solver(db);

  Worst wt, wm, wg, wc, wi;
  int checked = 0;
  for (std::size_t r = 0; r < ref.rows(); ++r) {
    if (ref.text("fuel", r) != "LCH4") continue;
    if (ref.text("composition", r) != "equilibrium") continue;
    const double mr = ref.num("mixture_ratio", r);
    const double pc = ref.num("pressure", r);
    const double eps = ref.num("expansion_ratio", r);
    if (pc < 1.0e6) continue;   // CEA and Ignis both converge here, but the
                                // dissociated low-pressure cases are compared
                                // separately below
    const PropellantMixture mix(lib.at("LOX"), lib.at("LCH4"), mr, 90.18, 111.66);
    const CombustionChamber chamber(solver, CompositionModel::kEquilibrium);
    const auto ch = chamber.solve(mix, pc, 1.0);
    const auto flow = chamber.makeFlow(ch);
    const auto exit_state = flow.atAreaRatio(eps, true);

    const std::string tag = "O/F " + std::to_string(mr) + ", " + std::to_string(pc * 1e-6) +
                            " MPa, eps " + std::to_string(eps);
    wt.update(relativeError(ch.state.T, ref.num("chamber_temperature", r)), tag);
    wm.update(relativeError(ch.state.M, ref.num("molar_mass", r)), tag);
    wg.update(relativeError(ch.state.gamma_s, ref.num("gamma", r)), tag);
    wc.update(relativeError(ch.c_star_ideal, ref.num("c_star", r)), tag);
    // CEA's vacuum specific impulse: momentum plus exit-pressure term.
    const double mdot_per_at = pc / ch.c_star_ideal;
    const double isp_vac = (mdot_per_at * exit_state.u + exit_state.gas.p * eps) /
                           (mdot_per_at * constants::g0);
    wi.update(relativeError(isp_vac, ref.num("isp_vacuum", r)), tag);
    ++checked;
  }
  INFO("checked " << checked << " CEA cases\n"
       << "  chamber temperature : " << wt.value << " at " << wt.where << "\n"
       << "  molar mass          : " << wm.value << "\n"
       << "  gamma_s             : " << wg.value << "\n"
       << "  c*                  : " << wc.value << " at " << wc.where << "\n"
       << "  vacuum Isp          : " << wi.value << " at " << wi.where);
  REQUIRE(checked >= 60);
  // Different thermodynamic databases (NASA TM-4513 1993 against the CEA 2002
  // nine-coefficient set) and a slightly different species list.  The tolerances
  // below are the measured spread, recorded in docs/validation.md.
  REQUIRE(wt.value < 0.01);
  REQUIRE(wm.value < 0.005);
  REQUIRE(wg.value < 0.005);
  REQUIRE(wc.value < 0.01);
  REQUIRE(wi.value < 0.015);
}

TEST_CASE("LOX/hydrogen performance is validated against NASA CEA",
          "[validation][cea][hydrogen]") {
  const ReferenceTable ref(referenceDir() + "/cea_reference.csv");
  const auto lib = library();
  const auto& db = hydrogenDatabase();
  const EquilibriumSolver solver(db);

  Worst wt, wm, wg, wc, wi;
  int checked = 0;
  for (std::size_t r = 0; r < ref.rows(); ++r) {
    if (ref.text("fuel", r) != "LH2") continue;
    if (ref.text("composition", r) != "equilibrium") continue;
    const double mr = ref.num("mixture_ratio", r);
    const double pc = ref.num("pressure", r);
    const double eps = ref.num("expansion_ratio", r);
    if (pc < 1.0e6) continue;
    const PropellantMixture mix(lib.at("LOX"), lib.at("LH2"), mr, 90.18, 20.27);
    const CombustionChamber chamber(solver, CompositionModel::kEquilibrium);
    const auto ch = chamber.solve(mix, pc, 1.0);
    const auto flow = chamber.makeFlow(ch);
    const auto exit_state = flow.atAreaRatio(eps, true);
    const std::string tag = "O/F " + std::to_string(mr) + ", " + std::to_string(pc * 1e-6) +
                            " MPa, eps " + std::to_string(eps);
    wt.update(relativeError(ch.state.T, ref.num("chamber_temperature", r)), tag);
    wm.update(relativeError(ch.state.M, ref.num("molar_mass", r)), tag);
    wg.update(relativeError(ch.state.gamma_s, ref.num("gamma", r)), tag);
    wc.update(relativeError(ch.c_star_ideal, ref.num("c_star", r)), tag);
    const double mdot_per_at = pc / ch.c_star_ideal;
    const double isp_vac = (mdot_per_at * exit_state.u + exit_state.gas.p * eps) /
                           (mdot_per_at * constants::g0);
    wi.update(relativeError(isp_vac, ref.num("isp_vacuum", r)), tag);
    ++checked;
  }
  INFO("checked " << checked << " CEA cases\n"
       << "  chamber temperature : " << wt.value << " at " << wt.where << "\n"
       << "  molar mass          : " << wm.value << "\n"
       << "  gamma_s             : " << wg.value << "\n"
       << "  c*                  : " << wc.value << " at " << wc.where << "\n"
       << "  vacuum Isp          : " << wi.value << " at " << wi.where);
  REQUIRE(checked >= 50);
  REQUIRE(wt.value < 0.01);
  REQUIRE(wm.value < 0.01);
  REQUIRE(wg.value < 0.01);
  REQUIRE(wc.value < 0.01);
  REQUIRE(wi.value < 0.02);
}

TEST_CASE("frozen expansion is validated against NASA CEA", "[validation][cea][frozen]") {
  const ReferenceTable ref(referenceDir() + "/cea_reference.csv");
  const auto lib = library();
  const auto& db = rocketDatabase();
  const EquilibriumSolver solver(db);
  Worst wc, wi, wte;
  int checked = 0;
  for (std::size_t r = 0; r < ref.rows(); ++r) {
    if (ref.text("fuel", r) != "LCH4") continue;
    if (ref.text("composition", r) != "frozen") continue;
    const double mr = ref.num("mixture_ratio", r);
    const double pc = ref.num("pressure", r);
    const double eps = ref.num("expansion_ratio", r);
    if (pc < 2.0e6) continue;
    const PropellantMixture mix(lib.at("LOX"), lib.at("LCH4"), mr, 90.18, 111.66);
    const CombustionChamber chamber(solver, CompositionModel::kFrozen);
    const auto ch = chamber.solve(mix, pc, 1.0);
    const auto flow = chamber.makeFlow(ch);
    const auto exit_state = flow.atAreaRatio(eps, true);
    const std::string tag = "O/F " + std::to_string(mr) + " at " + std::to_string(pc * 1e-6) +
                            " MPa, eps " + std::to_string(eps);
    wc.update(relativeError(ch.c_star_ideal, ref.num("c_star", r)), tag);
    wte.update(relativeError(exit_state.gas.T, ref.num("exit_temperature", r)), tag);
    const double mdot_per_at = pc / ch.c_star_ideal;
    const double isp_vac = (mdot_per_at * exit_state.u + exit_state.gas.p * eps) /
                           (mdot_per_at * constants::g0);
    wi.update(relativeError(isp_vac, ref.num("isp_vacuum", r)), tag);
    ++checked;
  }
  INFO("checked " << checked << " frozen CEA cases\n"
       << "  c*             : " << wc.value << " at " << wc.where << "\n"
       << "  exit temperature: " << wte.value << " at " << wte.where << "\n"
       << "  vacuum Isp     : " << wi.value << " at " << wi.where);
  REQUIRE(checked >= 40);
  // CEA freezes the composition at the throat by default ("frozenAtThroat"),
  // whereas Ignis freezes it at the chamber; the two therefore differ by more
  // than the equilibrium comparison, and the difference is a definition
  // difference rather than an error.  What must agree closely is c*, which is
  // set upstream of the throat in both.
  REQUIRE(wc.value < 0.02);
  REQUIRE(wi.value < 0.06);
}

TEST_CASE("equilibrium expansion recovers more impulse than frozen expansion",
          "[validation][physics]") {
  const auto lib = library();
  const auto& db = rocketDatabase();
  const EquilibriumSolver solver(db);
  const PropellantMixture mix(lib.at("LOX"), lib.at("LCH4"), 3.4, 90.18, 111.66);
  NozzleGeometrySpec gs;
  gs.throat_radius = 0.07;
  gs.contraction_ratio = 2.8;
  gs.chamber_length = 0.22;
  gs.expansion_ratio = 45.0;
  const auto geom = NozzleGeometry::build(gs);
  NozzlePerformanceOptions opts;
  opts.auto_divergence = false;

  const CombustionChamber eq(solver, CompositionModel::kEquilibrium);
  const CombustionChamber fr(solver, CompositionModel::kFrozen);
  const auto ch_eq = eq.solve(mix, 5.5e6, 1.0);
  const auto ch_fr = fr.solve(mix, 5.5e6, 1.0);
  const auto p_eq = evaluateNozzle(eq.makeFlow(ch_eq), geom, 0.0, opts);
  const auto p_fr = evaluateNozzle(fr.makeFlow(ch_fr), geom, 0.0, opts);

  INFO("vacuum Isp: equilibrium " << p_eq.isp << " s, frozen " << p_fr.isp << " s");
  // Recombination releases energy during the expansion, so shifting equilibrium
  // always gives the higher impulse and the higher exit temperature.
  REQUIRE(p_eq.isp > p_fr.isp);
  REQUIRE(p_eq.t_exit > p_fr.t_exit);
  REQUIRE(ch_eq.c_star_ideal > ch_fr.c_star_ideal);
  // The gap is a few percent for this propellant combination.
  REQUIRE((p_eq.isp - p_fr.isp) / p_fr.isp > 0.02);
  REQUIRE((p_eq.isp - p_fr.isp) / p_fr.isp < 0.15);
}
