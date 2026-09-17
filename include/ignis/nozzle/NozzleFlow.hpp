// SPDX-License-Identifier: MIT
#pragma once
/// \file NozzleFlow.hpp
/// \brief Steady quasi-one-dimensional compressible nozzle flow.
///
/// FORMULATION
/// -----------
/// The expansion is treated as steady, adiabatic, inviscid and
/// one-dimensional, so along a streamtube
///
///     h(p) + u^2/2 = h0      (energy, h0 = chamber stagnation enthalpy)
///     s(p)         = s0      (isentropic)
///     rho u A      = mdot    (mass)
///
/// Rather than assuming a constant ratio of specific heats, Ignis marches on
/// *static pressure*: at each p the thermodynamic state is obtained either by
///
///   * frozen composition -- the chamber mole numbers are held fixed and T is
///     found from s(n_c, T, p) = s0, or
///   * shifting equilibrium -- a constant-entropy, constant-pressure Gibbs
///     minimisation is solved at every station,
///
/// and then u = sqrt(2 (h0 - h)).  The local area ratio follows from mass
/// conservation, A/At = G*/G with G = rho u.  Because G(p) rises from zero at
/// p = p0 to a maximum and falls again, the throat is located by solving
/// M(p) = 1, and each area ratio is inverted on the appropriate branch.  This
/// makes the solver exact for arbitrary variable properties: no gamma is ever
/// assumed constant.
///
/// The speed of sound uses the isentropic exponent gamma_s = -(dln p/dln v)_s,
/// which for a shifting-equilibrium mixture differs from cp/cv.
///
/// NORMAL SHOCKS
/// -------------
/// A normal shock is treated as a discontinuity of zero thickness across which
/// composition is frozen (the residence time inside a shock is far shorter
/// than any chemical time).  The jump is solved from the exact conservation
/// equations with the real caloric equation of state, not from the
/// constant-gamma relations.  Downstream of the shock the flow is subsonic and
/// carries the post-shock stagnation state.
///
/// If no shock position reproduces the commanded ambient pressure, the
/// condition is reported -- a supersonic solution is never returned as if it
/// were valid.
///
/// SEPARATION
/// ----------
/// Flow separation is an empirical, viscous phenomenon that a quasi-1D
/// inviscid model cannot predict from first principles.  Two published
/// criteria are offered and both are clearly labelled empirical; they are
/// evaluated as post-processing diagnostics and never silently alter the
/// inviscid solution.

#include <memory>
#include <string>
#include <vector>

#include "ignis/equilibrium/Equilibrium.hpp"
#include "ignis/nozzle/NozzleGeometry.hpp"

namespace ignis {

/// State at one point in the expansion.
struct ExpansionState {
  GasState gas;
  double u = 0.0;           ///< m/s
  double mach = 0.0;
  double mass_flux = 0.0;   ///< rho u, kg/(m^2 s)
  double area_ratio = 0.0;  ///< A / At
  bool supersonic = false;
};

/// Relation between exit and ambient pressure.
enum class ExpansionRegime { kUnderExpanded, kIdeallyExpanded, kOverExpanded, kSubsonicExit };
std::string toString(ExpansionRegime r);

/// Empirical separation criteria.
enum class SeparationCriterion {
  kNone,
  /// Summerfield: separation when the wall pressure falls below 0.4 p_ambient.
  /// M. Summerfield, C. Foster and W. Swan, "Flow Separation in Overexpanded
  /// Supersonic Exhaust Nozzles", Jet Propulsion 24, 319-321 (1954).
  kSummerfield,
  /// Schmucker: p_sep/p_a = (1.88 M - 1)^-0.64.
  /// R. H. Schmucker, "Flow Processes in Overexpanded Chemical Rocket
  /// Nozzles", NASA TM-77396 (1984), translation of DLR/TU-Muenchen report.
  kSchmucker,
};
std::string toString(SeparationCriterion c);
SeparationCriterion separationCriterionFromString(const std::string& s);

/// Chamber stagnation reference for an expansion.
struct ChamberReference {
  Eigen::VectorXd b;   ///< element moles per kg of mixture
  GasState stagnation; ///< chamber (stagnation) state
};

/// Quasi-1D expansion from a fixed chamber stagnation state.
///
/// Construction locates the throat, so the object is relatively expensive to
/// build and cheap to query.  Instances are immutable after construction.
class NozzleFlow {
 public:
  NozzleFlow(const EquilibriumSolver& solver, CompositionModel model, ChamberReference ref);

  CompositionModel model() const { return model_; }
  const ChamberReference& chamber() const { return ref_; }
  const GasMixture& mixture() const { return solver_->mixture(); }
  const EquilibriumSolver& solver() const { return *solver_; }

  /// State at a commanded static pressure (Pa), 0 < p < p_chamber.
  ExpansionState atPressure(double p) const;
  /// Sonic (throat) state.
  const ExpansionState& throat() const { return throat_; }
  /// State at a commanded area ratio on the requested branch.
  ExpansionState atAreaRatio(double area_ratio, bool supersonic) const;
  /// Throat mass flux rho* u*, kg/(m^2 s).
  double throatMassFlux() const { return throat_.mass_flux; }
  /// Ideal characteristic velocity c* = p_c / (rho* u*), m/s.
  double cStarIdeal() const { return ref_.stagnation.p / throat_.mass_flux; }

  /// Mass-flow error of a solution, |rho u A - mdot| / mdot.  Used by the test
  /// suite to confirm conservation.
  double massFlowResidual(const ExpansionState& st) const;
  /// Stagnation-enthalpy error, |h + u^2/2 - h0| / |h0|.
  double energyResidual(const ExpansionState& st) const;

  /// Lowest static pressure at which the species polynomials remain valid, Pa.
  double pressureFloor() const { return p_floor_; }
  /// Largest area ratio reachable before the property data runs out.
  double maxAreaRatio() const { return max_area_ratio_; }

 private:
  ExpansionState buildState(const GasState& gas) const;

  const EquilibriumSolver* solver_;
  CompositionModel model_;
  ChamberReference ref_;
  double h0_ = 0.0, s0_ = 0.0, p0_ = 0.0;
  double p_floor_ = 0.0;
  double max_area_ratio_ = 0.0;
  ExpansionState throat_;
};

/// Result of a normal-shock jump.
struct NormalShockResult {
  ExpansionState upstream;
  GasState downstream;      ///< static state behind the shock
  double u_downstream = 0.0;
  double mach_downstream = 0.0;
  double stagnation_pressure_ratio = 0.0;  ///< p02 / p01
  double entropy_rise = 0.0;               ///< J/(kg K)
};

/// Solve the normal-shock jump for a real gas at frozen composition.
/// \throws ConvergenceError if the downstream branch cannot be found,
///         InfeasibleError if the upstream state is not supersonic.
NormalShockResult solveNormalShock(const GasMixture& mix, const ExpansionState& upstream);

/// Everything one steady operating point produces.
struct NozzlePerformance {
  double p_chamber = 0.0;
  double p_ambient = 0.0;
  double p_exit = 0.0;
  double t_exit = 0.0;
  double u_exit = 0.0;
  double mach_exit = 0.0;
  double area_ratio = 0.0;
  double throat_area = 0.0;
  double exit_area = 0.0;

  double mdot = 0.0;                  ///< kg/s
  double c_star_ideal = 0.0;          ///< m/s
  double c_star = 0.0;                ///< m/s, after combustion efficiency
  double thrust_momentum = 0.0;       ///< N
  double thrust_pressure = 0.0;       ///< N
  double thrust_ideal = 0.0;          ///< N, inviscid, full-flowing, no losses
  double thrust = 0.0;                ///< N, after the declared loss factors
  double isp_ideal = 0.0;             ///< s
  double isp = 0.0;                   ///< s
  double c_effective = 0.0;           ///< m/s, = thrust / mdot
  double cf_ideal = 0.0;
  double cf = 0.0;
  double isp_vacuum = 0.0;            ///< s, same nozzle at p_ambient = 0

  /// Declared, separately visible loss factors (all default to 1).
  double eta_c_star = 1.0;            ///< combustion efficiency
  double lambda_divergence = 1.0;     ///< divergence (geometric) efficiency
  double eta_nozzle = 1.0;            ///< additional nozzle efficiency

  ExpansionRegime regime = ExpansionRegime::kIdeallyExpanded;
  bool shock_in_nozzle = false;
  double shock_area_ratio = 0.0;
  double shock_x = 0.0;
  /// Set when a shock is required but no interior position matches p_ambient.
  bool shock_unresolved = false;
  std::string shock_note;

  bool separation_predicted = false;
  double separation_area_ratio = 0.0;
  double separation_x = 0.0;
  SeparationCriterion separation_criterion = SeparationCriterion::kNone;

  double mass_flow_residual = 0.0;
  double energy_residual = 0.0;
  std::string composition_model;
};

/// Options for a steady nozzle performance evaluation.
struct NozzlePerformanceOptions {
  double eta_c_star = 1.0;
  /// Divergence efficiency.  If `auto_divergence` is true it is computed from
  /// the exit half-angle as (1 + cos theta)/2 (conical) or (1 + cos theta_e)/2
  /// (bell) -- an explicitly empirical geometric correction.
  bool auto_divergence = true;
  double lambda_divergence = 1.0;
  double eta_nozzle = 1.0;
  SeparationCriterion separation = SeparationCriterion::kSummerfield;
  bool resolve_internal_shocks = true;
  /// |p_e/p_a - 1| below which the nozzle is reported as ideally expanded.
  double ideal_tolerance = 0.01;
};

/// Evaluate one steady operating point of a nozzle.
NozzlePerformance evaluateNozzle(const NozzleFlow& flow, const NozzleGeometry& geom,
                                 double p_ambient, const NozzlePerformanceOptions& opts);

/// Flow properties sampled along the contour.
struct AxialProfile {
  std::vector<double> x, radius, area_ratio, p, T, rho, u, mach, a, gamma_s, cp, molar_mass;
  std::vector<int> supersonic;
};

/// Sample the expansion at every contour station.
AxialProfile sampleAxialProfile(const NozzleFlow& flow, const NozzleGeometry& geom);

}  // namespace ignis
