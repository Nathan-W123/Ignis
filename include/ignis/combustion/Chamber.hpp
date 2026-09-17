// SPDX-License-Identifier: MIT
#pragma once
/// \file Chamber.hpp
/// \brief Combustion-chamber thermochemistry and characteristic velocity.
///
/// The chamber is treated as an *infinite-area combustor*: the propellants
/// burn to equilibrium at the commanded chamber pressure with no inlet
/// momentum, so the chamber state is a stagnation state.  This is the same
/// assumption NASA CEA makes by default and it is exact in the limit of a
/// large contraction ratio.  The contraction ratio actually used by the
/// geometry only enters through residence time, L* and the heat-transfer
/// model, never through the chamber thermodynamics.
///
/// IDEAL vs. CORRECTED
/// -------------------
/// Everything computed here is *ideal*: complete, adiabatic, equilibrium
/// combustion.  The combustion efficiency eta_c* is applied separately and
/// explicitly, so the two never blur:
///
///     c*_actual = eta_c* c*_ideal
///     mdot      = p_c A_t / c*_actual
///
/// eta_c* is an empirical input (typically 0.92-0.99 for a well-designed
/// injector); Ignis never invents a value for it.

#include <string>

#include "ignis/combustion/Propellant.hpp"
#include "ignis/equilibrium/Equilibrium.hpp"
#include "ignis/nozzle/NozzleFlow.hpp"

namespace ignis {

/// Result of one chamber thermochemistry evaluation.
struct ChamberResult {
  GasState state;                  ///< equilibrium chamber (stagnation) state
  Eigen::VectorXd element_moles;   ///< mol/kg, the conserved element vector
  double reactant_enthalpy = 0.0;  ///< J/kg
  double mixture_ratio = 0.0;
  double equivalence_ratio = 0.0;
  double stoichiometric_mixture_ratio = 0.0;

  double c_star_ideal = 0.0;       ///< m/s
  double eta_c_star = 1.0;
  double c_star = 0.0;             ///< m/s, = eta_c* * c_star_ideal
  double throat_mass_flux = 0.0;   ///< rho* u*, kg/(m^2 s), ideal
  GasState throat_state;           ///< sonic state
  double throat_mach = 0.0;

  CompositionModel model = CompositionModel::kEquilibrium;
  EquilibriumDiagnostics diagnostics;

  /// Throat area needed to pass `mdot` at this chamber pressure, m^2.
  double throatAreaFor(double mdot) const;
  /// Mass flow through a throat of area `at`, kg/s.
  double massFlowFor(double at) const;
  /// Residence time V / (mdot / rho_c), s.
  double residenceTime(double chamber_volume, double mdot) const;
  /// Characteristic length L* = V / At, m.
  static double characteristicLength(double chamber_volume, double throat_area) {
    return chamber_volume / throat_area;
  }
  std::string summary() const;
};

/// Chamber thermochemistry driver.
class CombustionChamber {
 public:
  CombustionChamber(const EquilibriumSolver& solver, CompositionModel model);

  /// Solve the chamber state for a propellant mixture at a chamber pressure.
  /// \param eta_c_star combustion efficiency, applied only to c*.
  ChamberResult solve(const PropellantMixture& mix, double p_chamber,
                      double eta_c_star = 1.0) const;

  /// Build the expansion object that matches a solved chamber state.
  NozzleFlow makeFlow(const ChamberResult& res) const;

  const EquilibriumSolver& solver() const { return *solver_; }
  CompositionModel model() const { return model_; }

 private:
  const EquilibriumSolver* solver_;
  CompositionModel model_;
};

}  // namespace ignis
