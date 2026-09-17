// SPDX-License-Identifier: MIT
#pragma once
/// \file GasState.hpp
/// \brief Complete thermodynamic state of an ideal-gas mixture.

#include <string>
#include <vector>

#include <Eigen/Dense>

namespace ignis {

/// How the composition of a state is allowed to respond to changes in T and p.
enum class CompositionModel {
  kFrozen,       ///< composition held fixed; derivative terms are the ideal-gas ones
  kEquilibrium,  ///< composition re-minimises Gibbs energy at every state
};

std::string toString(CompositionModel m);
CompositionModel compositionModelFromString(const std::string& s);

/// Thermodynamic state of an ideal-gas mixture on a *per kilogram* basis.
///
/// The mole-number vector `n` holds mol of each species per kg of mixture (the
/// basis used by NASA CEA), so that sum_j n_j * M_j = 1 kg by construction.
struct GasState {
  double T = 0.0;                  ///< K
  double p = 0.0;                  ///< Pa
  Eigen::VectorXd n;               ///< mol/kg, one entry per species
  double n_total = 0.0;            ///< mol/kg
  double M = 0.0;                  ///< mean molar mass, kg/mol
  double R = 0.0;                  ///< specific gas constant, J/(kg K)
  double rho = 0.0;                ///< kg/m^3
  double v = 0.0;                  ///< specific volume, m^3/kg
  double h = 0.0;                  ///< specific enthalpy (absolute), J/kg
  double u = 0.0;                  ///< specific internal energy (absolute), J/kg
  double s = 0.0;                  ///< specific entropy, J/(kg K)
  double g = 0.0;                  ///< specific Gibbs energy, J/kg

  double cp_frozen = 0.0;          ///< J/(kg K), composition held fixed
  double cv_frozen = 0.0;          ///< J/(kg K)
  double gamma_frozen = 0.0;       ///< cp_frozen / cv_frozen

  /// Effective properties including composition shift.  For a frozen state
  /// these equal the frozen values and the log-derivatives are exactly 1/-1.
  double cp_eff = 0.0;             ///< J/(kg K)
  double cv_eff = 0.0;             ///< J/(kg K)
  double gamma_eff = 0.0;          ///< cp_eff / cv_eff (NOT the isentropic exponent)
  double dlnV_dlnT_p = 1.0;        ///< (d ln v / d ln T)_p
  double dlnV_dlnp_T = -1.0;       ///< (d ln v / d ln p)_T
  /// Isentropic exponent gamma_s = -(d ln p / d ln v)_s.  Equal to
  /// gamma_frozen for a frozen mixture; smaller when the mixture recombines.
  double gamma_s = 0.0;
  double a = 0.0;                  ///< speed of sound sqrt(gamma_s * p / rho), m/s

  CompositionModel model = CompositionModel::kFrozen;

  /// Mole fractions.
  Eigen::VectorXd moleFractions() const { return n / n_total; }
  /// Mass fractions.
  Eigen::VectorXd massFractions(const Eigen::VectorXd& molar_masses) const {
    return n.cwiseProduct(molar_masses);
  }
};

}  // namespace ignis
