// SPDX-License-Identifier: MIT
#pragma once
/// \file GasMixture.hpp
/// \brief Ideal-gas mixture property evaluation on a per-kilogram basis.

#include <Eigen/Dense>

#include "ignis/thermo/GasState.hpp"
#include "ignis/thermo/SpeciesDatabase.hpp"

namespace ignis {

/// Ideal-gas mixture thermodynamics over a fixed species set.
///
/// All routines take a mole-number vector `n` in mol per kg of mixture.  The
/// caller is responsible for the normalisation sum_j n_j M_j = 1; `normalize()`
/// enforces it and `massResidual()` measures it.
class GasMixture {
 public:
  explicit GasMixture(const SpeciesDatabase& db) : db_(&db) {}

  const SpeciesDatabase& database() const { return *db_; }
  std::size_t numSpecies() const { return db_->size(); }

  /// sum_j n_j M_j - 1, which must vanish for a consistent per-kg basis.
  double massResidual(const Eigen::VectorXd& n) const;
  /// Scale `n` so that sum_j n_j M_j == 1 exactly.
  Eigen::VectorXd normalize(const Eigen::VectorXd& n) const;
  /// Build a per-kg mole-number vector from mole fractions.
  Eigen::VectorXd fromMoleFractions(const Eigen::VectorXd& X) const;
  /// Build a per-kg mole-number vector from mass fractions.
  Eigen::VectorXd fromMassFractions(const Eigen::VectorXd& Y) const;

  /// Mean molar mass, kg/mol.
  double molarMass(const Eigen::VectorXd& n) const;
  /// Specific enthalpy, J/kg (absolute: includes formation enthalpies).
  double enthalpy(const Eigen::VectorXd& n, double T) const;
  /// Specific internal energy, J/kg.
  double internalEnergy(const Eigen::VectorXd& n, double T) const;
  /// Frozen specific heat at constant pressure, J/(kg K).
  double cpFrozen(const Eigen::VectorXd& n, double T) const;
  /// Specific entropy, J/(kg K).  Includes the -R ln(X p/p0) mixing terms.
  double entropy(const Eigen::VectorXd& n, double T, double p) const;
  /// Element mole numbers per kg, b_i = sum_j a_ij n_j.
  Eigen::VectorXd elementMoles(const Eigen::VectorXd& n) const;

  /// Assemble a complete frozen-composition state.
  GasState frozenState(const Eigen::VectorXd& n, double T, double p) const;

  /// Solve enthalpy(n, T) = h for T by safeguarded Newton iteration.
  /// \throws ConvergenceError if the bracket cannot be reduced.
  double temperatureFromEnthalpy(const Eigen::VectorXd& n, double h, double T_guess) const;
  /// Solve internalEnergy(n, T) = u for T.
  double temperatureFromInternalEnergy(const Eigen::VectorXd& n, double u, double T_guess) const;
  /// Solve entropy(n, T, p) = s for T at fixed p and composition.
  double temperatureFromEntropy(const Eigen::VectorXd& n, double s, double p, double T_guess) const;

  /// Check every species temperature range; throws RangeError on violation.
  void requireInRange(double T) const;

 private:
  const SpeciesDatabase* db_;
};

}  // namespace ignis
