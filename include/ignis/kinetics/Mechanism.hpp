// SPDX-License-Identifier: MIT
#pragma once
/// \file Mechanism.hpp
/// \brief A gas-phase reaction mechanism: forward rates from a mechanism
///        file, reverse rates by detailed balance, net production rates.
///
/// FORWARD RATES.  Modified Arrhenius, k_f = A T^b exp(-Ea / RT), in three
/// forms, the ones GRI-Mech 3.0 uses for the reactions Ignis keeps
/// (tools/build_kinetics.py):
///   elementary   q = k_f prod C_i^nu_i - k_r prod C_j^nu_j
///   three-body   the same, times [M] = sum_k eps_k C_k
///   falloff      k_f = k_inf Pr/(1 + Pr) F,  Pr = k_0 [M] / k_inf,
///                F = 1 (Lindemann) or Troe's broadening,
///                log F = log F_cent / (1 + ((log Pr + c)/(n - 0.14 (log Pr + c)))^2),
///                F_cent = (1-a) e^{-T/T3} + a e^{-T/T1} + e^{-T2/T},
///                c = -0.4 - 0.67 log F_cent,  n = 0.75 - 1.27 log F_cent.
/// An explicit collider ("H + O2 + O2 <=> HO2 + O2") is a three-body reaction
/// whose only efficiency is that species'.
///
/// REVERSE RATES.  k_r = k_f / K_c with
///   K_c = exp(-dG0 / RT) (p0 / RT)^dnu,   dG0 = sum_k nu_k g0_k(T),
/// from the same species thermodynamics as the equilibrium solver.  A
/// mechanism read this way therefore relaxes to exactly the composition
/// EquilibriumSolver computes; it cannot disagree with it about where
/// equilibrium is, only about how fast the gas gets there.
///
/// UNITS.  Mol-based SI throughout: concentrations mol/m^3, rates
/// mol/(m^3 s), A in (m^3/mol)^(n-1)/s for overall order n (third body
/// counted), Ea in J/mol.

#include <string>
#include <utility>
#include <vector>

#include <Eigen/Dense>

#include "ignis/thermo/SpeciesDatabase.hpp"

namespace ignis {

struct ArrheniusRate {
  double A = 0.0;   ///< mol-based SI
  double b = 0.0;
  double Ea = 0.0;  ///< J/mol
  double operator()(double T) const;
};

enum class ReactionType { kElementary, kThreeBody, kFalloff };

struct Reaction {
  std::string equation;
  ReactionType type = ReactionType::kElementary;
  std::vector<int> reactants, products;            ///< database species indices
  std::vector<double> reactant_nu, product_nu;     ///< stoichiometric coefficients
  ArrheniusRate rate;   ///< elementary and three-body; falloff: the high-pressure limit
  ArrheniusRate low;    ///< falloff: the low-pressure limit
  bool troe = false;
  double troe_a = 0.0, troe_t3 = 0.0, troe_t1 = 0.0, troe_t2 = 0.0;
  bool troe_has_t2 = false;
  /// Third-body efficiencies: [M] = default * sum C + sum (eps_k - default) C_k.
  double default_efficiency = 1.0;
  std::vector<std::pair<int, double>> efficiencies;
  double delta_nu = 0.0;  ///< sum of product minus reactant coefficients
};

class Mechanism {
 public:
  /// Read a mechanism file (tools/build_kinetics.py's format) and map it onto
  /// `db`.  A reaction naming a species `db` does not carry cannot occur in
  /// this gas and is left out, as is one whose only collider is absent; both
  /// are listed in dropped().  Third-body efficiencies naming absent species
  /// are skipped, since they would multiply a zero concentration.
  static Mechanism loadYaml(const std::string& path, const SpeciesDatabase& db);

  /// The shipped nozzle mechanism, data/kinetics/gri30_nozzle.yaml, found by
  /// the same rules as the other data files.
  static std::string findDefault();

  const SpeciesDatabase& database() const { return *db_; }
  const std::vector<Reaction>& reactions() const { return reactions_; }
  const std::vector<std::string>& dropped() const { return dropped_; }
  const std::string& path() const { return path_; }
  /// Database indices of the species any kept reaction touches.
  const std::vector<int>& activeSpecies() const { return active_; }

  /// Net molar production rates, mol/(m^3 s), one entry per database species.
  /// \param C  molar concentrations, mol/m^3 (negative entries are treated as 0)
  /// \param rate_multiplier scales every forward and reverse rate alike, so
  ///        it changes how fast the gas relaxes but not where to
  void productionRates(double T, const Eigen::VectorXd& C, Eigen::VectorXd& wdot,
                       double rate_multiplier = 1.0) const;

  /// Forward and reverse rate constants of reaction r at T and third-body
  /// concentration M (mol/m^3), without the [M] factor of a three-body
  /// reaction.  For tests.
  void rateConstants(std::size_t r, double T, double M, double& kf, double& kr) const;

 private:
  const SpeciesDatabase* db_ = nullptr;
  std::vector<Reaction> reactions_;
  std::vector<std::string> dropped_;
  std::vector<int> active_;
  std::string path_;
};

}  // namespace ignis
