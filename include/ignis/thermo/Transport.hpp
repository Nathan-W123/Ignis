// SPDX-License-Identifier: MIT
#pragma once
/// \file Transport.hpp
/// \brief Dilute-gas transport properties from Chapman-Enskog kinetic theory.
///
/// MODEL AND SOURCES
/// -----------------
/// Pure-species dynamic viscosity (Chapman-Enskog first approximation,
/// Hirschfelder, Curtiss & Bird, "Molecular Theory of Gases and Liquids",
/// Wiley 1954, Eq. 8.2-18):
///
///     mu_k = (5/16) sqrt(pi m_k k_B T) / (pi sigma_k^2 Omega22(T*))
///
/// with T* = k_B T / eps_k.  The reduced collision integral Omega^(2,2)* uses
/// the correlation of
///
///     P. D. Neufeld, A. R. Janzen and R. A. Aziz,
///     J. Chem. Phys. 57, 1100 (1972):
///     Omega22* = 1.16145 T*^-0.14874 + 0.52487 exp(-0.77320 T*)
///                                     + 2.16178 exp(-2.43787 T*)
///
/// valid for 0.3 <= T* <= 100 with a stated maximum error of 0.064 %.
///
/// Polar species use the Brokaw correction to the Lennard-Jones collision
/// integral (Poling, Prausnitz & O'Connell, "The Properties of Gases and
/// Liquids", 5th ed., McGraw-Hill 2001, Eq. 9-4.3):
///
///     Omega22*_polar = Omega22*_LJ(T*) + 0.2 delta*^2 / T*
///     delta* = mu_D^2 / (2 * 4 pi eps0 * eps * sigma^3)     [SI]
///
/// Pure-species thermal conductivity uses the modified Eucken correlation
/// (Poling et al., Eq. 10-3.3):
///
///     lambda_k = (mu_k / M_k) (1.32 cv_k + 1.77 R)
///
/// which reduces exactly to the monatomic Eucken result 2.5 mu cv / M.
///
/// Mixture viscosity uses Wilke's rule (C. R. Wilke, J. Chem. Phys. 18, 517,
/// 1950) and mixture conductivity the Wassiljewa form with the same
/// interaction coefficients (Mason & Saxena, Phys. Fluids 1, 361, 1958).
///
/// LIMITATIONS
/// -----------
/// * Dilute-gas theory: no density (pressure) correction.  At 20 MPa and
///   3500 K the reduced density of the combustion products is still low
///   enough for this to be a few-percent effect, but it is *not* valid for
///   dense supercritical coolant at low temperature -- see CoolantFluid,
///   which reports when it leaves the validated region.
/// * Species without Lennard-Jones data are excluded from the mixture
///   averages; `coveredMoleFraction` reports how much of the mixture was
///   actually represented.

#include <Eigen/Dense>

#include "ignis/thermo/SpeciesDatabase.hpp"

namespace ignis {

/// Result of a mixture transport evaluation.
struct TransportResult {
  double viscosity = 0.0;             ///< Pa s
  double conductivity = 0.0;          ///< W/(m K)
  double prandtl = 0.0;               ///< mu cp / lambda, dimensionless
  double covered_mole_fraction = 0.0; ///< fraction of the mixture with LJ data
};

/// Reduced collision integral Omega^(2,2)* for the Lennard-Jones 12-6
/// potential (Neufeld-Janzen-Aziz correlation).
double omega22LJ(double t_star);
/// Reduced collision integral Omega^(1,1)* (same reference).
double omega11LJ(double t_star);

/// Chapman-Enskog transport for an ideal-gas mixture.
class TransportModel {
 public:
  explicit TransportModel(const SpeciesDatabase& db);

  /// Pure-species viscosity, Pa s.
  double speciesViscosity(std::size_t j, double T) const;
  /// Pure-species thermal conductivity, W/(m K).
  double speciesConductivity(std::size_t j, double T) const;

  /// Mixture viscosity, conductivity and Prandtl number.
  /// \param X mole fractions (need not be normalised)
  /// \param cp_mass mixture specific heat, J/(kg K), used only for Pr
  TransportResult mixture(const Eigen::VectorXd& X, double T, double cp_mass) const;

  /// Lowest and highest temperature at which every species with LJ data stays
  /// inside the 0.3 <= T* <= 100 validity window of the Neufeld correlation.
  double tMinValid() const { return t_min_valid_; }
  double tMaxValid() const { return t_max_valid_; }

 private:
  const SpeciesDatabase* db_;
  std::vector<double> delta_star_;  ///< reduced dipole moment, per species
  std::vector<bool> has_data_;
  double t_min_valid_ = 0.0;
  double t_max_valid_ = 0.0;
};

}  // namespace ignis
