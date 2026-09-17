// SPDX-License-Identifier: MIT
#pragma once
/// \file Species.hpp
/// \brief NASA 7-coefficient polynomial species model.
///
/// The NASA-7 (a.k.a. CHEMKIN / NASA "old" or "lower-case a") form represents
/// the ideal-gas molar properties of a species on two temperature intervals
/// [Tmin, Tmid] and [Tmid, Tmax] as
///
///   cp(T)/R = a1 + a2 T + a3 T^2 + a4 T^3 + a5 T^4
///   h(T)/(R T) = a1 + a2 T/2 + a3 T^2/3 + a4 T^3/4 + a5 T^4/5 + a6/T
///   s(T)/R = a1 ln T + a2 T + a3 T^2/2 + a4 T^3/3 + a5 T^4/4 + a7
///
/// h is the *absolute* enthalpy: it already contains the enthalpy of formation
/// (a6 R is the formation enthalpy offset), so no separate heat-of-reaction
/// bookkeeping is needed anywhere in Ignis.
///
/// s(T) is the standard-state entropy at p0 = 1 bar.  The pressure dependence
/// of an ideal gas is applied by the mixture layer as -R ln(p_partial / p0).
///
/// Gibbs energy follows as g0(T) = h(T) - T s0(T), and the chemical potential
/// of species j in an ideal-gas mixture is
///   mu_j = g0_j(T) + R T ln(X_j p / p0).

#include <array>
#include <map>
#include <string>
#include <vector>

#include "ignis/core/Constants.hpp"
#include "ignis/core/Exceptions.hpp"

namespace ignis {

/// Lennard-Jones (12-6) collision parameters used by the transport model.
struct TransportData {
  enum class Geometry { kAtom, kLinear, kNonlinear };
  Geometry geometry = Geometry::kNonlinear;
  double well_depth = 0.0;          ///< eps/k_B, K
  double diameter = 0.0;            ///< sigma, Angstrom
  double dipole = 0.0;              ///< Debye
  double polarizability = 0.0;      ///< Angstrom^3
  double rotational_relaxation = 0.0;
  bool valid = false;               ///< false => no data available for this species
};

/// Molar thermodynamic state of one species at a temperature.
struct SpeciesState {
  double cp = 0.0;   ///< J/(mol K)
  double h = 0.0;    ///< J/mol  (absolute, includes enthalpy of formation)
  double s0 = 0.0;   ///< J/(mol K), standard state at p0 = 1 bar
  double g0 = 0.0;   ///< J/mol, = h - T s0
};

/// One chemical species with its NASA-7 fit.
class Species {
 public:
  Species() = default;

  Species(std::string name, std::map<std::string, int> composition, double molar_mass,
          double t_min, double t_mid, double t_max, std::array<double, 7> low,
          std::array<double, 7> high, TransportData transport, std::string source);

  const std::string& name() const { return name_; }
  /// Elemental composition, e.g. {"C":1, "H":4}.
  const std::map<std::string, int>& composition() const { return composition_; }
  /// Molar mass in kg/mol.
  double molarMass() const { return molar_mass_; }
  double tMin() const { return t_min_; }
  double tMid() const { return t_mid_; }
  double tMax() const { return t_max_; }
  const TransportData& transport() const { return transport_; }
  const std::string& source() const { return source_; }
  /// Number of atoms of `element` in one molecule (0 if absent).
  int atoms(const std::string& element) const;

  /// Molar heat capacity at constant pressure, J/(mol K).
  double cp(double T) const;
  /// Absolute molar enthalpy, J/mol.
  double h(double T) const;
  /// Standard-state molar entropy at p0 = 1 bar, J/(mol K).
  double s0(double T) const;
  /// Standard-state molar Gibbs energy, J/mol.
  double g0(double T) const { return h(T) - T * s0(T); }
  /// Chemical potential in an ideal-gas mixture, J/mol.
  /// \param X mole fraction (must be > 0), \param p pressure in Pa.
  double mu(double T, double X, double p) const;

  /// All four properties in one pass (cheaper than four separate calls).
  SpeciesState state(double T) const;

  /// Dimensionless cp/R, h/(RT), s/R -- used by the equilibrium solver.
  void reduced(double T, double& cp_R, double& h_RT, double& s_R) const;

  /// True when T lies inside the fit range [tMin, tMax].
  bool inRange(double T) const { return T >= t_min_ && T <= t_max_; }

  /// Throws RangeError when T is outside [tMin, tMax].
  void requireInRange(double T) const;

 private:
  const std::array<double, 7>& coeffs(double T) const {
    return (T < t_mid_) ? low_ : high_;
  }

  std::string name_;
  std::map<std::string, int> composition_;
  double molar_mass_ = 0.0;   // kg/mol
  double t_min_ = 0.0, t_mid_ = 0.0, t_max_ = 0.0;
  std::array<double, 7> low_{};
  std::array<double, 7> high_{};
  TransportData transport_{};
  std::string source_;
};

}  // namespace ignis
