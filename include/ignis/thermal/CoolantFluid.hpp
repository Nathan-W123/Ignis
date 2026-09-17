// SPDX-License-Identifier: MIT
#pragma once
/// \file CoolantFluid.hpp
/// \brief Real-fluid properties for regenerative coolants.
///
/// A regeneratively cooled channel carries methane or hydrogen at 5-30 MPa and
/// 100-900 K -- supercritical, strongly non-ideal, and exactly where cubic
/// equations of state and the classical residual-property correlations are
/// weakest.  Ignis therefore evaluates coolant properties by interpolating a
/// table built from the published *reference* equations of state (Setzmann &
/// Wagner 1991 for methane, Leachman et al. 2009 for hydrogen, Schmidt &
/// Wagner 1985 for oxygen); see tools/build_coolant_tables.py for the full
/// citations and the regeneration command.
///
/// The table is bilinear in (T, ln p).  `PengRobinsonFluid` provides an
/// independent analytic model so that the size of the cubic-EOS approximation
/// is measured by the test suite rather than assumed.
///
/// Boiling: for p below the critical pressure the table carries the saturation
/// curve, and `at()` reports how close the state is to it.  The cooling solver
/// treats a bulk temperature at or above T_sat(p) as a hard failure rather
/// than continuing with single-phase correlations that no longer apply.

#include <string>
#include <vector>

namespace ignis {

/// Phase classification of a coolant state.
enum class CoolantPhase {
  kSupercriticalPressure,  ///< p >= p_crit: no phase boundary can be crossed
  kLiquid,                 ///< p < p_crit and T < T_sat(p)
  kVapor,                  ///< p < p_crit and T_sat(p) <= T < T_crit
  kSupercriticalTemperature,  ///< p < p_crit and T >= T_crit: a dense gas
};
std::string toString(CoolantPhase p);

/// One coolant thermodynamic/transport state.
struct CoolantState {
  double T = 0.0;              ///< K
  double p = 0.0;              ///< Pa
  double rho = 0.0;            ///< kg/m^3
  double cp = 0.0;             ///< J/(kg K)
  double h = 0.0;              ///< J/kg (table reference)
  double mu = 0.0;             ///< Pa s
  double k = 0.0;              ///< W/(m K)
  double prandtl = 0.0;
  bool supercritical_pressure = false;
  CoolantPhase phase = CoolantPhase::kSupercriticalPressure;
  /// Saturation temperature at this pressure, K.  Negative when p >= p_crit.
  double t_saturation = -1.0;
  /// True when the state is a sub-critical liquid within `kSaturationMargin`
  /// of its boiling point -- i.e. at risk of boiling in the channel.
  bool near_saturation = false;

  /// Relative temperature margin below T_sat that counts as "near".
  static constexpr double kSaturationMargin = 0.05;
};

/// Tabulated real-fluid properties.
class CoolantFluid {
 public:
  /// Load a table produced by tools/build_coolant_tables.py.
  static CoolantFluid loadCsv(const std::string& path);
  /// Locate a shipped table by fluid name ("methane", "hydrogen", "oxygen").
  static CoolantFluid load(const std::string& fluid_name);

  const std::string& name() const { return name_; }
  double criticalTemperature() const { return t_crit_; }
  double criticalPressure() const { return p_crit_; }
  double molarMass() const { return molar_mass_; }
  double tMin() const { return t_.front(); }
  double tMax() const { return t_.back(); }
  double pMin() const { return p_.front(); }
  double pMax() const { return p_.back(); }
  const std::string& provenance() const { return provenance_; }

  /// Properties at (T, p).  Throws RangeError outside the tabulated box.
  CoolantState at(double T, double p) const;
  /// Saturation temperature, K.  Throws RangeError for p >= p_crit.
  double saturationTemperature(double p) const;
  /// Invert h(T, p) for T.  Used by the channel energy march.
  double temperatureFromEnthalpy(double h, double p, double T_guess) const;

 private:
  double interp(const std::vector<double>& field, double T, double p) const;

  std::string name_, provenance_;
  std::vector<double> t_, p_, lnp_;
  std::vector<double> rho_, cp_, h_, mu_, k_;
  std::vector<double> sat_t_;   ///< saturation temperature per pressure node
  double t_crit_ = 0.0, p_crit_ = 0.0, molar_mass_ = 0.0, rho_crit_ = 0.0, acentric_ = 0.0;
};

/// Peng-Robinson cubic equation of state, used as an independent check on the
/// tabulated density and as a fallback for fluids with no table.
///
/// D.-Y. Peng and D. B. Robinson, "A New Two-Constant Equation of State",
/// Ind. Eng. Chem. Fundam. 15, 59-64 (1976).
///
///     p = R T / (v - b) - a alpha(T) / (v^2 + 2 b v - b^2)
///     a = 0.45724 R^2 Tc^2 / pc ,  b = 0.07780 R Tc / pc
///     alpha = [1 + kappa (1 - sqrt(T/Tc))]^2
///     kappa = 0.37464 + 1.54226 w - 0.26992 w^2
///
/// Accuracy: for methane between 5 and 30 MPa this reproduces the reference
/// EOS density to about 1-2 % above 250 K, degrading to roughly 12 % in the
/// dense liquid-like region near 120 K.  The test suite measures this.
class PengRobinsonFluid {
 public:
  PengRobinsonFluid(double t_crit, double p_crit, double acentric, double molar_mass);

  /// Compressibility factor of the stable (lowest-Gibbs) root.
  double compressibility(double T, double p) const;
  /// Density, kg/m^3.
  double density(double T, double p) const;
  /// Departure enthalpy h - h_ideal, J/kg.
  double enthalpyDeparture(double T, double p) const;
  /// Departure isobaric heat capacity cp - cp_ideal, J/(kg K).
  double cpDeparture(double T, double p) const;

 private:
  double t_crit_, p_crit_, acentric_, molar_mass_, a_c_, b_, kappa_;
};

}  // namespace ignis
