// SPDX-License-Identifier: MIT
#pragma once
/// \file HeatTransfer.hpp
/// \brief Hot-gas-side heat transfer and through-wall conduction.
///
/// SCOPE AND HONESTY
/// -----------------
/// This is an *engineering estimate*, not a conjugate CFD solution.  The
/// hot-gas film coefficient comes from a 1957 correlation fitted to a handful
/// of nozzles; real engines differ from it by tens of percent, which is why
/// `bartz_multiplier` exists and is one of the dispersed Monte Carlo inputs.
/// There is no boundary-layer solution, no film or transpiration cooling, no
/// combustion-zone non-uniformity, no injector streak, and no soot layer.
///
/// BARTZ CORRELATION
/// -----------------
/// D. R. Bartz, "A Simple Equation for Rapid Estimation of Rocket Nozzle
/// Convective Heat Transfer Coefficients", Jet Propulsion 27, 49-51 (1957):
///
///   h_g = (0.026 / Dt^0.2) (mu^0.2 cp / Pr^0.6) (p_c / c*)^0.8
///         (Dt / R_curv)^0.1 (At / A)^0.9 sigma
///
///   sigma = 1 / { [0.5 (T_wg/T_c)(1 + (g-1)/2 M^2) + 0.5]^0.68
///                 [1 + (g-1)/2 M^2]^0.12 }
///
/// mu, cp and Pr are evaluated at the chamber *stagnation* condition, as in
/// Bartz's original derivation.  Written in SI the correlation is
/// dimensionally consistent and needs no unit conversion factor.
///
/// RECOVERY TEMPERATURE
/// --------------------
///   T_aw = T_static [1 + r (g-1)/2 M^2],  r = Pr^(1/3)  (turbulent)
///
/// WALL CONDUCTION
/// ---------------
/// One-dimensional radial conduction through a cylindrical shell with a
/// temperature-dependent conductivity evaluated at the mean wall temperature:
///
///   q_inner = k(T_mean) (T_wg - T_wc) / (r_i ln(r_o / r_i))
///
/// The planar limit q = k (T_wg - T_wc)/t is recovered as t/r_i -> 0 and is
/// available for comparison.
///
/// RADIATION
/// ---------
/// Optional and kept strictly separate from convection:
///
///   q_rad = sigma_SB eps_wall eps_gas (T_gas^4 - T_wg^4)
///
/// eps_gas is a user input with no default.  A gray-gas emissivity is a crude
/// stand-in for the H2O/CO2 band structure; a Leckner- or Hottel-style band
/// model would be required for quantitative radiative flux.

#include <string>
#include <vector>

#include "ignis/core/Constants.hpp"
#include "ignis/core/Exceptions.hpp"

namespace ignis {

/// Wall material with a linear conductivity model.
struct WallMaterial {
  std::string name = "CuCrZr";
  double conductivity_reference = 320.0;  ///< W/(m K) at reference_temperature
  double reference_temperature = 300.0;   ///< K
  double conductivity_slope = 0.0;        ///< W/(m K^2)
  double valid_min = 250.0, valid_max = 900.0;  ///< K
  double max_temperature = 800.0;         ///< K, continuous-service limit
  double emissivity = 0.3;
  double density = 8900.0;                ///< kg/m^3
  std::string source;

  /// Thermal conductivity at a temperature, W/(m K).
  double conductivity(double T) const;
  /// True when T lies inside the fitted range of the conductivity model.
  bool inValidRange(double T) const { return T >= valid_min && T <= valid_max; }
};

/// The shipped material database.
class MaterialLibrary {
 public:
  static MaterialLibrary loadYaml(const std::string& path);
  static MaterialLibrary loadDefault();
  const WallMaterial& at(const std::string& name) const;
  std::vector<std::string> names() const;

 private:
  std::vector<WallMaterial> materials_;
};

/// Inputs the Bartz correlation needs from the chamber.
struct BartzReference {
  double throat_diameter = 0.0;     ///< m
  double curvature_radius = 0.0;    ///< m, throat wall radius of curvature
  double chamber_pressure = 0.0;    ///< Pa
  double c_star = 0.0;              ///< m/s
  double chamber_temperature = 0.0; ///< K (stagnation)
  double viscosity = 0.0;           ///< Pa s, at chamber stagnation
  double cp = 0.0;                  ///< J/(kg K), at chamber stagnation
  double prandtl = 0.0;             ///< at chamber stagnation
  double multiplier = 1.0;          ///< explicit empirical scaling
};

/// Hot-gas film coefficient from the Bartz correlation, W/(m^2 K).
/// \param area_ratio  local A/At (always >= 1, both sides of the throat)
/// \param mach        local Mach number
/// \param gamma       local isentropic exponent
/// \param t_wall_hot  hot-side wall temperature, K
double bartzFilmCoefficient(const BartzReference& ref, double area_ratio, double mach,
                            double gamma, double t_wall_hot);

/// Adiabatic (recovery) wall temperature, K.
double recoveryTemperature(double t_static, double mach, double gamma, double prandtl);

/// Radial conduction resistance of a cylindrical shell referred to its inner
/// surface, (m^2 K)/W.
double cylindricalWallResistance(double inner_radius, double thickness, double conductivity);

/// Planar conduction resistance referred to either surface, (m^2 K)/W.
inline double planarWallResistance(double thickness, double conductivity) {
  if (!(conductivity > 0.0)) throw ConfigError("wall conduction: conductivity must be positive");
  return thickness / conductivity;
}

/// Gray-gas radiative flux from the combustion gas to the wall, W/m^2.
double grayGasRadiation(double t_gas, double t_wall, double gas_emissivity,
                        double wall_emissivity);

}  // namespace ignis
