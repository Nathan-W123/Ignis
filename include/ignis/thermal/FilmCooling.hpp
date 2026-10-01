// SPDX-License-Identifier: MIT
#pragma once
/// \file FilmCooling.hpp
/// \brief Gaseous film cooling by tangential slot injection.
///
/// SOURCE
/// ------
/// J. E. Hatch and S. S. Papell, "Use of a theoretical flow model to correlate
/// data for film cooling or heating an adiabatic wall by tangential injection
/// of gases of different fluid properties", NASA TN D-130 (1959).  Equation
/// numbers below are the report's, checked against the page images.
///
///   ln eta = -[ H - 0.04 ] (S V_g / alpha_c)^0.125  f(V_g / V_c)       (12)
///   eta    = 1                                    for H < 0.04
///   f(r)   = 1 + 0.4 atan(r - 1)                  for r = V_g/V_c >= 1 (10)
///          = (1/r)^(1.5 (1/r - 1))                for r <= 1           (11)
///
/// eta = (T_ad - t_w) / (T_ad - t_c) is the effectiveness of the film on an
/// adiabatic wall: T_ad the gas stream's recovery temperature, t_c the
/// coolant's static temperature at the slot exit, t_w the wall the film
/// leaves behind.  S is the slot height, L the slot length (the chamber
/// circumference here), V_c and V_g the coolant and gas velocities at the
/// slot, alpha_c = k/(rho c_p) the coolant's diffusivity at the slot exit.
/// H = h L x / (w c_p)_c with h the gas stream's own turbulent coefficient,
/// h = 0.0265 (k/D_h) Re^0.8 Pr^0.3, every property at the mean of the gas
/// static and the coolant temperatures (the report's assumption 6).
///
/// IN A NOZZLE h, L and the gas state change along the wall, so H is the
/// integral of h L dx / (w c_p)_c from the slot -- the form the report's own
/// derivation (its eq. 2) has before it takes h and L constant.  The 0.04,
/// the 0.125 power and f are applied as fitted, with the slot quantities.
/// That generalisation is this implementation's, not the report's.
///
/// RANGE.  The report's data: gas 279-1092 K at 32-317 m/s, air and helium
/// coolants, V_g/V_c from 0.45 to 33.3, slots 1.6-12.7 mm, no reaction.  Its
/// stated accuracy is about +-5 % in wall temperature for eta from 0.2 to 1,
/// and up to 60 % in the coolant flow needed.  Reproduced against the
/// report's own Table I (helium) in tests/validation: within the stated
/// velocity range, a mean wall-temperature error of 4 %; below it, where eq.
/// (11) is extrapolated, errors of 30-50 %.  Outside 0.45-33.3 the velocity
/// ratio used in f is therefore held at the nearer limit and the caller is
/// told.  A rocket film -- dense supercritical fuel next to a 3500 K burning
/// gas -- is far outside every one of these ranges; the correlation is the
/// classical first estimate, not a validated rocket model.

#include <string>

namespace ignis {

/// The coolant at the slot exit and the gas stream it is injected into.
struct FilmSlot {
  double mass_flow = 0.0;      ///< kg/s, w
  double slot_height = 0.0;    ///< m, S
  double perimeter = 0.0;      ///< m, L (slot length)
  double temperature = 0.0;    ///< K, t_c
  double density = 0.0;        ///< kg/m^3 at the slot exit
  double cp = 0.0;             ///< J/(kg K)
  double conductivity = 0.0;   ///< W/(m K)
  double gas_velocity = 0.0;   ///< m/s, V_g at the slot

  double coolantVelocity() const;  ///< V_c = w / (rho L S)
  double diffusivity() const;      ///< alpha_c = k / (rho c_p)
};

/// Velocity ratio range of the report's data; f is evaluated inside it.
constexpr double kHatchPapellMinVelocityRatio = 0.45;
constexpr double kHatchPapellMaxVelocityRatio = 33.3;

/// f(V_g/V_c), eqs. (10) and (11), at the given ratio (no clamping).
double hatchPapellVelocityFunction(double vg_over_vc);

/// The multiplier (S V_g/alpha_c)^0.125 f(V_g/V_c) of eq. (12), with the
/// velocity ratio held inside the fitted range.  `clamped` reports whether
/// it had to be.
double hatchPapellFactor(const FilmSlot& slot, bool* clamped = nullptr);

/// Effectiveness from the heat-capacity group H and the factor above.
double hatchPapellEffectiveness(double heat_group, double factor);

/// The report's gas-side coefficient, 0.0265 (k/D) Re^0.8 Pr^0.3, from
/// properties already evaluated at the mean film temperature.
double hatchPapellGasCoefficient(double density, double velocity, double diameter,
                                 double viscosity, double conductivity, double prandtl);

}  // namespace ignis
