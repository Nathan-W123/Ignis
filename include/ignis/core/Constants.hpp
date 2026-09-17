// SPDX-License-Identifier: MIT
#pragma once
/// \file Constants.hpp
/// \brief Physical constants and unit conventions for Ignis.
///
/// UNIT CONVENTION (applies to the entire library unless a symbol name says
/// otherwise):
///
///   length              m
///   mass                kg
///   time                s
///   temperature         K
///   pressure            Pa
///   energy              J
///   amount of substance mol
///   molar mass          kg/mol      (note: the *database file* stores kg/kmol
///                                    == g/mol; conversion happens once at load)
///   specific enthalpy   J/kg
///   molar enthalpy      J/mol
///   specific entropy    J/(kg K)
///   molar entropy       J/(mol K)
///   mole numbers n_j    mol per kg of mixture  (the "CEA basis")
///   velocity            m/s
///   area                m^2
///   heat flux           W/m^2
///   heat-transfer coeff W/(m^2 K)
///   thermal conductivity W/(m K)
///   dynamic viscosity   Pa s
///
/// SIGN CONVENTION:
///   * Heat flux q is positive from gas -> wall -> coolant.
///   * The nozzle axial coordinate x increases from the injector face
///     (x = 0) towards the nozzle exit.
///   * Thrust is positive in the direction opposite to the exhaust velocity.

namespace ignis {

/// Fundamental and derived constants.
/// Values: CODATA 2018 (exact SI definitions where applicable).
namespace constants {

/// Ratio of a circle's circumference to its diameter.  Defined here rather
/// than relying on the POSIX M_PI macro, which is not part of standard C++.
inline constexpr double pi = 3.14159265358979323846;

/// Universal gas constant, J/(mol K).  Exact by the 2019 SI redefinition
/// (R = N_A k_B with both defined exactly).
inline constexpr double R_universal = 8.31446261815324;

/// Avogadro constant, 1/mol (exact).
inline constexpr double N_A = 6.02214076e23;

/// Boltzmann constant, J/K (exact).
inline constexpr double k_B = 1.380649e-23;

/// Standard-state pressure used by the NASA polynomials, Pa (1 bar).
/// NASA TM-4513 / Gordon & McBride use p0 = 1 bar = 1e5 Pa.
inline constexpr double p_reference = 1.0e5;

/// Standard gravitational acceleration, m/s^2 (defined, used only to convert
/// effective exhaust velocity to specific impulse in seconds).
inline constexpr double g0 = 9.80665;

/// Stefan-Boltzmann constant, W/(m^2 K^4)  (= 2 pi^5 k_B^4 / (15 h^3 c^2)).
inline constexpr double sigma_SB = 5.670374419e-8;

/// One standard atmosphere, Pa (defined).
inline constexpr double atm = 101325.0;

}  // namespace constants
}  // namespace ignis
