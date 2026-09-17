// SPDX-License-Identifier: MIT
#pragma once
/// \file FeedSystem.hpp
/// \brief Pressure-fed propellant supply, with an optional pump comparison.
///
/// MODEL
/// -----
/// For each propellant leg the tank pressure must overcome the line losses and
/// the injector pressure drop before the flow reaches the chamber:
///
///     p_tank = p_chamber + dp_injector + dp_lines + dp_dynamic
///
/// with
///     dp_injector = mdot^2 / (2 rho (Cd A)^2)      (incompressible orifice)
///     dp_lines    = [f L/D + sum K] rho V^2 / 2    (Darcy-Weisbach + fittings)
///     dp_dynamic  = rho V_line^2 / 2               (velocity head at entry)
///
/// The friction factor uses Colebrook-White, the same correlation as the
/// cooling channels.  Liquid densities are taken at the storage state.
///
/// Two questions can be asked:
///   * sizing  -- given a target chamber pressure and mass flow, what tank
///     pressure and injector area are needed?
///   * capability -- given fixed tank pressures and injector areas, what
///     chamber pressure and mixture ratio does the system actually deliver?
/// The second closes a fixed point against the chamber model: the chamber
/// pressure sets the flow, the flow sets c*, and c* sets the chamber pressure.
///
/// A minimum injector stiffness (dp_injector / p_chamber) is enforced because
/// a soft injector invites combustion instability.  Ignis reports the stiffness
/// and flags designs below the configured floor; it does not model combustion
/// stability itself.
///
/// PUMP COMPARISON
/// ---------------
/// The isentropic pump power for the same duty is reported for comparison:
///
///     P_pump = mdot dp / (rho eta_pump)
///
/// This is a conceptual sizing number.  No turbine, gas generator, bearing,
/// seal, shaft dynamics, cavitation (NPSH) or start-transient behaviour is
/// represented, and no cycle balance is closed.
///
/// OMITTED
/// -------
/// Water-hammer, feed-system coupled instability (POGO), two-phase flow,
/// cavitating venturis, valve dynamics, tank pressurisation thermodynamics
/// (the pressurant mass estimate below is an isothermal ideal-gas bound),
/// heat soak, and any structural sizing.

#include <string>

#include "ignis/combustion/Propellant.hpp"

namespace ignis {

/// One propellant leg.
struct FeedLeg {
  double mass_flow = 0.0;          ///< kg/s
  double density = 0.0;            ///< kg/m^3
  double line_length = 0.0;        ///< m
  double line_diameter = 0.0;      ///< m
  double line_roughness = 5.0e-6;  ///< m
  double fitting_k = 2.0;          ///< sum of minor-loss coefficients
  double viscosity = 1.0e-4;       ///< Pa s, for the line Reynolds number
  double injector_cd = 0.75;       ///< discharge coefficient
  double injector_area = 0.0;      ///< m^2, total effective orifice area
  double tank_pressure = 0.0;      ///< Pa (capability mode)
};

/// Result for one leg.
struct FeedLegResult {
  double mass_flow = 0.0;
  double line_velocity = 0.0;      ///< m/s
  double reynolds = 0.0;
  double friction_factor = 0.0;
  double dp_lines = 0.0;           ///< Pa
  double dp_dynamic = 0.0;         ///< Pa
  double dp_injector = 0.0;        ///< Pa
  double injector_velocity = 0.0;  ///< m/s
  double required_tank_pressure = 0.0;  ///< Pa
  double injector_area = 0.0;      ///< m^2
  double stiffness = 0.0;          ///< dp_injector / p_chamber
  double pump_power = 0.0;         ///< W, isentropic-equivalent comparison
};

/// Whole feed-system result.
struct FeedSystemResult {
  FeedLegResult oxidizer, fuel;
  double chamber_pressure = 0.0;
  double mixture_ratio = 0.0;
  double total_mass_flow = 0.0;
  double total_pump_power = 0.0;   ///< W
  bool stiffness_ok = true;
  std::string notes;
  std::string summary() const;
};

/// Sizing inputs.
struct FeedSystemSpec {
  FeedLeg oxidizer, fuel;
  double minimum_stiffness = 0.15;  ///< dp_injector / p_chamber floor
  double pump_efficiency = 0.65;    ///< for the comparison figure only
  double pump_inlet_pressure = 3.0e5;  ///< Pa, tank head for the pump case
};

/// Sizing mode: given the chamber pressure and the two mass flows, compute the
/// injector areas that hit `target_stiffness` and the required tank pressures.
FeedSystemResult sizeFeedSystem(const FeedSystemSpec& spec, double p_chamber,
                                double target_stiffness);

/// Capability mode: given fixed injector areas and tank pressures, find the
/// mass flows each leg delivers at a commanded chamber pressure.
/// \throws InfeasibleError if a tank cannot supply the chamber pressure.
FeedSystemResult evaluateFeedSystem(const FeedSystemSpec& spec, double p_chamber);

/// Isothermal ideal-gas pressurant mass needed to expel a propellant volume.
///
///     m_press = p_tank V_propellant M / (R T_press)
///
/// This is the isothermal bound; a real blowdown or regulated system with a
/// cold pressurant needs more.  Documented as a bound, not a prediction.
double pressurantMass(double tank_pressure, double propellant_volume, double pressurant_molar_mass,
                      double pressurant_temperature);

}  // namespace ignis
