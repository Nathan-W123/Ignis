// SPDX-License-Identifier: MIT
#pragma once
/// \file Turbomachinery.hpp
/// \brief Pumps, turbines and the gas generators and preburners that drive them.
///
/// These are the components the cycle balance (Cycle.hpp) is assembled from,
/// and they are validated on their own against published turbopump data
/// (tests/validation/test_cycle_validation.cpp).
///
/// PUMP.  Isentropic work is the integral of dp/rho along the isentrope.  Where
/// the propellant has a real-fluid table (data/coolants) that integral is taken
/// along the table: on an isentrope dh = v dp, so marching h with dp/rho(T, p)
/// and recovering T from (h, p) at each step follows the isentrope without the
/// table needing an entropy column.  This matters for liquid hydrogen, which is
/// compressible enough at 40 MPa that dp/rho_inlet overstates the work by
/// several percent.  Without a table the liquid is taken as incompressible at
/// the library density.  The actual enthalpy rise is the isentropic one over
/// the pump efficiency.
///
/// TURBINE.  Two kinds of drive gas.  Combustion products (gas generator,
/// preburner) expand isentropically at frozen composition on the species
/// thermodynamics: s(T_out,s, p_out) = s(T_in, p_in).  A heated coolant
/// (expander cycle) expands along its real-fluid table the same way the pump
/// compresses.  Power is mass flow times efficiency times the isentropic drop.
///
/// COMBUSTOR.  A gas generator or preburner runs far from stoichiometric so
/// that its gas is cool enough for a turbine.  Ignis solves the equilibrium
/// mixture ratio that gives a requested turbine inlet temperature, on the
/// fuel-rich or the oxidiser-rich branch, as the energy balance it is: the
/// reactants' enthalpy equals that of their equilibrium products at the
/// target temperature and pressure.  For hydrogen that equilibrium is a
/// good description of the gas.  For a hydrocarbon it is not: a fuel-rich
/// kerosene or methane generator makes soot and cracked hydrocarbons far from
/// equilibrium, and Ignis carries no condensed carbon (docs/limitations.md),
/// so the turbine gas there is an approximation of the right temperature and a
/// rough guess at the right composition.

#include <string>

#include <Eigen/Dense>

#include "ignis/combustion/Propellant.hpp"
#include "ignis/equilibrium/Equilibrium.hpp"
#include "ignis/thermal/CoolantFluid.hpp"
#include "ignis/thermo/GasMixture.hpp"

namespace ignis {

struct PumpResult {
  double mass_flow = 0.0;        ///< kg/s
  double p_in = 0.0, p_out = 0.0;   ///< Pa
  double t_in = 0.0, t_out = 0.0;   ///< K (t_out from the actual enthalpy rise)
  double isentropic_head = 0.0;  ///< J/kg, dh_s
  double efficiency = 0.0;
  double power = 0.0;            ///< W, absorbed
  bool real_fluid = false;       ///< integrated along a property table
};

/// Pump a liquid from (p_in, t_in) to p_out.  `table` may be null, in which
/// case the liquid is incompressible at `density`.
PumpResult pumpLiquid(const CoolantFluid* table, double density, double mass_flow, double p_in,
                      double t_in, double p_out, double efficiency);

struct TurbineResult {
  double mass_flow = 0.0;        ///< kg/s
  double p_in = 0.0, p_out = 0.0;   ///< Pa
  double t_in = 0.0, t_out = 0.0;   ///< K, after the actual expansion
  double pressure_ratio = 0.0;   ///< p_in / p_out
  double isentropic_drop = 0.0;  ///< J/kg
  double efficiency = 0.0;
  double power = 0.0;            ///< W, delivered
};

/// Combustion products of composition n (mol/kg) expanding at frozen
/// composition from (t_in, p_in) through `pressure_ratio`.
TurbineResult gasTurbine(const GasMixture& mix, const Eigen::VectorXd& n, double mass_flow,
                         double t_in, double p_in, double pressure_ratio, double efficiency);

/// A heated coolant expanding along its real-fluid table.
TurbineResult fluidTurbine(const CoolantFluid& fluid, double mass_flow, double t_in, double p_in,
                           double pressure_ratio, double efficiency);

/// A gas generator or preburner: the equilibrium state at the mixture ratio
/// that gives `turbine_inlet_temperature`.
struct CombustorResult {
  double mixture_ratio = 0.0;
  double temperature = 0.0;      ///< K, as solved
  double pressure = 0.0;         ///< Pa
  Eigen::VectorXd n;             ///< equilibrium composition, mol/kg
  double molar_mass = 0.0;       ///< kg/mol
  double cp_frozen = 0.0;        ///< J/(kg K)
  /// h_reactants - h_products at the solved mixture ratio, J/kg: how closely
  /// the adiabatic flame temperature sits on the target.
  double enthalpy_residual = 0.0;
};

/// The propellants enter at their storage temperatures (the library's
/// reference states) plus an enthalpy gain per kg of each: pump work, and for
/// a fuel that has cooled the chamber first, the jacket's heat.  Hydrogen
/// warmed that way needs far less oxygen to reach a turbine temperature.
/// \throws InfeasibleError if no mixture ratio on the requested side burns at
///         the target temperature.
CombustorResult combustorAtTemperature(const EquilibriumSolver& solver, const Propellant& oxidizer,
                                       const Propellant& fuel, double t_oxidizer, double t_fuel,
                                       double pressure, double turbine_inlet_temperature,
                                       bool fuel_rich, double oxidizer_enthalpy_gain = 0.0,
                                       double fuel_enthalpy_gain = 0.0);

}  // namespace ignis
