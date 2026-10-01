// SPDX-License-Identifier: MIT
#pragma once
/// \file Cycle.hpp
/// \brief Turbopump power balance: gas-generator, staged-combustion and
///        expander cycles.
///
/// A pump-fed engine has to make the power its pumps absorb, and how it makes
/// it is its cycle.  Each closure below finds the operating point at which
///
///     eta_mech * P_turbine = P_pumps
///
/// with the pumps and turbines of Turbomachinery.hpp, and reports what the
/// cycle costs and what pressures it needs.
///
/// PRESSURES.  A pump's discharge pressure is whatever lies downstream of it
/// adds up to:
///     main injector     p_c (1 + s_inj)          s_inj = dp_inj / p_c
///     valves and lines  + f_line p_c
///     cooling jacket    + dp_jacket              fuel side, from the cooling solve
/// and in a staged-combustion or expander cycle the turbine as well.
///
/// GAS GENERATOR (open cycle).  A small flow burns in a separate generator at
/// the turbine inlet temperature, drives the turbine across a set pressure
/// ratio and is thrown away.  The generator is fed from the main pump
/// discharges and runs at the chamber pressure unless told otherwise.  The
/// pumps also lift the generator's own propellant, so with w_i the work per
/// kg each pump does (its isentropic head over its efficiency), x_i the
/// generator's propellant mass fractions and dh_T the turbine's work per kg,
///
///     m_gg = (m_ox w_ox + m_f w_f) / (eta_m dh_T - x_ox w_ox - x_f w_f)
///
/// in closed form.  The turbine exhaust leaves through its own duct.  Its
/// thrust is bracketed: nothing (dumped), or an ideal frozen expansion to the
/// main nozzle's exit pressure.  The delivered specific impulse, main chamber
/// and exhaust over main and generator flow, is reported for both.
///
/// STAGED COMBUSTION (closed cycle).  All of one propellant -- the fuel for a
/// fuel-rich preburner, the oxidiser for an oxidiser-rich one -- or a set
/// fraction of it burns in a preburner with a little of the other, drives
/// the turbine, and enters the main chamber; nothing is thrown away.  The
/// turbine pressure ratio PR is what the balance solves for.  The preburner
/// sits at
///     p_pb = PR p_c (1 + s_hg)                   s_hg: hot-gas injector stiffness
/// and is fed at p_pb (1 + s_pb) + f_line p_c, through the jacket first when
/// its propellant is the coolant.  The preburner's minor propellant is
/// raised to that pressure by a boost stage, as on the RS-25.  A larger PR
/// gives the turbine more work per kg and asks more of the pump feeding it, so
/// for a given turbine temperature there is a highest chamber pressure the
/// cycle can close at.  Past it there is no solution, and the result says so
/// and how short of power the best pressure ratio leaves it.
///
/// EXPANDER (closed cycle).  The fuel, heated in the jacket, drives the turbine
/// and then goes to the main injector.  The turbine works on the jacket's heat
/// alone, so the balance fixes the fuel pump's discharge pressure, and how hot
/// the jacket leaves the fuel decides whether it closes at all.  The turbine
/// expands the real fluid along its table (Turbomachinery.hpp).
///
/// ENTHALPY.  The pump work and the jacket heat stay in the propellant, so a
/// generator or preburner fed after them is fed warm, and needs less of the
/// other propellant to reach its temperature.  On the RS-25 that is the
/// difference between a fuel-preburner mixture ratio of 1.02 and the 0.87 its
/// published flows imply (tests/validation/test_cycle_validation.cpp).
///
/// WHAT IS NOT HERE.  One equivalent turbine drives all the pumps: a
/// twin-shaft engine balances each shaft separately, and this checks their
/// sum.  There are no pump or turbine maps or speeds, no NPSH or cavitation,
/// no bearing-coolant, tank-pressurisation or igniter bleeds, and no start
/// transient.  Efficiencies are inputs.  Generator and preburner gas is
/// chemical equilibrium at the turbine inlet temperature, which for a
/// fuel-rich hydrocarbon is not what a real generator makes (Turbomachinery.hpp).

#include <string>
#include <vector>

#include "ignis/combustion/Propellant.hpp"
#include "ignis/cycle/Turbomachinery.hpp"
#include "ignis/equilibrium/Equilibrium.hpp"
#include "ignis/thermal/CoolantFluid.hpp"

namespace ignis {

enum class CycleType { kGasGenerator, kStagedCombustion, kExpander };

std::string toString(CycleType type);
/// "gas_generator", "staged_combustion" or "expander".
/// \throws ConfigError otherwise.
CycleType parseCycleType(const std::string& name);

/// What the cycle is built from.  Efficiencies, stiffnesses and temperatures
/// are design inputs, not predictions.
struct CycleSpec {
  CycleType type = CycleType::kGasGenerator;

  double pump_efficiency_oxidizer = 0.70;
  double pump_efficiency_fuel = 0.70;
  double boost_pump_efficiency = 0.70;      ///< staged combustion's boost stage
  double turbine_efficiency = 0.60;
  double mechanical_efficiency = 0.98;      ///< bearings, seals, gearing
  double pump_inlet_pressure_oxidizer = 0.3e6;  ///< Pa, from the tank
  double pump_inlet_pressure_fuel = 0.3e6;      ///< Pa

  double injector_stiffness = 0.20;   ///< main injector dp / p_c
  double line_loss_fraction = 0.05;   ///< valves and lines, as a fraction of p_c

  /// Gas generator or preburner.
  double turbine_inlet_temperature = 900.0;  ///< K
  bool fuel_rich = true;

  /// Gas generator.
  double gas_generator_pressure = 0.0;   ///< Pa; 0 = the chamber pressure
  double turbine_pressure_ratio = 20.0;  ///< turbine inlet / outlet

  /// Staged combustion.
  double preburner_injector_stiffness = 0.15;  ///< dp / p_preburner
  double hot_gas_injector_stiffness = 0.10;    ///< (turbine outlet - p_c) / p_c
  double preburner_flow_fraction = 1.0;        ///< of the major propellant

  /// Expander.
  double turbine_bypass_fraction = 0.0;  ///< of the jacket flow, round the turbine
};

/// What the rest of the engine analysis supplies.
struct CycleInputs {
  double chamber_pressure = 0.0;    ///< Pa
  double oxidizer_mass_flow = 0.0;  ///< kg/s into the main chamber
  double fuel_mass_flow = 0.0;      ///< kg/s into the main chamber
  const Propellant* oxidizer = nullptr;
  const Propellant* fuel = nullptr;
  double oxidizer_temperature = 0.0;  ///< K at the pump inlet
  double fuel_temperature = 0.0;      ///< K at the pump inlet

  /// The regenerative jacket, cooled by the fuel straight from its pump.
  bool jacket = false;
  double jacket_mass_flow = 0.0;          ///< kg/s
  double jacket_pressure_drop = 0.0;      ///< Pa
  double jacket_outlet_temperature = 0.0; ///< K
  double jacket_heat = 0.0;               ///< W into the coolant
  const CoolantFluid* jacket_fluid = nullptr;  ///< the coolant table (expander turbine)
  /// A fuel pump discharge set elsewhere -- the jacket's configured inlet
  /// pressure.  Where the cycle fixes the fuel pump's discharge by what lies
  /// downstream (gas generator, oxidiser-rich staged combustion), the pump
  /// delivers at least this, and the excess is throttled before the injector.
  double fuel_discharge_floor = 0.0;      ///< Pa, 0 = none

  /// The main nozzle, for the generator exhaust and the delivered impulse.
  double exit_pressure = 0.0;     ///< Pa
  double ambient_pressure = 0.0;  ///< Pa
  double thrust = 0.0;            ///< N at the ambient pressure
  double thrust_vacuum = 0.0;     ///< N
};

struct CycleResult {
  CycleType type = CycleType::kGasGenerator;
  /// False when no operating point balances: the turbine cannot make the
  /// power at any admissible pressure ratio.  Everything else then describes
  /// the best attempt, and `infeasibility` says what fell short.
  bool feasible = true;
  std::string infeasibility;
  double chamber_pressure = 0.0;

  PumpResult oxidizer_pump, fuel_pump;
  bool has_boost_pump = false;
  PumpResult boost_pump;
  std::string boost_propellant;     ///< "oxidizer" or "fuel"
  std::string oxidizer_pump_model;  ///< real-fluid table or incompressible, and why
  std::string fuel_pump_model;
  double pump_power = 0.0;          ///< W, every pump

  bool has_combustor = false;       ///< gas generator or preburner
  CombustorResult combustor;
  TurbineResult turbine;
  double turbine_power = 0.0;       ///< W, at the shaft before mechanical losses
  /// (eta_m P_T - P_pumps) / P_pumps at the operating point.
  double power_balance_residual = 0.0;
  /// Staged and expander: the most eta_m P_T / P_pumps any admissible
  /// pressure ratio reaches.  Above one the cycle has headroom (it could close
  /// at a higher chamber pressure); below one it cannot close.
  double max_power_ratio = 0.0;

  double oxidizer_discharge_pressure = 0.0;  ///< Pa
  double fuel_discharge_pressure = 0.0;      ///< Pa, also the jacket's inlet pressure
  /// Fuel pump discharge beyond what the injector needs, dropped across a
  /// valve because `fuel_discharge_floor` asked for it, Pa.
  double fuel_throttle_loss = 0.0;
  double combustor_pressure = 0.0;           ///< Pa
  double turbine_pressure_ratio = 0.0;
  double turbine_inlet_temperature = 0.0;    ///< K

  /// Gas generator.
  double gg_mass_flow = 0.0;        ///< kg/s
  double gg_oxidizer_flow = 0.0, gg_fuel_flow = 0.0;
  double gg_flow_fraction = 0.0;    ///< of the engine's total flow
  double exhaust_exit_pressure = 0.0;  ///< Pa
  double exhaust_velocity = 0.0;       ///< m/s, ideal frozen expansion
  double exhaust_thrust_vacuum = 0.0;  ///< N
  double exhaust_thrust = 0.0;         ///< N at the ambient pressure

  /// Staged combustion.
  double preburner_major_flow = 0.0;  ///< kg/s
  double preburner_minor_flow = 0.0;  ///< kg/s

  /// The whole engine: for a closed cycle these are the main chamber's.
  double total_mass_flow = 0.0;           ///< kg/s
  double overall_mixture_ratio = 0.0;
  double thrust_vacuum_delivered = 0.0;   ///< N, with the exhaust's ideal thrust
  double thrust_delivered = 0.0;          ///< N at the ambient pressure
  double isp_vacuum_delivered = 0.0;      ///< s
  double isp_delivered = 0.0;             ///< s at the ambient pressure
  double isp_vacuum_dumped = 0.0;         ///< s, with the exhaust's thrust counted as nothing

  std::vector<std::string> warnings;
  std::string summary() const;
};

/// The temperature a pump delivers `prop` at, pumped from (t_in, p_in) to
/// p_out the way the cycle pumps it -- for starting a jacket solve before the
/// cycle has closed.  An incompressible liquid leaves at its inlet temperature.
double pumpDischargeTemperature(const Propellant& prop, double t_in, double p_in, double p_out,
                                double efficiency);

/// Close the cycle.
/// \throws ConfigError for inputs that cannot describe a cycle (an expander
///         without a fuel-cooled jacket, a generator above the chamber
///         pressure, efficiencies outside (0, 1]).  An operating point that
///         simply does not balance is reported, not thrown.
CycleResult solveCycle(const CycleSpec& spec, const CycleInputs& inputs,
                       const EquilibriumSolver& solver);

}  // namespace ignis
