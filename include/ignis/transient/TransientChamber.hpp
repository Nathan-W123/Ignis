// SPDX-License-Identifier: MIT
#pragma once
/// \file TransientChamber.hpp
/// \brief Zero-dimensional transient chamber mass and energy model.
///
/// STATE AND EQUATIONS
/// -------------------
/// The chamber is a single well-stirred control volume of fixed volume V.  The
/// state is the mass of oxidiser-derived and fuel-derived material it holds
/// plus its total internal energy:
///
///   d m_ox / dt = mdot_ox_in - mdot_out (m_ox / m)
///   d m_f  / dt = mdot_f_in  - mdot_out (m_f  / m)
///   d U    / dt = mdot_ox_in h_ox + mdot_f_in h_f - mdot_out h_c
///                 - (1 - eta_heat(t)) (mdot_ox_in + mdot_f_in) q_comb(O/F)
///                 + Qdot_wall
///
/// with m = m_ox + m_f.  Tracking the two masses separately gives the
/// instantaneous chamber mixture ratio exactly, including the transient
/// excursions caused by asymmetric valve timing.
///
/// Enthalpies are absolute (they include the enthalpies of formation), so
/// combustion needs no source term: injecting cold propellant and letting the
/// mixture equilibrate *is* the heat release.  The ignition ramp is therefore
/// modelled by explicitly *withholding* part of the chemical energy:
/// eta_heat(t) is the fraction of the heat of combustion released, and the
/// remainder leaves the chamber with the unburned propellant.
///
/// eta_heat is not the same quantity as the steady-state combustion efficiency
/// eta_c*.  Withholding a fraction (1 - eta_heat) of the heat of combustion
/// q_comb lowers the chamber temperature by roughly
///     dT = (1 - eta_heat) q_comb / cp
/// and, since c* scales as sqrt(T), the two are related by
///     eta_c* ~ sqrt(1 - (1 - eta_heat) q_comb / (cp T_adiabatic)) .
/// For the shipped LOX/methane case q_comb / (cp T_ad) is about 0.36, so
/// eta_c* = 0.96 corresponds to eta_heat near 0.78, not 0.92.  The shipped
/// transient configuration uses the calibrated value so that the steady part of
/// the run reproduces the steady-state analysis.
///
/// CLOSURE
/// -------
/// Given (m_ox, m_f, U) the chamber temperature and pressure follow from
///   u = U/m,  rho = m/V,  O/F = m_ox/m_f
/// by inverting the tabulated equilibrium equation of state (EquilibriumTable).
///
/// OUTFLOW
/// -------
/// Choked:   mdot_out = p A_t / c*(O/F, p)     -- exact for the ideal
///                                                equilibrium nozzle
/// Unchoked: the standard compressible orifice relation with the local
///           isentropic exponent.  If p_c <= p_ambient the outflow is set to
///           zero and the sample is flagged: Ignis does not model backflow.
///
/// INTEGRATORS
/// -----------
///   * classical fixed-step RK4
///   * adaptive embedded Cash-Karp RK4(5) with PI step-size control
/// Both are exercised by the verification suite, which measures the observed
/// order of convergence.
///
/// SAFEGUARDS
/// ----------
/// A chamber that starts nearly empty makes u = U/m stiff.  The model requires
/// a finite, physically meaningful initial charge (the chamber is pre-filled
/// with combustion products at the ambient pressure by default), enforces a
/// mass floor, and reports a clear error rather than integrating through a
/// singularity.

#include <string>
#include <vector>

#include "ignis/transient/EquilibriumTable.hpp"

namespace ignis {

/// Piecewise-linear schedule of a quantity against time.
class Schedule {
 public:
  Schedule() = default;
  Schedule(std::vector<double> times, std::vector<double> values);

  /// Trapezoidal open/close ramp.
  static Schedule ramp(double t_open, double ramp_up, double steady_value, double t_close,
                       double ramp_down);
  /// Constant value for all time.
  static Schedule constant(double value);

  double at(double t) const;
  const std::vector<double>& times() const { return t_; }
  const std::vector<double>& values() const { return v_; }

 private:
  std::vector<double> t_, v_;
};

/// Transient run definition.
struct TransientSpec {
  double chamber_volume = 0.0;       ///< m^3
  double throat_area = 0.0;          ///< m^2
  double exit_area = 0.0;            ///< m^2, used only for the thrust estimate
  double ambient_pressure = 101325.0;///< Pa
  Schedule oxidizer_flow;            ///< kg/s
  Schedule fuel_flow;                ///< kg/s
  Schedule combustion_efficiency;    ///< eta_heat(t), 0..1
  double oxidizer_inlet_enthalpy = 0.0;  ///< J/kg (absolute)
  double fuel_inlet_enthalpy = 0.0;      ///< J/kg (absolute)
  double initial_pressure = 101325.0;///< Pa
  double initial_temperature = 300.0;///< K
  double initial_mixture_ratio = 0.0;///< O/F of the initial charge; 0 => table midpoint
  double wall_heat_rate = 0.0;       ///< W removed from the gas (positive = loss)
  double t_end = 1.0;                ///< s
  double dt = 1.0e-5;                ///< s, fixed step or initial adaptive step
  double dt_output = 1.0e-4;         ///< s
  std::string integrator = "rk45";   ///< "rk4" or "rk45"
  double rtol = 1.0e-8;
  double atol = 1.0e-10;
  double dt_min = 1.0e-10;
  double dt_max = 1.0e-3;
  double mass_floor = 1.0e-9;        ///< kg
};

/// One recorded instant.
struct TransientSample {
  double t = 0.0;
  double mass = 0.0;              ///< kg
  double mass_oxidizer = 0.0;     ///< kg
  double mass_fuel = 0.0;         ///< kg
  double internal_energy = 0.0;   ///< J
  double pressure = 0.0;          ///< Pa
  double temperature = 0.0;       ///< K
  double density = 0.0;           ///< kg/m^3
  double mixture_ratio = 0.0;
  double molar_mass = 0.0;        ///< kg/mol
  double gamma_s = 0.0;
  double c_star = 0.0;            ///< m/s
  double mdot_ox_in = 0.0, mdot_fuel_in = 0.0, mdot_out = 0.0;  ///< kg/s
  double eta_heat = 0.0;
  double thrust = 0.0;            ///< N, quasi-steady estimate
  bool choked = false;
  /// True when the *inlet stream* ratio was outside the tabulated range, so the
  /// releasable chemical energy was computed from the nearest tabulated ratio
  /// with the remainder treated as inert.  The chamber state itself is never
  /// clamped: a chamber ratio outside the table stops the run.
  bool inlet_mixture_ratio_clamped = false;
};

/// Result of a transient run.
struct TransientResult {
  std::vector<TransientSample> samples;
  int steps = 0;
  int rejected_steps = 0;
  int rhs_evaluations = 0;
  double max_pressure = 0.0;
  double time_to_90_percent = -1.0;  ///< s, from first flow to 90 % of the peak
  /// |(m_end - m_start) - integral of net mass flux| / max mass, where the
  /// integral is accumulated with the same Runge-Kutta weights as the state.
  /// The conservative form makes this exact in exact arithmetic, so the value
  /// reports accumulated round-off and confirms the state update is consistent.
  /// Global integration accuracy is measured separately by grid refinement in
  /// the verification suite.
  double mass_conservation_error = 0.0;
  /// The same test on the energy equation.
  double energy_conservation_error = 0.0;
  /// Number of recorded samples whose *inlet stream* ratio lay outside the
  /// table (normal at the valve edges, where one propellant leads the other).
  int inlet_clamped_samples = 0;
  double min_mixture_ratio = 0.0;
  double max_mixture_ratio = 0.0;
  double min_temperature = 0.0;
  double max_temperature = 0.0;
  std::string integrator;
  std::string summary() const;
};

/// Integrate a startup/shutdown transient.
TransientResult simulateTransient(const EquilibriumTable& table, const TransientSpec& spec);

}  // namespace ignis
