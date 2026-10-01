// SPDX-License-Identifier: MIT
#pragma once
/// \file NozzleKinetics.hpp
/// \brief Finite-rate recombination through the supersonic nozzle: a
///        one-dimensional kinetics (ODK) calculation.
///
/// Shifting equilibrium assumes the recombination of H, O, OH and CO keeps up
/// with the expansion all the way to the exit; frozen flow assumes it stops at
/// the throat.  Neither is true, and the real exit state lies between them.
/// This module integrates the species along the nozzle at the rates of a
/// reaction mechanism (Mechanism.hpp) and gives the specific impulse the
/// recombination actually recovers.
///
/// EQUATIONS.  Steady quasi-one-dimensional flow of an ideal-gas mixture with
/// the wall A(x) prescribed, composition n_k in mol/kg:
///
///     rho u A = mdot                          mass
///     h(T, n) + u^2/2 = h0                    energy (T recovered from it)
///     dn_k/dx = wdot_k / (rho u)              species
///     du/dx = u (-A'/A + Psi) / (1 - M_f^2)   momentum, combined with the above
///     Psi = sum(dn_k) / sum(n_k) - sum(H_k dn_k) / (c_p,f T)
///
/// with p = rho R_u T sum(n_k) and M_f the frozen Mach number (frozen sound
/// speed).  Psi is the heat release and mole change of recombination: it
/// accelerates the gas, which is where the recovered impulse comes from.
///
/// START.  The momentum equation is singular where M_f = 1, which lies just
/// downstream of the equilibrium throat (the frozen sound speed exceeds the
/// equilibrium one).  The integration starts from the shifting-equilibrium
/// solution at the station where M_f reaches `start_frozen_mach`, 1.10 by
/// default.  Upstream of it, at throat temperature and pressure, recombination
/// is fast enough that the gas is close to equilibrium, and the result is
/// insensitive to the exact start (tests/unit/test_kinetics.cpp measures how
/// insensitive).
///
/// INTEGRATOR.  The species equations are stiff: near the throat the chemistry
/// relaxes in nanoseconds and the flow crosses the nozzle in a fraction of a
/// millisecond.  The step is the linearly implicit Euler method with one
/// Richardson extrapolation (second order and L-stable), with a
/// finite-difference Jacobian.  It conserves every element exactly (linear
/// invariants pass through the implicit solve), and the energy equation is
/// satisfied exactly because T is solved from it at every evaluation.  The
/// difference between the one-step and two-half-step results controls the
/// step.
///
/// LIMITS.  With every rate multiplied by zero the gas is frozen at its start
/// composition and the result must equal a frozen isentropic expansion from
/// there; with every rate multiplied by a large factor it must approach
/// shifting equilibrium.  Both are tested.

#include <functional>
#include <string>
#include <vector>

#include <Eigen/Dense>

#include "ignis/kinetics/Mechanism.hpp"
#include "ignis/nozzle/NozzleFlow.hpp"
#include "ignis/nozzle/NozzleGeometry.hpp"

namespace ignis {

struct KineticNozzleOptions {
  /// Start where the equilibrium flow's frozen Mach number reaches this.
  double start_frozen_mach = 1.10;
  /// Scales every forward and reverse rate (1 = the mechanism's own).
  double rate_multiplier = 1.0;
  double relative_tolerance = 1.0e-6;
  double absolute_tolerance = 1.0e-12;  ///< on n_k, mol/kg
  int max_steps = 200000;
  int profile_points = 201;             ///< stations kept in the profile
};

struct KineticNozzleResult {
  double start_x = 0.0;               ///< m
  double start_area_ratio = 0.0;
  double start_frozen_mach = 0.0;
  double exit_area_ratio = 0.0;
  double mass_flow = 0.0;             ///< kg/s through the geometry

  double t_exit = 0.0, p_exit = 0.0, rho_exit = 0.0, u_exit = 0.0;
  double mach_frozen_exit = 0.0;
  Eigen::VectorXd n_exit;             ///< mol/kg, database order

  /// Vacuum specific impulse of the inviscid core, s: finite rate, and the two
  /// limits from the same geometry -- shifting equilibrium all the way, and
  /// frozen from the start station.
  double isp_vacuum = 0.0;
  double isp_vacuum_shifting = 0.0;
  double isp_vacuum_frozen = 0.0;
  double kinetic_efficiency = 1.0;    ///< isp_vacuum / isp_vacuum_shifting
  /// Share of the shifting-minus-frozen impulse the recombination recovers.
  double recovered_fraction = 1.0;

  /// The march, resampled evenly in x.
  std::vector<double> x, area_ratio, temperature, pressure, velocity, mach_frozen;
  std::vector<std::string> species;                  ///< the mechanism's species
  std::vector<std::vector<double>> mole_fraction;    ///< [species][station]

  int steps = 0, rejected_steps = 0, rhs_evaluations = 0;
  double energy_residual = 0.0;       ///< max |h + u^2/2 - h0| / |h0|
  double element_residual = 0.0;      ///< max relative change of any element
  std::vector<std::string> warnings;

  std::string summary() const;
};

/// Integrate the finite-rate nozzle from the equilibrium flow `flow` (which
/// must be a shifting-equilibrium flow) over the divergent part of `geom`.
KineticNozzleResult integrateKineticNozzle(const NozzleFlow& flow, const NozzleGeometry& geom,
                                           const Mechanism& mechanism,
                                           const KineticNozzleOptions& options = {});

/// The march itself: from the state (n0, t0, p0, u0) at x0 to x_end through
/// the passage area(x), which must start frozen-supersonic.  Fills the march,
/// exit and diagnostic fields of the result (area ratios relative to the
/// start); the nozzle's two limits are integrateKineticNozzle's.
KineticNozzleResult marchKinetic(const GasMixture& mix, const Mechanism& mechanism,
                                 const std::function<double(double)>& area, double x0,
                                 double x_end, const Eigen::VectorXd& n0, double t0, double p0,
                                 double u0, const KineticNozzleOptions& options = {});

/// Frozen isentropic expansion of composition n from (t0, p0), with
/// stagnation enthalpy h0 and mass flow mdot, to the supersonic state at
/// passage area `area`.  The frozen limit of the integration, by an
/// independent route.
struct FrozenExpansion {
  double t = 0.0, p = 0.0, rho = 0.0, u = 0.0;
  double isp_vacuum = 0.0;  ///< (u + p/(rho u)) / g0, s
};
FrozenExpansion frozenExpansion(const GasMixture& mix, const Eigen::VectorXd& n, double h0,
                                double t0, double p0, double mdot, double area);

}  // namespace ignis
