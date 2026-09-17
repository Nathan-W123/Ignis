// SPDX-License-Identifier: MIT
#pragma once
/// \file Equilibrium.hpp
/// \brief Constrained Gibbs-energy minimisation for ideal-gas mixtures.
///
/// PROBLEM
/// -------
/// Find the mole numbers n_j >= 0 (mol per kg of mixture) that minimise the
/// total Gibbs energy
///
///     G(n, T, p) = sum_j n_j mu_j ,
///     mu_j = g0_j(T) + R T [ ln(n_j / n) + ln(p / p0) ] ,  n = sum_j n_j
///
/// subject to the linear element-conservation constraints
///
///     sum_j a_ij n_j = b_i^0        (i = 1 .. E)
///
/// where a_ij is the number of atoms of element i in species j.
///
/// ALGORITHM
/// ---------
/// The Lagrangian stationarity conditions are
///
///     mu_j / (R T) = sum_i a_ij pi_i         (pi_i = -lambda_i / R T)
///
/// Ignis solves these with the descent method of Gordon & McBride (NASA
/// RP-1311, 1994), i.e. Newton's method applied to the stationarity and
/// constraint equations in the variables (ln n_j, ln n, pi_i [, ln T]).
/// Working in ln n_j makes every iterate strictly positive, so the
/// nonnegativity constraint is satisfied by construction rather than by
/// clipping.  Because the species corrections can be eliminated analytically,
/// each Newton step costs only an (E+1) or (E+2) square solve, independent of
/// the number of species.
///
/// Step length is controlled by the RP-1311 rules (Eqs. 3.1-3.3): a primary
/// limiter that caps |Delta ln n_j|, |Delta ln n| and |Delta ln T|, and a
/// secondary limiter that prevents trace species from overshooting.  If the
/// resulting step still increases the merit function, an Armijo backtracking
/// line search on the appropriate thermodynamic potential is applied.
///
/// SUPPORTED PROBLEMS
/// ------------------
///   kTP  fixed temperature and pressure           (minimise G)
///   kHP  fixed enthalpy and pressure (adiabatic)  (T is an unknown)
///   kSP  fixed entropy and pressure (isentropic)  (T is an unknown)
///   kUV  fixed internal energy and specific volume (outer Newton on T)
///   kTV  fixed temperature and specific volume     (outer Newton on p)
///
/// FAILURE REPORTING
/// -----------------
/// Every solve returns diagnostics containing the achieved element-balance
/// residual, the Gibbs optimality residual, the energy/entropy closure and the
/// per-iteration correction history.  A solve that does not reach tolerance
/// throws ConvergenceError; it never returns a clipped or partially converged
/// state flagged as valid.

#include <string>
#include <vector>

#include <Eigen/Dense>

#include "ignis/thermo/GasMixture.hpp"
#include "ignis/thermo/GasState.hpp"
#include "ignis/thermo/SpeciesDatabase.hpp"

namespace ignis {

/// Which thermodynamic pair is held fixed.
enum class EquilibriumProblem { kTP, kHP, kSP, kUV, kTV };

std::string toString(EquilibriumProblem p);

/// Tunable numerics for the equilibrium solver.
struct EquilibriumOptions {
  /// Convergence test on the mass-weighted species correction,
  /// max_j n_j |Delta ln n_j| / n  (RP-1311 uses 0.5e-5).
  double species_tolerance = 1.0e-9;
  /// Convergence test on |Delta ln n|.
  double moles_tolerance = 1.0e-9;
  /// Convergence test on |Delta ln T|.
  double temperature_tolerance = 1.0e-9;
  int max_iterations = 300;
  /// Mole fraction below which a species is treated as a trace: it is placed
  /// directly on its stationarity condition instead of taking a damped Newton
  /// step, which keeps trace species from destabilising the step limiter.
  double trace_fraction = 1.0e-14;
  /// Mole fraction above which a species must satisfy the Gibbs optimality
  /// condition for the solve to be accepted.
  double optimality_check_fraction = 1.0e-12;
  /// Accepted Gibbs optimality residual, max |mu_j/RT - sum_i a_ij pi_i|.
  double optimality_tolerance = 1.0e-6;
  /// Accepted relative element-balance residual.
  double element_tolerance = 1.0e-10;
  /// Initial temperature guess for problems where T is an unknown, K.
  double initial_temperature = 3800.0;
  /// Extra initial temperatures tried if the first attempt fails, K.
  std::vector<double> restart_temperatures = {2500.0, 4500.0, 1800.0, 5200.0, 1000.0};
  /// Maximum outer iterations for the kUV / kTV outer loop.
  int max_outer_iterations = 100;
  bool record_history = true;
};

/// Everything needed to judge whether a solve can be trusted.
struct EquilibriumDiagnostics {
  bool converged = false;
  int iterations = 0;
  int restarts = 0;
  int outer_iterations = 0;
  /// max_j n_j |Delta ln n_j| / n at the final step.
  double final_species_correction = 0.0;
  double final_moles_correction = 0.0;
  double final_temperature_correction = 0.0;
  /// max_i |sum_j a_ij n_j - b_i| , mol/kg.
  double element_residual_abs = 0.0;
  /// The same, divided by max_i b_i.
  double element_residual_rel = 0.0;
  /// max |mu_j/(RT) - sum_i a_ij pi_i| over species above
  /// `optimality_check_fraction`.
  double gibbs_residual = 0.0;
  /// |h - h_target| / max(1, |h_target|)  (kHP),
  /// |s - s_target| / max(1, |s_target|)  (kSP), etc.
  double state_residual = 0.0;
  /// sum_j n_j M_j - 1 : the per-kilogram normalisation error.
  double mass_residual = 0.0;
  /// max correction per iteration, for plotting solver convergence.
  std::vector<double> history;
  std::string message;
};

/// Result of an equilibrium calculation.
struct EquilibriumResult {
  GasState state;                 ///< full state incl. equilibrium derivatives
  Eigen::VectorXd pi;             ///< Lagrange multipliers, one per element
  EquilibriumDiagnostics diagnostics;
};

/// Constrained Gibbs minimiser.  Instances are immutable and therefore safe to
/// share between threads.
class EquilibriumSolver {
 public:
  explicit EquilibriumSolver(const SpeciesDatabase& db, EquilibriumOptions opts = {});

  const SpeciesDatabase& database() const { return *db_; }
  const GasMixture& mixture() const { return mix_; }
  const EquilibriumOptions& options() const { return opts_; }

  /// Equilibrium at fixed temperature and pressure.
  /// \param b element mole numbers per kg of mixture (length = numElements)
  EquilibriumResult tp(const Eigen::VectorXd& b, double T, double p,
                       const Eigen::VectorXd* n_guess = nullptr) const;

  /// Adiabatic equilibrium at fixed enthalpy and pressure.
  /// \param h_target specific enthalpy of the reactants, J/kg
  EquilibriumResult hp(const Eigen::VectorXd& b, double h_target, double p,
                       double T_guess = 0.0,
                       const Eigen::VectorXd* n_guess = nullptr) const;

  /// Isentropic equilibrium at fixed entropy and pressure.
  EquilibriumResult sp(const Eigen::VectorXd& b, double s_target, double p,
                       double T_guess = 0.0,
                       const Eigen::VectorXd* n_guess = nullptr) const;

  /// Equilibrium at fixed internal energy and specific volume.
  EquilibriumResult uv(const Eigen::VectorXd& b, double u_target, double v_target,
                       double T_guess = 0.0) const;

  /// Equilibrium at fixed temperature and specific volume.
  EquilibriumResult tv(const Eigen::VectorXd& b, double T, double v_target,
                       double p_guess = 0.0) const;

  /// Convert a per-kg element vector into the mixture mass it represents.
  /// Should equal 1 kg; used to validate caller-supplied element vectors.
  double elementMass(const Eigen::VectorXd& b) const;

 private:
  struct Reduced;  // active element / species subset

  /// Drop elements with b_i == 0 and the species that depend on them.
  static void buildReduced(const SpeciesDatabase& db, const Eigen::VectorXd& b_full,
                           Reduced& red);

  EquilibriumResult solveConstantPressure(const Eigen::VectorXd& b, EquilibriumProblem problem,
                                          double p, double target, double T_fixed_or_guess,
                                          const Eigen::VectorXd* n_guess) const;

  /// One full Newton solve from a specified start; throws on failure.
  EquilibriumResult attempt(const Reduced& red, EquilibriumProblem problem, double p,
                            double target, double T_start, const Eigen::VectorXd& n_start) const;

  /// Fill the equilibrium derivative fields (cp_eff, gamma_s, a, ...) of a
  /// converged state by solving the two auxiliary linear systems.
  void computeDerivatives(const Reduced& red, const Eigen::VectorXd& n_red, double T, double p,
                          GasState& st) const;

  const SpeciesDatabase* db_;
  GasMixture mix_;
  EquilibriumOptions opts_;
};

}  // namespace ignis
