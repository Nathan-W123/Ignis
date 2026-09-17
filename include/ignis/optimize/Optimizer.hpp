// SPDX-License-Identifier: MIT
#pragma once
/// \file Optimizer.hpp
/// \brief Constrained design optimisation.
///
/// PROBLEM
/// -------
///     maximise (or minimise)  f(x)          a named result metric
///     subject to              g_i(x) <= 0   named metrics with bounds
///                             x_lo <= x <= x_hi
///
/// METHOD
/// ------
/// The objective is an engine analysis: expensive, not differentiable in closed
/// form, and capable of failing outright for infeasible geometry.  Ignis
/// therefore uses a derivative-free method wrapped in an augmented Lagrangian:
///
///   * outer loop -- augmented Lagrangian (Hestenes-Powell-Rockafellar)
///
///         L(x; lambda, mu) = -s f(x)
///                            + (1/(2 mu)) sum_i [ max(0, lambda_i + mu g_i(x))^2
///                                                 - lambda_i^2 ]
///
///     with multiplier update lambda_i <- max(0, lambda_i + mu g_i(x)) and a
///     penalty increase when the worst violation fails to fall by a factor of
///     four.  This drives the iterates to the constraint boundary instead of
///     stopping short of it the way a fixed quadratic penalty does.
///
///   * inner loop -- Nelder-Mead simplex with the dimension-adaptive
///     coefficients of F. Gao and L. Han, Comput. Optim. Appl. 51, 259-277
///     (2012), operating on variables scaled to the unit box and reflected at
///     the bounds.
///
///   * multi-start -- a Latin hypercube of starting points (plus the nominal
///     configuration) so that a single local basin cannot be mistaken for the
///     global optimum.
///
/// A failed analysis returns a large finite penalty and is counted, so the
/// optimiser walks around infeasible regions instead of crashing.
///
/// NO CLAIM OF GLOBAL OPTIMALITY IS MADE.  The result is the best feasible
/// point found by this procedure with this model; it is not a flight design.

#include <string>
#include <vector>

#include "ignis/engine/SteadyEngine.hpp"
#include "ignis/io/Table.hpp"

namespace ignis {

/// A bounded design variable.
struct DesignVariable {
  std::string parameter;
  double min = 0.0;
  double max = 0.0;
  double start = 0.0;   ///< 0 => take the configuration's current value
};

/// An inequality constraint on a result metric.
struct DesignConstraint {
  std::string metric;
  /// "<=" or ">=".
  std::string op = "<=";
  double bound = 0.0;
  /// Scale used to normalise the violation; 0 => |bound| or 1.
  double scale = 0.0;

  /// g(x) <= 0 form.
  double violation(double value) const;
};

/// Optimisation problem definition.
struct OptimizationSpec {
  std::vector<DesignVariable> variables;
  std::vector<DesignConstraint> constraints;
  std::string objective;
  bool maximize = true;
  int starts = 4;               ///< Latin-hypercube multi-start count
  int max_evaluations = 4000;
  double simplex_tolerance = 1.0e-7;
  int outer_iterations = 8;
  double initial_penalty = 10.0;
  double penalty_growth = 5.0;
  unsigned seed = 20260917u;
  int threads = 0;
};

/// One recorded objective evaluation.
struct OptimizationStep {
  int evaluation = 0;
  int outer = 0;
  double objective = 0.0;         ///< raw metric value (NaN when the run failed)
  double merit = 0.0;             ///< augmented-Lagrangian value
  double max_violation = 0.0;
  bool feasible = false;
  std::vector<double> x;
};

/// Optimisation outcome.
struct OptimizationResult {
  std::vector<double> x;
  std::vector<std::string> variable_names;
  double objective = 0.0;
  std::vector<double> constraint_values;
  std::vector<double> constraint_violations;
  bool feasible = false;
  int evaluations = 0;
  int failures = 0;
  std::vector<OptimizationStep> history;
  SteadyEngineResult best;
  std::string summary(const OptimizationSpec& spec) const;
  /// Convergence history as an exportable table.
  Table historyTable(const OptimizationSpec& spec) const;
};

/// Run the optimisation.
OptimizationResult optimize(const SteadyEngine& engine, const OptimizationSpec& spec);

/// Parse an `optimization:` section.
OptimizationSpec parseOptimization(const ConfigNode& root);

}  // namespace ignis
