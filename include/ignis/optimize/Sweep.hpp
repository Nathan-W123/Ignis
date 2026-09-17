// SPDX-License-Identifier: MIT
#pragma once
/// \file Sweep.hpp
/// \brief Multi-dimensional parameter sweeps.
///
/// A sweep is a full factorial grid over one or more named parameters.  Every
/// grid point runs the same steady analysis as `ignis_engine`, so a sweep and a
/// single-point run can never disagree.  Points that fail (an infeasible
/// geometry, a non-convergent equilibrium, a boiling coolant) are recorded with
/// their error message rather than dropped, so the failed region of a design
/// space is visible in the output instead of appearing as a gap.

#include <string>
#include <vector>

#include "ignis/engine/SteadyEngine.hpp"
#include "ignis/io/Table.hpp"

namespace ignis {

/// One swept axis.
struct SweepAxis {
  std::string parameter;
  std::vector<double> values;
  /// Build a linearly spaced axis.
  static SweepAxis linear(const std::string& parameter, double lo, double hi, int points);
  /// Build a logarithmically spaced axis.
  static SweepAxis logarithmic(const std::string& parameter, double lo, double hi, int points);
};

/// Sweep definition.
struct SweepSpec {
  std::vector<SweepAxis> axes;
  std::vector<std::string> metrics;
  int threads = 0;   ///< 0 => hardware concurrency
};

/// Sweep outcome: a table with one row per grid point.
struct SweepResult {
  Table table{"sweep"};
  int points = 0;
  int failures = 0;
  std::vector<std::string> failure_messages;   ///< one per failed point
  std::string summary() const;
};

/// Run a sweep.  Deterministic and independent of the thread count.
SweepResult runSweep(const SteadyEngine& engine, const SweepSpec& spec);

/// Parse a `sweep:` section.
SweepSpec parseSweep(const ConfigNode& root);

}  // namespace ignis
