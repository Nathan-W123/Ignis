// SPDX-License-Identifier: MIT
#pragma once
/// \file MonteCarlo.hpp
/// \brief Deterministic multithreaded Monte Carlo and sensitivity analysis.
///
/// DETERMINISM
/// -----------
/// Each sample draws from its own generator, seeded by hashing the campaign
/// seed with the sample index through SplitMix64:
///
///     state = splitmix64(seed XOR splitmix64(index + GOLDEN))
///
/// Sample k therefore always sees the same random numbers no matter which
/// thread evaluates it, how many threads there are, or in what order the work
/// is scheduled.  `ignis_mc --threads N` reproduces `--threads 1` bit for bit,
/// and the test suite checks exactly that.
///
/// DISTRIBUTIONS
/// -------------
///   normal      (mean, sigma)
///   lognormal   (median, sigma_log)   -- sigma_log is the standard deviation
///                                        of ln(x); the median is exp(mu)
///   uniform     (low, high)
///   triangular  (low, mode, high)
/// Each input may also declare hard bounds; a draw outside them is re-drawn
/// (up to a limit) rather than clipped, so the realised distribution stays a
/// genuine truncated distribution.
///
/// SENSITIVITY
/// -----------
/// Two independent measures are reported:
///   * standardised regression coefficients (SRC) from a least-squares fit of
///     the standardised outputs on the standardised inputs, together with the
///     model R^2 so that the linearity assumption is visible;
///   * Spearman rank correlation, which survives monotone nonlinearity;
///   * local finite-difference elasticities at the nominal point,
///     (dy/dx)(x/y), computed by central differences with a relative step.

#include <cstdint>
#include <limits>
#include <map>
#include <string>
#include <vector>

#include "ignis/engine/SteadyEngine.hpp"
#include "ignis/io/Table.hpp"

namespace ignis {

enum class DistributionType { kNormal, kLogNormal, kUniform, kTriangular };
std::string toString(DistributionType d);
DistributionType distributionFromString(const std::string& s);

/// One dispersed input.
struct UncertainInput {
  std::string parameter;
  DistributionType distribution = DistributionType::kNormal;
  double a = 0.0;   ///< mean / median / low
  double b = 0.0;   ///< sigma / sigma_log / high / mode
  double c = 0.0;   ///< high (triangular only)
  double lower_bound = -std::numeric_limits<double>::infinity();
  double upper_bound = std::numeric_limits<double>::infinity();
  /// When true, a, b, c are multipliers applied to the nominal value.
  bool relative = false;
};

/// Campaign definition.
struct MonteCarloSpec {
  std::vector<UncertainInput> inputs;
  std::vector<std::string> outputs;
  int samples = 1000;
  std::uint64_t seed = 20260917ull;
  int threads = 0;              ///< 0 => hardware concurrency
  bool record_samples = true;
  /// Relative step for the local finite-difference sensitivities.
  double fd_step = 1.0e-3;
};

/// Summary statistics of one output.
struct Statistics {
  int count = 0;
  double mean = 0.0, stddev = 0.0, minimum = 0.0, maximum = 0.0;
  double p5 = 0.0, p50 = 0.0, p95 = 0.0, p99 = 0.0;
  double skewness = 0.0;
};

/// Sensitivity of one output to the dispersed inputs.
struct OutputSensitivity {
  std::string output;
  std::vector<double> src;        ///< standardised regression coefficients
  std::vector<double> spearman;   ///< rank correlation
  std::vector<double> elasticity; ///< local (dy/dx)(x/y) at nominal
  double r_squared = 0.0;
};

/// Campaign outcome.
struct MonteCarloResult {
  Table samples{"monte_carlo_samples"};
  std::vector<std::string> input_names, output_names;
  std::map<std::string, Statistics> statistics;
  std::vector<OutputSensitivity> sensitivity;
  int requested = 0, succeeded = 0, failed = 0;
  /// Failure message -> count.
  std::map<std::string, int> failure_modes;
  /// Counts of engineering constraint violations seen across the campaign.
  std::map<std::string, int> flags;
  std::string summary() const;
  Table sensitivityTable() const;
  Table statisticsTable() const;
};

/// Run a Monte Carlo campaign.
MonteCarloResult runMonteCarlo(const SteadyEngine& engine, const MonteCarloSpec& spec);

/// Parse a `monte_carlo:` section.
MonteCarloSpec parseMonteCarlo(const ConfigNode& root);

/// Statistics of a sample vector (sorted internally).
Statistics computeStatistics(std::vector<double> v);

}  // namespace ignis
