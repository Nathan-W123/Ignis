// SPDX-License-Identifier: MIT
#pragma once
/// \file Exceptions.hpp
/// \brief Ignis exception hierarchy.
///
/// Design rule enforced throughout the library: a computation that cannot
/// produce a physically valid answer throws.  It never clamps, never returns a
/// silently-fabricated state, and never reports a converged flag it did not
/// earn.  Callers that want tolerant behaviour (sweeps, Monte Carlo) catch
/// these and record the failure mode explicitly.

#include <stdexcept>
#include <string>

namespace ignis {

/// Base class for every error raised by Ignis.
class IgnisError : public std::runtime_error {
 public:
  explicit IgnisError(const std::string& what) : std::runtime_error(what) {}
};

/// Invalid, missing or self-inconsistent user configuration.
class ConfigError : public IgnisError {
 public:
  explicit ConfigError(const std::string& what) : IgnisError("config error: " + what) {}
};

/// A property was requested outside the validated range of its correlation or
/// polynomial fit.
class RangeError : public IgnisError {
 public:
  explicit RangeError(const std::string& what) : IgnisError("range error: " + what) {}
};

/// An iterative solver failed to reach its tolerance.
class ConvergenceError : public IgnisError {
 public:
  explicit ConvergenceError(const std::string& what)
      : IgnisError("convergence failure: " + what) {}
};

/// The requested state is physically impossible (e.g. negative area, a shock
/// solution that does not exist, a feed system that cannot supply the demanded
/// chamber pressure).
class InfeasibleError : public IgnisError {
 public:
  explicit InfeasibleError(const std::string& what)
      : IgnisError("infeasible: " + what) {}
};

}  // namespace ignis
