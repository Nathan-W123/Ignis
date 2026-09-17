// SPDX-License-Identifier: MIT
#pragma once
/// \file Parameters.hpp
/// \brief Named access to configuration inputs and result outputs.
///
/// Sweeps, optimisation and Monte Carlo all need to poke a named scalar into a
/// configuration and read a named scalar out of a result.  Centralising the
/// two maps here means a new study only has to name parameters in YAML, and an
/// unknown name produces an error that lists every valid choice instead of
/// silently doing nothing.

#include <string>
#include <vector>

#include "ignis/engine/SteadyEngine.hpp"

namespace ignis {

/// Every settable input parameter, as "section.field".
std::vector<std::string> parameterNames();
/// Every readable output metric, as "section.field".
std::vector<std::string> metricNames();

/// Units of a parameter or metric, for labelling exported tables and plots.
std::string parameterUnits(const std::string& name);
std::string metricUnits(const std::string& name);

/// Write `value` into `cfg`.  Throws ConfigError for an unknown name or an
/// unusable value.
void applyParameter(EngineConfig& cfg, const std::string& name, double value);
/// Read the current value of a parameter from a configuration.
double readParameter(const EngineConfig& cfg, const std::string& name);
/// Read a metric from a completed analysis.
double readMetric(const SteadyEngineResult& result, const std::string& name);

}  // namespace ignis
