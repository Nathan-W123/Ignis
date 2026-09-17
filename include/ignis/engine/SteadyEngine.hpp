// SPDX-License-Identifier: MIT
#pragma once
/// \file SteadyEngine.hpp
/// \brief One complete steady-state engine analysis.
///
/// This module owns the assembly order that the individual physics modules
/// deliberately do not know about:
///
///   1. species set  ->  equilibrium solver  ->  transport model
///   2. propellants  ->  chamber equilibrium  ->  c* and throat state
///   3. contour      ->  quasi-1D expansion   ->  thrust, Isp, Cf, regime
///   4. expansion + contour + chamber -> Bartz -> wall -> coolant channels
///   5. mass flows   ->  feed system
///
/// Nothing here introduces new physics; it only wires the modules together and
/// carries the results, so an application, a sweep, an optimiser and a Monte
/// Carlo sample all run exactly the same analysis.

#include <memory>
#include <string>

#include "ignis/combustion/Chamber.hpp"
#include "ignis/cycle/FeedSystem.hpp"
#include "ignis/io/Config.hpp"
#include "ignis/io/Json.hpp"
#include "ignis/io/Table.hpp"
#include "ignis/nozzle/Atmosphere.hpp"
#include "ignis/thermal/RegenerativeCooling.hpp"

namespace ignis {

/// Everything one steady operating point produces.
struct SteadyEngineResult {
  std::string name;
  ChamberResult chamber;
  NozzlePerformance performance;
  AxialProfile profile;

  double mdot = 0.0;          ///< kg/s total
  double mdot_oxidizer = 0.0; ///< kg/s
  double mdot_fuel = 0.0;     ///< kg/s
  double residence_time = 0.0;///< s
  double l_star = 0.0;        ///< m
  double chamber_volume = 0.0;///< m^3
  double throat_area = 0.0;   ///< m^2
  double exit_area = 0.0;     ///< m^2
  double ambient_pressure = 0.0;
  double altitude = 0.0;
  bool altitude_known = false;

  bool has_cooling = false;
  CoolingResult cooling;
  bool has_feed = false;
  FeedSystemResult feed;

  std::string summary() const;
  Json toJson() const;
  /// Axial profile as an exportable table.
  Table profileTable() const;
  /// Thermal stations as an exportable table (empty when cooling is off).
  Table coolingTable() const;
};

/// Steady engine driver.  Construction is the expensive part (loading the
/// species database and building the solver); `run()` is cheap enough to call
/// thousands of times from a sweep or a Monte Carlo campaign.
class SteadyEngine {
 public:
  explicit SteadyEngine(EngineConfig config);

  /// Run at the ambient pressure implied by the configuration.
  SteadyEngineResult run() const;
  /// Run at a specified ambient pressure (Pa).
  SteadyEngineResult runAt(double ambient_pressure) const;
  /// Run with the configuration overridden -- used by sweeps and Monte Carlo.
  SteadyEngineResult runWith(const EngineConfig& override_config) const;

  const EngineConfig& config() const { return config_; }
  const SpeciesDatabase& database() const { return *db_; }
  const EquilibriumSolver& solver() const { return *solver_; }
  const TransportModel& transport() const { return *transport_; }
  const PropellantLibrary& propellants() const { return *propellants_; }

  /// Build the propellant mixture described by a configuration.
  PropellantMixture mixture(const EngineConfig& cfg) const;

 private:
  SteadyEngineResult analyse(const EngineConfig& cfg, double ambient_pressure) const;

  EngineConfig config_;
  std::shared_ptr<SpeciesDatabase> db_;
  std::shared_ptr<EquilibriumSolver> solver_;
  std::shared_ptr<TransportModel> transport_;
  std::shared_ptr<PropellantLibrary> propellants_;
};

}  // namespace ignis
