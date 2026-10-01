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
///   5. mass flows   ->  turbopump cycle (which sets the jacket inlet pressure)
///                   ->  feed system
///
/// Nothing here introduces new physics; it only wires the modules together and
/// carries the results, so an application, a sweep, an optimiser and a Monte
/// Carlo sample all run exactly the same analysis.

#include <memory>
#include <string>

#include "ignis/combustion/Chamber.hpp"
#include "ignis/cycle/Cycle.hpp"
#include "ignis/cycle/FeedSystem.hpp"
#include "ignis/io/Config.hpp"
#include "ignis/io/Json.hpp"
#include "ignis/io/Table.hpp"
#include "ignis/nozzle/Atmosphere.hpp"
#include "ignis/thermal/BoundaryLayerLosses.hpp"
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
  /// Installed-envelope geometry.  A nozzle that does not fit the vehicle is
  /// not a design, however good its specific impulse, so these are exposed as
  /// optimisable metrics rather than left inside NozzleGeometry.
  double exit_radius = 0.0;       ///< m
  double exit_diameter = 0.0;     ///< m
  double total_length = 0.0;      ///< m, injector face to exit plane
  double divergent_length = 0.0;  ///< m, throat to exit plane
  double ambient_pressure = 0.0;
  double altitude = 0.0;
  bool altitude_known = false;

  bool has_cooling = false;
  CoolingResult cooling;
  /// The wall boundary layer and what it cost (applied to `performance`).
  bool has_boundary_layer = false;
  BoundaryLayerLossResult boundary_layer;

  /// What a wall film costs in specific impulse, as a bracket.  The headline
  /// `performance` assumes the film mixes and burns completely, i.e. the
  /// engine at its overall mixture ratio.  The other limit keeps the film as a
  /// separate unburnt stream: the core burns at the mixture ratio the film
  /// leaves it, the film expands as a calorically perfect gas from its
  /// injection state to the core's exit pressure, and the two are
  /// mass-weighted in vacuum.  A real film partly mixes and burns, so the
  /// truth lies between; neither limit predicts where.
  struct FilmPerformance {
    bool present = false;
    double film_mass_flow = 0.0;          ///< kg/s
    double film_fraction_of_total = 0.0;  ///< of the engine's mass flow
    double core_mixture_ratio = 0.0;      ///< O/F of the core without the film
    double film_exhaust_velocity = 0.0;   ///< m/s, the unburnt film expanded alone
    double isp_vacuum_mixed = 0.0;        ///< s, the headline (fully mixed) value
    double isp_vacuum_unmixed = 0.0;      ///< s, the two-stream limit
    double isp_vacuum_penalty = 0.0;      ///< 1 - unmixed / mixed
  } film;
  bool has_feed = false;
  FeedSystemResult feed;
  /// The turbopump cycle, and how it set the jacket's inlet pressure.
  bool has_cycle = false;
  CycleResult cycle;
  struct CycleCoupling {
    bool jacket_inlet_from_cycle = false;  ///< the cycle, not the config, set it
    int passes = 0;                        ///< cooling solves while iterating
    double jacket_inlet_pressure = 0.0;    ///< Pa, the jacket's final solve
    double mismatch = 0.0;                 ///< final fuel discharge / that - 1
  } cycle_coupling;

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
