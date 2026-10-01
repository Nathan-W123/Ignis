// SPDX-License-Identifier: MIT
#pragma once
/// \file Config.hpp
/// \brief YAML configuration parsing with path-aware validation.
///
/// Every Ignis application is driven by a YAML file.  The parser is strict:
/// unknown keys are an error (which catches typos before they silently change
/// nothing), required keys are named when missing, and every value is range
/// checked at load time rather than deep inside a solver.
///
/// All physical quantities in configuration files are SI, matching the library
/// convention documented in core/Constants.hpp.  The one exception is angles,
/// which are given in degrees and named accordingly.

#include <memory>
#include <string>
#include <vector>

#include <yaml-cpp/yaml.h>

#include "ignis/combustion/Propellant.hpp"
#include "ignis/cycle/Cycle.hpp"
#include "ignis/cycle/FeedSystem.hpp"
#include "ignis/kinetics/NozzleKinetics.hpp"
#include "ignis/nozzle/NozzleFlow.hpp"
#include "ignis/nozzle/NozzleGeometry.hpp"
#include "ignis/thermal/RegenerativeCooling.hpp"
#include "ignis/thermo/GasState.hpp"
#include "ignis/transient/TransientChamber.hpp"

namespace ignis {

/// Wrapper over a YAML node that keeps track of its path for error messages.
class ConfigNode {
 public:
  ConfigNode() = default;
  ConfigNode(YAML::Node node, std::string path, std::string file)
      : node_(std::move(node)), path_(std::move(path)), file_(std::move(file)) {}

  bool exists() const { return node_ && !node_.IsNull(); }
  const std::string& path() const { return path_; }

  /// Child node; throws ConfigError when absent.
  ConfigNode operator[](const std::string& key) const;
  /// Child node; returns a non-existing node when absent.
  ConfigNode optional(const std::string& key) const;
  bool has(const std::string& key) const;

  double number() const;
  double number(double fallback) const;
  /// Number restricted to [lo, hi]; the error names the offending key.
  double number(double lo, double hi, double fallback) const;
  int integer() const;
  int integer(int fallback) const;
  bool boolean(bool fallback) const;
  std::string text() const;
  std::string text(const std::string& fallback) const;
  std::vector<double> numbers() const;
  std::vector<std::string> strings() const;
  std::vector<std::string> strings(const std::vector<std::string>& fallback) const;

  /// Reject any key not in `allowed`, listing the valid ones.  Call this on
  /// every mapping so that a misspelt option is loud rather than ignored.
  void requireOnly(const std::vector<std::string>& allowed) const;

  const YAML::Node& raw() const { return node_; }

 private:
  [[noreturn]] void fail(const std::string& what) const;

  YAML::Node node_;
  std::string path_;
  std::string file_;
};

/// Everything a steady engine analysis needs.
struct EngineConfig {
  std::string file;
  std::string name = "engine";
  std::string description;

  // species selection
  std::string species_database;             ///< empty => shipped default
  std::vector<std::string> elements = {"C", "H", "O"};
  std::vector<std::string> species_include; ///< explicit list wins over elements

  // propellants
  std::string propellant_library;
  std::string oxidizer = "LOX";
  std::string fuel = "LCH4";
  double oxidizer_temperature = 0.0;        ///< 0 => the library reference value
  double fuel_temperature = 0.0;
  double mixture_ratio = 3.4;

  // chamber
  double chamber_pressure = 1.0e7;
  double eta_c_star = 1.0;
  CompositionModel composition = CompositionModel::kEquilibrium;

  // geometry and performance
  NozzleGeometrySpec nozzle;
  NozzlePerformanceOptions performance;
  /// March the wall boundary layer and correct the nozzle for it (discharge
  /// coefficient, displaced exit, momentum deficit).  On by default.
  bool boundary_layer_losses = true;
  /// Wall temperature the boundary layer assumes where no cooling jacket
  /// supplies one, K.
  double uncooled_wall_temperature = 1000.0;
  /// Finite-rate recombination through the supersonic nozzle
  /// (performance.kinetics).  With `kinetics_apply` the delivered thrust
  /// carries the kinetic efficiency; otherwise it is reported alongside.
  bool kinetics_enabled = false;
  bool kinetics_apply = true;
  std::string kinetics_mechanism;   ///< empty = the shipped GRI-Mech 3.0 subset
  KineticNozzleOptions kinetics;
  bool use_altitude = false;
  double altitude = 0.0;
  double ambient_pressure = 101325.0;

  // cooling
  bool cooling_enabled = false;
  CoolingSpec cooling;
  std::string wall_material = "CuCrZr";
  double coolant_fuel_fraction = 1.0;       ///< fraction of the fuel flow used
  /// Wall film as a fraction of the fuel flow (cooling.film.fuel_fraction);
  /// the engine turns it into cooling.film_mass_flow.  0 = no film unless
  /// cooling.film.mass_flow is given directly.
  double film_fuel_fraction = 0.0;

  // feed system
  bool feed_enabled = false;
  FeedSystemSpec feed;
  double injector_stiffness = 0.2;
  bool feed_capability_mode = false;

  /// Turbopump cycle (pump-fed engines).  When the jacket's inlet pressure
  /// is not configured, the cycle sets it: it is the fuel pump's discharge.
  bool cycle_enabled = false;
  CycleSpec cycle;

  /// Sample the axial flow profile during a steady analysis.  Sweeps, Monte
  /// Carlo and the optimiser turn this off: the profile costs a nozzle solve at
  /// every contour station and none of them read it.
  bool sample_profile = true;

  // output
  std::string output_directory = "results";
  std::string output_prefix;

  /// The raw document, so that sweep / Monte Carlo / optimisation sections can
  /// be read by their own modules.
  std::shared_ptr<YAML::Node> document;
  ConfigNode root() const;

  static EngineConfig load(const std::string& path);
  /// Resolve the ambient pressure, applying the atmosphere model if an
  /// altitude was given.
  double resolveAmbientPressure() const;
};

/// Parse a transient section (returns false when the section is absent).
struct TransientConfig {
  TransientSpec spec;
  TableGrid grid;
  int interpolation_samples = 200;
  unsigned interpolation_seed = 20260917u;
  bool enabled = false;
};
TransientConfig parseTransient(const ConfigNode& root, const EngineConfig& engine);

}  // namespace ignis
