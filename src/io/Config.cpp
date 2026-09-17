// SPDX-License-Identifier: MIT
#include "ignis/io/Config.hpp"

#include <algorithm>
#include <cmath>
#include <sstream>

#include "ignis/nozzle/Atmosphere.hpp"
#include "ignis/thermal/HeatTransfer.hpp"

namespace ignis {

void ConfigNode::fail(const std::string& what) const {
  std::ostringstream os;
  os << file_ << ": at '" << (path_.empty() ? "<root>" : path_) << "': " << what;
  throw ConfigError(os.str());
}

ConfigNode ConfigNode::operator[](const std::string& key) const {
  if (!node_ || !node_.IsMap()) fail("expected a mapping while looking up '" + key + "'");
  const auto child = node_[key];
  if (!child || child.IsNull()) fail("required key '" + key + "' is missing");
  return ConfigNode(child, path_.empty() ? key : path_ + "." + key, file_);
}

ConfigNode ConfigNode::optional(const std::string& key) const {
  if (!node_ || !node_.IsMap()) return ConfigNode({}, path_.empty() ? key : path_ + "." + key, file_);
  return ConfigNode(node_[key], path_.empty() ? key : path_ + "." + key, file_);
}

bool ConfigNode::has(const std::string& key) const {
  return node_ && node_.IsMap() && node_[key] && !node_[key].IsNull();
}

double ConfigNode::number() const {
  if (!exists()) fail("expected a number");
  try {
    return node_.as<double>();
  } catch (const YAML::Exception&) {
    fail("expected a number, found '" + node_.as<std::string>("?") + "'");
  }
}

double ConfigNode::number(double fallback) const { return exists() ? number() : fallback; }

double ConfigNode::number(double lo, double hi, double fallback) const {
  const double v = number(fallback);
  if (!(v >= lo && v <= hi)) {
    std::ostringstream os;
    os << "value " << v << " is outside the allowed range [" << lo << ", " << hi << "]";
    fail(os.str());
  }
  return v;
}

int ConfigNode::integer() const {
  if (!exists()) fail("expected an integer");
  try {
    return node_.as<int>();
  } catch (const YAML::Exception&) {
    fail("expected an integer");
  }
}

int ConfigNode::integer(int fallback) const { return exists() ? integer() : fallback; }

bool ConfigNode::boolean(bool fallback) const {
  if (!exists()) return fallback;
  try {
    return node_.as<bool>();
  } catch (const YAML::Exception&) {
    fail("expected true or false");
  }
}

std::string ConfigNode::text() const {
  if (!exists()) fail("expected a string");
  return node_.as<std::string>();
}

std::string ConfigNode::text(const std::string& fallback) const {
  return exists() ? text() : fallback;
}

std::vector<double> ConfigNode::numbers() const {
  if (!exists() || !node_.IsSequence()) fail("expected a sequence of numbers");
  std::vector<double> out;
  for (const auto& e : node_) out.push_back(e.as<double>());
  return out;
}

std::vector<std::string> ConfigNode::strings() const {
  if (!exists() || !node_.IsSequence()) fail("expected a sequence of strings");
  std::vector<std::string> out;
  for (const auto& e : node_) out.push_back(e.as<std::string>());
  return out;
}

std::vector<std::string> ConfigNode::strings(const std::vector<std::string>& fallback) const {
  return exists() ? strings() : fallback;
}

void ConfigNode::requireOnly(const std::vector<std::string>& allowed) const {
  if (!exists() || !node_.IsMap()) return;
  std::vector<std::string> unknown;
  for (const auto& kv : node_) {
    const auto key = kv.first.as<std::string>();
    if (std::find(allowed.begin(), allowed.end(), key) == allowed.end()) unknown.push_back(key);
  }
  if (unknown.empty()) return;
  std::ostringstream os;
  os << "unknown key" << (unknown.size() > 1 ? "s" : "") << ":";
  for (const auto& u : unknown) os << " '" << u << "'";
  os << ". Recognised keys are:";
  for (const auto& a : allowed) os << " " << a;
  fail(os.str());
}

ConfigNode EngineConfig::root() const {
  if (!document) throw ConfigError("configuration document was not loaded");
  return ConfigNode(*document, "", file);
}

double EngineConfig::resolveAmbientPressure() const {
  if (use_altitude) return Atmosphere::at(altitude).pressure;
  return ambient_pressure;
}

EngineConfig EngineConfig::load(const std::string& path) {
  YAML::Node doc;
  try {
    doc = YAML::LoadFile(path);
  } catch (const YAML::BadFile&) {
    throw ConfigError("cannot open configuration file '" + path + "'");
  } catch (const YAML::Exception& e) {
    throw ConfigError("cannot parse configuration file '" + path + "': " + e.what());
  }

  EngineConfig cfg;
  cfg.file = path;
  cfg.document = std::make_shared<YAML::Node>(doc);
  ConfigNode root(doc, "", path);
  root.requireOnly({"name", "description", "species", "propellants", "chamber", "nozzle",
                    "performance", "cooling", "feed", "transient", "sweep", "optimization",
                    "monte_carlo", "output"});

  cfg.name = root.optional("name").text("engine");
  cfg.description = root.optional("description").text("");

  // --- species ---------------------------------------------------------
  if (root.has("species")) {
    const auto s = root["species"];
    s.requireOnly({"database", "elements", "include"});
    cfg.species_database = s.optional("database").text("");
    cfg.elements = s.optional("elements").strings(cfg.elements);
    cfg.species_include = s.optional("include").strings({});
  }

  // --- propellants -----------------------------------------------------
  {
    const auto p = root["propellants"];
    p.requireOnly({"library", "oxidizer", "fuel", "oxidizer_temperature", "fuel_temperature",
                   "mixture_ratio"});
    cfg.propellant_library = p.optional("library").text("");
    cfg.oxidizer = p["oxidizer"].text();
    cfg.fuel = p["fuel"].text();
    cfg.oxidizer_temperature = p.optional("oxidizer_temperature").number(0.0);
    cfg.fuel_temperature = p.optional("fuel_temperature").number(0.0);
    cfg.mixture_ratio = p["mixture_ratio"].number(1e-3, 1e3, 0.0);
  }

  // --- chamber ---------------------------------------------------------
  {
    const auto c = root["chamber"];
    c.requireOnly({"pressure", "eta_c_star", "composition"});
    cfg.chamber_pressure = c["pressure"].number(1.0e3, 1.0e9, 0.0);
    cfg.eta_c_star = c.optional("eta_c_star").number(0.1, 1.0, 1.0);
    cfg.composition = compositionModelFromString(c.optional("composition").text("equilibrium"));
  }

  // --- nozzle ----------------------------------------------------------
  {
    const auto n = root["nozzle"];
    n.requireOnly({"throat_radius", "throat_area", "chamber_radius", "contraction_ratio",
                   "chamber_length", "converging_half_angle", "chamber_fillet_ratio",
                   "throat_upstream_ratio", "throat_downstream_ratio", "expansion_ratio",
                   "type", "cone_half_angle", "bell_length_fraction", "bell_initial_angle",
                   "bell_exit_angle", "stations", "bartz_curvature"});
    auto& g = cfg.nozzle;
    g.throat_radius = n.optional("throat_radius").number(0.0);
    g.throat_area = n.optional("throat_area").number(0.0);
    g.chamber_radius = n.optional("chamber_radius").number(0.0);
    g.contraction_ratio = n.optional("contraction_ratio").number(0.0);
    g.chamber_length = n["chamber_length"].number(0.0, 100.0, 0.0);
    g.converging_half_angle = n.optional("converging_half_angle").number(2.0, 80.0, 30.0);
    g.chamber_fillet_ratio = n.optional("chamber_fillet_ratio").number(0.05, 5.0, 0.5);
    g.throat_upstream_ratio = n.optional("throat_upstream_ratio").number(0.1, 10.0, 1.5);
    g.throat_downstream_ratio = n.optional("throat_downstream_ratio").number(0.05, 10.0, 0.382);
    g.expansion_ratio = n["expansion_ratio"].number(1.0001, 1000.0, 0.0);
    g.divergent = divergentTypeFromString(n.optional("type").text("bell"));
    g.cone_half_angle = n.optional("cone_half_angle").number(1.0, 45.0, 15.0);
    g.bell_length_fraction = n.optional("bell_length_fraction").number(0.3, 1.2, 0.8);
    g.bell_initial_angle = n.optional("bell_initial_angle").number(0.0, 60.0, 0.0);
    g.bell_exit_angle = n.optional("bell_exit_angle").number(0.0, 45.0, 0.0);
    g.num_stations = n.optional("stations").integer(400);
    g.bartz_curvature = n.optional("bartz_curvature").text("mean");
  }

  // --- performance -----------------------------------------------------
  if (root.has("performance")) {
    const auto p = root["performance"];
    p.requireOnly({"ambient_pressure", "altitude", "auto_divergence", "lambda_divergence",
                   "eta_nozzle", "separation", "resolve_internal_shocks", "ideal_tolerance",
                   "ascent_profile"});
    if (p.has("altitude") && p.has("ambient_pressure"))
      throw ConfigError(path + ": performance: give either 'altitude' or 'ambient_pressure'");
    if (p.has("altitude")) {
      cfg.use_altitude = true;
      cfg.altitude = p["altitude"].number(-5000.0, 1.0e6, 0.0);
    } else {
      cfg.ambient_pressure = p.optional("ambient_pressure").number(0.0, 1.0e7, 101325.0);
    }
    cfg.performance.auto_divergence = p.optional("auto_divergence").boolean(true);
    cfg.performance.lambda_divergence = p.optional("lambda_divergence").number(0.5, 1.0, 1.0);
    cfg.performance.eta_nozzle = p.optional("eta_nozzle").number(0.5, 1.0, 1.0);
    cfg.performance.separation =
        separationCriterionFromString(p.optional("separation").text("summerfield"));
    cfg.performance.resolve_internal_shocks =
        p.optional("resolve_internal_shocks").boolean(true);
    cfg.performance.ideal_tolerance = p.optional("ideal_tolerance").number(1e-6, 0.5, 0.01);
    if (p.has("ascent_profile")) {
      const auto a = p["ascent_profile"];
      a.requireOnly({"altitudes", "weights"});
      cfg.performance.ascent_altitudes = a["altitudes"].numbers();
      cfg.performance.ascent_weights = a["weights"].numbers();
      if (cfg.performance.ascent_altitudes.empty())
        throw ConfigError(a.path() + ": the ascent profile needs at least one altitude");
      if (cfg.performance.ascent_weights.size() != cfg.performance.ascent_altitudes.size())
        throw ConfigError(a.path() + ": 'altitudes' and 'weights' must be the same length (" +
                          std::to_string(cfg.performance.ascent_altitudes.size()) + " vs " +
                          std::to_string(cfg.performance.ascent_weights.size()) + ")");
      double w_sum = 0.0;
      for (std::size_t i = 0; i < cfg.performance.ascent_weights.size(); ++i) {
        if (cfg.performance.ascent_weights[i] < 0.0)
          throw ConfigError(a.path() + ": weights must not be negative");
        if (cfg.performance.ascent_altitudes[i] < -5000.0 ||
            cfg.performance.ascent_altitudes[i] > 86000.0)
          throw ConfigError(a.path() + ": altitudes must lie within the U.S. Standard "
                                       "Atmosphere 1976 range [-5000, 86000] m");
        w_sum += cfg.performance.ascent_weights[i];
      }
      if (w_sum <= 0.0) throw ConfigError(a.path() + ": the weights sum to zero");
    }
  }
  cfg.performance.eta_c_star = cfg.eta_c_star;

  // --- cooling ---------------------------------------------------------
  if (root.has("cooling")) {
    const auto c = root["cooling"];
    c.requireOnly({"enabled", "coolant", "num_channels", "width_mode", "width_fraction",
                   "channel_width", "channel_height", "min_land_width", "wall_thickness",
                   "roughness", "coolant_fuel_fraction", "inlet_temperature", "inlet_pressure",
                   "counterflow", "x_start", "x_end", "x_end_area_ratio", "nusselt_correlation",
                   "nusselt_multiplier", "bartz_multiplier", "gas_emissivity", "material",
                   "segments"});
    cfg.cooling_enabled = c.optional("enabled").boolean(true);
    auto& s = cfg.cooling;
    s.coolant = c.optional("coolant").text("methane");
    s.num_channels = c.optional("num_channels").integer(100);
    s.width_mode = channelWidthModeFromString(c.optional("width_mode").text("fraction_of_pitch"));
    s.width_fraction = c.optional("width_fraction").number(0.05, 0.95, 0.5);
    s.channel_width = c.optional("channel_width").number(0.0);
    s.channel_height = c.optional("channel_height").number(1e-5, 0.1, 3.0e-3);
    s.min_land_width = c.optional("min_land_width").number(1e-5, 0.05, 5.0e-4);
    s.wall_thickness = c.optional("wall_thickness").number(1e-5, 0.05, 1.0e-3);
    s.roughness = c.optional("roughness").number(0.0, 1e-3, 5.0e-6);
    s.inlet_temperature = c.optional("inlet_temperature").number(1.0, 2000.0, 0.0);
    s.inlet_pressure = c.optional("inlet_pressure").number(1e4, 1e9, 0.0);
    s.counterflow = c.optional("counterflow").boolean(true);
    s.x_start = c.optional("x_start").number(0.0);
    s.x_end = c.optional("x_end").number(-1.0);
    s.x_end_area_ratio = c.optional("x_end_area_ratio").number(0.0, 1000.0, 0.0);
    s.nusselt_correlation = c.optional("nusselt_correlation").text("dittus-boelter");
    s.nusselt_multiplier = c.optional("nusselt_multiplier").number(0.1, 10.0, 1.0);
    s.bartz_multiplier = c.optional("bartz_multiplier").number(0.1, 10.0, 1.0);
    s.gas_emissivity = c.optional("gas_emissivity").number(0.0, 1.0, 0.0);
    s.num_segments = c.optional("segments").integer(200);
    cfg.wall_material = c.optional("material").text("CuCrZr");
    cfg.coolant_fuel_fraction = c.optional("coolant_fuel_fraction").number(0.01, 1.0, 1.0);
    s.wall = MaterialLibrary::loadDefault().at(cfg.wall_material);
  }

  // --- feed system -----------------------------------------------------
  if (root.has("feed")) {
    const auto f = root["feed"];
    f.requireOnly({"enabled", "mode", "injector_stiffness", "minimum_stiffness",
                   "pump_efficiency", "pump_inlet_pressure", "oxidizer", "fuel"});
    cfg.feed_enabled = f.optional("enabled").boolean(true);
    cfg.injector_stiffness = f.optional("injector_stiffness").number(0.01, 0.9, 0.2);
    cfg.feed.minimum_stiffness = f.optional("minimum_stiffness").number(0.0, 0.9, 0.15);
    cfg.feed.pump_efficiency = f.optional("pump_efficiency").number(0.05, 1.0, 0.65);
    cfg.feed.pump_inlet_pressure = f.optional("pump_inlet_pressure").number(1e3, 1e8, 3.0e5);
    cfg.feed_capability_mode = (f.optional("mode").text("sizing") == "capability");
    auto leg = [&](const ConfigNode& n, FeedLeg& l) {
      n.requireOnly({"line_length", "line_diameter", "line_roughness", "fitting_k", "viscosity",
                     "injector_cd", "injector_area", "tank_pressure"});
      l.line_length = n.optional("line_length").number(0.0, 100.0, 1.0);
      l.line_diameter = n.optional("line_diameter").number(1e-4, 1.0, 0.05);
      l.line_roughness = n.optional("line_roughness").number(0.0, 1e-3, 5.0e-6);
      l.fitting_k = n.optional("fitting_k").number(0.0, 100.0, 2.0);
      l.viscosity = n.optional("viscosity").number(1e-7, 1.0, 1.0e-4);
      l.injector_cd = n.optional("injector_cd").number(0.05, 1.0, 0.75);
      l.injector_area = n.optional("injector_area").number(0.0);
      l.tank_pressure = n.optional("tank_pressure").number(0.0);
    };
    if (f.has("oxidizer")) leg(f["oxidizer"], cfg.feed.oxidizer);
    if (f.has("fuel")) leg(f["fuel"], cfg.feed.fuel);
  }

  // --- output ----------------------------------------------------------
  if (root.has("output")) {
    const auto o = root["output"];
    o.requireOnly({"directory", "prefix"});
    cfg.output_directory = o.optional("directory").text("results");
    cfg.output_prefix = o.optional("prefix").text("");
  }
  if (cfg.output_prefix.empty()) {
    cfg.output_prefix = cfg.name;
    for (auto& ch : cfg.output_prefix)
      if (!std::isalnum(static_cast<unsigned char>(ch))) ch = '_';
  }
  return cfg;
}

TransientConfig parseTransient(const ConfigNode& root, const EngineConfig& engine) {
  TransientConfig out;
  if (!root.has("transient")) return out;
  const auto t = root["transient"];
  t.requireOnly({"enabled", "chamber_volume", "throat_area", "exit_area", "ambient_pressure",
                 "altitude", "oxidizer_schedule", "fuel_schedule", "ignition", "initial_pressure",
                 "initial_temperature", "initial_mixture_ratio", "wall_heat_rate", "t_end", "dt",
                 "dt_output", "integrator", "rtol", "atol", "dt_min", "dt_max", "table",
                 "interpolation_samples", "interpolation_seed"});
  out.enabled = t.optional("enabled").boolean(true);
  auto& s = out.spec;
  s.chamber_volume = t.optional("chamber_volume").number(0.0);
  s.throat_area = t.optional("throat_area").number(0.0);
  s.exit_area = t.optional("exit_area").number(0.0);
  if (t.has("altitude")) {
    s.ambient_pressure = Atmosphere::at(t["altitude"].number(-5000.0, 1e6, 0.0)).pressure;
  } else {
    s.ambient_pressure = t.optional("ambient_pressure").number(0.0, 1e7, 101325.0);
  }
  auto schedule = [&](const ConfigNode& n, const char* what) {
    n.requireOnly({"open", "ramp_up", "value", "close", "ramp_down", "times", "values"});
    if (n.has("times")) {
      const auto times = n["times"].numbers();
      const auto values = n["values"].numbers();
      if (times.size() != values.size())
        throw ConfigError(std::string("transient: ") + what +
                          " schedule needs equal-length 'times' and 'values'");
      return Schedule(times, values);
    }
    return Schedule::ramp(n["open"].number(0.0, 1e4, 0.0), n["ramp_up"].number(1e-9, 1e4, 0.0),
                          n["value"].number(0.0, 1e6, 0.0), n["close"].number(0.0, 1e4, 0.0),
                          n["ramp_down"].number(1e-9, 1e4, 0.0));
  };
  s.oxidizer_flow = schedule(t["oxidizer_schedule"], "oxidizer");
  s.fuel_flow = schedule(t["fuel_schedule"], "fuel");
  {
    const auto ig = t["ignition"];
    ig.requireOnly({"times", "values"});
    s.combustion_efficiency = Schedule(ig["times"].numbers(), ig["values"].numbers());
  }
  s.initial_pressure = t.optional("initial_pressure").number(1.0, 1e8, 101325.0);
  s.initial_temperature = t.optional("initial_temperature").number(100.0, 4000.0, 300.0);
  s.initial_mixture_ratio =
      t.optional("initial_mixture_ratio").number(0.0, 100.0, engine.mixture_ratio);
  s.wall_heat_rate = t.optional("wall_heat_rate").number(-1e12, 1e12, 0.0);
  s.t_end = t["t_end"].number(1e-6, 1e4, 0.0);
  s.dt = t.optional("dt").number(1e-12, 1.0, 1.0e-5);
  s.dt_output = t.optional("dt_output").number(1e-9, 1.0, 1.0e-4);
  s.integrator = t.optional("integrator").text("rk45");
  s.rtol = t.optional("rtol").number(1e-14, 1e-2, 1.0e-8);
  s.atol = t.optional("atol").number(1e-18, 1e-2, 1.0e-10);
  s.dt_min = t.optional("dt_min").number(1e-14, 1.0, 1.0e-10);
  s.dt_max = t.optional("dt_max").number(1e-9, 1.0, 1.0e-3);

  if (t.has("table")) {
    const auto g = t["table"];
    g.requireOnly({"mr_min", "mr_max", "mr_points", "t_min", "t_max", "t_points", "p_min",
                   "p_max", "p_points"});
    out.grid.mr_min = g.optional("mr_min").number(0.05, 100.0, 1.0);
    out.grid.mr_max = g.optional("mr_max").number(0.05, 200.0, 8.0);
    out.grid.mr_points = g.optional("mr_points").integer(41);
    out.grid.t_min = g.optional("t_min").number(200.0, 3000.0, 200.0);
    out.grid.t_max = g.optional("t_max").number(1000.0, 6000.0, 4000.0);
    out.grid.t_points = g.optional("t_points").integer(97);
    out.grid.p_min = g.optional("p_min").number(1.0, 1e9, 5.0e4);
    out.grid.p_max = g.optional("p_max").number(1e3, 1e9, 2.0e7);
    out.grid.p_points = g.optional("p_points").integer(13);
  }
  out.interpolation_samples = t.optional("interpolation_samples").integer(200);
  out.interpolation_seed =
      static_cast<unsigned>(t.optional("interpolation_seed").integer(20260917));
  return out;
}

}  // namespace ignis
