// SPDX-License-Identifier: MIT
#include "ignis/thermal/HeatTransfer.hpp"

#include <yaml-cpp/yaml.h>

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <sstream>

#include "ignis/core/Version.hpp"

namespace ignis {

double WallMaterial::conductivity(double T) const {
  const double k = conductivity_reference + conductivity_slope * (T - reference_temperature);
  if (!(k > 0.0)) {
    std::ostringstream os;
    os << "material " << name << ": the linear conductivity model gives k = " << k
       << " W/(m K) at " << T << " K, which is unphysical";
    throw RangeError(os.str());
  }
  return k;
}

MaterialLibrary MaterialLibrary::loadYaml(const std::string& path) {
  YAML::Node root;
  try {
    root = YAML::LoadFile(path);
  } catch (const YAML::Exception& e) {
    throw ConfigError("cannot parse material library '" + path + "': " + e.what());
  }
  MaterialLibrary lib;
  for (const auto& node : root["materials"]) {
    WallMaterial m;
    m.name = node["name"].as<std::string>();
    m.conductivity_reference = node["conductivity-reference"].as<double>();
    m.reference_temperature = node["reference-temperature"].as<double>();
    m.conductivity_slope = node["conductivity-slope"].as<double>(0.0);
    const auto vr = node["valid-range"].as<std::vector<double>>(std::vector<double>{250.0, 1000.0});
    m.valid_min = vr.at(0);
    m.valid_max = vr.at(1);
    m.max_temperature = node["max-temperature"].as<double>();
    m.emissivity = node["emissivity"].as<double>(0.3);
    m.density = node["density"].as<double>(8000.0);
    m.source = node["source"].as<std::string>("");
    lib.materials_.push_back(m);
  }
  if (lib.materials_.empty()) throw ConfigError("material library '" + path + "' is empty");
  return lib;
}

MaterialLibrary MaterialLibrary::loadDefault() {
  namespace fs = std::filesystem;
  std::vector<std::string> tried;
  auto probe = [&](const fs::path& p) {
    tried.push_back(p.string());
    std::error_code ec;
    return fs::exists(p, ec);
  };
  const char* env = std::getenv("IGNIS_DATA_DIR");
  if (env != nullptr) {
    const fs::path p = fs::path(env) / "materials" / "ignis_materials.yaml";
    if (probe(p)) return loadYaml(p.string());
  }
  {
    const fs::path p = fs::path(kDefaultDataDir) / "materials" / "ignis_materials.yaml";
    if (probe(p)) return loadYaml(p.string());
  }
  for (const char* rel : {"data/materials/ignis_materials.yaml",
                          "../data/materials/ignis_materials.yaml",
                          "../../data/materials/ignis_materials.yaml"}) {
    if (probe(rel)) return loadYaml(rel);
  }
  std::ostringstream os;
  os << "cannot locate the Ignis material library. Tried:";
  for (const auto& t : tried) os << "\n  " << t;
  throw ConfigError(os.str());
}

const WallMaterial& MaterialLibrary::at(const std::string& name) const {
  for (const auto& m : materials_)
    if (m.name == name) return m;
  std::ostringstream os;
  os << "unknown wall material '" << name << "'. Available:";
  for (const auto& m : materials_) os << " " << m.name;
  throw ConfigError(os.str());
}

std::vector<std::string> MaterialLibrary::names() const {
  std::vector<std::string> n;
  for (const auto& m : materials_) n.push_back(m.name);
  return n;
}

double bartzFilmCoefficient(const BartzReference& ref, double area_ratio, double mach,
                            double gamma, double t_wall_hot) {
  if (!(ref.throat_diameter > 0.0)) throw ConfigError("Bartz: throat diameter must be positive");
  if (!(ref.curvature_radius > 0.0)) throw ConfigError("Bartz: curvature radius must be positive");
  if (!(ref.c_star > 0.0)) throw ConfigError("Bartz: c* must be positive");
  if (!(ref.prandtl > 0.0)) throw ConfigError("Bartz: Prandtl number must be positive");
  if (!(ref.chamber_temperature > 0.0)) throw ConfigError("Bartz: chamber temperature must be positive");
  if (!(area_ratio >= 1.0)) throw ConfigError("Bartz: local area ratio must be at least 1");
  if (!(t_wall_hot > 0.0)) throw ConfigError("Bartz: wall temperature must be positive");

  const double m2 = 0.5 * (gamma - 1.0) * mach * mach;
  const double sigma = 1.0 / (std::pow(0.5 * (t_wall_hot / ref.chamber_temperature) *
                                           (1.0 + m2) + 0.5, 0.68) *
                              std::pow(1.0 + m2, 0.12));
  return ref.multiplier * (0.026 / std::pow(ref.throat_diameter, 0.2)) *
         (std::pow(ref.viscosity, 0.2) * ref.cp / std::pow(ref.prandtl, 0.6)) *
         std::pow(ref.chamber_pressure / ref.c_star, 0.8) *
         std::pow(ref.throat_diameter / ref.curvature_radius, 0.1) *
         std::pow(1.0 / area_ratio, 0.9) * sigma;
}

double recoveryTemperature(double t_static, double mach, double gamma, double prandtl) {
  if (!(prandtl > 0.0)) throw ConfigError("recovery temperature: Prandtl number must be positive");
  const double r = std::cbrt(prandtl);
  return t_static * (1.0 + r * 0.5 * (gamma - 1.0) * mach * mach);
}

double cylindricalWallResistance(double inner_radius, double thickness, double conductivity) {
  if (!(inner_radius > 0.0)) throw ConfigError("wall conduction: inner radius must be positive");
  if (!(thickness > 0.0)) throw ConfigError("wall conduction: thickness must be positive");
  if (!(conductivity > 0.0)) throw ConfigError("wall conduction: conductivity must be positive");
  return inner_radius * std::log((inner_radius + thickness) / inner_radius) / conductivity;
}

double grayGasRadiation(double t_gas, double t_wall, double gas_emissivity,
                        double wall_emissivity) {
  if (gas_emissivity <= 0.0) return 0.0;
  if (!(gas_emissivity <= 1.0) || !(wall_emissivity > 0.0 && wall_emissivity <= 1.0))
    throw ConfigError("radiation: emissivities must lie in (0, 1]");
  const double tg2 = t_gas * t_gas, tw2 = t_wall * t_wall;
  return constants::sigma_SB * gas_emissivity * wall_emissivity * (tg2 * tg2 - tw2 * tw2);
}

}  // namespace ignis
