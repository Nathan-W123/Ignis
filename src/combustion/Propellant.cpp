// SPDX-License-Identifier: MIT
#include "ignis/combustion/Propellant.hpp"

#include <yaml-cpp/yaml.h>

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <sstream>

#include "ignis/core/Version.hpp"

namespace ignis {

double Propellant::molarEnthalpy(double T_inlet, const SpeciesDatabase& db) const {
  if (!liquid) {
    if (species.empty())
      throw ConfigError("propellant " + name + ": gaseous entry has no 'species' field");
    return db.at(species).h(T_inlet);
  }
  const double dT = T_inlet - reference_temperature;
  if (std::abs(dT) > max_temperature_offset) {
    std::ostringstream os;
    os << "propellant " << name << ": inlet temperature " << T_inlet << " K is "
       << std::abs(dT) << " K from the tabulated reference state (" << reference_temperature
       << " K); the linear sensible correction is only accepted within "
       << max_temperature_offset << " K";
    throw RangeError(os.str());
  }
  return reference_enthalpy + cp_liquid * dT;
}

PropellantLibrary PropellantLibrary::loadYaml(const std::string& path) {
  YAML::Node root;
  try {
    root = YAML::LoadFile(path);
  } catch (const YAML::Exception& e) {
    throw ConfigError("cannot parse propellant library '" + path + "': " + e.what());
  }
  if (!root["propellants"] || !root["propellants"].IsSequence())
    throw ConfigError("propellant library '" + path + "' has no 'propellants' sequence");

  PropellantLibrary lib;
  lib.path_ = path;
  for (const auto& node : root["propellants"]) {
    Propellant p;
    p.name = node["name"].as<std::string>();
    p.description = node["description"].as<std::string>("");
    p.role = node["role"].as<std::string>("fuel");
    const auto phase = node["phase"].as<std::string>("liquid");
    if (phase != "liquid" && phase != "gas")
      throw ConfigError("propellant " + p.name + ": phase must be 'liquid' or 'gas'");
    p.liquid = (phase == "liquid");
    for (const auto& kv : node["composition"])
      p.composition[kv.first.as<std::string>()] = kv.second.as<int>();
    p.molar_mass = node["molar-mass"].as<double>() * 1.0e-3;
    p.reference_temperature = node["reference-temperature"].as<double>();
    p.reference_enthalpy = node["reference-enthalpy"].as<double>(0.0);
    p.cp_liquid = node["cp-liquid"].as<double>(0.0);
    p.density = node["density"].as<double>(0.0);
    p.species = node["species"].as<std::string>("");
    p.coolant_table = node["coolant-table"].as<std::string>("");
    p.max_temperature_offset = node["max-temperature-offset"].as<double>(20.0);
    p.source = node["source"].as<std::string>("");
    if (p.liquid && !(p.density > 0.0))
      throw ConfigError("propellant " + p.name + ": liquid entries need a positive density");
    if (p.composition.empty())
      throw ConfigError("propellant " + p.name + ": empty composition");
    lib.by_name_[p.name] = p;
  }
  if (lib.by_name_.empty()) throw ConfigError("propellant library '" + path + "' is empty");
  return lib;
}

const Propellant& PropellantLibrary::at(const std::string& name) const {
  auto it = by_name_.find(name);
  if (it == by_name_.end()) {
    std::ostringstream os;
    os << "unknown propellant '" << name << "'. Available:";
    for (const auto& kv : by_name_) os << " " << kv.first;
    throw ConfigError(os.str());
  }
  return it->second;
}

std::vector<std::string> PropellantLibrary::names() const {
  std::vector<std::string> n;
  for (const auto& kv : by_name_) n.push_back(kv.first);
  return n;
}

std::string findDefaultPropellantLibrary() {
  namespace fs = std::filesystem;
  std::vector<std::string> tried;
  auto probe = [&](const fs::path& p) {
    tried.push_back(p.string());
    std::error_code ec;
    return fs::exists(p, ec) && fs::is_regular_file(p, ec);
  };
  const char* env = std::getenv("IGNIS_DATA_DIR");
  if (env != nullptr) {
    const fs::path p = fs::path(env) / "propellants" / "ignis_propellants.yaml";
    if (probe(p)) return p.string();
  }
  {
    const fs::path p = fs::path(kDefaultDataDir) / "propellants" / "ignis_propellants.yaml";
    if (probe(p)) return p.string();
  }
  for (const char* rel : {"data/propellants/ignis_propellants.yaml",
                          "../data/propellants/ignis_propellants.yaml",
                          "../../data/propellants/ignis_propellants.yaml"}) {
    if (probe(rel)) return std::string(rel);
  }
  std::ostringstream os;
  os << "cannot locate the Ignis propellant library. Tried:";
  for (const auto& t : tried) os << "\n  " << t;
  throw ConfigError(os.str());
}

PropellantMixture::PropellantMixture(Propellant oxidizer, Propellant fuel, double mixture_ratio,
                                     double T_oxidizer, double T_fuel)
    : ox_(std::move(oxidizer)),
      fuel_(std::move(fuel)),
      mr_(mixture_ratio),
      T_ox_(T_oxidizer),
      T_fuel_(T_fuel) {
  if (!(mr_ > 0.0) || !std::isfinite(mr_))
    throw ConfigError("mixture ratio must be a finite positive number, got " +
                      std::to_string(mr_));
  if (!(T_ox_ > 0.0) || !(T_fuel_ > 0.0))
    throw ConfigError("propellant inlet temperatures must be positive");
}

Eigen::VectorXd PropellantMixture::elementMoles(const SpeciesDatabase& db) const {
  Eigen::VectorXd b = Eigen::VectorXd::Zero(static_cast<Eigen::Index>(db.numElements()));
  const double n_ox = oxidizerMassFraction() / ox_.molar_mass;    // mol/kg of mixture
  const double n_fu = fuelMassFraction() / fuel_.molar_mass;
  auto add = [&](const Propellant& p, double n) {
    for (const auto& kv : p.composition) {
      const int i = db.elementIndex(kv.first);
      if (i < 0)
        throw ConfigError("propellant " + p.name + " contains element '" + kv.first +
                          "' which no species in the active set can carry");
      b(i) += n * kv.second;
    }
  };
  add(ox_, n_ox);
  add(fuel_, n_fu);
  return b;
}

double PropellantMixture::enthalpy(const SpeciesDatabase& db) const {
  const double n_ox = oxidizerMassFraction() / ox_.molar_mass;
  const double n_fu = fuelMassFraction() / fuel_.molar_mass;
  return n_ox * ox_.molarEnthalpy(T_ox_, db) + n_fu * fuel_.molarEnthalpy(T_fuel_, db);
}

double PropellantMixture::stoichiometricMixtureRatio() const {
  // Oxygen atoms required to take all C to CO2 and all H to H2O, minus the
  // oxygen the fuel already carries.
  auto count = [](const std::map<std::string, int>& c, const char* e) {
    auto it = c.find(e);
    return it == c.end() ? 0 : it->second;
  };
  const int nC = count(fuel_.composition, "C");
  const int nH = count(fuel_.composition, "H");
  const int nO_fuel = count(fuel_.composition, "O");
  const double o_atoms_needed = 2.0 * nC + 0.5 * nH - nO_fuel;
  if (!(o_atoms_needed > 0.0))
    throw ConfigError("propellant " + fuel_.name + " carries enough oxygen to burn itself; "
                      "the stoichiometric mixture ratio is undefined");
  const int nO_ox = count(ox_.composition, "O");
  if (nO_ox <= 0)
    throw ConfigError("oxidiser " + ox_.name + " contains no oxygen");
  // moles of oxidiser per mole of fuel, converted to a mass ratio
  const double mol_ox_per_mol_fuel = o_atoms_needed / static_cast<double>(nO_ox);
  return mol_ox_per_mol_fuel * ox_.molar_mass / fuel_.molar_mass;
}

double PropellantMixture::bulkDensity() const {
  if (!(ox_.density > 0.0) || !(fuel_.density > 0.0))
    throw ConfigError("bulk density needs a density for both propellants");
  const double vo = oxidizerMassFraction() / ox_.density;
  const double vf = fuelMassFraction() / fuel_.density;
  return 1.0 / (vo + vf);
}

PropellantMixture PropellantMixture::withMixtureRatio(double mr) const {
  return PropellantMixture(ox_, fuel_, mr, T_ox_, T_fuel_);
}

}  // namespace ignis
