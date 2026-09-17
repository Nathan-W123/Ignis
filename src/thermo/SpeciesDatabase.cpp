// SPDX-License-Identifier: MIT
#include "ignis/thermo/SpeciesDatabase.hpp"

#include <yaml-cpp/yaml.h>

#include <algorithm>
#include <cstdlib>
#include <filesystem>
#include <set>
#include <sstream>

#include "ignis/core/Version.hpp"

namespace ignis {
namespace {

TransportData::Geometry parseGeometry(const std::string& g) {
  if (g == "atom") return TransportData::Geometry::kAtom;
  if (g == "linear") return TransportData::Geometry::kLinear;
  if (g == "nonlinear") return TransportData::Geometry::kNonlinear;
  throw ConfigError("unknown transport geometry '" + g + "'");
}

}  // namespace

SpeciesDatabase SpeciesDatabase::loadYaml(const std::string& path) {
  YAML::Node root;
  try {
    root = YAML::LoadFile(path);
  } catch (const YAML::Exception& e) {
    throw ConfigError("cannot parse species database '" + path + "': " + e.what());
  }
  if (!root["species"] || !root["species"].IsSequence())
    throw ConfigError("species database '" + path + "' has no 'species' sequence");

  SpeciesDatabase db;
  db.provenance_.path = path;
  if (root["generated"]) db.provenance_.generated = root["generated"].as<std::string>();
  if (root["sources"]) {
    for (const auto& s : root["sources"]) {
      std::ostringstream os;
      os << s["name"].as<std::string>("(unnamed)");
      if (s["url"]) os << " <" << s["url"].as<std::string>() << ">";
      if (s["sha256"]) os << " sha256=" << s["sha256"].as<std::string>().substr(0, 16) << "...";
      db.provenance_.sources.push_back(os.str());
    }
  }
  if (root["atomic-weights"]) {
    for (const auto& kv : root["atomic-weights"])
      db.atomic_weight_[kv.first.as<std::string>()] = kv.second.as<double>() * 1.0e-3;  // g/mol -> kg/mol
  }

  for (const auto& node : root["species"]) {
    const auto name = node["name"].as<std::string>();
    std::map<std::string, int> comp;
    for (const auto& kv : node["composition"]) comp[kv.first.as<std::string>()] = kv.second.as<int>();

    const auto tr = node["temperature-ranges"].as<std::vector<double>>();
    if (tr.size() != 3)
      throw ConfigError("species " + name + ": expected 3 temperature-range entries");

    auto readCoeffs = [&](const char* key) {
      const auto v = node[key].as<std::vector<double>>();
      if (v.size() != 7) throw ConfigError("species " + name + ": " + key + " needs 7 values");
      std::array<double, 7> a{};
      std::copy(v.begin(), v.end(), a.begin());
      return a;
    };

    TransportData trans;
    if (node["transport"]) {
      const auto& t = node["transport"];
      trans.geometry = parseGeometry(t["geometry"].as<std::string>());
      trans.well_depth = t["well-depth"].as<double>();
      trans.diameter = t["diameter"].as<double>();
      trans.dipole = t["dipole"].as<double>(0.0);
      trans.polarizability = t["polarizability"].as<double>(0.0);
      trans.rotational_relaxation = t["rotational-relaxation"].as<double>(0.0);
      trans.valid = true;
    }

    db.species_.emplace_back(name, comp,
                             node["molar-mass"].as<double>() * 1.0e-3,  // kg/kmol -> kg/mol
                             tr[0], tr[1], tr[2], readCoeffs("low-coeffs"),
                             readCoeffs("high-coeffs"), trans,
                             node["source"].as<std::string>(""));
  }
  if (db.species_.empty()) throw ConfigError("species database '" + path + "' is empty");
  db.rebuild();
  return db;
}

SpeciesDatabase SpeciesDatabase::fromSpecies(std::vector<Species> species,
                                             std::map<std::string, double> atomic_weights) {
  SpeciesDatabase db;
  db.species_ = std::move(species);
  db.atomic_weight_ = std::move(atomic_weights);
  db.provenance_.path = "<in-memory>";
  db.provenance_.sources.push_back("constructed programmatically");
  if (db.species_.empty()) throw ConfigError("fromSpecies: no species supplied");
  db.rebuild();
  return db;
}

void SpeciesDatabase::rebuild() {
  index_.clear();
  for (std::size_t j = 0; j < species_.size(); ++j) {
    if (index_.count(species_[j].name()))
      throw ConfigError("duplicate species '" + species_[j].name() + "' in database");
    index_[species_[j].name()] = static_cast<int>(j);
  }

  std::set<std::string> els;
  for (const auto& s : species_)
    for (const auto& kv : s.composition()) els.insert(kv.first);
  elements_.assign(els.begin(), els.end());
  element_index_.clear();
  for (std::size_t i = 0; i < elements_.size(); ++i)
    element_index_[elements_[i]] = static_cast<int>(i);

  a_.setZero(static_cast<Eigen::Index>(elements_.size()),
             static_cast<Eigen::Index>(species_.size()));
  mw_.setZero(static_cast<Eigen::Index>(species_.size()));
  for (std::size_t j = 0; j < species_.size(); ++j) {
    mw_(static_cast<Eigen::Index>(j)) = species_[j].molarMass();
    for (const auto& kv : species_[j].composition())
      a_(element_index_.at(kv.first), static_cast<Eigen::Index>(j)) = kv.second;
  }
}

SpeciesDatabase SpeciesDatabase::subset(const std::vector<std::string>& names) const {
  SpeciesDatabase out;
  out.provenance_ = provenance_;
  out.atomic_weight_ = atomic_weight_;
  for (const auto& n : names) {
    const int j = index(n);
    if (j < 0) throw ConfigError("species '" + n + "' is not in database " + provenance_.path);
    out.species_.push_back(species_[static_cast<std::size_t>(j)]);
  }
  if (out.species_.empty()) throw ConfigError("requested species subset is empty");
  out.rebuild();
  return out;
}

SpeciesDatabase SpeciesDatabase::restrictToElements(const std::vector<std::string>& elements) const {
  const std::set<std::string> allowed(elements.begin(), elements.end());
  std::vector<std::string> keep;
  for (const auto& s : species_) {
    bool ok = true;
    for (const auto& kv : s.composition())
      if (!allowed.count(kv.first)) { ok = false; break; }
    if (ok) keep.push_back(s.name());
  }
  if (keep.empty())
    throw ConfigError("no species are composed solely of the requested elements");
  return subset(keep);
}

const Species& SpeciesDatabase::at(const std::string& name) const {
  const int j = index(name);
  if (j < 0) throw ConfigError("species '" + name + "' not present");
  return species_[static_cast<std::size_t>(j)];
}

int SpeciesDatabase::index(const std::string& name) const {
  auto it = index_.find(name);
  return it == index_.end() ? -1 : it->second;
}

int SpeciesDatabase::elementIndex(const std::string& symbol) const {
  auto it = element_index_.find(symbol);
  return it == element_index_.end() ? -1 : it->second;
}

double SpeciesDatabase::atomicWeight(const std::string& symbol) const {
  auto it = atomic_weight_.find(symbol);
  if (it == atomic_weight_.end())
    throw ConfigError("no atomic weight for element '" + symbol + "'");
  return it->second;
}

double SpeciesDatabase::tMaxCommon() const {
  double t = 1e30;
  for (const auto& s : species_) t = std::min(t, s.tMax());
  return t;
}

double SpeciesDatabase::tMinCommon() const {
  double t = 0.0;
  for (const auto& s : species_) t = std::max(t, s.tMin());
  return t;
}

std::vector<std::string> SpeciesDatabase::names() const {
  std::vector<std::string> n;
  n.reserve(species_.size());
  for (const auto& s : species_) n.push_back(s.name());
  return n;
}

std::string findDefaultSpeciesDatabase() {
  namespace fs = std::filesystem;
  std::vector<std::string> tried;
  auto probe = [&](const fs::path& p) -> bool {
    tried.push_back(p.string());
    std::error_code ec;
    return fs::exists(p, ec) && fs::is_regular_file(p, ec);
  };

  const char* env = std::getenv("IGNIS_DATA_DIR");
  if (env != nullptr) {
    const fs::path p = fs::path(env) / "thermo" / "ignis_nasa7.yaml";
    if (probe(p)) return p.string();
    const fs::path p2 = fs::path(env) / "ignis_nasa7.yaml";
    if (probe(p2)) return p2.string();
  }
  {
    const fs::path p = fs::path(kDefaultDataDir) / "thermo" / "ignis_nasa7.yaml";
    if (probe(p)) return p.string();
  }
  for (const char* rel : {"data/thermo/ignis_nasa7.yaml", "../data/thermo/ignis_nasa7.yaml",
                          "../../data/thermo/ignis_nasa7.yaml",
                          "../../../data/thermo/ignis_nasa7.yaml"}) {
    if (probe(rel)) return std::string(rel);
  }
  std::ostringstream os;
  os << "cannot locate the Ignis species database. Tried:";
  for (const auto& t : tried) os << "\n  " << t;
  os << "\nSet IGNIS_DATA_DIR to the directory containing thermo/ignis_nasa7.yaml.";
  throw ConfigError(os.str());
}

}  // namespace ignis
