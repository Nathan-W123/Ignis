// SPDX-License-Identifier: MIT
#pragma once
/// \file SpeciesDatabase.hpp
/// \brief Container for the species set used by a calculation.

#include <map>
#include <string>
#include <vector>

#include <Eigen/Dense>

#include "ignis/thermo/Species.hpp"

namespace ignis {

/// Provenance record for the on-disk database.
struct DatabaseProvenance {
  std::string generated;
  std::vector<std::string> sources;   ///< human-readable citation lines
  std::string path;
};

/// An indexed collection of species plus the element/species stoichiometric
/// matrix used by the equilibrium solver.
///
/// Row i, column j of the element matrix `a` holds the number of atoms of
/// element i in one molecule of species j.
class SpeciesDatabase {
 public:
  /// Load every species in the YAML file at `path`.
  static SpeciesDatabase loadYaml(const std::string& path);

  /// Build a database directly from species objects.  Used by the verification
  /// suite to construct synthetic gases (for example a calorically perfect one)
  /// without going through a file.
  static SpeciesDatabase fromSpecies(std::vector<Species> species,
                                     std::map<std::string, double> atomic_weights);

  /// Return the subset containing exactly `names`, preserving that order.
  /// Throws ConfigError if a name is unknown.
  SpeciesDatabase subset(const std::vector<std::string>& names) const;

  /// Return the subset whose species are composed only of `elements`.
  SpeciesDatabase restrictToElements(const std::vector<std::string>& elements) const;

  std::size_t size() const { return species_.size(); }
  const std::vector<Species>& species() const { return species_; }
  const Species& operator[](std::size_t j) const { return species_[j]; }
  const Species& at(const std::string& name) const;
  /// Index of `name`, or -1 if absent.
  int index(const std::string& name) const;
  bool has(const std::string& name) const { return index(name) >= 0; }

  /// Element symbols, sorted, as they index the rows of elementMatrix().
  const std::vector<std::string>& elements() const { return elements_; }
  std::size_t numElements() const { return elements_.size(); }
  int elementIndex(const std::string& symbol) const;

  /// (numElements x size()) matrix of atom counts.
  const Eigen::MatrixXd& elementMatrix() const { return a_; }
  /// Molar masses, kg/mol, one per species.
  const Eigen::VectorXd& molarMasses() const { return mw_; }
  /// Atomic weight of an element, kg/mol.
  double atomicWeight(const std::string& symbol) const;

  const DatabaseProvenance& provenance() const { return provenance_; }

  /// Hottest temperature valid for *every* species in the set, K.
  double tMaxCommon() const;
  /// Coldest temperature valid for *every* species in the set, K.
  double tMinCommon() const;

  /// Names, in index order.
  std::vector<std::string> names() const;

 private:
  void rebuild();

  std::vector<Species> species_;
  std::map<std::string, int> index_;
  std::vector<std::string> elements_;
  std::map<std::string, int> element_index_;
  std::map<std::string, double> atomic_weight_;  // kg/mol
  Eigen::MatrixXd a_;
  Eigen::VectorXd mw_;
  DatabaseProvenance provenance_;
};

/// Locate the shipped species database.
///
/// Search order: `IGNIS_DATA_DIR` environment variable, the compiled-in install
/// path, then a handful of paths relative to the current directory so that the
/// tools work from a build tree.  Throws ConfigError with the full list of
/// attempted locations if nothing is found.
std::string findDefaultSpeciesDatabase();

}  // namespace ignis
