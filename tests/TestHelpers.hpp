// SPDX-License-Identifier: MIT
#pragma once
/// Shared fixtures and helpers for the Ignis test suite.
#include <map>
#include <string>
#include <vector>

#include "ignis/thermo/SpeciesDatabase.hpp"

namespace ignis_test {

/// Repository root, injected by CMake.
inline std::string sourceDir() { return IGNIS_SOURCE_DIR; }
inline std::string speciesDbPath() { return sourceDir() + "/data/thermo/ignis_nasa7.yaml"; }
inline std::string configDir() { return sourceDir() + "/configs"; }
inline std::string referenceDir() { return sourceDir() + "/validation/reference"; }

/// Process-wide species database, loaded once.
const ignis::SpeciesDatabase& fullDatabase();
/// The C/H/O subset used for LOX/hydrocarbon work.
const ignis::SpeciesDatabase& rocketDatabase();
/// The H/O subset used for LOX/hydrogen work.
const ignis::SpeciesDatabase& hydrogenDatabase();

/// A CSV table read from validation/reference: `#` lines are comments, the
/// first non-comment line is the header.
class ReferenceTable {
 public:
  explicit ReferenceTable(const std::string& path);
  std::size_t rows() const { return rows_; }
  bool has(const std::string& column) const { return index_.count(column) > 0; }
  /// Numeric cell; NaN for an unparseable entry.
  double num(const std::string& column, std::size_t row) const;
  const std::string& text(const std::string& column, std::size_t row) const;
  const std::vector<std::string>& headers() const { return headers_; }
  const std::string& path() const { return path_; }

 private:
  std::string path_;
  std::vector<std::string> headers_;
  std::map<std::string, std::size_t> index_;
  std::vector<std::vector<std::string>> cells_;  // [column][row]
  std::size_t rows_ = 0;
};

/// Relative difference |a - b| / max(|b|, floor).
double relativeError(double a, double b, double floor_value = 1.0e-30);

}  // namespace ignis_test
