// SPDX-License-Identifier: MIT
#include "TestHelpers.hpp"

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace ignis_test {

const ignis::SpeciesDatabase& fullDatabase() {
  static const ignis::SpeciesDatabase db = ignis::SpeciesDatabase::loadYaml(speciesDbPath());
  return db;
}

const ignis::SpeciesDatabase& rocketDatabase() {
  static const ignis::SpeciesDatabase db = fullDatabase().restrictToElements({"C", "H", "O"});
  return db;
}

const ignis::SpeciesDatabase& hydrogenDatabase() {
  static const ignis::SpeciesDatabase db = fullDatabase().restrictToElements({"H", "O"});
  return db;
}

ReferenceTable::ReferenceTable(const std::string& path) : path_(path) {
  std::ifstream in(path);
  if (!in) throw std::runtime_error("cannot open reference table '" + path + "'");
  std::string line;
  bool have_header = false;
  std::vector<std::vector<std::string>> rows;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#') continue;
    std::vector<std::string> fields;
    std::string cell;
    std::istringstream ss(line);
    while (std::getline(ss, cell, ',')) fields.push_back(cell);
    if (!have_header) {
      headers_ = fields;
      for (std::size_t i = 0; i < headers_.size(); ++i) index_[headers_[i]] = i;
      have_header = true;
    } else {
      rows.push_back(std::move(fields));
    }
  }
  rows_ = rows.size();
  cells_.assign(headers_.size(), std::vector<std::string>(rows_));
  for (std::size_t r = 0; r < rows_; ++r)
    for (std::size_t c = 0; c < headers_.size() && c < rows[r].size(); ++c)
      cells_[c][r] = rows[r][c];
}

double ReferenceTable::num(const std::string& column, std::size_t row) const {
  const auto it = index_.find(column);
  if (it == index_.end()) throw std::runtime_error("reference table has no column " + column);
  const std::string& s = cells_[it->second].at(row);
  if (s.empty()) return std::nan("");
  return std::atof(s.c_str());
}

const std::string& ReferenceTable::text(const std::string& column, std::size_t row) const {
  const auto it = index_.find(column);
  if (it == index_.end()) throw std::runtime_error("reference table has no column " + column);
  return cells_[it->second].at(row);
}

double relativeError(double a, double b, double floor_value) {
  return std::abs(a - b) / std::max(std::abs(b), floor_value);
}

}  // namespace ignis_test
