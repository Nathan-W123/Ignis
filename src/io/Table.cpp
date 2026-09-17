// SPDX-License-Identifier: MIT
#include "ignis/io/Table.hpp"

#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <sstream>

#include "ignis/core/Exceptions.hpp"

namespace ignis {

void Table::checkLength(std::size_t n, const std::string& header) const {
  if (!headers_.empty() && n != rows()) {
    std::ostringstream os;
    os << "table '" << name_ << "': column '" << header << "' has " << n
       << " rows but the table has " << rows();
    throw ConfigError(os.str());
  }
}

std::size_t Table::rows() const {
  if (headers_.empty()) return 0;
  for (std::size_t c = 0; c < headers_.size(); ++c) {
    if (is_string_[c] >= 0) return strings_[static_cast<std::size_t>(is_string_[c])].size();
  }
  return numeric_.empty() ? 0 : numeric_.front().size();
}

void Table::addColumn(const std::string& header, std::vector<double> values,
                      const std::string& units) {
  checkLength(values.size(), header);
  headers_.push_back(header);
  units_.push_back(units);
  is_string_.push_back(-1);
  numeric_.push_back(std::move(values));
}

void Table::addColumn(const std::string& header, std::vector<std::string> values) {
  checkLength(values.size(), header);
  headers_.push_back(header);
  units_.push_back("");
  is_string_.push_back(static_cast<int>(strings_.size()));
  strings_.push_back(std::move(values));
}

const std::vector<double>& Table::numeric(std::size_t c) const {
  if (c >= headers_.size()) throw ConfigError("table: column index out of range");
  if (is_string_[c] >= 0) throw ConfigError("table: column '" + headers_[c] + "' is textual");
  std::size_t k = 0;
  for (std::size_t i = 0; i < c; ++i)
    if (is_string_[i] < 0) ++k;
  return numeric_[k];
}

void Table::writeCsv(const std::string& path) const {
  namespace fs = std::filesystem;
  const fs::path p(path);
  if (p.has_parent_path()) {
    std::error_code ec;
    fs::create_directories(p.parent_path(), ec);
  }
  std::ofstream out(path);
  if (!out) throw ConfigError("cannot write CSV file '" + path + "'");

  bool any_units = false;
  for (const auto& u : units_) if (!u.empty()) any_units = true;
  if (any_units) {
    out << "# units:";
    for (std::size_t c = 0; c < headers_.size(); ++c)
      out << ' ' << headers_[c] << '[' << (units_[c].empty() ? "-" : units_[c]) << ']';
    out << '\n';
  }
  for (std::size_t c = 0; c < headers_.size(); ++c) {
    if (c) out << ',';
    out << headers_[c];
  }
  out << '\n';

  const std::size_t n = rows();
  std::size_t num_idx = 0;
  std::vector<std::size_t> numeric_index(headers_.size(), 0);
  for (std::size_t c = 0; c < headers_.size(); ++c)
    if (is_string_[c] < 0) numeric_index[c] = num_idx++;

  char buf[40];
  for (std::size_t r = 0; r < n; ++r) {
    for (std::size_t c = 0; c < headers_.size(); ++c) {
      if (c) out << ',';
      if (is_string_[c] >= 0) {
        out << strings_[static_cast<std::size_t>(is_string_[c])][r];
      } else {
        const double v = numeric_[numeric_index[c]][r];
        if (!std::isfinite(v)) {
          out << (std::isnan(v) ? "nan" : (v > 0 ? "inf" : "-inf"));
        } else {
          std::snprintf(buf, sizeof(buf), "%.10g", v);
          out << buf;
        }
      }
    }
    out << '\n';
  }
}

Json Table::toJson() const {
  Json j = Json::object();
  j["name"] = Json(name_);
  Json cols = Json::object();
  Json units = Json::object();
  std::size_t num_idx = 0;
  for (std::size_t c = 0; c < headers_.size(); ++c) {
    if (is_string_[c] >= 0) {
      cols[headers_[c]] = Json::of(strings_[static_cast<std::size_t>(is_string_[c])]);
    } else {
      cols[headers_[c]] = Json::of(numeric_[num_idx++]);
    }
    units[headers_[c]] = Json(units_[c]);
  }
  j["columns"] = cols;
  j["units"] = units;
  j["rows"] = Json(rows());
  return j;
}

}  // namespace ignis
