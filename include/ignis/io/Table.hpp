// SPDX-License-Identifier: MIT
#pragma once
/// \file Table.hpp
/// \brief Column-oriented result table with CSV and JSON export.

#include <string>
#include <vector>

#include "ignis/io/Json.hpp"

namespace ignis {

/// A rectangular table of named numeric (or string) columns.
class Table {
 public:
  explicit Table(std::string name = "table") : name_(std::move(name)) {}

  /// Append a numeric column.  Every column must have the same length.
  void addColumn(const std::string& header, std::vector<double> values,
                 const std::string& units = "");
  /// Append a string column.
  void addColumn(const std::string& header, std::vector<std::string> values);

  std::size_t rows() const;
  std::size_t columns() const { return headers_.size(); }
  const std::vector<std::string>& headers() const { return headers_; }
  const std::vector<std::string>& units() const { return units_; }
  /// Numeric column by index; throws if that column is a string column.
  const std::vector<double>& numeric(std::size_t c) const;

  /// Write CSV, creating parent directories.  The first line is a comment
  /// carrying the units, the second is the header row.
  void writeCsv(const std::string& path) const;
  /// Represent as a JSON object of arrays.
  Json toJson() const;

  const std::string& name() const { return name_; }

 private:
  void checkLength(std::size_t n, const std::string& header) const;

  std::string name_;
  std::vector<std::string> headers_, units_;
  std::vector<std::vector<double>> numeric_;
  std::vector<std::vector<std::string>> strings_;
  std::vector<int> is_string_;  // index into strings_ or -1
};

}  // namespace ignis
