// SPDX-License-Identifier: MIT
#pragma once
/// \file Json.hpp
/// \brief Minimal JSON value type and writer.
///
/// Ignis only ever *writes* JSON, so a small self-contained emitter avoids a
/// dependency.  Doubles are written with 17 significant digits (round-trip
/// exact) and non-finite values are written as null with a sibling note, since
/// JSON has no representation for them.

#include <map>
#include <memory>
#include <string>
#include <vector>

namespace ignis {

class Json {
 public:
  enum class Type { kNull, kBool, kNumber, kString, kArray, kObject };

  Json() = default;
  Json(bool v) : type_(Type::kBool), bool_(v) {}
  Json(double v) : type_(Type::kNumber), num_(v) {}
  Json(int v) : type_(Type::kNumber), num_(v) {}
  Json(std::size_t v) : type_(Type::kNumber), num_(static_cast<double>(v)) {}
  Json(const char* v) : type_(Type::kString), str_(v) {}
  Json(std::string v) : type_(Type::kString), str_(std::move(v)) {}

  static Json array() { Json j; j.type_ = Type::kArray; return j; }
  static Json object() { Json j; j.type_ = Type::kObject; return j; }
  static Json of(const std::vector<double>& v);
  static Json of(const std::vector<int>& v);
  static Json of(const std::vector<std::string>& v);

  Json& operator[](const std::string& key);
  void push(Json v);

  Type type() const { return type_; }
  /// Serialise with two-space indentation.
  std::string dump(int indent = 2) const;
  /// Write to a file, creating parent directories as needed.
  void writeFile(const std::string& path, int indent = 2) const;

 private:
  void dumpTo(std::string& out, int indent, int level) const;

  Type type_ = Type::kNull;
  bool bool_ = false;
  double num_ = 0.0;
  std::string str_;
  std::vector<Json> arr_;
  std::vector<std::pair<std::string, Json>> obj_;  // insertion-ordered
};

}  // namespace ignis
