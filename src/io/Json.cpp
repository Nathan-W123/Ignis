// SPDX-License-Identifier: MIT
#include "ignis/io/Json.hpp"

#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>

#include "ignis/core/Exceptions.hpp"

namespace ignis {
namespace {

std::string escape(const std::string& s) {
  std::string out;
  out.reserve(s.size() + 8);
  for (char c : s) {
    switch (c) {
      case '"': out += "\\\""; break;
      case '\\': out += "\\\\"; break;
      case '\n': out += "\\n"; break;
      case '\r': out += "\\r"; break;
      case '\t': out += "\\t"; break;
      default:
        if (static_cast<unsigned char>(c) < 0x20) {
          char buf[8];
          std::snprintf(buf, sizeof(buf), "\\u%04x", c);
          out += buf;
        } else {
          out += c;
        }
    }
  }
  return out;
}

std::string number(double v) {
  if (!std::isfinite(v)) return "null";
  char buf[40];
  std::snprintf(buf, sizeof(buf), "%.17g", v);
  return buf;
}

}  // namespace

Json Json::of(const std::vector<double>& v) {
  Json j = array();
  for (double x : v) j.push(Json(x));
  return j;
}
Json Json::of(const std::vector<int>& v) {
  Json j = array();
  for (int x : v) j.push(Json(x));
  return j;
}
Json Json::of(const std::vector<std::string>& v) {
  Json j = array();
  for (const auto& x : v) j.push(Json(x));
  return j;
}

Json& Json::operator[](const std::string& key) {
  if (type_ == Type::kNull) type_ = Type::kObject;
  if (type_ != Type::kObject) throw IgnisError("JSON: not an object");
  for (auto& kv : obj_)
    if (kv.first == key) return kv.second;
  obj_.emplace_back(key, Json());
  return obj_.back().second;
}

void Json::push(Json v) {
  if (type_ == Type::kNull) type_ = Type::kArray;
  if (type_ != Type::kArray) throw IgnisError("JSON: not an array");
  arr_.push_back(std::move(v));
}

void Json::dumpTo(std::string& out, int indent, int level) const {
  const std::string pad(static_cast<std::size_t>(indent * level), ' ');
  const std::string pad1(static_cast<std::size_t>(indent * (level + 1)), ' ');
  const char* nl = indent > 0 ? "\n" : "";
  switch (type_) {
    case Type::kNull: out += "null"; break;
    case Type::kBool: out += bool_ ? "true" : "false"; break;
    case Type::kNumber: out += number(num_); break;
    case Type::kString: out += '"' + escape(str_) + '"'; break;
    case Type::kArray: {
      if (arr_.empty()) { out += "[]"; break; }
      bool scalars = true;
      for (const auto& e : arr_)
        if (e.type_ == Type::kArray || e.type_ == Type::kObject) { scalars = false; break; }
      if (scalars) {
        // Keep numeric columns on one line so exported tables stay readable.
        out += '[';
        for (std::size_t i = 0; i < arr_.size(); ++i) {
          if (i) out += ", ";
          arr_[i].dumpTo(out, 0, 0);
        }
        out += ']';
        break;
      }
      out += '[';
      out += nl;
      for (std::size_t i = 0; i < arr_.size(); ++i) {
        out += pad1;
        arr_[i].dumpTo(out, indent, level + 1);
        if (i + 1 < arr_.size()) out += ',';
        out += nl;
      }
      out += pad;
      out += ']';
      break;
    }
    case Type::kObject: {
      if (obj_.empty()) { out += "{}"; break; }
      out += '{';
      out += nl;
      for (std::size_t i = 0; i < obj_.size(); ++i) {
        out += pad1;
        out += '"' + escape(obj_[i].first) + "\": ";
        obj_[i].second.dumpTo(out, indent, level + 1);
        if (i + 1 < obj_.size()) out += ',';
        out += nl;
      }
      out += pad;
      out += '}';
      break;
    }
  }
}

std::string Json::dump(int indent) const {
  std::string out;
  dumpTo(out, indent, 0);
  out += '\n';
  return out;
}

void Json::writeFile(const std::string& path, int indent) const {
  namespace fs = std::filesystem;
  const fs::path p(path);
  if (p.has_parent_path()) {
    std::error_code ec;
    fs::create_directories(p.parent_path(), ec);
  }
  std::ofstream out(path);
  if (!out) throw ConfigError("cannot write JSON file '" + path + "'");
  out << dump(indent);
}

}  // namespace ignis
