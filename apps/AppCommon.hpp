// SPDX-License-Identifier: MIT
#pragma once
/// Helpers shared by the Ignis command-line applications.
#include <filesystem>
#include <iostream>
#include <string>

#include "ignis/core/Version.hpp"
#include "ignis/io/Cli.hpp"
#include "ignis/io/Config.hpp"
#include "ignis/io/Json.hpp"
#include "ignis/optimize/Parameters.hpp"

namespace ignis {
namespace app {

/// Apply the command-line overrides and output settings to a configuration.
inline void applyCommandLine(EngineConfig& cfg, const CommandLine& cli) {
  for (const auto& kv : cli.overrides) applyParameter(cfg, kv.first, kv.second);
  if (!cli.output_dir.empty()) cfg.output_directory = cli.output_dir;
  if (!cli.prefix.empty()) cfg.output_prefix = cli.prefix;
}

/// Build "<dir>/<prefix>_<name>", creating the directory.
inline std::string outputPath(const EngineConfig& cfg, const std::string& name) {
  namespace fs = std::filesystem;
  std::error_code ec;
  fs::create_directories(cfg.output_directory, ec);
  return (fs::path(cfg.output_directory) / (cfg.output_prefix + "_" + name)).string();
}

/// Provenance block attached to every JSON output.
inline Json provenance(const EngineConfig& cfg) {
  Json j = Json::object();
  j["ignis_version"] = Json(std::string(kVersion));
  j["git"] = Json(std::string(kGitDescribe));
  j["build_type"] = Json(std::string(kBuildType));
  j["compiler"] = Json(std::string(kCompilerId) + " " + kCompilerVersion);
  j["config_file"] = Json(cfg.file);
  j["case_name"] = Json(cfg.name);
  if (!cfg.description.empty()) j["description"] = Json(cfg.description);
  return j;
}

/// Standard main() wrapper: turns any Ignis or std exception into a clear
/// message on stderr and a non-zero exit status.
template <typename Fn>
int run(Fn&& fn) {
  try {
    return fn();
  } catch (const IgnisError& e) {
    std::cerr << "ignis: " << e.what() << "\n";
    return 2;
  } catch (const std::exception& e) {
    std::cerr << "ignis: unexpected error: " << e.what() << "\n";
    return 3;
  }
}

}  // namespace app
}  // namespace ignis
