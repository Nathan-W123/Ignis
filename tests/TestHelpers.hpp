// SPDX-License-Identifier: MIT
#pragma once
/// Shared fixtures for the Ignis test suite.
#include <memory>
#include <string>
#include <vector>

#include "ignis/thermo/SpeciesDatabase.hpp"

namespace ignis_test {

/// Path to the repository root, injected by CMake.
inline std::string sourceDir() { return IGNIS_SOURCE_DIR; }
inline std::string speciesDbPath() { return sourceDir() + "/data/thermo/ignis_nasa7.yaml"; }
inline std::string configDir() { return sourceDir() + "/configs"; }

/// Process-wide species database, loaded once.
const ignis::SpeciesDatabase& fullDatabase();

/// The CHON species subset used for LOX/hydrocarbon work.
const ignis::SpeciesDatabase& rocketDatabase();

}  // namespace ignis_test
