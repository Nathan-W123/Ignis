// SPDX-License-Identifier: MIT
#include "TestHelpers.hpp"

namespace ignis_test {

const ignis::SpeciesDatabase& fullDatabase() {
  static const ignis::SpeciesDatabase db = ignis::SpeciesDatabase::loadYaml(speciesDbPath());
  return db;
}

const ignis::SpeciesDatabase& rocketDatabase() {
  static const ignis::SpeciesDatabase db =
      fullDatabase().restrictToElements({"C", "H", "O", "N"});
  return db;
}

}  // namespace ignis_test
