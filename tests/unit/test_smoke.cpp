// SPDX-License-Identifier: MIT
#include <catch2/catch_test_macros.hpp>
#include "TestHelpers.hpp"
TEST_CASE("species database loads", "[thermo]") {
  const auto& db = ignis_test::fullDatabase();
  REQUIRE(db.size() == 40);
}
