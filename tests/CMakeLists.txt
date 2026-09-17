# SPDX-License-Identifier: MIT
file(GLOB_RECURSE IGNIS_TEST_SOURCES CONFIGURE_DEPENDS
     ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)
add_executable(ignis_tests ${IGNIS_TEST_SOURCES})
target_link_libraries(ignis_tests PRIVATE ignis Catch2::Catch2WithMain ${IGNIS_YAML_TARGET})
target_include_directories(ignis_tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_definitions(ignis_tests PRIVATE
  IGNIS_SOURCE_DIR="${PROJECT_SOURCE_DIR}")
set_target_properties(ignis_tests PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

include(Catch)
catch_discover_tests(ignis_tests
  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
  PROPERTIES ENVIRONMENT "IGNIS_DATA_DIR=${IGNIS_DATA_DIR}")
