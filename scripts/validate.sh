#!/usr/bin/env bash
# Run only the verification and validation cases and print their measured errors.
set -euo pipefail
cd "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD_DIR="${BUILD_DIR:-build}"
exec "$BUILD_DIR/bin/ignis_tests" "[validation],[verification]" --success --reporter console 2>&1 |
  grep -vE '^\s*$|with expansion|^-+$|^\.\.\.|PASSED|REQUIRE'
