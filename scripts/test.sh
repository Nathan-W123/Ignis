#!/usr/bin/env bash
# Run the whole test suite, or a tagged subset:  ./scripts/test.sh '[equilibrium]'
set -euo pipefail
cd "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD_DIR="${BUILD_DIR:-build}"
if [[ $# -gt 0 ]]; then
  exec "$BUILD_DIR/bin/ignis_tests" "$@"
fi
cd "$BUILD_DIR" && ctest --output-on-failure -j "${JOBS:-$( (command -v nproc >/dev/null && nproc) || echo 4)}"
