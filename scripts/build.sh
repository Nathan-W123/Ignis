#!/usr/bin/env bash
# Configure and build Ignis in Release mode.
set -euo pipefail
cd "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BUILD_DIR="${BUILD_DIR:-build}"
cmake -S . -B "$BUILD_DIR" -DCMAKE_BUILD_TYPE="${BUILD_TYPE:-Release}" "$@"
cmake --build "$BUILD_DIR" -j "${JOBS:-$( (command -v nproc >/dev/null && nproc) || echo 4)}"
echo "binaries in $BUILD_DIR/bin"
