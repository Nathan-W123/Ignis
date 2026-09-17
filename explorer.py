#!/usr/bin/env python3
"""Launcher for the Ignis Engine Explorer.

    python3 explorer.py

The Explorer needs the Ignis binaries.  Build them first with

    ./scripts/build.sh

or point it at an existing build directory:

    IGNIS_BUILD_DIR=/path/to/build python3 explorer.py
"""
from __future__ import annotations

import os
import sys

_ROOT = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(_ROOT, "python"))

# Keep matplotlib's font cache inside the repo rather than in the user's home.
os.environ.setdefault("MPLCONFIGDIR", os.path.join(_ROOT, ".cache", "matplotlib"))
os.makedirs(os.environ["MPLCONFIGDIR"], exist_ok=True)


def main() -> int:
    try:
        from ignis_explorer.app import launch
    except ImportError as exc:  # pragma: no cover - dependency guidance
        missing = getattr(exc, "name", "") or str(exc)
        print(f"Ignis Engine Explorer needs a dependency that is not installed: {missing}\n\n"
              "Install the Explorer's requirements:\n"
              "    pip install -r python/requirements-explorer.txt\n", file=sys.stderr)
        return 1
    return launch(_ROOT)


if __name__ == "__main__":
    raise SystemExit(main())
