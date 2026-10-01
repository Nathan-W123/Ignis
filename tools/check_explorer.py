#!/usr/bin/env python3
"""Check the Explorer's offline path without a compiler or a display.

What the Explorer shows when the binaries are missing is only as good as
three promises, and each has been broken before:

* a preset loaded into the design panel comes back out as itself -- the panel
  once read the wall thickness in the wrong unit and dropped every field it
  had no control for, so no preset but the first could be replayed;
* every preset has a saved solve, and the saved solve is of that preset;
* the shipped engine model is the file its provenance names, is drawn only
  around its own geometry, and sits on the solved wall rather than through it.

Runs in CI with QT_QPA_PLATFORM=offscreen.  Exits non-zero on the first
broken promise, naming it.
"""
from __future__ import annotations

import hashlib
import json
import os
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "python"))
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import numpy as np  # noqa: E402
from PySide6 import QtWidgets  # noqa: E402

from ignis_explorer import engines  # noqa: E402
from ignis_explorer.app import DesignPanel  # noqa: E402
from ignis_explorer.solver import PRESETS, FrozenSolver, _same_design  # noqa: E402

failures = []


def check(ok: bool, what: str) -> None:
    print(("  ok    " if ok else "  FAIL  ") + what)
    if not ok:
        failures.append(what)


def main() -> int:
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication([])  # noqa: F841

    print("presets through the design panel")
    panel = DesignPanel("A", next(iter(PRESETS.values())))
    check(_same_design(panel.design(), next(iter(PRESETS.values()))),
          "the panel's first design is the first preset")
    for name, design in PRESETS.items():
        panel.load(design)
        check(_same_design(panel.design(), design), f"{name!r} round-trips")

    print("saved solves")
    frozen = FrozenSolver()
    for name, design in PRESETS.items():
        result = frozen.run(design)
        check(result.ok and result.replayed, f"{name!r} replays from data/results")
        if result.ok:
            check(bool(result.profile) and bool(result.altitude) and bool(result.composition),
                  f"{name!r} carries the profile, altitude and composition every tab draws")
            if design.cooling:
                check(bool(result.thermal), f"{name!r} carries the thermal profile")

    print("engine models")
    for model in engines.available():
        path = os.path.join(model.directory, model.spec["model"]["file"])
        with open(path, "rb") as fh:
            digest = hashlib.sha256(fh.read()).hexdigest()
        check(digest == model.spec["model"]["sha256"],
              f"{model.key}: the model file is the one its provenance names")
        check("frame" in model.spec, f"{model.key}: has a fitted frame")
        preset = model.preset
        check(preset is not None, f"{model.key}: its preset exists")
        if preset is None:
            continue
        check(engines.for_design(preset) is model, f"{model.key}: drawn around its preset")
        for name, other in PRESETS.items():
            if other is not preset:
                check(engines.for_design(other) is not model,
                      f"{model.key}: not drawn around {name!r}")
        result = frozen.run(preset)
        if not result.ok:
            continue
        x = np.asarray(result.profile["x"], float)
        r = np.asarray(result.profile["radius"], float)
        verts = np.vstack([m.vertices for m in model.meshes(float(x[-1]))])
        radial = np.hypot(verts[:, 1], verts[:, 2])
        wall = np.interp(verts[:, 0], x, r)
        inside = float(np.mean(radial < wall))
        check(inside < 0.02, f"{model.key}: {100 * inside:.1f} % of the model's vertices "
                             "lie inside the solved gas-side wall (want < 2 %)")
        exit_model = model.spec["frame"]["exit_radius"] * model.spec["model"]["units_to_metres"]
        check(abs(exit_model / r[-1] - 1.0) < 0.01,
              f"{model.key}: model exit radius {exit_model:.4f} m is within 1 % of "
              f"the solved {r[-1]:.4f} m")
        check(model.is_published_point(preset), f"{model.key}: its preset is the published point")

    print()
    if failures:
        print(f"{len(failures)} check(s) failed")
        return 1
    print("all checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
