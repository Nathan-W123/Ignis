"""Solved engines frozen to disk, so the plume runs without a compiler.

The Explorer drives the compiled solver, which means a C++ toolchain stands
between someone and the first thing on screen.  The plume march does not need
one: `ignis_viz.plume` is NumPy, and all it wants from the nozzle is eight
numbers at the exit plane.

`tools/make_case_library.py` writes those numbers out, with the contour for
drawing and the altitude-to-ambient table, and this module reads them back.

WHAT IS REPLAYED AND WHAT IS LIVE
---------------------------------
Replayed: chamber equilibrium, the nozzle expansion, the exit state, the
performance figures.  Those were solved by the binaries at the commit each
file records and are not recomputed here.

Live: the plume.  Ambient pressure is a free parameter, and the axisymmetric
Euler march runs from scratch every time it is asked for.  So the part that
moves on screen is a solution computed on the machine showing it, not a
recording -- and the UI says which is which rather than presenting the whole
thing as one live calculation.
"""
from __future__ import annotations

import json
import os
from dataclasses import dataclass
from typing import List, Optional, Tuple

import numpy as np

CASE_DIR = os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
    "data", "cases")


@dataclass
class Case:
    """One frozen engine: enough to draw it and to march its plume."""

    stem: str
    name: str
    description: str
    note: str
    solved_by: str
    exit: dict
    contour: np.ndarray                 # (n, 2) of (x, r) in metres
    performance: dict
    altitudes: np.ndarray               # m
    ambient: np.ndarray                 # Pa
    separated: np.ndarray               # 1 where Ignis predicts separation

    @classmethod
    def load(cls, path: str) -> "Case":
        with open(path) as fh:
            doc = json.load(fh)
        if doc.get("schema") != 1:
            raise ValueError(f"{os.path.basename(path)}: unknown case schema "
                             f"{doc.get('schema')!r}")
        table = doc["altitude_table"]
        return cls(
            stem=os.path.splitext(os.path.basename(path))[0],
            name=doc["name"], description=doc.get("description", ""),
            note=doc.get("note", ""), solved_by=doc.get("solved_by", "?"),
            exit=doc["exit"],
            contour=np.asarray(doc["contour"], dtype=float),
            performance=doc.get("performance", {}),
            altitudes=np.asarray(table["altitude"], dtype=float),
            ambient=np.asarray(table["ambient_pressure"], dtype=float),
            separated=np.asarray(table["separation_predicted"], dtype=int))

    def exit_state(self):
        """The exit plane, as `ignis_viz.plume` wants it."""
        from ignis_viz.plume import ExitState
        e = self.exit
        return ExitState(radius=e["radius"], pressure=e["pressure"],
                         temperature=e["temperature"], density=e["density"],
                         velocity=e["velocity"], mach=e["mach"],
                         gamma=e["gamma_s"], molar_mass=e["molar_mass"])

    def at_altitude(self, altitude_m: float) -> Tuple[float, bool]:
        """Ambient pressure there, and whether Ignis predicts separation.

        The pressure is interpolated because the atmosphere is smooth.  The
        separation flag is NOT: it is a threshold the solver crossed at a
        particular altitude, so it is read from the nearest tabulated row
        rather than blended into a fraction of a flag.
        """
        p = float(np.interp(altitude_m, self.altitudes, self.ambient))
        i = int(np.argmin(np.abs(self.altitudes - altitude_m)))
        return p, bool(self.separated[i])

    def lowest_attached_altitude(self) -> float:
        """The first altitude in the table where the nozzle flows full."""
        attached = np.flatnonzero(self.separated == 0)
        if attached.size == 0:
            return float(self.altitudes[-1])
        return float(self.altitudes[attached[0]])

    def drawn_contour(self, exit_radius: float) -> Optional[np.ndarray]:
        """The last stretch of wall, shifted so the exit plane is at x = 0."""
        if self.contour.shape[0] < 2:
            return None
        x = self.contour[:, 0] - self.contour[-1, 0]
        keep = x >= -6.0 * exit_radius
        if keep.sum() < 2:
            return None
        return np.column_stack([x[keep], self.contour[keep, 1]])


def available(directory: str = CASE_DIR) -> List[Case]:
    """Every case on disk, by name.  Unreadable ones are skipped, not raised.

    A case library that has one bad file should still offer the others: this
    is a fallback path, and failing it entirely would leave someone with a
    window that does nothing and no way to tell why.
    """
    if not os.path.isdir(directory):
        return []
    out = []
    for name in sorted(os.listdir(directory)):
        if not name.endswith(".json"):
            continue
        try:
            out.append(Case.load(os.path.join(directory, name)))
        except (OSError, KeyError, ValueError, json.JSONDecodeError):
            continue
    return sorted(out, key=lambda c: c.name)


__all__ = ["CASE_DIR", "Case", "available"]
