"""Real engine hardware, drawn around a solve of that engine.

An engine model here is a published 3-D model plus the facts needed to stand
it on an Ignis solution: its units, its axis (measured from the mesh by
`tools/fit_engine_model.py`, not eyeballed), and the design it belongs to.
Each lives in `data/engines/<name>/` with an `engine.json` and a README that
says where every number came from.

The model is display only.  The physics of a design comes from the design --
published dimensions and operating point, solved by Ignis -- and the model is
drawn only when the solved geometry is that engine's geometry.  Change the
throat or the expansion ratio and the hardware no longer describes the
nozzle being solved, so it is not drawn rather than being drawn wrongly.
"""
from __future__ import annotations

import json
import os
from dataclasses import dataclass, replace
from typing import Dict, List, Optional

import numpy as np

from .meshio import read_glb
from .solver import PRESETS, Design, _same_design
from .viewport3d import Mesh

ENGINES_DIR = os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
    "data", "engines")

# The fields that make two designs the same piece of hardware.  Propellant,
# pressure and altitude are how it is run, not what it is.
_GEOMETRY = ("throat_radius", "expansion_ratio", "contraction_ratio",
             "chamber_length", "bell_length_fraction", "contour_file")


def srgb(linear: float) -> float:
    """glTF base colours are linear; the viewport's vertex colours are not."""
    if linear <= 0.0031308:
        return 12.92 * linear
    return 1.055 * linear ** (1.0 / 2.4) - 0.055


@dataclass
class EngineModel:
    key: str
    directory: str
    spec: dict

    @property
    def title(self) -> str:
        return self.spec.get("title", self.key)

    @property
    def credit(self) -> str:
        return self.spec["model"].get("credit", "")

    @property
    def preset(self) -> Optional[Design]:
        return PRESETS.get(self.spec.get("preset", ""))

    @property
    def published(self) -> Dict[str, dict]:
        return self.spec.get("published", {})

    def matches(self, design: Optional[Design]) -> bool:
        """True when `design` is this engine's hardware, however it is run."""
        ref = self.preset
        if design is None or ref is None:
            return False
        for key in _GEOMETRY:
            a, b = getattr(design, key), getattr(ref, key)
            if isinstance(a, float):
                if abs(a - b) > 1e-9 * max(1.0, abs(b)):
                    return False
            elif a != b:
                return False
        return True

    def is_published_point(self, design: Optional[Design]) -> bool:
        """True when `design` is the operating point the published data are for.

        Altitude is left out on purpose: the published numbers this compares
        against are vacuum performance and the chamber's own cooling, neither
        of which depends on where the engine is.
        """
        ref = self.preset
        if design is None or ref is None:
            return False
        return _same_design(replace(design, altitude=ref.altitude), ref)

    def meshes(self, x_exit: float) -> List[Mesh]:
        """The model in the solver's frame, its exit plane on `x_exit`.

        The solver's x runs downstream from the injector face and the engine
        is axisymmetric about it, so the model is moved, scaled to metres and
        rotated until its fitted axis is the x axis pointing upstream from
        the exit.  Faces are cut away in the shader, not here, so the cut edge
        is straight however coarse the model is.
        """
        m, frame = self.spec["model"], self.spec["frame"]
        scale = float(m["units_to_metres"])
        origin = np.array(frame["exit_centre"], float)
        up = np.array(frame["upstream"], float)
        up /= np.linalg.norm(up)
        # Any two directions square to the axis complete the frame; roll_deg
        # turns the model about its axis so its plumbing faces where it reads.
        helper = np.array([0.0, 0.0, 1.0]) if abs(up[2]) < 0.9 else np.array([0.0, 1.0, 0.0])
        e1 = np.cross(helper, up)
        e1 /= np.linalg.norm(e1)
        e2 = np.cross(up, e1)
        roll = np.radians(float(m.get("roll_deg", 0.0)))
        e1, e2 = (np.cos(roll) * e1 + np.sin(roll) * e2,
                  -np.sin(roll) * e1 + np.cos(roll) * e2)
        # Rows map model directions to solver x, y, z.  Upstream is -x; keep
        # the map a rotation, never a reflection, so the model is not drawn
        # as its own mirror image.
        rot = np.stack([-up, e1, e2])
        if np.linalg.det(rot) < 0.0:
            rot[2] = -rot[2]

        out: List[Mesh] = []
        for prim in read_glb(os.path.join(self.directory, m["file"])):
            if prim.material not in m.get("draw", []):
                continue
            verts = (prim.vertices - origin) @ rot.T * scale
            verts[:, 0] += x_exit
            normals = None if prim.normals is None else prim.normals @ rot.T
            colour = np.array([srgb(c) for c in prim.colour])
            out.append(Mesh(vertices=verts, faces=prim.faces,
                            colors=np.tile(colour, (len(verts), 1)),
                            normals=normals, cut_in_shader=True,
                            name=f"{self.spec.get('name', self.key)} ({prim.material})"))
        return out

    def skin_colour(self):
        """The model's own colour, for the parts of the engine Ignis draws."""
        m = self.spec["model"]
        for prim in read_glb(os.path.join(self.directory, m["file"])):
            if prim.material in m.get("draw", []):
                return tuple(srgb(c) for c in prim.colour)
        return None


_CACHE: Optional[List[EngineModel]] = None


def available(directory: str = ENGINES_DIR) -> List[EngineModel]:
    """Every engine model shipped in `data/engines`, read once."""
    global _CACHE
    if _CACHE is None:
        found = []
        if os.path.isdir(directory):
            for key in sorted(os.listdir(directory)):
                path = os.path.join(directory, key, "engine.json")
                if os.path.exists(path):
                    with open(path) as fh:
                        found.append(EngineModel(key, os.path.dirname(path), json.load(fh)))
        _CACHE = found
    return _CACHE


def for_design(design: Optional[Design]) -> Optional[EngineModel]:
    """The engine model whose hardware `design` is, if any."""
    for model in available():
        if "frame" in model.spec and model.matches(design):
            return model
    return None


__all__ = ["EngineModel", "available", "for_design", "ENGINES_DIR"]
