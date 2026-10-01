#!/usr/bin/env python3
"""Find the axis of an engine model, so the Explorer can stand it on a solve.

A display model arrives in whatever frame its author used: NASA's Shuttle
engine bell is in inches, in orbiter coordinates, with the engine canted.  To
draw it around an Ignis solution the Explorer needs three things, and this
tool measures them from the mesh rather than leaving them to be eyeballed:

* the exit plane's centre,
* the direction of the axis (pointing upstream, from the exit to the throat),
* the exit radius -- the model's own, to be compared with the solved one.

Method: take the named reference surface (a body of revolution in the model),
pick the vertices within a sliver of each end along a rough axis, fit a circle
to each ring in 3-D (plane by SVD, then a linear least-squares circle in that
plane), and run the axis through the two centres.  The ring fits report their
worst residual, so a model that is not round where it should be says so.

Writes the result into the engine's engine.json under "frame".

    python3 tools/fit_engine_model.py data/engines/rs25
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "python"))
from ignis_explorer.meshio import read_glb  # noqa: E402


def _circle(points: np.ndarray):
    centre = points.mean(axis=0)
    _, _, vt = np.linalg.svd(points - centre)
    e1, e2, normal = vt[0], vt[1], vt[2]
    q = np.c_[(points - centre) @ e1, (points - centre) @ e2]
    a = np.c_[2.0 * q, np.ones(len(q))]
    sol, *_ = np.linalg.lstsq(a, (q ** 2).sum(axis=1), rcond=None)
    radius = float(np.sqrt(sol[2] + sol[0] ** 2 + sol[1] ** 2))
    residual = np.linalg.norm(q - sol[:2], axis=1) - radius
    return centre + sol[0] * e1 + sol[1] * e2, normal, radius, residual


def fit_ring(points: np.ndarray):
    """Centre, normal, radius and worst residual of a circle through 3-D points.

    An end ring often carries both edges of a lip -- the gas side and the
    outside of the wall -- and one circle through both lands between them,
    matching neither.  So when the first fit leaves points more than 1 % of
    the radius off it, the points are split at the first circle into the two
    edges, each is fitted on its own, and the inner (gas-side) edge is the
    one returned.  How many points were set aside is returned too, so a ring
    that needed this is visible rather than silently cleaned.
    """
    centre, normal, radius, residual = _circle(points)
    off = np.abs(residual) > 0.01 * radius
    if not off.any():
        return centre, normal, radius, float(np.abs(residual).max()), 0
    inner = residual < 0.0
    if inner.sum() < 6:
        raise ValueError("the ring is not round and has too few points on its "
                         "inner edge to fit that edge alone")
    centre, normal, radius, residual = _circle(points[inner])
    return centre, normal, radius, float(np.abs(residual).max()), int((~inner).sum())


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("engine_dir", help="directory holding engine.json")
    ap.add_argument("--sliver", type=float, default=0.02,
                    help="fraction of the length taken as each end ring (default 0.02)")
    args = ap.parse_args()

    spec_path = os.path.join(args.engine_dir, "engine.json")
    with open(spec_path) as fh:
        spec = json.load(fh)
    model = spec["model"]
    prims = read_glb(os.path.join(args.engine_dir, model["file"]))
    ref = [p for p in prims if p.material == model["reference_surface"]]
    if not ref:
        print(f"no primitive with material {model['reference_surface']!r}; "
              f"have {[p.material for p in prims]}", file=sys.stderr)
        return 1
    v = np.vstack([np.unique(p.vertices.round(6), axis=0) for p in ref])

    rough = np.zeros(3)
    rough["xyz".index(model["rough_axis"])] = 1.0
    s = v @ rough
    lo, hi = s.min(), s.max()
    sliver = args.sliver * (hi - lo)
    rings = [fit_ring(v[s <= lo + sliver]), fit_ring(v[s >= hi - sliver])]
    # The exit is the wider end.
    exit_ring, top_ring = sorted(rings, key=lambda r: -r[2])
    upstream = top_ring[0] - exit_ring[0]
    length = float(np.linalg.norm(upstream))
    upstream /= length
    tilt = float(np.degrees(np.arccos(min(1.0, abs(upstream @ exit_ring[1])))))

    frame = {
        "method": "tools/fit_engine_model.py: circles fitted to the end rings "
                  f"of the {model['reference_surface']!r} surface",
        "exit_centre": [round(float(c), 6) for c in exit_ring[0]],
        "upstream": [round(float(c), 8) for c in upstream],
        "exit_radius": round(exit_ring[2], 4),
        "top_radius": round(top_ring[2], 4),
        "length": round(length, 4),
        "exit_ring_worst_residual": round(exit_ring[3], 4),
        "exit_ring_points_dropped": exit_ring[4],
        "top_ring_worst_residual": round(top_ring[3], 4),
        "top_ring_points_dropped": top_ring[4],
        "exit_plane_tilt_deg": round(tilt, 3),
    }
    spec["frame"] = frame
    with open(spec_path, "w") as fh:
        json.dump(spec, fh, indent=2)
        fh.write("\n")
    scale = model["units_to_metres"]
    print(f"exit radius {frame['exit_radius']} model units = {frame['exit_radius'] * scale:.4f} m")
    print(f"top radius  {frame['top_radius']} model units = {frame['top_radius'] * scale:.4f} m")
    print(f"length      {frame['length']} model units = {frame['length'] * scale:.4f} m")
    print(f"exit plane tilted {tilt:.2f} deg to the axis; ring residuals "
          f"{frame['exit_ring_worst_residual']} / {frame['top_ring_worst_residual']}")
    print(f"wrote the frame into {spec_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
