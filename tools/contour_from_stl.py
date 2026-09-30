#!/usr/bin/env python3
"""Extract an axisymmetric wall contour from a CAD model of an engine.

WHAT THIS DOES, AND WHAT IT CANNOT
----------------------------------
Ignis solves a quasi-1D internal flow, so the only thing it can take from a
three-dimensional model is the wall radius as a function of axial position.
That is not a limitation of this script -- it is what the solver is.  A real
engine's manifolds, injector face, channel routing and mounting hardware are
not represented and cannot be: there is nowhere for them to go.

So this reads a triangle mesh and reduces it to (x, r) pairs, which is exactly
the revolve profile the chamber and nozzle were drawn from in the first place.
Feed the result to Ignis with

    nozzle:
      contour_file: my_engine_contour.csv

STL only.  STEP and IGES are boundary representations that need a real CAD
kernel (OpenCASCADE) to evaluate; exporting STL from whatever drew the model is
a better use of your time than making this script depend on one.

THE INNER WALL, NOT THE OUTER ONE
---------------------------------
A solid engine model has two surfaces at every station: the gas-side wall and
the outside of the structure.  The flow sees the inner one.  By default this
takes the *minimum* radius in each axial bin, which is the gas side for a
model of the solid; `--outer` takes the maximum instead, which is what a model
of the flow volume (a "negative" or fluid domain) needs.  Getting this backwards
produces a contour that looks plausible and is wrong, so the script reports
which convention it used and what the other one would have given.

THE AXIS
--------
The model has to be axisymmetric about a coordinate axis and this script does
not hunt for one: pass `--axis x|y|z`.  It checks how axisymmetric the mesh
actually is about that axis and reports the spread, because a model that is not
a body of revolution -- one that still has its bolt flanges, say -- will
produce a contour with no error message and no meaning.
"""
from __future__ import annotations

import argparse
import struct
import sys
from typing import Tuple

import numpy as np

AXES = {"x": 0, "y": 1, "z": 2}


def read_stl(path: str) -> np.ndarray:
    """Vertices of every triangle, (n, 3, 3), from a binary or ASCII STL."""
    with open(path, "rb") as fh:
        head = fh.read(5)
        fh.seek(0)
        if head[:5].lower() == b"solid":
            text = fh.read().decode("utf-8", errors="replace")
            if "facet" in text:
                return _read_ascii(text)
            fh.seek(0)
        return _read_binary(fh)


def _read_ascii(text: str) -> np.ndarray:
    verts = []
    for line in text.splitlines():
        parts = line.split()
        if len(parts) == 4 and parts[0] == "vertex":
            verts.append([float(v) for v in parts[1:]])
    if len(verts) < 3 or len(verts) % 3:
        raise ValueError(
            f"the ASCII STL yielded {len(verts)} vertices, which is not a whole "
            f"number of triangles")
    return np.asarray(verts, dtype=float).reshape(-1, 3, 3)


def _read_binary(fh) -> np.ndarray:
    fh.seek(80)
    raw = fh.read(4)
    if len(raw) != 4:
        raise ValueError("the file is too short to be a binary STL")
    count = struct.unpack("<I", raw)[0]
    body = fh.read(50 * count)
    if len(body) != 50 * count:
        raise ValueError(
            f"the header claims {count} triangles but only {len(body) // 50} "
            f"are present -- the file looks truncated")
    data = np.frombuffer(body, dtype=np.dtype([
        ("normal", "<3f4"), ("v", "<3,3f4"), ("attr", "<u2")]), count=count)
    return np.asarray(data["v"], dtype=float)


def axisymmetry_report(x: np.ndarray, r: np.ndarray, theta: np.ndarray,
                       bins: int, outer: bool,
                       sectors: int = 24) -> Tuple[float, float]:
    """How much the extracted wall radius varies with azimuth.

    The obvious test -- the spread of all radii in an axial bin -- is wrong,
    and wrong in a way that looks right: a model of a *solid* has two surfaces
    at every station, so that spread is dominated by the wall thickness and
    a perfect body of revolution scores as badly asymmetric.

    What matters instead is whether the surface being extracted is the same
    radius all the way round.  So the bin is cut into azimuthal sectors, the
    same min-or-max rule is applied inside each one, and the spread of *those*
    is reported.  That is zero for any body of revolution however thick its
    wall, and large exactly when the model still has a flange, a boss or a
    feed port on it.
    """
    edges = np.linspace(x.min(), x.max(), bins + 1)
    index = np.clip(np.digitize(x, edges) - 1, 0, bins - 1)
    sector = np.clip(((theta + np.pi) / (2 * np.pi) * sectors).astype(int),
                     0, sectors - 1)
    spreads = []
    for b in range(bins):
        in_bin = index == b
        if not np.any(in_bin):
            continue
        per_sector = []
        for s in range(sectors):
            sel = r[in_bin & (sector == s)]
            if sel.size:
                per_sector.append(sel.max() if outer else sel.min())
        if len(per_sector) < sectors // 2:
            continue          # too sparse here to say anything
        arr = np.asarray(per_sector)
        hi = arr.max()
        if hi > 0:
            spreads.append((hi - arr.min()) / hi)
    if not spreads:
        return float("nan"), float("nan")
    return float(np.max(spreads)), float(np.median(spreads))


def extract(vertices: np.ndarray, axis: int, bins: int,
            outer: bool) -> Tuple[np.ndarray, np.ndarray, dict]:
    """Reduce a mesh to (x, r), one radius per axial bin."""
    points = vertices.reshape(-1, 3)
    x = points[:, axis]
    other = [i for i in range(3) if i != axis]
    # Centre the radial coordinates on the mesh's own axis rather than assuming
    # the model was drawn on the origin.
    centre = [float(np.mean(points[:, i])) for i in other]
    dy = points[:, other[0]] - centre[0]
    dz = points[:, other[1]] - centre[1]
    r = np.hypot(dy, dz)
    theta = np.arctan2(dz, dy)

    worst, median = axisymmetry_report(x, r, theta, bins, outer)

    edges = np.linspace(x.min(), x.max(), bins + 1)
    index = np.clip(np.digitize(x, edges) - 1, 0, bins - 1)
    xs, rs, dropped = [], [], 0
    for b in range(bins):
        sel = r[index == b]
        if sel.size == 0:
            dropped += 1
            continue
        xs.append(0.5 * (edges[b] + edges[b + 1]))
        rs.append(sel.max() if outer else sel.min())
    info = {"axisymmetry_worst": worst, "axisymmetry_median": median,
            "centre": centre, "empty_bins": dropped,
            "triangles": len(vertices)}
    return np.asarray(xs), np.asarray(rs), info


def trim_to_nozzle(x: np.ndarray, r: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Cut back to the largest run that contracts to the throat then expands.

    A model usually carries a little of something else at each end -- a flange
    face, the injector plate, a lip.  Those show up as a radius that turns the
    wrong way, which Ignis rejects outright.  Trimming to the monotone run on
    each side of the minimum is the smallest thing that makes the contour a
    nozzle, and the script says how much it removed.
    """
    throat = int(np.argmin(r))
    lo = throat
    while lo > 0 and r[lo - 1] >= r[lo]:
        lo -= 1
    hi = throat
    while hi < len(r) - 1 and r[hi + 1] >= r[hi]:
        hi += 1
    return x[lo:hi + 1], r[lo:hi + 1]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("stl", help="input STL file")
    ap.add_argument("-o", "--output", default="contour.csv", help="output CSV")
    ap.add_argument("--axis", choices=sorted(AXES), default="x",
                    help="axis of revolution in the model (default: x)")
    ap.add_argument("--bins", type=int, default=400,
                    help="axial stations to produce (default: 400)")
    ap.add_argument("--scale", type=float, default=1.0,
                    help="multiply model units by this to get metres; use 0.001 "
                         "for a model drawn in millimetres")
    ap.add_argument("--outer", action="store_true",
                    help="take the outer surface instead of the gas-side wall; "
                         "use for a model of the flow volume, not of the solid")
    ap.add_argument("--no-trim", action="store_true",
                    help="keep the contour exactly as extracted, including any "
                         "non-monotone ends Ignis will then reject")
    args = ap.parse_args(argv)

    try:
        vertices = read_stl(args.stl)
    except (OSError, ValueError) as exc:
        print(f"could not read {args.stl}: {exc}", file=sys.stderr)
        return 1

    x, r, info = extract(vertices, AXES[args.axis], args.bins, args.outer)
    if x.size < 9:
        print(f"only {x.size} usable stations came out of the mesh; Ignis needs "
              f"at least 9", file=sys.stderr)
        return 1

    print(f"{info['triangles']} triangles, axis {args.axis}, "
          f"radial centre ({info['centre'][0]:.6g}, {info['centre'][1]:.6g})")
    worst = info["axisymmetry_worst"]
    print(f"axisymmetry: worst bin spread {worst * 100:.2f} %, "
          f"median {info['axisymmetry_median'] * 100:.2f} %")
    if not (worst < 0.05):
        print("  WARNING: this model is not a body of revolution about that axis. "
              "A contour taken from it is not the engine's wall, and nothing "
              "computed from it will mean anything. Check --axis first.",
              file=sys.stderr)
    if info["empty_bins"]:
        print(f"  {info['empty_bins']} of {args.bins} bins were empty and skipped")

    surface = "outer" if args.outer else "inner (gas-side)"
    print(f"took the {surface} surface in each bin")

    n_before = x.size
    if not args.no_trim:
        x, r = trim_to_nozzle(x, r)
        if x.size != n_before:
            print(f"trimmed {n_before - x.size} of {n_before} stations from the "
                  f"ends to leave a contour that contracts then expands")
    if x.size < 9:
        print(f"trimming left only {x.size} stations; the model may not contain "
              f"a convergent-divergent passage", file=sys.stderr)
        return 1

    x = (x - x[0]) * args.scale
    r = r * args.scale
    throat = int(np.argmin(r))
    print(f"contour: {x.size} stations, length {x[-1] * 1e3:.2f} mm, "
          f"throat radius {r[throat] * 1e3:.3f} mm at x = {x[throat] * 1e3:.2f} mm")
    print(f"         contraction {(r[0] / r[throat]) ** 2:.3f}, "
          f"expansion {(r[-1] / r[throat]) ** 2:.3f}")

    with open(args.output, "w") as fh:
        fh.write(f"# wall contour extracted from {args.stl}\n")
        fh.write(f"# axis {args.axis}, {surface} surface, scale {args.scale}\n")
        fh.write(f"# axisymmetry: worst bin spread {worst * 100:.2f} %\n")
        fh.write("# x and r are in metres, measured from the first station\n")
        fh.write("x_m,r_m\n")
        for xi, ri in zip(x, r):
            fh.write(f"{xi:.9f},{ri:.9f}\n")
    print(f"wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
