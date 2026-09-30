#!/usr/bin/env python3
"""Freeze solved engines into a small library the Explorer can open offline.

WHY THIS EXISTS
---------------
The Explorer shells out to the compiled solver, so it needs a C++ toolchain
before it will show anything.  The plume march does not: `ignis_viz.plume` is
NumPy, and everything it needs from the nozzle is eight numbers at the exit
plane.  So a library of solved exit states lets the Flow tab run -- with the
plume genuinely computed, live, on whatever machine opens it -- without a
compiler anywhere.

WHAT IS AND IS NOT LIVE
-----------------------
What is frozen here is the *upstream* answer: chamber equilibrium, the nozzle
expansion, the exit state, the performance figures.  Those were solved by the
binaries at the commit recorded in each file and are reproduced, not recomputed,
when the Explorer opens one.

What stays live is the plume: the ambient pressure is a free parameter, and the
axisymmetric Euler march is run from scratch every time.  That is the part that
moves on screen, and it is real physics each time, not a recording.

The Explorer says which is which rather than presenting the whole thing as one
live calculation.

USAGE
-----
    python3 tools/make_case_library.py            # rebuild every case
    python3 tools/make_case_library.py --list     # show what would be built
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import subprocess
import sys
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT_DIR = os.path.join(REPO, "data", "cases")

# (file stem, config, human name, one-line description)
CASES = [
    ("m1", "methane_nominal.yaml", "Ignis-M1",
     "LOX/CH4, 5.5 MPa, O/F 3.4, expansion 20 - the reference engine"),
    ("m1_vacuum", "methane_vacuum.yaml", "Ignis-M1 vacuum",
     "The same chamber with a vacuum-optimised bell"),
    ("h1", "hydrogen_nominal.yaml", "Ignis-H1",
     "LOX/H2, the hydrogen engine"),
]

# Exit-plane quantities the plume solver needs, and nothing else.
EXIT_KEYS = ("radius", "pressure", "temperature", "density", "velocity",
             "mach", "gamma_s", "molar_mass")


def read_csv(path: str) -> dict:
    with open(path) as fh:
        rows = [line for line in fh if not line.startswith("#")]
    out: dict = {}
    for row in csv.DictReader(rows):
        for key, value in row.items():
            try:
                out.setdefault(key, []).append(float(value))
            except (TypeError, ValueError):
                out.setdefault(key, []).append(value)
    return out


def build_case(stem: str, config: str, name: str, description: str,
               binary: str) -> dict:
    with tempfile.TemporaryDirectory() as tmp:
        subprocess.run(
            [binary, "-c", os.path.join(REPO, "configs", config),
             "-o", tmp, "--prefix", "case", "--quiet"],
            check=True, capture_output=True, text=True)
        profile = read_csv(os.path.join(tmp, "case_profile.csv"))
        engine = json.load(open(os.path.join(tmp, "case_engine.json")))

        # The altitude table comes from ignis_nozzle, which sweeps ambient.
        subprocess.run(
            [binary.replace("ignis_engine", "ignis_nozzle"),
             "-c", os.path.join(REPO, "configs", config), "-o", tmp,
             "--prefix", "alt", "--quiet",
             "--altitude-min", "0", "--altitude-max", "80000", "--points", "81"],
            check=True, capture_output=True, text=True)
        altitude = read_csv(os.path.join(tmp, "alt_altitude.csv"))

    result = engine["result"]
    perf, chamber = result["performance"], result["chamber"]
    exit_state = {key: profile[key][-1] for key in EXIT_KEYS}

    # The contour is only drawn, so it is thinned: the viewer cannot resolve
    # 400 stations across a few hundred pixels of engine silhouette.
    step = max(1, len(profile["x"]) // 120)
    contour = [[round(x, 6), round(r, 6)]
               for x, r in zip(profile["x"][::step], profile["radius"][::step])]
    if contour[-1][0] != round(profile["x"][-1], 6):
        contour.append([round(profile["x"][-1], 6),
                        round(profile["radius"][-1], 6)])

    return {
        "schema": 1,
        "name": name,
        "description": description,
        "config": config,
        # The binary's own stamp, not a `git describe` of the tree now: this
        # tool writes several files in turn, so by the second one the tree is
        # dirty with the first, and a live describe would blame the solve on
        # changes that did not exist when it ran.
        "solved_by": engine.get("git", "unknown"),
        "note": ("The exit state and performance below were solved by the Ignis "
                 "binaries and are replayed, not recomputed. The plume is "
                 "solved live from this exit state every time it is run."),
        "exit": exit_state,
        "performance": {
            "thrust": perf.get("thrust"),
            "isp": perf.get("isp"),
            "chamber_pressure": chamber.get("pressure"),
            "chamber_temperature": chamber.get("temperature"),
            "mixture_ratio": result.get("propellants", {}).get("mixture_ratio"),
            "expansion_ratio": result.get("geometry", {}).get("expansion_ratio"),
            "ambient_pressure": perf.get("ambient_pressure"),
        },
        "contour": contour,
        "altitude_table": {
            "altitude": [round(a, 1) for a in altitude["altitude"]],
            "ambient_pressure": [round(p, 4) for p in altitude["ambient_pressure"]],
            "separation_predicted": [int(s) for s in altitude["separation_predicted"]],
        },
    }



def freeze_presets(build_dir: str) -> int:
    """Solve every Explorer preset and write the whole Result to data/results.

    This goes through the Explorer's own `Solver`, not a parallel reading of
    the binaries' output, so what `FrozenSolver` later replays is the object a
    live solve of the same design returns -- the property grid, every chart
    tab and every tile read identical numbers either way.
    """
    sys.path.insert(0, os.path.join(REPO, "python"))
    from ignis_explorer.solver import (PRESETS, RESULTS_DIR, Solver,
                                       result_to_dict)

    solver = Solver(REPO, build_dir)
    os.makedirs(RESULTS_DIR, exist_ok=True)
    for name, design in PRESETS.items():
        result = solver.run(design)
        if not result.ok:
            print(f"{name}: the solver refused it\n{result.error}", file=sys.stderr)
            return 1
        stem = "".join(c if c.isalnum() else "_" for c in name.lower()).strip("_")
        stem = "_".join(part for part in stem.split("_") if part)
        path = os.path.join(RESULTS_DIR, f"{stem}.json")
        doc = result_to_dict(result)
        doc["preset"] = name
        with open(path, "w") as fh:
            json.dump(doc, fh, separators=(",", ":"))
            fh.write("\n")
        print(f"{name:24s} {result.get('performance.thrust') / 1e3:7.1f} kN  "
              f"{os.path.getsize(path) / 1024:6.1f} kB  -> "
              f"{os.path.relpath(path, REPO)}")
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--build-dir", default=os.path.join(REPO, "build"))
    ap.add_argument("--list", action="store_true", help="show the cases and exit")
    args = ap.parse_args(argv)

    if args.list:
        for stem, config, name, description in CASES:
            print(f"{stem:12s} {name:20s} {config:24s} {description}")
        return 0

    binary = os.path.join(args.build_dir, "bin", "ignis_engine")
    if not os.path.exists(binary):
        print(f"the solver is not built: {binary} does not exist", file=sys.stderr)
        return 1

    if freeze_presets(args.build_dir) != 0:
        return 1

    os.makedirs(OUT_DIR, exist_ok=True)
    for stem, config, name, description in CASES:
        try:
            case = build_case(stem, config, name, description, binary)
        except subprocess.CalledProcessError as exc:
            print(f"{stem}: solver failed\n{exc.stderr}", file=sys.stderr)
            return 1
        path = os.path.join(OUT_DIR, f"{stem}.json")
        with open(path, "w") as fh:
            json.dump(case, fh, indent=1)
            fh.write("\n")
        size = os.path.getsize(path)
        print(f"{stem:12s} {name:20s} exit M {case['exit']['mach']:.2f}  "
              f"{case['performance']['thrust'] / 1e3:7.1f} kN  "
              f"{size / 1024:5.1f} kB  -> {os.path.relpath(path, REPO)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
