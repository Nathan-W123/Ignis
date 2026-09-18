#!/usr/bin/env python3
"""Render the Ignis cover images and animations.

WHAT THIS MAKES
---------------
Three things, from one Ignis case:

  hero      a still of the engine firing -- the chamber, the throat, the bell
            and the exhaust plume, with the numbers that produced it
  flow      the engine running, with the flow made visible by tracers
  ignition  the engine lighting up, driven by the start-up transient

All three are emission renders of the solver's own field: the hue at every
point is Planck's law at the computed temperature and the brightness is the
computed density, integrated along the camera rays as an optically thin gray
emitter.  Nothing in the picture is painted on.

THE PLUME
---------
Ignis solves the nozzle and stops at the exit plane.  Everything downstream is
computed here by a separate axisymmetric Euler solver (`ignis_viz.plume`),
started from the Ignis exit state -- it is not part of the validated core, and
its assumptions and omissions are listed in that module.  The solve takes tens
of minutes, so it is cached; pass --solve to recompute it.

The case must be run first: see scripts/run_all.sh, or

    ignis_engine --config configs/methane_nominal.yaml --output results/...
"""
from __future__ import annotations

import argparse
import json
import os
import pickle
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "python"))

import numpy as np
import pandas as pd
from PIL import Image

from ignis_viz import animation, plume as plume_mod, render as R

# The nominal M1 is over-expanded at sea level and Ignis predicts the flow
# separates inside the bell, so a full-flowing sea-level plume would contradict
# the solver.  Three kilometres is the lowest round altitude where the altitude
# sweep reports the nozzle running full.
DEFAULT_ALTITUDE = 3000.0


def _ambient_pressure(altitude_csv: str, altitude: float):
    table = pd.read_csv(altitude_csv, comment="#")
    row = table.iloc[(table["altitude"] - altitude).abs().idxmin()]
    if bool(row["separation_predicted"]):
        raise SystemExit(
            f"Ignis predicts flow separation inside the nozzle at {altitude:.0f} m "
            f"(exit/ambient {row['exit_ambient_ratio']:.3f}).  Rendering an attached "
            f"plume there would contradict the solver; pick a higher altitude.")
    return float(row["ambient_pressure"]), row


def _load_plume(args, profile: pd.DataFrame, ambient: float):
    if not args.solve and os.path.exists(args.plume_cache):
        with open(args.plume_cache, "rb") as fh:
            return pickle.load(fh)
    exit_state = plume_mod.ExitState.from_profile(profile)
    started = time.time()
    field = plume_mod.solve(
        exit_state, ambient, nx=args.plume_nx, nr=args.plume_nr,
        length=args.plume_length, width=args.plume_width,
        max_steps=args.plume_steps, tolerance=args.plume_tolerance,
        report=lambda s, r, dt: print(f"  plume {s:6d}  residual {r:.3e}  "
                                      f"{time.time() - started:6.0f} s", flush=True))
    with open(args.plume_cache, "wb") as fh:
        pickle.dump(field, fh)
    return field


def _validate(exit_state, ambient: float, field) -> str:
    """Check the lip shock against the exact oblique-shock solution."""
    theory = plume_mod.lip_shock_theory(exit_state, ambient)
    if theory is None:
        return "plume is under-expanded; no attached lip shock to check"
    measured = plume_mod.measure_lip_shock(field)
    if measured is None:
        return "lip shock not located in the solution"
    beta, _, _ = theory
    return (f"lip shock  theory {np.degrees(beta):.2f} deg   "
            f"computed {np.degrees(measured):.2f} deg   "
            f"difference {abs(np.degrees(measured - beta)):.2f} deg")


def _jacket(case_prefix: str):
    """Cooling-channel count and land fraction, from the case's own config."""
    try:
        import yaml
        meta = json.load(open(f"{case_prefix}_engine.json"))
        path = meta.get("config_file", "")
        if not os.path.isabs(path):
            path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", path)
        cooling = yaml.safe_load(open(path)).get("cooling", {})
        channels = int(cooling.get("num_channels", 0))
        land = 1.0 - float(cooling.get("width_fraction", 0.5))
        return channels, land
    except Exception:
        return 0, 0.5


def _shell(args, field: R.Field, channels: int, land: float) -> R.Shell:
    """The nozzle wall as Ignis sizes it, drawn as a half section.

    The thickness is the hot-gas wall plus the coolant channel, which is what
    the cooling solver actually sizes; there is no structural close-out, no
    dome and no plumbing, because none of those are in the model.  The surfaces
    are kept very dark so the wall never veils the gas behind it, and the
    brightness goes to the sectioned face, which traces the contour and is what
    makes the bell read as a bell.
    """
    return R.Shell(transmittance=0.0, thickness=args.wall_thickness,
                   cut=True, ambient=0.002, diffuse=0.014, rim=0.045,
                   rim_power=4.0, section=0.60, glow=0.08)


def main() -> int:
    here = os.path.dirname(os.path.abspath(__file__))
    root = os.path.normpath(os.path.join(here, ".."))
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--case", default=os.path.join(root, "results/methane_nominal/m1"),
                    help="result prefix, e.g. results/methane_nominal/m1")
    ap.add_argument("--out", default=os.path.join(root, "results/figures"))
    ap.add_argument("--altitude", type=float, default=DEFAULT_ALTITUDE)
    ap.add_argument("--plume-cache", default=os.path.join(root, "results/plume_cache.pkl"))
    ap.add_argument("--solve", action="store_true", help="recompute the plume")
    ap.add_argument("--no-plume", action="store_true")
    ap.add_argument("--plume-nx", type=int, default=600)
    ap.add_argument("--plume-nr", type=int, default=200)
    ap.add_argument("--plume-length", type=float, default=10.0, help="exit radii")
    ap.add_argument("--plume-width", type=float, default=3.0, help="exit radii")
    ap.add_argument("--plume-steps", type=int, default=14000)
    ap.add_argument("--plume-tolerance", type=float, default=3.0e-3)
    ap.add_argument("--size", default="2400x1350")
    ap.add_argument("--max-step", type=float, default=2.2e-3,
                    help="ray-march sample spacing, metres")
    ap.add_argument("--wall-thickness", type=float, default=5.6e-3,
                    help="drawn wall: hot-gas wall + coolant channel, metres")
    ap.add_argument("--channels", type=int, default=0,
                    help="draw N cooling channels; 0 leaves the jacket "
                         "smooth, which is what 300 of them look like "
                         "at any sane image size")
    ap.add_argument("--head-depth", type=float, default=0.085,
                    help="injector dome depth, metres (a drawing choice)")
    ap.add_argument("--reach", type=float, default=3.2, help="plume in frame, exit radii")
    ap.add_argument("--exposure", type=float, default=0.52)
    ap.add_argument("--latitude", type=float, default=900.0)
    ap.add_argument("--tracer-gain", type=float, default=1.1)
    ap.add_argument("--tracers", type=int, default=90000)
    ap.add_argument("--margin", type=float, default=0.10)
    ap.add_argument("--dx", type=float, default=0.02)
    ap.add_argument("--dy", type=float, default=0.0)
    ap.add_argument("--move-width", type=int, default=1920)
    ap.add_argument("--move-height", type=int, default=1080)
    ap.add_argument("--move-frames", type=int, default=180)
    ap.add_argument("--only", choices=["hero", "flow", "ignition", "move"],
                    action="append")
    args = ap.parse_args()

    width, height = (int(v) for v in args.size.split("x"))
    want = set(args.only or ["hero", "flow", "move", "ignition"])
    os.makedirs(args.out, exist_ok=True)

    profile_csv = f"{args.case}_profile.csv"
    profile = pd.read_csv(profile_csv, comment="#")
    ambient, alt_row = _ambient_pressure(f"{args.case}_altitude.csv", args.altitude)
    print(f"altitude {args.altitude:.0f} m   ambient {ambient / 1e3:.2f} kPa   "
          f"{alt_row['regime']}   thrust {alt_row['thrust'] / 1e3:.1f} kN   "
          f"Isp {alt_row['isp']:.1f} s")

    field = R.Field.from_profile(profile, bins=320)
    channels, land = _jacket(args.case)
    shell = _shell(args, field, channels, land)
    print(f"jacket: {args.wall_thickness * 1e3:.1f} mm wall; "
          f"{channels} channels at {100 * land:.0f} % land, "
          f"{'drawn' if args.channels else 'not drawn'}")
    exit_state = plume_mod.ExitState.from_profile(profile)
    plume_field = None
    if not args.no_plume:
        plume_field = _load_plume(args, profile, ambient)
        print(" ", _validate(exit_state, ambient, plume_field))

    length, r_exit = field.length, float(field.radius[-1])
    extra = None
    if plume_field is not None:
        extra = R.plume_hull(plume_field, length)
        extra = extra[extra[:, 0] <= length + args.reach * r_exit]

    if "hero" in want:
        camera = R.shift_camera(
            R.framed_camera(field, width, height, pitch_deg=13.0, margin=args.margin,
                            extra_points=extra),
            ndc_dx=args.dx, ndc_dy=args.dy)
        started = time.time()
        scene = R.Scene(field, shell, plume_field,
                        plume_reach=(args.reach + 2.4) * r_exit)
        gas, wall, depth = R.march_scene(scene, camera, max_step=args.max_step)
        reference = R.exposure_reference(gas)

        from ignis_viz import tracers as T
        tr = T.Tracers(field, plume_field, count=args.tracers, lifetime=5.0e-4)
        dt = 1.0 / 30.0 / 3000.0
        for _ in range(40):
            tr.advance(dt / 6.0)
        gas = gas + T.streak(tr, camera, (height, width), dt,
                             gain=args.tracer_gain,
                             depth=depth).reshape(-1, 3) * reference

        image = R.to_uint8(R.compose(gas, wall, (height, width), exposure=args.exposure,
                                     latitude=args.latitude, reference=reference))
        plain = os.path.join(args.out, "17_engine_render.png")
        Image.fromarray(image).save(plain)
        print(f"  hero {time.time() - started:.0f} s -> {plain}")

    if "flow" in want:
        started = time.time()
        out = animation.flow_animation(
            profile_csv, plume_field, os.path.join(args.out, "19_engine_flow.mp4"),
            width=1920, height=1080, max_step=args.max_step, gain=args.tracer_gain,
            shell=shell, plume_reach=args.reach, exposure=args.exposure,
            latitude=args.latitude)
        print(f"  flow {time.time() - started:.0f} s -> {out}")

    if "move" in want:
        started = time.time()
        out = animation.camera_move_animation(
            profile_csv, plume_field, os.path.join(args.out, "21_engine_orbit.mp4"),
            width=args.move_width, height=args.move_height, frames=args.move_frames,
            max_step=args.max_step, gain=args.tracer_gain, shell=shell,
            exposure=args.exposure, latitude=args.latitude,
            progress=lambda k, n: (
                print(f"    frame {k + 1}/{n}  {time.time() - started:5.0f} s",
                      flush=True) if k % 10 == 0 else None))
        print(f"  move {time.time() - started:.0f} s -> {out}")

    if "ignition" in want:
        started = time.time()
        out = animation.firing_animation(
            os.path.join(root, "results/startup_transient/startup_transient_transient.csv"),
            profile_csv, os.path.join(args.out, "20_engine_ignition.mp4"),
            width=1920, height=1080, steps=1800, bins=320, shell=shell)
        print(f"  ignition {time.time() - started:.0f} s -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
