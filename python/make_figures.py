#!/usr/bin/env python3
"""Generate every Ignis figure from solver output.

The figures are drawn only from files the C++ applications wrote.  Run
``scripts/run_all.sh`` first (or the individual ``ignis_*`` commands) so the
inputs exist, then::

    python3 python/make_figures.py --results results --figures results/figures

Missing inputs are reported by name and skipped rather than fabricated, and the
exit status is non-zero if any requested figure could not be made.
"""
from __future__ import annotations

import argparse
import os
import sys
import traceback

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from ignis_viz import animation, figures, style  # noqa: E402


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", default="results", help="directory holding solver output")
    ap.add_argument("--figures", default="results/figures", help="where to write the PNGs")
    ap.add_argument("--no-animation", action="store_true", help="skip the animation")
    ap.add_argument("--only", nargs="*", default=None, help="only these figure names")
    args = ap.parse_args(argv)

    style.apply()
    R, F = args.results, args.figures
    os.makedirs(F, exist_ok=True)

    jobs = {
        "composition_vs_mixture_ratio": lambda: figures.mixture_ratio_composition(
            f"{R}/mixture_ratio_sweep/mixture_ratio_sweep_sweep.csv",
            f"{R}/mixture_ratio_sweep/mixture_ratio_sweep_composition_vs_mr.csv",
            f"{F}/01_composition_vs_mixture_ratio.png"),
        "performance_vs_mixture_ratio": lambda: figures.mixture_ratio_performance(
            f"{R}/mixture_ratio_sweep/mixture_ratio_sweep_sweep.csv",
            f"{F}/02_performance_vs_mixture_ratio.png"),
        "nozzle_contour_mach": lambda: figures.nozzle_contour_mach(
            f"{R}/methane_nominal/m1_profile.csv", f"{F}/03_nozzle_contour_mach.png"),
        "axial_profiles": lambda: figures.axial_profiles(
            f"{R}/methane_nominal/m1_profile.csv", f"{F}/04_axial_profiles.png"),
        "altitude_performance": lambda: figures.altitude_performance(
            f"{R}/methane_vacuum/m1vac_altitude.csv", f"{F}/05_altitude_performance.png"),
        "expansion_trade": lambda: figures.expansion_trade(
            f"{R}/expansion_ratio_sweep/expansion_ratio_sweep_sweep.csv",
            f"{F}/06_expansion_ratio_trade.png"),
        "thermal_profiles": lambda: figures.thermal_profiles(
            f"{R}/methane_nominal/m1_thermal.csv", f"{F}/07_thermal_profiles.png"),
        "cooling_map": lambda: figures.cooling_map(
            f"{R}/cooling_sweep/cooling_sweep_sweep.csv", f"{F}/08_cooling_design_space.png"),
        "transient_history": lambda: figures.transient_history(
            f"{R}/startup_transient/startup_transient_transient.csv",
            f"{F}/09_startup_transient.png"),
        "optimization_convergence": lambda: figures.optimization_convergence(
            f"{R}/optimization/optimization_optimization_history.csv",
            f"{R}/optimization/optimization_sweep.json",
            f"{F}/10_optimization_convergence.png"),
        "monte_carlo_histograms": lambda: figures.monte_carlo_histograms(
            f"{R}/monte_carlo/monte_carlo_mc_samples.csv",
            f"{R}/monte_carlo/monte_carlo_mc.json", f"{F}/11_monte_carlo_histograms.png"),
        "sensitivity_tornado": lambda: figures.sensitivity_tornado(
            f"{R}/monte_carlo/monte_carlo_mc_sensitivity.csv",
            f"{F}/12_sensitivity_ranking.png",
            ["performance.isp", "performance.thrust", "cooling.max_wall_temperature",
             "cooling.pressure_drop"]),
        "monte_carlo_scatter": lambda: figures.monte_carlo_scatter(
            f"{R}/monte_carlo/monte_carlo_mc_samples.csv",
            f"{F}/13_monte_carlo_scatter.png"),
        "hydrogen_axial_profiles": lambda: figures.axial_profiles(
            f"{R}/hydrogen_nominal/h1_profile.csv", f"{F}/14_hydrogen_axial_profiles.png"),
        "hydrogen_thermal": lambda: figures.thermal_profiles(
            f"{R}/hydrogen_nominal/h1_thermal.csv", f"{F}/15_hydrogen_thermal.png"),
    }
    if not args.no_animation:
        jobs["startup_animation"] = lambda: animation.startup_animation(
            f"{R}/startup_transient/startup_transient_transient.csv",
            f"{R}/methane_nominal/m1_profile.csv", f"{F}/00_startup_animation.gif")

    wanted = args.only or list(jobs)
    unknown = [n for n in wanted if n not in jobs]
    if unknown:
        print(f"unknown figure name(s): {', '.join(unknown)}", file=sys.stderr)
        print("available: " + ", ".join(jobs), file=sys.stderr)
        return 2

    made, failed = [], []
    for name in wanted:
        try:
            path = jobs[name]()
            made.append(path)
            print(f"  wrote {path}")
        except FileNotFoundError as exc:
            failed.append((name, str(exc)))
            print(f"  SKIPPED {name}: {exc}", file=sys.stderr)
        except Exception:
            failed.append((name, traceback.format_exc(limit=2)))
            print(f"  FAILED {name}:\n{traceback.format_exc(limit=3)}", file=sys.stderr)

    print(f"\n{len(made)} figure(s) written to {F}")
    if failed:
        print(f"{len(failed)} figure(s) could not be made:", file=sys.stderr)
        for name, why in failed:
            print(f"  {name}: {why.splitlines()[-1] if why else ''}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
