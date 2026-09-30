"""Run the Ignis solver for the Explorer and collect everything it produced.

The Explorer does not reimplement any physics.  It writes a configuration,
invokes the same `ignis_engine` and `ignis_equilibrium` binaries the command
line uses, and reads back their JSON and CSV.  A design that the Explorer shows
is a design the solver actually solved; if the solver refuses, the Explorer
shows the refusal rather than a plausible-looking number.
"""
from __future__ import annotations

import csv
import json
import os
import subprocess
import tempfile
from dataclasses import asdict, dataclass, field, replace
from typing import Dict, List, Optional

PROPELLANTS = {
    "LOX / CH4": {
        "oxidizer": "LOX",
        "fuel": "LCH4",
        "oxidizer_temperature": 90.18,
        "fuel_temperature": 111.66,
        "coolant": "methane",
        "elements": ["C", "H", "O"],
        "mr_range": (1.5, 6.0),
        "mr_default": 3.4,
        "coolant_inlet_temperature": 111.66,
    },
    "LOX / H2": {
        "oxidizer": "LOX",
        "fuel": "LH2",
        "oxidizer_temperature": 90.18,
        "fuel_temperature": 20.27,
        "coolant": "hydrogen",
        "elements": ["H", "O"],
        "mr_range": (3.0, 9.0),
        "mr_default": 5.5,
        "coolant_inlet_temperature": 25.0,
    },
}


@dataclass
class Design:
    """Everything the Explorer lets the user change."""

    propellant: str = "LOX / CH4"
    chamber_pressure: float = 5.5e6     # Pa
    mixture_ratio: float = 3.4
    throat_radius: float = 0.070        # m
    expansion_ratio: float = 20.0
    altitude: float = 0.0               # m
    # A tabulated wall contour to use instead of the analytic bell.  When set,
    # the throat radius and expansion ratio controls stop having any effect --
    # both are measured from the contour -- and the UI says so rather than
    # leaving live-looking controls that do nothing.
    contour_file: str = ""
    composition: str = "equilibrium"    # equilibrium (shifting) | frozen
    eta_c_star: float = 0.96
    cooling: bool = True
    # Jacket geometry, kept out of the main panel but part of the design.
    num_channels: int = 300
    channel_height: float = 5.0e-3      # m
    wall_thickness: float = 0.6e-3      # m
    coolant_inlet_pressure: float = 15.0e6  # Pa

    def with_propellant(self, name: str) -> "Design":
        """Switch propellant, moving the mixture ratio to that pair's default."""
        p = PROPELLANTS[name]
        return replace(self, propellant=name, mixture_ratio=p["mr_default"],
                       coolant_inlet_pressure=15.0e6 if name == "LOX / CH4" else 12.0e6)


# The shipped designs.  They live here rather than in the window because the
# case-freezing tool has to see them without importing a GUI toolkit.
#
# "10 km" is first because it is where the M1 is worth looking at: at sea level
# Ignis predicts the nozzle separates, so the Flow tab rightly refuses to march
# an attached plume there -- correct, but a poor first thing to be shown.
PRESETS: Dict[str, Design] = {
    "Ignis-M1  10 km": Design(altitude=10_000.0),
    "Ignis-M1  sea level": Design(),
    "Ignis-M1  vacuum": Design(expansion_ratio=45.0, altitude=80_000.0),
    "Ignis-H1  upper stage": Design(propellant="LOX / H2", mixture_ratio=5.5,
                                    throat_radius=0.060, expansion_ratio=60.0,
                                    altitude=80_000.0, eta_c_star=0.97,
                                    num_channels=240, channel_height=6.0e-3,
                                    wall_thickness=0.7e-3,
                                    coolant_inlet_pressure=12.0e6),
}


@dataclass
class Result:
    """A solved design, or the reason it could not be solved."""

    ok: bool = False
    error: str = ""
    design: Optional[Design] = None
    scalars: Dict[str, float] = field(default_factory=dict)
    strings: Dict[str, str] = field(default_factory=dict)
    warnings: List[str] = field(default_factory=list)
    profile: Dict[str, List[float]] = field(default_factory=dict)
    thermal: Dict[str, List[float]] = field(default_factory=dict)
    altitude: Dict[str, List[float]] = field(default_factory=dict)
    composition: List[tuple] = field(default_factory=list)   # (species, mole fraction)
    version: str = ""
    replayed: bool = False     # True when FrozenSolver answered, not the binaries

    def get(self, key: str, default: float = float("nan")) -> float:
        v = self.scalars.get(key, default)
        return default if v is None else v


def _read_csv(path: str) -> Dict[str, List[float]]:
    """Read an Ignis CSV, skipping its `# units:` header line."""
    if not os.path.exists(path):
        return {}
    with open(path, newline="") as fh:
        rows = [r for r in fh if not r.startswith("#")]
    reader = csv.DictReader(rows)
    out: Dict[str, List[float]] = {k: [] for k in (reader.fieldnames or [])}
    for row in reader:
        for k, v in row.items():
            if k is None:
                continue
            try:
                out[k].append(float(v))
            except (TypeError, ValueError):
                out[k].append(float("nan"))
    return out


def _config_text(d: Design) -> str:
    p = PROPELLANTS[d.propellant]
    elements = ", ".join(p["elements"])
    lines = [
        "# Written by the Ignis Engine Explorer.  Every number here came from a",
        "# control in the UI; nothing is hidden from the solver.",
        "name: explorer",
        f"species: {{ elements: [{elements}] }}",
        "propellants:",
        f"  oxidizer: {p['oxidizer']}",
        f"  fuel: {p['fuel']}",
        f"  oxidizer_temperature: {p['oxidizer_temperature']}",
        f"  fuel_temperature: {p['fuel_temperature']}",
        f"  mixture_ratio: {d.mixture_ratio!r}",
        "chamber:",
        f"  pressure: {d.chamber_pressure!r}",
        f"  eta_c_star: {d.eta_c_star!r}",
        f"  composition: {d.composition}",
        "nozzle:",
    ]
    if d.contour_file:
        # Everything the analytic path takes as input, the tabulated path
        # measures, so emitting the analytic keys here would be writing down
        # numbers the solver is about to ignore.  Only the fallback curvature
        # ratio and the station count still mean anything.
        lines += [
            f"  contour_file: {d.contour_file}",
            "  throat_upstream_ratio: 1.5",
            "  stations: 240",
        ]
    else:
        lines += [
            f"  throat_radius: {d.throat_radius!r}",
            "  contraction_ratio: 2.8",
            "  chamber_length: 0.22",
            "  converging_half_angle: 30.0",
            "  chamber_fillet_ratio: 0.5",
            "  throat_upstream_ratio: 1.5",
            "  throat_downstream_ratio: 0.382",
            f"  expansion_ratio: {d.expansion_ratio!r}",
            "  type: bell",
            "  bell_length_fraction: 0.8",
            "  bell_initial_angle: 33.0",
            "  bell_exit_angle: 8.0",
            "  stations: 240",
            "  bartz_curvature: mean",
        ]
    lines += [
        "performance:",
        f"  altitude: {d.altitude!r}",
        "  auto_divergence: true",
        "  eta_nozzle: 1.0",
        "  separation: summerfield",
        "  resolve_internal_shocks: true",
        "  ascent_profile:",
        "    altitudes: [0, 5000, 10000, 20000, 40000]",
        "    weights:   [0.30, 0.25, 0.20, 0.15, 0.10]",
    ]
    if d.cooling:
        lines += [
            "cooling:",
            "  enabled: true",
            f"  coolant: {p['coolant']}",
            "  material: CuCrZr",
            f"  num_channels: {int(d.num_channels)}",
            "  width_mode: fraction_of_pitch",
            "  width_fraction: 0.50",
            f"  channel_height: {d.channel_height!r}",
            f"  wall_thickness: {d.wall_thickness!r}",
            "  roughness: 5.0e-6",
            "  coolant_fuel_fraction: 1.0",
            f"  inlet_temperature: {p['coolant_inlet_temperature']!r}",
            f"  inlet_pressure: {d.coolant_inlet_pressure!r}",
            "  counterflow: true",
            "  x_end_area_ratio: 8.0",
            "  nusselt_correlation: dittus-boelter",
            "  segments: 160",
        ]
    lines += [
        "feed:",
        "  enabled: true",
        "  mode: sizing",
        "  injector_stiffness: 0.20",
        "  pump_efficiency: 0.65",
        "  pump_inlet_pressure: 3.0e5",
        "  oxidizer: { line_length: 1.5, line_diameter: 0.10, fitting_k: 2.5,",
        "              viscosity: 1.9e-4, injector_cd: 0.78 }",
        "  fuel:     { line_length: 1.8, line_diameter: 0.07, fitting_k: 3.0,",
        "              viscosity: 1.2e-4, injector_cd: 0.75 }",
    ]
    return "\n".join(lines) + "\n"


class Solver:
    """Locates the Ignis binaries and runs a design through them."""

    live = True

    def __init__(self, repo_root: str, build_dir: Optional[str] = None) -> None:
        self.repo_root = os.path.abspath(repo_root)
        self.build_dir = build_dir or os.environ.get(
            "IGNIS_BUILD_DIR", os.path.join(self.repo_root, "build"))
        self.bin_dir = os.path.join(self.build_dir, "bin")

    def binary(self, name: str) -> str:
        """Path to one solver binary, wherever this platform's build put it.

        `build/bin` is where every generator is told to write, but a build made
        before that was pinned -- or one driven by hand -- leaves Visual Studio
        and Xcode output under a per-configuration subdirectory instead.  Those
        are searched too, so an existing build keeps working rather than
        producing a "binaries not found" dialog next to a directory full of
        binaries.
        """
        exe = name + (".exe" if os.name == "nt" else "")
        candidates = [os.path.join(self.bin_dir, exe)]
        candidates += [os.path.join(self.bin_dir, cfg, exe)
                       for cfg in ("Release", "RelWithDebInfo", "Debug", "MinSizeRel")]
        for path in candidates:
            if os.path.exists(path):
                return path
        return candidates[0]

    def available(self) -> bool:
        return os.path.exists(self.binary("ignis_engine"))

    def missing_message(self) -> str:
        if os.name == "nt":
            how = ("Build them first, from the repository root:\n"
                   "    cmake -S . -B build -DCMAKE_BUILD_TYPE=Release\n"
                   "    cmake --build build --config Release")
        else:
            how = "Build them first:\n    ./scripts/build.sh"
        return ("Could not find the Ignis binaries.\n\n"
                f"Looked in: {self.bin_dir}\n"
                "(and in Release/ RelWithDebInfo/ Debug/ below it)\n\n"
                f"{how}\n\n"
                "or point the Explorer at an existing build with the "
                "IGNIS_BUILD_DIR environment variable.")

    def run(self, d: Design) -> Result:
        if not self.available():
            return Result(ok=False, design=d, error=self.missing_message())
        with tempfile.TemporaryDirectory(prefix="ignis-explorer-") as tmp:
            cfg = os.path.join(tmp, "explorer.yaml")
            with open(cfg, "w") as fh:
                fh.write(_config_text(d))
            try:
                self._exec([self.binary("ignis_engine"), "--config", cfg,
                            "--output", tmp, "--prefix", "x", "--contour", "--quiet"])
                self._exec([self.binary("ignis_equilibrium"), "--config", cfg,
                            "--output", tmp, "--prefix", "x", "--quiet"])
                # The altitude curve comes from the same validated atmosphere
                # the command line uses, not from a second implementation here.
                self._exec([self.binary("ignis_nozzle"), "--config", cfg,
                            "--output", tmp, "--prefix", "x", "--quiet",
                            "--altitude-min", "0", "--altitude-max", "80000",
                            "--points", "81"])
            except RuntimeError as exc:
                return Result(ok=False, design=d, error=str(exc))
            return self._collect(d, tmp)

    @staticmethod
    def _exec(cmd: List[str]) -> None:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=180)
        if proc.returncode != 0:
            # Ignis reports failures as a single actionable line on stderr.
            msg = (proc.stderr or proc.stdout or "").strip()
            raise RuntimeError(msg or f"{os.path.basename(cmd[0])} exited {proc.returncode}")

    def _collect(self, d: Design, tmp: str) -> Result:
        with open(os.path.join(tmp, "x_engine.json")) as fh:
            doc = json.load(fh)
        r = doc["result"]
        res = Result(ok=True, design=d, version=doc.get("git", ""))
        ch, pf, geo = r["chamber"], r["performance"], r.get("geometry", {})
        cool = r.get("cooling", {})
        feed = r.get("feed", {})

        def take(prefix: str, src: dict) -> None:
            for k, v in src.items():
                if isinstance(v, (int, float)):
                    res.scalars[f"{prefix}.{k}"] = float(v)
                elif isinstance(v, str):
                    res.strings[f"{prefix}.{k}"] = v
                elif isinstance(v, bool):
                    res.scalars[f"{prefix}.{k}"] = 1.0 if v else 0.0

        take("chamber", ch)
        take("performance", pf)
        take("geometry", geo)
        take("cooling", cool)
        take("feed", feed)
        res.warnings = list(cool.get("warnings", []) or [])

        res.profile = _read_csv(os.path.join(tmp, "x_profile.csv"))
        res.thermal = _read_csv(os.path.join(tmp, "x_thermal.csv"))
        res.altitude = _read_csv(os.path.join(tmp, "x_altitude.csv"))

        eq = _read_csv(os.path.join(tmp, "x_equilibrium.csv"))
        comp = []
        if eq:
            names = _read_csv_strings(os.path.join(tmp, "x_equilibrium.csv"), "species")
            fracs = eq.get("mole_fraction", [])
            for name, x in zip(names, fracs):
                if x == x and x > 0.0:
                    comp.append((name, x))
            comp.sort(key=lambda t: t[1], reverse=True)
        res.composition = comp
        return res


def _read_csv_strings(path: str, column: str) -> List[str]:
    if not os.path.exists(path):
        return []
    with open(path, newline="") as fh:
        rows = [r for r in fh if not r.startswith("#")]
    return [row[column] for row in csv.DictReader(rows) if column in row]


# ---------------------------------------------------------------------------
# Frozen results: every pane without a compiler
# ---------------------------------------------------------------------------

RESULTS_DIR = os.path.join(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
    "data", "results")


def result_to_dict(r: Result) -> dict:
    """Everything a Result carries, as plain JSON types."""
    return {
        "schema": 1,
        "ok": r.ok, "error": r.error, "version": r.version,
        "design": asdict(r.design) if r.design is not None else None,
        "scalars": r.scalars, "strings": r.strings, "warnings": r.warnings,
        "profile": r.profile, "thermal": r.thermal, "altitude": r.altitude,
        "composition": [list(c) for c in r.composition],
    }


def result_from_dict(d: dict) -> Result:
    if d.get("schema") != 1:
        raise ValueError(f"unknown frozen-result schema {d.get('schema')!r}")
    return Result(
        ok=d["ok"], error=d.get("error", ""), version=d.get("version", ""),
        design=Design(**d["design"]) if d.get("design") else None,
        scalars=d.get("scalars", {}), strings=d.get("strings", {}),
        warnings=list(d.get("warnings", [])),
        profile=d.get("profile", {}), thermal=d.get("thermal", {}),
        altitude=d.get("altitude", {}),
        composition=[tuple(c) for c in d.get("composition", [])])


def _same_design(a: Design, b: Design) -> bool:
    """Field-by-field, with a tolerance on floats.

    Exact equality is wrong here: the design panel rebuilds a Design from spin
    boxes, so 5.0 mm arrives as 5.0 * 1e-3, which is not bit-for-bit the
    literal 5.0e-3 a preset was written with.  One part in a billion is far
    tighter than any control can resolve and far looser than that rounding.
    """
    for key, va in asdict(a).items():
        vb = getattr(b, key)
        if isinstance(va, float) or isinstance(vb, float):
            if abs(float(va) - float(vb)) > 1e-9 * max(1.0, abs(float(va))):
                return False
        elif va != vb:
            return False
    return True


class FrozenSolver:
    """Stands in for `Solver` when the compiled binaries are missing.

    It answers from results that `tools/make_case_library.py` produced by
    running the real solver on each preset -- through this module's own
    `Solver`, so a replayed Result is identical to what a live solve of that
    design returns.  It answers ONLY for those designs.  Anything else comes
    back as a refusal naming what would be needed, never as the nearest
    preset's numbers dressed up as the answer to a different question.
    """

    live = False

    def __init__(self, directory: str = RESULTS_DIR) -> None:
        self.directory = directory
        self._frozen: List[Result] = []
        if os.path.isdir(directory):
            for name in sorted(os.listdir(directory)):
                if name.endswith(".json"):
                    try:
                        with open(os.path.join(directory, name)) as fh:
                            self._frozen.append(result_from_dict(json.load(fh)))
                    except (OSError, ValueError, KeyError, TypeError):
                        continue

    def available(self) -> bool:
        return bool(self._frozen)

    def missing_message(self) -> str:
        return f"no frozen results in {self.directory}"

    def run(self, d: Design) -> Result:
        for r in self._frozen:
            if r.design is not None and _same_design(r.design, d):
                return replace(r, replayed=True)
        return Result(ok=False, design=d, error=(
            "this design is not one of the saved presets, and solving a new "
            "one needs the compiled solver.\n"
            "Pick a preset from the ribbon, or build the binaries."))
