#!/usr/bin/env python3
"""Generate Ignis coolant / propellant real-fluid property tables.

WHY A TABLE
-----------
Regeneratively cooled channels run methane or hydrogen at 5-30 MPa and 100-900
K -- deep in the supercritical region where cubic equations of state and the
classical residual-property correlations (Jossi-Stiel-Thodos, Stiel-Thodos)
are only good to 10-50 %.  Ignis therefore ships a tabulated property database
built from the published *reference* equations of state and transport
correlations, evaluated with CoolProp:

  Methane
    EOS         U. Setzmann and W. Wagner, "A New Equation of State and Tables
                of Thermodynamic Properties for Methane Covering the Range from
                the Melting Line to 625 K at Pressures up to 1000 MPa",
                J. Phys. Chem. Ref. Data 20, 1061-1155 (1991).
    viscosity   E. W. Lemmon and R. T. Jacobsen, Int. J. Thermophys. 25, 21-69
                (2004) / CoolProp default methane correlation.
    conductivity  same reference.

  Hydrogen (normal)
    EOS         J. W. Leachman, R. T. Jacobsen, S. G. Penoncello and
                E. W. Lemmon, "Fundamental Equations of State for Parahydrogen,
                Normal Hydrogen, and Orthohydrogen",
                J. Phys. Chem. Ref. Data 38, 721-748 (2009).
    viscosity / conductivity  CoolProp defaults (Muzny et al. 2013,
                Assael et al. 2011).

  Oxygen
    EOS         R. Schmidt and W. Wagner, Fluid Phase Equilibria 19, 175-200
                (1985).

CoolProp itself: I. H. Bell, J. Wronski, S. Quoilin and V. Lemort,
Ind. Eng. Chem. Res. 53, 2498-2508 (2014).  BSD licensed.

The table is generated once and committed, so the C++ core has no runtime
dependency on CoolProp and the results are reproducible from a clean checkout.
`ignis::PengRobinsonFluid` provides an independent analytic model that the test
suite compares against these tables, so the size of the approximation is
measured rather than assumed.

USAGE
-----
    python3 tools/build_coolant_tables.py -o data/coolants
"""
from __future__ import annotations

import argparse
import datetime as _dt
import os
import sys

import numpy as np

try:
    import CoolProp.CoolProp as CP
except ImportError:  # pragma: no cover
    sys.exit("CoolProp is required to regenerate the tables: pip install CoolProp")

# fluid -> (CoolProp name, lowest T, highest T, pressure grid)
#
# The temperature axis is deliberately non-uniform: transport properties and cp
# vary steeply through the near-critical and pseudo-critical region and are
# smooth well above it, so the grid is refined below 1.5 T_crit and coarsened
# at high temperature.  This keeps the interpolation error small without
# inflating the committed file.
FLUIDS = {
    "methane":  dict(cp_name="Methane",  t_min=95.0, t_max=1100.0, p=(0.5e6, 40.0e6, 41)),
    "hydrogen": dict(cp_name="Hydrogen", t_min=22.0, t_max=1100.0, p=(0.5e6, 40.0e6, 41)),
    "oxygen":   dict(cp_name="Oxygen",   t_min=60.0, t_max=1000.0, p=(0.5e6, 40.0e6, 41)),
}


def temperature_axis(t_min, t_max, t_crit):
    """Refined below 1.5 T_crit, moderate to 2.5 T_crit, coarse above."""
    knees = [(t_min, min(1.5 * t_crit, t_max), 2.0),
             (min(1.5 * t_crit, t_max), min(2.5 * t_crit, t_max), 8.0),
             (min(2.5 * t_crit, t_max), t_max, 25.0)]
    pts = []
    for lo, hi, step in knees:
        if hi <= lo:
            continue
        n = max(1, int(round((hi - lo) / step)))
        pts.extend(lo + (hi - lo) * i / n for i in range(n))
    pts.append(t_max)
    return np.array(sorted(set(round(v, 6) for v in pts)))


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-o", "--output-dir", default="data/coolants")
    args = ap.parse_args(argv)
    os.makedirs(args.output_dir, exist_ok=True)
    stamp = _dt.datetime.now(_dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")

    for name, spec in FLUIDS.items():
        cpn = spec["cp_name"]
        Tc = CP.PropsSI("Tcrit", cpn)
        T = temperature_axis(spec["t_min"], spec["t_max"], Tc)
        p = np.geomspace(spec["p"][0], spec["p"][1], spec["p"][2])

        pc = CP.PropsSI("pcrit", cpn)
        Tt = CP.PropsSI("Ttriple", cpn)
        M = CP.PropsSI("M", cpn)
        acentric = CP.PropsSI("acentric", cpn)
        rhoc = CP.PropsSI("rhocrit", cpn)

        # Saturation curve on the pressure grid (only below the critical point).
        sat = []
        for pi in p:
            if pi < pc:
                try:
                    sat.append(CP.PropsSI("T", "P", float(pi), "Q", 0, cpn))
                except Exception:
                    sat.append(float("nan"))
            else:
                sat.append(float("nan"))

        rows, n_valid = [], 0
        for Ti in T:
            for pi in p:
                try:
                    rho = CP.PropsSI("D", "T", float(Ti), "P", float(pi), cpn)
                    cp_ = CP.PropsSI("C", "T", float(Ti), "P", float(pi), cpn)
                    h = CP.PropsSI("H", "T", float(Ti), "P", float(pi), cpn)
                    mu = CP.PropsSI("V", "T", float(Ti), "P", float(pi), cpn)
                    k = CP.PropsSI("L", "T", float(Ti), "P", float(pi), cpn)
                    ok = 1 if all(np.isfinite(v) and v > 0 for v in (rho, cp_, mu, k)) else 0
                except Exception:
                    rho = cp_ = h = mu = k = float("nan")
                    ok = 0
                if ok:
                    n_valid += 1
                rows.append((Ti, pi, rho, cp_, h, mu, k, ok))

        path = os.path.join(args.output_dir, f"{name}.csv")
        with open(path, "w") as fh:
            fh.write(f"# Ignis real-fluid property table: {name}\n")
            fh.write(f"# generated {stamp} by tools/build_coolant_tables.py\n")
            fh.write(f"# source: CoolProp {CP.get_global_param_string('version')}"
                     f" reference EOS for {cpn}\n")
            fh.write("# units: T[K] p[Pa] rho[kg/m3] cp[J/(kg K)] h[J/kg]"
                     " mu[Pa s] k[W/(m K)] valid[0/1]\n")
            fh.write(f"# critical: Tc={Tc:.6f} pc={pc:.6f} rhoc={rhoc:.6f}"
                     f" Ttriple={Tt:.6f} M={M:.9f} acentric={acentric:.6f}\n")
            fh.write("# saturation curve (p[Pa],Tsat[K]); NaN above the critical pressure\n")
            fh.write("#sat " + " ".join(f"{pi:.8g}:{ti:.8g}" for pi, ti in zip(p, sat)) + "\n")
            fh.write(f"#grid nT={len(T)} np={len(p)}\n")
            fh.write("T,p,rho,cp,h,mu,k,valid\n")
            for r in rows:
                fh.write("%.6g,%.8g,%.8g,%.8g,%.10g,%.8g,%.8g,%d\n" % r)
        print(f"wrote {path}: {len(rows)} points ({n_valid} valid), "
              f"T {T[0]:g}-{T[-1]:g} K, p {p[0]/1e6:g}-{p[-1]/1e6:g} MPa, "
              f"{os.path.getsize(path)/1024:.0f} kB")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
