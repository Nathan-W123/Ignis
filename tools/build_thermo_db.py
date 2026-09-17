#!/usr/bin/env python3
"""Build the Ignis NASA-7 thermodynamic + transport species database.

Sources
-------
1. Thermodynamic polynomials (NASA 7-coefficient form, 200-6000 K):
     Cantera's ``nasa_gas.yaml``, which is a machine-readable transcription of

       B. J. McBride, S. Gordon and M. A. Reno,
       "Coefficients for Calculating Thermodynamic and Transport Properties of
        Individual Species", NASA TM-4513, October 1993.
       https://ntrs.nasa.gov/citations/19940013151

     Upstream file:
       https://raw.githubusercontent.com/Cantera/cantera/main/data/nasa_gas.yaml

2. Lennard-Jones transport parameters:
     GRI-Mech 3.0 transport database (``gri30_tran.dat``), as transcribed in
     Cantera's ``gri30.yaml``.

       G. P. Smith, D. M. Golden, M. Frenklach, N. W. Moriarty, B. Eiteneer,
       M. Goldenberg, C. T. Bowman, R. K. Hanson, S. Song, W. C. Gardiner Jr.,
       V. V. Lissianski and Z. Qin, "GRI-Mech 3.0",
       http://combustion.berkeley.edu/gri-mech/

     Upstream file:
       https://raw.githubusercontent.com/Cantera/cantera/main/data/gri30.yaml

No coefficient is modified: values are copied verbatim.  The only processing
performed here is (a) selection of the species subset relevant to LOX/CH4 and
LOX/H2 propulsion, (b) renaming of NASA compound names to plain chemical
formulas, and (c) emission in the Ignis database schema.

Usage
-----
    python3 tools/build_thermo_db.py --nasa nasa_gas.yaml --gri gri30.yaml \
        -o data/thermo/ignis_nasa7.yaml
"""
from __future__ import annotations

import argparse
import datetime as _dt
import hashlib
import sys

import yaml


class PlainLoader(yaml.SafeLoader):
    """SafeLoader with YAML-1.1 bool/null resolution disabled.

    Required because species names such as ``NO``, ``N`` and ``ON`` are
    otherwise silently converted to booleans / nulls.
    """


for _ch in list(PlainLoader.yaml_implicit_resolvers):
    PlainLoader.yaml_implicit_resolvers[_ch] = [
        (t, r)
        for (t, r) in PlainLoader.yaml_implicit_resolvers[_ch]
        if t not in ("tag:yaml.org,2002:bool", "tag:yaml.org,2002:null")
    ]


# (NASA database name, Ignis name, GRI-Mech transport name or None)
SPECIES = [
    # --- H/O system -------------------------------------------------------
    ("H2", "H2", "H2"),
    ("H", "H", "H"),
    ("O", "O", "O"),
    ("O2", "O2", "O2"),
    ("OH", "OH", "OH"),
    ("H2O", "H2O", "H2O"),
    ("HO2", "HO2", "HO2"),
    ("H2O2", "H2O2", "H2O2"),
    # --- C/H/O system -----------------------------------------------------
    ("CO", "CO", "CO"),
    ("CO2", "CO2", "CO2"),
    ("C", "C", "C"),
    ("CH", "CH", "CH"),
    ("CH2", "CH2", "CH2"),
    ("CH3", "CH3", "CH3"),
    ("CH4", "CH4", "CH4"),
    ("HCO", "HCO", "HCO"),
    ("HCHO,formaldehy", "CH2O", "CH2O"),
    ("CH2OH", "CH2OH", "CH2OH"),
    ("CH3O", "CH3O", "CH3O"),
    ("CH3OH", "CH3OH", "CH3OH"),
    ("C2H", "C2H", "C2H"),
    ("C2H2,acetylene", "C2H2", "C2H2"),
    ("C2H4", "C2H4", "C2H4"),
    ("C2H6", "C2H6", "C2H6"),
    ("CH2CO,ketene", "CH2CO", "CH2CO"),
    ("C2", "C2", None),
    # --- N system ---------------------------------------------------------
    ("N2", "N2", "N2"),
    ("N", "N", "N"),
    ("NO", "NO", "NO"),
    ("NO2", "NO2", "NO2"),
    ("N2O", "N2O", "N2O"),
    ("NH", "NH", "NH"),
    ("NH2", "NH2", "NH2"),
    ("NH3", "NH3", "NH3"),
    ("HNO", "HNO", "HNO"),
    ("HCN", "HCN", "HCN"),
    ("CN", "CN", "CN"),
    ("NCO", "NCO", "NCO"),
    ("HNCO", "HNCO", "HNCO"),
    # --- inert ------------------------------------------------------------
    ("Ar", "AR", "AR"),
]

# Standard atomic weights, IUPAC 2021 (Prohaska et al., Pure Appl. Chem. 94,
# 573-600, 2022).  Conventional single values for the CHON/Ar elements.
ATOMIC_WEIGHT = {
    "C": 12.011,
    "H": 1.008,
    "O": 15.999,
    "N": 14.007,
    "Ar": 39.95,
}


def sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 16), b""):
            h.update(chunk)
    return h.hexdigest()


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--nasa", required=True, help="path to nasa_gas.yaml")
    ap.add_argument("--gri", required=True, help="path to gri30.yaml")
    ap.add_argument("-o", "--output", required=True)
    args = ap.parse_args(argv)

    nasa = yaml.load(open(args.nasa), Loader=PlainLoader)
    gri = yaml.load(open(args.gri), Loader=PlainLoader)

    nasa_by_name = {s["name"]: s for s in nasa["species"]}
    gri_by_name = {s["name"]: s for s in gri["species"]}

    out_species = []
    missing = []
    for nasa_name, ignis_name, tran_name in SPECIES:
        src = nasa_by_name.get(nasa_name)
        if src is None:
            missing.append(nasa_name)
            continue
        th = src["thermo"]
        if th["model"] != "NASA7":
            raise SystemExit(f"{nasa_name}: unexpected thermo model {th['model']}")
        tr = list(th["temperature-ranges"])
        data = [list(map(float, row)) for row in th["data"]]
        if len(tr) == 2:                      # single-range species (e.g. none here)
            data = [data[0], data[0]]
            tr = [tr[0], 0.5 * (tr[0] + tr[1]), tr[1]]
        if len(tr) != 3 or len(data) != 2:
            raise SystemExit(f"{nasa_name}: unsupported {len(tr)} temperature ranges")

        comp = {k: int(v) for k, v in src["composition"].items()}
        mw = sum(ATOMIC_WEIGHT[e] * n for e, n in comp.items())

        rec = {
            "name": ignis_name,
            "composition": comp,
            "molar-mass": round(mw, 6),          # kg/kmol == g/mol
            "temperature-ranges": [float(tr[0]), float(tr[1]), float(tr[2])],
            "low-coeffs": data[0],
            "high-coeffs": data[1],
            "source": f"NASA TM-4513 via Cantera nasa_gas.yaml (note: {th.get('note','')})",
        }
        if tran_name is not None and tran_name in gri_by_name:
            gt = gri_by_name[tran_name].get("transport")
            if gt is not None:
                rec["transport"] = {
                    "geometry": gt["geometry"],
                    "well-depth": float(gt["well-depth"]),      # K  (eps/kB)
                    "diameter": float(gt["diameter"]),          # Angstrom
                    "dipole": float(gt.get("dipole", 0.0)),     # Debye
                    "polarizability": float(gt.get("polarizability", 0.0)),  # A^3
                    "rotational-relaxation": float(gt.get("rotational-relaxation", 0.0)),
                    "source": "GRI-Mech 3.0 transport database",
                }
        out_species.append(rec)

    if missing:
        print(f"ERROR: species not found in NASA database: {missing}", file=sys.stderr)
        return 1

    header = {
        "description": (
            "Ignis species database: NASA 7-coefficient thermodynamic polynomials "
            "plus Lennard-Jones transport parameters.\n"
            "Generated by tools/build_thermo_db.py -- do not edit by hand.\n"
            "Coefficients are verbatim copies of the cited sources."
        ),
        "generated": _dt.datetime.now(_dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "units": {
            "molar-mass": "kg/kmol",
            "temperature": "K",
            "well-depth": "K (epsilon/k_B)",
            "diameter": "Angstrom",
            "dipole": "Debye",
            "polarizability": "Angstrom^3",
            "coefficients": "dimensionless NASA-7 (cp/R form)",
        },
        "sources": [
            {
                "name": "NASA TM-4513 (McBride, Gordon & Reno, 1993)",
                "url": "https://ntrs.nasa.gov/citations/19940013151",
                "file": "nasa_gas.yaml",
                "sha256": sha256(args.nasa),
            },
            {
                "name": "GRI-Mech 3.0 transport database (Smith et al.)",
                "url": "http://combustion.berkeley.edu/gri-mech/",
                "file": "gri30.yaml",
                "sha256": sha256(args.gri),
            },
            {
                "name": "IUPAC 2021 standard atomic weights (conventional values)",
                "url": "https://doi.org/10.1515/pac-2019-0603",
            },
        ],
        "atomic-weights": ATOMIC_WEIGHT,
    }

    with open(args.output, "w") as fh:
        fh.write("# " + "\n# ".join(header.pop("description").splitlines()) + "\n")
        yaml.dump(header, fh, sort_keys=False, default_flow_style=None, width=100)
        fh.write("\nspecies:\n")
        for rec in out_species:
            fh.write("- name: \"%s\"\n" % rec["name"])
            fh.write("  composition: {%s}\n" %
                     ", ".join(f"{k}: {v}" for k, v in rec["composition"].items()))
            fh.write("  molar-mass: %.6f\n" % rec["molar-mass"])
            fh.write("  temperature-ranges: [%.1f, %.1f, %.1f]\n" % tuple(rec["temperature-ranges"]))
            fh.write("  low-coeffs: [%s]\n" % ", ".join("%.9g" % c for c in rec["low-coeffs"]))
            fh.write("  high-coeffs: [%s]\n" % ", ".join("%.9g" % c for c in rec["high-coeffs"]))
            fh.write("  source: \"%s\"\n" % rec["source"].replace('"', "'"))
            if "transport" in rec:
                t = rec["transport"]
                fh.write("  transport: {geometry: %s, well-depth: %g, diameter: %g, "
                         "dipole: %g, polarizability: %g, rotational-relaxation: %g}\n"
                         % (t["geometry"], t["well-depth"], t["diameter"], t["dipole"],
                            t["polarizability"], t["rotational-relaxation"]))
        fh.write("\n")

    n_tr = sum(1 for r in out_species if "transport" in r)
    print(f"wrote {args.output}: {len(out_species)} species "
          f"({n_tr} with transport data)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
