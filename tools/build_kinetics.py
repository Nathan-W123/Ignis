#!/usr/bin/env python3
"""Extract Ignis's nozzle-recombination mechanism from GRI-Mech 3.0.

WHAT AND WHY
------------
In a rocket nozzle the combustion is over; what is left is recombination.
H, O and OH recombine into H2, O2 and H2O, and CO is oxidised to CO2, as the
gas cools.  Fully shifting equilibrium assumes those reactions are infinitely
fast and frozen flow assumes they stop at the throat; the real gas does
neither.  Ignis integrates the recombination at finite rate along the
divergent nozzle (ignis/kinetics/NozzleKinetics.hpp), and this script supplies
the rates.

SOURCE
------
GRI-Mech 3.0: G. P. Smith, D. M. Golden, M. Frenklach, N. W. Moriarty,
B. Eiteneer, M. Goldenberg, C. T. Bowman, R. K. Hanson, S. Song,
W. C. Gardiner, Jr., V. V. Lissianski and Z. Qin,
http://combustion.berkeley.edu/gri-mech/ -- through the copy Cantera ships
(gri30.yaml).  Only the forward rate parameters are taken.  Ignis computes
each reverse rate from its own species thermodynamics (detailed balance), so
the mechanism and the equilibrium solver can never disagree about where
equilibrium is.

THE SUBSET
----------
Every GRI-Mech 3.0 reaction whose reactants and products all lie in
{H, H2, O, O2, OH, H2O, HO2, H2O2, CO, CO2, HCO}: the hydrogen-oxygen system,
CO oxidation, and the formyl radical through which CO and H exchange.  Those
are the species of a hydrogen or hydrocarbon rocket's exhaust above the ppm
level.  GRI's third-body efficiencies are kept whole, and so are its
reactions with an explicit N2 or Ar collider: Ignis's reader skips any
species the run's database does not carry, so with carbon-hydrogen-oxygen
propellants those reactions simply have no collider.  Hydrocarbon species
present in a fuel-rich chamber at trace level are carried through the nozzle
unreacted.

UNITS
-----
Cantera stores A in kmol-based SI units.  The file written here is mol-based
SI: A in (m^3/mol)^(n-1)/s for a reaction of overall order n (counting the
third body), b dimensionless, Ea in J/mol.

USAGE
-----
    python3 tools/build_kinetics.py -o data/kinetics/gri30_nozzle.yaml
"""
from __future__ import annotations

import argparse
import datetime as _dt
import sys

SPECIES = ["H2", "H", "O", "O2", "OH", "H2O", "HO2", "H2O2", "CO", "CO2", "HCO"]


def _cantera():
    try:
        import cantera as ct
    except ImportError:  # pragma: no cover
        sys.exit("Cantera is required to regenerate the mechanism: pip install cantera")
    return ct


def arrhenius(rate, order):
    """Cantera Arrhenius (kmol units) -> mol-based SI dict."""
    return {"A": rate.pre_exponential_factor * 1e-3 ** (order - 1),
            "b": rate.temperature_exponent,
            "Ea": rate.activation_energy * 1e-3}


def fmt(x):
    return "%.10g" % x


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-o", "--output", default="data/kinetics/gri30_nozzle.yaml")
    args = ap.parse_args(argv)
    ct = _cantera()
    gas = ct.Solution("gri30.yaml")
    keep = set(SPECIES)
    out = []
    for i, r in enumerate(gas.reactions()):
        species = set(r.reactants) | set(r.products)
        if not species <= keep:
            continue
        tb = r.third_body
        order = sum(r.reactants.values())
        entry = {"gri": i + 1, "equation": r.equation, "reactants": dict(r.reactants),
                 "products": dict(r.products), "duplicate": bool(r.duplicate)}
        if not r.reversible:
            sys.exit(f"reaction {r.equation} is irreversible; the reader expects reversible ones")
        rtype = r.reaction_type
        if rtype == "Arrhenius":
            entry["type"] = "elementary"
            entry["rate"] = arrhenius(r.rate, order)
        elif rtype == "three-body-Arrhenius":
            entry["type"] = "three-body"
            entry["rate"] = arrhenius(r.rate, order + 1)
        elif rtype.startswith("falloff"):
            entry["type"] = "falloff"
            entry["high"] = arrhenius(r.rate.high_rate, order)
            entry["low"] = arrhenius(r.rate.low_rate, order + 1)
            coeffs = list(r.rate.falloff_coeffs)
            if rtype == "falloff-Troe":
                entry["troe"] = coeffs
            elif rtype != "falloff-Lindemann":
                sys.exit(f"unsupported falloff form {rtype} in {r.equation}")
        else:
            sys.exit(f"unsupported reaction type {rtype} in {r.equation}")
        if tb is not None:
            if tb.name == "M":
                entry["default_efficiency"] = tb.default_efficiency
                entry["efficiencies"] = dict(tb.efficiencies)
            else:   # an explicit collider, e.g. H + O2 + O2 <=> HO2 + O2
                entry["default_efficiency"] = 0.0
                entry["efficiencies"] = {tb.name: 1.0}
        out.append(entry)

    stamp = _dt.datetime.now(_dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    with open(args.output, "w") as fh:
        fh.write("# Ignis nozzle-recombination mechanism: the H/O/CO/HCO subset of GRI-Mech 3.0\n")
        fh.write(f"# generated {stamp} by tools/build_kinetics.py from Cantera "
                 f"{ct.__version__}'s gri30.yaml\n")
        fh.write("# source: G. P. Smith et al., GRI-Mech 3.0, "
                 "http://combustion.berkeley.edu/gri-mech/\n")
        fh.write("# forward rates only; reverse rates follow from Ignis's own species\n"
                 "# thermodynamics by detailed balance\n")
        fh.write("# units: A in (m^3/mol)^(n-1)/s for overall order n (third body counted),"
                 " b dimensionless, Ea in J/mol\n")
        fh.write("species: [" + ", ".join(SPECIES) + "]\n")
        fh.write("reactions:\n")
        for e in out:
            fh.write(f"- equation: \"{e['equation']}\"\n")
            fh.write(f"  gri: {e['gri']}\n")
            fh.write(f"  type: {e['type']}\n")
            fh.write("  reactants: {" + ", ".join(f"{k}: {fmt(v)}" for k, v in e["reactants"].items()) + "}\n")
            fh.write("  products: {" + ", ".join(f"{k}: {fmt(v)}" for k, v in e["products"].items()) + "}\n")
            for key in ("rate", "high", "low"):
                if key in e:
                    r = e[key]
                    fh.write(f"  {key}: {{A: {fmt(r['A'])}, b: {fmt(r['b'])}, Ea: {fmt(r['Ea'])}}}\n")
            if "troe" in e:
                fh.write("  troe: [" + ", ".join(fmt(v) for v in e["troe"]) + "]\n")
            if "default_efficiency" in e:
                fh.write(f"  default_efficiency: {fmt(e['default_efficiency'])}\n")
                fh.write("  efficiencies: {" + ", ".join(f"{k}: {fmt(v)}" for k, v in e["efficiencies"].items()) + "}\n")
            if e["duplicate"]:
                fh.write("  duplicate: true\n")
    print(f"wrote {args.output}: {len(out)} reactions among {len(SPECIES)} species")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
