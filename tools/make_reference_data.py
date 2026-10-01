#!/usr/bin/env python3
"""Generate the external reference data that Ignis is verified and validated against.

The outputs are committed under validation/reference/ so the C++ test suite and
the validation report can run from a clean checkout without Cantera, CoolProp or
NASA CEA installed.  Re-run this script to refresh them.

WHAT EACH FILE IS FOR
---------------------
thermo_reference.csv
    Species cp, h and s at a range of temperatures, computed twice:
      * "nasa"  -- Cantera reading the same NASA TM-4513 coefficients Ignis uses.
                   Agreement here is VERIFICATION: two independent
                   implementations of the same polynomial must agree to
                   round-off.
      * "gri30" -- Cantera reading GRI-Mech 3.0, an independent fit of the same
                   species from a different evaluation of the underlying
                   spectroscopic and calorimetric data.  Disagreement here
                   measures the spread between published data sets, not a code
                   error.

transport_reference.csv
    Mixture-averaged viscosity and thermal conductivity from Cantera's
    kinetic-theory transport, for verifying the Chapman-Enskog implementation.

equilibrium_reference.csv
    Adiabatic equilibrium states from Cantera's Gibbs minimiser over mixture
    ratio and pressure, restricted to exactly the species Ignis uses.  Cantera's
    CHEMKIN-derived files default to a 1 atm standard state, whereas NASA
    TM-4513 is referenced to 1 bar; the reference pressure is therefore reset to
    1 bar so the two solve the same problem.  Agreement is VERIFICATION of the
    minimiser against an independent implementation on identical data.

cea_reference.csv
    NASA CEA results (the actual NASA Glenn Fortran program, through the
    rocketcea wrapper) for both propellant combinations.  CEA uses its own 2002
    nine-coefficient thermodynamic database and its own species list, so
    agreement here is VALIDATION against an authoritative independent tool, and
    the residual differences are a real measure of data-set spread.

coolant_reference.csv
    Reference-equation-of-state properties for methane, hydrogen and oxygen from
    CoolProp, for checking the shipped coolant tables and the Peng-Robinson
    implementation.

CITATIONS
---------
Cantera   I. H. Bell, J. Wronski, S. Quoilin, V. Lemort (CoolProp),
          Ind. Eng. Chem. Res. 53, 2498 (2014);
          D. G. Goodwin et al., Cantera, https://cantera.org
CEA       B. J. McBride and S. Gordon, NASA RP-1311 Parts I (1994), II (1996);
          rocketcea, https://rocketcea.readthedocs.io
"""
from __future__ import annotations

import argparse
import datetime as _dt
import os
import sys

# The species Ignis carries, and the name each one has in the upstream files.
IGNIS_TO_NASA = {
    "H2": "H2", "H": "H", "O": "O", "O2": "O2", "OH": "OH", "H2O": "H2O",
    "HO2": "HO2", "H2O2": "H2O2", "CO": "CO", "CO2": "CO2", "C": "C", "CH": "CH",
    "CH2": "CH2", "CH3": "CH3", "CH4": "CH4", "HCO": "HCO",
    "CH2O": "HCHO,formaldehy", "CH2OH": "CH2OH", "CH3O": "CH3O", "CH3OH": "CH3OH",
    "C2H": "C2H", "C2H2": "C2H2,acetylene", "C2H4": "C2H4", "C2H6": "C2H6",
    "CH2CO": "CH2CO,ketene", "C2": "C2", "N2": "N2", "N": "N", "NO": "NO",
    "NO2": "NO2", "N2O": "N2O", "NH": "NH", "NH2": "NH2", "NH3": "NH3",
    "HNO": "HNO", "HCN": "HCN", "CN": "CN", "NCO": "NCO", "HNCO": "HNCO",
    "AR": "Ar",
}
IGNIS_TO_GRI = {k: k for k in IGNIS_TO_NASA}
IGNIS_TO_GRI["AR"] = "AR"
IGNIS_TO_GRI.pop("C2", None)     # not in GRI-Mech 3.0

CHO_SPECIES = ["H2", "H", "O", "O2", "OH", "H2O", "HO2", "H2O2", "CO", "CO2", "C",
               "CH", "CH2", "CH3", "CH4", "HCO", "CH2O", "CH2OH", "CH3O", "CH3OH",
               "C2H", "C2H2", "C2H4", "C2H6", "CH2CO", "C2"]
HO_SPECIES = ["H2", "H", "O", "O2", "OH", "H2O", "HO2", "H2O2"]

CAL = 4.184          # thermochemical calorie, J
# NASA CEA propellant library values, as printed in the CEA input deck.
CEA_PROPELLANT = {
    "LOX":  dict(h_cal=-3102.0,  T=90.18,  atoms={"O": 2}, M=31.998),
    "LCH4": dict(h_cal=-21390.0, T=111.66, atoms={"C": 1, "H": 4}, M=16.043),
    "LH2":  dict(h_cal=-2154.0,  T=20.27,  atoms={"H": 2}, M=2.016),
}


def header(fh, what, sources):
    stamp = _dt.datetime.now(_dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    fh.write(f"# Ignis reference data: {what}\n")
    fh.write(f"# generated {stamp} by tools/make_reference_data.py\n")
    for s in sources:
        fh.write(f"# source: {s}\n")


def write_thermo(out_dir):
    import cantera as ct
    nasa = {s.name: s for s in ct.Species.list_from_file("nasa_gas.yaml")}
    gri = {s.name: s for s in ct.Solution("gri30.yaml").species()}
    R = ct.gas_constant  # J/(kmol K)

    def props(sp, T):
        g = ct.Solution(thermo="ideal-gas", species=[sp])
        g.TP = T, 1e5
        return g.cp_mole / 1000.0, g.enthalpy_mole / 1000.0, g.entropy_mole / 1000.0

    temps = [200.0, 298.15, 500.0, 999.9, 1000.0, 1000.1, 1500.0, 2000.0, 3000.0,
             4000.0, 5000.0, 6000.0]
    path = os.path.join(out_dir, "thermo_reference.csv")
    with open(path, "w") as fh:
        header(fh, "species thermodynamic properties",
               ["Cantera %s reading nasa_gas.yaml (NASA TM-4513) with the standard "
                "state reset to 1 bar" % ct.__version__,
                "Cantera reading gri30.yaml (GRI-Mech 3.0), an independent fit"])
        fh.write("# units: T[K] cp[J/(mol K)] h[J/mol] s[J/(mol K)]\n")
        fh.write("species,T,cp_nasa,h_nasa,s_nasa,cp_gri,h_gri,s_gri\n")
        n = 0
        for ignis_name, nasa_name in IGNIS_TO_NASA.items():
            src = nasa.get(nasa_name)
            if src is None:
                sys.exit(f"missing NASA species {nasa_name}")
            d = src.input_data
            d["thermo"]["reference-pressure"] = "1 bar"
            sp_nasa = ct.Species.from_dict(d)
            gri_name = IGNIS_TO_GRI.get(ignis_name)
            sp_gri = None
            if gri_name and gri_name in gri:
                dg = gri[gri_name].input_data
                dg["thermo"]["reference-pressure"] = "1 bar"
                sp_gri = ct.Species.from_dict(dg)
            for T in temps:
                lo, hi = sp_nasa.thermo.min_temp, sp_nasa.thermo.max_temp
                if T < lo or T > hi:
                    continue
                cp1, h1, s1 = props(sp_nasa, T)
                if sp_gri is not None and sp_gri.thermo.min_temp <= T <= sp_gri.thermo.max_temp:
                    cp2, h2, s2 = props(sp_gri, T)
                else:
                    cp2 = h2 = s2 = float("nan")
                fh.write("%s,%.6f,%.12g,%.12g,%.12g,%.12g,%.12g,%.12g\n"
                         % (ignis_name, T, cp1, h1, s1, cp2, h2, s2))
                n += 1
    print(f"wrote {path}: {n} rows")
    return R


def _cho_solution(names):
    import cantera as ct
    wanted = {IGNIS_TO_NASA[n] for n in names}
    sp = []
    for s in ct.Species.list_from_file("nasa_gas.yaml"):
        if s.name in wanted:
            d = s.input_data
            d["thermo"]["reference-pressure"] = "1 bar"
            sp.append(ct.Species.from_dict(d))
    return ct.Solution(thermo="ideal-gas", species=sp)


def write_equilibrium(out_dir):
    import cantera as ct
    path = os.path.join(out_dir, "equilibrium_reference.csv")
    rows = []
    cases = [("LOX", "LCH4", CHO_SPECIES, [2.4, 2.8, 3.2, 3.4, 3.6, 4.0, 4.4]),
             ("LOX", "LH2", HO_SPECIES, [3.5, 4.5, 5.5, 6.0, 7.0, 8.0])]
    for ox, fu, names, ratios in cases:
        g = _cho_solution(names)
        O, F = CEA_PROPELLANT[ox], CEA_PROPELLANT[fu]
        for p in (1.0e5, 2.0e6, 5.5e6, 1.0e7, 2.0e7):
            for mr in ratios:
                mo, mf = mr / (1 + mr), 1 / (1 + mr)
                # Reactant enthalpy on the absolute NASA scale, per kg.
                h0 = (mo / (O["M"] * 1e-3)) * O["h_cal"] * CAL + \
                     (mf / (F["M"] * 1e-3)) * F["h_cal"] * CAL
                # The unburned mixture cannot hold the liquid-propellant
                # enthalpy at any temperature inside the polynomial range, so
                # burn it first at a nominal temperature and only then impose
                # the reactant enthalpy.  The element vector is what is
                # conserved, so the intermediate state does not matter.
                feed = {"CH4": mf, "O2": mo} if fu == "LCH4" else {"H2": mf, "O2": mo}
                g.TPY = 3000.0, p, feed
                g.equilibrate("TP")
                g.HP = h0, p
                g.equilibrate("HP")
                row = dict(oxidizer=ox, fuel=fu, mixture_ratio=mr, pressure=p,
                           reactant_enthalpy=h0, temperature=g.T,
                           molar_mass=g.mean_molecular_weight * 1e-3,
                           density=g.density, entropy=g.entropy_mass,
                           cp_frozen=g.cp_mass)
                for n in names:
                    cn = IGNIS_TO_NASA[n]
                    row["X_" + n] = g[cn].X[0] if cn in g.species_names else 0.0
                rows.append(row)
    keys = ["oxidizer", "fuel", "mixture_ratio", "pressure", "reactant_enthalpy",
            "temperature", "molar_mass", "density", "entropy", "cp_frozen"]
    species_cols = sorted({k for r in rows for k in r if k.startswith("X_")})
    with open(path, "w") as fh:
        header(fh, "adiabatic equilibrium states",
               ["Cantera %s Gibbs minimiser on NASA TM-4513 data with a 1 bar "
                "standard state, restricted to the Ignis species set" % ct.__version__,
                "reactant enthalpies from the NASA CEA propellant library"])
        fh.write("# units: pressure[Pa] enthalpy[J/kg] temperature[K] molar_mass[kg/mol] "
                 "density[kg/m^3] entropy[J/(kg K)] cp[J/(kg K)] X[-]\n")
        fh.write(",".join(keys + species_cols) + "\n")
        for r in rows:
            fh.write(",".join([str(r[k]) if isinstance(r[k], str) else "%.12g" % r[k]
                               for k in keys] +
                              ["%.12g" % r.get(c, 0.0) for c in species_cols]) + "\n")
    print(f"wrote {path}: {len(rows)} rows")


def write_transport(out_dir):
    import cantera as ct
    g = ct.Solution("gri30.yaml")
    path = os.path.join(out_dir, "transport_reference.csv")
    mixtures = [
        ("H2O", {"H2O": 1.0}),
        ("CO2", {"CO2": 1.0}),
        ("N2", {"N2": 1.0}),
        ("H2", {"H2": 1.0}),
        ("O2", {"O2": 1.0}),
        ("CO", {"CO": 1.0}),
        ("products_ch4", {"H2O": 0.4950, "CO": 0.1781, "CO2": 0.1305, "H2": 0.0801,
                          "OH": 0.0630, "O2": 0.0226, "H": 0.0210, "O": 0.0095}),
        ("products_h2", {"H2O": 0.7300, "H2": 0.1400, "OH": 0.0600, "H": 0.0350,
                         "O": 0.0150, "O2": 0.0200}),
    ]
    with open(path, "w") as fh:
        header(fh, "mixture-averaged transport properties",
               ["Cantera %s kinetic-theory transport with GRI-Mech 3.0 "
                "Lennard-Jones parameters" % ct.__version__])
        fh.write("# units: T[K] p[Pa] viscosity[Pa s] conductivity[W/(m K)] cp[J/(kg K)]\n")
        fh.write("mixture,T,p,viscosity,conductivity,cp_mass,composition\n")
        n = 0
        for name, X in mixtures:
            for T in (300.0, 600.0, 1000.0, 2000.0, 3000.0, 3500.0):
                g.TPX = T, 1e5, X
                comp = " ".join("%s:%.6g" % kv for kv in X.items())
                fh.write("%s,%.4f,%.6g,%.10g,%.10g,%.10g,%s\n"
                         % (name, T, 1e5, g.viscosity, g.thermal_conductivity,
                            g.cp_mass, comp))
                n += 1
    print(f"wrote {path}: {n} rows")


def write_cea(out_dir):
    from rocketcea.cea_obj import CEA_Obj
    path = os.path.join(out_dir, "cea_reference.csv")
    PSIA = 1.0 / 6894.757293168361
    FT_S = 0.3048
    # RP1_NASA is rocketcea's card for CEA's own library RP-1 (CH1.95,
    # h = -5907.672 cal/mol at 298.15 K), which is what Ignis's RP1 entry holds.
    cases = [("LOX", "CH4", "LCH4", [2.4, 2.8, 3.2, 3.4, 3.6, 4.0, 4.4]),
             ("LOX", "LH2", "LH2", [3.5, 4.5, 5.5, 6.0, 7.0, 8.0]),
             ("LOX", "RP1_NASA", "RP1", [2.0, 2.2, 2.4, 2.6, 2.8, 3.0])]
    rows = []
    for ox, cea_fuel, ignis_fuel, ratios in cases:
        for frozen in (0, 1):
            C = CEA_Obj(oxName=ox, fuelName=cea_fuel)
            for pc in (1.0e5, 2.0e6, 5.5e6, 1.0e7, 2.0e7):
                for mr in ratios:
                    for eps in (20.0, 45.0, 60.0):
                        pc_psia = pc * PSIA
                        try:
                            # frozenAtThroat=0 freezes the composition at the
                            # chamber, which is what Ignis's frozen model does.
                            # Leaving it at 1 would compare a chamber-frozen
                            # expansion against a throat-frozen one and the c*
                            # values would differ by definition.
                            isp, cstar, tc = C.get_IvacCstrTc(Pc=pc_psia, MR=mr, eps=eps,
                                                             frozen=frozen, frozenAtThroat=0)
                            mw, gam = C.get_Chamber_MolWt_gamma(Pc=pc_psia, MR=mr, eps=eps)
                            # get_PambCf divides by the ambient pressure, so the
                            # vacuum thrust coefficient is taken from its
                            # definition CF = Isp g0 / c* instead.
                            cf = isp * 9.80665 / (cstar * FT_S)
                            pcpe = C.get_PcOvPe(Pc=pc_psia, MR=mr, eps=eps, frozen=frozen,
                                                frozenAtThroat=0)
                            te = C.get_Temperatures(Pc=pc_psia, MR=mr, eps=eps, frozen=frozen,
                                                    frozenAtThroat=0)
                        except Exception as exc:  # pragma: no cover
                            print(f"  CEA failed for {ox}/{cea_fuel} MR={mr} Pc={pc}: {exc}")
                            continue
                        rows.append(dict(oxidizer=ox, fuel=ignis_fuel,
                                         composition="frozen" if frozen else "equilibrium",
                                         mixture_ratio=mr, pressure=pc, expansion_ratio=eps,
                                         chamber_temperature=tc / 1.8,
                                         molar_mass=mw * 1e-3, gamma=gam,
                                         c_star=cstar * FT_S, isp_vacuum=isp, cf_vacuum=cf,
                                         pc_over_pe=pcpe,
                                         throat_temperature=te[1] / 1.8,
                                         exit_temperature=te[2] / 1.8))
    keys = list(rows[0].keys())
    with open(path, "w") as fh:
        header(fh, "NASA CEA rocket performance",
               ["NASA Glenn CEA (McBride & Gordon, NASA RP-1311) through rocketcea",
                "liquid propellants at their normal boiling points, infinite-area combustor"])
        fh.write("# units: pressure[Pa] temperature[K] molar_mass[kg/mol] c_star[m/s] isp[s]\n")
        fh.write(",".join(keys) + "\n")
        for r in rows:
            fh.write(",".join(str(r[k]) if isinstance(r[k], str) else "%.10g" % r[k]
                              for k in keys) + "\n")
    print(f"wrote {path}: {len(rows)} rows")


def write_coolant(out_dir):
    import CoolProp.CoolProp as CP
    path = os.path.join(out_dir, "coolant_reference.csv")
    # The last points of each cryogen are turbopump states: a pump inlet at a
    # few bar on the cold liquid, and discharges up to the 41 MPa of a
    # hydrogen pump and the 48 MPa of an oxygen preburner boost pump, where
    # tests/validation/test_cycle_validation.cpp integrates the pump work.
    points = {
        "methane": [(120.0, 1.0e7), (150.0, 1.5e7), (200.0, 1.5e7), (300.0, 1.5e7),
                    (400.0, 1.0e7), (600.0, 1.5e7), (800.0, 2.0e7), (500.0, 5.0e6),
                    (112.0, 3.0e5), (115.0, 1.0e7)],
        "hydrogen": [(30.0, 1.2e7), (60.0, 1.2e7), (100.0, 1.2e7), (200.0, 1.0e7),
                     (400.0, 1.0e7), (700.0, 1.5e7),
                     (21.0, 2.5e5), (21.0, 2.0e6), (25.0, 1.0e7), (32.0, 3.0e7), (45.0, 4.1e7)],
        "oxygen": [(100.0, 1.0e7), (150.0, 1.0e7), (300.0, 1.0e7), (500.0, 1.5e7),
                   (91.0, 7.0e5), (92.0, 2.8e7), (95.0, 4.8e7)],
        # The RP-1 surrogate: a cold inlet, the jacket's working range, and
        # near the critical temperature at a jacket pressure.
        "dodecane": [(300.0, 1.0e7), (400.0, 1.0e7), (500.0, 1.5e7), (600.0, 1.0e7),
                     (650.0, 2.0e7)],
    }
    names = {"methane": "Methane", "hydrogen": "Hydrogen", "oxygen": "Oxygen",
             "dodecane": "n-Dodecane"}
    with open(path, "w") as fh:
        header(fh, "coolant real-fluid properties",
               ["CoolProp %s reference equations of state (Setzmann & Wagner 1991 for "
                "methane, Leachman et al. 2009 for hydrogen, Schmidt & Wagner 1985 for "
                "oxygen, Lemmon & Huber 2004 for n-dodecane)"
                % CP.get_global_param_string("version")])
        fh.write("# units: T[K] p[Pa] rho[kg/m^3] cp[J/(kg K)] mu[Pa s] k[W/(m K)]\n")
        fh.write("fluid,T,p,rho,cp,mu,k,Tcrit,pcrit,acentric,molar_mass\n")
        n = 0
        for fluid, pts in points.items():
            cn = names[fluid]
            for T, p in pts:
                fh.write("%s,%.4f,%.6g,%.10g,%.10g,%.10g,%.10g,%.6f,%.6f,%.6f,%.9f\n"
                         % (fluid, T, p, CP.PropsSI("D", "T", T, "P", p, cn),
                            CP.PropsSI("C", "T", T, "P", p, cn),
                            CP.PropsSI("V", "T", T, "P", p, cn),
                            CP.PropsSI("L", "T", T, "P", p, cn),
                            CP.PropsSI("Tcrit", cn), CP.PropsSI("pcrit", cn),
                            CP.PropsSI("acentric", cn), CP.PropsSI("M", cn)))
                n += 1
    print(f"wrote {path}: {n} rows")


def _plain(o):
    """Cantera's input_data mappings as plain Python, for yaml.safe_dump."""
    if hasattr(o, "items"):
        return {str(k): _plain(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_plain(v) for v in o]
    if hasattr(o, "item"):          # numpy scalars
        return o.item()
    return o


def _kinetic_gas(names):
    """Cantera gas on Ignis's own species data (NASA TM-4513, 1 bar standard
    state, Ignis's names) with the GRI-Mech 3.0 reactions of Ignis's nozzle
    mechanism -- identical data, an independent implementation."""
    import cantera as ct
    import yaml
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from build_kinetics import SPECIES as REACTING
    by_nasa = {IGNIS_TO_NASA[n]: n for n in names}
    species = []
    for s in ct.Species.list_from_file("nasa_gas.yaml"):
        if s.name in by_nasa:
            d = _plain(s.input_data)
            d["thermo"]["reference-pressure"] = "1 bar"
            d["name"] = by_nasa[s.name]
            species.append(d)
    present = {d["name"] for d in species}
    gri = ct.Solution("gri30.yaml")
    reactions, ids = [], []
    for i, r in enumerate(gri.reactions()):
        sp = set(r.reactants) | set(r.products)
        if not (sp <= set(REACTING) and sp <= present):
            continue
        tb = r.third_body
        if tb is not None and tb.name != "M" and tb.name not in present:
            continue
        reactions.append(_plain(r.input_data))
        ids.append(i + 1)
    doc = {"phases": [{"name": "gas", "thermo": "ideal-gas", "elements": ["O", "H", "C"],
                       "species": [d["name"] for d in species], "kinetics": "gas",
                       "reactions": "all", "skip-undeclared-third-bodies": True}],
           "species": species, "reactions": reactions}
    return ct.Solution(yaml=yaml.safe_dump(doc)), ids


def write_kinetics(out_dir):
    """Forward rate constants, and a finite-rate nozzle march through a cone,
    both from Cantera on identical species data."""
    import cantera as ct
    import numpy as np
    from scipy.integrate import solve_ivp
    gas, ids = _kinetic_gas(CHO_SPECIES)
    names = gas.species_names
    src = "Cantera %s on NASA TM-4513 data (1 bar) with GRI-Mech 3.0's rates" % ct.__version__

    # --- the burnt gas: LOX/CH4 at O/F 3.4 and 5.5 MPa, as the equilibrium
    # reference makes it, then a shifting-equilibrium expansion. -----------
    O, F = CEA_PROPELLANT["LOX"], CEA_PROPELLANT["LCH4"]
    mr, pc = 3.4, 5.5e6
    mo, mf = mr / (1 + mr), 1 / (1 + mr)
    h0 = (mo / (O["M"] * 1e-3)) * O["h_cal"] * CAL + (mf / (F["M"] * 1e-3)) * F["h_cal"] * CAL
    gas.TPY = 3000.0, pc, {"CH4": mf, "O2": mo}
    gas.equilibrate("TP")
    gas.HP = h0, pc
    gas.equilibrate("HP")
    s0 = gas.entropy_mass

    def shifting(p):
        gas.SP = s0, p
        gas.equilibrate("SP")
        u = np.sqrt(max(2.0 * (h0 - gas.enthalpy_mass), 0.0))
        return u, gas.density * u

    # Rate constants at four temperatures, on this gas's equilibrium mixture.
    rates_path = os.path.join(out_dir, "kinetics_rates_reference.csv")
    with open(rates_path, "w") as fh:
        header(fh, "forward rate constants", [src])
        fh.write("# units: T[K] p[Pa] M[mol/m^3] kf[(m^3/mol)^(n-1)/s, the third body excluded]\n")
        fh.write("gri,T,p,M,kf\n")
        n = 0
        for T in (800.0, 1500.0, 2500.0, 3500.0):
            gas.TP = T, 2.0e6
            gas.equilibrate("TP")
            kf = gas.forward_rate_constants
            M = gas.third_body_concentrations
            for j, r in enumerate(gas.reactions()):
                order = sum(r.reactants.values())
                extra = order if r.reaction_type.startswith("three-body") else order - 1
                m = M[j] * 1e3 if np.isfinite(M[j]) else 0.0
                fh.write("%d,%.6f,%.6g,%.12g,%.12g\n" % (ids[j], T, 2.0e6, m, kf[j] * 1e-3 ** extra))
                n += 1
    print(f"wrote {rates_path}: {n} rows")

    # Throat: the largest mass flux; start: frozen Mach 1.10 downstream.
    from scipy.optimize import brentq, minimize_scalar
    res = minimize_scalar(lambda p: -shifting(p)[1], bounds=(0.4 * pc, 0.75 * pc),
                          method="bounded", options={"xatol": 1e-3})
    p_star = res.x
    g_star = shifting(p_star)[1]

    def frozen_mach(p):
        u, _ = shifting(p)
        cp, cv = gas.cp_mass, gas.cv_mass
        return u / np.sqrt(cp / cv * gas.P / gas.density)

    p_start = brentq(lambda p: frozen_mach(p) - 1.10, 0.05 * pc, p_star, xtol=1e-6, rtol=1e-14)
    u0, g0 = shifting(p_start)
    eps_start = g_star / g0
    T0, X0, Y0 = gas.T, gas.X.copy(), gas.Y.copy()
    rt = 0.070
    rs = rt * np.sqrt(eps_start)
    tan = np.tan(np.radians(15.0))
    x_end = (rt * np.sqrt(20.0) - rs) / tan
    area = lambda x: np.pi * (rs + x * tan) ** 2
    darea = lambda x: 2.0 * np.pi * (rs + x * tan) * tan
    mdot = gas.density * u0 * area(0.0)
    W = gas.molecular_weights * 1e-3                  # kg/mol
    react = [names.index(k) for k in names if k in set(["H2", "H", "O", "O2", "OH", "H2O", "HO2",
                                                         "H2O2", "CO", "CO2", "HCO"])]
    n_fixed = Y0 / W                                  # mol/kg, every species
    Ru = ct.gas_constant * 1e-3                       # J/(mol K)

    def unpack(y):
        n = n_fixed.copy()
        n[react] = y[:-1]
        return n, y[-1]

    last_T = [T0]

    def temperature(Y, h):
        # Newton on h(T) at fixed composition, to round-off: Cantera's own
        # HP setter stops near 1e-8, and that noise would stall a 1e-10 march.
        T = last_T[0]
        for _ in range(50):
            gas.TPY = T, 1e5, Y
            dT = (gas.enthalpy_mass - h) / gas.cp_mass
            T -= dT
            if abs(dT) < 1e-13 * T:
                break
        last_T[0] = T
        return T

    def rhs(x, y):
        n, u = unpack(y)
        Y = n * W
        T = temperature(Y / Y.sum(), h0 - 0.5 * u * u)
        rho = mdot / (u * area(x))
        gas.TDY = T, rho, Y / Y.sum()
        wdot = gas.net_production_rates * 1e3          # mol/(m^3 s)
        dn = wdot[react] / (rho * u)
        H = gas.partial_molar_enthalpies[react] * 1e-3 # J/mol
        N = n.sum()
        cp = gas.cp_mass
        R = Ru * N
        psi = dn.sum() / N - (H * dn).sum() / (cp * T)
        m2 = u * u / (cp / (cp - R) * R * T)
        du = u * (-darea(x) / area(x) + psi) / (1.0 - m2)
        return np.append(dn, du)

    y0 = np.append(n_fixed[react], u0)
    stations = [0.0, 0.002, 0.01, 0.05, 0.2, 0.5, 1.0]
    sol = solve_ivp(rhs, (0.0, x_end), y0, method="Radau", rtol=1e-9,
                    atol=np.append(np.full(len(react), 1e-16), 1e-8),
                    t_eval=[f * x_end for f in stations])
    if not sol.success:
        sys.exit("kinetic march failed: " + sol.message)

    path = os.path.join(out_dir, "kinetics_nozzle_reference.csv")
    with open(path, "w") as fh:
        header(fh, "finite-rate nozzle march through a 15-degree cone",
               [src, "LOX/CH4 at O/F 3.4 and 5.5 MPa with CEA's propellant enthalpies, "
                "shifting equilibrium to frozen Mach 1.10, then the march by scipy's Radau "
                "integrator at rtol 1e-9"])
        fh.write("# geometry: r(x) = r_start + x tan(15 deg) from x = 0 to x_end; throat radius "
                 "%.4f m\n" % rt)
        fh.write("# units: x[m] r[m] T[K] p[Pa] u[m/s] X[-]; row 'start' is the march's input\n")
        fh.write("kind,x,r,T,p,u," + ",".join("X_" + k for k in names) + "\n")
        fh.write("start,0,%.12g,%.12g,%.12g,%.12g," % (rs, T0, p_start, u0) +
                 ",".join("%.12g" % v for v in X0) + "\n")
        fh.write("end,%.12g,%.12g,0,0,0," % (x_end, rs + x_end * tan) +
                 ",".join("0" for _ in names) + "\n")
        for x, y in zip(sol.t, sol.y.T):
            n, u = unpack(y)
            Y = n * W
            T = temperature(Y / Y.sum(), h0 - 0.5 * u * u)
            rho = mdot / (u * area(x))
            gas.TDY = T, rho, Y / Y.sum()
            fh.write("station,%.12g,%.12g,%.12g,%.12g,%.12g," % (x, rs + x * tan, T, gas.P, u) +
                     ",".join("%.12g" % v for v in gas.X) + "\n")
    print(f"wrote {path}: start at A/A* {eps_start:.5f}, {len(sol.t)} stations, "
          f"exit {sol.y[-1, -1]:.3f} m/s")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-o", "--output-dir", default="validation/reference")
    ap.add_argument("--skip", nargs="*", default=[],
                    help="sections to skip: thermo transport equilibrium cea coolant kinetics")
    args = ap.parse_args(argv)
    os.makedirs(args.output_dir, exist_ok=True)
    for name, fn in (("thermo", write_thermo), ("transport", write_transport),
                     ("equilibrium", write_equilibrium), ("cea", write_cea),
                     ("coolant", write_coolant), ("kinetics", write_kinetics)):
        if name in args.skip:
            print(f"skipping {name}")
            continue
        try:
            fn(args.output_dir)
        except ImportError as exc:
            print(f"cannot generate {name}: {exc}", file=sys.stderr)
            return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
