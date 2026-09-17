"""Every Ignis figure.

Each function takes the paths of the result files it needs and writes one PNG.
Nothing here recomputes physics: the figures show exactly what the C++ solver
wrote, so a figure and the numbers in the report can never disagree.
"""
from __future__ import annotations

import os

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection

from . import style
from .data import Table, read_json, read_table, require

MPA = 1e-6
KN = 1e-3
MM = 1e3


def _save(fig, out: str) -> str:
    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)
    fig.savefig(out)
    plt.close(fig)
    return out


# ---------------------------------------------------------------------------
# Chamber thermochemistry
# ---------------------------------------------------------------------------

def mixture_ratio_composition(sweep_csv: str, equilibrium_csv: str, out: str,
                              max_species: int = 8) -> str:
    """Equilibrium species mole fractions against mixture ratio.

    The sweep table carries only scalar metrics, so the composition is rebuilt
    by reading the per-mixture-ratio equilibrium exports written alongside it.
    """
    t = read_table(require(sweep_csv))
    comp = read_table(require(equilibrium_csv))
    # Rank species by their peak mole fraction and keep the leading ones; the
    # tail is summed into a single "other" series rather than given new hues.
    species_cols = [c for c in comp.columns if c.startswith("X_")]
    peak = {c: comp[c].max() for c in species_cols}
    ranked = sorted(peak, key=peak.get, reverse=True)
    keep, tail = ranked[:max_species], ranked[max_species:]

    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    mr = comp["mixture_ratio"]
    for i, col in enumerate(keep):
        name = col[2:]
        # Eight series is past the direct-label budget: on a log axis their peaks
        # sit on top of each other, so identity is carried by the legend alone.
        ax.plot(mr, comp[col], color=style.CATEGORICAL[i], label=name)
    if tail:
        ax.plot(mr, comp[tail].sum(axis=1), color=style.INK_MUTED, linewidth=1.4,
                label=f"other ({len(tail)} species)")
    ax.set_yscale("log")
    ax.set_ylim(1e-5, 1.2)
    ax.set_xlabel("mixture ratio O/F [-]")
    ax.set_ylabel("mole fraction [-]")
    ax.set_title("Equilibrium chamber composition against mixture ratio")
    ax.legend(ncols=3, loc="lower center", bbox_to_anchor=(0.5, -0.40))
    style.caption(fig, "LOX/CH4, shifting equilibrium at the design chamber pressure. "
                       "Species ranked by peak mole fraction; the remainder summed.")
    return _save(fig, out)


def mixture_ratio_performance(sweep_csv: str, out: str) -> str:
    """Flame temperature, c* and specific impulse against mixture ratio."""
    t = read_table(require(sweep_csv))
    mr = t["propellants.mixture_ratio"]
    fig, axes = plt.subplots(3, 1, figsize=(6.6, 7.6), sharex=True)

    ax = axes[0]
    ax.plot(mr, t["chamber.temperature"], color=style.CATEGORICAL[0])
    k = int(np.argmax(t["chamber.temperature"].to_numpy()))
    ax.plot([mr.iloc[k]], [t["chamber.temperature"].iloc[k]], "o",
            color=style.CATEGORICAL[0], zorder=5)
    style.annotate(ax, mr.iloc[k], t["chamber.temperature"].iloc[k],
                   f"peak {t['chamber.temperature'].iloc[k]:.0f} K at O/F {mr.iloc[k]:.2f}",
                   dx=0, dy=-12, ha="center", va="top")
    ax.set_ylabel(t.label("chamber.temperature", "flame temperature"))
    ax.set_title("Adiabatic flame temperature")

    ax = axes[1]
    ax.plot(mr, t["chamber.c_star_ideal"], color=style.CATEGORICAL[0], label="ideal")
    ax.plot(mr, t["chamber.c_star"], color=style.CATEGORICAL[1], label="after combustion efficiency")
    k = int(np.argmax(t["chamber.c_star_ideal"].to_numpy()))
    ax.plot([mr.iloc[k]], [t["chamber.c_star_ideal"].iloc[k]], "o",
            color=style.CATEGORICAL[0], zorder=5)
    style.annotate(ax, mr.iloc[k], t["chamber.c_star_ideal"].iloc[k],
                   f"peak {t['chamber.c_star_ideal'].iloc[k]:.0f} m/s at O/F {mr.iloc[k]:.2f}",
                   dx=14, dy=-14, va="top")
    ax.set_ylabel(t.label("chamber.c_star_ideal", "characteristic velocity"))
    ax.set_title("Characteristic velocity")
    ax.legend(loc="lower center", ncols=2)

    ax = axes[2]
    ax.plot(mr, t["performance.isp_vacuum"], color=style.CATEGORICAL[0], label="vacuum")
    ax.plot(mr, t["performance.isp"], color=style.CATEGORICAL[1], label="sea level")
    for col, colour in (("performance.isp_vacuum", style.CATEGORICAL[0]),
                        ("performance.isp", style.CATEGORICAL[1])):
        k = int(np.argmax(t[col].to_numpy()))
        ax.plot([mr.iloc[k]], [t[col].iloc[k]], "o", color=colour, zorder=5)
        style.annotate(ax, mr.iloc[k], t[col].iloc[k],
                       f"{t[col].iloc[k]:.1f} s at O/F {mr.iloc[k]:.2f}",
                       dx=8, dy=3, color=colour)
    ax.set_ylabel(t.label("performance.isp", "specific impulse"))
    ax.set_xlabel("mixture ratio O/F [-]")
    ax.set_title("Specific impulse")
    ax.legend(loc="lower center", ncols=2)

    style.caption(fig, "LOX/CH4, shifting equilibrium. The impulse optimum sits fuel-rich of the "
                       "temperature optimum because a lighter exhaust beats a hotter one.")
    return _save(fig, out)


# ---------------------------------------------------------------------------
# Nozzle
# ---------------------------------------------------------------------------

def nozzle_contour_mach(profile_csv: str, out: str) -> str:
    """Nozzle contour coloured by local Mach number."""
    t = read_table(require(profile_csv))
    x = t["x"].to_numpy() * MM
    r = t["radius"].to_numpy() * MM
    m = t["mach"].to_numpy()

    fig, ax = plt.subplots(figsize=(8.6, 4.2))
    norm = plt.Normalize(m.min(), m.max())

    # The quasi-1D field is a function of x alone, so the interior is filled
    # with a vertical-band gradient: that is exactly what the solver computed,
    # with no interpolation across the radius implied.
    ny = 220
    yy = np.linspace(-r.max(), r.max(), ny)
    field = np.tile(m, (ny, 1))
    inside = np.abs(yy)[:, None] <= r[None, :]
    im = ax.pcolormesh(x, yy, np.ma.masked_where(~inside, field), cmap=style.SEQUENTIAL,
                       norm=norm, shading="nearest", rasterized=True)
    for sign in (1.0, -1.0):
        ax.plot(x, sign * r, color=style.INK_PRIMARY, linewidth=1.4, solid_joinstyle="round")

    i_throat = int(np.argmin(r))
    ax.plot([x[i_throat], x[i_throat]], [-r[i_throat], r[i_throat]],
            color=style.SURFACE, linewidth=1.0)
    style.annotate(ax, x[i_throat], r.max() * 0.98, "throat, M = 1", dy=2, ha="center")
    style.annotate(ax, x[-1], r[-1], f"exit  M = {m[-1]:.2f}", dx=-6, ha="right", dy=8)

    ax.set_xlim(x.min(), x.max() * 1.02)
    ax.set_ylim(-r.max() * 1.25, r.max() * 1.25)
    ax.set_aspect("equal")
    ax.set_xlabel("axial position from the injector face [mm]")
    ax.set_ylabel("radius [mm]")
    ax.set_title("Nozzle contour coloured by Mach number")
    ax.grid(False)
    cb = fig.colorbar(im, ax=ax, pad=0.015, fraction=0.035)
    cb.set_label("Mach number [-]", color=style.INK_SECONDARY, fontsize=9)
    cb.outline.set_visible(False)
    style.caption(fig, "Quasi-1D shifting-equilibrium expansion, so the Mach number varies only "
                       "along the axis. The wall is the analytic chamber/fillet/cone/"
                       "throat-arc/bell construction.")
    return _save(fig, out)


def axial_profiles(profile_csv: str, out: str, case: str = "Ignis-M1") -> str:
    """Pressure, temperature, density, velocity and Mach number along the axis."""
    t = read_table(require(profile_csv))
    x = t["x"] * MM
    xt = float(t["x"][t["mach"].sub(1.0).abs().idxmin()] * MM)

    panels = [("pressure", "static pressure [MPa]", MPA),
              ("temperature", "static temperature [K]", 1.0),
              ("density", "density [kg/m^3]", 1.0),
              ("velocity", "velocity [m/s]", 1.0),
              ("mach", "Mach number [-]", 1.0),
              ("gamma_s", "isentropic exponent [-]", 1.0)]
    fig, axes = plt.subplots(3, 2, figsize=(9.6, 7.4), sharex=True)
    for ax, (col, label, scale) in zip(axes.ravel(), panels):
        ax.plot(x, t[col] * scale, color=style.CATEGORICAL[0])
        ax.axvline(xt, color=style.INK_MUTED, linewidth=0.9)
        ax.set_ylabel(label)
        if col == "pressure":
            ax.set_yscale("log")
        if col == "mach":
            ax.axhline(1.0, color=style.GRID, linewidth=0.9)
    axes[0, 0].annotate("throat", (xt, 1.0), xycoords=("data", "axes fraction"),
                        textcoords="offset points", xytext=(4, -12), fontsize=8.5,
                        color=style.INK_MUTED)
    for ax in axes[-1]:
        ax.set_xlabel("axial position [mm]")
    style.suptitle(fig, f"{case}: flow properties along the nozzle axis")
    style.caption(fig, "The vertical rule marks the sonic throat. Every panel comes from the "
                       "same converged station march.")
    return _save(fig, out)


def altitude_performance(altitude_csv: str, out: str) -> str:
    """Thrust, specific impulse and the thrust split against altitude."""
    t = read_table(require(altitude_csv))
    frame = t.frame[np.isfinite(t.frame["altitude"])]
    z = frame["altitude"] / 1000.0

    fig, axes = plt.subplots(3, 1, figsize=(6.8, 8.0), sharex=True)

    ax = axes[0]
    ax.plot(z, frame["thrust"] * KN, color=style.CATEGORICAL[0], label="total")
    ax.plot(z, frame["thrust_momentum"] * KN, color=style.CATEGORICAL[1], label="momentum")
    ax.plot(z, frame["thrust_pressure"] * KN, color=style.CATEGORICAL[2], label="pressure")
    ax.axhline(0.0, color=style.GRID, linewidth=0.9)
    style.annotate(ax, z.iloc[-1], frame["thrust"].iloc[-1] * KN,
                   f"{frame['thrust'].iloc[-1] * KN:.0f} kN", dx=-4, ha="right",
                   color=style.CATEGORICAL[0])
    ax.set_ylabel("thrust [kN]")
    ax.set_title("Thrust and its two contributions")
    ax.legend(ncols=3, loc="lower right")

    ax = axes[1]
    ax.plot(z, frame["isp"], color=style.CATEGORICAL[0])
    ax.set_ylabel("specific impulse [s]")
    ax.set_title("Specific impulse")

    ax = axes[2]
    ax.plot(z, frame["exit_ambient_ratio"], color=style.CATEGORICAL[0])
    ax.axhline(1.0, color=style.INK_MUTED, linewidth=0.9)
    style.annotate(ax, z.iloc[0], 1.0, "ideally expanded", dy=4, color=style.INK_MUTED)
    ax.set_yscale("log")
    ax.set_ylabel("exit / ambient pressure [-]")
    ax.set_xlabel("geometric altitude [km]")
    ax.set_title("Expansion regime")

    # Shade the regimes using the status palette, which is reserved for state.
    regimes = frame["regime"].to_numpy()
    edges = [0] + [i for i in range(1, len(regimes)) if regimes[i] != regimes[i - 1]] + [len(regimes)]
    shades = {"over-expanded": style.STATUS["warning"], "ideally-expanded": style.STATUS["good"],
              "under-expanded": style.CATEGORICAL[0], "subsonic-exit": style.STATUS["critical"]}
    for a, b in zip(edges[:-1], edges[1:]):
        ax.axvspan(z.iloc[a], z.iloc[min(b, len(z) - 1)],
                   color=shades.get(regimes[a], style.GRID), alpha=0.08, linewidth=0)
        if b - a > len(regimes) // 12:
            ax.text(0.5 * (z.iloc[a] + z.iloc[min(b, len(z) - 1)]), ax.get_ylim()[1] * 0.35,
                    regimes[a].replace("-", "\n"), ha="center", fontsize=8,
                    color=style.INK_MUTED)
    style.caption(fig, "U.S. Standard Atmosphere 1976. Separation, where predicted, is an "
                       "empirical criterion and does not alter the inviscid solution.")
    return _save(fig, out)


def expansion_trade(sweep_csv: str, out: str) -> str:
    """Specific impulse against expansion ratio at several altitudes."""
    t = read_table(require(sweep_csv))
    fig, ax = plt.subplots(figsize=(7.0, 4.6))
    alts = sorted(t["performance.altitude"].unique())
    colours = style.ordinal(len(alts))
    for colour, z in zip(colours, alts):
        sub = t.frame[t.frame["performance.altitude"] == z].sort_values("nozzle.expansion_ratio")
        ax.plot(sub["nozzle.expansion_ratio"], sub["performance.isp"], color=colour,
                label=f"{z / 1000:.0f} km")
        k = int(np.argmax(sub["performance.isp"].to_numpy()))
        ax.plot([sub["nozzle.expansion_ratio"].iloc[k]], [sub["performance.isp"].iloc[k]], "o",
                color=colour, zorder=5,
                markeredgecolor=style.SURFACE, markeredgewidth=1.4)
    ax.set_xscale("log")
    ax.set_xlabel("expansion ratio Ae/At [-]")
    ax.set_ylabel("specific impulse [s]")
    ax.set_title("Specific impulse against expansion ratio")
    ax.legend(title="altitude", ncols=4, loc="lower center",
              title_fontproperties={"size": 9})
    ax.set_ylim(bottom=max(0.0, ax.get_ylim()[0]))
    style.annotate(ax, *ax.get_xlim()[:1], 0, "")
    style.caption(fig, "Markers show the optimum expansion ratio at each altitude. The ladder "
                       "uses one hue because altitude is an ordered magnitude, not an identity.")
    return _save(fig, out)


# ---------------------------------------------------------------------------
# Thermal
# ---------------------------------------------------------------------------

def thermal_profiles(thermal_csv: str, out: str, case: str = "Ignis-M1") -> str:
    """Film coefficient, heat flux, wall and coolant conditions along the jacket."""
    t = read_table(require(thermal_csv))
    x = t["x"] * MM
    xq = float(t["x"][t["q_total"].idxmax()] * MM)

    fig, axes = plt.subplots(4, 1, figsize=(7.2, 9.4), sharex=True)

    ax = axes[0]
    ax.plot(x, t["h_gas"] * 1e-3, color=style.CATEGORICAL[0], label="hot gas")
    ax.plot(x, t["h_coolant"] * 1e-3, color=style.CATEGORICAL[1], label="coolant")
    ax.set_yscale("log")
    ax.set_ylabel("film coefficient [kW/(m^2 K)]")
    ax.set_title("Convective film coefficients")
    ax.legend(ncols=2, loc="upper right")

    ax = axes[1]
    ax.plot(x, t["q_total"] * 1e-6, color=style.CATEGORICAL[0], label="total")
    if t["q_radiative"].max() > 0:
        ax.plot(x, t["q_convective"] * 1e-6, color=style.CATEGORICAL[1], label="convective")
        ax.plot(x, t["q_radiative"] * 1e-6, color=style.CATEGORICAL[2], label="radiative")
        ax.legend(ncols=3, loc="upper right")
    k = int(t["q_total"].idxmax())
    ax.plot([x.iloc[k]], [t["q_total"].iloc[k] * 1e-6], "o", color=style.CATEGORICAL[0], zorder=5)
    style.annotate(ax, x.iloc[k], t["q_total"].iloc[k] * 1e-6,
                   f"peak {t['q_total'].iloc[k] * 1e-6:.1f} MW/m^2")
    ax.set_ylabel("wall heat flux [MW/m^2]")
    ax.set_title("Wall heat flux")

    ax = axes[2]
    ax.plot(x, t["t_adiabatic_wall"], color=style.INK_MUTED, linewidth=1.3,
            label="adiabatic wall")
    ax.plot(x, t["t_wall_hot"], color=style.CATEGORICAL[0], label="hot-side wall")
    ax.plot(x, t["t_wall_cold"], color=style.CATEGORICAL[1], label="cold-side wall")
    ax.plot(x, t["coolant_T"], color=style.CATEGORICAL[2], label="coolant bulk")
    k = int(t["t_wall_hot"].idxmax())
    ax.plot([x.iloc[k]], [t["t_wall_hot"].iloc[k]], "o", color=style.CATEGORICAL[0], zorder=5)
    style.annotate(ax, x.iloc[k], t["t_wall_hot"].iloc[k],
                   f"peak {t['t_wall_hot'].iloc[k]:.0f} K")
    ax.set_ylabel("temperature [K]")
    ax.set_title("Temperatures through the wall")
    # The adiabatic-wall line runs across the top of this panel, so the legend
    # goes into the empty band between the wall temperatures and it.
    ax.legend(ncols=2, loc="center left", bbox_to_anchor=(0.02, 0.55))

    ax = axes[3]
    ax.plot(x, t["coolant_p"] * MPA, color=style.CATEGORICAL[0])
    ax.set_ylabel("coolant pressure [MPa]")
    ax.set_xlabel("axial position [mm]")
    ax.set_title("Coolant pressure")
    ax.text(0.03, 0.92, f"{(t['coolant_p'].iloc[-1] - t['coolant_p'].iloc[0]) * MPA:+.2f} MPa "
            f"across the jacket\n(counter-flow: inlet at the nozzle end)",
            transform=ax.transAxes, va="top", fontsize=8.5, color=style.INK_SECONDARY)

    for ax in axes:
        ax.axvline(xq, color=style.INK_MUTED, linewidth=0.9)
    style.suptitle(fig, f"{case}: regenerative cooling along the chamber and nozzle")
    style.caption(fig, "Bartz hot-gas correlation with a coupled wall and coolant balance. "
                       "Engineering estimate, not a conjugate CFD solution.")
    return _save(fig, out)


def cooling_map(sweep_csv: str, out: str) -> str:
    """Peak wall temperature and jacket pressure drop over the channel design space."""
    t = read_table(require(sweep_csv))
    h = np.sort(t["cooling.channel_height"].unique())
    w = np.sort(t["cooling.wall_thickness"].unique())
    def grid(col):
        g = np.full((len(w), len(h)), np.nan)
        for _, row in t.frame.iterrows():
            i = int(np.argmin(np.abs(w - row["cooling.wall_thickness"])))
            j = int(np.argmin(np.abs(h - row["cooling.channel_height"])))
            g[i, j] = row[col]
        return g

    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.3))
    for ax, col, label, scale in (
            (axes[0], "cooling.max_wall_temperature", "peak hot-wall temperature [K]", 1.0),
            (axes[1], "cooling.pressure_drop", "jacket pressure drop [MPa]", MPA)):
        g = grid(col) * scale
        im = ax.pcolormesh(h * MM, w * MM, g, cmap=style.SEQUENTIAL, shading="nearest")
        cs = ax.contour(h * MM, w * MM, g, colors=[style.INK_SECONDARY], linewidths=0.8,
                        levels=6)
        ax.clabel(cs, inline=True, fontsize=7.5, fmt="%.0f" if scale == 1.0 else "%.1f")
        ax.set_xlabel("channel height [mm]")
        ax.set_ylabel("wall thickness [mm]")
        ax.set_title(label)
        ax.grid(False)
        cb = fig.colorbar(im, ax=ax, pad=0.015, fraction=0.045)
        cb.outline.set_visible(False)
        # Points the jacket cannot pass are left blank rather than filled with a
        # clamped value, so the blank band is the answer and is named as such.
        blank = np.isnan(g)
        if blank.any():
            x_blank = h[np.where(blank.any(axis=0))[0]].max() * MM
            ax.text(0.03, 0.5, "infeasible:\nthe jacket cannot\npass the flow",
                    transform=ax.transAxes, va="center", fontsize=8,
                    color=style.INK_MUTED)
            ax.axvline(0.5 * (x_blank + h[np.where(~blank.any(axis=0))[0]].min() * MM),
                       color=style.GRID, linewidth=1.0)
    style.suptitle(fig, "Cooling-channel design space")
    style.caption(fig, "Deeper channels drop less pressure but cool less well; a thinner wall "
                       "runs cooler but has less structural margin. The trade is visible here.")
    return _save(fig, out)


# ---------------------------------------------------------------------------
# Transient
# ---------------------------------------------------------------------------

def transient_history(transient_csv: str, out: str) -> str:
    """Chamber pressure, temperature, mixture ratio and flows through a start."""
    t = read_table(require(transient_csv))
    ms = t["time"] * 1e3

    fig, axes = plt.subplots(4, 1, figsize=(7.0, 8.8), sharex=True)

    ax = axes[0]
    ax.plot(ms, t["pressure"] * MPA, color=style.CATEGORICAL[0])
    k = int(t["pressure"].idxmax())
    style.annotate(ax, ms.iloc[k], t["pressure"].iloc[k] * MPA,
                   f"peak {t['pressure'].iloc[k] * MPA:.2f} MPa", dx=14, dy=-4, va="top")
    ax.set_ylabel("chamber pressure [MPa]")
    ax.set_title("Chamber pressure")
    # Everything interesting happens in the first 60 ms; the rest is a plateau.
    zoom = ax.inset_axes([0.40, 0.14, 0.34, 0.62])
    early = ms <= 60.0
    zoom.plot(ms[early], t["pressure"][early] * MPA, color=style.CATEGORICAL[0])
    zoom.set_xlim(0.0, 60.0)
    zoom.tick_params(labelsize=7.5)
    zoom.set_title("first 60 ms", fontsize=8.5, color=style.INK_SECONDARY,
                   fontweight="normal", pad=3)

    ax = axes[1]
    ax.plot(ms, t["temperature"], color=style.CATEGORICAL[0])
    ax.set_ylabel("chamber temperature [K]")
    ax.set_title("Chamber temperature")

    ax = axes[2]
    ax.plot(ms, t["mixture_ratio"], color=style.CATEGORICAL[0])
    ax.set_ylabel("chamber O/F [-]")
    ax.set_title("Chamber mixture ratio")
    # The fuel-lead dip is a start-up feature: the run-wide minimum belongs to
    # the shutdown blow-down, which is a different thing entirely.
    start = t["time"] <= 0.10
    j = int(np.argmin(t["mixture_ratio"][start].to_numpy()))
    style.annotate(ax, ms[start].iloc[j], t["mixture_ratio"][start].iloc[j],
                   "fuel-lead dip at ignition", dx=12, dy=-2, va="center")

    ax = axes[3]
    ax.plot(ms, t["mdot_ox_in"], color=style.CATEGORICAL[0], label="oxidiser in")
    ax.plot(ms, t["mdot_fuel_in"], color=style.CATEGORICAL[1], label="fuel in")
    ax.plot(ms, t["mdot_out"], color=style.CATEGORICAL[2], label="nozzle out")
    ax.set_ylabel("mass flow [kg/s]")
    ax.set_xlabel("time [ms]")
    ax.set_title("Mass flows")
    ax.legend(ncols=3, loc="center right")

    style.suptitle(fig, "Startup and shutdown transient")
    style.caption(fig, "Zero-dimensional chamber mass and energy balance with a tabulated "
                       "equilibrium equation of state; valid from ignition onwards.")
    return _save(fig, out)


# ---------------------------------------------------------------------------
# Optimisation and uncertainty
# ---------------------------------------------------------------------------

def optimization_convergence(history_csv: str, json_path: str, out: str) -> str:
    """Objective and constraint violation against evaluation count."""
    t = read_table(require(history_csv))
    meta = read_json(require(json_path))["optimization"]
    objective = meta["objective_metric"]
    feasible = t.frame["feasible"] > 0.5

    fig, axes = plt.subplots(2, 1, figsize=(7.2, 5.8), sharex=True)

    ax = axes[0]
    ax.scatter(t["evaluation"][~feasible], t[objective][~feasible], s=7,
               color=style.INK_MUTED, label="infeasible", alpha=0.55, linewidths=0)
    ax.scatter(t["evaluation"][feasible], t[objective][feasible], s=9,
               color=style.CATEGORICAL[0], label="feasible", linewidths=0)
    best = t.frame[feasible][objective].cummax()
    ax.plot(t["evaluation"][feasible], best, color=style.CATEGORICAL[1],
            label="best so far", linewidth=1.8, zorder=6)
    ax.axhline(meta["objective"], color=style.INK_MUTED, linewidth=0.9)
    style.annotate(ax, t["evaluation"].min(), meta["objective"],
                   f"best feasible {meta['objective']:.2f} s", dx=2, dy=4)
    ax.set_ylabel(t.label(objective, objective.split(".")[-1].replace("_", " ")))
    ax.set_title("Objective")
    ax.legend(ncols=3, loc="lower right")
    lo = np.nanpercentile(t.frame[objective][np.isfinite(t.frame[objective])], 5)
    ax.set_ylim(bottom=lo)

    ax = axes[1]
    # The violation jumps between feasible and infeasible on consecutive
    # evaluations: joining them with a line paints solid blocks and hides the
    # trend, so the samples are drawn as marks and the running best as a line.
    v = np.maximum(t["max_violation"].to_numpy(), 1e-12)
    ax.scatter(t["evaluation"], v, s=5, linewidths=0, alpha=0.40,
               color=style.CATEGORICAL[0])
    ax.set_yscale("log")
    ax.set_ylim(5e-13, max(2.0, v.max() * 2.0))
    ax.set_ylabel("worst normalised violation [-]")
    ax.set_xlabel("objective evaluation [-]")
    ax.set_title("Constraint violation")
    ax.text(0.01, 0.06, "marks on the floor are feasible designs; the augmented Lagrangian "
                        "walks outside the feasible set and is pulled back",
            transform=ax.transAxes, fontsize=8, color=style.INK_MUTED)

    style.caption(fig, "Augmented Lagrangian around a Nelder-Mead simplex with "
                       f"Latin-hypercube multi-start. {meta['evaluations']} evaluations, "
                       f"{meta['failed_evaluations']} failed analyses.")
    return _save(fig, out)


def monte_carlo_histograms(samples_csv: str, json_path: str, out: str) -> str:
    """Output distributions with their 5th, 50th and 95th percentiles."""
    t = read_table(require(samples_csv))
    meta = read_json(require(json_path))["monte_carlo"]
    # The sample matrix carries the dispersed INPUTS first and the requested
    # outputs after them, and both families share the same name prefixes.
    # Selecting on the prefix alone would histogram the inputs, so the outputs
    # are taken from the statistics table and the input names removed.
    dispersed = {d["parameter"] for d in meta["inputs"]}
    metrics = [m for m in meta["statistics"]["columns"]["metric"]
               if m not in dispersed and "residual" not in m and m in t.columns
               and t.frame[m].notna().sum() > 10]
    preferred = ["performance.thrust", "performance.isp", "performance.isp_vacuum",
                 "cooling.max_wall_temperature", "cooling.max_heat_flux",
                 "cooling.pressure_drop"]
    outputs = [m for m in preferred if m in metrics]
    outputs += [m for m in metrics if m not in outputs]
    outputs = outputs[:6]

    fig, axes = plt.subplots(2, 3, figsize=(10.6, 5.8))
    for ax, col in zip(axes.ravel(), outputs):
        v = t.frame[col].dropna().to_numpy()
        ax.hist(v, bins=45, color=style.CATEGORICAL[0], alpha=0.9, linewidth=0)
        for q, lw in ((5, 1.0), (50, 1.6), (95, 1.0)):
            ax.axvline(np.percentile(v, q), color=style.INK_SECONDARY, linewidth=lw)
        ax.text(0.97, 0.94, f"p50 {np.percentile(v, 50):.4g}", transform=ax.transAxes,
                ha="right", va="top", fontsize=8.5, color=style.INK_SECONDARY)
        ax.set_xlabel(t.label(col, col.split(".")[-1].replace("_", " ")))
        ax.set_ylabel("samples [-]")
        ax.set_title(col, fontsize=10)
    for ax in axes.ravel()[len(outputs):]:
        ax.set_visible(False)
    style.suptitle(fig, "Monte Carlo output distributions")
    style.caption(fig, f"{meta['samples_succeeded']} of {meta['samples_requested']} samples "
                       "succeeded. Vertical rules mark the 5th, 50th and 95th percentiles.")
    return _save(fig, out)


def sensitivity_tornado(sensitivity_csv: str, out: str, outputs: list[str] | None = None) -> str:
    """Standardised regression coefficients, ranked, for a few outputs."""
    t = read_table(require(sensitivity_csv))
    available = list(dict.fromkeys(t.frame["output"]))
    chosen = outputs or available[:4]
    chosen = [c for c in chosen if c in available][:4]

    fig, axes = plt.subplots(1, len(chosen), figsize=(3.5 * len(chosen), 4.6), sharey=False)
    if len(chosen) == 1:
        axes = [axes]
    for ax, out_name in zip(axes, chosen):
        sub = t.frame[t.frame["output"] == out_name].copy()
        sub["mag"] = sub["src"].abs()
        sub = sub.sort_values("mag").tail(8)
        colours = [style.POLE_POSITIVE if v >= 0 else style.POLE_NEGATIVE for v in sub["src"]]
        ax.barh(range(len(sub)), sub["src"], color=colours, height=0.62, linewidth=0)
        ax.set_yticks(range(len(sub)))
        ax.set_yticklabels([n.split(".")[-1].replace("_", " ") for n in sub["input"]],
                           fontsize=8.5)
        ax.axvline(0.0, color=style.INK_SECONDARY, linewidth=0.9)
        ax.set_xlabel("standardised regression coefficient [-]")
        r2 = float(sub["r_squared"].iloc[0])
        ax.set_title(f"{out_name}\nlinear model R^2 = {r2:.3f}", fontsize=9.5)
        ax.grid(axis="y", visible=False)
        lim = max(0.15, 1.15 * sub["mag"].max())
        ax.set_xlim(-lim, lim)
    style.suptitle(fig, "Which uncertain inputs drive which outputs")
    style.caption(fig, "Sign is polarity, not identity: blue raises the output, red lowers it. "
                       "Magnitudes are comparable within a panel, not across panels.")
    return _save(fig, out)


def monte_carlo_scatter(samples_csv: str, out: str) -> str:
    """The dominant input/output relationship, with the failure flags shown."""
    t = read_table(require(samples_csv))
    ok = t.frame["status"] == "ok"
    fig, axes = plt.subplots(1, 2, figsize=(9.4, 4.2))

    ax = axes[0]
    ax.scatter(t.frame.loc[ok, "cooling.bartz_multiplier"],
               t.frame.loc[ok, "cooling.max_wall_temperature"],
               s=7, linewidths=0, alpha=0.55, color=style.CATEGORICAL[0])
    hot = ok & (t.frame["wall_limit_exceeded"] > 0.5)
    if hot.any():
        ax.scatter(t.frame.loc[hot, "cooling.bartz_multiplier"],
                   t.frame.loc[hot, "cooling.max_wall_temperature"],
                   s=9, linewidths=0, color=style.STATUS["critical"],
                   label=f"above the material limit ({int(hot.sum())} of {int(ok.sum())})")
        ax.legend(loc="upper left")
    ax.set_xlabel("Bartz correlation multiplier [-]")
    ax.set_ylabel("peak hot-wall temperature [K]")
    ax.set_title("The correlation drives the wall temperature")

    ax = axes[1]
    ax.scatter(t.frame.loc[ok, "performance.isp"], t.frame.loc[ok, "performance.thrust"] * KN,
               s=7, linewidths=0, alpha=0.55, color=style.CATEGORICAL[0])
    ax.set_xlabel("specific impulse [s]")
    ax.set_ylabel("thrust [kN]")
    ax.set_title("Thrust against specific impulse")
    style.caption(fig, "Left: red marks the samples whose peak wall temperature exceeds the "
                       "CuCrZr limit. At most three categorical hues are used in scatter form, "
                       "where every pair of colours is on screen at once.")
    return _save(fig, out)
