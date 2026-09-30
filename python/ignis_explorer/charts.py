"""Matplotlib viewports for the Ignis Engine Explorer.

Every chart takes one or two solved designs.  Colour is assigned by the job it
does (see styles.py): the two design slots are an identity and keep fixed hues,
Mach number along the contour is a magnitude and gets the sequential ramp, and
constraint state is status colour that is never reused for a series.
"""
from __future__ import annotations

import textwrap
from typing import List, Optional, Sequence

import numpy as np
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.collections import LineCollection
from matplotlib.figure import Figure

from . import styles
from .solver import Result

MM = 1e3
MPA = 1e-6
KN = 1e-3


class Chart(FigureCanvasQTAgg):
    """A figure that redraws itself from the current design slots."""

    title = ""
    caption = ""

    def __init__(self, parent=None, figsize=(6.0, 4.0)) -> None:
        self.figure = Figure(figsize=figsize)
        super().__init__(self.figure)
        self.setParent(parent)
        self.setMinimumHeight(240)

    def show_results(self, results: Sequence[Optional[Result]], labels: Sequence[str]) -> None:
        self.figure.clear()
        live = [(r, lab) for r, lab in zip(results, labels) if r is not None and r.ok]
        if not live:
            ax = self.figure.add_subplot(111)
            ax.text(0.5, 0.5, "run a design to see this chart", ha="center", va="center",
                    transform=ax.transAxes, color=styles.TEXT_DIM, fontsize=10)
            ax.set_xticks([]); ax.set_yticks([]); ax.grid(False)
            for sp in ax.spines.values():
                sp.set_visible(False)
        else:
            self.draw_chart([r for r, _ in live], [lab for _, lab in live])
            if self.caption:
                # Wrap to the canvas: an unwrapped caption runs off a narrow
                # panel, and the panel is resizable.
                width = max(48, int(self.figure.get_figwidth() * 15))
                self.figure.supxlabel("\n".join(textwrap.wrap(self.caption, width)),
                                      color=styles.TEXT_DIM, fontsize=7.5,
                                      ha="left", x=0.005)
        self.draw_idle()

    def draw_chart(self, results: List[Result], labels: List[str]) -> None:  # pragma: no cover
        raise NotImplementedError

    @staticmethod
    def slot_colour(label: str) -> str:
        return styles.DESIGN_A if label.startswith("A") else styles.DESIGN_B


class ContourChart(Chart):
    """Nozzle half-section coloured by the local Mach number."""

    title = "Nozzle contour"
    caption = ("Axisymmetric half-section, true to scale. Colour is the quasi-1D Mach number, "
               "which varies only along the axis.")

    def draw_chart(self, results, labels):
        n = len(results)
        axes = self.figure.subplots(n, 1, squeeze=False)[:, 0]
        vmax = max(np.nanmax(r.profile.get("mach", [1.0])) for r in results)
        norm = __import__("matplotlib").colors.Normalize(0.0, vmax)
        lc = None
        for ax, r, lab in zip(axes, results, labels):
            x = np.asarray(r.profile["x"]) * MM
            rad = np.asarray(r.profile["radius"]) * MM
            mach = np.asarray(r.profile["mach"])
            pts = np.array([x, rad]).T.reshape(-1, 1, 2)
            segs = np.concatenate([pts[:-1], pts[1:]], axis=1)
            lc = LineCollection(segs, cmap=styles.SEQUENTIAL, norm=norm, linewidth=3.2)
            lc.set_array(0.5 * (mach[:-1] + mach[1:]))
            ax.add_collection(lc)
            # The interior is a flat, recessive fill rather than a ramp step: the
            # ramp is the Mach encoding, and a filled area in the same family
            # would compete with the wall line that carries it.
            ax.fill_between(x, 0.0, rad, color="#243447", zorder=0)
            ax.axhline(0.0, color=styles.BORDER, linewidth=1.0, zorder=1)
            xt = x[int(np.argmin(rad))]
            ax.axvline(xt, color=styles.BORDER, linewidth=0.9, zorder=1)
            ax.set_xlim(x.min(), x.max())
            ax.set_ylim(-0.04 * rad.max(), rad.max() * 1.12)
            ax.set_aspect("equal")
            ax.set_anchor("C")
            ax.set_ylabel("radius [mm]")
            eps = r.get("performance.expansion_ratio", r.design.expansion_ratio)
            ax.set_title(f"{lab}   ε = {eps:.1f},  exit ⌀ "
                         f"{2 * rad.max():.0f} mm,  Mₑ = "
                         f"{r.get('performance.exit_mach'):.2f}",
                         color=self.slot_colour(lab))
        axes[-1].set_xlabel("axial position from the injector face [mm]")
        if lc is not None:
            cb = self.figure.colorbar(lc, ax=list(axes), pad=0.012, fraction=0.03)
            cb.set_label("Mach number [-]", color=styles.TEXT_MUTED, fontsize=8)
            cb.outline.set_visible(False)


class AxialChart(Chart):
    """Pressure, temperature, velocity and Mach number along the axis."""

    title = "Axial profiles"
    caption = "The vertical rule marks the sonic throat. Every panel comes from the same station march."

    def draw_chart(self, results, labels):
        axes = self.figure.subplots(2, 2).ravel()
        panels = [("pressure", MPA, "static pressure [MPa]", True),
                  ("temperature", 1.0, "static temperature [K]", False),
                  ("velocity", 1.0, "velocity [m/s]", False),
                  ("mach", 1.0, "Mach number [-]", False)]
        for ax, (col, scale, label, logy) in zip(axes, panels):
            for r, lab in zip(results, labels):
                x = np.asarray(r.profile["x"]) * MM
                y = np.asarray(r.profile[col]) * scale
                ax.plot(x, y, color=self.slot_colour(lab), label=lab)
                rad = np.asarray(r.profile["radius"])
                ax.axvline(x[int(np.argmin(rad))], color=styles.BORDER, linewidth=0.9, zorder=0)
            if logy:
                ax.set_yscale("log")
            ax.set_ylabel(label)
        for ax in axes[2:]:
            ax.set_xlabel("axial position [mm]")
        if len(results) > 1:
            axes[0].legend(ncols=2, loc="lower left")


class ThermalChart(Chart):
    """Heat flux and the wall / coolant temperatures along the jacket."""

    title = "Wall and coolant"
    caption = ("Bartz hot-gas correlation with a coupled wall and coolant balance. "
               "Engineering estimate, not a conjugate CFD solution.")

    def draw_chart(self, results, labels):
        have = [(r, lab) for r, lab in zip(results, labels) if r.thermal]
        if not have:
            ax = self.figure.add_subplot(111)
            ax.text(0.5, 0.5, "cooling is switched off for this design",
                    ha="center", va="center", transform=ax.transAxes,
                    color=styles.TEXT_DIM, fontsize=10)
            ax.set_xticks([]); ax.set_yticks([]); ax.grid(False)
            for sp in ax.spines.values():
                sp.set_visible(False)
            return
        ax_q, ax_t = self.figure.subplots(2, 1, sharex=True)
        for r, lab in have:
            x = np.asarray(r.thermal["x"]) * MM
            colour = self.slot_colour(lab)
            ax_q.plot(x, np.asarray(r.thermal["q_total"]) * 1e-6, color=colour, label=lab)
            # Hot-side wall solid, coolant bulk lighter: two surfaces of one
            # design share its hue, separated by weight rather than by a new hue.
            ax_t.plot(x, r.thermal["t_wall_hot"], color=colour, label=f"{lab} hot wall")
            ax_t.plot(x, r.thermal["coolant_T"], color=colour, linewidth=1.0, alpha=0.55,
                      label=f"{lab} coolant bulk")
            limit = r.get("cooling.wall_limit_temperature", float("nan"))
            if limit == limit and limit > 0:
                ax_t.axhline(limit, color=styles.STATUS["critical"], linewidth=0.9, zorder=0)
                material = r.strings.get("cooling.wall_material", "material")
                ax_t.annotate(f"{material} limit {limit:.0f} K", (0.995, limit),
                              xycoords=("axes fraction", "data"),
                              xytext=(0, 3), textcoords="offset points",
                              ha="right", va="bottom", fontsize=7.5,
                              color=styles.STATUS["critical"])
        ax_q.set_ylabel("wall heat flux [MW/m²]")
        ax_q.set_title("Wall heat flux")
        ax_t.set_ylabel("temperature [K]")
        ax_t.set_title("Hot-side wall and coolant bulk")
        ax_t.set_xlabel("axial position [mm]")
        ax_t.legend(ncols=2, loc="lower left")
        if len(have) > 1:
            ax_q.legend(ncols=2, loc="upper right")


class CompositionChart(Chart):
    """Chamber equilibrium composition, ranked."""

    title = "Chamber composition"
    caption = "Equilibrium mole fractions at the chamber stagnation state; species below 1e-5 are omitted."

    def draw_chart(self, results, labels):
        ax = self.figure.add_subplot(111)
        keep = 10
        names: List[str] = []
        for r in results:
            for n, _ in r.composition[:keep]:
                if n not in names:
                    names.append(n)
        names = names[:keep]
        idx = np.arange(len(names))
        width = 0.8 / max(1, len(results))
        for k, (r, lab) in enumerate(zip(results, labels)):
            table = dict(r.composition)
            vals = [max(table.get(n, 0.0), 1e-6) for n in names]
            ax.barh(idx + k * width - 0.4 + width / 2, vals, height=width * 0.88,
                    color=self.slot_colour(lab), label=lab, linewidth=0)
        ax.set_yticks(idx)
        ax.set_yticklabels(names, fontsize=8)
        ax.invert_yaxis()
        ax.set_xscale("log")
        ax.set_xlim(1e-5, 1.3)
        ax.set_xlabel("mole fraction [-]")
        ax.grid(axis="y", visible=False)
        if len(results) > 1:
            ax.legend(ncols=2, loc="lower right")


class AltitudeChart(Chart):
    """Thrust and specific impulse against altitude, from the converged exit state."""

    title = "Altitude sweep"
    caption = ("Thrust is exactly linear in ambient pressure for a full-flowing nozzle, so these "
               "curves follow from the converged exit state. A shaded band marks where the "
               "Summerfield criterion predicts separation; the curve there is optimistic because "
               "the inviscid solution is not modified.")

    def draw_chart(self, results, labels):
        # The curve comes from ignis_nozzle, not from a Python atmosphere: there
        # is one implementation of the U.S. Standard Atmosphere in this
        # repository and it is the validated C++ one.
        have = [(r, lab) for r, lab in zip(results, labels) if r.altitude]
        if not have:
            ax = self.figure.add_subplot(111)
            ax.text(0.5, 0.5, "no altitude sweep available for this design",
                    ha="center", va="center", transform=ax.transAxes,
                    color=styles.TEXT_DIM, fontsize=10)
            ax.set_xticks([]); ax.set_yticks([]); ax.grid(False)
            for sp in ax.spines.values():
                sp.set_visible(False)
            return
        ax_f, ax_i = self.figure.subplots(2, 1, sharex=True)
        for k, (r, lab) in enumerate(have):
            z = np.asarray(r.altitude["altitude"])
            keep = np.isfinite(z)
            z = z[keep]
            thrust = np.asarray(r.altitude["thrust"])[keep]
            isp = np.asarray(r.altitude["isp"])[keep]
            sep = np.asarray(r.altitude["separation_predicted"])[keep]
            ax_f.plot(z * 1e-3, thrust * KN, color=self.slot_colour(lab), label=lab)
            ax_i.plot(z * 1e-3, isp, color=self.slot_colour(lab), label=lab)
            if sep.any():
                # Each design gets its band in its own hue and its own text
                # line, so two overlapping bands stay readable.
                z_sep = z[sep > 0.5].max()
                colour = self.slot_colour(lab)
                for ax in (ax_f, ax_i):
                    ax.axvspan(0.0, z_sep * 1e-3, color=colour, alpha=0.12, zorder=0, lw=0)
                    ax.axvline(z_sep * 1e-3, color=colour, alpha=0.5, lw=0.9, zorder=0)
                ax_f.annotate(f"{lab}: separation predicted below {z_sep * 1e-3:.1f} km",
                              (0.5, 0.06 + 0.09 * k), xycoords="axes fraction",
                              ha="center", va="bottom", fontsize=7.5, color=colour)
        ax_f.set_ylabel("thrust [kN]")
        ax_i.set_ylabel("specific impulse [s]")
        ax_i.set_xlabel("geometric altitude [km]")
        ax_f.set_title("Thrust")
        ax_i.set_title("Specific impulse")
        if len(have) > 1:
            ax_f.legend(ncols=2, loc="lower right")
