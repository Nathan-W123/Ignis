"""Shared plotting style for every Ignis figure.

COLOUR CONTRACT
---------------
One validated palette is used everywhere, and colour is assigned by the job it
does rather than by taste:

* categorical (identity - species, propellant, correlation)  -> CATEGORICAL,
  taken in fixed slot order and never cycled.  A ninth series folds into
  "other" rather than inventing a hue.
* sequential (magnitude - Mach number along a contour, an altitude ladder)
  -> one hue, light to dark.
* diverging (polarity - a signed sensitivity coefficient) -> blue/red with a
  neutral grey midpoint.

The palette is the validated reference instance: worst adjacent colour-vision
separation 9.1 (OKLab dE x100, 8 is the target) and worst adjacent
normal-vision separation 19.6 (15 is the floor).  Forms that put every pair on
screen at once (scatter, field plots) use at most the first three slots, which
clear the all-pairs gate.

Marks are thin, grids are solid hairlines one shade off the surface, and text
always wears an ink colour rather than a series colour, so identity is never
carried by colour alone.
"""
from __future__ import annotations

import matplotlib as mpl
import matplotlib.colors as mcolors

# --- categorical: fixed slot order, never cycled ---------------------------
CATEGORICAL = [
    "#2a78d6",  # 1 blue
    "#eb6834",  # 2 orange
    "#1baf7a",  # 3 aqua
    "#eda100",  # 4 yellow
    "#e87ba4",  # 5 magenta
    "#008300",  # 6 green
    "#4a3aa7",  # 7 violet
    "#e34948",  # 8 red
]
# Forms that show every pair simultaneously are capped at these three.
CATEGORICAL_ALLPAIRS = CATEGORICAL[:3]

# --- sequential: one hue, light to dark ------------------------------------
SEQUENTIAL_STEPS = [
    "#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec", "#5598e7",
    "#3987e5", "#2a78d6", "#256abf", "#1c5cab", "#184f95", "#104281", "#0d366b",
]
SEQUENTIAL = mcolors.LinearSegmentedColormap.from_list("ignis_seq", SEQUENTIAL_STEPS)
# Ordinal use (discrete ordered marks) must stay clear of the surface.
ORDINAL_STEPS = SEQUENTIAL_STEPS[3:]

# --- diverging: warm/cool poles, neutral midpoint --------------------------
DIVERGING = mcolors.LinearSegmentedColormap.from_list(
    "ignis_div", ["#104281", "#2a78d6", "#9ec5f4", "#f0efec", "#f3a9a8", "#e34948", "#a01f1e"])
POLE_POSITIVE = "#2a78d6"
POLE_NEGATIVE = "#e34948"

# --- status: reserved, never reused for a series ---------------------------
STATUS = {"good": "#008300", "warning": "#eda100", "serious": "#eb6834", "critical": "#e34948"}

# --- surfaces and ink ------------------------------------------------------
SURFACE = "#fcfcfb"
INK_PRIMARY = "#0b0b0b"
INK_SECONDARY = "#52514e"
INK_MUTED = "#8a8880"
GRID = "#e3e2dd"


def ordinal(n: int) -> list[str]:
    """`n` steps of the sequential ramp, spaced for discrete ordered series."""
    if n <= 1:
        return [ORDINAL_STEPS[len(ORDINAL_STEPS) // 2]]
    idx = [round(i * (len(ORDINAL_STEPS) - 1) / (n - 1)) for i in range(n)]
    return [ORDINAL_STEPS[i] for i in idx]


def apply() -> None:
    """Install the Ignis style into matplotlib's global defaults."""
    mpl.rcParams.update({
        "figure.facecolor": SURFACE,
        "figure.dpi": 130,
        "savefig.dpi": 160,
        "savefig.facecolor": SURFACE,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.25,
        "axes.facecolor": SURFACE,
        "axes.edgecolor": GRID,
        "axes.linewidth": 0.8,
        "axes.labelcolor": INK_SECONDARY,
        "axes.titlecolor": INK_PRIMARY,
        "axes.titlesize": 11.5,
        "axes.titleweight": "bold",
        "axes.titlelocation": "left",
        "axes.titlepad": 10,
        "axes.labelsize": 9.5,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "axes.axisbelow": True,
        "axes.prop_cycle": mpl.cycler(color=CATEGORICAL),
        "grid.color": GRID,
        "grid.linewidth": 0.7,
        "grid.linestyle": "-",      # solid hairlines; dashing reads as a threshold
        "grid.alpha": 1.0,
        "lines.linewidth": 2.0,
        "lines.markersize": 4.5,
        "lines.solid_capstyle": "round",
        "xtick.color": INK_MUTED,
        "ytick.color": INK_MUTED,
        "xtick.labelcolor": INK_SECONDARY,
        "ytick.labelcolor": INK_SECONDARY,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "xtick.major.size": 3.5,
        "ytick.major.size": 3.5,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        "legend.frameon": False,
        "legend.fontsize": 9,
        "legend.labelcolor": INK_SECONDARY,
        "legend.handlelength": 1.6,
        "legend.borderaxespad": 0.0,
        "font.size": 10,
        "font.family": "sans-serif",
        "figure.constrained_layout.use": True,
        "figure.titlesize": 12.5,
        "figure.titleweight": "bold",
    })


def caption(fig, text: str) -> None:
    """A single muted source/assumption line under a figure."""
    fig.supxlabel(text, color=INK_MUTED, fontsize=8, ha="left", x=0.005)


def annotate(ax, x, y, text, *, dx=6, dy=6, color=None, ha="left", va="bottom"):
    """A selective direct label: used for an endpoint or an extreme, never for
    every point."""
    ax.annotate(text, (x, y), textcoords="offset points", xytext=(dx, dy),
                fontsize=8.5, color=color or INK_SECONDARY, ha=ha, va=va)
