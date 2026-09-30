"""Dark theme for the Ignis Engine Explorer.

The palette follows the Aero CFD Studio desktop GUI so the two tools read as
one family: a deep blue-grey chrome with a teal accent.

COLOUR CONTRACT
---------------
Chart colour is assigned by the job it does, never by taste:

* categorical (identity - a species, a design slot, a wall surface)
  -> ``CATEGORICAL``, taken in fixed slot order and never cycled.
* sequential (magnitude - Mach number along a contour) -> one hue, light to
  dark.
* status (state - a constraint that passes, warns or fails) -> ``STATUS``,
  reserved and never reused for a series.

``CATEGORICAL`` was validated against this module's chart surface (#1a2332) in
dark mode rather than chosen by eye: worst adjacent colour-vision separation
14.1 (OKLab dE x100, 8 is the target), worst adjacent normal-vision separation
16.2 (15 is the floor), every slot inside the dark lightness band and above
3:1 contrast against the surface.  Forms that put every pair on screen at once
use ``CATEGORICAL_ALLPAIRS``, the first three slots, which clear the stricter
all-pairs gate (11.5 / 17.4).
"""
from __future__ import annotations

import matplotlib.colors as mcolors

# --- chrome ---------------------------------------------------------------
BG_DEEP = "#0a0e14"
BG_MAIN = "#0f1419"
BG_PANEL = "#151c24"
BG_ELEVATED = "#1a2332"
BG_INPUT = "#1e2a3a"
BORDER = "#2d3a4f"
BORDER_FOCUS = "#3d7ea6"
TEXT = "#e2e8f0"
TEXT_MUTED = "#94a3b8"
TEXT_DIM = "#64748b"
ACCENT = "#2dd4bf"
ACCENT_BRIGHT = "#4de8ff"
ACCENT_DEEP = "#14b8a6"
ACCENT_HOVER = "#5eead4"

# --- categorical: fixed slot order, never cycled --------------------------
CATEGORICAL = [
    "#0891b2",  # 1 cyan
    "#16a34a",  # 2 green
    "#a855f7",  # 3 violet
    "#f43f5e",  # 4 rose
    "#0284c7",  # 5 blue
    "#b45309",  # 6 bronze
    "#ec4899",  # 7 pink
    "#ea580c",  # 8 orange
]
CATEGORICAL_ALLPAIRS = CATEGORICAL[:3]

# The two design slots are an identity, so they get fixed hues that never move.
DESIGN_A = "#4de8ff"
DESIGN_B = "#f0a020"

# --- sequential: one hue, dark to light on a dark surface -----------------
# The darkest step still has to clear the chart surface (#1a2332): on a dark
# background a ramp that starts at the surface makes its own low end invisible.
SEQUENTIAL_STEPS = [
    "#3580a2", "#3a8cb0", "#3f98be", "#44a4cc", "#4ab0d9", "#55bbe2",
    "#57bce3", "#6cc7ea", "#85d3f0", "#9fdef5", "#bae9fa", "#d6f3fd",
]
SEQUENTIAL = mcolors.LinearSegmentedColormap.from_list("ignis_dark_seq", SEQUENTIAL_STEPS)

# --- status: reserved, never reused for a series --------------------------
STATUS = {
    "good": "#22c55e",
    "warning": "#eab308",
    "serious": "#f97316",
    "critical": "#ef4444",
}

CHART = {
    "figure_bg": BG_PANEL,
    "axes_bg": BG_ELEVATED,
    "text": TEXT_MUTED,
    "grid": BORDER,
}


def apply_matplotlib() -> None:
    """Install the Explorer's chart style into matplotlib's defaults."""
    import matplotlib as mpl

    mpl.rcParams.update({
        "figure.facecolor": CHART["figure_bg"],
        "savefig.facecolor": CHART["figure_bg"],
        "axes.facecolor": CHART["axes_bg"],
        "axes.edgecolor": CHART["grid"],
        "axes.linewidth": 0.8,
        "axes.labelcolor": CHART["text"],
        "axes.titlecolor": TEXT,
        "axes.titlesize": 10.0,
        "axes.titleweight": "bold",
        "axes.titlelocation": "left",
        "axes.titlepad": 7,
        "axes.labelsize": 8.5,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "axes.axisbelow": True,
        "axes.prop_cycle": mpl.cycler(color=CATEGORICAL),
        "grid.color": CHART["grid"],
        "grid.linewidth": 0.6,
        "grid.linestyle": "-",     # solid hairlines; dashes read as a threshold
        "grid.alpha": 1.0,
        "lines.linewidth": 1.8,
        "lines.solid_capstyle": "round",
        "xtick.color": TEXT_DIM,
        "ytick.color": TEXT_DIM,
        "xtick.labelcolor": CHART["text"],
        "ytick.labelcolor": CHART["text"],
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "xtick.major.size": 3.0,
        "ytick.major.size": 3.0,
        "legend.frameon": False,
        "legend.fontsize": 8,
        "legend.labelcolor": CHART["text"],
        "legend.handlelength": 1.5,
        "font.size": 9,
        "figure.constrained_layout.use": True,
    })


STYLESHEET = f"""
QMainWindow, QWidget#centralRoot {{
    background-color: {BG_DEEP};
    color: {TEXT};
}}
QWidget {{ font-size: 12px; }}
QSplitter::handle {{ background-color: {BORDER}; }}
QSplitter::handle:horizontal {{ width: 1px; }}
QSplitter::handle:vertical {{ height: 1px; }}
QScrollArea {{ background-color: {BG_MAIN}; border: none; }}
QStatusBar {{
    background-color: {BG_PANEL};
    color: {TEXT_MUTED};
    border-top: 1px solid {BORDER};
}}
QFrame#ribbonBar {{
    background-color: {BG_PANEL};
    border-bottom: 1px solid {BORDER};
}}
QFrame#ribbonGroup {{
    background-color: transparent;
    border-right: 1px solid {BORDER};
}}
QLabel#ribbonGroupTitle {{
    color: {TEXT_DIM};
    font-size: 10px;
    font-weight: 600;
    letter-spacing: 1px;
}}
QLabel#appTitle {{
    color: {ACCENT_BRIGHT};
    font-size: 15px;
    font-weight: 700;
}}
QLabel#appSubtitle {{ color: {TEXT_DIM}; font-size: 11px; }}
QFrame#panelFrame {{
    background-color: {BG_PANEL};
    border: 1px solid {BORDER};
    border-radius: 6px;
}}
QLabel#panelTitle {{
    color: {ACCENT_BRIGHT};
    font-size: 12px;
    font-weight: 700;
    letter-spacing: 0.5px;
}}
QLabel#fieldLabel {{ color: {TEXT_MUTED}; font-size: 11px; }}
QLabel#unitLabel {{ color: {TEXT_DIM}; font-size: 10px; }}
QLabel#cardValue {{ color: {TEXT}; font-size: 19px; font-weight: 700; }}
QLabel#cardLabel {{ color: {TEXT_DIM}; font-size: 10px; letter-spacing: 0.6px; }}
QLabel#cardDelta {{ font-size: 11px; font-weight: 600; }}
QFrame#card {{
    background-color: {BG_ELEVATED};
    border: 1px solid {BORDER};
    border-radius: 6px;
}}
QDoubleSpinBox, QSpinBox, QComboBox, QLineEdit {{
    background-color: {BG_INPUT};
    color: {TEXT};
    border: 1px solid {BORDER};
    border-radius: 4px;
    padding: 3px 6px;
    selection-background-color: {ACCENT_DEEP};
}}
QDoubleSpinBox:focus, QSpinBox:focus, QComboBox:focus {{
    border: 1px solid {BORDER_FOCUS};
}}
QComboBox::drop-down {{ border: none; width: 16px; }}
QComboBox QAbstractItemView {{
    background-color: {BG_ELEVATED};
    color: {TEXT};
    border: 1px solid {BORDER};
    selection-background-color: {BG_INPUT};
    selection-color: {ACCENT_BRIGHT};
}}
QPushButton {{
    background-color: {BG_INPUT};
    color: {TEXT};
    border: 1px solid {BORDER};
    border-radius: 4px;
    padding: 5px 12px;
}}
QPushButton:hover {{ background-color: {BG_ELEVATED}; border-color: {BORDER_FOCUS}; }}
QPushButton:disabled {{ color: {TEXT_DIM}; border-color: {BORDER}; }}
QPushButton#runButton {{
    background-color: {ACCENT_DEEP};
    color: #04201c;
    font-weight: 700;
    border: none;
    padding: 6px 18px;
}}
QPushButton#runButton:hover {{ background-color: {ACCENT_HOVER}; }}
QPushButton#runButton:disabled {{ background-color: {BORDER}; color: {TEXT_DIM}; }}
QCheckBox {{ color: {TEXT_MUTED}; }}
QCheckBox::indicator {{
    width: 13px; height: 13px;
    border: 1px solid {BORDER};
    border-radius: 3px;
    background-color: {BG_INPUT};
}}
QCheckBox::indicator:checked {{ background-color: {ACCENT_DEEP}; border-color: {ACCENT_DEEP}; }}
QTabWidget::pane {{ border: 1px solid {BORDER}; border-radius: 6px; background-color: {BG_PANEL}; }}
QTabBar::tab {{
    background-color: {BG_MAIN};
    color: {TEXT_MUTED};
    border: 1px solid {BORDER};
    border-bottom: none;
    border-top-left-radius: 5px;
    border-top-right-radius: 5px;
    padding: 5px 14px;
    margin-right: 2px;
}}
QTabBar::tab:selected {{ background-color: {BG_PANEL}; color: {ACCENT_BRIGHT}; }}
QTextEdit, QPlainTextEdit {{
    background-color: {BG_ELEVATED};
    color: {TEXT_MUTED};
    border: 1px solid {BORDER};
    border-radius: 5px;
    font-size: 11px;
}}
QScrollBar:vertical {{ background: {BG_MAIN}; width: 10px; margin: 0; }}
QScrollBar::handle:vertical {{ background: {BORDER}; border-radius: 5px; min-height: 24px; }}
QScrollBar::handle:vertical:hover {{ background: {BORDER_FOCUS}; }}
QScrollBar::add-line, QScrollBar::sub-line {{ height: 0; }}
QToolTip {{
    background-color: {BG_ELEVATED};
    color: {TEXT};
    border: 1px solid {BORDER_FOCUS};
    padding: 4px;
}}
"""
