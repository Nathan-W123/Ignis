"""Animations built from Ignis solver output."""
from __future__ import annotations

import os

import matplotlib.animation as manim
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection

from . import style
from .data import read_table, require

MM = 1e3
MPA = 1e-6


def startup_animation(transient_csv: str, profile_csv: str, out: str, fps: int = 25,
                      frames: int = 200) -> str:
    """Chamber pressure building through a start, with the nozzle filling.

    The left panel replays the zero-dimensional transient.  The right panel
    shows the steady axial pressure profile scaled by the instantaneous chamber
    pressure, which is the quasi-steady approximation the transient's thrust
    estimate already makes; the caption says so rather than implying the nozzle
    flow was solved at every instant.
    """
    tr = read_table(require(transient_csv))
    pr = read_table(require(profile_csv))

    t_ms = tr["time"].to_numpy() * 1e3
    p = tr["pressure"].to_numpy()
    T = tr["temperature"].to_numpy()
    mr = tr["mixture_ratio"].to_numpy()
    idx = np.linspace(0, len(t_ms) - 1, min(frames, len(t_ms))).astype(int)

    x = pr["x"].to_numpy() * MM
    r = pr["radius"].to_numpy() * MM
    p_shape = pr["pressure"].to_numpy() / pr["pressure"].to_numpy()[0]
    mach = pr["mach"].to_numpy()

    fig, (ax_hist, ax_noz) = plt.subplots(
        1, 2, figsize=(11.0, 4.3), gridspec_kw={"width_ratios": [1.0, 1.25]})

    ax_hist.plot(t_ms, p * MPA, color=style.GRID, linewidth=1.4)
    (live,) = ax_hist.plot([], [], color=style.CATEGORICAL[0], linewidth=2.2)
    (head,) = ax_hist.plot([], [], "o", color=style.CATEGORICAL[0], markersize=6,
                           markeredgecolor=style.SURFACE, markeredgewidth=1.5)
    ax_hist.set_xlim(t_ms.min(), t_ms.max())
    ax_hist.set_ylim(0, p.max() * MPA * 1.12)
    ax_hist.set_xlabel("time [ms]")
    ax_hist.set_ylabel("chamber pressure [MPa]")
    ax_hist.set_title("Chamber pressure")
    readout = ax_hist.text(0.03, 0.94, "", transform=ax_hist.transAxes, va="top",
                           fontsize=9.5, color=style.INK_SECONDARY, family="monospace")

    norm = plt.Normalize(0.0, mach.max())
    pts = np.array([x, r]).T.reshape(-1, 1, 2)
    segs = np.concatenate([pts[:-1], pts[1:]], axis=1)
    lc = LineCollection(segs, cmap=style.SEQUENTIAL, norm=norm, linewidth=2.6)
    lc.set_array(0.5 * (mach[:-1] + mach[1:]))
    ax_noz.add_collection(lc)
    pts2 = np.array([x, -r]).T.reshape(-1, 1, 2)
    segs2 = np.concatenate([pts2[:-1], pts2[1:]], axis=1)
    lc2 = LineCollection(segs2, cmap=style.SEQUENTIAL, norm=norm, linewidth=2.6)
    lc2.set_array(0.5 * (mach[:-1] + mach[1:]))
    ax_noz.add_collection(lc2)
    fill = ax_noz.fill_between(x, -r, r, color=style.SEQUENTIAL(0.05), zorder=0)

    ax_p = ax_noz.twiny()  # a second *x* scale is fine: it shares the y axis with nothing
    ax_p.set_visible(False)
    ax_press = ax_noz.inset_axes([0.0, 0.0, 1.0, 1.0], zorder=3)
    ax_press.patch.set_alpha(0.0)
    (pline,) = ax_press.plot([], [], color=style.CATEGORICAL[1], linewidth=2.0)
    ax_press.set_xlim(x.min(), x.max())
    ax_press.set_ylim(0.0, p.max() * MPA * 1.12)
    ax_press.set_yticks([])
    ax_press.set_xticks([])
    for spine in ax_press.spines.values():
        spine.set_visible(False)
    ax_press.grid(False)
    ax_press.text(0.985, 0.94, "axial static pressure", transform=ax_press.transAxes,
                  ha="right", va="top", fontsize=8.5, color=style.CATEGORICAL[1])

    ax_noz.set_xlim(x.min(), x.max())
    ax_noz.set_ylim(-r.max() * 1.2, r.max() * 1.2)
    ax_noz.set_aspect("equal")
    ax_noz.set_xlabel("axial position [mm]")
    ax_noz.set_ylabel("radius [mm]")
    ax_noz.set_title("Nozzle, coloured by steady Mach number")
    cb = fig.colorbar(lc, ax=ax_noz, pad=0.015, fraction=0.035)
    cb.set_label("Mach number [-]", color=style.INK_SECONDARY, fontsize=9)
    cb.outline.set_visible(False)

    style.caption(fig, "Left: the zero-dimensional transient. Right: the converged steady axial "
                       "pressure profile scaled by the instantaneous chamber pressure, i.e. the "
                       "same quasi-steady approximation the transient thrust estimate makes.")

    def update(k):
        i = idx[k]
        live.set_data(t_ms[:i + 1], p[:i + 1] * MPA)
        head.set_data([t_ms[i]], [p[i] * MPA])
        readout.set_text(f"t   {t_ms[i]:7.1f} ms\n"
                         f"pc  {p[i] * MPA:7.3f} MPa\n"
                         f"Tc  {T[i]:7.0f} K\n"
                         f"O/F {mr[i]:7.3f}")
        pline.set_data(x, p_shape * p[i] * MPA)
        return live, head, readout, pline

    anim = manim.FuncAnimation(fig, update, frames=len(idx), blit=False, interval=1000 / fps)
    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)
    writer = manim.PillowWriter(fps=fps) if out.endswith(".gif") else manim.FFMpegWriter(fps=fps)
    anim.save(out, writer=writer, dpi=110)
    plt.close(fig)
    return out
