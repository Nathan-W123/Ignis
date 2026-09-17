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


def startup_animation(transient_csv: str, profile_csv: str, out: str, fps: int = 20,
                      frames: int = 130, dpi: int = 96) -> str:
    """Chamber pressure building through a start, with the nozzle filling.

    Three panels, no shared scales: the left panel replays the zero-dimensional
    transient; the upper right panel is the nozzle contour coloured by the
    converged steady Mach number; the lower right panel carries the axial static
    pressure on its own axis, scaled by the instantaneous chamber pressure.
    That last panel is a quasi-steady approximation -- the same one the
    transient's own thrust estimate makes -- and the caption says so rather than
    implying the nozzle flow was solved at every instant.
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

    fig = plt.figure(figsize=(11.0, 6.4))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 0.95], hspace=0.06)
    ax_noz = fig.add_subplot(gs[0, :])
    ax_hist = fig.add_subplot(gs[1, 0])
    ax_press = fig.add_subplot(gs[1, 1])

    # --- top: the contour, coloured by magnitude --------------------------
    # An axisymmetric half-section: the contour is mirror-symmetric, and
    # drawing one half lets the panel stay true to scale without going tall.
    norm = plt.Normalize(0.0, mach.max())
    pts = np.array([x, r]).T.reshape(-1, 1, 2)
    segs = np.concatenate([pts[:-1], pts[1:]], axis=1)
    lc = LineCollection(segs, cmap=style.SEQUENTIAL, norm=norm, linewidth=2.8)
    lc.set_array(0.5 * (mach[:-1] + mach[1:]))
    ax_noz.add_collection(lc)
    ax_noz.fill_between(x, 0.0, r, color=style.SEQUENTIAL(0.05), zorder=0)
    ax_noz.axhline(0.0, color=style.GRID, linewidth=1.0, zorder=1)
    x_throat = x[int(np.argmin(r))]
    ax_noz.set_xlim(x.min(), x.max())
    ax_noz.set_ylim(-18.0, r.max() * 1.12)
    ax_noz.set_aspect("equal")
    ax_noz.set_anchor("C")
    ax_noz.set_xlabel("axial position [mm]")
    ax_noz.set_ylabel("radius [mm]")
    ax_noz.set_title("Ignis-M1 nozzle, coloured by steady Mach number")
    cb = fig.colorbar(lc, ax=ax_noz, pad=0.012, fraction=0.020)
    cb.set_label("Mach number [-]", color=style.INK_SECONDARY, fontsize=9)
    cb.outline.set_visible(False)
    style.annotate(ax_noz, x_throat, float(r.min()), "throat", dx=0, dy=-13,
                   ha="center", va="top", color=style.INK_MUTED)
    ax_noz.text(0.004, 0.06, "centreline (axisymmetric half-section)",
                transform=ax_noz.transAxes, fontsize=8, color=style.INK_MUTED)

    # --- bottom left: the transient itself --------------------------------
    ax_hist.plot(t_ms, p * MPA, color=style.GRID, linewidth=1.4)
    (live,) = ax_hist.plot([], [], color=style.CATEGORICAL[0], linewidth=2.2)
    (head,) = ax_hist.plot([], [], "o", color=style.CATEGORICAL[0], markersize=6,
                           markeredgecolor=style.SURFACE, markeredgewidth=1.5)
    ax_hist.set_xlim(t_ms.min(), t_ms.max())
    ax_hist.set_ylim(0, p.max() * MPA * 1.12)
    ax_hist.set_xlabel("time [ms]")
    ax_hist.set_ylabel("chamber pressure [MPa]")
    ax_hist.set_title("Chamber pressure (0-D transient)")
    # The readout sits in the empty quadrant under the pressure plateau, so it
    # never collides with the rising edge or with the plateau itself.
    readout = ax_hist.text(0.40, 0.62, "", transform=ax_hist.transAxes, va="top",
                           fontsize=9.5, color=style.INK_SECONDARY, family="monospace")

    # --- bottom right: axial pressure, on its own scale -------------------
    (pline,) = ax_press.plot([], [], color=style.CATEGORICAL[1], linewidth=2.0)
    ax_press.set_xlim(x.min(), x.max())
    ax_press.set_ylim(0.0, p.max() * MPA * 1.12)
    ax_press.set_xlabel("axial position [mm]")
    ax_press.set_ylabel("static pressure [MPa]")
    ax_press.set_title("Axial static pressure (quasi-steady)")
    ax_press.axvline(x_throat, color=style.GRID, linewidth=1.2, zorder=0)

    style.caption(fig, "Top and bottom right: the converged steady solution, scaled by the "
                       "instantaneous chamber pressure -- the same quasi-steady\napproximation the "
                       "transient's own thrust estimate makes. Bottom left: the transient that was "
                       "actually integrated.")

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
    if out.endswith(".gif"):
        # A GIF is a palette format: quantising to 64 colours costs nothing
        # visible on a two-hue figure and roughly halves the committed file.
        writer = manim.PillowWriter(fps=fps)
        anim.save(out, writer=writer, dpi=dpi,
                  savefig_kwargs={"facecolor": style.SURFACE})
        _quantise_gif(out, colors=128)
    else:
        anim.save(out, writer=manim.FFMpegWriter(fps=fps), dpi=dpi)
    plt.close(fig)
    return out


def _quantise_gif(path: str, colors: int = 64) -> None:
    """Re-encode a GIF with a shared adaptive palette, in place."""
    from PIL import Image, ImageSequence

    with Image.open(path) as im:
        duration = im.info.get("duration", 50)
        frames = [f.convert("RGB") for f in ImageSequence.Iterator(im)]
    palette = frames[len(frames) // 2].quantize(colors=colors, method=Image.MEDIANCUT)
    quantised = [f.quantize(palette=palette, dither=Image.NONE) for f in frames]
    # disposal=1 ("leave in place") lets the encoder store only the changed
    # bounding box of each frame; most of this figure never moves.
    quantised[0].save(path, save_all=True, append_images=quantised[1:],
                      duration=duration, loop=0, optimize=True, disposal=1)
