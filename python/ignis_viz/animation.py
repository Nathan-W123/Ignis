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


def firing_animation(transient_csv: str, profile_csv: str, out: str, *,
                     width: int = 1280, height: int = 720, fps: int = 30,
                     steps: int = 1600, bins: int = 240,
                     pitch_deg: float = 13.0, shell=None) -> str:
    """The engine lighting up and shutting down, as a volumetric render.

    The brightness of every frame is the solver's: the zero-dimensional start
    transient gives chamber pressure and temperature against time, and the
    nozzle field is scaled onto that state before its thermal emission is
    integrated along the camera rays.

    Two approximations are worth naming.  The nozzle is scaled quasi-steadily
    -- the same approximation the transient's own thrust estimate makes -- so
    the normalised shape of the axial profile is held fixed and only its level
    follows the chamber.  And the transient is driven onto the steady profile
    through *normalised* pressure and temperature, because the zero-dimensional
    chamber model carries a heat-loss efficiency the steady equilibrium
    solution does not, and the two would otherwise disagree at the same instant.

    Exposure is fixed across the whole sequence from the steady-state frame, so
    a frame that looks brighter is an engine that got brighter.
    """
    import numpy as np
    from PIL import Image, ImageDraw, ImageFont
    import imageio.v2 as imageio
    import pandas as pd

    from . import render as R

    tr = pd.read_csv(transient_csv, comment="#")
    profile = pd.read_csv(profile_csv, comment="#")
    field = R.Field.from_profile(profile, bins=bins)

    steady = tr.loc[tr["pressure"].idxmax()]
    p_ref, t_ref = float(steady["pressure"]), float(steady["temperature"])

    camera = R.shift_camera(
        R.framed_camera(field, width, height, pitch_deg=pitch_deg, margin=0.16),
        ndc_dy=0.08)
    if shell is None:
        shell = R.Shell(transmittance=0.0, thickness=5.6e-3, ambient=0.002,
                        diffuse=0.014, rim=0.045, rim_power=4.0, section=0.60,
                        glow=0.08, cut=True)
    geometry = R.build_geometry(field, camera, steps=steps, shell=shell)

    # Dwell on the start and the shutdown; skip across the steady plateau.
    def ramp(a, b, n):
        return np.linspace(a, b, n, endpoint=False)
    times = np.concatenate([ramp(0.0, 0.052, 150), ramp(0.052, 0.596, 24),
                            ramp(0.596, 0.634, 72), np.full(8, 0.634)])

    def state(t):
        p = float(np.interp(t, tr["time"], tr["pressure"]))
        temp = float(np.interp(t, tr["time"], tr["temperature"]))
        thrust = float(np.interp(t, tr["time"], tr["thrust"]))
        return p, temp, thrust

    reference = R.exposure_reference(
        geometry.paths @ R.emission_per_slice(field))

    font_path = _font_file("Liberation Sans")
    mono_path = _font_file("Liberation Mono")
    big = ImageFont.truetype(font_path, max(14, height // 30))
    small = ImageFont.truetype(mono_path, max(11, height // 46))

    writer = imageio.get_writer(out, fps=fps, codec="libx264", quality=8,
                                macro_block_size=8,
                                ffmpeg_params=["-pix_fmt", "yuv420p"])
    try:
        for t in times:
            p, temp, thrust = state(t)
            frame_field = field.scaled(field.temperature[0] * temp / t_ref, p / p_ref)
            img = R.render(geometry, frame_field, latitude=900.0, exposure=0.78,
                           reference=reference)
            pic = Image.fromarray(R.to_uint8(img))
            draw = ImageDraw.Draw(pic)
            pad = height // 24
            draw.text((pad, height - pad - big.size - small.size * 3.4),
                      f"t = {t * 1e3:6.1f} ms", font=big, fill="#f2f4f7")
            for i, line in enumerate((f"chamber   {p / 1e6:6.2f} MPa",
                                      f"          {temp:6.0f} K",
                                      f"thrust    {thrust / 1e3:6.1f} kN")):
                draw.text((pad, height - pad - small.size * (3.2 - i * 1.1)),
                          line, font=small, fill="#9aa4b2")
            writer.append_data(np.asarray(pic))
    finally:
        writer.close()
    return out


def _font_file(name: str) -> str:
    from matplotlib import font_manager
    return font_manager.findfont(font_manager.FontProperties(family=name))


def flow_animation(profile_csv: str, plume, out: str, *, width: int = 1920,
                   height: int = 1080, fps: int = 30, seconds: float = 10.0,
                   slowdown: float = 3000.0, max_step: float = 2.2e-3,
                   bins: int = 320,
                   tracer_count: int = 60000, gain: float = 1.0,
                   lifetime: float = 5.0e-4, pitch_deg: float = 13.0,
                   margin: float = 0.06, plume_reach: float = 3.2,
                   exposure: float = 0.45, latitude: float = 900.0,
                   shell=None, caption: str = "") -> str:
    """The engine running: steady emission, with the flow made visible.

    The field is steady, so the ray march is done once and every frame reuses
    it.  What moves is a population of tracers advected by the solver's own
    velocity -- see `ignis_viz.tracers` for exactly what they are and are not.

    Playback is `slowdown` times slower than reality: the gas crosses the whole
    picture in about a millisecond, which at real speed would be a third of a
    single frame.
    """
    import numpy as np
    from PIL import Image, ImageDraw, ImageFont
    import imageio.v2 as imageio
    import pandas as pd

    from . import render as R
    from . import tracers as T

    profile = pd.read_csv(profile_csv, comment="#")
    field = R.Field.from_profile(profile, bins=bins)
    length, r_exit = field.length, float(field.radius[-1])

    extra = None
    if plume is not None:
        extra = R.plume_hull(plume, length)
        extra = extra[extra[:, 0] <= length + plume_reach * r_exit]
    camera = R.framed_camera(field, width, height, pitch_deg=pitch_deg,
                             margin=margin, extra_points=extra)
    if shell is None:
        shell = R.Shell(transmittance=0.0, thickness=5.6e-3, ambient=0.002,
                        diffuse=0.014, rim=0.045, rim_power=4.0, section=0.60,
                        glow=0.08, cut=True)
    scene = R.Scene(field, shell, plume,
                    plume_reach=(plume_reach + 2.4) * r_exit)
    gas, wall, depth = R.march_scene(scene, camera, max_step=max_step)
    reference = R.exposure_reference(gas)

    dt = 1.0 / fps / slowdown
    tr = T.Tracers(field, plume, count=tracer_count, lifetime=lifetime)
    for _ in range(40):                       # let the population settle
        tr.advance(dt / 6.0)

    mono = ImageFont.truetype(_font_file("Liberation Mono"), max(11, height // 64))
    writer = imageio.get_writer(out, fps=fps, codec="libx264", quality=9,
                                macro_block_size=8,
                                ffmpeg_params=["-pix_fmt", "yuv420p"])
    try:
        for _ in range(int(round(seconds * fps))):
            streaks = T.streak(tr, camera, (height, width), dt, gain=gain,
                               depth=depth)
            frame = gas + streaks.reshape(-1, 3) * reference
            img = R.compose(frame, wall, (height, width), latitude=latitude,
                            exposure=exposure, reference=reference)
            pic = Image.fromarray(R.to_uint8(img))
            if caption:
                draw = ImageDraw.Draw(pic)
                pad = height // 26
                for i, line in enumerate(caption.split("\n")):
                    draw.text((pad, height - pad - mono.size * (len(caption.split("\n")) - i) * 1.5),
                              line, font=mono, fill="#9aa4b2")
            writer.append_data(np.asarray(pic))
    finally:
        writer.close()
    return out


def camera_move_animation(profile_csv: str, plume, out: str, *, width: int = 1280,
                          height: int = 720, fps: int = 24, frames: int = 168,
                          max_step: float = 2.2e-3, bins: int = 320,
                          yaw_from: float = 0.0, yaw_to: float = 46.0,
                          pitch_from: float = 8.0, pitch_to: float = 26.0,
                          reach_from: float = 1.1, reach_to: float = 3.4,
                          centre_from: float = 0.34, centre_to: float = 0.5,
                          tracer_count: int = 60000, gain: float = 1.0,
                          lifetime: float = 5.0e-4, slowdown: float = 3000.0,
                          exposure: float = 0.45, latitude: float = 900.0,
                          shell=None, progress=None,
                          transient_csv: str = "", ignition_fraction: float = 0.42,
                          ignition_end: float = 0.048) -> str:
    """A moving camera over the running engine.

    The scene is axisymmetric, so swinging the camera *around* the nozzle axis
    produces an identical picture frame after frame -- the one direction that
    buys nothing.  What changes the view is the angle between the line of sight
    and the axis, so the move swings downstream from a square side elevation to
    an oblique three-quarter, rising as it goes and pulling back from the
    throat to the whole plume.

    Every frame is a fresh march -- the geometry changes when the camera does,
    so none of it can be reused -- and the tracers keep advancing throughout, so
    the flow runs while the camera moves.
    """
    import numpy as np
    import imageio.v2 as imageio
    import pandas as pd

    from . import render as R
    from . import tracers as T

    profile = pd.read_csv(profile_csv, comment="#")
    field = R.Field.from_profile(profile, bins=bins)
    length, r_exit = field.length, float(field.radius[-1])
    if shell is None:
        shell = R.Shell(transmittance=0.0, thickness=5.6e-3, ambient=0.002,
                        diffuse=0.014, rim=0.045, rim_power=4.0, section=0.60,
                        glow=0.08, cut=True)

    hull = R.plume_hull(plume, length) if plume is not None else None

    # Optionally run the start-up transient while the camera moves.  The plume
    # is a steady solution at full chamber pressure, and during the ramp the
    # nozzle is so far over-expanded that Ignis predicts the flow separates
    # inside it -- so there is genuinely no attached plume until the chamber
    # comes up, and none is drawn until it does.
    ramp = None
    if transient_csv:
        tr_table = pd.read_csv(transient_csv, comment="#")
        steady = tr_table.loc[tr_table["pressure"].idxmax()]
        p_ref, t_ref = float(steady["pressure"]), float(steady["temperature"])
        attach = _attach_pressure(field, plume, p_ref)
        ramp = (tr_table, p_ref, t_ref, attach)

    dt = 1.0 / fps / slowdown
    tr = T.Tracers(field, plume, count=tracer_count, lifetime=lifetime)
    for _ in range(40):
        tr.advance(dt / 6.0)

    def ease(u):                       # smoothstep, so the move starts and ends at rest
        return u * u * (3.0 - 2.0 * u)

    def frame_camera(u):
        reach = reach_from + (reach_to - reach_from) * u
        extra = None
        if hull is not None:
            extra = hull[hull[:, 0] <= length + reach * r_exit]
        return R.framed_camera(
            field, width, height,
            yaw_deg=yaw_from + (yaw_to - yaw_from) * u,
            pitch_deg=pitch_from + (pitch_to - pitch_from) * u,
            margin=0.10, centre=centre_from + (centre_to - centre_from) * u,
            extra_points=extra), reach

    # Expose once, on the engine at full thrust, halfway through the move.  Set
    # from the opening frame instead, the exposure would be keyed to a chamber
    # that has not lit yet and everything after ignition would be white.
    mid_camera, mid_reach = frame_camera(0.5)
    reference = R.exposure_reference(R.march_scene(
        R.Scene(field, shell, plume, plume_reach=(mid_reach + 2.4) * r_exit),
        mid_camera, max_step=max_step)[0])

    writer = imageio.get_writer(out, fps=fps, codec="libx264", quality=9,
                                macro_block_size=8,
                                ffmpeg_params=["-pix_fmt", "yuv420p"])
    try:
        for k in range(frames):
            camera, reach = frame_camera(ease(k / max(frames - 1, 1)))
            # Frame on `reach` but march well past it: auto-framing leaves a
            # margin, and an oblique view sees further down the plume, so a
            # march that stopped at the framing limit would print the end of
            # the marched domain as a hard edge inside the picture.
            frame_field, frame_plume = field, plume
            if ramp is not None:
                tr_table, p_ref, t_ref, attach = ramp
                span = max(frames * ignition_fraction, 1.0)
                t_now = min(k / span, 1.0) * ignition_end
                p_now = float(np.interp(t_now, tr_table["time"], tr_table["pressure"]))
                temp = float(np.interp(t_now, tr_table["time"], tr_table["temperature"]))
                frame_field = field.scaled(field.temperature[0] * temp / t_ref,
                                           p_now / p_ref)
                frame_plume = plume if p_now >= attach else None
                tr.field = frame_field
            scene = R.Scene(frame_field, shell, frame_plume,
                            plume_reach=(reach + 2.4) * r_exit)
            gas, wall, depth = R.march_scene(scene, camera, max_step=max_step)
            gas = gas + T.streak(tr, camera, (height, width), dt, gain=gain,
                                 depth=depth,
                                 max_x=None if frame_plume is not None else length
                                 ).reshape(-1, 3) * reference
            img = R.compose(gas, wall, (height, width), latitude=latitude,
                            exposure=exposure, reference=reference)
            writer.append_data(R.to_uint8(img))
            if progress is not None:
                progress(k, frames)
    finally:
        writer.close()
    return out


def _attach_pressure(field, plume, p_steady: float) -> float:
    """Chamber pressure at which the nozzle stops being separated.

    Ignis's Summerfield criterion puts separation at an exit-to-ambient
    pressure ratio near 0.35.  The exit pressure of a choked nozzle is a fixed
    fraction of the chamber pressure, so that maps straight onto a chamber
    pressure -- above it the nozzle flows full and the computed plume applies,
    below it the flow is separated and no plume is drawn.
    """
    if plume is None:
        return float("inf")
    exit_fraction = float(plume.exit.pressure) / p_steady
    return 0.35 * float(plume.ambient_pressure) / max(exit_fraction, 1e-12)
