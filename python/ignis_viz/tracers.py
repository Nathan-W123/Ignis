"""Massless tracers, advected by the computed velocity field.

A steady solution has nothing moving in it, and a render of one is a still life
however bright it is.  What actually conveys a rocket engine is the flow: gas
crawling through the chamber, accelerating hard through the throat, and leaving
at three kilometres a second.

So this module carries tracer particles through the solved field.  They are
markers, not physics -- they have no mass, no volume and no influence on
anything -- and every one of them moves at the local velocity the solver
computed: the quasi-1D axial velocity inside the nozzle, and the two-dimensional
velocity of the Euler plume outside it.  Inside the nozzle a tracer keeps its
fractional radius r/R(x), which is the streamline of a quasi-1D flow.

Because each particle moves at the local speed, the streak it leaves in a frame
is proportional to that speed: the chamber barely smears while the exhaust draws
long lines.  That contrast is the solver's, not a stylistic choice.

Where the tracers are is a visualisation choice and is not.  Injected at a
steady rate from the injector face they would settle at a number density
proportional to the gas density, which falls by a factor of eighty-five between
the chamber and the exit plane, leaving the exhaust -- the part worth looking at
-- almost empty.  Instead each tracer is given a finite life and respawns
somewhere in the flow, which spreads them evenly over the picture.  Their
density therefore means nothing; their direction, their length and their colour
are the solution.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Tuple

import numpy as np
from scipy.ndimage import gaussian_filter

from .blackbody import blackbody_rgb
from .render import Camera, Field, pixel_coords


@dataclass
class Tracers:
    """A population of tracers moving through the engine and its plume."""

    field: Field
    plume: object = None
    count: int = 24000
    lifetime: float = 4.0e-4        # s, mean tracer life before it respawns
    seed: int = 20260917

    def __post_init__(self) -> None:
        self._rng = np.random.default_rng(self.seed)
        self.length = self.field.length
        if self.plume is not None:
            self._p_dx = float(self.plume.x[1] - self.plume.x[0])
            self._p_dr = float(self.plume.r[1] - self.plume.r[0])
            self._p_nr, self._p_nx = self.plume.temperature.shape
            self.reach = self.length + self._p_nx * self._p_dx
            self._p_rmax = self._p_nr * self._p_dr
        else:
            self.reach = self.length
            self._p_rmax = 0.0
        # Visible boundary of the flow against x, for seeding.
        self._edge_x = self.field.x.copy()
        self._edge_r = self.field.radius.copy()
        if self.plume is not None:
            jet = self.plume.jet_fraction
            r_edge = np.array([(self.plume.r[np.flatnonzero(c >= 0.5)[-1]]
                                if np.any(c >= 0.5) else 0.0) for c in jet.T])
            self._edge_x = np.concatenate([self._edge_x, self.plume.x + self.length])
            self._edge_r = np.concatenate([self._edge_r, r_edge])

        self.x = np.zeros(self.count)
        self.r = np.zeros(self.count)
        self.theta = np.zeros(self.count)
        self.age = np.zeros(self.count)
        self.life = np.ones(self.count)
        self._seed(np.ones(self.count, dtype=bool))
        # Stagger the ages so the whole population does not blink together.
        self.age = self._rng.uniform(0.0, 1.0, self.count) * self.life

    # -- seeding ----------------------------------------------------------
    def _wall(self, x: np.ndarray) -> np.ndarray:
        return np.interp(x, self.field.x, self.field.radius)

    def _edge(self, x: np.ndarray) -> np.ndarray:
        return np.interp(x, self._edge_x, self._edge_r)

    def _seed(self, mask: np.ndarray) -> None:
        """Respawn tracers anywhere in the flow, with a fresh lifetime."""
        n = int(mask.sum())
        if n == 0:
            return
        x = self._rng.uniform(0.0, self.reach * 0.995, n)
        self.x[mask] = x
        # Uniform over the cross-section, not over the radius.
        self.r[mask] = self._edge(x) * np.sqrt(self._rng.uniform(0.0, 0.97, n))
        self.theta[mask] = self._rng.uniform(0.0, 2.0 * np.pi, n)
        self.age[mask] = 0.0
        self.life[mask] = self.lifetime * self._rng.uniform(0.55, 1.45, n)

    # -- advection --------------------------------------------------------
    def _plume_sample(self, arrays, x: np.ndarray, r: np.ndarray):
        """Bilinear sample of plume arrays at (x measured from the exit, r)."""
        fx = np.clip((x - self.plume.x[0]) / self._p_dx, 0.0, self._p_nx - 1.001)
        fr = np.clip((r - self.plume.r[0]) / self._p_dr, 0.0, self._p_nr - 1.001)
        i0, j0 = fx.astype(np.int32), fr.astype(np.int32)
        tx, tr = fx - i0, fr - j0
        out = []
        for a in arrays:
            out.append((a[j0, i0] * (1 - tx) + a[j0, i0 + 1] * tx) * (1 - tr)
                       + (a[j0 + 1, i0] * (1 - tx) + a[j0 + 1, i0 + 1] * tx) * tr)
        return out

    def advance(self, dt: float) -> None:
        """Move every tracer by one step, re-seeding the ones that leave."""
        inside = self.x < self.length
        if inside.any():
            x = self.x[inside]
            frac = self.r[inside] / np.maximum(self._wall(x), 1e-9)
            u = np.interp(x, self.field.x, self.field.velocity)
            x_new = x + u * dt
            self.x[inside] = x_new
            # A quasi-1D streamline keeps its fraction of the local radius.
            self.r[inside] = frac * self._wall(np.minimum(x_new, self.length))

        out = ~inside
        if out.any() and self.plume is not None:
            xp = self.x[out] - self.length
            u, v = self._plume_sample(
                (self.plume.velocity_x, self.plume.velocity_r), xp, self.r[out])
            self.x[out] = self.x[out] + u * dt
            self.r[out] = np.abs(self.r[out] + v * dt)

        self.age += dt
        if self.plume is None:
            gone = self.x >= self.length
        else:
            gone = (self.x >= self.reach) | (self.r >= self._p_rmax)
        gone |= self.age >= self.life
        if self.plume is not None:
            beyond = self.x > self.length
            if beyond.any():
                y, = self._plume_sample((self.plume.jet_fraction,),
                                        self.x[beyond] - self.length, self.r[beyond])
                left_jet = np.zeros_like(gone)
                left_jet[beyond] = y < 0.35
                gone |= left_jet
        self._seed(gone)

    # -- rendering --------------------------------------------------------
    def positions(self) -> np.ndarray:
        return np.stack([self.x, self.r * np.cos(self.theta),
                         self.r * np.sin(self.theta)], axis=-1)

    def _local(self) -> Tuple[np.ndarray, np.ndarray]:
        """Local temperature and density at every tracer."""
        t = np.interp(self.x, self.field.x, self.field.temperature)
        rho = np.interp(self.x, self.field.x, self.field.density)
        if self.plume is not None:
            beyond = self.x > self.length
            if beyond.any():
                tp, rp = self._plume_sample(
                    (self.plume.temperature, self.plume.density),
                    self.x[beyond] - self.length, self.r[beyond])
                t[beyond], rho[beyond] = tp, rp
        return t, rho

    def colours(self) -> np.ndarray:
        """Per-tracer colour: blackbody hue, weighted to stay visible throughout.

        The gas spans more than two decades of brightness between the chamber
        and the exhaust, so a tracer of fixed brightness would be lost in one
        and blinding in the other.  Weighting each by the local emission keeps
        every tracer at a similar ratio to the gas around it, which is what
        makes them legible over the whole field.  A smooth fade over the
        tracer's life keeps them from appearing and vanishing abruptly.
        """
        t, rho = self._local()
        emission = rho * np.clip(t, 1.0, None) ** 4
        weight = emission / max(float(emission.max()), 1e-30)
        fade = np.sin(np.pi * np.clip(self.age / np.maximum(self.life, 1e-12), 0.0, 1.0))
        return (blackbody_rgb(np.clip(t, 300.0, None))
                * (weight * fade)[:, None])


def streak(tracers: Tracers, camera: Camera, shape: Tuple[int, int], dt: float, *,
           substeps: int = 6, sigma: float = 0.9, gain: float = 1.0,
           depth: Optional[np.ndarray] = None,
           max_x: Optional[float] = None) -> np.ndarray:
    """Advance the tracers one frame and splat the path each one swept.

    The frame is built from `substeps` intermediate positions so a fast tracer
    draws a streak rather than a dot; the streak length is then literally
    speed times exposure, as it would be on a camera.

    `depth` is the distance to the nearest opaque surface on each ray, from the
    renderer.  Splatting is a screen-space operation and knows nothing about the
    scene, so without it a tracer inside the chamber draws straight through the
    engine wall in front of it.
    """
    h, w = shape
    buffer = np.zeros((h, w, 3), dtype=np.float64)
    sub = dt / substeps
    for k in range(substeps):
        tracers.advance(sub)
        points = tracers.positions()
        px = pixel_coords(camera, points)
        col = tracers.colours() * (gain * (k + 1) / substeps)
        ix, iy = px[:, 0].astype(np.int64), px[:, 1].astype(np.int64)
        ok = (ix >= 0) & (ix < w) & (iy >= 0) & (iy < h)
        if max_x is not None:
            # Nothing is drawn downstream of `max_x`, for the frames of a start
            # where the nozzle is still separated and there is no plume behind
            # the tracers to belong to.
            ok &= tracers.x <= max_x
        if depth is not None:
            # The ray through a tracer's pixel passes through the tracer, so its
            # distance from the eye is directly comparable with the wall depth.
            flat = np.clip(iy, 0, h - 1) * w + np.clip(ix, 0, w - 1)
            ok &= np.linalg.norm(points - camera.eye, axis=1) < depth[flat]
        np.add.at(buffer, (iy[ok], ix[ok]), col[ok])
    return gaussian_filter(buffer, sigma=(sigma, sigma, 0))
