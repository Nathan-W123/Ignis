"""Axisymmetric Euler solution for the exhaust plume.

WHY THIS EXISTS
---------------
Ignis solves the nozzle and stops at the exit plane, which is the right place
for a quasi-1D internal-flow model to stop: everything downstream is a genuinely
two-dimensional external flow.  But the plume is the part of a rocket engine
everyone recognises, and its shock-cell structure is the visible consequence of
the exit state the solver computed.  Drawing a decorative plume onto a render
would misrepresent the project, so this module computes one instead.

WHAT IT SOLVES
--------------
The axisymmetric Euler equations,

    d(U)/dt + d(F)/dx + d(G)/dr = -(1/r) (rho v, rho u v, rho v^2, v(E+p), rho v Y),

for U = (rho, rho u, rho v, E, rho Y), by a second-order finite-volume method:
MUSCL reconstruction with a minmod limiter on primitive variables, an HLLC
approximate Riemann solver at the faces, and a two-stage SSP Runge-Kutta step.
The jet enters through the exit plane at the state Ignis computed and expands
into still ambient air.  `Y` is a passive scalar, 1 in jet gas and 0 in ambient
air, carried so that thermal emission is never attributed to the surrounding
atmosphere.

This resolves what sets the plume's appearance: the Prandtl-Meyer fan or lip
shock, the barrel shock, the Mach disc where one forms, and the shock-cell
train downstream.

WHAT IT DOES NOT SOLVE
----------------------
This is an inviscid model and it is NOT part of the validated Ignis core:

  * No viscosity and no turbulence model, so the jet shear layer spreads only
    by numerical diffusion and the shock cells persist much further downstream
    than they would in a real plume.
  * One ratio of specific heats for the whole domain, frozen at the nozzle exit
    value.  The ambient air really has gamma = 1.4; the jet core structure is
    set by the jet's own gamma and by p_e/p_a, and the ambient contributes
    nothing to the emission, so a single gamma is used and declared here rather
    than a second species equation being carried.
  * Frozen chemistry.  The products leave the nozzle fuel-rich and a real plume
    at low altitude afterburns with entrained air, which is bright.  None of
    that is modelled, so the plume computed here is dimmer and cooler
    downstream than a real one.
  * No radiation transport; emission is treated as optically thin at render
    time, as in the nozzle.

The one quantitative check available is the Mach-disc standoff distance of an
underexpanded jet, which has a well-established correlation
(Ashkenas & Sherman 1966): x_M/D_e = 0.67 sqrt(p_e/p_a).  See
`validate_mach_disc` and `tools/validate_plume.py`.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Tuple

import numpy as np

R_UNIVERSAL = 8.31446261815324
M_AIR = 0.0289647          # kg/mol, dry air


@dataclass
class ExitState:
    """The nozzle exit plane, as Ignis computed it."""

    radius: float           # m
    pressure: float         # Pa
    temperature: float      # K
    density: float          # kg/m^3
    velocity: float         # m/s, axial
    mach: float
    gamma: float
    molar_mass: float       # kg/mol

    @classmethod
    def from_profile(cls, profile) -> "ExitState":
        e = profile.iloc[-1]
        return cls(radius=float(e["radius"]), pressure=float(e["pressure"]),
                   temperature=float(e["temperature"]), density=float(e["density"]),
                   velocity=float(e["velocity"]), mach=float(e["mach"]),
                   gamma=float(e["gamma_s"]), molar_mass=float(e["molar_mass"]))


@dataclass
class PlumeField:
    """A converged plume solution on the (x, r) grid, cell centres."""

    x: np.ndarray            # m, (nx,), measured from the exit plane
    r: np.ndarray            # m, (nr,)
    density: np.ndarray      # kg/m^3, (nr, nx)
    velocity_x: np.ndarray
    velocity_r: np.ndarray
    pressure: np.ndarray
    temperature: np.ndarray
    mach: np.ndarray
    jet_fraction: np.ndarray # 1 in jet gas, 0 in ambient air
    gamma: float
    residual_history: np.ndarray
    ambient_pressure: float
    exit: ExitState


def _minmod(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    return np.where(a * b <= 0.0, 0.0, np.where(np.abs(a) < np.abs(b), a, b))


def _slope(w: np.ndarray, axis: int) -> np.ndarray:
    """Minmod-limited slopes; zero in the outermost ghost layer."""
    lo = [slice(None)] * w.ndim
    mid = [slice(None)] * w.ndim
    hi = [slice(None)] * w.ndim
    lo[axis], mid[axis], hi[axis] = slice(None, -2), slice(1, -1), slice(2, None)
    s = np.zeros_like(w)
    s[tuple(mid)] = _minmod(w[tuple(mid)] - w[tuple(lo)], w[tuple(hi)] - w[tuple(mid)])
    return s


def _faces(w: np.ndarray, axis: int) -> Tuple[np.ndarray, np.ndarray]:
    """Left and right reconstructed states on every interior face."""
    s = _slope(w, axis)
    a = [slice(None)] * w.ndim
    b = [slice(None)] * w.ndim
    a[axis], b[axis] = slice(1, -2), slice(2, -1)
    return (w[tuple(a)] + 0.5 * s[tuple(a)], w[tuple(b)] - 0.5 * s[tuple(b)])


def _hllc(rl, unl, utl, pl, yl, rr, unr, utr, pr, yr, g):
    """HLLC flux in the face-normal direction (Toro, 3rd ed., ch. 10).

    Returns (mass, normal momentum, tangential momentum, energy, scalar).
    """
    rl = np.maximum(rl, 1e-12)
    rr = np.maximum(rr, 1e-12)
    pl = np.maximum(pl, 1e-9)
    pr = np.maximum(pr, 1e-9)
    al = np.sqrt(g * pl / rl)
    ar = np.sqrt(g * pr / rr)
    el = pl / (g - 1.0) + 0.5 * rl * (unl * unl + utl * utl)
    er = pr / (g - 1.0) + 0.5 * rr * (unr * unr + utr * utr)

    # Pressure-based (PVRS) wave-speed estimate with the two-rarefaction guard.
    rbar = 0.5 * (rl + rr)
    abar = 0.5 * (al + ar)
    pstar = np.maximum(0.0, 0.5 * (pl + pr) - 0.5 * (unr - unl) * rbar * abar)
    ql = np.where(pstar > pl, np.sqrt(1.0 + (g + 1.0) / (2.0 * g) * (pstar / pl - 1.0)), 1.0)
    qr = np.where(pstar > pr, np.sqrt(1.0 + (g + 1.0) / (2.0 * g) * (pstar / pr - 1.0)), 1.0)
    sl = unl - al * ql
    sr = unr + ar * qr
    denom = rl * (sl - unl) - rr * (sr - unr)
    sm = (pr - pl + rl * unl * (sl - unl) - rr * unr * (sr - unr)) / np.where(
        np.abs(denom) < 1e-30, 1e-30, denom)

    def flux(r_, un, ut, p_, e_, y_):
        return np.stack([r_ * un, r_ * un * un + p_, r_ * un * ut,
                         un * (e_ + p_), r_ * un * y_])

    fl = flux(rl, unl, utl, pl, el, yl)
    fr = flux(rr, unr, utr, pr, er, yr)

    def star(r_, un, ut, p_, e_, y_, s_):
        coef = r_ * (s_ - un) / np.where(np.abs(s_ - sm) < 1e-30, 1e-30, s_ - sm)
        u_star = np.stack([
            coef,
            coef * sm,
            coef * ut,
            coef * (e_ / r_ + (sm - un) * (sm + p_ / (r_ * np.where(
                np.abs(s_ - un) < 1e-30, 1e-30, s_ - un)))),
            coef * y_,
        ])
        u_old = np.stack([r_, r_ * un, r_ * ut, e_, r_ * y_])
        return u_star, u_old

    usl, uol = star(rl, unl, utl, pl, el, yl, sl)
    usr, uor = star(rr, unr, utr, pr, er, yr, sr)
    fsl = fl + sl * (usl - uol)
    fsr = fr + sr * (usr - uor)

    return np.where(sl >= 0.0, fl,
                    np.where(sm >= 0.0, fsl,
                             np.where(sr >= 0.0, fsr, fr)))


def _apply_bc(w: np.ndarray, r_centres: np.ndarray, inflow: np.ndarray,
              jet_mask: np.ndarray) -> None:
    """Ghost cells, in place.  `w` is (5, nr+4, nx+4) primitives."""
    # Left boundary, x = 0.  Inside the exit radius the jet enters at a fixed
    # supersonic state, so every characteristic points inwards and the whole
    # state is prescribed.  Outside it is the nozzle's own exit face: a solid
    # wall, mirrored with the axial velocity reversed.
    for ghost, source in ((1, 2), (0, 3)):
        w[:, :, ghost] = np.where(jet_mask[None, :], inflow[:, None], w[:, :, source])
        w[1, :, ghost] = np.where(jet_mask, inflow[1], -w[1, :, source])

    # Right boundary: outflow.  The core is supersonic; zero-gradient is the
    # standard choice and the domain is long enough that the subsonic pockets
    # behind Mach discs are well upstream of it.
    w[:, :, -2] = w[:, :, -3]
    w[:, :, -1] = w[:, :, -3]

    # Axis of symmetry, r = 0: mirror with the radial velocity reversed.
    for ghost, source in ((1, 2), (0, 3)):
        w[:, ghost, :] = w[:, source, :]
        w[2, ghost, :] = -w[2, source, :]

    # Far field: zero-gradient.
    w[:, -2, :] = w[:, -3, :]
    w[:, -1, :] = w[:, -3, :]


def _to_conserved(w: np.ndarray, g: float) -> np.ndarray:
    rho, u, v, p, y = w
    e = p / (g - 1.0) + 0.5 * rho * (u * u + v * v)
    return np.stack([rho, rho * u, rho * v, e, rho * y])


def _to_primitive(u: np.ndarray, g: float) -> np.ndarray:
    rho = np.maximum(u[0], 1e-12)
    vx, vr = u[1] / rho, u[2] / rho
    p = np.maximum((g - 1.0) * (u[3] - 0.5 * rho * (vx * vx + vr * vr)), 1e-9)
    return np.stack([rho, vx, vr, p, np.clip(u[4] / rho, 0.0, 1.0)])


def _residual(w: np.ndarray, g: float, dx: float, dr: float,
              r_centres: np.ndarray) -> np.ndarray:
    """-(dF/dx + dG/dr) - source, on interior cells."""
    # Axial faces.
    wl, wr = _faces(w, axis=2)
    fx = _hllc(wl[0], wl[1], wl[2], wl[3], wl[4],
               wr[0], wr[1], wr[2], wr[3], wr[4], g)
    # Radial faces: the normal velocity is v, the tangential one u, so the two
    # momentum components swap places on the way in and on the way out.
    wb, wt = _faces(w, axis=1)
    fr = _hllc(wb[0], wb[2], wb[1], wb[3], wb[4],
               wt[0], wt[2], wt[1], wt[3], wt[4], g)
    fr = np.stack([fr[0], fr[2], fr[1], fr[3], fr[4]])

    div = ((fx[:, 2:-2, 1:] - fx[:, 2:-2, :-1]) / dx
           + (fr[:, 1:, 2:-2] - fr[:, :-1, 2:-2]) / dr)

    rho, u, v, p, y = w[:, 2:-2, 2:-2]
    e = p / (g - 1.0) + 0.5 * rho * (u * u + v * v)
    src = np.stack([rho * v, rho * u * v, rho * v * v, v * (e + p), rho * v * y])
    return -div - src / r_centres[None, :, None]


def solve(exit_state: ExitState, ambient_pressure: float, *,
          ambient_temperature: float = 288.15, nx: int = 720, nr: int = 280,
          length: float = 14.0, width: float = 3.2, cfl: float = 0.35,
          max_steps: int = 40000, tolerance: float = 2.0e-4,
          report: Optional[callable] = None) -> PlumeField:
    """March the plume to steady state.

    `length` and `width` are in exit radii.  The gamma used throughout is the
    nozzle exit value -- see the module docstring for why, and for everything
    this model leaves out.
    """
    g = exit_state.gamma
    r_jet = R_UNIVERSAL / exit_state.molar_mass
    r_air = R_UNIVERSAL / M_AIR

    x_max = length * exit_state.radius
    r_max = width * exit_state.radius
    dx, dr = x_max / nx, r_max / nr
    x_centres = (np.arange(nx) + 0.5) * dx
    r_centres = (np.arange(nr) + 0.5) * dr

    # Padded with two ghost layers on every side.
    r_pad = (np.arange(-2, nr + 2) + 0.5) * dr
    jet_mask = np.abs(r_pad) < exit_state.radius

    inflow = np.array([exit_state.density, exit_state.velocity, 0.0,
                       exit_state.pressure, 1.0])
    ambient = np.array([ambient_pressure / (r_air * ambient_temperature),
                        0.0, 0.0, ambient_pressure, 0.0])

    # Start with the jet column already established rather than with the whole
    # domain at rest: the same steady state, reached without first paying for
    # the jet to blast a cavity through still air.
    w = np.repeat(np.repeat(ambient[:, None, None], nr + 4, axis=1), nx + 4, axis=2)
    w[:, jet_mask, :] = inflow[:, None, None]
    _apply_bc(w, r_pad, inflow, jet_mask)
    u = _to_conserved(w, g)

    history = []
    peak = 0.0
    for step in range(max_steps):
        w = _to_primitive(u, g)
        _apply_bc(w, r_pad, inflow, jet_mask)
        a = np.sqrt(g * w[3] / w[0])
        dt = cfl / np.max((np.abs(w[1]) + a) / dx + (np.abs(w[2]) + a) / dr)

        k1 = _residual(w, g, dx, dr, r_centres)
        u1 = u.copy()
        u1[:, 2:-2, 2:-2] = u[:, 2:-2, 2:-2] + dt * k1

        w1 = _to_primitive(u1, g)
        _apply_bc(w1, r_pad, inflow, jet_mask)
        k2 = _residual(w1, g, dx, dr, r_centres)
        u[:, 2:-2, 2:-2] = 0.5 * (u[:, 2:-2, 2:-2] + u1[:, 2:-2, 2:-2] + dt * k2)

        # L2 of the density residual against the largest value seen so far, so
        # the criterion is "the flow has stopped changing relative to how hard
        # it was changing at its worst" and does not depend on the units.
        res = float(np.sqrt(np.mean(k1[0] ** 2)))
        peak = max(peak, res)
        history.append(res / peak if peak > 0.0 else 0.0)
        if report is not None and step % 500 == 0:
            report(step, history[-1], dt)
        if step > 400 and history[-1] < tolerance:
            break

    w = _to_primitive(u, g)
    _apply_bc(w, r_pad, inflow, jet_mask)
    rho, vx, vr, p, y = w[:, 2:-2, 2:-2]
    r_mix = y * r_jet + (1.0 - y) * r_air
    temperature = p / (rho * r_mix)
    mach = np.sqrt(vx * vx + vr * vr) / np.sqrt(g * p / rho)
    return PlumeField(x=x_centres, r=r_centres, density=rho, velocity_x=vx,
                      velocity_r=vr, pressure=p, temperature=temperature,
                      mach=mach, jet_fraction=y, gamma=g,
                      residual_history=np.asarray(history),
                      ambient_pressure=ambient_pressure, exit=exit_state)


def mach_disc_location(field: PlumeField) -> Optional[float]:
    """Axial position of the first Mach disc, or None if the axis stays supersonic.

    The Mach disc is a normal shock standing on the axis, so it shows up as the
    first downstream crossing of M = 1 in the axis row.
    """
    m = field.mach[0]
    below = np.flatnonzero(m < 1.0)
    if below.size == 0:
        return None
    i = int(below[0])
    if i == 0:
        return float(field.x[0])
    # Linear interpolation onto M = 1 between the bracketing cells.
    m0, m1 = m[i - 1], m[i]
    t = (m0 - 1.0) / (m0 - m1)
    return float(field.x[i - 1] + t * (field.x[i] - field.x[i - 1]))


def lip_shock_theory(exit_state: ExitState, ambient_pressure: float):
    """Exact oblique-shock solution at the lip of an over-expanded jet.

    An over-expanded jet leaves the nozzle below ambient pressure, so it must
    be compressed back up to it.  That happens through an oblique shock
    attached to the nozzle lip, and its wave angle is not a correlation: given
    the upstream Mach number and the required pressure ratio, the normal-shock
    relation fixes the normal Mach number,

        p2/p1 = 1 + 2 gamma / (gamma + 1) (Mn1^2 - 1),

    the wave angle follows from Mn1 = M1 sin(beta), and the flow deflection
    from the theta-beta-M relation.  Comparing that angle against the one the
    solver actually produces is a check with an exact answer, which a
    correlation would not give.

    Returns (beta, theta, mach_normal) in radians, or None if the jet is not
    over-expanded or the shock cannot attach.
    """
    g, m1 = exit_state.gamma, exit_state.mach
    ratio = ambient_pressure / exit_state.pressure
    if ratio <= 1.0:
        return None
    mn1 = np.sqrt((ratio - 1.0) * (g + 1.0) / (2.0 * g) + 1.0)
    if mn1 >= m1:
        return None
    beta = np.arcsin(mn1 / m1)
    theta = np.arctan(2.0 / np.tan(beta) * (m1 * m1 * np.sin(beta) ** 2 - 1.0)
                      / (m1 * m1 * (g + np.cos(2.0 * beta)) + 2.0))
    return float(beta), float(theta), float(mn1)


def measure_lip_shock(field: PlumeField, r_lo: float = 0.45,
                      r_hi: float = 0.92) -> Optional[float]:
    """Wave angle of the lip shock in the computed field, in radians.

    The shock is located at each radial station as the first axial crossing of
    the midpoint between the exit and ambient pressures, over a band of radii
    clear of both the lip singularity and the axis.  A straight-line fit
    through those points gives the angle.
    """
    p_mid = 0.5 * (field.exit.pressure + field.ambient_pressure)
    r_e = field.exit.radius
    rows = np.flatnonzero((field.r > r_lo * r_e) & (field.r < r_hi * r_e))
    xs, rs = [], []
    for j in rows:
        above = np.flatnonzero(field.pressure[j] > p_mid)
        if above.size == 0:
            continue
        i = int(above[0])
        if i == 0:
            continue
        p0, p1 = field.pressure[j, i - 1], field.pressure[j, i]
        t = (p_mid - p0) / (p1 - p0)
        xs.append(field.x[i - 1] + t * (field.x[i] - field.x[i - 1]))
        rs.append(field.r[j])
    if len(xs) < 4:
        return None
    # r decreases as x increases along the shock; the wave angle is measured
    # from the upstream flow direction, which is the axis.
    slope = np.polyfit(np.asarray(xs), np.asarray(rs), 1)[0]   # dr/dx along the shock
    return float(np.arctan(abs(slope)))
