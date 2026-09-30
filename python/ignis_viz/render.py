"""Volumetric render of the Ignis flow field.

WHAT THIS IS
------------
An emission-only volume render of the axisymmetric quasi-1D solution.  The
solver gives temperature, density and velocity as functions of axial position;
revolving that about the axis gives a three-dimensional field, and the image is
the line integral of its thermal emission along each camera ray.

For an optically thin gray emitter the emission coefficient is

    j = kappa rho B(T),

so radiance along a ray is the integral of `rho * sigma T^4` weighted by the
blackbody chromaticity at the local temperature.  That is what is integrated
here: the hue comes from Planck's law at the computed temperature (see
blackbody.py) and the brightness from the computed density and temperature.
The grey absorption coefficient is unknown, so the absolute scale is arbitrary
and the image is normalised for display -- the *relative* structure is the
solver's.

WHAT THIS IS NOT
----------------
It is not a radiative-transfer solution: there is no absorption, no scattering
and no band structure, and the products' real emissivity is dominated by the
H2O and CO2 bands rather than being gray.  It is not a 2-D or 3-D flow
solution either: the quasi-1D field is uniform across each cross-section, and
revolving it says exactly that and no more.

PERFORMANCE
-----------
The geometry is fixed while the flow changes, so the expensive part -- ray
marching to find how much path each pixel spends in each axial slice -- is done
once and reused.  A frame is then a matrix multiply of that path matrix against
the per-slice emission, which makes an animation cheap.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Tuple

import numpy as np
from scipy.ndimage import gaussian_filter

from .blackbody import blackbody_rgb

SIGMA_SB = 5.670374419e-8


@dataclass
class Field:
    """The revolved quasi-1D solution, sampled onto a uniform axial grid."""

    x: np.ndarray           # m, bin centres
    radius: np.ndarray      # m, wall radius at each bin
    temperature: np.ndarray # K
    density: np.ndarray     # kg/m^3
    mach: np.ndarray
    velocity: Optional[np.ndarray] = None   # m/s, axial

    @property
    def length(self) -> float:
        return float(self.x[-1] + 0.5 * (self.x[1] - self.x[0]))

    @property
    def max_radius(self) -> float:
        return float(self.radius.max())

    @classmethod
    def from_profile(cls, profile, bins: int = 220) -> "Field":
        """Resample a solver profile onto `bins` uniform axial slices."""
        x = np.asarray(profile["x"], dtype=np.float64)
        edges = np.linspace(x.min(), x.max(), bins + 1)
        centres = 0.5 * (edges[:-1] + edges[1:])

        def grid(name):
            return np.interp(centres, x, np.asarray(profile[name], dtype=np.float64))

        return cls(x=centres, radius=grid("radius"), temperature=grid("temperature"),
                   density=grid("density"), mach=grid("mach"),
                   velocity=grid("velocity") if "velocity" in profile else None)

    def scaled(self, t_chamber: float, p_ratio: float) -> "Field":
        """The same normalised profile at a different chamber state.

        Quasi-steady scaling, the same approximation the transient's own thrust
        estimate makes: the shape of the normalised temperature and density
        distributions is held fixed and the level follows the chamber.
        """
        ratio = t_chamber / self.temperature[0]
        t = self.temperature * ratio
        rho = self.density * p_ratio / max(ratio, 1e-9)
        # At a fixed Mach distribution the speed follows the local sound speed,
        # so it scales with the square root of the temperature.
        v = None if self.velocity is None else self.velocity * np.sqrt(max(ratio, 0.0))
        return Field(self.x, self.radius, t, rho, self.mach, v)


@dataclass
class Camera:
    eye: np.ndarray
    target: np.ndarray
    up: np.ndarray
    fov_deg: float
    width: int
    height: int

    def rays(self) -> Tuple[np.ndarray, np.ndarray]:
        """Origin and unit direction for every pixel, row-major."""
        forward = self.target - self.eye
        forward = forward / np.linalg.norm(forward)
        right = np.cross(forward, self.up)
        right = right / np.linalg.norm(right)
        up = np.cross(right, forward)

        aspect = self.width / self.height
        half_h = np.tan(np.radians(self.fov_deg) * 0.5)
        half_w = half_h * aspect
        # Pixel centres, y downwards in image space.
        px = (np.arange(self.width) + 0.5) / self.width * 2.0 - 1.0
        py = 1.0 - (np.arange(self.height) + 0.5) / self.height * 2.0
        gx, gy = np.meshgrid(px, py)
        d = (forward[None, None, :]
             + (gx * half_w)[..., None] * right[None, None, :]
             + (gy * half_h)[..., None] * up[None, None, :])
        d /= np.linalg.norm(d, axis=-1, keepdims=True)
        return self.eye.astype(np.float32), d.reshape(-1, 3).astype(np.float32)

    def basis(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Forward, right and up of the view frame."""
        forward = self.target - self.eye
        forward = forward / np.linalg.norm(forward)
        right = np.cross(forward, self.up)
        right = right / np.linalg.norm(right)
        return forward, right, np.cross(right, forward)

    def project(self, points: np.ndarray) -> np.ndarray:
        """Normalised device coordinates of world points; (N, 2), +-1 at the edges."""
        forward, right, up = self.basis()
        half_h = np.tan(np.radians(self.fov_deg) * 0.5)
        half_w = half_h * self.width / self.height
        v = np.asarray(points, dtype=np.float64) - self.eye
        depth = np.maximum(v @ forward, 1e-9)
        return np.stack([(v @ right) / depth / half_w,
                         (v @ up) / depth / half_h], axis=-1)


def default_camera(field: Field, width: int, height: int) -> Camera:
    """A three-quarter view that shows the throat and looks into the bell."""
    L = field.length
    return Camera(eye=np.array([L * 0.30, L * 0.30, L * 1.28]),
                  target=np.array([L * 0.56, L * 0.005, 0.0]),
                  up=np.array([0.0, 1.0, 0.0]),
                  fov_deg=34.0, width=width, height=height)


def aft_camera(field: Field, width: int, height: int, *,
               azimuth_deg: float = 0.0, elevation_deg: float = 18.0,
               distance: float = 1.85, fov_deg: float = 30.0) -> Camera:
    """Behind the exit plane, looking forward up the bore.

    This is the view that shows what the solver actually computed.  The bell is
    a dark silhouette; inside it the gas column recedes towards the throat,
    which is both the narrowest and by far the brightest part of the field, so
    the image has a bright structured core against black -- rather than the
    side elevation, where a bell nozzle reads as a featureless cone.

    `azimuth_deg` swings the eye around the axis, which is what an orbit
    animation varies; everything else is held.
    """
    L = field.length
    a = np.radians(azimuth_deg)
    e = np.radians(elevation_deg)
    target = np.array([L * 0.30, 0.0, 0.0])
    # `elevation_deg` is the angle off the nozzle axis; `azimuth_deg` rolls the
    # eye around it.  +x is downstream, so the eye sits aft of the exit plane.
    direction = np.array([np.cos(e), np.sin(e) * np.cos(a), np.sin(e) * np.sin(a)])
    return Camera(eye=target + distance * L * direction, target=target,
                  up=np.array([0.0, 1.0, 0.0]),
                  fov_deg=fov_deg, width=width, height=height)


def _hull(field: Field, n_theta: int = 24) -> np.ndarray:
    """Points on the outer surface of revolution, for framing tests."""
    theta = np.linspace(0.0, 2.0 * np.pi, n_theta, endpoint=False)
    r = field.radius[:, None]
    return np.stack([np.repeat(field.x[:, None], n_theta, axis=1).ravel(),
                     (r * np.cos(theta)[None, :]).ravel(),
                     (r * np.sin(theta)[None, :]).ravel()], axis=-1)


def plume_hull(plume, axial_offset: float, n_theta: int = 24,
               threshold: float = 0.5) -> np.ndarray:
    """Points on the jet boundary of a plume, for framing.

    The boundary is where the passive jet-fraction scalar falls through
    `threshold`, which is the visible edge of the plume -- framing on the whole
    computational domain instead would include a lot of dark ambient air.
    """
    r_edge = np.array([
        (plume.r[np.flatnonzero(col >= threshold)[-1]]
         if np.any(col >= threshold) else 0.0)
        for col in plume.jet_fraction.T])
    theta = np.linspace(0.0, 2.0 * np.pi, n_theta, endpoint=False)
    x = np.repeat((plume.x + axial_offset)[:, None], n_theta, axis=1).ravel()
    return np.stack([x, (r_edge[:, None] * np.cos(theta)).ravel(),
                     (r_edge[:, None] * np.sin(theta)).ravel()], axis=-1)


def framed_camera(field: Field, width: int, height: int, *, yaw_deg: float = 0.0,
                  pitch_deg: float = 8.0, fov_deg: float = 26.0,
                  margin: float = 0.10, centre: float = 0.5,
                  extra_points: Optional[np.ndarray] = None) -> Camera:
    """A camera at the given angle, pulled back until the whole engine fits.

    Framing by hand means re-guessing a distance every time the angle or the
    aspect ratio changes, and an animation that orbits the engine would breathe
    in and out as the silhouette turns.  Here the distance is solved for: the
    outer surface of revolution is projected into the image and the eye is
    pushed back until every point sits inside `1 - margin` of the frame.

    `yaw_deg` is 0 for a side elevation and swings towards the nozzle axis;
    `pitch_deg` raises the eye above the axis.
    """
    hull = _hull(field)
    if extra_points is not None and len(extra_points):
        hull = np.concatenate([hull, np.asarray(extra_points, dtype=float)])
    # Centre on what is actually being framed, which with a plume attached is
    # well downstream of the middle of the engine.
    span = float(hull[:, 0].max() - hull[:, 0].min())
    target = np.array([hull[:, 0].min() + span * centre, 0.0, 0.0])
    L = max(field.length, span)
    yaw, pitch = np.radians(yaw_deg), np.radians(pitch_deg)
    direction = np.array([np.sin(yaw) * np.cos(pitch), np.sin(pitch),
                          np.cos(yaw) * np.cos(pitch)])
    limit = 1.0 - margin

    def fits(distance: float) -> bool:
        cam = Camera(eye=target + distance * direction, target=target,
                     up=np.array([0.0, 1.0, 0.0]), fov_deg=fov_deg,
                     width=width, height=height)
        return bool(np.abs(cam.project(hull)).max() <= limit)

    lo, hi = 0.25 * L, 0.5 * L
    for _ in range(64):                       # bracket: grow until it fits
        if fits(hi):
            break
        lo, hi = hi, hi * 1.5
    else:
        raise RuntimeError("could not frame the field")
    for _ in range(40):                       # then bisect onto the tightest fit
        mid = 0.5 * (lo + hi)
        if fits(mid):
            hi = mid
        else:
            lo = mid
    return Camera(eye=target + hi * direction, target=target,
                  up=np.array([0.0, 1.0, 0.0]), fov_deg=fov_deg,
                  width=width, height=height)


@dataclass
class Shell:
    """How the nozzle wall is drawn.

    This is a lighting choice, not physics: the solver says nothing about what
    the outside of the engine looks like.  `transmittance` is the fraction of a
    ray that survives one crossing of the wall -- 1.0 makes the shell invisible,
    0.0 makes it solid.  A value near zero gives the honest opaque engine; a
    value near 0.3 gives a cutaway that shows the whole gas column at the cost
    of pretending the wall is glass.
    """

    transmittance: float = 0.30
    thickness: float = 0.0       # m, radial thickness of the drawn wall
    ambient: float = 0.004       # light reaching the wall from nowhere in particular
    diffuse: float = 0.055       # lambertian response to the key light
    rim: float = 0.42            # grazing-angle (fresnel) response
    rim_power: float = 3.0       # how tightly the rim hugs the silhouette
    section: float = 0.0         # brightness of the sectioned (cut) face
    glow: float = 0.55           # how much of the gas emission the inner face picks up
    cut: bool = False            # remove the shell nearest the camera
    quarter: bool = False        # remove only the near QUARTER, above the axis
    head: float = 0.0            # m, depth of the injector dome closing x = 0
    channels: int = 0            # cooling channels to cut into the jacket
    land_fraction: float = 0.5   # share of the channel pitch that is solid land
    interior_dim: float = 1.0    # how much the far inner face is dimmed
    manifolds: Tuple[Tuple[float, float], ...] = ()   # (x, tube radius), metres
    key: Tuple[float, float, float] = (0.34, 0.42, 0.55)
    light: Tuple[float, float, float] = (-0.35, 0.80, 0.49)
    specular: float = 0.0        # Blinn-Phong highlight, which is what makes a
    shininess: float = 40.0      # curved metal surface read as curved metal


@dataclass
class Geometry:
    """Everything about the scene that does not change when the flow does."""

    paths: np.ndarray        # (Npix, Nbins) float32, metres of gas per slice
    wall_rgb: np.ndarray     # (Npix, 3) float32, shaded wall radiance
    wall_depth: np.ndarray   # (Npix,) float32, distance to the wall hit
    shape: Tuple[int, int]


def build_geometry(field: Field, camera: Camera, *, steps: int = 520,
                   shell: Optional["Shell"] = None,
                   seed: int = 20260917) -> Geometry:
    """Ray-march the fixed geometry once.

    The nozzle shell is drawn TRANSLUCENT rather than opaque.  A real engine
    wall is not, of course; the point of this image is the interior field, and a
    solid shell would hide exactly the thing the solver computed.  Each crossing
    of the wall attenuates the ray and adds a metallic surface term, so the
    structure still reads as structure while the flow shows through it.

    Everything here depends on geometry alone, never on the flow, which is what
    lets a whole animation reuse one call: the returned `paths` matrix is the
    transmittance-weighted path length each pixel spends in each axial slice.
    """
    shell = shell if shell is not None else Shell()
    origin, dirs = camera.rays()
    n_pix = dirs.shape[0]
    n_bins = field.x.size

    L = field.length
    r_max = field.max_radius * 1.08

    # March between the two planes that bound the nozzle along the view ray,
    # padded so grazing rays are not clipped.
    span = np.linalg.norm(camera.eye - camera.target) + L + 2.0 * r_max
    t0 = max(0.0, float(np.linalg.norm(camera.eye - camera.target) - L - 2.0 * r_max))
    ds = np.float32((span - t0) / steps)
    if shell.thickness > 0.0 and ds > 0.5 * shell.thickness:
        raise ValueError(
            f"march step {ds * 1e3:.2f} mm is too coarse for a "
            f"{shell.thickness * 1e3:.2f} mm wall; rays would tunnel through it. "
            f"Use at least {int(np.ceil((span - t0) / (0.5 * shell.thickness)))} steps.")

    bin_width = L / n_bins

    paths = np.zeros((n_pix, n_bins), dtype=np.float32)
    wall_rgb = np.zeros((n_pix, 3), dtype=np.float32)
    wall_depth = np.full(n_pix, np.inf, dtype=np.float32)
    trans = np.ones(n_pix, dtype=np.float32)

    radius32 = field.radius.astype(np.float32)
    wall = _Wall(shell, camera, field, dirs)

    # Stratified jitter: a fixed step lands the sample points on the same
    # surfaces for neighbouring rays and the slice boundaries print as rings.
    # One random offset per ray turns that structure into noise, which the
    # eye ignores and the bloom filter smooths away.
    rng = np.random.default_rng(seed)
    jitter = rng.random(n_pix, dtype=np.float32)

    for i in range(steps):
        t = (np.float32(t0) + (np.float32(i) + jitter) * ds)[:, None]
        p = origin[None, :] + dirs * t
        px, py, pz = p[:, 0], p[:, 1], p[:, 2]
        r = np.sqrt(py * py + pz * pz)

        inside_x = (px >= 0.0) & (px < L) & (trans > 0.004)
        idx = np.clip((px / bin_width).astype(np.int32), 0, n_bins - 1)
        in_gas = inside_x & (r < radius32[idx])
        # Emission is attenuated by whatever the ray has already passed through.
        if in_gas.any():
            paths[in_gas, idx[in_gas]] += ds * trans[in_gas]

        wall.step(inside_x, idx, px, r, py, pz, t, trans, wall_rgb, wall_depth)

    return Geometry(paths=paths, wall_rgb=wall_rgb, wall_depth=wall_depth,
                    shape=(camera.height, camera.width))


def emission_per_slice(field: Field, t_grid: Optional[np.ndarray] = None) -> np.ndarray:
    """Radiance emitted per metre of path in each axial slice, linear sRGB."""
    t = np.clip(field.temperature, 1.0, None)
    chroma = blackbody_rgb(t)                       # unit-luminance hue
    power = field.density * SIGMA_SB * t ** 4       # optically thin gray emitter
    return (chroma * power[:, None]).astype(np.float32)


def render(geometry: Geometry, field: Field, *, exposure: float = 1.0,
           bloom: float = 0.7, background: float = 0.004,
           latitude: float = 900.0, reference: Optional[float] = None) -> np.ndarray:
    """Composite one frame: gas emission, the shell, bloom and a tone curve.

    The scene's linear dynamic range is enormous -- emission goes as rho T^4,
    and between the chamber and the exit plane rho falls 85x while T^4 falls 8x,
    so the chamber is some seven hundred times brighter than the exhaust.  A
    linear mapping would show the chamber and nothing else.  `latitude` sets a
    logarithmic response, which is what a camera does with a scene like this;
    the *ordering* and the hue are the solver's, the compression is display.

    `reference` fixes the exposure across an animation so that a frame getting
    brighter means the engine got brighter, not that the tone curve moved.
    """
    radiance = geometry.paths @ emission_per_slice(field)   # (Npix, 3)
    return compose(radiance, geometry.wall_rgb, geometry.shape, exposure=exposure,
                   bloom=bloom, background=background, latitude=latitude,
                   reference=reference)


def compose(radiance: np.ndarray, wall_rgb: np.ndarray, shape: Tuple[int, int], *,
            exposure: float = 1.0, bloom: float = 0.7, background: float = 0.004,
            latitude: float = 900.0, reference: Optional[float] = None,
            bloom_threshold: float = 0.55) -> np.ndarray:
    """Gas radiance plus the shell, through the camera response, to display sRGB."""
    h, w = shape
    ref = reference if reference is not None else float(
        np.percentile(radiance.sum(axis=1), 99.95))
    img = np.clip(radiance / ref if ref > 0 else radiance, 0.0, None) * exposure
    img = np.log1p(latitude * img) / np.log1p(latitude)
    img = (img + wall_rgb).reshape(h, w, 3) + background

    if bloom > 0.0:
        bright = np.clip(img - bloom_threshold, 0.0, None)
        glow = sum(gaussian_filter(bright, sigma=(sg, sg, 0)) * wgt
                   for sg, wgt in ((3.0, 0.50), (9.0, 0.34), (24.0, 0.22)))
        img = img + bloom * glow

    return _tonemap(img)


def exposure_reference(radiance: np.ndarray, percentile: float = 99.95) -> float:
    """The level `compose` maps to 1; pass it back to hold exposure across frames."""
    return float(np.percentile(radiance.sum(axis=1), percentile))


_LUMA = np.array([0.2126, 0.7152, 0.0722])


def _tonemap(img: np.ndarray) -> np.ndarray:
    """Filmic shoulder applied to luminance, then the sRGB transfer.

    Running the shoulder on each channel separately is what turned the chamber
    flat white: red saturates first, then green, and by the time the curve has
    compressed the highlight the hue has been thrown away.  A 3515 K chamber is
    not white, it is deep orange, and that is the one colour in the frame the
    solver actually determined -- so the curve is applied to luminance and the
    chromaticity carried through it unchanged.

    Channels can then land above 1.  Those are bleached smoothly towards
    neutral rather than clipped, which is the way film runs out of headroom:
    the very centre of the chamber goes white because it is brighter than the
    medium can hold, not because the arithmetic lost the hue on the way.
    """
    img = np.clip(img, 0.0, None)
    lum = img @ _LUMA
    tone = (lum * (2.51 * lum + 0.03)) / (lum * (2.43 * lum + 0.59) + 0.14)
    img = img * (tone / np.maximum(lum, 1e-9))[..., None]

    peak = img.max(axis=-1, keepdims=True)
    bleach = np.clip(peak - 1.0, 0.0, 1.0)
    img = img * (1.0 - bleach) + peak * bleach
    img = img / np.maximum(peak, 1.0)
    return np.clip(np.clip(img, 0.0, 1.0) ** (1.0 / 2.2), 0.0, 1.0)


def to_uint8(img: np.ndarray) -> np.ndarray:
    return (np.clip(img, 0.0, 1.0) * 255.0 + 0.5).astype(np.uint8)


# ---------------------------------------------------------------------------
# Scene marching: nozzle interior plus the computed plume
# ---------------------------------------------------------------------------


@dataclass
class Scene:
    """What the camera looks at: the engine, its shell, and optionally a plume.

    The nozzle field is quasi-1D, so its emission depends on x alone; the plume
    is a genuine (x, r) field and is interpolated.  Both are revolved about the
    same axis and share the exit plane, which is where the plume's own x
    coordinate starts.
    """

    field: Field
    shell: Shell
    plume: object = None          # ignis_viz.plume.PlumeField, or None
    plume_reach: float = 0.0      # m past the exit plane to march; 0 = all of it


def _plume_emission_table(plume) -> np.ndarray:
    """Linear-sRGB radiance per metre of path, on the plume's own (r, x) grid.

    Emission is masked by the passive jet-fraction scalar so that the cold
    ambient air the jet is expanding into never contributes light.
    """
    t = np.clip(plume.temperature, 1.0, None)
    chroma = blackbody_rgb(t.ravel()).reshape(t.shape + (3,))
    power = plume.jet_fraction * plume.density * SIGMA_SB * t ** 4
    return (chroma * power[..., None]).astype(np.float32)


def _cylinder_range(origin: np.ndarray, dirs: np.ndarray, radius: float,
                    x_lo: float, x_hi: float):
    """Entry and exit parameters of each ray through the scene's bounding cylinder.

    Returns (t_in, t_out, hit).  The cylinder is coaxial with the nozzle, so
    this is a quadratic in the two transverse components clipped by the axial
    slab -- cheap, and it both drops the rays that see nothing and tightens the
    march to the part of each ray that can contribute.
    """
    ox, oy, oz = (float(origin[0]), float(origin[1]), float(origin[2]))
    dx, dy, dz = dirs[:, 0], dirs[:, 1], dirs[:, 2]

    a = dy * dy + dz * dz
    b = 2.0 * (oy * dy + oz * dz)
    c = oy * oy + oz * oz - radius * radius
    disc = b * b - 4.0 * a * c
    ok = (disc > 0.0) & (a > 1e-20)
    root = np.sqrt(np.maximum(disc, 0.0))
    denom = np.where(a > 1e-20, 2.0 * a, 1.0)
    t_near = (-b - root) / denom
    t_far = (-b + root) / denom

    # Axial slab.  A ray parallel to the axis is inside it or misses entirely.
    parallel = np.abs(dx) < 1e-12
    safe_dx = np.where(parallel, 1.0, dx)
    ta = (x_lo - ox) / safe_dx
    tb = (x_hi - ox) / safe_dx
    slab_lo = np.where(parallel, -np.inf, np.minimum(ta, tb))
    slab_hi = np.where(parallel, np.inf, np.maximum(ta, tb))
    ok &= ~(parallel & ((ox < x_lo) | (ox > x_hi)))

    t_in = np.maximum(np.maximum(t_near, slab_lo), 0.0)
    t_out = np.minimum(t_far, slab_hi)
    ok &= t_out > t_in
    return t_in.astype(np.float32), t_out.astype(np.float32), ok


class _Wall:
    """The engine hardware, as an implicit solid, marched along with the gas.

    A bare shell of revolution does not read as a rocket engine, because a
    rocket engine is not one.  What is drawn here is an implicit solid built
    from pieces, `inside` returning how deep a point lies within it:

      * the cooling jacket -- the hot-gas wall plus the coolant passage, at the
        thickness Ignis's cooling solver sizes, with the channels themselves cut
        into it at the computed count and land fraction;
      * an injector dome closing the chamber at x = 0;
      * coolant manifolds, rings around the engine where the jacket is fed and
        collected.

    Of those, the wall thickness, the channel count and the land fraction are
    the solver's, and the manifold positions follow from its counterflow
    circuit -- coolant in at the nozzle exit, out at the injector.  The dome
    shape, the manifold sizes and the surface finish are drawing choices: the
    model has nothing to say about them, and they carry no result.

    Surface normals come from the gradient of that same implicit function, so
    every piece is shaded consistently and adding another piece needs no new
    normal maths.  The gradient is evaluated only where a ray first enters the
    solid, which is a few thousand points rather than the whole march.

    With `cut` set the solid nearest the camera is removed; with `quarter` only
    the part of it above the axis is, so one half of the picture shows the
    outside of the engine and the other the flow through it.
    """

    def __init__(self, shell: Shell, camera: Camera, field: Field, dirs: np.ndarray):
        self.shell = shell
        self.dirs = dirs
        self.fx = field.x.astype(np.float64)
        self.fr = field.radius.astype(np.float64)
        self.length = float(field.length)
        self.key = np.array(shell.key, dtype=np.float32)
        light = np.array(shell.light, dtype=np.float32)
        self.light = light / np.linalg.norm(light)
        power = (field.density * field.temperature ** 4).astype(np.float32)
        self.slice_power = power / max(float(power.max()), 1e-30)
        self.slice_chroma = blackbody_rgb(
            np.clip(field.temperature, 1.0, None)).astype(np.float32)

        n = camera.eye - camera.target
        n = np.array([0.0, n[1], n[2]])
        norm = float(np.linalg.norm(n))
        self.cut_normal = (n / norm).astype(np.float32) if norm > 1e-12 else None
        # The second cut plane of a quarter section: the camera's own up
        # direction, flattened perpendicular to the axis, so "above" in the
        # picture is "above" in the cut.
        _, _, up = camera.basis()
        u = np.array([0.0, up[1], up[2]])
        u_norm = float(np.linalg.norm(u))
        self.up_normal = (u / u_norm).astype(np.float32) if u_norm > 1e-12 else None

        self.was_solid = np.zeros(dirs.shape[0], dtype=bool)

    # -- the solid --------------------------------------------------------
    def inside(self, px, py, pz):
        """Depth inside the sectioned solid; positive inside, in metres."""
        hardware = self._hardware(px, py, pz)
        keep = self._keep(py, pz)
        return hardware if keep is None else np.minimum(hardware, keep)

    def cut_active(self, px, py, pz):
        """True where the section, not the hardware, is the binding surface.

        A sectioned face has to be found by asking which constraint the surface
        came from.  Testing whether the normal happens to point along the cut
        plane catches every other face that shares that direction too -- on a
        surface of revolution that is a whole wedge of the bell, and it prints
        as one.
        """
        keep = self._keep(py, pz)
        if keep is None:
            return np.zeros(np.shape(px), dtype=bool)
        return keep < self._hardware(px, py, pz)

    def _keep(self, py, pz):
        """Depth into the part of space the section keeps, or None if uncut."""
        if not self.shell.cut or self.cut_normal is None:
            return None
        keep = -(py * self.cut_normal[1] + pz * self.cut_normal[2])
        if self.shell.quarter and self.up_normal is not None:
            keep = np.maximum(keep, -(py * self.up_normal[1] + pz * self.up_normal[2]))
        return keep

    def _hardware(self, px, py, pz):
        """Depth inside the hardware alone, before the section is applied."""
        shell = self.shell
        r = np.sqrt(py * py + pz * pz)
        r_in = np.interp(px, self.fx, self.fr)
        r_wall = r_in + shell.thickness

        if shell.channels > 0:
            # Channels are milled into the jacket, so the jacket is at full
            # thickness over a land and cut back over a passage.  The groove is
            # drawn as a smooth corrugation rather than a square slot: at three
            # hundred channels the pitch falls below a pixel near the
            # silhouette, and a hard edge there samples into speckle while a
            # smooth one only softens.
            phase = np.arctan2(pz, py) * shell.channels
            ripple = 0.5 * (1.0 + np.cos(phase))
            groove = shell.thickness * 0.45 * (1.0 - shell.land_fraction) * 2.0
            r_out = r_wall - groove * (1.0 - ripple)
        else:
            r_out = r_wall

        depth = np.minimum(r - r_in, r_out - r)
        depth = np.where((px >= 0.0) & (px <= self.length), depth, -1.0)

        if shell.head > 0.0:
            # A torispherical cap closing the chamber at the injector end.
            a, b = shell.head, float(self.fr[0] + shell.thickness)
            q = np.sqrt(np.square(np.minimum(px, 0.0) / a) + np.square(r / b))
            depth = np.maximum(depth, np.where(px <= 0.0, b * (1.0 - q), -1.0))

        for x_m, tube in shell.manifolds:
            ring = np.interp(x_m, self.fx, self.fr) + shell.thickness
            depth = np.maximum(
                depth, tube - np.sqrt(np.square(px - x_m) + np.square(r - ring)))

        return depth

    def _normal(self, px, py, pz, h=2.0e-4):
        """Outward surface normal, from the gradient of `inside`."""
        gx = self.inside(px + h, py, pz) - self.inside(px - h, py, pz)
        gy = self.inside(px, py + h, pz) - self.inside(px, py - h, pz)
        gz = self.inside(px, py, pz + h) - self.inside(px, py, pz - h)
        n = np.stack([gx, gy, gz], axis=-1).astype(np.float32)
        # `inside` grows towards the interior, so its gradient points inwards.
        n = -n / np.maximum(np.linalg.norm(n, axis=-1, keepdims=True), 1e-12)
        return n

    # -- marching ---------------------------------------------------------
    def step(self, inside_x, idx, px, r, py, pz, t_col, trans, wall_rgb, wall_depth):
        """Shade whatever surface each ray entered on this step."""
        shell = self.shell
        solid = (self.inside(px, py, pz) > 0.0) & (trans > 0.004)
        entering = solid & ~self.was_solid
        self.was_solid = solid
        if not entering.any():
            return

        hit = entering
        normals = self._normal(px[hit], py[hit], pz[hit])
        nx, ny, nz = normals[:, 0], normals[:, 1], normals[:, 2]
        on_cut = self.cut_active(px[hit], py[hit], pz[hit])
        view = -self.dirs[hit]
        # A face pointing away from the camera is one we are seeing from inside.
        from_outside = (nx * view[:, 0] + ny * view[:, 1] + nz * view[:, 2]) > 0.0
        flip = np.where(from_outside, 1.0, -1.0).astype(np.float32)
        nx, ny, nz = nx * flip, ny * flip, nz * flip

        facing = np.abs(nx * view[:, 0] + ny * view[:, 1] + nz * view[:, 2])
        lam = np.clip(nx * self.light[0] + ny * self.light[1] + nz * self.light[2],
                      0.0, 1.0)
        fres = np.clip(1.0 - facing, 0.0, 1.0) ** shell.rim_power
        spec = 0.0
        if shell.specular > 0.0:
            half = self.light[None, :] + view
            half /= np.maximum(np.linalg.norm(half, axis=1, keepdims=True), 1e-9)
            spec = (shell.specular
                    * np.clip(nx * half[:, 0] + ny * half[:, 1] + nz * half[:, 2],
                              0.0, 1.0) ** shell.shininess)[:, None]

        b = idx[hit]
        # Only a face looking back into the gas is lit by it.
        interior = (~from_outside & ~on_cut).astype(np.float32)[:, None]
        glow = self.slice_power[b][:, None] * self.slice_chroma[b] * interior * shell.glow
        # The far side of the bell, seen from inside, sits behind the gas -- at
        # full brightness it washes out the flow it is meant to contain.
        dim = 1.0 - interior * (1.0 - shell.interior_dim)
        body = ((shell.ambient + shell.diffuse * lam[:, None]) * self.key[None, :]
                + shell.rim * fres[:, None] * self.key[None, :]) * dim
        shade = (body + spec * (1.0 - interior)
                 + shell.section * on_cut.astype(np.float32)[:, None] * self.key[None, :]
                 + glow)
        wall_rgb[hit] += (shade * trans[hit][:, None]).astype(np.float32)
        wall_depth[hit] = np.minimum(wall_depth[hit], t_col[hit, 0])
        trans[hit] *= np.float32(shell.transmittance)


class _PlumeLookup:
    """Bilinear sampling of the plume's emission table on its own (x, r) grid."""

    def __init__(self, plume, reach: float):
        self.table = _plume_emission_table(plume)          # (nr, nx, 3)
        self.dx = float(plume.x[1] - plume.x[0])
        self.dr = float(plume.r[1] - plume.r[0])
        self.x0, self.r0 = float(plume.x[0]), float(plume.r[0])
        self.n_r, self.n_x = plume.temperature.shape
        full = self.n_x * self.dx
        # Marching plume that is off the side of the frame costs steps and
        # changes nothing, and with a long domain that is most of it.
        self.length = min(full, reach) if reach > 0.0 else full
        self.r_max = self.n_r * self.dr
        # How wide the jet actually is, as against how wide the grid is.  The
        # emission is masked by the jet-fraction scalar, so the ambient air out
        # to the grid edge contributes nothing -- but a bounding volume drawn
        # around the grid makes every ray march through all of it.
        columns = int(np.ceil(self.length / self.dx))
        jet = plume.jet_fraction[:, :max(columns, 1)]
        rows = np.flatnonzero((jet >= 0.5).any(axis=1))
        self.r_visible = (float(plume.r[rows[-1]]) + 2.0 * self.dr
                          if rows.size else 0.0)

    @classmethod
    def build(cls, scene: "Scene"):
        return None if scene.plume is None else cls(scene.plume, scene.plume_reach)

    def sample(self, x: np.ndarray, r: np.ndarray) -> np.ndarray:
        fx = (x - self.x0) / self.dx
        fr = (r - self.r0) / self.dr
        i0 = np.clip(fx.astype(np.int32), 0, self.n_x - 2)
        j0 = np.clip(fr.astype(np.int32), 0, self.n_r - 2)
        tx = np.clip(fx - i0, 0.0, 1.0)[:, None]
        tr = np.clip(fr - j0, 0.0, 1.0)[:, None]
        t = self.table
        return ((t[j0, i0] * (1.0 - tx) + t[j0, i0 + 1] * tx) * (1.0 - tr)
                + (t[j0 + 1, i0] * (1.0 - tx) + t[j0 + 1, i0 + 1] * tx) * tr)


def march_scene(scene: Scene, camera: Camera, *, steps: int = 900,
                max_step: Optional[float] = None,
                chunk: int = 48000, seed: int = 20260917) -> Tuple[np.ndarray, np.ndarray]:
    """Integrate emission along every camera ray.

    Returns (gas radiance, shell radiance, wall depth): the first two are
    (Npix, 3), the third the distance to the nearest opaque surface on each ray
    (infinite where there is none), which is what lets screen-space overlays
    such as the flow tracers be hidden behind the engine.

    Unlike `build_geometry` this accumulates colour directly instead of a
    per-slice path matrix: with a two-dimensional plume the matrix would need a
    column per (x, r) cell, which is a hundred thousand columns and no longer a
    useful thing to precompute.  The trade is that a frame costs a full march,
    so this path is for stills and the path-matrix one is for animation.
    """
    field, shell = scene.field, scene.shell
    origin, dirs_all = camera.rays()
    n_all = dirs_all.shape[0]
    n_bins = field.x.size

    L = field.length
    bin_width = L / n_bins
    radius32 = field.radius.astype(np.float32)
    nozzle_emission = emission_per_slice(field)

    plume = _PlumeLookup.build(scene)
    p_len = plume.length if plume is not None else 0.0
    p_visible = plume.r_visible if plume is not None else 0.0

    # Every ray is clipped to the scene's bounding cylinder before marching, and
    # rays that miss it are dropped.  Stepping uniformly from the eye across a
    # bounding sphere spends most of the work on empty space -- typically two
    # thirds of the steps and two thirds of the pixels -- and the step has to
    # stay smaller than the wall is thick or rays tunnel straight through it,
    # so that waste is what decides whether an animated camera is affordable.
    bulge = max((tube for _, tube in shell.manifolds), default=0.0)
    r_scene = max(field.max_radius + shell.thickness + 2.0 * bulge,
                  p_visible) * 1.001
    t_in, t_out, visible = _cylinder_range(origin, dirs_all, r_scene,
                                           -shell.head, L + p_len)
    active = np.flatnonzero(visible)
    dirs = dirs_all[active]
    n_pix = active.size
    if n_pix == 0:
        return (np.zeros((n_all, 3), dtype=np.float32),
                np.zeros((n_all, 3), dtype=np.float32),
                np.full(n_all, np.inf, dtype=np.float32))
    t_in, t_out = t_in[active], t_out[active]
    if max_step is not None:
        # Hold the sample spacing fixed instead of the step count, so a moving
        # camera keeps the same fidelity as the view changes shape.
        steps = max(8, int(np.ceil(float((t_out - t_in).max()) / max_step)))
    ds = ((t_out - t_in) / steps).astype(np.float32)
    if shell.thickness > 0.0:
        coarsest = float(ds.max())
        if coarsest > 0.5 * shell.thickness:
            raise ValueError(
                f"march step {coarsest * 1e3:.2f} mm is too coarse for a "
                f"{shell.thickness * 1e3:.2f} mm wall; rays would tunnel through "
                f"it.  Use at least "
                f"{int(np.ceil(steps * coarsest / (0.5 * shell.thickness)))} steps.")

    gas_all = np.zeros((n_pix, 3), dtype=np.float32)
    wall_all = np.zeros((n_pix, 3), dtype=np.float32)
    depth_all = np.full(n_pix, np.inf, dtype=np.float32)
    rng = np.random.default_rng(seed)
    jitter_all = rng.random(n_pix, dtype=np.float32)

    # One pass over all rays touches tens of megabytes per step and spends the
    # whole march missing cache.  Marching a slice of rays all the way through
    # before moving to the next keeps the working set resident, which is worth
    # more here than any amount of arithmetic tuning.
    for lo in range(0, n_pix, chunk):
        hi = min(lo + chunk, n_pix)
        dirs_c = dirs[lo:hi]
        t_in_c, ds_c = t_in[lo:hi], ds[lo:hi]
        jitter = jitter_all[lo:hi]
        gas = gas_all[lo:hi]
        wall_rgb = wall_all[lo:hi]
        trans = np.ones(hi - lo, dtype=np.float32)
        wall = _Wall(shell, camera, field, dirs_c)
        wall_depth = depth_all[lo:hi]
        _march_chunk(steps, origin, dirs_c, t_in_c, ds_c, jitter, trans, gas,
                     wall_rgb, wall, wall_depth, L, n_bins, bin_width, radius32,
                     nozzle_emission, plume)

    out_gas = np.zeros((n_all, 3), dtype=np.float32)
    out_wall = np.zeros((n_all, 3), dtype=np.float32)
    out_depth = np.full(n_all, np.inf, dtype=np.float32)
    out_gas[active] = gas_all
    out_wall[active] = wall_all
    out_depth[active] = depth_all
    return out_gas, out_wall, out_depth


def _march_chunk(steps, origin, dirs, t_in, ds, jitter, trans, gas, wall_rgb,
                 wall, wall_depth, L, n_bins, bin_width, radius32,
                 nozzle_emission, plume):
    """March one slice of rays all the way through the scene."""
    for i in range(steps):
        t = (t_in + (np.float32(i) + jitter) * ds)[:, None]
        p = origin[None, :] + dirs * t
        px, py, pz = p[:, 0], p[:, 1], p[:, 2]
        r = np.sqrt(py * py + pz * pz)
        live = trans > 0.004

        # --- nozzle interior -------------------------------------------------
        inside_x = (px >= 0.0) & (px < L) & live
        idx = np.clip((px / bin_width).astype(np.int32), 0, n_bins - 1)
        in_gas = inside_x & (r < radius32[idx])
        if in_gas.any():
            gas[in_gas] += (nozzle_emission[idx[in_gas]]
                            * (ds[in_gas] * trans[in_gas])[:, None])
        wall.step(inside_x, idx, px, r, py, pz, t, trans, wall_rgb, wall_depth)

        # --- plume -----------------------------------------------------------
        if plume is not None:
            xp = px - L
            # From the exit plane, not from the plume's first cell centre: the
            # half-cell between them is a grid detail, and leaving it empty
            # prints a dark seam right where the flow leaves the engine.
            inside_p = (xp >= 0.0) & (xp < plume.length) & (r < plume.r_max) & live
            if inside_p.any():
                gas[inside_p] += (plume.sample(xp[inside_p], r[inside_p])
                                  * (ds[inside_p] * trans[inside_p])[:, None])


def shift_camera(camera: Camera, ndc_dx: float = 0.0, ndc_dy: float = 0.0) -> Camera:
    """Slide the subject within the frame without changing the viewing angle.

    Auto-framing centres the engine, but a cover wants it off to one side with
    the plume running out of the frame.  Translating the eye and the target
    together keeps the view direction, the distance and therefore the whole
    perspective; only where the subject lands in the image changes.
    """
    forward, right, up = camera.basis()
    depth = float(np.linalg.norm(camera.target - camera.eye))
    half_h = np.tan(np.radians(camera.fov_deg) * 0.5)
    half_w = half_h * camera.width / camera.height
    delta = (-ndc_dx * half_w * depth) * right + (-ndc_dy * half_h * depth) * up
    return Camera(eye=camera.eye + delta, target=camera.target + delta, up=camera.up,
                  fov_deg=camera.fov_deg, width=camera.width, height=camera.height)


def pixel_coords(camera: Camera, points: np.ndarray) -> np.ndarray:
    """Image pixel coordinates of world points, (N, 2), origin top-left."""
    ndc = camera.project(points)
    return np.stack([(ndc[:, 0] + 1.0) * 0.5 * camera.width,
                     (1.0 - ndc[:, 1]) * 0.5 * camera.height], axis=-1)
