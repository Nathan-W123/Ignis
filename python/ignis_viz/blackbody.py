"""Blackbody colour, computed rather than chosen.

The Ignis-M1 chamber sits at 3524 K and its exhaust leaves at 2099 K.  Gas at
those temperatures is genuinely incandescent, so the colour of a rendered
Ignis flow field is not a design decision: it is the chromaticity of a Planck
radiator at the temperature the solver computed.

    B(lambda, T) = 2 h c^2 / lambda^5 / (exp(h c / (lambda k T)) - 1)

is integrated against the CIE 1931 2-degree colour-matching functions, and the
resulting tristimulus value is converted to linear sRGB with the standard
primaries.  The CMFs use the analytic multi-lobe Gaussian fits of

    Wyman, Sloan and Shirley, "Simple Analytic Approximations to the CIE XYZ
    Color Matching Functions", Journal of Computer Graphics Techniques 2(2),
    1-11 (2013),

which reproduce the tabulated functions to well under a percent -- far below
anything a screen can show.

The returned colour is normalised to unit luminance, so hue carries
temperature and brightness is left to carry whatever the caller is encoding
with it (in the renderer: emitted power).
"""
from __future__ import annotations

import numpy as np

# Physical constants (CODATA), SI.
_H = 6.62607015e-34      # J s
_C = 2.99792458e8        # m/s
_KB = 1.380649e-23       # J/K

# Linear sRGB primaries, IEC 61966-2-1.
_XYZ_TO_SRGB = np.array([
    [3.2406255, -1.5372080, -0.4986286],
    [-0.9689307, 1.8757561, 0.0415175],
    [0.0557101, -0.2040211, 1.0569959],
])


def _gauss(x: np.ndarray, mu: float, s1: float, s2: float) -> np.ndarray:
    """Piecewise-Gaussian lobe: sigma differs either side of the peak."""
    s = np.where(x < mu, s1, s2)
    return np.exp(-0.5 * ((x - mu) / s) ** 2)


def cie_cmf(lam_nm: np.ndarray) -> tuple:
    """CIE 1931 2-degree colour-matching functions (Wyman et al. 2013)."""
    x = (1.056 * _gauss(lam_nm, 599.8, 37.9, 31.0)
         + 0.362 * _gauss(lam_nm, 442.0, 16.0, 26.7)
         - 0.065 * _gauss(lam_nm, 501.1, 20.4, 26.2))
    y = (0.821 * _gauss(lam_nm, 568.8, 46.9, 40.5)
         + 0.286 * _gauss(lam_nm, 530.9, 16.3, 31.1))
    z = (1.217 * _gauss(lam_nm, 437.0, 11.8, 36.0)
         + 0.681 * _gauss(lam_nm, 459.0, 26.0, 13.8))
    return x, y, z


def planck(lam_m: np.ndarray, T: np.ndarray) -> np.ndarray:
    """Spectral radiance of a blackbody, W / (m^2 sr m)."""
    lam = lam_m[None, :]
    t = np.asarray(T, dtype=np.float64)[:, None]
    # exp argument is large at short wavelength / low temperature; clip so the
    # exponential cannot overflow before the reciprocal makes it irrelevant.
    arg = np.clip(_H * _C / (lam * _KB * t), 0.0, 700.0)
    return 2.0 * _H * _C ** 2 / lam ** 5 / np.expm1(arg)


def blackbody_rgb(T, samples: int = 256) -> np.ndarray:
    """Linear sRGB chromaticity of a blackbody at temperature `T` (K).

    Returns an (N, 3) array normalised to unit luminance, gamut-clipped.
    """
    T = np.atleast_1d(np.asarray(T, dtype=np.float64))
    # Below roughly the Draper point a body emits no visible light at all, so
    # the integrals over the visible band underflow and normalising them to
    # unit luminance divides by nothing -- which comes out as arbitrary colour
    # rather than as no colour.  The chromaticity is clamped at that floor;
    # what makes a cold body dark is its radiated power, which is applied
    # separately and does go to zero.
    T = np.maximum(T, 700.0)
    lam_nm = np.linspace(380.0, 780.0, samples)
    dlam = lam_nm[1] - lam_nm[0]
    xb, yb, zb = cie_cmf(lam_nm)
    rad = planck(lam_nm * 1e-9, T)                     # (N, samples)

    xyz = np.stack([(rad * xb).sum(axis=1),
                    (rad * yb).sum(axis=1),
                    (rad * zb).sum(axis=1)], axis=1) * dlam
    rgb = xyz @ _XYZ_TO_SRGB.T
    # Out-of-gamut blues at high temperature go slightly negative; desaturate
    # towards white rather than clipping to black, which would shift the hue.
    deficit = np.maximum(0.0, -rgb.min(axis=1, keepdims=True))
    rgb = rgb + deficit
    peak = rgb.max(axis=1, keepdims=True)
    rgb = np.divide(rgb, peak, out=np.zeros_like(rgb), where=peak > 0)
    return np.clip(rgb, 0.0, 1.0)


def blackbody_lookup(t_min: float, t_max: float, n: int = 512) -> tuple:
    """A (temperature grid, RGB table) pair for fast lookup during rendering."""
    grid = np.linspace(t_min, t_max, n)
    return grid, blackbody_rgb(grid).astype(np.float32)
