"""A 3-D viewport for the engine and its flow, rasterised in NumPy.

WHY NOT OPENGL
--------------
The obvious way to draw a 3-D scene in Qt is QOpenGLWidget.  It is not used
here for two reasons.  A GPU path cannot be exercised in the headless
environment this is developed and tested in, so it would ship unverified; and
it would put a driver between the project and its output when the scene is a
few thousand triangles -- small enough that a z-buffered rasteriser in NumPy
keeps up with an orbiting camera.  Nothing here needs a GPU, so nothing here
asks for one.

WHAT IS DRAWN
-------------
Two things, which is how a CFD viewport usually reads:

  * The engine, as a surface of revolution built from the contour Ignis
    solved, swept through less than a full turn so the wedge that is missing
    lets you see the gas-side wall on the far side.  The inner surface is
    coloured by the flow field; the outer surface is left as metal, because
    colouring the outside of a structure by a gas property would be inventing
    a quantity that does not exist there.

  * A slice plane through the axis, carrying the plume field as a texture.
    The axisymmetric solution lives on (x, r), so the slice is not a
    reconstruction -- it is the solution, shown where it was computed.  The
    shock cells sit on the axis and a surface rendering would hide them.

COORDINATES
-----------
Model space is x along the axis (downstream positive), y and z radial.  The
camera orbits the axis; because the flow is axisymmetric the azimuth only
changes which side the cut wedge faces, while elevation genuinely changes the
view.
"""
from __future__ import annotations

from dataclasses import dataclass, field as dc_field
from typing import Optional, Tuple

import numpy as np
from matplotlib import colormaps
from PySide6 import QtCore, QtGui, QtWidgets

from . import styles


@dataclass
class Mesh:
    """Triangles with per-vertex colour, and optionally a texture."""

    vertices: np.ndarray                      # (n, 3) float
    faces: np.ndarray                         # (m, 3) int
    colors: Optional[np.ndarray] = None       # (n, 3) float 0..1, per vertex
    uv: Optional[np.ndarray] = None           # (n, 2) float 0..1
    texture: Optional[np.ndarray] = None      # (h, w, 3) float 0..1
    # Per-texel coverage, same (h, w) as the texture.  Fragments below
    # `alpha_cutoff` are discarded outright rather than blended: the plume
    # slice is mostly ambient air, and drawing that as an opaque black
    # rectangle would hide everything behind it.  A cutout needs no sorting
    # and no blend state, which is the whole reason to prefer it here.
    alpha: Optional[np.ndarray] = None
    alpha_cutoff: float = 0.02
    emissive: bool = False                    # skip the lighting term
    name: str = ""
    # Per-vertex normals from the source file.  Without them they are worked
    # out from the faces, which is right for a revolved surface and wrong for
    # a model whose author chose where its edges are sharp.
    normals: Optional[np.ndarray] = None
    # Cut the cutaway wedge out of this mesh when drawing it, instead of the
    # wedge being left out of its geometry.  For models that are not built by
    # revolving a profile, so have no seam at the cut to stop on.
    cut_in_shader: bool = False


def revolve(x: np.ndarray, r: np.ndarray, *, theta0: float, theta1: float,
            n_theta: int, scalar: Optional[np.ndarray] = None,
            ramp: str = "inferno", lo: float = 0.0, hi: float = 1.0,
            flat_colour: Optional[Tuple[float, float, float]] = None,
            flip_normals: bool = False) -> Mesh:
    """Sweep a profile around the x axis into a triangle mesh.

    `scalar` colours the surface through `ramp`; `flat_colour` overrides it for
    a surface that carries no field, such as the outside of the structure.
    """
    theta = np.linspace(theta0, theta1, n_theta)
    ct, st = np.cos(theta), np.sin(theta)
    nx = x.size

    verts = np.empty((n_theta * nx, 3))
    verts[:, 0] = np.tile(x, n_theta)
    verts[:, 1] = np.repeat(ct, nx) * np.tile(r, n_theta)
    verts[:, 2] = np.repeat(st, nx) * np.tile(r, n_theta)

    if flat_colour is not None:
        colors = np.tile(np.asarray(flat_colour, dtype=float), (verts.shape[0], 1))
    else:
        values = np.zeros(nx) if scalar is None else np.asarray(scalar, dtype=float)
        span = hi - lo
        norm = np.clip((values - lo) / span, 0.0, 1.0) if span > 0 else np.zeros(nx)
        row = colormaps[ramp](norm)[:, :3]
        colors = np.tile(row, (n_theta, 1))

    # Two triangles per quad between adjacent profile rows.
    i = np.arange(nx - 1)
    j = np.arange(n_theta - 1)
    jj, ii = np.meshgrid(j, i, indexing="ij")
    a = (jj * nx + ii).ravel()
    b = a + 1
    c = a + nx
    d = c + 1
    if flip_normals:
        faces = np.concatenate([np.stack([a, c, b], axis=1),
                                np.stack([b, c, d], axis=1)])
    else:
        faces = np.concatenate([np.stack([a, b, c], axis=1),
                                np.stack([b, d, c], axis=1)])
    return Mesh(vertices=verts, faces=faces, colors=colors)


def slice_plane(x0: float, x1: float, r_max: float, texture: np.ndarray,
                *, azimuth: float = np.pi / 2) -> Mesh:
    """A quad through the axis carrying a field image as a texture.

    The default azimuth puts the plane's radial direction along z -- upright.
    A plane laid flat in x-y is a floor when seen from above and vanishes to a
    line when seen from the side; upright it reads as the cross-section it is.
    """
    c, s = np.cos(azimuth), np.sin(azimuth)
    verts = np.array([
        [x0, -r_max * c, -r_max * s],
        [x1, -r_max * c, -r_max * s],
        [x1, r_max * c, r_max * s],
        [x0, r_max * c, r_max * s],
    ])
    uv = np.array([[0.0, 1.0], [1.0, 1.0], [1.0, 0.0], [0.0, 0.0]])
    faces = np.array([[0, 1, 2], [0, 2, 3]])
    return Mesh(vertices=verts, faces=faces, uv=uv, texture=texture,
                emissive=True, name="slice")


@dataclass
class Camera3D:
    """An orbit camera: azimuth and elevation about a target, plus distance."""

    target: np.ndarray = dc_field(default_factory=lambda: np.zeros(3))
    distance: float = 3.0
    # Side-on and a little above, which is how an engine is drawn and the only
    # view where the whole contour reads at once.  Looking down the axis puts
    # the camera inside the bell's own shadow and shows a disc.
    #
    # The sign matters: at +pi/2 the screen-right vector comes out as -x, so
    # the engine faces left and the exhaust runs backwards across the frame.
    azimuth: float = -np.pi / 2
    elevation: float = 0.28     # rad
    fov: float = 0.9            # rad, vertical

    def basis(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        ce, se = np.cos(self.elevation), np.sin(self.elevation)
        ca, sa = np.cos(self.azimuth), np.sin(self.azimuth)
        forward = np.array([-ce * ca, -ce * sa, -se])
        eye = self.target - forward * self.distance
        world_up = np.array([0.0, 0.0, 1.0])
        right = np.cross(forward, world_up)
        n = np.linalg.norm(right)
        right = np.array([0.0, 1.0, 0.0]) if n < 1e-9 else right / n
        up = np.cross(right, forward)
        return eye, forward, right, up

    def view(self, points: np.ndarray) -> np.ndarray:
        """World points into camera space: +x right, +y up, +z forward."""
        eye, forward, right, up = self.basis()
        rel = points - eye
        return np.stack([rel @ right, rel @ up, rel @ forward], axis=1)


class Rasterizer:
    """Z-buffered triangle rasteriser, batched over triangles.

    The obvious implementation -- one Python-level pass per triangle with the
    per-pixel work vectorised inside it -- is what this started as, and it ran
    a 115 000-triangle engine at 0.2 fps.  The per-triangle interpreter
    overhead dominated completely.

    So triangles are rasterised *together* instead.  Every triangle small
    enough gets the same fixed TILE x TILE stamp, which makes the barycentric
    test one array operation over (n_triangles, TILE, TILE); the few that are
    larger are split off and handled on their own.  The depth test is then a
    scatter-minimum over packed keys -- the quantised depth in the high bits
    and the fragment index in the low bits -- so a single `np.minimum.at`
    resolves the whole buffer and the winning fragment index falls out of the
    low bits.  That is the standard trick for a z-buffer without a GPU, and it
    is what makes the camera orbit instead of crawl.
    """

    TILE = 12

    def __init__(self, width: int, height: int, background: Tuple[float, float, float]):
        self.width = width
        self.height = height
        self.color = np.tile(np.asarray(background, dtype=np.float32),
                             (height, width, 1))
        # Packed (depth, fragment) keys; the sentinel is "nothing here yet".
        self._key = np.full(height * width, np.iinfo(np.int64).max, dtype=np.int64)
        self._frag_rgb: list = []
        self._frag_offset = 0

    def draw(self, mesh: Mesh, camera: Camera3D,
             light: np.ndarray = np.array([0.4, 0.7, 0.55])) -> None:
        cam = camera.view(mesh.vertices)
        z = cam[:, 2]
        near = 1e-4
        scale = (self.height * 0.5) / np.tan(camera.fov * 0.5)
        with np.errstate(divide="ignore", invalid="ignore"):
            sx = self.width * 0.5 + cam[:, 0] / z * scale
            sy = self.height * 0.5 - cam[:, 1] / z * scale

        faces = mesh.faces
        # Triangles with any vertex behind the eye are dropped rather than
        # clipped: the scene is always viewed from outside its own bounds, so
        # a correct near-plane clip would cost more than it can ever save.
        faces = faces[(z[faces] > near).all(axis=1)]
        if faces.size == 0:
            return
        v0, v1, v2 = faces[:, 0], faces[:, 1], faces[:, 2]

        area = ((sx[v1] - sx[v0]) * (sy[v2] - sy[v0]) -
                (sx[v2] - sx[v0]) * (sy[v1] - sy[v0]))
        faces = faces[np.abs(area) > 1e-9]
        if faces.size == 0:
            return
        v0, v1, v2 = faces[:, 0], faces[:, 1], faces[:, 2]

        shade = np.ones(faces.shape[0], dtype=np.float32)
        if not mesh.emissive:
            e0 = mesh.vertices[v1] - mesh.vertices[v0]
            e1 = mesh.vertices[v2] - mesh.vertices[v0]
            normal = np.cross(e0, e1)
            normal /= np.maximum(np.linalg.norm(normal, axis=1, keepdims=True), 1e-12)
            lam = np.abs(normal @ (light / np.linalg.norm(light)))
            shade = (0.30 + 0.70 * lam).astype(np.float32)

        x0, x1, x2 = sx[v0], sx[v1], sx[v2]
        y0, y1, y2 = sy[v0], sy[v1], sy[v2]
        min_x = np.floor(np.minimum(np.minimum(x0, x1), x2))
        min_y = np.floor(np.minimum(np.minimum(y0, y1), y2))
        max_x = np.ceil(np.maximum(np.maximum(x0, x1), x2))
        max_y = np.ceil(np.maximum(np.maximum(y0, y1), y2))
        small = ((max_x - min_x) <= self.TILE) & ((max_y - min_y) <= self.TILE)

        if small.any():
            self._batch(mesh, faces[small], sx, sy, z, shade[small],
                        min_x[small], min_y[small], self.TILE)
        if (~small).any():
            # A handful of large triangles -- the slice plane, mostly.  Each
            # gets its own stamp, sized to it.
            big = np.flatnonzero(~small)
            for t in big:
                w = int(max_x[t] - min_x[t]) + 1
                h = int(max_y[t] - min_y[t]) + 1
                if w <= 0 or h <= 0 or w * h > 4_000_000:
                    continue
                self._batch(mesh, faces[t:t + 1] if False else faces[[t]],
                            sx, sy, z, shade[[t]], min_x[[t]], min_y[[t]],
                            max(w, h))

    def _batch(self, mesh: Mesh, faces, sx, sy, z, shade, min_x, min_y,
               tile: int) -> None:
        n = faces.shape[0]
        i0, i1, i2 = faces[:, 0], faces[:, 1], faces[:, 2]
        off = np.arange(tile, dtype=np.float32)
        gx = min_x[:, None, None] + off[None, None, :] + 0.5
        gy = min_y[:, None, None] + off[None, :, None] + 0.5

        x0, x1, x2 = sx[i0][:, None, None], sx[i1][:, None, None], sx[i2][:, None, None]
        y0, y1, y2 = sy[i0][:, None, None], sy[i1][:, None, None], sy[i2][:, None, None]
        denom = (y1 - y2) * (x0 - x2) + (x2 - x1) * (y0 - y2)
        denom = np.where(np.abs(denom) < 1e-12, np.nan, denom)
        w0 = ((y1 - y2) * (gx - x2) + (x2 - x1) * (gy - y2)) / denom
        w1 = ((y2 - y0) * (gx - x2) + (x0 - x2) * (gy - y2)) / denom
        w2 = 1.0 - w0 - w1

        # gx is (n, 1, tile) and gy is (n, tile, 1); the barycentrics broadcast
        # to (n, tile, tile), so the pixel indices have to as well or the flat
        # selection below indexes the wrong array.
        shape = (n, tile, tile)
        px = np.broadcast_to(np.rint(gx - 0.5).astype(np.int64), shape)
        py = np.broadcast_to(np.rint(gy - 0.5).astype(np.int64), shape)
        inside = ((w0 >= -1e-6) & (w1 >= -1e-6) & (w2 >= -1e-6) &
                  (px >= 0) & (px < self.width) & (py >= 0) & (py < self.height))
        if not inside.any():
            return

        z0, z1, z2 = z[i0][:, None, None], z[i1][:, None, None], z[i2][:, None, None]
        iz = w0 / z0 + w1 / z1 + w2 / z2
        with np.errstate(divide="ignore", invalid="ignore"):
            depth = 1.0 / iz
        inside &= np.isfinite(depth) & (depth > 0)
        if not inside.any():
            return

        sel = np.flatnonzero(inside.ravel())
        d = np.broadcast_to(depth, shape).ravel()[sel]
        pix = (py.ravel()[sel] * self.width + px.ravel()[sel])

        b0 = np.broadcast_to((w0 / z0) / iz, shape).ravel()[sel]
        b1 = np.broadcast_to((w1 / z1) / iz, shape).ravel()[sel]
        b2 = 1.0 - b0 - b1
        tri = np.repeat(np.arange(n), tile * tile)[sel]

        if mesh.texture is not None and mesh.uv is not None:
            u = b0 * mesh.uv[i0[tri], 0] + b1 * mesh.uv[i1[tri], 0] + b2 * mesh.uv[i2[tri], 0]
            v = b0 * mesh.uv[i0[tri], 1] + b1 * mesh.uv[i1[tri], 1] + b2 * mesh.uv[i2[tri], 1]
            th, tw = mesh.texture.shape[:2]
            ui = np.clip((u * (tw - 1)).astype(np.int32), 0, tw - 1)
            vi = np.clip((v * (th - 1)).astype(np.int32), 0, th - 1)
            rgb = mesh.texture[vi, ui]
            if mesh.alpha is not None:
                visible = mesh.alpha[vi, ui] >= mesh.alpha_cutoff
                if not visible.any():
                    return
                sel, d, pix, rgb = (sel[visible], d[visible],
                                    pix[visible], rgb[visible])
        else:
            c = mesh.colors
            rgb = (b0[:, None] * c[i0[tri]] + b1[:, None] * c[i1[tri]] +
                   b2[:, None] * c[i2[tri]])
        if not mesh.emissive:
            rgb = rgb * shade[tri][:, None]

        # Pack depth into the high bits so a scatter-minimum over the key
        # resolves both "which fragment is nearest" and "where is its colour"
        # in one pass.  24 bits of depth is ~16 million levels across the
        # scene, far finer than anything the eye or the geometry needs.
        far = float(np.max(d)) * 1.0000001 + 1e-9
        q = np.clip((d / far) * ((1 << 24) - 1), 0, (1 << 24) - 1).astype(np.int64)
        idx = np.arange(sel.size, dtype=np.int64) + self._frag_offset
        key = (q << 32) | idx
        np.minimum.at(self._key, pix, key)

        self._frag_rgb.append(np.ascontiguousarray(rgb, dtype=np.float32))
        self._frag_offset += sel.size

    def image(self) -> np.ndarray:
        """Resolve the key buffer into pixels.

        Nothing is written to the colour buffer until here, so a fragment that
        loses the depth test costs a comparison and never a write.
        """
        if self._frag_rgb:
            rgb = np.concatenate(self._frag_rgb, axis=0)
            hit = self._key != np.iinfo(np.int64).max
            if hit.any():
                frag = (self._key[hit] & 0xFFFFFFFF).astype(np.int64)
                flat = self.color.reshape(-1, 3)
                flat[hit] = rgb[frag]
        return np.clip(self.color, 0.0, 1.0)


class Viewport3D(QtWidgets.QWidget):
    """The 3-D pane: orbit with the left button, zoom with the wheel."""

    def __init__(self) -> None:
        super().__init__()
        self.setObjectName("viewport3d")
        self.setMinimumSize(360, 260)
        self.setSizePolicy(QtWidgets.QSizePolicy.Expanding,
                           QtWidgets.QSizePolicy.Expanding)
        self.setMouseTracking(True)
        self.camera = Camera3D()
        self._meshes: list = []
        self._image: Optional[QtGui.QImage] = None
        self._buffer: Optional[np.ndarray] = None
        self._drag: Optional[QtCore.QPoint] = None
        self._scale = 2          # render at 1/2 resolution, then upscale
        self._dirty = True
        self._caption = ""

    # -------------------------------------------------------------- contents
    def set_scene(self, meshes, *, target: np.ndarray, span: float,
                  caption: str = "") -> None:
        self._meshes = list(meshes)
        self.camera.target = np.asarray(target, dtype=float)
        if self._dirty:
            self.camera.distance = span * 1.9
            self._dirty = False
        self._caption = caption
        self._render()

    def frame(self, span: Optional[float] = None) -> None:
        if span is not None:
            self.camera.distance = span * 1.9
        self._render()

    # ---------------------------------------------------------------- render
    def _render(self) -> None:
        w = max(80, self.width() // self._scale)
        h = max(60, self.height() // self._scale)
        bg = QtGui.QColor(styles.BG_DEEP)
        raster = Rasterizer(w, h, (bg.redF(), bg.greenF(), bg.blueF()))
        cam = Camera3D(target=self.camera.target, distance=self.camera.distance,
                       azimuth=self.camera.azimuth,
                       elevation=self.camera.elevation, fov=self.camera.fov)
        cam.fov = self.camera.fov
        for mesh in self._meshes:
            raster.draw(mesh, cam)
        img = (raster.image() * 255.0 + 0.5).astype(np.uint8)
        self._buffer = np.ascontiguousarray(img)
        self._image = QtGui.QImage(self._buffer.data, w, h, 3 * w,
                                   QtGui.QImage.Format_RGB888)
        self.update()

    def paintEvent(self, event) -> None:      # noqa: N802  (Qt override)
        p = QtGui.QPainter(self)
        p.fillRect(self.rect(), QtGui.QColor(styles.BG_DEEP))
        if self._image is not None:
            p.setRenderHint(QtGui.QPainter.SmoothPixmapTransform, True)
            p.drawImage(self.rect(), self._image)
        if self._caption:
            p.setPen(QtGui.QColor(styles.TEXT_MUTED))
            p.drawText(self.rect().adjusted(10, 8, -10, -8),
                       QtCore.Qt.AlignTop | QtCore.Qt.AlignLeft, self._caption)
        else:
            p.setPen(QtGui.QColor(styles.TEXT_DIM))
            p.drawText(self.rect(), QtCore.Qt.AlignCenter,
                       "Solve a design to build the model.")
        p.end()

    def resizeEvent(self, event) -> None:     # noqa: N802  (Qt override)
        super().resizeEvent(event)
        if self._meshes:
            self._render()

    # ----------------------------------------------------------- interaction
    def mousePressEvent(self, event) -> None:     # noqa: N802
        if event.button() == QtCore.Qt.LeftButton:
            self._drag = event.position().toPoint()

    def mouseMoveEvent(self, event) -> None:      # noqa: N802
        if self._drag is None or not self._meshes:
            return
        pos = event.position().toPoint()
        dx = pos.x() - self._drag.x()
        dy = pos.y() - self._drag.y()
        self._drag = pos
        self.camera.azimuth += dx * 0.01
        limit = np.pi / 2 - 0.02
        self.camera.elevation = float(
            np.clip(self.camera.elevation + dy * 0.01, -limit, limit))
        self._render()

    def mouseReleaseEvent(self, event) -> None:   # noqa: N802
        self._drag = None

    def wheelEvent(self, event) -> None:          # noqa: N802
        if not self._meshes:
            return
        step = event.angleDelta().y() / 120.0
        self.camera.distance *= float(np.exp(-step * 0.12))
        self._render()


def engine_meshes(contour: np.ndarray, *, wall: float,
                  scalar: Optional[np.ndarray] = None, ramp: str = "inferno",
                  lo: float = 0.0, hi: float = 1.0,
                  cut_deg: float = 80.0, cut_azimuth: float = np.pi,
                  wall_colour: Optional[Tuple[float, float, float]] = None,
                  outer_colour: Tuple[float, float, float] = (0.30, 0.33, 0.38),
                  n_theta: int = 64):
    """The engine as an inner gas-side surface, an outer shell, and end caps.

    A wedge of `cut_deg` is left out so the camera can see the gas-side wall
    on the far side -- the standard cutaway, and the only way a surface
    rendering of a duct shows anything but its outside.
    """
    x, r = contour[:, 0], contour[:, 1]
    half = np.radians(cut_deg) / 2.0
    # The gap has to face the camera or the cutaway shows nothing: the default
    # view sits on -y, and y = r cos(theta), so the missing wedge is centred on
    # theta = pi.
    theta0 = cut_azimuth + half
    theta1 = cut_azimuth + 2.0 * np.pi - half

    inner = revolve(x, r, theta0=theta0, theta1=theta1, n_theta=n_theta,
                    scalar=scalar, ramp=ramp, lo=lo, hi=hi,
                    flat_colour=wall_colour, flip_normals=True)
    inner.name = "wall (gas side)"
    outer = revolve(x, r + wall, theta0=theta0, theta1=theta1, n_theta=n_theta,
                    flat_colour=outer_colour)
    outer.name = "structure"

    # The two faces exposed by the cut, so the wall reads as having thickness
    # rather than as an infinitely thin sheet.
    caps = []
    for ang in (theta0, theta1):
        # theta0 and theta1 are the two edges the sweep stops at.
        c, s = np.cos(ang), np.sin(ang)
        n = x.size
        verts = np.empty((2 * n, 3))
        verts[:n, 0] = x
        verts[:n, 1] = r * c
        verts[:n, 2] = r * s
        verts[n:, 0] = x
        verts[n:, 1] = (r + wall) * c
        verts[n:, 2] = (r + wall) * s
        i = np.arange(n - 1)
        faces = np.concatenate([
            np.stack([i, i + 1, i + n], axis=1),
            np.stack([i + 1, i + 1 + n, i + n], axis=1)])
        colors = np.tile(np.array([0.46, 0.50, 0.57]), (2 * n, 1))
        caps.append(Mesh(vertices=verts, faces=faces, colors=colors,
                         name="cut face"))
    return [outer, inner] + caps


def plume_slice(field, texture: np.ndarray, x_offset: float,
                alpha: Optional[np.ndarray] = None) -> Mesh:
    """The plume field as a slice plane, starting at the nozzle exit.

    `alpha` is normally the jet fraction, so the ambient air the solver also
    carries is simply not drawn.  Without it the slice is a rectangle of
    background colour that hides whatever is behind it.
    """
    dx = float(field.x[1] - field.x[0]) if field.x.size > 1 else 0.0
    dr = float(field.r[1] - field.r[0]) if field.r.size > 1 else 0.0
    length = float(field.x[-1]) + 0.5 * dx
    r_max = float(field.r[-1]) + 0.5 * dr
    mesh = slice_plane(x_offset, x_offset + length, r_max, texture)
    mesh.alpha = alpha
    return mesh


__all__ = ["Camera3D", "Mesh", "Rasterizer", "Viewport3D", "engine_meshes",
           "plume_slice", "revolve", "slice_plane"]
