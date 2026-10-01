"""A GPU viewport: the engine as lit geometry, the flow as a true 3-D volume.

WHY THIS REPLACED THE NUMPY RASTERISER
--------------------------------------
The NumPy path drew a slice through the plume -- a flat plane -- and did it on
the CPU, at a fraction of the window's resolution and a few frames a second.
It looked like a picture of a slice because it was one.

This draws the flow as a volume.  That is not a different simulation, and it
is worth being exact about why it is legitimate.  The plume solution is
axisymmetric: it lives on (x, r), and the full three-dimensional field is
exactly that field swept around the axis.  So a ray fired through the scene,
sampled at points (x, y, z), reads the solution at (x, sqrt(y^2 + z^2)) -- the
real field at that point in space, not an extrusion or an illustration of it.
The same holds inside the engine, where the quasi-1-D solution is by
definition uniform across each station: the gas flowing out of the chamber is
drawn from the axial profile the solver produced.

WHAT THE VOLUME SHOWS
---------------------
Colour is the selected field through the selected ramp, exactly as on the
colour bar.  Density -- how much a stretch of ray contributes -- is the jet
fraction the plume solver carries, so ambient air is transparent and exhaust
is not.  That is a visualisation choice, the standard emission-absorption
model a CFD post-processor uses for a scalar field, and it is not a radiation
calculation; the hero render in `ignis_viz.render` is the place for that.

HOW IT IS DRAWN
---------------
Two passes.  The engine is drawn into an offscreen buffer that keeps its depth
as a texture.  A full-screen pass then fires a ray per pixel, intersects the
bounding cylinder of the flow, marches through it, and stops at the engine's
depth -- so the metal occludes gas behind it and gas in front of the metal
glows over it, without the two ever being sorted by hand.

Everything uses OpenGL 3.3 core, which any GPU from the last decade and Mesa's
software renderer both provide.  If a context cannot be had, the widget says
so and the Explorer falls back to the NumPy viewport rather than showing a
black rectangle.
"""
from __future__ import annotations

import ctypes
from typing import Optional

import numpy as np
from matplotlib import colormaps
from PySide6 import QtCore, QtGui
from PySide6.QtOpenGLWidgets import QOpenGLWidget

try:
    from OpenGL import GL
    HAVE_PYOPENGL = True
except Exception:                                   # noqa: BLE001
    GL = None
    HAVE_PYOPENGL = False

from . import styles

# --- shaders ----------------------------------------------------------------

MESH_VS = """
#version 330 core
layout(location = 0) in vec3 aPos;
layout(location = 1) in vec3 aNrm;
layout(location = 2) in vec3 aCol;
uniform mat4 uViewProj;
out vec3 vPos;
out vec3 vNrm;
out vec3 vCol;
void main() {
    vPos = aPos;
    vNrm = aNrm;
    vCol = aCol;
    gl_Position = uViewProj * vec4(aPos, 1.0);
}
"""

MESH_FS = """
#version 330 core
in vec3 vPos;
in vec3 vNrm;
in vec3 vCol;
uniform vec3 uEye;
uniform vec3 uKey;
uniform vec3 uFill;
out vec4 fragColor;
void main() {
    vec3 n = normalize(vNrm);
    vec3 v = normalize(uEye - vPos);
    // The cutaway shows inner and outer faces alike, so light both sides.
    if (dot(n, v) < 0.0) n = -n;
    float key  = max(dot(n, uKey), 0.0);
    float fill = max(dot(n, uFill), 0.0);
    vec3 h = normalize(uKey + v);
    float spec = pow(max(dot(n, h), 0.0), 60.0);
    float rim = pow(1.0 - max(dot(n, v), 0.0), 3.0);
    vec3 c = vCol * (0.22 + 0.62 * key + 0.22 * fill) + vec3(0.45) * spec
           + vec3(0.10) * rim;
    fragColor = vec4(c, 1.0);
}
"""

QUAD_VS = """
#version 330 core
out vec2 vNdc;
void main() {
    // One triangle that covers the screen; no vertex buffer needed.
    vec2 p = vec2((gl_VertexID << 1) & 2, gl_VertexID & 2);
    vNdc = p * 2.0 - 1.0;
    gl_Position = vec4(vNdc, 0.0, 1.0);
}
"""

VOLUME_FS = """
#version 330 core
in vec2 vNdc;
uniform sampler2D uSceneColor;
uniform sampler2D uSceneDepth;
uniform sampler2D uPlume;     // (x, r) -> R: normalised field, G: jet fraction
uniform sampler2D uInside;    // x -> R: normalised field, G: wall radius
uniform sampler2D uRamp;      // 256 x 1 colour ramp
uniform mat4 uInvViewProj;
uniform vec3 uEye;
uniform vec3 uBackground;
uniform float uPlumeX0;
uniform float uPlumeX1;
uniform float uPlumeR;
uniform float uInsideX0;
uniform float uInsideX1;
uniform float uInsideR;
uniform int uHasPlume;
uniform int uHasInside;
uniform float uDensity;        // plume, per unit length
uniform float uInsideDensity;  // inside the engine, per unit length
uniform float uGain;
uniform int uSteps;
out vec4 fragColor;

bool cylinder(vec3 o, vec3 d, float R, float xa, float xb,
              out float t0, out float t1) {
    float tx0 = -1e30, tx1 = 1e30;
    if (abs(d.x) > 1e-9) {
        float a = (xa - o.x) / d.x, b = (xb - o.x) / d.x;
        tx0 = min(a, b); tx1 = max(a, b);
    } else if (o.x < xa || o.x > xb) {
        return false;
    }
    float A = d.y * d.y + d.z * d.z;
    float B = 2.0 * (o.y * d.y + o.z * d.z);
    float C = o.y * o.y + o.z * o.z - R * R;
    float tc0 = -1e30, tc1 = 1e30;
    if (A > 1e-12) {
        float disc = B * B - 4.0 * A * C;
        if (disc < 0.0) return false;
        float s = sqrt(disc);
        tc0 = (-B - s) / (2.0 * A);
        tc1 = (-B + s) / (2.0 * A);
    } else if (C > 0.0) {
        return false;
    }
    t0 = max(max(tx0, tc0), 0.0);
    t1 = min(tx1, tc1);
    return t1 > t0;
}

void main() {
    vec2 uv = vNdc * 0.5 + 0.5;
    vec3 base = texture(uSceneColor, uv).rgb;
    float depth = texture(uSceneDepth, uv).r;

    vec4 far = uInvViewProj * vec4(vNdc, 1.0, 1.0);
    far /= far.w;
    vec3 dir = normalize(far.xyz - uEye);
    float tScene = 1e30;
    if (depth < 1.0) {
        vec4 hit = uInvViewProj * vec4(vNdc, depth * 2.0 - 1.0, 1.0);
        hit /= hit.w;
        tScene = dot(hit.xyz - uEye, dir);
    }

    float xa = (uHasInside == 1) ? uInsideX0 : uPlumeX0;
    float xb = (uHasPlume == 1) ? uPlumeX1 : uInsideX1;
    float R = max(uHasPlume == 1 ? uPlumeR : 0.0, uHasInside == 1 ? uInsideR : 0.0);
    float t0, t1;
    if ((uHasPlume == 0 && uHasInside == 0) || !cylinder(uEye, dir, R, xa, xb, t0, t1)) {
        fragColor = vec4(base, 1.0);
        return;
    }
    t1 = min(t1, tScene);
    if (t1 <= t0) {
        fragColor = vec4(base, 1.0);
        return;
    }

    float dt = (t1 - t0) / float(uSteps);
    vec3 acc = vec3(0.0);
    float trans = 1.0;
    for (int i = 0; i < uSteps; ++i) {
        float t = t0 + (float(i) + 0.5) * dt;
        vec3 p = uEye + dir * t;
        float r = length(p.yz);
        float value = 0.0, cover = 0.0, rho = uDensity;
        if (uHasPlume == 1 && p.x >= uPlumeX0 && p.x <= uPlumeX1 && r <= uPlumeR) {
            vec2 s = texture(uPlume, vec2((p.x - uPlumeX0) / (uPlumeX1 - uPlumeX0),
                                          r / uPlumeR)).rg;
            value = s.r;
            cover = s.g;
        } else if (uHasInside == 1 && p.x >= uInsideX0 && p.x < uInsideX1) {
            vec2 s = texture(uInside, vec2((p.x - uInsideX0) / (uInsideX1 - uInsideX0),
                                           0.5)).rg;
            // Quasi-1-D: the station's state fills the duct out to the wall.
            cover = (r < s.g) ? 1.0 : 0.0;
            value = s.r;
            rho = uInsideDensity;
        }
        if (cover < 0.004) continue;
        vec3 c = texture(uRamp, vec2(clamp(value, 0.0, 1.0), 0.5)).rgb;
        // Hot gas is made more opaque than cool gas.  With a flat opacity the
        // outer layers of the jet hide its core, and the shock cells -- which
        // sit on the axis -- vanish behind a uniform shell.
        float a = 1.0 - exp(-rho * cover * (0.12 + 0.88 * value) * dt);
        acc += trans * a * c * uGain;
        trans *= 1.0 - a;
        if (trans < 0.01) break;
    }
    fragColor = vec4(base * trans + acc, 1.0);
}
"""


# --- small matrix helpers ----------------------------------------------------

def perspective(fovy: float, aspect: float, near: float, far: float) -> np.ndarray:
    f = 1.0 / np.tan(fovy / 2.0)
    m = np.zeros((4, 4))
    m[0, 0] = f / aspect
    m[1, 1] = f
    m[2, 2] = (far + near) / (near - far)
    m[2, 3] = 2.0 * far * near / (near - far)
    m[3, 2] = -1.0
    return m


def look_at(eye: np.ndarray, target: np.ndarray, up: np.ndarray) -> np.ndarray:
    f = target - eye
    f = f / np.linalg.norm(f)
    s = np.cross(f, up)
    n = np.linalg.norm(s)
    s = np.array([0.0, 1.0, 0.0]) if n < 1e-9 else s / n
    u = np.cross(s, f)
    m = np.eye(4)
    m[0, :3], m[1, :3], m[2, :3] = s, u, -f
    m[:3, 3] = -m[:3, :3] @ eye
    return m


def vertex_normals(vertices: np.ndarray, faces: np.ndarray) -> np.ndarray:
    """Area-weighted vertex normals, so a smooth body shades smooth."""
    v0, v1, v2 = (vertices[faces[:, i]] for i in range(3))
    fn = np.cross(v1 - v0, v2 - v0)
    out = np.zeros_like(vertices)
    for i in range(3):
        np.add.at(out, faces[:, i], fn)
    length = np.linalg.norm(out, axis=1, keepdims=True)
    return out / np.maximum(length, 1e-12)


class GLViewport(QOpenGLWidget):
    """Orbit with the left button, pan with the right, zoom with the wheel."""

    #: emitted once, after the first attempt to initialise GL: True if usable
    ready = QtCore.Signal(bool)

    def __init__(self) -> None:
        super().__init__()
        fmt = QtGui.QSurfaceFormat()
        fmt.setVersion(3, 3)
        fmt.setProfile(QtGui.QSurfaceFormat.CoreProfile)
        fmt.setDepthBufferSize(24)
        fmt.setSamples(4)
        self.setFormat(fmt)
        self.setMinimumSize(360, 260)
        self.ok = False
        self.error = ""

        # Camera: orbit about a target on the axis.
        self.target = np.zeros(3)
        self.distance = 3.0
        # Upstream of the scene and to the side the cutaway faces, so the
        # first view is the engine with its exhaust running away from it.
        self.azimuth = -np.pi / 2 - 0.62
        self.elevation = 0.30
        self._drag: Optional[QtCore.QPointF] = None
        self._button = None

        # Scene inputs, held on the CPU and uploaded when they change.
        self._mesh: Optional[np.ndarray] = None
        self._mesh_dirty = False
        self._plume: Optional[np.ndarray] = None
        self._plume_box = (0.0, 1.0, 1.0)
        self._plume_dirty = False
        self._inside: Optional[np.ndarray] = None
        self._inside_box = (0.0, 1.0, 1.0)
        self._inside_dirty = False
        self._ramp_name = "inferno"
        self._ramp_dirty = True
        self.opacity = 5.0          # optical depth across the plume diameter
        self.inside_opacity = 5.0   # ... and across the chamber
        self.gain = 1.0
        self.steps = 160
        self.caption = ""

        self._fbo = None
        self._size = (0, 0)

    # ------------------------------------------------------------ scene API
    def set_engine(self, meshes) -> None:
        """Engine meshes (viewport3d.Mesh) with per-vertex colour."""
        parts = []
        for m in meshes:
            if m.colors is None:
                continue
            n = vertex_normals(m.vertices, m.faces)
            idx = m.faces.reshape(-1)
            parts.append(np.hstack([m.vertices[idx], n[idx], m.colors[idx]]))
        self._mesh = np.vstack(parts).astype(np.float32) if parts else None
        self._mesh_dirty = True
        self.update()

    def set_plume(self, value: Optional[np.ndarray], cover: Optional[np.ndarray],
                  x0: float, length: float, r_max: float) -> None:
        """Plume field on (r, x): value normalised to 0..1, cover the jet fraction."""
        if value is None:
            self._plume = None
        else:
            self._plume = np.ascontiguousarray(
                np.stack([value, cover], axis=-1), dtype=np.float32)
        self._plume_box = (float(x0), float(x0 + length), float(r_max))
        self._plume_dirty = True
        self.update()

    def set_inside(self, x: Optional[np.ndarray], value: Optional[np.ndarray],
                   r_wall: Optional[np.ndarray]) -> None:
        """The gas inside the engine, from the quasi-1-D axial profile."""
        if x is None:
            self._inside = None
        else:
            # Resample onto a uniform x grid: the texture is sampled linearly
            # in x, and the solver's stations are not uniformly spaced.
            grid = np.linspace(float(x[0]), float(x[-1]), 512)
            v = np.interp(grid, x, value)
            w = np.interp(grid, x, r_wall)
            self._inside = np.ascontiguousarray(
                np.stack([v, w], axis=-1)[None, :, :], dtype=np.float32)
            self._inside_box = (float(x[0]), float(x[-1]), float(np.max(r_wall)))
        self._inside_dirty = True
        self.update()

    def set_ramp(self, name: str) -> None:
        if name != self._ramp_name:
            self._ramp_name = name
            self._ramp_dirty = True
            self.update()

    def frame_scene(self, x_min: float, x_max: float, radius: float) -> None:
        self.target = np.array([x_min + 0.42 * (x_max - x_min), 0.0, 0.0])
        self.distance = 0.55 * (x_max - x_min) + 1.6 * radius
        self.update()

    # ----------------------------------------------------------- GL plumbing
    def initializeGL(self) -> None:          # noqa: N802  (Qt override)
        if not HAVE_PYOPENGL:
            self.error = "PyOpenGL is not installed"
            self.ready.emit(False)
            return
        try:
            self._mesh_prog = self._program(MESH_VS, MESH_FS)
            self._vol_prog = self._program(QUAD_VS, VOLUME_FS)
            self._mesh_vao = GL.glGenVertexArrays(1)
            self._mesh_vbo = GL.glGenBuffers(1)
            self._quad_vao = GL.glGenVertexArrays(1)
            self._tex_plume = self._texture()
            self._tex_inside = self._texture()
            self._tex_ramp = self._texture()
            self._mesh_count = 0
            self.ok = True
        except Exception as exc:                        # noqa: BLE001
            self.error = str(exc)
            self.ok = False
        self.ready.emit(self.ok)

    def resizeGL(self, w: int, h: int) -> None:     # noqa: N802  (Qt override)
        if self.ok:
            self._make_fbo(int(w * self.devicePixelRatioF()),
                           int(h * self.devicePixelRatioF()))

    def paintGL(self) -> None:               # noqa: N802  (Qt override)
        if not self.ok:
            return
        w = int(self.width() * self.devicePixelRatioF())
        h = int(self.height() * self.devicePixelRatioF())
        if (w, h) != self._size:
            self._make_fbo(w, h)
        self._upload()
        # The caption is drawn with QPainter, and Qt's text path leaves
        # blending on with glBlendFunc(GL_CONSTANT_COLOR, GL_ONE_MINUS_SRC_COLOR)
        # and the pen as the blend colour.  Inherited, that multiplies every
        # later frame into the caption's grey -- the first frame is clean and
        # every one after it is fogged.  So set what the passes rely on.
        GL.glDisable(GL.GL_BLEND)
        GL.glDisable(GL.GL_SCISSOR_TEST)
        GL.glDisable(GL.GL_STENCIL_TEST)
        GL.glDisable(GL.GL_CULL_FACE)
        GL.glDepthMask(GL.GL_TRUE)
        GL.glColorMask(GL.GL_TRUE, GL.GL_TRUE, GL.GL_TRUE, GL.GL_TRUE)

        eye, view_proj = self._matrices(w, h)
        bg = QtGui.QColor(styles.BG_DEEP)
        background = (bg.redF(), bg.greenF(), bg.blueF())

        # Pass 1: the engine, into a buffer whose depth is a texture.
        GL.glBindFramebuffer(GL.GL_FRAMEBUFFER, self._fbo)
        GL.glViewport(0, 0, w, h)
        GL.glClearColor(*background, 1.0)
        GL.glClear(GL.GL_COLOR_BUFFER_BIT | GL.GL_DEPTH_BUFFER_BIT)
        GL.glEnable(GL.GL_DEPTH_TEST)
        if self._mesh_count:
            GL.glUseProgram(self._mesh_prog)
            self._set_mat(self._mesh_prog, "uViewProj", view_proj)
            self._set_vec(self._mesh_prog, "uEye", eye)
            key = np.array([0.35, -0.55, 0.75])
            fill = np.array([-0.6, 0.4, 0.2])
            self._set_vec(self._mesh_prog, "uKey", key / np.linalg.norm(key))
            self._set_vec(self._mesh_prog, "uFill", fill / np.linalg.norm(fill))
            GL.glBindVertexArray(self._mesh_vao)
            GL.glDrawArrays(GL.GL_TRIANGLES, 0, self._mesh_count)

        # Pass 2: march the flow over it, into the widget's own framebuffer.
        GL.glBindFramebuffer(GL.GL_FRAMEBUFFER, self.defaultFramebufferObject())
        GL.glViewport(0, 0, w, h)
        GL.glDisable(GL.GL_DEPTH_TEST)
        GL.glUseProgram(self._vol_prog)
        p = self._vol_prog
        for unit, (name, tex) in enumerate([("uSceneColor", self._tex_color),
                                            ("uSceneDepth", self._tex_depth),
                                            ("uPlume", self._tex_plume),
                                            ("uInside", self._tex_inside),
                                            ("uRamp", self._tex_ramp)]):
            GL.glActiveTexture(GL.GL_TEXTURE0 + unit)
            GL.glBindTexture(GL.GL_TEXTURE_2D, tex)
            GL.glUniform1i(GL.glGetUniformLocation(p, name), unit)
        self._set_mat(p, "uInvViewProj", np.linalg.inv(view_proj))
        self._set_vec(p, "uEye", eye)
        self._set_vec(p, "uBackground", np.array(background))
        x0, x1, rp = self._plume_box
        xi0, xi1, ri = self._inside_box
        for name, value in [("uPlumeX0", x0), ("uPlumeX1", x1), ("uPlumeR", rp),
                            ("uInsideX0", xi0), ("uInsideX1", xi1), ("uInsideR", ri),
                            ("uDensity", self.opacity / max(2.0 * rp, 1e-6)),
                            ("uInsideDensity", self.inside_opacity / max(2.0 * ri, 1e-6)),
                            ("uGain", self.gain)]:
            GL.glUniform1f(GL.glGetUniformLocation(p, name), float(value))
        GL.glUniform1i(GL.glGetUniformLocation(p, "uHasPlume"),
                       int(self._plume is not None))
        GL.glUniform1i(GL.glGetUniformLocation(p, "uHasInside"),
                       int(self._inside is not None))
        GL.glUniform1i(GL.glGetUniformLocation(p, "uSteps"), int(self.steps))
        GL.glBindVertexArray(self._quad_vao)
        GL.glDrawArrays(GL.GL_TRIANGLES, 0, 3)
        GL.glActiveTexture(GL.GL_TEXTURE0)

        if self.caption:
            painter = QtGui.QPainter(self)
            painter.setPen(QtGui.QColor(styles.TEXT_MUTED))
            painter.drawText(self.rect().adjusted(10, 8, -10, -8),
                             QtCore.Qt.AlignTop | QtCore.Qt.AlignLeft, self.caption)
            painter.end()

    # -------------------------------------------------------------- helpers
    def _program(self, vs_src: str, fs_src: str) -> int:
        def compile_one(kind, src):
            shader = GL.glCreateShader(kind)
            GL.glShaderSource(shader, src)
            GL.glCompileShader(shader)
            if not GL.glGetShaderiv(shader, GL.GL_COMPILE_STATUS):
                raise RuntimeError(GL.glGetShaderInfoLog(shader).decode())
            return shader
        prog = GL.glCreateProgram()
        for kind, src in ((GL.GL_VERTEX_SHADER, vs_src), (GL.GL_FRAGMENT_SHADER, fs_src)):
            GL.glAttachShader(prog, compile_one(kind, src))
        GL.glLinkProgram(prog)
        if not GL.glGetProgramiv(prog, GL.GL_LINK_STATUS):
            raise RuntimeError(GL.glGetProgramInfoLog(prog).decode())
        return prog

    @staticmethod
    def _texture() -> int:
        tex = GL.glGenTextures(1)
        GL.glBindTexture(GL.GL_TEXTURE_2D, tex)
        GL.glTexParameteri(GL.GL_TEXTURE_2D, GL.GL_TEXTURE_MIN_FILTER, GL.GL_LINEAR)
        GL.glTexParameteri(GL.GL_TEXTURE_2D, GL.GL_TEXTURE_MAG_FILTER, GL.GL_LINEAR)
        GL.glTexParameteri(GL.GL_TEXTURE_2D, GL.GL_TEXTURE_WRAP_S, GL.GL_CLAMP_TO_EDGE)
        GL.glTexParameteri(GL.GL_TEXTURE_2D, GL.GL_TEXTURE_WRAP_T, GL.GL_CLAMP_TO_EDGE)
        return tex

    def _make_fbo(self, w: int, h: int) -> None:
        w, h = max(w, 1), max(h, 1)
        if self._fbo is not None:
            GL.glDeleteFramebuffers(1, [self._fbo])
            GL.glDeleteTextures(2, [self._tex_color, self._tex_depth])
        self._tex_color = GL.glGenTextures(1)
        GL.glBindTexture(GL.GL_TEXTURE_2D, self._tex_color)
        GL.glTexImage2D(GL.GL_TEXTURE_2D, 0, GL.GL_RGBA8, w, h, 0,
                        GL.GL_RGBA, GL.GL_UNSIGNED_BYTE, None)
        for pname in (GL.GL_TEXTURE_MIN_FILTER, GL.GL_TEXTURE_MAG_FILTER):
            GL.glTexParameteri(GL.GL_TEXTURE_2D, pname, GL.GL_NEAREST)
        self._tex_depth = GL.glGenTextures(1)
        GL.glBindTexture(GL.GL_TEXTURE_2D, self._tex_depth)
        GL.glTexImage2D(GL.GL_TEXTURE_2D, 0, GL.GL_DEPTH_COMPONENT24, w, h, 0,
                        GL.GL_DEPTH_COMPONENT, GL.GL_UNSIGNED_INT, None)
        for pname in (GL.GL_TEXTURE_MIN_FILTER, GL.GL_TEXTURE_MAG_FILTER):
            GL.glTexParameteri(GL.GL_TEXTURE_2D, pname, GL.GL_NEAREST)
        self._fbo = GL.glGenFramebuffers(1)
        GL.glBindFramebuffer(GL.GL_FRAMEBUFFER, self._fbo)
        GL.glFramebufferTexture2D(GL.GL_FRAMEBUFFER, GL.GL_COLOR_ATTACHMENT0,
                                  GL.GL_TEXTURE_2D, self._tex_color, 0)
        GL.glFramebufferTexture2D(GL.GL_FRAMEBUFFER, GL.GL_DEPTH_ATTACHMENT,
                                  GL.GL_TEXTURE_2D, self._tex_depth, 0)
        status = GL.glCheckFramebufferStatus(GL.GL_FRAMEBUFFER)
        if status != GL.GL_FRAMEBUFFER_COMPLETE:
            raise RuntimeError(f"offscreen buffer incomplete (status {status:#x})")
        GL.glBindFramebuffer(GL.GL_FRAMEBUFFER, self.defaultFramebufferObject())
        self._size = (w, h)

    def _upload(self) -> None:
        if self._mesh_dirty:
            self._mesh_dirty = False
            GL.glBindVertexArray(self._mesh_vao)
            GL.glBindBuffer(GL.GL_ARRAY_BUFFER, self._mesh_vbo)
            data = self._mesh if self._mesh is not None else np.zeros((0, 9), np.float32)
            GL.glBufferData(GL.GL_ARRAY_BUFFER, data.nbytes, data, GL.GL_STATIC_DRAW)
            stride = 9 * 4
            for loc in range(3):
                GL.glEnableVertexAttribArray(loc)
                GL.glVertexAttribPointer(loc, 3, GL.GL_FLOAT, False, stride,
                                         ctypes.c_void_p(loc * 12))
            self._mesh_count = int(data.shape[0])
        if self._plume_dirty and self._plume is not None:
            self._plume_dirty = False
            h, w = self._plume.shape[:2]
            GL.glBindTexture(GL.GL_TEXTURE_2D, self._tex_plume)
            GL.glTexImage2D(GL.GL_TEXTURE_2D, 0, GL.GL_RG32F, w, h, 0,
                            GL.GL_RG, GL.GL_FLOAT, self._plume)
        if self._inside_dirty and self._inside is not None:
            self._inside_dirty = False
            h, w = self._inside.shape[:2]
            GL.glBindTexture(GL.GL_TEXTURE_2D, self._tex_inside)
            GL.glTexImage2D(GL.GL_TEXTURE_2D, 0, GL.GL_RG32F, w, h, 0,
                            GL.GL_RG, GL.GL_FLOAT, self._inside)
        if self._ramp_dirty:
            self._ramp_dirty = False
            try:
                table = colormaps[self._ramp_name](np.linspace(0, 1, 256))[:, :3]
            except KeyError:
                table = colormaps["inferno"](np.linspace(0, 1, 256))[:, :3]
            data = np.ascontiguousarray(table[None, :, :], dtype=np.float32)
            GL.glBindTexture(GL.GL_TEXTURE_2D, self._tex_ramp)
            GL.glTexImage2D(GL.GL_TEXTURE_2D, 0, GL.GL_RGB32F, 256, 1, 0,
                            GL.GL_RGB, GL.GL_FLOAT, data)

    def _matrices(self, w: int, h: int):
        ce, se = np.cos(self.elevation), np.sin(self.elevation)
        ca, sa = np.cos(self.azimuth), np.sin(self.azimuth)
        eye = self.target + self.distance * np.array([ce * ca, ce * sa, se])
        view = look_at(eye, self.target, np.array([0.0, 0.0, 1.0]))
        proj = perspective(np.radians(40.0), w / max(h, 1),
                           self.distance * 0.01, self.distance * 20.0)
        return eye, proj @ view

    @staticmethod
    def _set_mat(prog, name, m):
        GL.glUniformMatrix4fv(GL.glGetUniformLocation(prog, name), 1, True,
                              np.ascontiguousarray(m, dtype=np.float32))

    @staticmethod
    def _set_vec(prog, name, v):
        GL.glUniform3f(GL.glGetUniformLocation(prog, name), *[float(c) for c in v])

    # ------------------------------------------------------------ interaction
    def mousePressEvent(self, event) -> None:      # noqa: N802
        self._drag = event.position()
        self._button = event.button()

    def mouseMoveEvent(self, event) -> None:       # noqa: N802
        if self._drag is None:
            return
        pos = event.position()
        dx, dy = pos.x() - self._drag.x(), pos.y() - self._drag.y()
        self._drag = pos
        if self._button == QtCore.Qt.RightButton:
            # Pan along the axis, which is the only direction worth panning.
            self.target = self.target + np.array([-dx * self.distance * 0.0015, 0, 0])
        else:
            self.azimuth -= dx * 0.008
            self.elevation = float(np.clip(self.elevation + dy * 0.008, -1.45, 1.45))
        self.update()

    def mouseReleaseEvent(self, event) -> None:    # noqa: N802
        self._drag = None

    def wheelEvent(self, event) -> None:           # noqa: N802
        self.distance *= float(np.exp(-event.angleDelta().y() / 120.0 * 0.12))
        self.update()


__all__ = ["GLViewport", "HAVE_PYOPENGL"]
