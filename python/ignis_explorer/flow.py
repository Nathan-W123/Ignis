"""Time-resolved plume view for the Ignis Engine Explorer.

WHAT THIS SHOWS
---------------
The nozzle solution stops at the exit plane, which is where a quasi-1D internal
model should stop.  This tab takes that exit state and marches the axisymmetric
Euler equations outward from it (`ignis_viz.plume`), starting from rest, so what
you watch is the plume establishing itself: the jet front driving into still
air, the starting vortex, the shock cells forming from the exit plane outward,
and the Mach disc settling into place.

Every frame is a solution of the equations at that instant.  Nothing is
interpolated between frames and nothing is drawn on.

WHAT IT IS NOT
--------------
The plume model is inviscid and axisymmetric.  It has no turbulence model, so
the shear layer spreads only by numerical diffusion and you will not see the
fine-grained turbulent breakup of a real jet -- that is three-dimensional and
needs a different code entirely.  What is real here is the shock structure and
the wave propagation.  `ignis_viz/plume.py` lists everything else the model
leaves out, and `docs/rendering.md` explains why the plume is computed rather
than painted.

COLOUR
------
Every field shown is a magnitude, so it gets a sequential ramp: one hue, dark
to light, perceptually ordered.  `inferno` is the default because it matches
the thermal aesthetic of the project's renders and is perceptually uniform, so
equal steps in the data are equal steps in apparent brightness.  `turbo` is
offered because it is what most published CFD figures use and people expect it,
but it is not the default: a rainbow ramp invents contours the data does not
have.
"""
from __future__ import annotations

import contextlib
import io
import os
import sys
from dataclasses import dataclass
from typing import List, Optional

import numpy as np
from matplotlib import colormaps
from PySide6 import QtCore, QtGui, QtWidgets

from . import shell, styles, viewport3d
from .solver import Result
from .widgets import panel

# Fields the viewer offers, in the order they appear in the selector.
FIELDS = {
    "Temperature": ("temperature", "K", "inferno"),
    "Mach number": ("mach", "-", "viridis"),
    "Pressure": ("pressure", "Pa", "magma"),
    "Jet fraction": ("jet_fraction", "-", "cividis"),
    "Axial velocity": ("velocity_x", "m/s", "inferno"),
}

RAMPS = ["inferno", "viridis", "magma", "cividis", "plasma", "turbo"]

# (nx, nr, length, width, max steps, snapshot cadence, rough wall time).
# Length and width are in exit radii.  The domain is long and narrow because a
# plume is: a square window spends most of its pixels on ambient air.
#
# The first row's time is measured -- 9000 steps on a 300x90 grid took 561 s on
# the machine this was written on.  The other two are that measurement scaled
# by the product of cell count and step count, which is what an explicit
# finite-volume march costs.  They are labelled "~" because neither the
# measurement nor the scaling transfers to another machine.
#
# Nothing here is a wait in the usual sense: frames start arriving a few
# seconds in and the plume is watchable while it builds, so the number is how
# long until it stops changing, not how long until there is something to see.
GRIDS = [
    (300, 90, 12.0, 4.5, 9000, 60, "9 min"),
    (440, 130, 14.0, 4.5, 13000, 90, "28 min"),
    (640, 190, 16.0, 5.0, 18000, 130, "85 min"),
]


@dataclass
class Frame:
    """One instant of the march, kept small enough to hold a few hundred."""

    step: int
    time: float                      # s since the nozzle started flowing
    data: dict                       # field name -> float32 (nr, nx)
    extent: tuple                    # (length, full mirrored height), m

    @classmethod
    def capture(cls, step: int, time: float, field) -> "Frame":
        dx = float(field.x[1] - field.x[0]) if field.x.size > 1 else 0.0
        dr = float(field.r[1] - field.r[0]) if field.r.size > 1 else 0.0
        return cls(step=step, time=time,
                   extent=(float(field.x[-1]) + 0.5 * dx,
                           2.0 * (float(field.r[-1]) + 0.5 * dr)),
                   data={key: np.asarray(getattr(field, key), dtype=np.float32)
                         for key, _, _ in FIELDS.values()})


class PlumeWorker(QtCore.QObject):
    """Marches the plume on a worker thread, handing back frames as they land."""

    frame = QtCore.Signal(object)
    finished = QtCore.Signal(object)     # final PlumeField, or None on failure
    failed = QtCore.Signal(str)
    progress = QtCore.Signal(str)

    def __init__(self, exit_state, ambient: float, *, nx: int, nr: int,
                 length: float, width: float, max_steps: int,
                 every: int) -> None:
        super().__init__()
        self._exit = exit_state
        self._ambient = ambient
        self._nx, self._nr = nx, nr
        self._length, self._width = length, width
        self._max_steps = max_steps
        self._every = every
        self._stop = False

    @QtCore.Slot()
    def stop(self) -> None:
        self._stop = True

    @QtCore.Slot()
    def run(self) -> None:
        from ignis_viz import plume as plume_mod

        captured: List[Frame] = []

        def snapshot(step, time, field):
            if self._stop:
                raise _Cancelled()
            frame = Frame.capture(step, time, field)
            captured.append(frame)
            self.frame.emit(frame)
            self.progress.emit(
                f"step {step}   t = {time * 1e3:.3f} ms   "
                f"{len(captured)} frames")

        try:
            final = plume_mod.solve(
                self._exit, self._ambient, nx=self._nx, nr=self._nr,
                length=self._length, width=self._width,
                max_steps=self._max_steps, start="rest",
                snapshot=snapshot, snapshot_every=self._every)
        except _Cancelled:
            self.progress.emit(f"stopped after {len(captured)} frames")
            self.finished.emit(None)
            return
        except Exception as exc:                       # noqa: BLE001
            self.failed.emit(str(exc))
            self.finished.emit(None)
            return
        self.finished.emit(final)


class _Cancelled(Exception):
    """Raised inside the snapshot callback to unwind a cancelled march."""


class FieldView(QtWidgets.QWidget):
    """Paints one scalar field as a colour image, mirrored about the axis.

    The solver carries only r >= 0 because the flow is axisymmetric.  Mirroring
    is therefore not decoration: the lower half is the same solution, and
    showing it is what makes the picture read as a jet rather than as half of
    one.
    """

    def __init__(self) -> None:
        super().__init__()
        self._image: Optional[QtGui.QImage] = None
        self._caption = ""
        self._extent = (1.0, 1.0)          # physical (length, full height), m
        self._contour: Optional[np.ndarray] = None   # (x, r) in m, exit at x=0
        self.setMinimumHeight(260)
        self.setSizePolicy(QtWidgets.QSizePolicy.Expanding,
                           QtWidgets.QSizePolicy.Expanding)

    def set_contour(self, contour: Optional[np.ndarray]) -> None:
        """Engine wall profile to draw upstream of the exit plane.

        Given as (n, 2) of (x, r) in metres with the exit plane at x = 0 and
        x negative inside the engine.  This is the contour Ignis actually
        solved, not a sketch.
        """
        self._contour = contour
        self.update()

    def show_field(self, values: np.ndarray, ramp: str, lo: float, hi: float,
                   caption: str, extent) -> None:
        span = hi - lo
        norm = (values - lo) / span if span > 0 else np.zeros_like(values)
        rgba = colormaps[ramp](np.clip(norm, 0.0, 1.0), bytes=True)
        mirrored = np.ascontiguousarray(np.concatenate([rgba[::-1], rgba], axis=0))
        h, w = mirrored.shape[:2]
        # Keep a reference: QImage does not copy the buffer it is handed.
        self._buffer = mirrored
        self._image = QtGui.QImage(mirrored.data, w, h, 4 * w,
                                   QtGui.QImage.Format_RGBA8888)
        self._caption = caption
        self._extent = extent
        self.update()

    def clear(self) -> None:
        self._image = None
        self._caption = ""
        self.update()

    def paintEvent(self, event) -> None:      # noqa: N802  (Qt override)
        painter = QtGui.QPainter(self)
        painter.fillRect(self.rect(), QtGui.QColor(styles.BG_DEEP))
        if self._image is None:
            painter.setPen(QtGui.QColor(styles.TEXT_DIM))
            painter.drawText(self.rect(), QtCore.Qt.AlignCenter,
                             "Solve a design, then press Run flow.")
            painter.end()
            return
        painter.setRenderHint(QtGui.QPainter.SmoothPixmapTransform, True)
        painter.setRenderHint(QtGui.QPainter.Antialiasing, True)
        length, height = self._extent
        engine_len = 0.0 if self._contour is None else float(-self._contour[:, 0].min())
        target = self._fit(length + engine_len, height)
        # Pixels per metre, shared by the field image and the engine, so the
        # engine is drawn to the same scale as the flow rather than to taste.
        ppm = target.width() / (length + engine_len)
        field_rect = QtCore.QRect(
            target.x() + int(round(engine_len * ppm)), target.y(),
            target.width() - int(round(engine_len * ppm)), target.height())
        painter.drawImage(field_rect, self._image)
        if self._contour is not None:
            self._draw_engine(painter, target, ppm, height)
        if self._caption:
            painter.setPen(QtGui.QColor(styles.TEXT_MUTED))
            painter.drawText(self.rect().adjusted(10, 6, -10, -6),
                             QtCore.Qt.AlignTop | QtCore.Qt.AlignLeft,
                             self._caption)
        painter.end()

    def _fit(self, width_m: float, height_m: float) -> QtCore.QRect:
        """Largest rectangle with the PHYSICAL aspect ratio that fits the widget.

        The grid is not square -- nx and nr are chosen for accuracy, not for
        shape -- so scaling by pixel counts would stretch the plume.
        """
        available = self.rect().adjusted(8, 8, -8, -8)
        scale = min(available.width() / width_m, available.height() / height_m)
        w, h = int(width_m * scale), int(height_m * scale)
        return QtCore.QRect(available.x() + (available.width() - w) // 2,
                            available.y() + (available.height() - h) // 2, w, h)

    def _draw_engine(self, painter, target, ppm: float, height_m: float) -> None:
        """Fill the engine wall, mirrored, immediately upstream of the field."""
        x, r = self._contour[:, 0], self._contour[:, 1]
        x0 = target.x() + float(-x.min()) * ppm
        mid = target.y() + target.height() / 2.0

        def poly(sign: float) -> QtGui.QPolygonF:
            pts = [QtCore.QPointF(x0 + xi * ppm, mid + sign * ri * ppm)
                   for xi, ri in zip(x, r)]
            wall = max(target.height() * 0.012, 2.0)
            pts += [QtCore.QPointF(x0 + xi * ppm,
                                   mid + sign * (ri * ppm + wall))
                    for xi, ri in zip(x[::-1], r[::-1])]
            return QtGui.QPolygonF(pts)

        painter.setPen(QtGui.QPen(QtGui.QColor("#7e8ba3"), 1.2))
        painter.setBrush(QtGui.QColor("#39424f"))
        for sign in (-1.0, 1.0):
            painter.drawPolygon(poly(sign))


class FlowTab(QtWidgets.QWidget):
    """Run the plume from a solved design and scrub through its time evolution."""

    #: Emitted with a contour path (or "" to go back to the analytic bell) when
    #: the user imports geometry.  The Explorer owns solving, so the tab asks
    #: rather than running the solver itself.
    contourChosen = QtCore.Signal(str)

    def __init__(self) -> None:
        super().__init__()
        self._frames: List[Frame] = []
        self._result: Optional[Result] = None
        self._thread: Optional[QtCore.QThread] = None
        self._worker: Optional[PlumeWorker] = None
        self._ranges: dict = {}
        self._engine_meshes = None
        self._engine_length = 0.0
        self._scene_dirty = True
        self._build()

        self._timer = QtCore.QTimer(self)
        self._timer.setInterval(60)
        self._timer.timeout.connect(self._advance)

    # ------------------------------------------------------------------ build
    def _build(self) -> None:
        outer = QtWidgets.QVBoxLayout(self)
        outer.setContentsMargins(6, 6, 6, 6)
        outer.setSpacing(8)

        geometry, glay = panel("Engine")
        grow = QtWidgets.QHBoxLayout()
        grow.setSpacing(8)
        self.geometry_label = QtWidgets.QLabel("Ignis-M1 (analytic bell)")
        self.geometry_label.setStyleSheet(f"color: {styles.TEXT};")
        grow.addWidget(self.geometry_label, 1)
        self.import_button = QtWidgets.QPushButton("Import engine...")
        self.import_button.clicked.connect(self._import)
        grow.addWidget(self.import_button)
        self.revert_button = QtWidgets.QPushButton("Back to M1")
        self.revert_button.clicked.connect(lambda: self._use_contour(""))
        self.revert_button.setEnabled(False)
        grow.addWidget(self.revert_button)
        glay.addLayout(grow)
        outer.addWidget(geometry)

        self.view = FieldView()
        self.view3d = viewport3d.Viewport3D()
        self.stack = QtWidgets.QStackedWidget()
        self.stack.addWidget(self.view)
        self.stack.addWidget(self.view3d)
        outer.addWidget(self.stack, 1)

        self.colorbar = shell.ColorBar(list(FIELDS), RAMPS)
        self.colorbar.fieldChanged.connect(self._field_changed)
        self.colorbar.rampChanged.connect(lambda _: self._draw())
        outer.addWidget(self.colorbar)

        controls, lay = panel("Flow")
        row = QtWidgets.QHBoxLayout()
        row.setSpacing(8)

        self.run_button = QtWidgets.QPushButton("Run flow")
        self.run_button.clicked.connect(self._toggle)
        self.run_button.setEnabled(False)
        row.addWidget(self.run_button)

        self.play_button = QtWidgets.QPushButton("Play")
        self.play_button.clicked.connect(self._toggle_play)
        self.play_button.setEnabled(False)
        row.addWidget(self.play_button)

        self.view_box = QtWidgets.QComboBox()
        self.view_box.addItems(["2-D slice", "3-D model"])
        self.view_box.currentIndexChanged.connect(self._view_changed)
        row.addWidget(QtWidgets.QLabel("View"))
        row.addWidget(self.view_box)

        self.quality_box = QtWidgets.QComboBox()
        self.quality_box.addItems([f"{g[0]}x{g[1]}  ~{g[6]}" for g in GRIDS])
        row.addWidget(QtWidgets.QLabel("Grid"))
        row.addWidget(self.quality_box)

        self.save_button = QtWidgets.QPushButton("Save video...")
        self.save_button.clicked.connect(self._save_video)
        self.save_button.setEnabled(False)
        row.addWidget(self.save_button)

        row.addStretch(1)
        self.status = QtWidgets.QLabel("")
        self.status.setStyleSheet(f"color: {styles.TEXT_DIM};")
        row.addWidget(self.status)
        lay.addLayout(row)

        self.slider = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        self.slider.setEnabled(False)
        self.slider.valueChanged.connect(lambda _: self._draw())
        lay.addWidget(self.slider)

        note = QtWidgets.QLabel(
            "Inviscid axisymmetric Euler, marched from rest. Shock structure "
            "and wave propagation are solved; turbulent breakup is not "
            "modelled.")
        note.setWordWrap(True)
        note.setStyleSheet(f"color: {styles.TEXT_DIM}; font-size: 11px;")
        lay.addWidget(note)

        outer.addWidget(controls)

    # ------------------------------------------------------------- public API
    def set_result(self, result: Optional[Result]) -> None:
        """Hand in the design the Explorer has solved, or None to disable."""
        self._result = result
        ready = result is not None and bool(result.profile)
        self.run_button.setEnabled(ready and self._thread is None)
        if not ready:
            self.status.setText("Solve a design first.")

    # --------------------------------------------------------------- the march
    def _toggle(self) -> None:
        if self._thread is not None:
            self._worker.stop()
            self.run_button.setText("Stopping...")
            self.run_button.setEnabled(False)
            return
        self._start()

    def _start(self) -> None:
        try:
            exit_state, ambient = self._exit_and_ambient()
        except Exception as exc:                        # noqa: BLE001
            self.status.setText(str(exc))
            return

        nx, nr, length, width, max_steps, every, _ = GRIDS[
            self.quality_box.currentIndex()]

        self._frames = []
        self._ranges = {}
        self.slider.setEnabled(False)
        self.play_button.setEnabled(False)
        self.view.clear()

        self.view.set_contour(self._contour_points(exit_state.radius))
        self._worker = PlumeWorker(exit_state, ambient, nx=nx, nr=nr,
                                   length=length, width=width,
                                   max_steps=max_steps, every=every)
        self._thread = QtCore.QThread(self)
        self._worker.moveToThread(self._thread)
        self._thread.started.connect(self._worker.run)
        self._worker.frame.connect(self._got_frame)
        self._worker.progress.connect(self.status.setText)
        self._worker.failed.connect(lambda m: self.status.setText(f"failed: {m}"))
        self._worker.finished.connect(self._march_finished)
        self._thread.start()
        self.run_button.setText("Stop")

    def _contour_points(self, exit_radius: float) -> Optional[np.ndarray]:
        """The solved contour, shifted so the exit plane sits at x = 0.

        Only the last stretch is drawn: the whole chamber would dwarf the
        plume window and push the interesting part off the left edge.
        """
        profile = self._result.profile if self._result else None
        if not profile or "x" not in profile or "radius" not in profile:
            return None
        x = np.asarray(profile["x"], dtype=float)
        r = np.asarray(profile["radius"], dtype=float)
        x = x - x[-1]
        keep = x >= -6.0 * exit_radius
        if keep.sum() < 2:
            return None
        return np.column_stack([x[keep], r[keep]])

    def _exit_and_ambient(self):
        """Exit state and ambient pressure for the solved design."""
        from ignis_viz.plume import ExitState

        profile = self._result.profile
        if not profile:
            raise RuntimeError("the solved design carries no nozzle profile")

        def last(column: str) -> float:
            return float(profile[column][-1])

        exit_state = ExitState(
            radius=last("radius"), pressure=last("pressure"),
            temperature=last("temperature"), density=last("density"),
            velocity=last("velocity"), mach=last("mach"),
            gamma=last("gamma_s"), molar_mass=last("molar_mass"))

        ambient = self._result.get("performance.ambient_pressure")
        if not np.isfinite(ambient) or ambient <= 0.0:
            raise RuntimeError(
                "this design has no positive ambient pressure, so there is no "
                "plume to solve; a vacuum jet never forms a shock cell")
        if exit_state.pressure / ambient < 0.4:
            raise RuntimeError(
                f"Ignis predicts the nozzle separates at this ambient "
                f"(exit/ambient {exit_state.pressure / ambient:.2f}); marching "
                f"an attached plume would contradict the solver -- raise the "
                f"altitude")
        return exit_state, ambient

    @QtCore.Slot(object)
    def _got_frame(self, frame: Frame) -> None:
        first = not self._frames
        self._frames.append(frame)
        for key, values in frame.data.items():
            # 1st/99.5th percentile rather than min/max: the starting shock
            # carries a handful of cells far hotter than anything that
            # survives, and scaling to them leaves the rest of the plume
            # black for the whole run.  The caption says the true extremes.
            lo, hi = (float(v) for v in np.percentile(values, (1.0, 99.5)))
            true_lo, true_hi = float(values.min()), float(values.max())
            if key in self._ranges:
                p_lo, p_hi, t_lo, t_hi = self._ranges[key]
                self._ranges[key] = (min(p_lo, lo), max(p_hi, hi),
                                     min(t_lo, true_lo), max(t_hi, true_hi))
            else:
                self._ranges[key] = (lo, hi, true_lo, true_hi)
        self.slider.setMaximum(len(self._frames) - 1)
        follow = first or self.slider.value() == len(self._frames) - 2
        if follow:
            self.slider.setValue(len(self._frames) - 1)
        self._draw()

    @QtCore.Slot(object)
    def _march_finished(self, final) -> None:
        if self._thread is not None:
            self._thread.quit()
            self._thread.wait()
        self._thread = None
        self._worker = None
        self.run_button.setText("Run flow")
        self.run_button.setEnabled(self._result is not None)
        if self._frames:
            self.slider.setEnabled(True)
            self.play_button.setEnabled(True)
            self.save_button.setEnabled(True)
            last = self._frames[-1]
            self.status.setText(
                f"{len(self._frames)} frames   "
                f"{last.time * 1e3:.3f} ms of flow")

    # ----------------------------------------------------------------- display
    def _field_changed(self, _: str) -> None:
        # Each field carries the ramp that suits it, but a ramp the user picked
        # by hand should not be overwritten on every field change, so this only
        # moves the selector, silently, when the field itself changes.
        default_ramp = FIELDS[self.colorbar.field()][2]
        if self.colorbar.ramp() != default_ramp:
            self.colorbar.ramp_box.blockSignals(True)
            self.colorbar.ramp_box.setCurrentText(default_ramp)
            self.colorbar.ramp_box.blockSignals(False)
            self.colorbar._ramp = default_ramp
            self.colorbar.strip.set_ramp(default_ramp)
        self._draw()

    def _draw(self) -> None:
        if not self._frames:
            return
        index = min(self.slider.value(), len(self._frames) - 1)
        frame = self._frames[index]
        name = self.colorbar.field()
        key, unit, _ = FIELDS[name]
        values = frame.data[key]
        lo, hi, true_lo, true_hi = self._ranges[key]
        self.colorbar.set_range(lo, hi, unit)
        caption = (f"{name}   full range {true_lo:.4g} to {true_hi:.4g} {unit}\n"
                   f"t = {frame.time * 1e3:.3f} ms   step {frame.step}   "
                   f"frame {index + 1} of {len(self._frames)}")
        if self.stack.currentIndex() == 0:
            self.view.show_field(values, self.colorbar.ramp(), lo, hi,
                                 caption, frame.extent)
        else:
            self._draw3d(frame, values, lo, hi, caption)

    def _view_changed(self, index: int) -> None:
        self.stack.setCurrentIndex(index)
        self._scene_dirty = True
        self._draw()

    def _draw3d(self, frame: Frame, values: np.ndarray, lo: float, hi: float,
                caption: str) -> None:
        """Rebuild the 3-D scene for this frame.

        The engine mesh does not change between frames, so it is built once
        and kept; only the slice texture is remade, which is a colormap lookup
        and nothing more.
        """
        if self._engine_meshes is None and self._result is not None:
            profile = self._result.profile
            x = np.asarray(profile["x"], dtype=float)
            r = np.asarray(profile["radius"], dtype=float)
            temperature = np.asarray(profile["temperature"], dtype=float)
            # A smooth body of revolution does not need 400 stations to look
            # smooth, and every one of them costs 2 * n_theta triangles.
            step = max(1, x.size // 70)
            self._engine_meshes = viewport3d.engine_meshes(
                np.column_stack([x[::step], r[::step]]),
                wall=max(0.006, 0.02 * float(r.min())),
                scalar=temperature[::step], ramp=self.colorbar.ramp(),
                lo=float(temperature.min()), hi=float(temperature.max()),
                n_theta=48)
            self._engine_length = float(x[-1])

        span = self._engine_length or 1.0
        norm = np.clip((values - lo) / (hi - lo), 0.0, 1.0) if hi > lo else values * 0
        rgb = colormaps[self.colorbar.ramp()](norm)[:, :, :3].astype(np.float32)
        cover = frame.data["jet_fraction"]
        texture = np.concatenate([rgb[::-1], rgb], axis=0)
        alpha = np.concatenate([cover[::-1], cover], axis=0).astype(np.float32)

        class _F:
            pass

        f = _F()
        f.x = np.linspace(0.0, frame.extent[0], values.shape[1])
        f.r = np.linspace(0.0, frame.extent[1] * 0.5, values.shape[0])
        meshes = [viewport3d.plume_slice(f, texture, span, alpha=alpha)]
        meshes += self._engine_meshes or []
        # Frame the whole scene, not the engine: the plume runs another ten
        # exit radii past the exit plane, so framing on the engine alone puts
        # most of what was computed outside the window.
        total = span + frame.extent[0]
        self.view3d.set_scene(meshes, target=np.array([total * 0.44, 0.0, 0.0]),
                              span=total * 0.52, caption=caption)

    def _toggle_play(self) -> None:
        if self._timer.isActive():
            self._timer.stop()
            self.play_button.setText("Play")
        else:
            if self.slider.value() >= self.slider.maximum():
                self.slider.setValue(0)
            self._timer.start()
            self.play_button.setText("Pause")

    def _advance(self) -> None:
        if self.slider.value() >= self.slider.maximum():
            self.slider.setValue(0)
        else:
            self.slider.setValue(self.slider.value() + 1)

    # ---------------------------------------------------------------- geometry
    def _import(self) -> None:
        path, _ = QtWidgets.QFileDialog.getOpenFileName(
            self, "Import an engine", "",
            "Engine geometry (*.stl *.csv *.txt *.dat);;"
            "Triangle mesh (*.stl);;Contour table (*.csv *.txt *.dat)")
        if not path:
            return
        if path.lower().endswith(".stl"):
            path = self._contour_from_stl(path)
            if not path:
                return
        self._use_contour(path)

    def _contour_from_stl(self, stl: str) -> str:
        """Reduce a mesh to a contour, asking only what cannot be guessed."""
        axis, ok = QtWidgets.QInputDialog.getItem(
            self, "Axis of revolution",
            "Which model axis is the engine's centreline?",
            ["x", "y", "z"], 0, False)
        if not ok:
            return ""
        units, ok = QtWidgets.QInputDialog.getItem(
            self, "Model units", "The model is drawn in:",
            ["millimetres", "metres", "inches"], 0, False)
        if not ok:
            return ""
        scale = {"millimetres": 1e-3, "metres": 1.0, "inches": 0.0254}[units]

        sys.path.insert(0, os.path.join(
            os.path.dirname(os.path.dirname(os.path.dirname(
                os.path.abspath(__file__)))), "tools"))
        try:
            import contour_from_stl as extractor
        except ImportError as exc:
            self.status.setText(f"cannot load the STL reader: {exc}")
            return ""

        out = os.path.splitext(stl)[0] + "_contour.csv"
        argv = [stl, "-o", out, "--axis", axis, "--scale", str(scale)]
        buf = io.StringIO()
        err = io.StringIO()
        try:
            with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(err):
                code = extractor.main(argv)
        except Exception as exc:                        # noqa: BLE001
            self.status.setText(f"the STL could not be read: {exc}")
            return ""
        report = (buf.getvalue() + err.getvalue()).strip()
        if code != 0:
            QtWidgets.QMessageBox.warning(self, "Import failed", report)
            return ""
        # The axisymmetry warning is the one that matters: a model that is not
        # a body of revolution yields a contour that looks fine and means
        # nothing, so it is put in front of the user rather than logged.
        if "WARNING" in report:
            answer = QtWidgets.QMessageBox.warning(
                self, "This model may not be a body of revolution",
                report + "\n\nUse it anyway?",
                QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No,
                QtWidgets.QMessageBox.No)
            if answer != QtWidgets.QMessageBox.Yes:
                return ""
        self.status.setText(report.splitlines()[-2] if report else "")
        return out

    def _use_contour(self, path: str) -> None:
        self.geometry_label.setText(
            os.path.basename(path) if path else "Ignis-M1 (analytic bell)")
        self.revert_button.setEnabled(bool(path))
        self._frames = []
        self._ranges = {}
        self.slider.setEnabled(False)
        self.play_button.setEnabled(False)
        self.save_button.setEnabled(False)
        self.view.clear()
        self.view.set_contour(None)
        self.contourChosen.emit(path)

    # ------------------------------------------------------------------ export
    def frames(self) -> List[Frame]:
        return list(self._frames)

    def _save_video(self) -> None:
        """Write the captured march out as an mp4, at the size on screen."""
        if not self._frames:
            return
        path, _ = QtWidgets.QFileDialog.getSaveFileName(
            self, "Save flow video", "ignis_flow.mp4", "Video (*.mp4)")
        if not path:
            return
        try:
            import imageio.v2 as imageio
        except ImportError:
            self.status.setText(
                "imageio is not installed -- pip install -r "
                "python/requirements-render.txt")
            return

        was = self.slider.value()
        playing = self._timer.isActive()
        self._timer.stop()
        try:
            with imageio.get_writer(path, fps=24, macro_block_size=8) as writer:
                for index in range(len(self._frames)):
                    self.slider.setValue(index)
                    self._draw()
                    self.view.repaint()
                    image = self.view.grab().toImage().convertToFormat(
                        QtGui.QImage.Format_RGB888)
                    width, height = image.width(), image.height()
                    buf = image.constBits()
                    array = np.frombuffer(buf, np.uint8, count=height * image.bytesPerLine())
                    array = array.reshape(height, image.bytesPerLine())[:, : width * 3]
                    writer.append_data(array.reshape(height, width, 3).copy())
        except Exception as exc:                        # noqa: BLE001
            self.status.setText(f"could not write the video: {exc}")
            return
        finally:
            self.slider.setValue(was)
            if playing:
                self._timer.start()
        self.status.setText(
            f"wrote {len(self._frames)} frames to {path} "
            f"({len(self._frames) / 24.0:.1f} s at 24 fps)")

    def closeEvent(self, event) -> None:      # noqa: N802  (Qt override)
        if self._worker is not None:
            self._worker.stop()
        if self._thread is not None:
            self._thread.quit()
            self._thread.wait(2000)
        super().closeEvent(event)


__all__ = ["FlowTab", "FieldView", "PlumeWorker", "Frame", "FIELDS"]
