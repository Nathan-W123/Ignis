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

from dataclasses import dataclass
from typing import List, Optional

import numpy as np
from matplotlib import colormaps
from PySide6 import QtCore, QtGui, QtWidgets

from . import styles
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
# plume is: a square window spends most of its pixels on ambient air.  The
# times are measured on the machine this was written on and will not hold
# everywhere, which is why they are labelled "~".
GRIDS = [
    (300, 90, 12.0, 4.5, 9000, 60, "2 min"),
    (440, 130, 14.0, 4.5, 13000, 90, "8 min"),
    (640, 190, 16.0, 5.0, 18000, 130, "30 min"),
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

    def __init__(self) -> None:
        super().__init__()
        self._frames: List[Frame] = []
        self._result: Optional[Result] = None
        self._thread: Optional[QtCore.QThread] = None
        self._worker: Optional[PlumeWorker] = None
        self._ranges: dict = {}
        self._build()

        self._timer = QtCore.QTimer(self)
        self._timer.setInterval(60)
        self._timer.timeout.connect(self._advance)

    # ------------------------------------------------------------------ build
    def _build(self) -> None:
        outer = QtWidgets.QVBoxLayout(self)
        outer.setContentsMargins(6, 6, 6, 6)
        outer.setSpacing(8)

        self.view = FieldView()
        outer.addWidget(self.view, 1)

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

        self.field_box = QtWidgets.QComboBox()
        self.field_box.addItems(FIELDS.keys())
        self.field_box.currentTextChanged.connect(self._field_changed)
        row.addWidget(QtWidgets.QLabel("Field"))
        row.addWidget(self.field_box)

        self.ramp_box = QtWidgets.QComboBox()
        self.ramp_box.addItems(RAMPS)
        self.ramp_box.currentTextChanged.connect(lambda _: self._draw())
        row.addWidget(QtWidgets.QLabel("Colour"))
        row.addWidget(self.ramp_box)

        self.quality_box = QtWidgets.QComboBox()
        self.quality_box.addItems([f"{g[0]}x{g[1]}  ~{g[5]}" for g in GRIDS])
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
        default_ramp = FIELDS[self.field_box.currentText()][2]
        if self.ramp_box.currentText() != default_ramp:
            self.ramp_box.blockSignals(True)
            self.ramp_box.setCurrentText(default_ramp)
            self.ramp_box.blockSignals(False)
        self._draw()

    def _draw(self) -> None:
        if not self._frames:
            return
        index = min(self.slider.value(), len(self._frames) - 1)
        frame = self._frames[index]
        name = self.field_box.currentText()
        key, unit, _ = FIELDS[name]
        values = frame.data[key]
        lo, hi, true_lo, true_hi = self._ranges[key]
        caption = (f"{name}   colour {lo:.4g} to {hi:.4g} {unit}   "
                   f"(full range {true_lo:.4g} to {true_hi:.4g})\n"
                   f"t = {frame.time * 1e3:.3f} ms   step {frame.step}   "
                   f"frame {index + 1} of {len(self._frames)}")
        self.view.show_field(values, self.ramp_box.currentText(), lo, hi,
                             caption, frame.extent)

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
