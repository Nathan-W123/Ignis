"""Ignis Engine Explorer - an interactive front end for the Ignis solver.

The Explorer owns no physics.  Every number on screen came out of the same
`ignis_engine`, `ignis_equilibrium` and `ignis_nozzle` binaries the command
line drives, run on a configuration the UI writes.  When the solver refuses a
design the Explorer shows the refusal; it never fills the gap with an estimate.

Two design slots, A and B, are solved and displayed together so a change can be
read as a difference rather than remembered.
"""
from __future__ import annotations

import os
import sys
from typing import Dict, List, Optional

import matplotlib
matplotlib.use("QtAgg")

from PySide6 import QtCore, QtGui, QtWidgets  # noqa: E402

from . import styles  # noqa: E402
from .charts import (AltitudeChart, AxialChart, CompositionChart,  # noqa: E402
                     ContourChart, ThermalChart)
from .solver import PROPELLANTS, Design, Result, Solver  # noqa: E402
from .widgets import Card, Choice, ConstraintRow, Field, panel  # noqa: E402

ASSUMPTIONS = """\
Ignis-M1 and Ignis-H1 are CONCEPTUAL engines invented for this project. Nothing \
here has been compared with a test stand and nothing here is flight-ready.

• Gas-phase equilibrium only — no condensed carbon, no finite-rate kinetics. \
Frozen and shifting expansion bracket the truth; neither is it.
• Quasi-1D and inviscid. No boundary layer. Separation is predicted by an \
empirical criterion and flagged, but the inviscid solution is not modified, so \
a separated nozzle's reported thrust is optimistic.
• The thermal model is an engineering estimate. Bartz carries ±20–30 % scatter \
against measured rocket heat flux. 1-D wall, no axial conduction, no thermal \
stress, no life analysis.
• η_c* is an assumed input, not a prediction: there is no injector or mixing \
model. Ideal and corrected values are shown separately so the assumption stays \
visible.
• Ideal thermal equation of state (p = ρRT, Z ≡ 1) with fully variable caloric \
properties. No compressibility factor or fugacity model.
• The feed panel sizes pressures; it does not close an engine cycle.

Full detail: docs/limitations.md"""


class SolveWorker(QtCore.QObject):
    """Runs the solver off the UI thread."""

    finished = QtCore.Signal(str, object)

    def __init__(self, solver: Solver) -> None:
        super().__init__()
        self.solver = solver

    @QtCore.Slot(str, object)
    def solve(self, slot: str, design: Design) -> None:
        self.finished.emit(slot, self.solver.run(design))


class DesignPanel(QtWidgets.QWidget):
    """The controls for one design slot."""

    changed = QtCore.Signal()

    def __init__(self, slot: str, design: Design) -> None:
        super().__init__()
        self.slot = slot
        lay = QtWidgets.QVBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.setSpacing(5)

        self.propellant = Choice("Propellant", list(PROPELLANTS), design.propellant)
        mr_lo, mr_hi = PROPELLANTS[design.propellant]["mr_range"]
        self.pressure = Field("Chamber pressure", "MPa", 0.5, 25.0,
                              design.chamber_pressure * 1e-6, 0.25, 3)
        self.mixture = Field("Mixture ratio O/F", "–", mr_lo, mr_hi,
                             design.mixture_ratio, 0.05, 3)
        self.throat = Field("Throat radius", "mm", 10.0, 250.0,
                            design.throat_radius * 1e3, 1.0, 2)
        self.expansion = Field("Expansion ratio", "–", 2.0, 200.0,
                               design.expansion_ratio, 1.0, 2)
        self.altitude = Field("Altitude", "km", 0.0, 80.0,
                              design.altitude * 1e-3, 1.0, 2)
        self.eta = Field("η c*", "–", 0.7, 1.0, design.eta_c_star, 0.01, 3)
        self.composition = Choice("Expansion", ["shifting equilibrium", "frozen"],
                                  "shifting equilibrium")
        self.cooling = QtWidgets.QCheckBox("Regenerative cooling")
        self.cooling.setChecked(design.cooling)
        self.channel = Field("Channel height", "mm", 1.0, 12.0,
                             design.channel_height * 1e3, 0.25, 3)
        self.wall = Field("Wall thickness", "mm", 0.2, 2.5,
                          design.wall_thickness * 1e3, 0.05, 3)

        for w in (self.propellant, self.pressure, self.mixture, self.throat,
                  self.expansion, self.altitude, self.eta, self.composition,
                  self.channel, self.wall):
            lay.addWidget(w)
            w.changed.connect(self.changed)
        lay.addWidget(self.cooling)
        self.cooling.toggled.connect(self.changed)
        self.propellant.changed.connect(self._propellant_changed)
        lay.addStretch(1)

    def _propellant_changed(self) -> None:
        p = PROPELLANTS[self.propellant.value()]
        lo, hi = p["mr_range"]
        self.mixture.set_range(lo, hi)
        self.mixture.set_value(p["mr_default"])

    def design(self) -> Design:
        name = self.propellant.value()
        return Design(
            propellant=name,
            chamber_pressure=self.pressure.value() * 1e6,
            mixture_ratio=self.mixture.value(),
            throat_radius=self.throat.value() * 1e-3,
            expansion_ratio=self.expansion.value(),
            altitude=self.altitude.value() * 1e3,
            composition=("equilibrium" if self.composition.value().startswith("shifting")
                         else "frozen"),
            eta_c_star=self.eta.value(),
            cooling=self.cooling.isChecked(),
            channel_height=self.channel.value() * 1e-3,
            wall_thickness=self.wall.value() * 1e-3,
            coolant_inlet_pressure=15.0e6 if name == "LOX / CH4" else 12.0e6,
        )

    def load(self, d: Design) -> None:
        self.propellant.combo.setCurrentText(d.propellant)
        self.pressure.set_value(d.chamber_pressure * 1e-6)
        self.mixture.set_value(d.mixture_ratio)
        self.throat.set_value(d.throat_radius * 1e3)
        self.expansion.set_value(d.expansion_ratio)
        self.altitude.set_value(d.altitude * 1e-3)
        self.eta.set_value(d.eta_c_star)
        self.channel.set_value(d.channel_height * 1e3)
        self.wall.set_value(d.wall_thickness * 1e-3 * 1e3)
        self.cooling.setChecked(d.cooling)


PRESETS: Dict[str, Design] = {
    "Ignis-M1  sea level": Design(),
    "Ignis-M1  vacuum": Design(expansion_ratio=45.0, altitude=80_000.0),
    "Ignis-H1  upper stage": Design(propellant="LOX / H2", mixture_ratio=5.5,
                                    throat_radius=0.060, expansion_ratio=60.0,
                                    altitude=80_000.0, eta_c_star=0.97,
                                    num_channels=240, channel_height=6.0e-3,
                                    wall_thickness=0.7e-3,
                                    coolant_inlet_pressure=12.0e6),
}


class Explorer(QtWidgets.QMainWindow):
    request = QtCore.Signal(str, object)

    def __init__(self, repo_root: str) -> None:
        super().__init__()
        self.setWindowTitle("Ignis Engine Explorer")
        self.resize(1560, 950)
        self.solver = Solver(repo_root)
        self.results: Dict[str, Optional[Result]] = {"A": None, "B": None}
        self.pending: Dict[str, bool] = {"A": False, "B": False}
        self.compare = False

        styles.apply_matplotlib()
        self.setStyleSheet(styles.STYLESHEET)

        root = QtWidgets.QWidget()
        root.setObjectName("centralRoot")
        self.setCentralWidget(root)
        outer = QtWidgets.QVBoxLayout(root)
        outer.setContentsMargins(0, 0, 0, 0)
        outer.setSpacing(0)
        outer.addWidget(self._ribbon())

        body = QtWidgets.QSplitter(QtCore.Qt.Horizontal)
        body.addWidget(self._left_column())
        body.addWidget(self._charts())
        body.addWidget(self._right_column())
        body.setStretchFactor(0, 0)
        body.setStretchFactor(1, 1)
        body.setStretchFactor(2, 0)
        body.setSizes([300, 790, 460])
        inner = QtWidgets.QWidget()
        inner_lay = QtWidgets.QVBoxLayout(inner)
        inner_lay.setContentsMargins(10, 8, 10, 8)
        inner_lay.addWidget(body)
        outer.addWidget(inner, 1)

        self.status = self.statusBar()
        self._set_status("ready")

        self.thread = QtCore.QThread(self)
        self.worker = SolveWorker(self.solver)
        self.worker.moveToThread(self.thread)
        self.request.connect(self.worker.solve)
        self.worker.finished.connect(self._solved)
        self.thread.start()

        if not self.solver.available():
            QtWidgets.QMessageBox.warning(self, "Ignis binaries not found",
                                          self.solver.missing_message())
        else:
            QtCore.QTimer.singleShot(60, lambda: self.run_slot("A"))

    # --- construction ----------------------------------------------------
    def _ribbon(self) -> QtWidgets.QWidget:
        bar = QtWidgets.QFrame()
        bar.setObjectName("ribbonBar")
        lay = QtWidgets.QHBoxLayout(bar)
        lay.setContentsMargins(14, 8, 14, 8)
        lay.setSpacing(14)

        titles = QtWidgets.QVBoxLayout()
        titles.setSpacing(0)
        t = QtWidgets.QLabel("IGNIS  ·  Engine Explorer")
        t.setObjectName("appTitle")
        sub = QtWidgets.QLabel("thermochemical liquid-rocket propulsion — conceptual designs only")
        sub.setObjectName("appSubtitle")
        titles.addWidget(t)
        titles.addWidget(sub)
        lay.addLayout(titles)
        lay.addSpacing(18)

        self.preset = QtWidgets.QComboBox()
        self.preset.addItems(list(PRESETS))
        self.preset.setMinimumWidth(180)
        self.preset.currentTextChanged.connect(self._load_preset)
        lay.addWidget(QtWidgets.QLabel("Preset"))
        lay.addWidget(self.preset)

        self.compare_box = QtWidgets.QCheckBox("Compare A / B")
        self.compare_box.toggled.connect(self._toggle_compare)
        lay.addWidget(self.compare_box)

        lay.addStretch(1)
        self.run_button = QtWidgets.QPushButton("Run")
        self.run_button.setObjectName("runButton")
        self.run_button.clicked.connect(self.run_all)
        self.run_button.setShortcut(QtGui.QKeySequence("Ctrl+Return"))
        self.run_button.setToolTip("Solve every visible design slot   (Ctrl+Enter)")
        lay.addWidget(self.run_button)
        return bar

    def _left_column(self) -> QtWidgets.QWidget:
        col = QtWidgets.QWidget()
        lay = QtWidgets.QVBoxLayout(col)
        lay.setContentsMargins(0, 0, 6, 0)
        lay.setSpacing(8)

        frame_a, lay_a = panel("Design A")
        self.panel_a = DesignPanel("A", Design())
        self.panel_a.changed.connect(lambda: self._mark_stale("A"))
        lay_a.addWidget(self.panel_a)
        lay.addWidget(frame_a)

        self.frame_b, lay_b = panel("Design B")
        self.panel_b = DesignPanel("B", Design(expansion_ratio=45.0, altitude=80_000.0))
        self.panel_b.changed.connect(lambda: self._mark_stale("B"))
        lay_b.addWidget(self.panel_b)
        self.frame_b.setVisible(False)
        lay.addWidget(self.frame_b)
        lay.addStretch(1)

        scroll = QtWidgets.QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setWidget(col)
        scroll.setMinimumWidth(300)
        scroll.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        return scroll

    def _charts(self) -> QtWidgets.QWidget:
        self.tabs = QtWidgets.QTabWidget()
        self.charts = {
            "Contour": ContourChart(),
            "Axial": AxialChart(),
            "Thermal": ThermalChart(),
            "Composition": CompositionChart(),
            "Altitude": AltitudeChart(),
        }
        for name, chart in self.charts.items():
            holder = QtWidgets.QWidget()
            hl = QtWidgets.QVBoxLayout(holder)
            hl.setContentsMargins(6, 6, 6, 6)
            hl.addWidget(chart)
            self.tabs.addTab(holder, name)
        return self.tabs

    def _right_column(self) -> QtWidgets.QWidget:
        col = QtWidgets.QWidget()
        lay = QtWidgets.QVBoxLayout(col)
        lay.setContentsMargins(6, 0, 0, 0)
        lay.setSpacing(8)

        frame_cards, lay_cards = panel("Performance")
        grid = QtWidgets.QGridLayout()
        grid.setSpacing(6)
        self.cards = {
            "performance.thrust": Card("THRUST", "kN", "{:.1f}", "up"),
            "performance.isp": Card("SPECIFIC IMPULSE", "s", "{:.2f}", "up"),
            "performance.isp_vacuum": Card("VACUUM Isp", "s", "{:.2f}", "up"),
            "chamber.c_star_ideal": Card("C* IDEAL", "m/s", "{:.1f}", "up"),
            "chamber.temperature": Card("FLAME TEMPERATURE", "K", "{:.1f}"),
            "performance.mdot": Card("MASS FLOW", "kg/s", "{:.2f}"),
            "cooling.max_wall_temperature": Card("PEAK WALL TEMPERATURE", "K", "{:.1f}", "down"),
            "cooling.coolant_pressure_drop": Card("COOLANT Δp", "MPa", "{:.3f}", "down"),
        }
        self._card_scale = {"performance.thrust": 1e-3,
                            "cooling.coolant_pressure_drop": 1e-6}
        for i, card in enumerate(self.cards.values()):
            grid.addWidget(card, i // 2, i % 2)
        lay_cards.addLayout(grid)
        lay.addWidget(frame_cards)

        frame_con, lay_con = panel("Constraints and warnings")
        self.constraint_rows = [ConstraintRow() for _ in range(6)]
        for row in self.constraint_rows:
            lay_con.addWidget(row)
        self.warning_box = QtWidgets.QTextEdit()
        self.warning_box.setReadOnly(True)
        self.warning_box.setMaximumHeight(96)
        lay_con.addWidget(self.warning_box)
        lay.addWidget(frame_con)

        frame_note, lay_note = panel("Model assumptions and limitations")
        note = QtWidgets.QTextEdit()
        note.setReadOnly(True)
        note.setPlainText(ASSUMPTIONS)
        note.setMinimumHeight(190)
        lay_note.addWidget(note)
        lay.addWidget(frame_note, 1)

        scroll = QtWidgets.QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setWidget(col)
        scroll.setMinimumWidth(436)
        scroll.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        return scroll

    # --- behaviour -------------------------------------------------------
    def _load_preset(self, name: str) -> None:
        if name in PRESETS:
            self.panel_a.load(PRESETS[name])
            self.run_slot("A")

    def _toggle_compare(self, on: bool) -> None:
        self.compare = on
        self.frame_b.setVisible(on)
        if on and self.results["B"] is None:
            self.run_slot("B")
        else:
            self._refresh()

    def _mark_stale(self, slot: str) -> None:
        self._set_status(f"design {slot} changed — press Run (Ctrl+Enter)")

    def run_all(self) -> None:
        self.run_slot("A")
        if self.compare:
            self.run_slot("B")

    def run_slot(self, slot: str) -> None:
        if not self.solver.available():
            return
        panel_widget = self.panel_a if slot == "A" else self.panel_b
        self.pending[slot] = True
        self.run_button.setEnabled(False)
        self._set_status(f"solving design {slot} …")
        self.request.emit(slot, panel_widget.design())

    @QtCore.Slot(str, object)
    def _solved(self, slot: str, result: Result) -> None:
        self.results[slot] = result
        self.pending[slot] = False
        if not any(self.pending.values()):
            self.run_button.setEnabled(True)
        if not result.ok:
            self._set_status(f"design {slot}: {result.error.splitlines()[0]}", "critical")
        else:
            self._set_status(f"design {slot} solved — Ignis {result.version}")
        self._refresh()

    def _slots(self) -> List[str]:
        return ["A", "B"] if self.compare else ["A"]

    def _refresh(self) -> None:
        slots = self._slots()
        results = [self.results.get(s) for s in slots]
        labels = [f"{s}" for s in slots]
        for chart in self.charts.values():
            chart.show_results(results, labels)

        a = self.results.get("A")
        b = self.results.get("B") if self.compare else None
        for key, card in self.cards.items():
            scale = self._card_scale.get(key, 1.0)
            va = a.get(key) * scale if (a and a.ok) else None
            vb = b.get(key) * scale if (b and b.ok) else None
            card.set_values(va, vb)
        self._refresh_constraints(a)

    def _refresh_constraints(self, r: Optional[Result]) -> None:
        rows: List[tuple] = []
        if r is not None and r.ok:
            regime = r.strings.get("performance.regime", "?")
            sep = r.get("performance.separation_margin")
            if sep >= 1e5:
                rows.append(("Flow separation", "not applicable", "vacuum / criterion off", "good"))
            else:
                state = "good" if sep >= 0.02 else ("warning" if sep >= 0.0 else "critical")
                rows.append(("Separation margin", f"{sep * 100:+.1f} %", "want > 0", state))
            rows.append(("Expansion regime", regime, "", 
                         "good" if regime == "ideally-expanded" else "warning"))
            if r.design is not None and r.design.cooling:
                t = r.get("cooling.max_wall_temperature")
                lim = r.get("cooling.wall_limit_temperature")
                frac = t / lim if lim == lim and lim > 0 else 0.0
                state = "good" if frac < 0.9 else ("warning" if frac <= 1.0 else "critical")
                rows.append(("Peak wall temperature", f"{t:.0f} K",
                             f"limit {lim:.0f} K", state))
                dp = r.get("cooling.coolant_pressure_drop") * 1e-6
                pin = r.design.coolant_inlet_pressure * 1e-6
                state = "good" if dp < 0.35 * pin else ("warning" if dp < 0.6 * pin else "critical")
                rows.append(("Jacket pressure drop", f"{dp:.2f} MPa",
                             f"inlet {pin:.1f} MPa", state))
                boil = r.get("cooling.boiling_detected", 0.0)
                rows.append(("Coolant boiling", "detected" if boil > 0.5 else "none",
                             "", "critical" if boil > 0.5 else "good"))
            else:
                rows.append(("Regenerative cooling", "switched off", "", "warning"))
            res = max(r.get("performance.mass_flow_residual", 0.0),
                      r.get("chamber.element_residual", 0.0))
            rows.append(("Solver residuals", f"{res:.1e}", "want < 1e-10",
                         "good" if res < 1e-10 else "warning"))
        for row, data in zip(self.constraint_rows, rows):
            row.setVisible(True)
            row.set_state(*data)
        for row in self.constraint_rows[len(rows):]:
            row.setVisible(False)

        if r is not None and not r.ok:
            self.warning_box.setPlainText(r.error)
            self.warning_box.setStyleSheet(f"color: {styles.STATUS['critical']};")
        elif r is not None and r.warnings:
            self.warning_box.setPlainText("\n\n".join("⚠  " + w for w in r.warnings))
            self.warning_box.setStyleSheet(f"color: {styles.STATUS['warning']};")
        else:
            self.warning_box.setPlainText("No solver warnings for this design.")
            self.warning_box.setStyleSheet(f"color: {styles.TEXT_DIM};")

    def _set_status(self, text: str, level: str = "") -> None:
        colour = styles.STATUS[level] if level else styles.TEXT_MUTED
        self.status.setStyleSheet(f"color: {colour};")
        self.status.showMessage(text)

    def closeEvent(self, event: QtGui.QCloseEvent) -> None:
        self.thread.quit()
        self.thread.wait(3000)
        super().closeEvent(event)


def launch(repo_root: str = ".") -> int:
    app = QtWidgets.QApplication.instance() or QtWidgets.QApplication(sys.argv)
    app.setApplicationName("Ignis Engine Explorer")
    window = Explorer(os.path.abspath(repo_root))
    window.show()
    return app.exec()
