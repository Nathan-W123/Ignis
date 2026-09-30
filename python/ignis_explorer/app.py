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
from . import shell  # noqa: E402
from .flow import FlowTab  # noqa: E402
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

    #: Path of an imported wall contour, or "" for the analytic bell.  Owned
    #: here rather than in the Flow tab because it is part of the design the
    #: solver is asked for, not part of how the plume is displayed.
    contour_file: str = ""

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
            contour_file=self.contour_file,
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

        self.setStyleSheet(styles.STYLESHEET + shell.SHELL_STYLE)

        # The window is a simulation application, not a dashboard: a ribbon of
        # verbs on top, the model on the left, the viewport in the middle with
        # its legend beneath it, the read-outs docked right.  Panes are real
        # QDockWidgets so they can be torn off, re-tabbed or closed, which is
        # what anyone who drives one of these every day will try first.
        self.addToolBar(QtCore.Qt.TopToolBarArea, self._ribbon_bar())

        self.viewport = self._charts()
        self.setCentralWidget(self.viewport)

        self.tree = shell.ModelTree()
        self.tree.selected.connect(self._tree_selected)
        self.properties = shell.PropertyGrid()

        browser = shell.dock("Model", self.tree)
        details = shell.dock("Properties", self.properties)
        inputs = shell.dock("Design", self._left_column())
        readouts = shell.dock("Results", self._right_column())

        self.addDockWidget(QtCore.Qt.LeftDockWidgetArea, browser)
        self.addDockWidget(QtCore.Qt.LeftDockWidgetArea, details)
        self.addDockWidget(QtCore.Qt.LeftDockWidgetArea, inputs)
        self.addDockWidget(QtCore.Qt.RightDockWidgetArea, readouts)
        self.splitDockWidget(browser, details, QtCore.Qt.Vertical)
        self.tabifyDockWidget(details, inputs)
        details.raise_()
        self.resizeDocks([browser, details], [270, 330], QtCore.Qt.Vertical)
        self.resizeDocks([browser, readouts], [330, 430], QtCore.Qt.Horizontal)
        self._docks = {"Model": browser, "Properties": details,
                       "Design": inputs, "Results": readouts}
        self._build_tree()

        self.status = self.statusBar()
        self.solver_light = QtWidgets.QLabel()
        self.solver_light.setObjectName("solverLight")
        self.status.addPermanentWidget(self.solver_light)
        self._set_status("ready")

        self.thread = QtCore.QThread(self)
        self.worker = SolveWorker(self.solver)
        self.worker.moveToThread(self.thread)
        self.request.connect(self.worker.solve)
        self.worker.finished.connect(self._solved)
        self.thread.start()

        if not self.solver.available():
            # Without the binaries the design panel cannot solve anything, but
            # the plume can still be marched from a frozen case -- so the
            # window opens on the Flow tab instead of opening on a dead
            # dashboard behind a modal apology.
            self._no_solver()
        else:
            QtCore.QTimer.singleShot(60, lambda: self.run_slot("A"))

    def _no_solver(self) -> None:
        saved = len(self.flow._cases)
        self.run_button.setEnabled(False)
        self.preset.setEnabled(False)
        self.compare_box.setEnabled(False)
        if saved:
            self._show_tab("Flow")
            self._set_status(
                f"solver not built - {saved} saved engines loaded; the plume "
                f"still solves live", "warning")
            self.properties.show_rows(
                "No solver", [
                    ("state", "the Ignis binaries are not built"),
                    ("looked in", self.solver.bin_dir),
                    ("saved engines", str(saved)),
                    ("what still works", "the Flow tab: the plume is marched "
                                         "here, in Python"),
                    ("what does not", "solving a new design, and every other "
                                      "tab"),
                ])
        else:
            QtWidgets.QMessageBox.warning(self, "Ignis binaries not found",
                                          self.solver.missing_message())

    # --- construction ----------------------------------------------------
    def _ribbon_bar(self) -> QtWidgets.QToolBar:
        bar = QtWidgets.QToolBar("Ribbon")
        bar.setObjectName("ribbonToolBar")
        bar.setMovable(False)
        bar.setFloatable(False)

        ribbon = shell.Ribbon()

        case = ribbon.group("Case")
        self.preset = QtWidgets.QComboBox()
        self.preset.addItems(list(PRESETS))
        self.preset.setMinimumWidth(190)
        self.preset.currentTextChanged.connect(self._load_preset)
        case.button("open", "Preset", self.preset.showPopup,
                    tip="Load one of the shipped designs")
        case.stack([self.preset])

        sim = ribbon.group("Simulation")
        self.run_button = sim.button("run", "Solve", self.run_all,
                                     tip="Solve every visible design  (Ctrl+Enter)")
        self.run_button.setShortcut(QtGui.QKeySequence("Ctrl+Return"))
        self.compare_box = QtWidgets.QCheckBox("Compare A / B")
        self.compare_box.toggled.connect(self._toggle_compare)
        sim.stack([self.compare_box])

        geom = ribbon.group("Geometry")
        geom.button("import", "Import", lambda: self._go_to_flow_and("import"),
                    tip="Load an engine wall contour, or a CAD mesh to take one from")
        geom.button("engine", "Contour", lambda: self._show_tab("Contour"),
                    tip="Show the wall contour of the solved design")

        flow = ribbon.group("Flow")
        flow.button("flow", "Run", lambda: self._go_to_flow_and("run"),
                    tip="March the plume from the solved exit state")
        flow.button("play", "Play", lambda: self._go_to_flow_and("play"),
                    tip="Play the captured march")
        flow.button("video", "Export", lambda: self._go_to_flow_and("save"),
                    tip="Write the march out as an mp4")

        view = ribbon.group("View")
        view.button("chart", "Results", lambda: self._toggle_dock("Results"),
                    tip="Show or hide the results pane")
        view.button("grid", "Model", lambda: self._toggle_dock("Model"),
                    tip="Show or hide the model browser")
        view.button("thermal", "Thermal", lambda: self._show_tab("Thermal"),
                    tip="Show the wall and coolant temperatures")

        titles = QtWidgets.QVBoxLayout()
        titles.setSpacing(0)
        title = QtWidgets.QLabel("IGNIS")
        title.setObjectName("appTitle")
        title.setAlignment(QtCore.Qt.AlignRight)
        sub = QtWidgets.QLabel("conceptual designs only")
        sub.setObjectName("appSubtitle")
        sub.setAlignment(QtCore.Qt.AlignRight)
        titles.addWidget(title)
        titles.addWidget(sub)
        holder = QtWidgets.QWidget()
        holder.setLayout(titles)
        ribbon.finish(holder)

        bar.addWidget(ribbon)
        return bar

    # --- the model browser ------------------------------------------------
    def _build_tree(self) -> None:
        self.tree.build([
            ("engine", "Engine", "engine", [
                ("propellants", "Propellants", "flow"),
                ("chamber", "Chamber", "thermal"),
                ("nozzle", "Nozzle", "engine"),
                ("cooling", "Cooling jacket", "grid"),
            ]),
            ("results", "Results", "chart", [
                ("performance", "Performance", "chart"),
                ("geometry", "Geometry", "grid"),
                ("thermal", "Thermal", "thermal"),
                ("plume", "Plume", "flow"),
            ]),
        ])

    def _tree_selected(self, key: str) -> None:
        r = self.results.get("A")
        if r is None or not r.ok:
            self.properties.show_rows(
                key.title(), [("state", "not solved yet")])
            return
        d = r.design

        def num(k: str, scale: float = 1.0, unit: str = "", fmt: str = "{:.4g}") -> str:
            v = r.get(k)
            return "-" if v != v else f"{fmt.format(v * scale)} {unit}".strip()

        pages = {
            "engine": ("Engine", [
                ("name", d.propellant if d else "-"),
                ("solved by", r.version or "-"),
                ("thrust", num("performance.thrust", 1e-3, "kN")),
                ("specific impulse", num("performance.isp", 1.0, "s")),
            ]),
            "propellants": ("Propellants", [
                ("oxidiser", "LOX"),
                ("fuel", "LCH4" if d and d.propellant.endswith("CH4") else "LH2"),
                ("mixture ratio", f"{d.mixture_ratio:.3f}" if d else "-"),
                ("composition", d.composition if d else "-"),
            ]),
            "chamber": ("Chamber", [
                ("pressure", num("chamber.pressure", 1e-6, "MPa")),
                ("temperature", num("chamber.temperature", 1.0, "K")),
                ("c* ideal", num("chamber.c_star_ideal", 1.0, "m/s")),
                ("eta c*", f"{d.eta_c_star:.3f}" if d else "-"),
                ("mass flow", num("performance.mdot", 1.0, "kg/s")),
            ]),
            "nozzle": ("Nozzle", [
                ("contour", "imported" if (d and d.contour_file) else "analytic bell"),
                ("expansion ratio", num("performance.area_ratio")),
                ("exit radius", num("geometry.exit_radius", 1e3, "mm")),
                ("exit Mach", num("performance.exit_mach")),
                ("exit pressure", num("performance.exit_pressure", 1e-3, "kPa")),
                ("regime", r.strings.get("performance.regime", "-")),
            ]),
            "cooling": ("Cooling jacket", [
                ("enabled", "yes" if (d and d.cooling) else "no"),
                ("peak heat flux", num("cooling.max_heat_flux", 1e-6, "MW/m^2")),
                ("peak wall", num("cooling.max_wall_temperature", 1.0, "K")),
                ("coolant rise", num("cooling.coolant_temperature_rise", 1.0, "K")),
                ("pressure drop", num("cooling.coolant_pressure_drop", 1e-6, "MPa")),
            ]),
            "performance": ("Performance", [
                ("thrust", num("performance.thrust", 1e-3, "kN")),
                ("specific impulse", num("performance.isp", 1.0, "s")),
                ("thrust coefficient", num("performance.cf")),
                ("ambient", num("performance.ambient_pressure", 1e-3, "kPa")),
                ("altitude", f"{d.altitude / 1e3:.1f} km" if d else "-"),
            ]),
            "geometry": ("Geometry", [
                ("exit radius", num("geometry.exit_radius", 1e3, "mm")),
                ("total length", num("geometry.total_length", 1e3, "mm")),
                ("L*", num("geometry.l_star", 1e3, "mm")),
                ("chamber volume", num("geometry.chamber_volume", 1e6, "cm^3")),
                ("residence time", num("geometry.residence_time", 1e3, "ms")),
            ]),
            "thermal": ("Thermal", [
                ("peak heat flux", num("cooling.max_heat_flux", 1e-6, "MW/m^2")),
                ("peak wall", num("cooling.max_wall_temperature", 1.0, "K")),
                ("wall limit", num("cooling.wall_limit_temperature", 1.0, "K")),
                ("coolant out", num("cooling.coolant_outlet_temperature", 1.0, "K")),
            ]),
            "plume": ("Plume", [
                ("solver", "axisymmetric Euler, inviscid"),
                ("frames captured", str(len(self.flow.frames()))),
                ("turbulence model", "none - not modelled"),
                ("exit Mach", num("performance.exit_mach")),
            ]),
            "results": ("Results", [
                ("thrust", num("performance.thrust", 1e-3, "kN")),
                ("specific impulse", num("performance.isp", 1.0, "s")),
                ("warnings", str(len(r.warnings))),
            ]),
        }
        title, rows = pages.get(key, (key.title(), [("", "")]))
        self.properties.show_rows(title, rows)

    # --- ribbon actions ---------------------------------------------------
    def _show_tab(self, name: str) -> None:
        for i in range(self.tabs.count()):
            if self.tabs.tabText(i) == name:
                self.tabs.setCurrentIndex(i)
                return

    def _go_to_flow_and(self, what: str) -> None:
        self._show_tab("Flow")
        button = {"run": self.flow.run_button, "play": self.flow.play_button,
                  "save": self.flow.save_button,
                  "import": self.flow.import_button}.get(what)
        if button is not None and button.isEnabled():
            button.click()
        elif button is not None:
            self._set_status(
                {"run": "solve a design first, then Flow > Run",
                 "play": "run the flow first - there are no frames to play",
                 "save": "run the flow first - there is nothing to export",
                 "import": "import is busy"}.get(what, ""), "warning")

    def _toggle_dock(self, name: str) -> None:
        d = self._docks.get(name)
        if d is not None:
            d.setVisible(not d.isVisible())

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

        # The Flow tab is not a chart: it runs its own solver on design A's
        # exit state and owns a worker thread, so it is added by hand and fed
        # separately in _refresh.
        self.flow = FlowTab()
        self.flow.contourChosen.connect(self._contour_chosen)
        self.tabs.addTab(self.flow, "Flow")
        return self.tabs

    def _contour_chosen(self, path: str) -> None:
        """Re-solve design A on an imported contour (or back on the bell)."""
        self.panel_a.contour_file = path
        # Throat radius and expansion ratio are measured from a contour, so the
        # controls that set them stop meaning anything; dimming them is the
        # honest signal that they are no longer in the loop.
        for widget in (self.panel_a.throat, self.panel_a.expansion):
            widget.setEnabled(not path)
        self._set_status(
            f"solving on {os.path.basename(path)} …" if path
            else "solving on the analytic bell …")
        self.run_slot("A")

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
        self.flow.set_result(a if (a is not None and a.ok) else None)
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
