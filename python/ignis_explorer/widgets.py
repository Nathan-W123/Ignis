"""Small reusable widgets for the Ignis Engine Explorer."""
from __future__ import annotations

from typing import List, Optional

from PySide6 import QtCore, QtWidgets

from . import styles


def panel(title: str) -> tuple:
    """A titled panel frame, returned with the layout to fill."""
    frame = QtWidgets.QFrame()
    frame.setObjectName("panelFrame")
    outer = QtWidgets.QVBoxLayout(frame)
    outer.setContentsMargins(10, 8, 10, 10)
    outer.setSpacing(6)
    label = QtWidgets.QLabel(title.upper())
    label.setObjectName("panelTitle")
    outer.addWidget(label)
    return frame, outer


class Card(QtWidgets.QFrame):
    """One headline number, with an optional comparison against design B."""

    def __init__(self, label: str, unit: str, fmt: str = "{:.4g}",
                 better: str = "") -> None:
        super().__init__()
        self.setObjectName("card")
        self._fmt = fmt
        self._better = better          # "up", "down" or "" (no preference)
        lay = QtWidgets.QVBoxLayout(self)
        lay.setContentsMargins(10, 7, 10, 8)
        lay.setSpacing(1)
        top = QtWidgets.QLabel(f"{label}  [{unit}]" if unit else label)
        top.setObjectName("cardLabel")
        self._value = QtWidgets.QLabel("—")
        self._value.setObjectName("cardValue")
        self._delta = QtWidgets.QLabel("")
        self._delta.setObjectName("cardDelta")
        lay.addWidget(top)
        lay.addWidget(self._value)
        lay.addWidget(self._delta)

    def set_values(self, a: Optional[float], b: Optional[float]) -> None:
        if a is None or a != a:
            self._value.setText("—")
            self._delta.setText("")
            return
        self._value.setText(self._fmt.format(a))
        if b is None or b != b:
            self._delta.setText("")
            return
        d = a - b
        if abs(b) > 0:
            pct = 100.0 * d / abs(b)
            text = f"Δ {d:+.3g}  ({pct:+.1f} %)"
        else:
            text = f"Δ {d:+.3g}"
        colour = styles.TEXT_DIM
        if self._better and abs(d) > 1e-12:
            good = (d > 0) if self._better == "up" else (d < 0)
            colour = styles.STATUS["good"] if good else styles.STATUS["serious"]
        self._delta.setText(text)
        self._delta.setStyleSheet(f"color: {colour};")


class Field(QtWidgets.QWidget):
    """A labelled numeric control with its unit."""

    changed = QtCore.Signal()

    def __init__(self, label: str, unit: str, lo: float, hi: float, value: float,
                 step: float, decimals: int = 3) -> None:
        super().__init__()
        lay = QtWidgets.QHBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.setSpacing(6)
        name = QtWidgets.QLabel(label)
        name.setObjectName("fieldLabel")
        name.setMinimumWidth(104)
        self.spin = QtWidgets.QDoubleSpinBox()
        self.spin.setRange(lo, hi)
        self.spin.setValue(value)
        self.spin.setSingleStep(step)
        self.spin.setDecimals(decimals)
        self.spin.setMinimumWidth(92)
        self.spin.setKeyboardTracking(False)
        self.spin.valueChanged.connect(lambda _: self.changed.emit())
        unit_label = QtWidgets.QLabel(unit)
        unit_label.setObjectName("unitLabel")
        unit_label.setMinimumWidth(42)
        lay.addWidget(name)
        lay.addWidget(self.spin, 1)
        lay.addWidget(unit_label)

    def value(self) -> float:
        return self.spin.value()

    def set_value(self, v: float) -> None:
        blocked = self.spin.blockSignals(True)
        self.spin.setValue(v)
        self.spin.blockSignals(blocked)

    def set_range(self, lo: float, hi: float) -> None:
        self.spin.setRange(lo, hi)


class Choice(QtWidgets.QWidget):
    """A labelled combo box."""

    changed = QtCore.Signal()

    def __init__(self, label: str, options: List[str], current: str = "") -> None:
        super().__init__()
        lay = QtWidgets.QHBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.setSpacing(6)
        name = QtWidgets.QLabel(label)
        name.setObjectName("fieldLabel")
        name.setMinimumWidth(104)
        self.combo = QtWidgets.QComboBox()
        self.combo.addItems(options)
        if current:
            self.combo.setCurrentText(current)
        self.combo.currentTextChanged.connect(lambda _: self.changed.emit())
        lay.addWidget(name)
        lay.addWidget(self.combo, 1)
        lay.addSpacing(42)

    def value(self) -> str:
        return self.combo.currentText()


class ConstraintRow(QtWidgets.QWidget):
    """One constraint, its value, its bound and whether it is met."""

    def __init__(self) -> None:
        super().__init__()
        lay = QtWidgets.QHBoxLayout(self)
        lay.setContentsMargins(0, 1, 0, 1)
        lay.setSpacing(8)
        self.dot = QtWidgets.QLabel("●")
        self.dot.setFixedWidth(12)
        self.name = QtWidgets.QLabel()
        self.name.setObjectName("fieldLabel")
        self.name.setMinimumWidth(118)
        self.value = QtWidgets.QLabel()
        self.value.setStyleSheet(f"color: {styles.TEXT};")
        self.value.setAlignment(QtCore.Qt.AlignRight | QtCore.Qt.AlignVCenter)
        self.value.setMinimumWidth(78)
        self.bound = QtWidgets.QLabel()
        self.bound.setObjectName("unitLabel")
        self.bound.setMinimumWidth(86)
        lay.addWidget(self.dot)
        lay.addWidget(self.name, 1)
        lay.addWidget(self.value)
        lay.addWidget(self.bound)

    def set_state(self, name: str, value: str, bound: str, status: str) -> None:
        self.name.setText(name)
        self.value.setText(value)
        self.bound.setText(bound)
        self.dot.setStyleSheet(f"color: {styles.STATUS[status]};")
        self.setToolTip(f"{name}: {value}  {bound}")
