"""Application shell for the Explorer: ribbon, model tree, properties, colour bar.

This is the furniture a simulation application is expected to have, and the
reason for it is not decoration.  A ribbon of grouped verbs, a tree of the
thing being modelled, a property grid for whatever the tree has selected, a
viewport with a legend under it and dockable panes around it -- that layout is
what someone who uses CFD every day already knows how to drive, so the window
stops needing to be explained.

Nothing here computes anything.  Every widget in this module is a way of
showing or collecting a value; the physics lives behind `solver.py` and
`ignis_viz`.

ICONS
-----
Drawn with QPainter rather than shipped as files.  A handful of glyphs at one
size does not justify an asset pipeline, and vector drawing means they stay
crisp on a high-DPI display and recolour with the theme instead of being
baked.
"""
from __future__ import annotations

from typing import Callable, Iterable, List, Optional, Sequence, Tuple

import numpy as np
from matplotlib import colormaps
from PySide6 import QtCore, QtGui, QtWidgets

from . import styles

ICON_SIZE = 22

# Icons carry their own colour, the way a simulation ribbon's do.  A monochrome
# strip makes every verb look the same weight and forces the label to do all
# the work; a green Solve, a red Stop and an amber folder are recognised before
# they are read.
#
# One set serves both themes, which constrains them more than it looks: to
# clear 3:1 against a near-white panel AND against a near-black one, a colour's
# relative luminance has to sit inside roughly 0.14 to 0.26.  The amber and the
# orange were originally outside it and washed out on light; they were solved
# back into the band rather than nudged by eye.  Two separate sets would drift
# apart, and a colour that only works on one theme is wrong half the time.
ICON_COLOURS = {
    "run": "#2e9e4f",        # go
    "play": "#2e9e4f",
    "stop": "#d1434a",       # halt
    "pause": "#b9801c",
    "open": "#b9801c",       # a folder is amber everywhere
    "save": "#3b82c4",
    "import": "#3b82c4",
    "video": "#8b5cd6",
    "fit": "#6b7f99",
    "chart": "#3b82c4",
    "engine": "#7c8b9e",     # metal
    "flow": "#e56426",       # the plume is hot
    "grid": "#6b7f99",
    "thermal": "#d1434a",
}


def _pen(colour: str, width: float = 1.8) -> QtGui.QPen:
    pen = QtGui.QPen(QtGui.QColor(colour), width)
    pen.setCapStyle(QtCore.Qt.RoundCap)
    pen.setJoinStyle(QtCore.Qt.RoundJoin)
    return pen


def icon(name: str, colour: Optional[str] = None) -> QtGui.QIcon:
    """One of the shell's glyphs, drawn at `ICON_SIZE`.

    Without an explicit colour the glyph takes its semantic one from
    `ICON_COLOURS`, falling back to body text for anything unlisted.
    """
    if colour is None:
        colour = ICON_COLOURS.get(name, styles.TEXT)
    size = ICON_SIZE * 2                     # drawn oversize, scaled down
    pix = QtGui.QPixmap(size, size)
    pix.fill(QtCore.Qt.transparent)
    p = QtGui.QPainter(pix)
    p.setRenderHint(QtGui.QPainter.Antialiasing, True)
    p.setPen(_pen(colour, 3.2))
    c = QtGui.QColor(colour)
    m = size * 0.18                          # margin
    s = size - 2 * m

    if name == "run":
        path = QtGui.QPainterPath()
        path.moveTo(m + s * 0.15, m)
        path.lineTo(m + s, m + s * 0.5)
        path.lineTo(m + s * 0.15, m + s)
        path.closeSubpath()
        p.fillPath(path, c)
    elif name == "stop":
        p.fillRect(QtCore.QRectF(m, m, s, s), c)
    elif name == "play":
        path = QtGui.QPainterPath()
        path.moveTo(m + s * 0.2, m)
        path.lineTo(m + s * 0.95, m + s * 0.5)
        path.lineTo(m + s * 0.2, m + s)
        path.closeSubpath()
        p.fillPath(path, c)
    elif name == "pause":
        w = s * 0.28
        p.fillRect(QtCore.QRectF(m + s * 0.08, m, w, s), c)
        p.fillRect(QtCore.QRectF(m + s * 0.64, m, w, s), c)
    elif name == "open":
        p.drawLine(QtCore.QPointF(m, m + s * 0.25), QtCore.QPointF(m + s * 0.42, m + s * 0.25))
        p.drawRect(QtCore.QRectF(m, m + s * 0.25, s, s * 0.62))
        p.drawLine(QtCore.QPointF(m, m + s * 0.25), QtCore.QPointF(m + s * 0.18, m + s * 0.08))
        p.drawLine(QtCore.QPointF(m + s * 0.18, m + s * 0.08),
                   QtCore.QPointF(m + s * 0.5, m + s * 0.08))
    elif name == "save":
        p.drawRect(QtCore.QRectF(m, m, s, s))
        p.drawRect(QtCore.QRectF(m + s * 0.25, m + s * 0.52, s * 0.5, s * 0.48))
        p.drawRect(QtCore.QRectF(m + s * 0.28, m, s * 0.44, s * 0.3))
    elif name == "import":
        p.drawLine(QtCore.QPointF(m + s * 0.5, m), QtCore.QPointF(m + s * 0.5, m + s * 0.62))
        p.drawLine(QtCore.QPointF(m + s * 0.22, m + s * 0.36),
                   QtCore.QPointF(m + s * 0.5, m + s * 0.64))
        p.drawLine(QtCore.QPointF(m + s * 0.78, m + s * 0.36),
                   QtCore.QPointF(m + s * 0.5, m + s * 0.64))
        p.drawLine(QtCore.QPointF(m, m + s), QtCore.QPointF(m + s, m + s))
    elif name == "video":
        p.drawRoundedRect(QtCore.QRectF(m, m + s * 0.2, s * 0.66, s * 0.6),
                          s * 0.08, s * 0.08)
        path = QtGui.QPainterPath()
        path.moveTo(m + s * 0.72, m + s * 0.5)
        path.lineTo(m + s, m + s * 0.24)
        path.lineTo(m + s, m + s * 0.76)
        path.closeSubpath()
        p.fillPath(path, c)
    elif name == "fit":
        for sx, sy in ((0, 0), (1, 0), (0, 1), (1, 1)):
            x = m + sx * s
            y = m + sy * s
            dx = s * 0.26 * (1 if sx == 0 else -1)
            dy = s * 0.26 * (1 if sy == 0 else -1)
            p.drawLine(QtCore.QPointF(x, y), QtCore.QPointF(x + dx, y))
            p.drawLine(QtCore.QPointF(x, y), QtCore.QPointF(x, y + dy))
    elif name == "chart":
        p.drawLine(QtCore.QPointF(m, m), QtCore.QPointF(m, m + s))
        p.drawLine(QtCore.QPointF(m, m + s), QtCore.QPointF(m + s, m + s))
        path = QtGui.QPainterPath()
        path.moveTo(m + s * 0.15, m + s * 0.75)
        path.cubicTo(m + s * 0.45, m + s * 0.72, m + s * 0.5, m + s * 0.2,
                     m + s * 0.92, m + s * 0.12)
        p.drawPath(path)
    elif name == "engine":
        p.drawLine(QtCore.QPointF(m, m + s * 0.3), QtCore.QPointF(m + s * 0.42, m + s * 0.3))
        p.drawLine(QtCore.QPointF(m + s * 0.42, m + s * 0.3),
                   QtCore.QPointF(m + s * 0.58, m + s * 0.5))
        p.drawLine(QtCore.QPointF(m + s * 0.58, m + s * 0.5), QtCore.QPointF(m + s, m + s * 0.12))
        p.drawLine(QtCore.QPointF(m, m + s * 0.7), QtCore.QPointF(m + s * 0.42, m + s * 0.7))
        p.drawLine(QtCore.QPointF(m + s * 0.42, m + s * 0.7),
                   QtCore.QPointF(m + s * 0.58, m + s * 0.5))
        p.drawLine(QtCore.QPointF(m + s * 0.58, m + s * 0.5), QtCore.QPointF(m + s, m + s * 0.88))
    elif name == "flow":
        for k, y in enumerate((0.22, 0.5, 0.78)):
            path = QtGui.QPainterPath()
            path.moveTo(m, m + s * y)
            path.cubicTo(m + s * 0.35, m + s * (y - 0.14 + 0.09 * k),
                         m + s * 0.65, m + s * (y + 0.14 - 0.09 * k),
                         m + s, m + s * y)
            p.drawPath(path)
    elif name == "grid":
        p.drawRect(QtCore.QRectF(m, m, s, s))
        p.drawLine(QtCore.QPointF(m + s / 3, m), QtCore.QPointF(m + s / 3, m + s))
        p.drawLine(QtCore.QPointF(m + 2 * s / 3, m), QtCore.QPointF(m + 2 * s / 3, m + s))
        p.drawLine(QtCore.QPointF(m, m + s / 3), QtCore.QPointF(m + s, m + s / 3))
        p.drawLine(QtCore.QPointF(m, m + 2 * s / 3), QtCore.QPointF(m + s, m + 2 * s / 3))
    elif name == "thermal":
        p.drawEllipse(QtCore.QRectF(m + s * 0.28, m + s * 0.62, s * 0.44, s * 0.38))
        p.drawLine(QtCore.QPointF(m + s * 0.5, m), QtCore.QPointF(m + s * 0.5, m + s * 0.66))
        p.drawLine(QtCore.QPointF(m + s * 0.72, m + s * 0.18), QtCore.QPointF(m + s, m + s * 0.18))
        p.drawLine(QtCore.QPointF(m + s * 0.72, m + s * 0.42), QtCore.QPointF(m + s, m + s * 0.42))
    else:
        p.drawEllipse(QtCore.QRectF(m, m, s, s))
    p.end()
    return QtGui.QIcon(pix.scaled(ICON_SIZE, ICON_SIZE, QtCore.Qt.KeepAspectRatio,
                                  QtCore.Qt.SmoothTransformation))


class RibbonButton(QtWidgets.QToolButton):
    """A ribbon verb: glyph above its label, the way a CFD ribbon reads."""

    def __init__(self, name: str, label: str, slot: Optional[Callable] = None,
                 *, tip: str = "", big: bool = True) -> None:
        super().__init__()
        self.setToolButtonStyle(QtCore.Qt.ToolButtonTextUnderIcon if big
                                else QtCore.Qt.ToolButtonTextBesideIcon)
        self.setIcon(icon(name))
        self.setIconSize(QtCore.QSize(ICON_SIZE, ICON_SIZE))
        self.setText(label)
        self.setToolTip(tip or label)
        self.setAutoRaise(True)
        self.setObjectName("ribbonButton")
        self.setMinimumWidth(60)
        if slot is not None:
            self.clicked.connect(slot)


class RibbonGroup(QtWidgets.QWidget):
    """A titled cluster of ribbon buttons with a rule down its right edge."""

    def __init__(self, title: str) -> None:
        super().__init__()
        self.setObjectName("ribbonGroup")
        outer = QtWidgets.QVBoxLayout(self)
        outer.setContentsMargins(8, 4, 8, 2)
        outer.setSpacing(2)
        self.row = QtWidgets.QHBoxLayout()
        self.row.setSpacing(2)
        self.row.setContentsMargins(0, 0, 0, 0)
        outer.addLayout(self.row, 1)
        caption = QtWidgets.QLabel(title)
        caption.setObjectName("ribbonGroupTitle")
        caption.setAlignment(QtCore.Qt.AlignHCenter)
        outer.addWidget(caption)

    def add(self, widget: QtWidgets.QWidget) -> QtWidgets.QWidget:
        self.row.addWidget(widget)
        return widget

    def button(self, name: str, label: str, slot: Optional[Callable] = None,
               *, tip: str = "") -> RibbonButton:
        return self.add(RibbonButton(name, label, slot, tip=tip))

    def stack(self, widgets: Sequence[QtWidgets.QWidget]) -> None:
        """Put small controls in a column, as a ribbon does with combo boxes."""
        column = QtWidgets.QVBoxLayout()
        column.setSpacing(3)
        column.setContentsMargins(2, 2, 2, 2)
        for w in widgets:
            column.addWidget(w)
        column.addStretch(1)
        holder = QtWidgets.QWidget()
        holder.setLayout(column)
        self.row.addWidget(holder)


class Ribbon(QtWidgets.QFrame):
    """The strip of grouped verbs across the top of the window."""

    def __init__(self) -> None:
        super().__init__()
        self.setObjectName("ribbonBar")
        self._lay = QtWidgets.QHBoxLayout(self)
        self._lay.setContentsMargins(6, 2, 10, 0)
        self._lay.setSpacing(0)

    def group(self, title: str) -> RibbonGroup:
        if self._lay.count():
            rule = QtWidgets.QFrame()
            rule.setObjectName("ribbonRule")
            rule.setFrameShape(QtWidgets.QFrame.VLine)
            self._lay.addWidget(rule)
        g = RibbonGroup(title)
        self._lay.addWidget(g)
        return g

    def finish(self, trailing: Optional[QtWidgets.QWidget] = None) -> None:
        self._lay.addStretch(1)
        if trailing is not None:
            self._lay.addWidget(trailing)


class ModelTree(QtWidgets.QTreeWidget):
    """The model browser: what the case is made of, as a tree.

    Selecting a node is what fills the property grid, so the tree is the
    navigation and the grid is the detail -- the same split every simulation
    front end uses, and the reason neither needs a label explaining it.
    """

    #: node key of whatever is selected
    selected = QtCore.Signal(str)

    def __init__(self) -> None:
        super().__init__()
        self.setObjectName("modelTree")
        self.setHeaderHidden(True)
        self.setIndentation(14)
        self.setAnimated(True)
        self.setRootIsDecorated(True)
        self.currentItemChanged.connect(self._changed)

    def build(self, spec: Iterable[Tuple[str, str, str, Sequence[Tuple[str, str, str]]]]) -> None:
        """`spec` is (key, label, icon, children) with children the same shape."""
        self.clear()
        for key, label, glyph, children in spec:
            parent = QtWidgets.QTreeWidgetItem([label])
            parent.setIcon(0, icon(glyph))
            parent.setData(0, QtCore.Qt.UserRole, key)
            font = parent.font(0)
            font.setBold(True)
            parent.setFont(0, font)
            for ckey, clabel, cglyph in children:
                child = QtWidgets.QTreeWidgetItem([clabel])
                child.setIcon(0, icon(cglyph))
                child.setData(0, QtCore.Qt.UserRole, ckey)
                parent.addChild(child)
            self.addTopLevelItem(parent)
            parent.setExpanded(True)

    def _changed(self, current, _previous) -> None:
        if current is not None:
            key = current.data(0, QtCore.Qt.UserRole)
            if key:
                self.selected.emit(key)


class PropertyGrid(QtWidgets.QTreeWidget):
    """Two-column key/value detail for whatever the tree has selected.

    Values are shown, not edited: the controls that change a design live in
    the design panel, and a grid that looked editable but was not would be
    worse than one that plainly is not.
    """

    def __init__(self) -> None:
        super().__init__()
        self.setObjectName("propertyGrid")
        self.setColumnCount(2)
        self.setHeaderLabels(["Property", "Value"])
        self.setRootIsDecorated(True)
        self.setIndentation(12)
        self.setAlternatingRowColors(True)
        self.header().setSectionResizeMode(0, QtWidgets.QHeaderView.Stretch)
        self.header().setSectionResizeMode(1, QtWidgets.QHeaderView.ResizeToContents)

    def show_rows(self, title: str,
                  rows: Sequence[Tuple[str, str]],
                  groups: Optional[Sequence[Tuple[str, Sequence[Tuple[str, str]]]]] = None) -> None:
        self.clear()
        head = QtWidgets.QTreeWidgetItem([title, ""])
        font = head.font(0)
        font.setBold(True)
        head.setFont(0, font)
        head.setForeground(0, QtGui.QColor(styles.ACCENT))
        self.addTopLevelItem(head)
        for key, value in rows:
            # No explicit foreground: the stylesheet owns text colour, so the
            # grid follows a theme switch instead of keeping the colour it was
            # built with.
            self.addTopLevelItem(QtWidgets.QTreeWidgetItem([key, value]))
        for name, sub in groups or ():
            parent = QtWidgets.QTreeWidgetItem([name, ""])
            pf = parent.font(0)
            pf.setBold(True)
            parent.setFont(0, pf)
            for key, value in sub:
                parent.addChild(QtWidgets.QTreeWidgetItem([key, value]))
            self.addTopLevelItem(parent)
            parent.setExpanded(True)


class ColorBar(QtWidgets.QWidget):
    """The legend strip under a viewport: gradient, its limits, and the map.

    A field plot without one is a picture.  With one it is a measurement, so
    this is not chrome -- it is the difference between the viewport showing
    something and the viewport meaning something.
    """

    rampChanged = QtCore.Signal(str)
    fieldChanged = QtCore.Signal(str)

    def __init__(self, fields: Sequence[str], ramps: Sequence[str]) -> None:
        super().__init__()
        self.setObjectName("colorBar")
        self._ramp = ramps[0]
        self._lo = 0.0
        self._hi = 1.0
        self._unit = ""
        self._gradient: Optional[QtGui.QImage] = None

        lay = QtWidgets.QHBoxLayout(self)
        lay.setContentsMargins(8, 4, 8, 4)
        lay.setSpacing(10)

        self.strip = _GradientStrip(self)
        lay.addWidget(self.strip, 1)

        lay.addWidget(QtWidgets.QLabel("Field"))
        self.field_box = QtWidgets.QComboBox()
        self.field_box.addItems(list(fields))
        self.field_box.currentTextChanged.connect(self.fieldChanged)
        lay.addWidget(self.field_box)

        lay.addWidget(QtWidgets.QLabel("Colour Map"))
        self.ramp_box = QtWidgets.QComboBox()
        self.ramp_box.addItems(list(ramps))
        self.ramp_box.currentTextChanged.connect(self._ramp_changed)
        lay.addWidget(self.ramp_box)

        self.set_range(0.0, 1.0, "")

    def _ramp_changed(self, name: str) -> None:
        self._ramp = name
        self.strip.set_ramp(name)
        self.rampChanged.emit(name)

    def ramp(self) -> str:
        return self._ramp

    def field(self) -> str:
        return self.field_box.currentText()

    def set_range(self, lo: float, hi: float, unit: str) -> None:
        self._lo, self._hi, self._unit = lo, hi, unit
        self.strip.set_labels(lo, hi, unit)
        self.strip.set_ramp(self._ramp)


class _GradientStrip(QtWidgets.QWidget):
    """The gradient itself, with min / mid / max written on it."""

    def __init__(self, parent=None) -> None:
        super().__init__(parent)
        self._ramp = "inferno"
        self._labels = ("", "", "")
        self.setMinimumHeight(26)
        self.setMinimumWidth(240)
        self.setSizePolicy(QtWidgets.QSizePolicy.Expanding,
                           QtWidgets.QSizePolicy.Fixed)

    def set_ramp(self, name: str) -> None:
        self._ramp = name
        self.update()

    def set_labels(self, lo: float, hi: float, unit: str) -> None:
        mid = 0.5 * (lo + hi)
        suffix = f" {unit}" if unit and unit != "-" else ""
        self._labels = (f"{lo:.4g}{suffix}", f"{mid:.4g}", f"{hi:.4g}{suffix}")
        self.update()

    def paintEvent(self, event) -> None:      # noqa: N802  (Qt override)
        p = QtGui.QPainter(self)
        rect = self.rect().adjusted(0, 0, -1, -1)
        try:
            table = colormaps[self._ramp](np.linspace(0.0, 1.0, 256), bytes=True)
        except KeyError:
            table = colormaps["inferno"](np.linspace(0.0, 1.0, 256), bytes=True)
        grad = QtGui.QLinearGradient(rect.left(), 0, rect.right(), 0)
        for i in range(0, 256, 8):
            r, g, b, _ = table[i]
            grad.setColorAt(i / 255.0, QtGui.QColor(int(r), int(g), int(b)))
        p.fillRect(rect, QtGui.QBrush(grad))
        p.setPen(_pen(styles.BORDER, 1.0))
        p.drawRect(rect)

        # The labels sit on the gradient, so each needs a contrasting halo.
        font = p.font()
        font.setPointSizeF(max(7.5, font.pointSizeF() - 0.5))
        p.setFont(font)
        flags = (QtCore.Qt.AlignLeft, QtCore.Qt.AlignHCenter, QtCore.Qt.AlignRight)
        inner = rect.adjusted(6, 0, -6, 0)
        for text, flag in zip(self._labels, flags):
            p.setPen(QtGui.QColor(0, 0, 0, 170))
            p.drawText(inner.adjusted(1, 1, 1, 1), flag | QtCore.Qt.AlignVCenter, text)
            p.setPen(QtGui.QColor("#ffffff"))
            p.drawText(inner, flag | QtCore.Qt.AlignVCenter, text)
        p.end()


def dock(title: str, widget: QtWidgets.QWidget, *,
         areas=QtCore.Qt.AllDockWidgetAreas) -> QtWidgets.QDockWidget:
    """Wrap a widget as a dock pane that can be moved, tabbed or closed."""
    d = QtWidgets.QDockWidget(title)
    d.setObjectName(f"dock_{title.replace(' ', '_').lower()}")
    d.setAllowedAreas(areas)
    d.setFeatures(QtWidgets.QDockWidget.DockWidgetMovable |
                  QtWidgets.QDockWidget.DockWidgetFloatable)
    d.setWidget(widget)
    return d


def shell_style() -> str:
    return f"""
QFrame#ribbonBar {{
    background: {styles.BG_PANEL};
    border-bottom: 1px solid {styles.BORDER};
}}
QFrame#ribbonRule {{
    color: {styles.BORDER};
    max-width: 1px;
    margin: 4px 2px 16px 2px;
}}
QLabel#ribbonGroupTitle {{
    color: {styles.TEXT_DIM};
    font-size: 10px;
    padding-top: 1px;
}}
QToolButton#ribbonButton {{
    color: {styles.TEXT};
    border: 1px solid transparent;
    border-radius: 4px;
    padding: 4px 6px 3px 6px;
    font-size: 11px;
}}
QToolButton#ribbonButton:hover {{
    background: {styles.BG_ELEVATED};
    border-color: {styles.BORDER};
}}
QToolButton#ribbonButton:pressed {{ background: {styles.BG_INPUT}; }}
QToolButton#ribbonButton:disabled {{ color: {styles.TEXT_DIM}; }}

QTreeWidget#modelTree, QTreeWidget#propertyGrid {{
    background: {styles.BG_PANEL};
    border: 1px solid {styles.BORDER};
    border-radius: 4px;
    color: {styles.TEXT};
    outline: 0;
}}
QTreeWidget#modelTree::item, QTreeWidget#propertyGrid::item {{
    padding: 3px 2px;
}}
QTreeWidget#modelTree::item:selected, QTreeWidget#propertyGrid::item:selected {{
    background: {styles.BG_ELEVATED};
    color: {styles.ACCENT_BRIGHT};
}}
QTreeWidget#propertyGrid {{ alternate-background-color: {styles.BG_MAIN}; }}
QHeaderView::section {{
    background: {styles.BG_ELEVATED};
    color: {styles.TEXT_MUTED};
    border: 0;
    border-bottom: 1px solid {styles.BORDER};
    padding: 4px 6px;
    font-size: 11px;
}}

QWidget#colorBar {{
    background: {styles.BG_PANEL};
    border-top: 1px solid {styles.BORDER};
}}

QDockWidget {{
    color: {styles.TEXT_MUTED};
    font-size: 11px;
    titlebar-close-icon: none;
    titlebar-normal-icon: none;
}}
QDockWidget::title {{
    background: {styles.BG_ELEVATED};
    padding: 5px 8px;
    border-bottom: 1px solid {styles.BORDER};
}}
QStatusBar {{
    background: {styles.BG_PANEL};
    color: {styles.TEXT_MUTED};
    border-top: 1px solid {styles.BORDER};
}}
QStatusBar::item {{ border: 0; }}
"""

__all__ = ["ColorBar", "ModelTree", "PropertyGrid", "Ribbon", "RibbonButton",
           "RibbonGroup", "dock", "icon", "shell_style"]
