"""Loading Ignis result files.

Every Ignis application writes a CSV table and a JSON document.  The CSV files
carry a ``# units:`` comment line ahead of the header; :func:`read_table`
returns both the data frame and that mapping so axis labels can be built from
the solver's own unit declarations instead of being retyped by hand.
"""
from __future__ import annotations

import json
import os
import re
from dataclasses import dataclass, field

import pandas as pd

_UNIT_RE = re.compile(r"(\S+?)\[([^\]]*)\]")


@dataclass
class Table:
    """A result table plus the units the solver declared for each column."""
    frame: pd.DataFrame
    units: dict[str, str] = field(default_factory=dict)
    path: str = ""

    def __getitem__(self, key):
        return self.frame[key]

    def label(self, column: str, pretty: str | None = None) -> str:
        """Axis label built from the column name and its declared unit."""
        name = pretty if pretty is not None else column.replace("_", " ")
        unit = self.units.get(column, "")
        return f"{name} [{unit}]" if unit and unit != "-" else name

    @property
    def columns(self):
        return self.frame.columns


def read_table(path: str) -> Table:
    units: dict[str, str] = {}
    with open(path) as fh:
        for line in fh:
            if not line.startswith("#"):
                break
            if line.startswith("# units:"):
                units = {m.group(1): m.group(2) for m in _UNIT_RE.finditer(line)}
    frame = pd.read_csv(path, comment="#")
    return Table(frame=frame, units=units, path=path)


def read_json(path: str) -> dict:
    with open(path) as fh:
        return json.load(fh)


def require(path: str) -> str:
    """Return `path`, or raise with a message naming the command that makes it."""
    if os.path.exists(path):
        return path
    raise FileNotFoundError(
        f"{path} is missing. Run scripts/run_all.sh (or the matching ignis_* "
        f"command) to produce it before generating figures.")
