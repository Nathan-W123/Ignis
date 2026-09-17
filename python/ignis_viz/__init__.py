"""Ignis visualisation and analysis package.

Figures are drawn only from files the C++ solver wrote, so a figure can never
disagree with the numbers in the report.  See `style` for the colour contract.
"""
from . import animation, data, figures, style  # noqa: F401

__all__ = ["animation", "data", "figures", "style"]
__version__ = "1.0.0"
