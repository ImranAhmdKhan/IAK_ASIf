"""IAK - Intelligent Automated Kluster-generator Pipeline."""
from __future__ import annotations

from .pipeline import Pipeline
from .models import Config, RunMode, Molecule, Atom
from .constants import EH2KCAL, KCAL2KJ

__all__ = ["IAKApp", "Pipeline", "Config", "RunMode", "Molecule", "Atom", "EH2KCAL", "KCAL2KJ"]


def __getattr__(name: str):
    """Lazy import of GUI-dependent names (tkinter may not be available)."""
    if name == "IAKApp":
        from .app import IAKApp  # noqa: PLC0415
        return IAKApp
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
