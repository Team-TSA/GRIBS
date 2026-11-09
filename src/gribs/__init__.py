"""GRIBS package (A-6 stage 1, ADR_A6_PACKAGE_SPLIT).

Public API (U-39 (a), unit 10): run_simulation, Result and __version__ (= PROGRAM_VERSION, never the pyproject
version). run_simulation and Result are resolved on first access by importing their owner module (gribs.api,
gribs.simulation.result); `import gribs` itself imports only gribs._version, so it does not import numpy, scipy
or matplotlib. Any other unknown name raises AttributeError.

Code was moved here from the single-file program unit by unit; that program file, at the repository root,
remains the entry point and imports the moved names back. The package never imports the entry.
"""
from __future__ import annotations

import importlib

from gribs._version import PROGRAM_VERSION as __version__

__all__ = ["run_simulation", "Result", "__version__"]

#: public name -> owner module, imported on first access (no from-import, so guard T3 holds)
_LAZY = {"run_simulation": "gribs.api", "Result": "gribs.simulation.result"}


def __getattr__(name):
    owner = _LAZY.get(name)
    if owner is None:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    return getattr(importlib.import_module(owner), name)


def __dir__():
    return sorted(set(globals()) | set(__all__))
