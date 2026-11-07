"""Low-level helpers with no GRIBS dependencies (A-6 unit 1)."""
from __future__ import annotations

from typing import Any


def _report_value_text(raw: Any, limit: int = 200) -> str:
    """Text for a reported value that never raises.

    Normal values keep their ``repr`` (truncated to ``limit`` characters).
    ``repr`` can fail, e.g. for an int beyond ``sys.get_int_max_str_digits()``
    (ValueError) or an object whose ``__repr__`` raises; such values are
    described without converting them to text.
    """
    try:
        return repr(raw)[:limit]
    except Exception as exc:  # noqa: BLE001 - reporting must not raise
        if isinstance(raw, int) and not isinstance(raw, bool):
            return f"<int with {raw.bit_length()} bits; decimal repr unavailable ({type(exc).__name__})>"
        return f"<{type(raw).__name__} value; repr failed ({type(exc).__name__})>"
