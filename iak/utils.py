from __future__ import annotations

import math
import os
import time
from typing import Optional


def _now():
    return time.strftime("%Y-%m-%d %H:%M:%S")


def _fmt_duration(seconds: float):
    if seconds is None or seconds < 0 or math.isinf(seconds) or math.isnan(seconds):
        return "estimating"
    seconds = int(seconds)
    h, rem = divmod(seconds, 3600)
    m, s = divmod(rem, 60)
    return f"{h:02d}:{m:02d}:{s:02d}"


def _safe_float(value):
    try:
        if value in (None, "", "N/A"):
            return None
        return float(value)
    except Exception:
        return None


def _safe_int(value):
    try:
        return int(value)
    except Exception:
        return None


def _validate_xyz(path: Optional[str], label: str, required: bool):
    record = {
        "label": label,
        "path": os.path.abspath(path) if path else "",
        "required": bool(required),
        "provided": bool(path),
        "valid": False,
        "result": "INCORRECT",
        "atoms_declared": 0,
        "atoms_parsed": 0,
        "message": "",
    }
    if not path:
        record["valid"] = not required
        record["result"] = "CORRECT" if record["valid"] else "INCORRECT"
        record["message"] = "Optional input not provided." if not required else "Required XYZ input was not provided."
        return record
    if not os.path.isfile(path):
        record["message"] = "File does not exist."
        return record
    try:
        with open(path, "r", encoding="utf-8", errors="replace") as handle:
            lines = [line.rstrip("\n") for line in handle]
        n_atoms = int(lines[0].strip())
        parsed = 0
        malformed = []
        for idx, line in enumerate(lines[2 : 2 + n_atoms], start=3):
            parts = line.split()
            if len(parts) < 4:
                malformed.append(idx)
                continue
            try:
                float(parts[1])
                float(parts[2])
                float(parts[3])
                parsed += 1
            except Exception:
                malformed.append(idx)
        record["atoms_declared"] = n_atoms
        record["atoms_parsed"] = parsed
        record["valid"] = n_atoms > 0 and parsed == n_atoms and not malformed
        record["result"] = "CORRECT" if record["valid"] else "INCORRECT"
        if record["valid"]:
            record["message"] = f"Valid XYZ with {n_atoms} atoms."
        elif malformed:
            record["message"] = f"Malformed coordinate lines: {', '.join(map(str, malformed[:8]))}."
        else:
            record["message"] = f"Declared {n_atoms} atoms but parsed {parsed}."
    except Exception as exc:
        record["message"] = f"Could not parse XYZ: {exc}"
    return record

