from __future__ import annotations

import dataclasses
import os
import re
from enum import Enum
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np

from .constants import REACTION_TYPE_CHOICES


class RunMode(Enum):
    FAST = "fast"
    BALANCED = "balanced"
    THOROUGH = "thorough"


@dataclasses.dataclass
class Config:
    n_generate: int = 200
    n_keep_scored: int = 50
    n_keep_clustered: int = 40
    n_run_xtb: int = 20
    n_run_crest: int = 5
    rmsd_cutoff: float = 0.5
    clash_cutoff: float = 1.2
    xtb_method: str = "--gfn2"
    crest_method: str = "--gfn2"
    orca_method: str = "B97-3c Opt Freq"
    charge: int = 0
    multiplicity: int = 1
    energy_a: Optional[float] = None
    energy_b: Optional[float] = None
    xtb_ewin_kcal: float = 5.0
    crest_ewin_kcal: float = 3.0
    random_seed: int = 42
    max_placement_attempts: int = 50
    preopt_inputs: bool = True
    cores: int = 4
    maxcore: int = 2000

    @classmethod
    def from_mode(cls, mode: RunMode):
        if mode == RunMode.FAST:
            return cls(n_generate=50, n_keep_scored=20, n_keep_clustered=10, n_run_xtb=5, n_run_crest=1)
        if mode == RunMode.THOROUGH:
            return cls(
                n_generate=1000,
                n_keep_scored=300,
                n_keep_clustered=100,
                n_run_xtb=50,
                n_run_crest=10,
                xtb_ewin_kcal=10.0,
                crest_ewin_kcal=5.0,
            )
        return cls()


class Atom:
    __slots__ = ("symbol", "x", "y", "z")

    def __init__(self, s, x, y, z):
        self.symbol, self.x, self.y, self.z = s, x, y, z

    @property
    def coords(self):
        return np.array([self.x, self.y, self.z], dtype=float)

    @coords.setter
    def coords(self, v):
        self.x, self.y, self.z = float(v[0]), float(v[1]), float(v[2])


class Molecule:
    def __init__(self, atoms: List[Atom], name: str = "mol"):
        self.atoms = atoms
        self.name = name
        self.score = 0.0
        self.energy_eh = 0.0
        self.gibbs_eh = 0.0
        self.imag_freqs = 0
        self.lineage = []

    @classmethod
    def from_xyz(cls, path: str):
        if not path or not os.path.exists(path):
            return cls([])
        with open(path, "r", encoding="utf-8", errors="replace") as f:
            raw = f.readlines()
        # Strip leading/trailing whitespace per line but do NOT discard empty
        # lines so the mandatory (possibly blank) comment line keeps its index.
        stripped = [l.strip() for l in raw]
        # Skip any leading blank lines before the atom count
        idx = 0
        while idx < len(stripped) and not stripped[idx]:
            idx += 1
        try:
            n = int(stripped[idx])
        except (ValueError, IndexError):
            raise ValueError(f"Invalid XYZ format in {path}")
        # The comment is always the very next line (even when empty)
        comment_line = stripped[idx + 1] if idx + 1 < len(stripped) else ""
        atoms = []
        for line in stripped[idx + 2 : idx + 2 + n]:
            p = line.split()
            if len(p) >= 4:
                try:
                    atoms.append(Atom(p[0], float(p[1]), float(p[2]), float(p[3])))
                except ValueError:
                    continue
        return cls(atoms, comment_line or "mol")

    def to_xyz(self, path, comment=""):
        with open(path, "w", encoding="utf-8", newline="\n") as f:
            f.write(f"{len(self.atoms)}\n{comment}\n")
            for a in self.atoms:
                f.write(f"{a.symbol:<4} {a.x:15.6f} {a.y:15.6f} {a.z:15.6f}\n")

    def coords_array(self):
        return np.array([a.coords for a in self.atoms])

    def centroid(self):
        return np.mean(self.coords_array(), axis=0) if len(self.atoms) > 0 else np.array([0, 0, 0])

    def translate(self, v):
        for a in self.atoms:
            a.coords += v

    def rotate(self, R):
        for a in self.atoms:
            a.coords = R @ a.coords

    def merge(self, other):
        self.atoms.extend([Atom(a.symbol, a.x, a.y, a.z) for a in other.atoms])

    def n_atoms(self):
        return len(self.atoms)


def read_multi_xyz(path: str) -> List[Molecule]:
    mols = []
    if not os.path.exists(path):
        return mols
    with open(path, "r", encoding="utf-8", errors="replace") as f:
        while True:
            line = f.readline()
            if not line:
                break
            line = line.strip()
            if not line:
                continue
            try:
                n = int(line)
            except ValueError:
                continue
            comment = f.readline().strip()
            atoms = []
            for _ in range(n):
                p = f.readline().split()
                if len(p) >= 4:
                    atoms.append(Atom(p[0], float(p[1]), float(p[2]), float(p[3])))
            m = Molecule(atoms, name=comment)
            try:
                m.energy_eh = float(comment.split()[0])
            except Exception:
                pass
            mols.append(m)
    return mols

