"""
prebiotic_chem.analysis
=======================
Structural analysis tools for prebiotic chemistry calculations.

Provides:
- RMSD-based conformer clustering (Kabsch algorithm)
- H-bond detection and scoring
- Nucleophilic / electrophilic site classification
- Interaction energy proxy scoring for clusters
- Diversity-based conformer selection
"""

from __future__ import annotations

import math
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from .constants import ATOMIC_NUMBERS


# ---------------------------------------------------------------------------
# Element / atom utilities
# ---------------------------------------------------------------------------

def _z(symbol: str) -> Optional[int]:
    return ATOMIC_NUMBERS.get(symbol.strip().capitalize())


def vdw_radius(symbol: str) -> float:
    """Return approximate van der Waals radius in Å (Bondi 1964 values)."""
    _VDW: Dict[str, float] = {
        "H": 1.20, "C": 1.70, "N": 1.55, "O": 1.52, "F": 1.47,
        "P": 1.80, "S": 1.80, "Cl": 1.75, "Br": 1.85, "I": 1.98,
        "Na": 2.27, "Mg": 1.73, "Ca": 2.31, "K": 2.75, "Fe": 2.05,
    }
    return _VDW.get(symbol.strip().capitalize(), 2.00)


def covalent_radius(symbol: str) -> float:
    """Return approximate single-bond covalent radius in Å (Alvarez 2008)."""
    _COV: Dict[str, float] = {
        "H": 0.31, "C": 0.76, "N": 0.71, "O": 0.66, "F": 0.57,
        "P": 1.07, "S": 1.05, "Cl": 1.02, "Br": 1.20, "I": 1.39,
        "Na": 1.66, "Mg": 1.41, "Ca": 1.76,
    }
    return _COV.get(symbol.strip().capitalize(), 1.50)


def are_bonded(
    sym_i: str, coord_i: np.ndarray,
    sym_j: str, coord_j: np.ndarray,
    tolerance: float = 0.45,
) -> bool:
    """Return ``True`` if the distance between two atoms is within covalent bond range."""
    dist = float(np.linalg.norm(coord_i - coord_j))
    threshold = covalent_radius(sym_i) + covalent_radius(sym_j) + tolerance
    return dist < threshold


# ---------------------------------------------------------------------------
# Kabsch RMSD
# ---------------------------------------------------------------------------

def kabsch_rmsd(
    coords1: np.ndarray,
    coords2: np.ndarray,
) -> float:
    """
    Compute the RMSD between two conformers after optimal superposition
    (Kabsch 1976 algorithm).

    Parameters
    ----------
    coords1, coords2:
        (N, 3) coordinate arrays.

    Returns
    -------
    RMSD in Å.
    """
    c1 = coords1 - coords1.mean(axis=0)
    c2 = coords2 - coords2.mean(axis=0)
    H = c1.T @ c2
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    rotated = c1 @ R.T
    return float(np.sqrt(np.mean(np.sum((rotated - c2) ** 2, axis=1))))


# ---------------------------------------------------------------------------
# RMSD matrix
# ---------------------------------------------------------------------------

def rmsd_matrix(conformers: List[np.ndarray]) -> np.ndarray:
    """
    Compute the pairwise RMSD matrix for a list of conformer coordinate arrays.

    Parameters
    ----------
    conformers:
        List of (N, 3) arrays.

    Returns
    -------
    Symmetric (M, M) RMSD matrix.
    """
    n = len(conformers)
    mat = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            r = kabsch_rmsd(conformers[i], conformers[j])
            mat[i, j] = mat[j, i] = r
    return mat


# ---------------------------------------------------------------------------
# RMSD-based greedy clustering
# ---------------------------------------------------------------------------

def rmsd_cluster(
    conformers: List[np.ndarray],
    cutoff: float = 0.5,
) -> List[int]:
    """
    Greedy RMSD clustering: pick one representative per cluster.

    The first conformer (lowest energy if pre-sorted) is always kept.
    Subsequent conformers are added only if their RMSD to all already-kept
    structures exceeds *cutoff*.

    Parameters
    ----------
    conformers:
        List of (N, 3) arrays, ideally sorted by ascending energy.
    cutoff:
        RMSD cutoff in Å below which two structures are considered the same.

    Returns
    -------
    List of indices of the selected representative structures.
    """
    kept: List[int] = []
    kept_coords: List[np.ndarray] = []
    for i, c in enumerate(conformers):
        if not kept_coords or all(
            kabsch_rmsd(c, c_ref) >= cutoff for c_ref in kept_coords
        ):
            kept.append(i)
            kept_coords.append(c)
    return kept


# ---------------------------------------------------------------------------
# H-bond analysis
# ---------------------------------------------------------------------------

_HBOND_DONOR_HEAVY = frozenset({"N", "O", "S", "F"})
_HBOND_ACCEPTOR    = frozenset({"N", "O", "F", "P", "S"})


def find_hbond_donors(
    symbols: Sequence[str],
    coords: np.ndarray,
) -> List[Tuple[int, int]]:
    """
    Identify hydrogen-bond donor pairs (H index, heavy-atom index).

    A hydrogen is considered a donor H if it is within bonding distance
    of a donor-capable heavy atom (N, O, S, F).

    Returns
    -------
    List of (H_index, heavy_index) tuples.
    """
    heavy = {i: s for i, s in enumerate(symbols)
             if s.capitalize() in _HBOND_DONOR_HEAVY}
    donors = []
    for i, s in enumerate(symbols):
        if s.capitalize() != "H":
            continue
        for j, hs in heavy.items():
            if are_bonded("H", coords[i], hs, coords[j], tolerance=0.3):
                donors.append((i, j))
                break
    return donors


def find_hbond_acceptors(symbols: Sequence[str]) -> List[int]:
    """Return indices of potential H-bond acceptor atoms (N, O, F, P, S)."""
    return [i for i, s in enumerate(symbols) if s.capitalize() in _HBOND_ACCEPTOR]


def detect_hbonds(
    symbols: Sequence[str],
    coords: np.ndarray,
    d_h_cutoff: float = 2.5,
    angle_cutoff_deg: float = 120.0,
) -> List[Dict]:
    """
    Detect hydrogen bonds using distance + angle criteria.

    Criteria (Jeffrey 1997):
      - D…A distance ≤ 3.5 Å  (donor heavy atom to acceptor)
      - H…A distance ≤ 2.5 Å  (default ``d_h_cutoff``)
      - D–H…A angle ≥ 120°    (default ``angle_cutoff_deg``)

    Parameters
    ----------
    symbols, coords:
        Molecular geometry.
    d_h_cutoff:
        Maximum H…A distance in Å.
    angle_cutoff_deg:
        Minimum D–H…A angle in degrees.

    Returns
    -------
    List of dicts: ``donor_H``, ``donor_heavy``, ``acceptor``,
    ``dist_HA_ang``, ``dist_DA_ang``, ``angle_DHA_deg``.
    """
    donors   = find_hbond_donors(symbols, coords)
    acceptors = find_hbond_acceptors(symbols)

    hbonds = []
    for (h_idx, d_idx) in donors:
        for a_idx in acceptors:
            if a_idx == d_idx:
                continue
            dist_HA = float(np.linalg.norm(coords[h_idx] - coords[a_idx]))
            if dist_HA > d_h_cutoff:
                continue
            dist_DA = float(np.linalg.norm(coords[d_idx] - coords[a_idx]))
            # Angle D–H…A
            v1 = coords[d_idx] - coords[h_idx]
            v2 = coords[a_idx] - coords[h_idx]
            cos_a = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2) + 1e-12)
            angle = math.degrees(math.acos(float(np.clip(cos_a, -1.0, 1.0))))
            if angle >= angle_cutoff_deg:
                hbonds.append({
                    "donor_H":      h_idx,
                    "donor_heavy":  d_idx,
                    "acceptor":     a_idx,
                    "dist_HA_ang":  dist_HA,
                    "dist_DA_ang":  dist_DA,
                    "angle_DHA_deg": angle,
                })
    return hbonds


# ---------------------------------------------------------------------------
# Prebiotic interaction scoring
# ---------------------------------------------------------------------------

def score_prebiotic_geometry(
    symbols: Sequence[str],
    coords: np.ndarray,
) -> float:
    """
    Compute a heuristic prebiotic interaction score for a cluster geometry.

    Scoring rules (inspired by ChemRefine's ``score_geometry``):
      +5 per optimal H-bond (H…A distance 1.6–2.3 Å, Gaussian peak at 1.9 Å)
      −25 per clashing H-bond (H…A < 1.3 Å)
      −0.5 × (max atom distance from centroid) — penalises extended structures

    Returns
    -------
    float score (higher is better).
    """
    hbonds = detect_hbonds(symbols, coords)
    score = 0.0
    for hb in hbonds:
        d = hb["dist_HA_ang"]
        if d < 1.3:
            score -= 25.0
        elif 1.6 <= d <= 2.3:
            score += 5.0 * math.exp(-0.5 * ((d - 1.9) / 0.2) ** 2)
    centroid = coords.mean(axis=0)
    max_r = float(np.max(np.linalg.norm(coords - centroid, axis=1)))
    score -= max_r * 0.5
    return score


# ---------------------------------------------------------------------------
# Nucleophilic / electrophilic site identification
# ---------------------------------------------------------------------------

_NUCLEOPHILIC = frozenset({"N", "O", "S", "P", "F", "Cl", "Br", "I"})
_ELECTROPHILIC = frozenset({"C", "B", "Si", "P", "S"})


def nucleophilic_sites(symbols: Sequence[str]) -> List[int]:
    """Return indices of atoms likely to act as nucleophiles."""
    sites = [i for i, s in enumerate(symbols) if s.capitalize() in _NUCLEOPHILIC]
    return sites if sites else [i for i, s in enumerate(symbols) if s.capitalize() != "H"]


def electrophilic_sites(symbols: Sequence[str]) -> List[int]:
    """Return indices of atoms likely to act as electrophiles."""
    sites = [i for i, s in enumerate(symbols) if s.capitalize() in _ELECTROPHILIC]
    return sites if sites else list(range(len(symbols)))


# ---------------------------------------------------------------------------
# Molecular formula
# ---------------------------------------------------------------------------

def molecular_formula(symbols: Sequence[str]) -> str:
    """Return Hill-ordered molecular formula string."""
    counts: Dict[str, int] = {}
    for s in symbols:
        sym = s.strip().capitalize()
        counts[sym] = counts.get(sym, 0) + 1

    def _key(sym: str) -> Tuple:
        if sym == "C":
            return (0, sym)
        if sym == "H":
            return (1, sym)
        return (2, sym)

    parts = []
    for sym in sorted(counts, key=_key):
        n = counts[sym]
        parts.append(f"{sym}{n if n > 1 else ''}")
    return "".join(parts)
