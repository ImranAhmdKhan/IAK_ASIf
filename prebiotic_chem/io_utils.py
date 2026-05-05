"""
prebiotic_chem.io_utils
=======================
File I/O utilities for the prebiotic computational chemistry workflow.

Provides:
- XYZ file reading / writing
- xTB input preparation and output parsing
- CREST output parsing
- ORCA input template generation and output parsing
- SMILES → 3-D coordinate generation (via RDKit when available)
"""

from __future__ import annotations

import logging
import os
import re
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np

from .constants import HARTREE_TO_KCAL, KCAL_TO_KJ

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Optional RDKit import
# ---------------------------------------------------------------------------
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdmolfiles
    _RDKIT_AVAILABLE = True
except ImportError:
    _RDKIT_AVAILABLE = False


# ---------------------------------------------------------------------------
# XYZ I/O
# ---------------------------------------------------------------------------

def read_xyz(path: str | Path) -> Tuple[List[str], np.ndarray, str]:
    """
    Parse an XYZ file.

    Parameters
    ----------
    path:
        Path to the ``.xyz`` file.

    Returns
    -------
    (symbols, coords, comment)
        - *symbols* : list of element symbols (length N)
        - *coords*  : (N, 3) float array in Å
        - *comment* : second-line comment string
    """
    lines = Path(path).read_text().splitlines()
    n_atoms = int(lines[0].strip())
    comment = lines[1].rstrip() if len(lines) > 1 else ""
    symbols: List[str] = []
    coords: List[List[float]] = []
    for line in lines[2: 2 + n_atoms]:
        parts = line.split()
        symbols.append(parts[0])
        coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
    return symbols, np.array(coords, dtype=float), comment


def write_xyz(
    path: str | Path,
    symbols: Sequence[str],
    coords: np.ndarray,
    comment: str = "",
) -> None:
    """
    Write an XYZ file.

    Parameters
    ----------
    path:
        Output file path.
    symbols:
        Element symbols.
    coords:
        (N, 3) coordinate array in Å.
    comment:
        Comment line (second line of XYZ format).
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [str(len(symbols)), comment]
    for sym, (x, y, z) in zip(symbols, coords):
        lines.append(f"{sym:<4s} {x:14.8f} {y:14.8f} {z:14.8f}")
    path.write_text("\n".join(lines) + "\n")


def read_xyz_trajectory(path: str | Path) -> List[Tuple[List[str], np.ndarray, str]]:
    """
    Parse a multi-structure XYZ trajectory file (ensemble / conformer file).

    Returns
    -------
    List of (symbols, coords, comment) tuples, one per frame.
    """
    text = Path(path).read_text().splitlines()
    frames = []
    i = 0
    while i < len(text):
        line = text[i].strip()
        if not line:
            i += 1
            continue
        try:
            n = int(line)
        except ValueError:
            i += 1
            continue
        comment = text[i + 1].rstrip() if i + 1 < len(text) else ""
        syms: List[str] = []
        coords: List[List[float]] = []
        for j in range(i + 2, i + 2 + n):
            parts = text[j].split()
            syms.append(parts[0])
            coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
        frames.append((syms, np.array(coords, dtype=float), comment))
        i += 2 + n
    return frames


# ---------------------------------------------------------------------------
# SMILES → 3-D coordinates
# ---------------------------------------------------------------------------

def smiles_to_xyz(
    smiles: str,
    output_path: Optional[str | Path] = None,
    n_confs: int = 10,
    comment: str = "",
) -> Tuple[List[str], np.ndarray]:
    """
    Generate 3-D coordinates from a SMILES string using RDKit.

    Parameters
    ----------
    smiles:
        SMILES representation.
    output_path:
        If given, writes the best conformer to this XYZ file.
    n_confs:
        Number of ETKDG conformers to generate; the lowest-energy one is kept.
    comment:
        Comment for the XYZ file.

    Returns
    -------
    (symbols, coords)

    Raises
    ------
    ImportError
        If RDKit is not installed.
    ValueError
        If RDKit cannot parse the SMILES or embed a conformer.
    """
    if not _RDKIT_AVAILABLE:
        raise ImportError(
            "RDKit is required for SMILES → XYZ conversion. "
            "Install it with: conda install -c conda-forge rdkit"
        )
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"RDKit could not parse SMILES: {smiles!r}")
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    ids = AllChem.EmbedMultipleConfs(mol, numConfs=n_confs, params=params)
    if not ids:
        raise ValueError(f"RDKit failed to embed conformer for SMILES: {smiles!r}")
    # Pick conformer with lowest MMFF energy
    props = AllChem.MMFFGetMoleculeProperties(mol)
    best_id, best_e = 0, float("inf")
    for conf_id in ids:
        ff = AllChem.MMFFGetMoleculeForceField(mol, props, confId=conf_id)
        if ff is None:
            continue
        e = ff.CalcEnergy()
        if e < best_e:
            best_e, best_id = e, conf_id

    conf = mol.GetConformer(best_id)
    symbols = [atom.GetSymbol() for atom in mol.GetAtoms()]
    coords = np.array(conf.GetPositions(), dtype=float)

    if output_path is not None:
        write_xyz(output_path, symbols, coords, comment or smiles)
    return symbols, coords


# ---------------------------------------------------------------------------
# xTB helpers
# ---------------------------------------------------------------------------

def run_xtb(
    xyz_path: str | Path,
    flags: List[str],
    work_dir: Optional[str | Path] = None,
    n_cores: int = 4,
) -> Dict[str, Any]:
    """
    Run xTB on a single XYZ file and parse key output values.

    Parameters
    ----------
    xyz_path:
        Path to the input XYZ file.
    flags:
        Command-line flags, e.g. ``["--gfn2", "--alpb", "water", "--opt"]``.
    work_dir:
        Working directory for the xTB run. Defaults to the XYZ file's directory.
    n_cores:
        Number of OpenMP threads.

    Returns
    -------
    dict with keys:
        ``returncode``, ``stdout``, ``energy_hartree`` (float | None),
        ``converged`` (bool), ``opt_xyz`` (path | None).
    """
    xyz_path = Path(xyz_path).resolve()
    if work_dir is None:
        work_dir = xyz_path.parent
    work_dir = Path(work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)

    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(n_cores)
    env.setdefault("MKL_NUM_THREADS", "1")

    cmd = ["xtb", str(xyz_path)] + flags
    logger.debug("xTB command: %s", " ".join(cmd))
    result = subprocess.run(
        cmd, cwd=str(work_dir), capture_output=True, text=True, env=env
    )
    stdout = result.stdout + result.stderr

    energy = _parse_xtb_energy(stdout)
    converged = "GEOMETRY OPTIMIZATION CONVERGED" in stdout or \
                "normal termination of xtb" in stdout.lower()

    opt_xyz: Optional[Path] = None
    candidate = work_dir / "xtbopt.xyz"
    if candidate.exists():
        opt_xyz = candidate

    return {
        "returncode": result.returncode,
        "stdout": stdout,
        "energy_hartree": energy,
        "converged": converged,
        "opt_xyz": opt_xyz,
    }


def _parse_xtb_energy(text: str) -> Optional[float]:
    """Extract the final total energy from xTB stdout."""
    # Matches: "          | TOTAL ENERGY             -X.XXXXXXX Eh |"
    pattern = re.compile(r"TOTAL ENERGY\s+([-\d.]+)\s+Eh")
    matches = pattern.findall(text)
    if matches:
        return float(matches[-1])
    # Fallback: "total energy  :"
    m = re.search(r"total energy\s*:\s*([-\d.]+)", text, re.IGNORECASE)
    if m:
        return float(m.group(1))
    return None


# ---------------------------------------------------------------------------
# CREST helpers
# ---------------------------------------------------------------------------

def run_crest(
    xyz_path: str | Path,
    flags: List[str],
    work_dir: Optional[str | Path] = None,
    n_cores: int = 4,
) -> Dict[str, Any]:
    """
    Run CREST conformer sampling and parse the output ensemble.

    Parameters
    ----------
    xyz_path:
        Input XYZ file.
    flags:
        CREST flags, e.g. ``["--gfn2", "--alpb", "water"]``.
    work_dir:
        Working directory.
    n_cores:
        Number of CPU cores (``-T`` flag added automatically).

    Returns
    -------
    dict with keys:
        ``returncode``, ``stdout``, ``n_conformers`` (int),
        ``ensemble_path`` (Path | None),
        ``conformers`` (list of (symbols, coords, energy_hartree)).
    """
    xyz_path = Path(xyz_path).resolve()
    if work_dir is None:
        work_dir = xyz_path.parent
    work_dir = Path(work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)

    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(n_cores)

    cmd = ["crest", str(xyz_path)] + flags + ["-T", str(n_cores)]
    logger.debug("CREST command: %s", " ".join(cmd))
    result = subprocess.run(
        cmd, cwd=str(work_dir), capture_output=True, text=True, env=env
    )
    stdout = result.stdout + result.stderr

    ensemble_path: Optional[Path] = None
    for candidate_name in ("crest_conformers.xyz", "crest_best.xyz"):
        p = work_dir / candidate_name
        if p.exists():
            ensemble_path = p
            break

    conformers = []
    if ensemble_path and ensemble_path.name == "crest_conformers.xyz":
        frames = read_xyz_trajectory(ensemble_path)
        for syms, coords, comment in frames:
            energy = _parse_energy_from_comment(comment)
            conformers.append((syms, coords, energy))

    return {
        "returncode": result.returncode,
        "stdout": stdout,
        "n_conformers": len(conformers),
        "ensemble_path": ensemble_path,
        "conformers": conformers,
    }


def _parse_energy_from_comment(comment: str) -> Optional[float]:
    """Try to extract an energy value from an XYZ comment line."""
    m = re.search(r"[-]?\d+\.\d+", comment)
    if m:
        return float(m.group())
    return None


# ---------------------------------------------------------------------------
# ORCA input generation
# ---------------------------------------------------------------------------

def write_orca_input(
    path: str | Path,
    symbols: Sequence[str],
    coords: np.ndarray,
    method: str = "B97-3c Opt Freq",
    charge: int = 0,
    multiplicity: int = 1,
    n_cores: int = 4,
    maxcore_mb: int = 2000,
    solvent_keyword: str = "",
    extra_blocks: str = "",
    comment: str = "",
) -> Path:
    """
    Write an ORCA ``.inp`` input file.

    Parameters
    ----------
    path:
        Output file path.
    symbols, coords:
        Molecular geometry.
    method:
        ORCA keyword line (e.g. ``"B97-3c Opt Freq"``).
    charge, multiplicity:
        System charge and multiplicity.
    n_cores:
        PAL cores.
    maxcore_mb:
        Memory per core in MB.
    solvent_keyword:
        If non-empty, appended to the keyword line (e.g. ``"CPCM(Water)"``).
    extra_blocks:
        Additional ``%block ... end`` text inserted before the coordinate block.
    comment:
        ``# comment`` header line.

    Returns
    -------
    Path to the written file.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    kw = method.strip()
    if solvent_keyword:
        kw += f" {solvent_keyword}"

    lines = []
    if comment:
        lines.append(f"# {comment}")
    lines.append(f"! {kw}")
    lines.append(f"%pal nprocs {n_cores} end")
    lines.append(f"%maxcore {maxcore_mb}")
    if extra_blocks:
        lines.append(extra_blocks.strip())
    lines.append(f"* xyz {charge} {multiplicity}")
    for sym, (x, y, z) in zip(symbols, coords):
        lines.append(f"  {sym:<4s} {x:14.8f} {y:14.8f} {z:14.8f}")
    lines.append("*")
    path.write_text("\n".join(lines) + "\n")
    return path


def parse_orca_energy(orca_out: str | Path) -> Optional[float]:
    """
    Parse the final SCF or DFT total energy from an ORCA output file.

    Returns
    -------
    Energy in Hartree, or ``None`` if not found.
    """
    text = Path(orca_out).read_text(errors="replace")
    # "FINAL SINGLE POINT ENERGY" is printed for SP and after geometry steps
    matches = re.findall(r"FINAL SINGLE POINT ENERGY\s+([-\d.]+)", text)
    if matches:
        return float(matches[-1])
    return None


def parse_orca_frequencies(orca_out: str | Path) -> Optional[List[float]]:
    """
    Parse vibrational frequencies (cm⁻¹) from an ORCA frequency output.

    Returns
    -------
    List of frequencies, or ``None`` if not found.
    """
    text = Path(orca_out).read_text(errors="replace")
    block_match = re.search(
        r"VIBRATIONAL FREQUENCIES\s*\n[-]+\n(.*?)\n\n",
        text,
        re.DOTALL,
    )
    if block_match is None:
        return None
    freqs = re.findall(r"\d+:\s+([-\d.]+)\s+cm\*\*-1", block_match.group(1))
    return [float(f) for f in freqs]


def parse_orca_thermochemistry(orca_out: str | Path) -> Dict[str, Optional[float]]:
    """
    Parse ORCA thermochemistry block (ZPE, H, G, S).

    Returns
    -------
    dict with keys: ``zero_point_correction_hartree``, ``enthalpy_hartree``,
    ``gibbs_free_energy_hartree``, ``entropy_cal_mol_K``.
    """
    text = Path(orca_out).read_text(errors="replace")
    result: Dict[str, Optional[float]] = {
        "zero_point_correction_hartree": None,
        "enthalpy_hartree":              None,
        "gibbs_free_energy_hartree":     None,
        "entropy_cal_mol_K":             None,
    }

    def _find(pattern: str) -> Optional[float]:
        m = re.search(pattern, text, re.IGNORECASE)
        return float(m.group(1)) if m else None

    result["zero_point_correction_hartree"] = _find(
        r"Zero point energy\s*\.\.\.\s*([-\d.]+)"
    )
    result["enthalpy_hartree"] = _find(
        r"Total enthalpy\s*\.\.\.\s*([-\d.]+)"
    )
    result["gibbs_free_energy_hartree"] = _find(
        r"Final Gibbs free energy\s*\.\.\.\s*([-\d.]+)"
    )
    # Entropy in cal/mol/K — ORCA prints in Hartree/K; convert
    s_hartree = _find(r"Total entropy correction\s*\.\.\.\s*([-\d.]+)")
    if s_hartree is not None:
        result["entropy_cal_mol_K"] = s_hartree * HARTREE_TO_KCAL * 1000.0
    return result
