from __future__ import annotations
from typing import Any

EH2KCAL = 627.509
KCAL2KJ = 4.184
ORCA_CMD = "orca"

REACTION_TYPE_CHOICES = [
    "Non-covalent",
    "Covalent",
    "Substitution",
    "Addition",
    "Nucleophilic",
    "Electrophilic",
]

REACTION_TYPE_ALIASES = {
    "nucleophillic": "Nucleophilic",
    "nucleophilic": "Nucleophilic",
    "electrophilic": "Electrophilic",
    "electrophillic": "Electrophilic",
    "electriphillic": "Electrophilic",
}

PERIODIC_SYMBOLS = [
    "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar",
    "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr",
    "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe",
    "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf",
    "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
    "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs",
    "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og",
]

ATOMIC_NUMBERS = {symbol: idx + 1 for idx, symbol in enumerate(PERIODIC_SYMBOLS)}

# ORCA input-file keyword sets — used when building the `!` keyword line.
# Task keywords are stripped from orca_method so only the level-of-theory tokens remain.
ORCA_TASK_KEYWORDS: frozenset = frozenset({
    "opt", "tightopt", "looseopt", "normalopt", "verytightopt",
    "freq", "anfreq", "numfreq", "nofreq",
})

# Substrings that indicate a genuine MPI / launcher failure in an ORCA output file.
# Deliberately narrow: normal ORCA parallel output contains "OpenMPI" header text which
# would false-positive on a plain `"mpi" in content` check.
ORCA_MPI_ERROR_TOKENS: tuple = (
    "mpirun: not found",
    "mpirun was unable",
    "aborting the run",
    "orte_init",
    "mpi_abort",
    "orte has lost",
    "hwloc",
)


def normalize_reaction_type(value: Any) -> str:
    text = str(value or "").strip()
    if not text:
        return "Non-covalent"
    lowered = text.lower()
    if lowered in REACTION_TYPE_ALIASES:
        return REACTION_TYPE_ALIASES[lowered]
    for candidate in REACTION_TYPE_CHOICES:
        if lowered == candidate.lower():
            return candidate
    return "Non-covalent"
