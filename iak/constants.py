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
