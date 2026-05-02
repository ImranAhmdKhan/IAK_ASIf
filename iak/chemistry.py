from __future__ import annotations

import math
import random
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

from .constants import ATOMIC_NUMBERS, normalize_reaction_type
from .models import Atom, Molecule


def _canonical_symbol(symbol: str) -> str:
    text = str(symbol or "").strip()
    if not text:
        return text
    return text[0].upper() + text[1:].lower()


def atom_counts_for_molecule(mol: Molecule) -> Dict[str, int]:
    counts: Dict[str, int] = {}
    if not mol:
        return counts
    for atom in mol.atoms:
        symbol = _canonical_symbol(atom.symbol)
        counts[symbol] = counts.get(symbol, 0) + 1
    return counts


def _merge_weighted_atom_counts(target: Dict[str, int], source: Dict[str, int], scale: int):
    if scale <= 0:
        return
    for symbol, count in source.items():
        target[symbol] = target.get(symbol, 0) + count * scale


def expected_cluster_atom_counts(anchor_mol: Molecule, guest_mols: List[Molecule], n_anchor: int, n_guests_list: List[int]) -> Dict[str, int]:
    counts: Dict[str, int] = {}
    if n_anchor > 0 and anchor_mol is not None:
        _merge_weighted_atom_counts(counts, atom_counts_for_molecule(anchor_mol), n_anchor)
    for idx, guest in enumerate(guest_mols or []):
        copies = int(n_guests_list[idx]) if idx < len(n_guests_list) else 0
        if copies > 0:
            _merge_weighted_atom_counts(counts, atom_counts_for_molecule(guest), copies)
    return counts


def _formula_sort_key(symbol: str):
    if symbol == "C":
        return (0, symbol)
    if symbol == "H":
        return (1, symbol)
    return (2, symbol)


def atom_counts_to_formula(counts: Dict[str, int]) -> str:
    if not counts:
        return "(empty)"
    parts = []
    for symbol in sorted(counts.keys(), key=_formula_sort_key):
        count = counts[symbol]
        parts.append(f"{symbol}{count if count != 1 else ''}")
    return " ".join(parts)


def compare_atom_counts(reference: Dict[str, int], observed: Dict[str, int]) -> Tuple[bool, str]:
    if reference == observed:
        return True, "Atom counts match expected stoichiometry."
    return (
        False,
        f"Atom composition mismatch. Expected [{atom_counts_to_formula(reference)}] but got [{atom_counts_to_formula(observed)}].",
    )


def total_electrons_for_system(anchor_mol: Molecule, guest_mols: List[Molecule], n_anchor: int, n_guests_list: List[int], total_charge: int) -> Tuple[Optional[int], List[str]]:
    unknown_symbols = set()
    nuclear_charge = 0

    def _accumulate(mol: Molecule, copies: int):
        nonlocal nuclear_charge
        if not mol or copies <= 0:
            return
        for atom in mol.atoms:
            symbol = _canonical_symbol(atom.symbol)
            z = ATOMIC_NUMBERS.get(symbol)
            if z is None:
                unknown_symbols.add(symbol)
                continue
            nuclear_charge += z * copies

    _accumulate(anchor_mol, n_anchor)
    for idx, guest in enumerate(guest_mols or []):
        copies = int(n_guests_list[idx]) if idx < len(n_guests_list) else 0
        _accumulate(guest, copies)

    if unknown_symbols:
        return None, sorted(unknown_symbols)
    return nuclear_charge - int(total_charge), []


def suggest_multiplicity(total_electrons: int, current_mult: int) -> int:
    required_parity = 1 - (total_electrons % 2)
    upper = max(3, min(total_electrons + 1, max(current_mult + 8, 25)))
    candidates = [m for m in range(1, upper + 1) if m % 2 == required_parity]
    if not candidates:
        return max(1, current_mult)
    return min(candidates, key=lambda m: abs(m - current_mult))


def validate_charge_multiplicity(
    anchor_mol: Molecule,
    guest_mols: List[Molecule],
    n_anchor: int,
    n_guests_list: List[int],
    total_charge: int,
    multiplicity: int,
) -> Dict[str, Any]:
    result = {
        "valid": True,
        "total_electrons": None,
        "unknown_symbols": [],
        "suggested_multiplicity": multiplicity,
        "message": "Charge/multiplicity accepted.",
    }
    if multiplicity < 1:
        result["valid"] = False
        result["message"] = "Multiplicity must be >= 1."
        return result

    total_electrons, unknown_symbols = total_electrons_for_system(anchor_mol, guest_mols, n_anchor, n_guests_list, total_charge)
    result["total_electrons"] = total_electrons
    result["unknown_symbols"] = unknown_symbols

    if unknown_symbols:
        result["message"] = (
            "Electron-based validation skipped because atomic numbers are unavailable for: "
            + ", ".join(unknown_symbols)
        )
        return result

    if total_electrons is None or total_electrons <= 0:
        result["valid"] = False
        result["message"] = (
            f"Invalid total electron count ({total_electrons}) from charge {total_charge}. "
            "Please verify stoichiometry and charge."
        )
        return result

    parity_ok = (total_electrons % 2) != (multiplicity % 2)
    upper_limit = total_electrons + 1
    if multiplicity > upper_limit:
        parity_ok = False

    suggested = suggest_multiplicity(total_electrons, multiplicity)
    result["suggested_multiplicity"] = suggested
    if parity_ok:
        result["message"] = f"Total electrons: {total_electrons}. Charge/multiplicity parity is consistent."
        return result

    result["valid"] = False
    result["message"] = (
        f"Charge/multiplicity mismatch for {total_electrons} electrons (charge={total_charge}, multiplicity={multiplicity}). "
        f"Suggested multiplicity: {suggested}."
    )
    return result


def validate_reactant_product_atom_balance(reactant_mol: Molecule, product_mol: Molecule, reaction_type: str) -> Tuple[bool, str]:
    r_count = reactant_mol.n_atoms()
    p_count = product_mol.n_atoms()
    if r_count != p_count:
        return (
            False,
            f"{reaction_type} path requires equal atom count in reactant/product. Found {r_count} vs {p_count}.",
        )

    react_counts = atom_counts_for_molecule(reactant_mol)
    prod_counts = atom_counts_for_molecule(product_mol)
    if react_counts != prod_counts:
        return (
            False,
            f"{reaction_type} path requires atom-type conservation. Reactant [{atom_counts_to_formula(react_counts)}], "
            f"Product [{atom_counts_to_formula(prod_counts)}].",
        )
    return True, "Reactant/product atom balance validated."


def find_hbond_acceptors(mol):
    _acc_set = {"O", "N", "F"}
    return [i for i, a in enumerate(mol.atoms) if _canonical_symbol(a.symbol) in _acc_set]


def find_hbond_donors(mol):
    donors = []
    _don_heavy = {"N", "O", "S", "F"}
    heavy = {i: a for i, a in enumerate(mol.atoms) if _canonical_symbol(a.symbol) in _don_heavy}
    for i, a in enumerate(mol.atoms):
        if _canonical_symbol(a.symbol) == "H":
            for hj, ha in heavy.items():
                if np.linalg.norm(a.coords - ha.coords) < 1.2:
                    donors.append((i, hj))
                    break
    return donors


def find_nucleophilic_sites(mol):
    preferred = {"N", "O", "S", "P", "F", "Cl", "Br", "I"}
    sites = [i for i, atom in enumerate(mol.atoms) if _canonical_symbol(atom.symbol) in preferred]
    if sites:
        return sites
    return [i for i, atom in enumerate(mol.atoms) if _canonical_symbol(atom.symbol) != "H"]


def find_electrophilic_sites(mol):
    preferred = {"C", "B", "Si", "P", "S"}
    sites = [i for i, atom in enumerate(mol.atoms) if _canonical_symbol(atom.symbol) in preferred]
    if sites:
        return sites
    return [i for i in range(mol.n_atoms())]


def kabsch_rmsd(c1, c2):
    c1_c, c2_c = c1 - np.mean(c1, axis=0), c2 - np.mean(c2, axis=0)
    H = c1_c.T @ c2_c
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    return float(np.sqrt(np.mean(np.sum((c1_c - c2_c @ R.T) ** 2, axis=1))))


def score_geometry(mol):
    score = 0.0
    c = mol.coords_array()
    acc = find_hbond_acceptors(mol)
    don = [d[0] for d in find_hbond_donors(mol)]
    for d_idx in don:
        for a_idx in acc:
            dist = np.linalg.norm(c[d_idx] - c[a_idx])
            if 1.6 < dist < 2.3:
                score += 5.0 * np.exp(-0.5 * ((dist - 1.9) / 0.2) ** 2)
            elif dist < 1.3:
                score -= 25.0
    score -= np.max(np.linalg.norm(c - mol.centroid(), axis=1)) * 0.5
    mol.score = score
    return score


def _axis_angle_matrix(axis, angle):
    axis = axis / (np.linalg.norm(axis) + 1e-12)
    c, s = math.cos(angle), math.sin(angle)
    x, y, z = axis
    return np.array(
        [
            [c + x * x * (1 - c), x * y * (1 - c) - z * s, x * z * (1 - c) + y * s],
            [y * x * (1 - c) + z * s, c + y * y * (1 - c), y * z * (1 - c) - x * s],
            [z * x * (1 - c) - y * s, z * y * (1 - c) + x * s, c + z * z * (1 - c)],
        ]
    )


def _align_vectors(v_from, v_to):
    v_from = v_from / (np.linalg.norm(v_from) + 1e-12)
    v_to = v_to / (np.linalg.norm(v_to) + 1e-12)
    cross, dot = np.cross(v_from, v_to), np.dot(v_from, v_to)
    if np.linalg.norm(cross) < 1e-6:
        perp = np.array([1, 0, 0]) if abs(v_from[0]) < 0.9 else np.array([0, 1, 0])
        return _axis_angle_matrix(np.cross(v_from, perp), math.pi) if dot < 0 else np.eye(3)
    return _axis_angle_matrix(cross, math.acos(np.clip(dot, -1.0, 1.0)))


def _reaction_profile(reaction_type: str) -> Dict[str, Any]:
    profile_map = {
        "Non-covalent": {"target_role": "acceptor", "guest_role": "donor_h", "distance": (1.8, 2.2)},
        "Covalent": {"target_role": "electrophile", "guest_role": "nucleophile", "distance": (1.55, 2.05)},
        "Substitution": {"target_role": "electrophile", "guest_role": "nucleophile", "distance": (1.65, 2.15)},
        "Addition": {"target_role": "electrophile", "guest_role": "nucleophile", "distance": (1.60, 2.10)},
        "Nucleophilic": {"target_role": "electrophile", "guest_role": "nucleophile", "distance": (1.55, 2.00)},
        "Electrophilic": {"target_role": "nucleophile", "guest_role": "electrophile", "distance": (1.55, 2.00)},
    }
    return profile_map.get(normalize_reaction_type(reaction_type), profile_map["Non-covalent"])


def _pick_sites_by_role(mol: Molecule, role: str) -> List[int]:
    if role == "acceptor":
        return find_hbond_acceptors(mol)
    if role == "donor_h":
        return [idx for idx, _ in find_hbond_donors(mol)]
    if role == "nucleophile":
        return find_nucleophilic_sites(mol)
    if role == "electrophile":
        return find_electrophilic_sites(mol)
    return list(range(mol.n_atoms()))


def generate_cluster(anchor_in, guests_in, n_guests_list, rng, max_att=50, n_anchor=1, reaction_type="Non-covalent"):
    def _clone(mol: Molecule) -> Molecule:
        return Molecule([Atom(a.symbol, a.x, a.y, a.z) for a in mol.atoms], name=mol.name)

    cluster = Molecule([])
    if anchor_in is not None and anchor_in.n_atoms() > 0 and n_anchor > 0:
        cluster = _clone(anchor_in)
        for extra_idx in range(1, n_anchor):
            placed_anchor = False
            for _ in range(max(5, max_att // 2)):
                extra = _clone(anchor_in)
                extra.translate(-extra.centroid())
                direction = np.array([rng.uniform(-1.0, 1.0) for _ in range(3)], dtype=float)
                if np.linalg.norm(direction) < 1e-6:
                    direction = np.array([0.0, 0.0, 1.0], dtype=float)
                direction = direction / np.linalg.norm(direction)
                distance = rng.uniform(3.8, 6.4 + 0.25 * extra_idx)
                extra.translate(cluster.centroid() + direction * distance)
                clash = False
                existing = cluster.coords_array()
                candidate = extra.coords_array()
                for p in existing:
                    if np.min(np.linalg.norm(candidate - p, axis=1)) < 1.25:
                        clash = True
                        break
                if not clash:
                    cluster.merge(extra)
                    placed_anchor = True
                    break
            if not placed_anchor:
                return None

    placed_guests = []
    profile = _reaction_profile(reaction_type)

    if sum(n_guests_list) == 0:
        return cluster

    def _pick_target_molecule():
        available = []
        if cluster.n_atoms() > 0:
            available.append(cluster)
        available.extend(placed_guests)
        if not available:
            return None
        if len(available) > 1 and rng.random() > 0.4:
            return rng.choice(available)
        return available[0]

    for g_idx, guest_in in enumerate(guests_in):
        n_copies = int(n_guests_list[g_idx]) if g_idx < len(n_guests_list) else 0
        for _ in range(n_copies):
            success = False
            if guest_in is None or guest_in.n_atoms() == 0:
                # Nothing physical to place; count as placed so we don't abort
                success = True
                break
            for _ in range(max_att):
                g = _clone(guest_in)
                g.translate(-g.centroid())
                target_mol = _pick_target_molecule()
                if target_mol is None:
                    g.translate(np.array([rng.uniform(-1.2, 1.2) for _ in range(3)]))
                    placed_guests.append(g)
                    success = True
                    break

                target_sites = _pick_sites_by_role(target_mol, profile["target_role"])
                if not target_sites:
                    target_sites = list(range(target_mol.n_atoms()))
                if not target_sites:
                    break
                target_idx = rng.choice(target_sites)

                donor_pair = None
                if profile["guest_role"] == "donor_h":
                    donors = find_hbond_donors(g)
                    if donors:
                        donor_pair = rng.choice(donors)
                        guest_site = donor_pair[0]
                    else:
                        fallback = _pick_sites_by_role(g, "nucleophile")
                        guest_site = rng.choice(fallback if fallback else list(range(g.n_atoms())))
                else:
                    guest_sites = _pick_sites_by_role(g, profile["guest_role"])
                    if not guest_sites:
                        guest_sites = list(range(g.n_atoms()))
                    guest_site = rng.choice(guest_sites)

                out_vec = target_mol.atoms[target_idx].coords - target_mol.centroid()
                if np.linalg.norm(out_vec) < 0.1:
                    out_vec = np.array([0.0, 0.0, 1.0])

                if donor_pair is not None:
                    d_vec = g.atoms[donor_pair[0]].coords - g.atoms[donor_pair[1]].coords
                else:
                    d_vec = g.atoms[guest_site].coords - g.centroid()
                    if np.linalg.norm(d_vec) < 0.1:
                        d_vec = np.array([1.0, 0.0, 0.0])

                g.rotate(_align_vectors(d_vec, -out_vec))
                g.rotate(_axis_angle_matrix(out_vec, rng.uniform(0, 2 * math.pi)))

                dmin, dmax = profile["distance"]
                target_pos = target_mol.atoms[target_idx].coords + (out_vec / np.linalg.norm(out_vec)) * rng.uniform(dmin, dmax)
                g.translate(target_pos - g.atoms[guest_site].coords)

                clash = False
                for existing in [cluster] + placed_guests:
                    if existing.n_atoms() == 0:
                        continue
                    c1, c2 = existing.coords_array(), g.coords_array()
                    for p in c1:
                        if np.min(np.linalg.norm(c2 - p, axis=1)) < 1.25:
                            clash = True
                            break
                    if clash:
                        break
                if not clash:
                    placed_guests.append(g)
                    success = True
                    break
            if not success:
                return None
    for pg in placed_guests:
        cluster.merge(pg)
    return cluster

