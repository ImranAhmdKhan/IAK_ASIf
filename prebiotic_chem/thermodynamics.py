"""
prebiotic_chem.thermodynamics
=============================
Thermodynamic analysis tools for prebiotic chemistry calculations.

Provides:
- Boltzmann population analysis of conformer ensembles
- Free-energy calculations from electronic energies + thermochemical corrections
- Reaction free-energy estimation
- Temperature-dependent equilibrium constant calculations
"""

from __future__ import annotations

import math
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from .constants import (
    GAS_CONST_KCAL,
    HARTREE_TO_KCAL,
    KCAL_TO_KJ,
    STANDARD_TEMPERATURE_K,
)


# ---------------------------------------------------------------------------
# Boltzmann analysis
# ---------------------------------------------------------------------------

def boltzmann_weights(
    energies_hartree: Sequence[float],
    temperature_K: float = STANDARD_TEMPERATURE_K,
) -> np.ndarray:
    """
    Compute Boltzmann weights for a set of electronic energies.

    Parameters
    ----------
    energies_hartree:
        Absolute electronic energies in Hartree.
    temperature_K:
        Temperature in Kelvin.

    Returns
    -------
    np.ndarray
        Normalised Boltzmann weights (sums to 1.0).
    """
    e = np.asarray(energies_hartree, dtype=float)
    delta_E_kcal = (e - e.min()) * HARTREE_TO_KCAL   # relative, ≥ 0
    RT = GAS_CONST_KCAL * temperature_K
    w = np.exp(-delta_E_kcal / RT)
    return w / w.sum()


def boltzmann_populations(
    energies_hartree: Sequence[float],
    ids: Optional[Sequence[str]] = None,
    temperature_K: float = STANDARD_TEMPERATURE_K,
) -> List[Dict]:
    """
    Return a list of dicts with energy, relative energy, and population data.

    Parameters
    ----------
    energies_hartree:
        Absolute energies in Hartree.
    ids:
        Optional structure identifiers.
    temperature_K:
        Temperature in Kelvin.

    Returns
    -------
    List of dicts sorted by ascending energy.
    Each dict contains:
        ``id``, ``energy_hartree``, ``rel_energy_kcal``,
        ``rel_energy_kJ``, ``boltzmann_weight``, ``population_pct``.
    """
    e = np.asarray(energies_hartree, dtype=float)
    if ids is None:
        ids = [str(i) for i in range(len(e))]

    order = np.argsort(e)
    e_sorted = e[order]
    ids_sorted = [ids[i] for i in order]

    delta_E_kcal = (e_sorted - e_sorted[0]) * HARTREE_TO_KCAL
    delta_E_kJ   = delta_E_kcal * KCAL_TO_KJ

    RT = GAS_CONST_KCAL * temperature_K
    w  = np.exp(-delta_E_kcal / RT)
    w /= w.sum()

    results = []
    for i, idx in enumerate(order):
        results.append({
            "id":                str(ids_sorted[i]),
            "energy_hartree":    float(e_sorted[i]),
            "rel_energy_kcal":   float(delta_E_kcal[i]),
            "rel_energy_kJ":     float(delta_E_kJ[i]),
            "boltzmann_weight":  float(w[i]),
            "population_pct":    float(w[i] * 100.0),
        })
    return results


def cumulative_population_cutoff(
    populations: List[Dict],
    threshold_pct: float = 99.0,
) -> List[Dict]:
    """
    Return structures up to a cumulative Boltzmann population threshold.

    Parameters
    ----------
    populations:
        Output of :func:`boltzmann_populations`.
    threshold_pct:
        Cumulative population percentage cutoff (default 99 %).

    Returns
    -------
    Filtered list (structures with highest populations first).
    """
    cumulative = 0.0
    selected = []
    for entry in populations:
        selected.append(entry)
        cumulative += entry["population_pct"]
        if cumulative >= threshold_pct:
            break
    return selected


# ---------------------------------------------------------------------------
# Free energy helpers
# ---------------------------------------------------------------------------

def gibbs_from_enthalpy_entropy(
    enthalpy_kcal: float,
    entropy_cal_mol_K: float,
    temperature_K: float = STANDARD_TEMPERATURE_K,
) -> float:
    """
    Compute Gibbs free energy: G = H − T·S.

    Parameters
    ----------
    enthalpy_kcal:
        Enthalpy in kcal mol⁻¹.
    entropy_cal_mol_K:
        Entropy in cal mol⁻¹ K⁻¹ (standard ORCA/Gaussian units).
    temperature_K:
        Temperature in Kelvin.

    Returns
    -------
    float
        Gibbs free energy in kcal mol⁻¹.
    """
    return enthalpy_kcal - temperature_K * (entropy_cal_mol_K / 1000.0)


def reaction_free_energy(
    delta_H_kcal: float,
    delta_S_cal_mol_K: float,
    temperature_K: float = STANDARD_TEMPERATURE_K,
) -> Tuple[float, float]:
    """
    Compute ΔG_rxn and equilibrium constant K_eq.

    Parameters
    ----------
    delta_H_kcal:
        Reaction enthalpy in kcal mol⁻¹.
    delta_S_cal_mol_K:
        Reaction entropy in cal mol⁻¹ K⁻¹.
    temperature_K:
        Temperature in Kelvin.

    Returns
    -------
    (delta_G_kcal, K_eq)
    """
    delta_G = gibbs_from_enthalpy_entropy(delta_H_kcal, delta_S_cal_mol_K, temperature_K)
    RT = GAS_CONST_KCAL * temperature_K
    K_eq = math.exp(-delta_G / RT) if RT > 0 else float("inf")
    return delta_G, K_eq


# ---------------------------------------------------------------------------
# Energy-window filtering (mirrors ChemRefine StructureRefiner.filter)
# ---------------------------------------------------------------------------

def filter_energy_window(
    energies_hartree: Sequence[float],
    ids: Optional[Sequence[str]] = None,
    window_kcal: float = 5.0,
) -> List[int]:
    """
    Return indices of structures within *window_kcal* of the global minimum.

    Parameters
    ----------
    energies_hartree:
        Energies in Hartree.
    ids:
        Optional structure identifiers (unused in filtering, kept for API symmetry).
    window_kcal:
        Energy window in kcal mol⁻¹ above the minimum.

    Returns
    -------
    Sorted list of indices satisfying the energy window.
    """
    e = np.asarray(energies_hartree, dtype=float)
    e_min = e.min()
    delta = (e - e_min) * HARTREE_TO_KCAL
    return sorted(int(i) for i in np.where(delta <= window_kcal)[0])


def filter_by_boltzmann_cutoff(
    energies_hartree: Sequence[float],
    temperature_K: float = STANDARD_TEMPERATURE_K,
    cumulative_pct: float = 99.0,
) -> List[int]:
    """
    Return indices of structures that together account for *cumulative_pct* of
    the Boltzmann-weighted ensemble population.

    Structures are sorted by energy (lowest first) before applying the cutoff.

    Parameters
    ----------
    energies_hartree:
        Energies in Hartree.
    temperature_K:
        Temperature in Kelvin.
    cumulative_pct:
        Cumulative population threshold in percent.

    Returns
    -------
    List of original indices (unsorted order matches ``energies_hartree``).
    """
    e = np.asarray(energies_hartree, dtype=float)
    order = np.argsort(e)

    w = boltzmann_weights(e, temperature_K)
    w_sorted = w[order]

    cumulative = 0.0
    selected_positions = []
    for pos, orig_idx in enumerate(order):
        selected_positions.append(int(orig_idx))
        cumulative += w_sorted[pos] * 100.0
        if cumulative >= cumulative_pct:
            break
    return selected_positions


# ---------------------------------------------------------------------------
# Temperature sweep
# ---------------------------------------------------------------------------

def temperature_sweep(
    delta_H_kcal: float,
    delta_S_cal_mol_K: float,
    temperatures_K: Optional[Sequence[float]] = None,
) -> List[Dict[str, float]]:
    """
    Compute ΔG and K_eq across a range of temperatures.

    Useful for assessing thermodynamic feasibility of a prebiotic reaction
    across different environmental temperature scenarios.

    Parameters
    ----------
    delta_H_kcal:
        Reaction enthalpy in kcal mol⁻¹.
    delta_S_cal_mol_K:
        Reaction entropy in cal mol⁻¹ K⁻¹.
    temperatures_K:
        Temperatures to evaluate. Defaults to 250–500 K in 25 K steps.

    Returns
    -------
    List of dicts with keys ``temperature_K``, ``temperature_C``,
    ``delta_G_kcal``, ``K_eq``, ``spontaneous``.
    """
    if temperatures_K is None:
        temperatures_K = list(range(250, 725, 25))

    rows = []
    for T in temperatures_K:
        dG, Keq = reaction_free_energy(delta_H_kcal, delta_S_cal_mol_K, float(T))
        rows.append({
            "temperature_K":  float(T),
            "temperature_C":  float(T) - 273.15,
            "delta_G_kcal":   dG,
            "K_eq":           Keq,
            "spontaneous":    dG < 0,
        })
    return rows
