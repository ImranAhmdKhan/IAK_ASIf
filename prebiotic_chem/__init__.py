"""
prebiotic_chem
==============
Professional computational chemistry app focused on prebiotic conditions.

Provides a complete workflow for studying prebiotic molecules and reactions
using state-of-the-art semi-empirical and DFT methods:

  - Curated library of prebiotic molecules (HCN, amino acids, nucleobases,
    sugars, nucleotides, lipids, and cofactor precursors)
  - Canonical early-Earth environment scenarios (warm little pond,
    hydrothermal vents, ice eutectic, Miller–Urey, and more)
  - Prebiotic reaction network definitions and pathway analysis
  - Multi-step computational workflow:
      xTB pre-opt → CREST conformer sampling → filtering →
      ORCA DFT geometry optimisation → ORCA high-level single point
  - Boltzmann population analysis, energy-window filtering, RMSD clustering
  - Temperature-sweep thermodynamics for reaction feasibility assessment
  - Plain-text and CSV report generation

Inspired by the ChemRefine framework (Sterling Group, University of Texas
at Dallas — https://github.com/sterling-group/ChemRefine).

Quick start::

    python -m prebiotic_chem --list-scenarios
    python -m prebiotic_chem --list-molecules
    python -m prebiotic_chem --scenario warm_little_pond --info
    python -m prebiotic_chem --xyz glycine.xyz --scenario warm_little_pond --run-dir run1
"""

from __future__ import annotations

__version__ = "1.0.0"
__author__  = "IAK Prebiotic Chemistry"

from .constants     import (
    HARTREE_TO_KCAL, GAS_CONST_KCAL, STANDARD_TEMPERATURE_K,
    PREBIOTIC_TEMP_RANGES, PREBIOTIC_PH_RANGES,
)
from .molecules     import (
    PrebioticMolecule, get_molecule, list_molecules, list_categories,
)
from .conditions    import (
    PrebioticScenario, SCENARIOS, get_scenario, list_scenarios,
)
from .reactions     import (
    ReactionStep, ReactionPathway, PATHWAYS, REACTION_STEPS,
    get_pathway, get_reaction_step, steps_for_scenario, pathways_for_molecule,
)
from .thermodynamics import (
    boltzmann_weights, boltzmann_populations, cumulative_population_cutoff,
    filter_energy_window, filter_by_boltzmann_cutoff,
    gibbs_from_enthalpy_entropy, reaction_free_energy, temperature_sweep,
)
from .analysis      import (
    kabsch_rmsd, rmsd_matrix, rmsd_cluster,
    detect_hbonds, find_hbond_donors, find_hbond_acceptors,
    nucleophilic_sites, electrophilic_sites,
    score_prebiotic_geometry, molecular_formula,
)
from .io_utils      import (
    read_xyz, write_xyz, read_xyz_trajectory,
    smiles_to_xyz, run_xtb, run_crest,
    write_orca_input, parse_orca_energy, parse_orca_thermochemistry,
)
from .pipeline      import PrebioticPipeline
from .report        import (
    PrebioticReport, hbond_report, temperature_sweep_report, write_population_csv,
)

__all__ = [
    # constants
    "HARTREE_TO_KCAL", "GAS_CONST_KCAL", "STANDARD_TEMPERATURE_K",
    "PREBIOTIC_TEMP_RANGES", "PREBIOTIC_PH_RANGES",
    # molecules
    "PrebioticMolecule", "get_molecule", "list_molecules", "list_categories",
    # conditions
    "PrebioticScenario", "SCENARIOS", "get_scenario", "list_scenarios",
    # reactions
    "ReactionStep", "ReactionPathway", "PATHWAYS", "REACTION_STEPS",
    "get_pathway", "get_reaction_step", "steps_for_scenario", "pathways_for_molecule",
    # thermodynamics
    "boltzmann_weights", "boltzmann_populations", "cumulative_population_cutoff",
    "filter_energy_window", "filter_by_boltzmann_cutoff",
    "gibbs_from_enthalpy_entropy", "reaction_free_energy", "temperature_sweep",
    # analysis
    "kabsch_rmsd", "rmsd_matrix", "rmsd_cluster",
    "detect_hbonds", "find_hbond_donors", "find_hbond_acceptors",
    "nucleophilic_sites", "electrophilic_sites",
    "score_prebiotic_geometry", "molecular_formula",
    # io
    "read_xyz", "write_xyz", "read_xyz_trajectory",
    "smiles_to_xyz", "run_xtb", "run_crest",
    "write_orca_input", "parse_orca_energy", "parse_orca_thermochemistry",
    # pipeline
    "PrebioticPipeline",
    # report
    "PrebioticReport", "hbond_report", "temperature_sweep_report",
    "write_population_csv",
]
