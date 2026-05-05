"""
prebiotic_chem.reactions
========================
Prebiotic reaction network definitions and pathway analysis tools.

Defines canonical prebiotic chemical reactions with:
  - Stoichiometry (reactants → products)
  - Prebiotic relevance description
  - Associated environmental scenario
  - Rough energy estimates when available

Also provides a simple directed graph representation for pathway tracing
(no external graph library required).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Set, Tuple


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class ReactionStep:
    """
    One elementary or net step in a prebiotic reaction pathway.

    Attributes
    ----------
    name:
        Unique identifier / short name.
    display_name:
        Human-readable reaction name.
    reactants:
        Molecule names (from ``molecules.py`` library or free text).
    products:
        Molecule names.
    conditions:
        Required environmental conditions (pH, T, catalyst, etc.).
    scenario_names:
        List of prebiotic scenario names where this step is relevant.
    description:
        Scientific narrative.
    delta_G_kcal:
        Approximate reaction free energy (kcal mol⁻¹); ``None`` if unknown.
    references:
        Key literature references.
    """

    name: str
    display_name: str
    reactants: Tuple[str, ...]
    products: Tuple[str, ...]
    conditions: str = ""
    scenario_names: Tuple[str, ...] = field(default_factory=tuple)
    description: str = ""
    delta_G_kcal: Optional[float] = None
    references: Tuple[str, ...] = field(default_factory=tuple)

    def __str__(self) -> str:  # pragma: no cover
        lhs = " + ".join(self.reactants)
        rhs = " + ".join(self.products)
        return f"{self.display_name}: {lhs} → {rhs}"


@dataclass
class ReactionPathway:
    """
    Ordered sequence of :class:`ReactionStep` objects forming a pathway.

    Parameters
    ----------
    name:
        Pathway identifier.
    display_name:
        Human-readable name.
    steps:
        Ordered list of reaction steps.
    description:
        Overview narrative.
    """

    name: str
    display_name: str
    steps: List[ReactionStep]
    description: str = ""

    def all_molecules(self) -> Set[str]:
        """Return the union of all reactant and product molecule names."""
        mols: Set[str] = set()
        for step in self.steps:
            mols.update(step.reactants)
            mols.update(step.products)
        return mols

    def summary(self) -> str:
        lines = [f"Pathway: {self.display_name}", "=" * 60]
        if self.description:
            lines.append(self.description)
            lines.append("")
        for i, step in enumerate(self.steps, 1):
            lhs = " + ".join(step.reactants)
            rhs = " + ".join(step.products)
            dG  = f"  ΔG ≈ {step.delta_G_kcal:.1f} kcal/mol" if step.delta_G_kcal is not None else ""
            lines.append(f"  Step {i}: {step.display_name}")
            lines.append(f"    {lhs} → {rhs}{dG}")
            if step.conditions:
                lines.append(f"    Conditions: {step.conditions}")
        return "\n".join(lines)


# ---------------------------------------------------------------------------
# Reaction step catalogue
# ---------------------------------------------------------------------------

# Strecker synthesis (amino acid from HCN, aldehyde, ammonia)
_S1 = ReactionStep(
    name="strecker_1_imine",
    display_name="Strecker Step 1 — imine formation",
    reactants=("formaldehyde", "ammonia"),
    products=("glycine_imine_intermediate",),
    conditions="pH 5–9, T = 25–80 °C",
    scenario_names=("warm_little_pond", "miller_urey", "tidal_pool"),
    description="Condensation of an aldehyde with ammonia to form an aldimine (Schiff base).",
    delta_G_kcal=-5.2,
    references=("Strecker 1850", "Miller 1955"),
)
_S2 = ReactionStep(
    name="strecker_2_aminonitrile",
    display_name="Strecker Step 2 — α-aminonitrile formation",
    reactants=("glycine_imine_intermediate", "hydrogen_cyanide"),
    products=("aminonitrile_intermediate",),
    conditions="pH 5–8, nucleophilic addition of HCN",
    scenario_names=("warm_little_pond", "miller_urey"),
    description="Addition of HCN to the imine to produce the α-aminonitrile.",
    delta_G_kcal=-8.5,
    references=("Miller 1955",),
)
_S3 = ReactionStep(
    name="strecker_3_hydrolysis",
    display_name="Strecker Step 3 — nitrile hydrolysis to amino acid",
    reactants=("aminonitrile_intermediate", "water"),
    products=("glycine",),
    conditions="pH 2–8, T = 50–100 °C, hours to days",
    scenario_names=("warm_little_pond", "miller_urey", "tidal_pool"),
    description="Acid/base-catalysed hydrolysis of the nitrile to the α-amino acid carboxylate.",
    delta_G_kcal=-12.0,
    references=("Miller & Orgel 1974",),
)

# Formose reaction (sugar synthesis)
_F1 = ReactionStep(
    name="formose_1_aldol",
    display_name="Formose Step 1 — glycolaldehyde formation",
    reactants=("formaldehyde", "formaldehyde"),
    products=("glycolaldehyde",),
    conditions="alkaline pH (≥ 9), T = 40–80 °C, Ca(OH)₂ catalyst",
    scenario_names=("warm_little_pond", "alkaline_vent"),
    description=(
        "The autocatalytic formose reaction initiates with aldol condensation of "
        "two formaldehyde molecules to give glycolaldehyde (C2 sugar)."
    ),
    delta_G_kcal=-4.0,
    references=("Butlerow 1861", "Breslow 1959"),
)
_F2 = ReactionStep(
    name="formose_2_glyceraldehyde",
    display_name="Formose Step 2 — glyceraldehyde formation",
    reactants=("glycolaldehyde", "formaldehyde"),
    products=("glyceraldehyde",),
    conditions="alkaline, autocatalytic",
    scenario_names=("warm_little_pond", "alkaline_vent"),
    description="Aldol addition extends the carbon chain to C3.",
    delta_G_kcal=-4.5,
    references=("Breslow 1959",),
)
_F3 = ReactionStep(
    name="formose_3_ribose",
    display_name="Formose Step 3 — ribose (C5 sugar) formation",
    reactants=("glyceraldehyde", "glycolaldehyde"),
    products=("ribose",),
    conditions="borate minerals (stabilises ribose), alkaline",
    scenario_names=("warm_little_pond",),
    description=(
        "Aldol condensation of C3 + C2 → C5; borate minerals selectively "
        "stabilise ribose over other pentoses (Ricardo et al. 2004, Science)."
    ),
    delta_G_kcal=-3.5,
    references=("Ricardo et al. 2004",),
)

# HCN polymerisation → adenine
_H1 = ReactionStep(
    name="hcn_poly_1_damn",
    display_name="HCN Polymerisation Step 1 — diaminomaleonitrile",
    reactants=("hydrogen_cyanide", "hydrogen_cyanide",
               "hydrogen_cyanide", "hydrogen_cyanide"),
    products=("diaminomaleonitrile",),
    conditions="pH 8–10, T = 0–25 °C, concentrated HCN (> 0.01 M)",
    scenario_names=("ice_eutectic", "warm_little_pond", "miller_urey"),
    description="HCN tetramer formation; key intermediate toward purine biosynthesis.",
    delta_G_kcal=-18.0,
    references=("Ferris & Orgel 1966",),
)
_H2 = ReactionStep(
    name="hcn_poly_2_adenine",
    display_name="HCN Polymerisation Step 2 — adenine synthesis",
    reactants=("diaminomaleonitrile", "hydrogen_cyanide"),
    products=("adenine",),
    conditions="UV irradiation or alkaline pH, T = 20–80 °C",
    scenario_names=("miller_urey", "warm_little_pond"),
    description="Pentamerisation of HCN (net) to give adenine (Oró 1960).",
    delta_G_kcal=-22.0,
    references=("Oró 1960", "Ferris & Orgel 1966"),
)

# Nucleotide condensation (nucleobase + ribose + phosphate)
_N1 = ReactionStep(
    name="nucleoside_formation",
    display_name="N-Glycosidic Bond Formation (nucleoside)",
    reactants=("adenine", "ribose"),
    products=("adenosine",),
    conditions="dry heating T = 100–160 °C, or formamide solvent, or mineral surface",
    scenario_names=("warm_little_pond", "volcanic_spring"),
    description=(
        "Formation of the N-glycosidic bond between purine/pyrimidine and ribose. "
        "Challenging prebiotically; facilitated by drying/heating or formamide."
    ),
    delta_G_kcal=4.5,
    references=("Fuller et al. 1972", "Becker et al. 2019"),
)
_N2 = ReactionStep(
    name="nucleotide_phosphorylation",
    display_name="Nucleotide Phosphorylation",
    reactants=("adenosine", "phosphoric_acid"),
    products=("adenosine_monophosphate",),
    conditions="dry heating, volcanic spring, or polyphosphate mineral source",
    scenario_names=("warm_little_pond", "volcanic_spring"),
    description=(
        "Phosphorylation of the 5′-hydroxyl of nucleoside to give the nucleotide. "
        "Driven by dry-heating or volcanic polyphosphate sources."
    ),
    delta_G_kcal=3.2,
    references=("Schwartz & Ponnamperuma 1968",),
)

# Iron–sulfur world (acetyl-CoA precursor)
_IS1 = ReactionStep(
    name="is_world_co2_reduction",
    display_name="Iron–Sulfur World — CO₂ reduction to acetate",
    reactants=("carbon_dioxide", "hydrogen_sulfide"),
    products=("acetic_acid",),
    conditions="FeS/Fe3S4 catalyst, T = 100–250 °C, acidic pH",
    scenario_names=("hydrothermal_vent",),
    description=(
        "FeS-catalysed reduction of CO₂ to acetic acid, the proposed precursor "
        "to acetyl-CoA in the iron–sulfur world (Wächtershäuser 1988)."
    ),
    delta_G_kcal=-7.0,
    references=("Wächtershäuser 1988", "Huber & Wächtershäuser 1997"),
)

# Peptide bond formation
_P1 = ReactionStep(
    name="peptide_condensation",
    display_name="Peptide Bond Formation (condensation)",
    reactants=("glycine", "alanine"),
    products=("gly_ala_dipeptide", "water"),
    conditions="dry-wet cycles, clay surface, or activated ester intermediate",
    scenario_names=("warm_little_pond", "tidal_pool", "volcanic_spring"),
    description=(
        "Condensation of two amino acids to form a dipeptide. Thermodynamically "
        "uphill in water (ΔG ≈ +3–4 kcal/mol); driven by dry-wet cycling or "
        "mineral surface catalysis."
    ),
    delta_G_kcal=3.5,
    references=("Fox 1958", "Lahav & Chang 1976"),
)

# Phospholipid self-assembly
_L1 = ReactionStep(
    name="vesicle_self_assembly",
    display_name="Fatty Acid Vesicle Formation",
    reactants=("decanoic_acid",),
    products=("fatty_acid_vesicle",),
    conditions="pH 7–9, T = 20–50 °C, dilute aqueous solution",
    scenario_names=("warm_little_pond", "tidal_pool"),
    description=(
        "Self-assembly of single-chain amphiphiles (fatty acids) into bilayer "
        "vesicles — protocell membrane precursors (Deamer 1985; Szostak 2001)."
    ),
    delta_G_kcal=-2.0,
    references=("Deamer 1985", "Szostak et al. 2001"),
)

# ---------------------------------------------------------------------------
# Reaction step catalogue dictionary
# ---------------------------------------------------------------------------

REACTION_STEPS: Dict[str, ReactionStep] = {
    step.name: step
    for step in [_S1, _S2, _S3, _F1, _F2, _F3, _H1, _H2, _N1, _N2, _IS1, _P1, _L1]
}

# ---------------------------------------------------------------------------
# Pathway catalogue
# ---------------------------------------------------------------------------

PATHWAYS: Dict[str, ReactionPathway] = {
    "strecker_amino_acid": ReactionPathway(
        name="strecker_amino_acid",
        display_name="Strecker Amino Acid Synthesis",
        steps=[_S1, _S2, _S3],
        description=(
            "Three-step synthesis of α-amino acids from simple feedstock molecules "
            "(aldehyde + NH₃ + HCN). Demonstrated in Miller-Urey sparking experiments."
        ),
    ),
    "formose_sugars": ReactionPathway(
        name="formose_sugars",
        display_name="Formose Reaction (Sugar Synthesis)",
        steps=[_F1, _F2, _F3],
        description=(
            "Autocatalytic condensation of formaldehyde to give pentose sugars "
            "including ribose. Mineral (borate) catalysis improves ribose selectivity."
        ),
    ),
    "hcn_to_adenine": ReactionPathway(
        name="hcn_to_adenine",
        display_name="HCN Polymerisation → Adenine",
        steps=[_H1, _H2],
        description=(
            "Prebiotic synthesis of adenine from HCN: tetramerisation → "
            "diaminomaleonitrile (DAMN) → pentamer → adenine."
        ),
    ),
    "nucleotide_synthesis": ReactionPathway(
        name="nucleotide_synthesis",
        display_name="Nucleotide Assembly (Nucleobase + Sugar + Phosphate)",
        steps=[_N1, _N2],
        description=(
            "Sequential formation of nucleoside (N-glycosidic bond) then "
            "nucleotide (phosphorylation) from prebiotic components."
        ),
    ),
    "iron_sulfur_acetate": ReactionPathway(
        name="iron_sulfur_acetate",
        display_name="Iron–Sulfur World: CO₂ → Acetate",
        steps=[_IS1],
        description=(
            "Proposed origin-of-metabolism pathway at hydrothermal vents via "
            "FeS-catalysed CO₂ reduction to acetic acid (Wächtershäuser)."
        ),
    ),
    "protocell_assembly": ReactionPathway(
        name="protocell_assembly",
        display_name="Protocell Assembly (Lipid Vesicle Formation)",
        steps=[_L1],
        description=(
            "Spontaneous self-assembly of prebiotic fatty acids into vesicles "
            "— the structural basis of early protocells."
        ),
    ),
}


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

def get_pathway(name: str) -> Optional[ReactionPathway]:
    """Return a pathway by name (case-insensitive).  Returns ``None`` if not found."""
    return PATHWAYS.get(name.strip().lower())


def get_reaction_step(name: str) -> Optional[ReactionStep]:
    """Return a reaction step by name (case-insensitive)."""
    return REACTION_STEPS.get(name.strip().lower())


def steps_for_scenario(scenario_name: str) -> List[ReactionStep]:
    """Return all reaction steps relevant to a given prebiotic scenario."""
    sn = scenario_name.strip().lower()
    return [s for s in REACTION_STEPS.values() if sn in s.scenario_names]


def pathways_for_molecule(molecule_name: str) -> List[ReactionPathway]:
    """Return pathways that involve a given molecule as reactant or product."""
    mn = molecule_name.strip().lower()
    result = []
    for pw in PATHWAYS.values():
        for step in pw.steps:
            if any(mn in r.lower() for r in step.reactants) or \
               any(mn in p.lower() for p in step.products):
                result.append(pw)
                break
    return result
