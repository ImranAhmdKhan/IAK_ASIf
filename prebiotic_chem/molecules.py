"""
prebiotic_chem.molecules
========================
Curated library of prebiotic molecules organised by chemical family.
Each entry provides a SMILES string, molecular name, molecular formula,
and contextual information about its prebiotic relevance.

A helper ``get_molecule`` function looks up by common name or alias, and
``list_molecules`` returns all entries for a given category.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional


@dataclass(frozen=True)
class PrebioticMolecule:
    """Metadata record for a prebiotic-chemistry-relevant molecule."""

    name: str
    smiles: str
    formula: str
    category: str
    charge: int = 0
    multiplicity: int = 1
    aliases: tuple[str, ...] = field(default_factory=tuple)
    relevance: str = ""

    def __str__(self) -> str:  # pragma: no cover
        return f"{self.name} ({self.formula})"


# ---------------------------------------------------------------------------
# Molecule library
# ---------------------------------------------------------------------------
_LIBRARY: List[PrebioticMolecule] = [

    # ── Simple gases / Miller-Urey feedstock ──────────────────────────────
    PrebioticMolecule(
        name="hydrogen_cyanide", smiles="C#N", formula="HCN",
        category="feedstock",
        aliases=("hcn", "hydrocyanic acid"),
        relevance=(
            "Central feedstock in prebiotic chemistry. Polymerises to nucleobases "
            "(Oró 1960), amino acids (Strecker synthesis), and purines. Present in "
            "comets, meteorites, and early-Earth volcanic/lightning discharges."
        ),
    ),
    PrebioticMolecule(
        name="formaldehyde", smiles="C=O", formula="CH2O",
        category="feedstock",
        aliases=("methanal", "hcho"),
        relevance=(
            "Key sugar precursor (formose reaction, Butlerow 1861). Produced by "
            "UV photolysis of CO2/CH4 in early atmosphere. Found in comets."
        ),
    ),
    PrebioticMolecule(
        name="ammonia", smiles="N", formula="NH3",
        category="feedstock",
        aliases=("nh3",),
        relevance=(
            "Nitrogen source in Strecker amino-acid synthesis and nucleobase "
            "formation. Present in Miller-Urey reducing atmosphere."
        ),
    ),
    PrebioticMolecule(
        name="water", smiles="O", formula="H2O",
        category="solvent",
        aliases=("h2o",),
        relevance="Universal prebiotic solvent.",
    ),
    PrebioticMolecule(
        name="carbon_dioxide", smiles="O=C=O", formula="CO2",
        category="feedstock",
        aliases=("co2",),
        relevance="Dominant carbon source in modern estimate of early atmosphere.",
    ),
    PrebioticMolecule(
        name="methane", smiles="C", formula="CH4",
        category="feedstock",
        aliases=("ch4",),
        relevance="Carbon source in Miller-Urey scenario.",
    ),
    PrebioticMolecule(
        name="hydrogen_sulfide", smiles="S", formula="H2S",
        category="feedstock",
        aliases=("h2s",),
        relevance=(
            "Sulfur source; involved in origin of metabolism at hydrothermal vents "
            "(Wächtershäuser 1990). Activates nucleotides in RNA world chemistry."
        ),
    ),
    PrebioticMolecule(
        name="cyanamide", smiles="N#CN", formula="CH2N2",
        category="feedstock",
        aliases=("h2ncn",),
        relevance=(
            "Condensing agent and activated amino-acid precursor. Formed from HCN "
            "hydrolysis and involved in peptide bond formation."
        ),
    ),
    PrebioticMolecule(
        name="phosphoric_acid", smiles="OP(=O)(O)O", formula="H3PO4",
        category="feedstock",
        aliases=("phosphate", "h3po4"),
        relevance=(
            "Phosphate source for nucleotide synthesis and energy coupling. "
            "Derived from volcanic activity and meteorite delivery."
        ),
    ),

    # ── Amino acids ──────────────────────────────────────────────────────
    PrebioticMolecule(
        name="glycine", smiles="NCC(=O)O", formula="C2H5NO2",
        category="amino_acid",
        aliases=("gly", "aminoacetic acid"),
        relevance=(
            "Simplest amino acid; detected in Murchison meteorite and cometary "
            "samples. Synthesised in Miller-Urey and Strecker reactions."
        ),
    ),
    PrebioticMolecule(
        name="alanine", smiles="C[C@@H](N)C(=O)O", formula="C3H7NO2",
        category="amino_acid",
        aliases=("ala",),
        relevance=(
            "Second simplest amino acid; abundant in Miller-Urey experiments "
            "and in carbonaceous chondrites."
        ),
    ),
    PrebioticMolecule(
        name="serine", smiles="N[C@@H](CO)C(=O)O", formula="C3H7NO3",
        category="amino_acid",
        aliases=("ser",),
        relevance="Hydroxyl amino acid; key intermediate in RNA world phosphorylation.",
    ),
    PrebioticMolecule(
        name="aspartic_acid", smiles="N[C@@H](CC(=O)O)C(=O)O", formula="C4H7NO4",
        category="amino_acid",
        aliases=("asp", "aspartate"),
        relevance=(
            "Acidic amino acid; found in meteorites; intermediate in purine biosynthesis."
        ),
    ),
    PrebioticMolecule(
        name="glutamic_acid", smiles="N[C@@H](CCC(=O)O)C(=O)O", formula="C5H9NO4",
        category="amino_acid",
        aliases=("glu", "glutamate"),
        relevance="Found in carbonaceous meteorites; important for early metabolism.",
    ),

    # ── Nucleobases ───────────────────────────────────────────────────────
    PrebioticMolecule(
        name="adenine", smiles="Nc1ncnc2[nH]cnc12", formula="C5H5N5",
        category="nucleobase",
        aliases=("ade", "6-aminopurine"),
        relevance=(
            "Purine base of ATP/ADP and DNA/RNA adenosine. Formed by pentamerisation "
            "of HCN (Oró 1960) and from interstellar ice photochemistry."
        ),
    ),
    PrebioticMolecule(
        name="guanine", smiles="Nc1nc2[nH]cnc2c(=O)[nH]1", formula="C5H5N5O",
        category="nucleobase",
        aliases=("gua",),
        relevance=(
            "Purine base; found in Murchison meteorite. Formed from HCN polymers "
            "and diaminomaleonitrile (DAMN) intermediates."
        ),
    ),
    PrebioticMolecule(
        name="cytosine", smiles="Nc1ccnc(=O)[nH]1", formula="C4H5N3O",
        category="nucleobase",
        aliases=("cyt",),
        relevance=(
            "Pyrimidine base; synthesised from cyanoacetylene + urea under prebiotic "
            "conditions (Ferris & Orgel 1966). Present in icy-body chemistry."
        ),
    ),
    PrebioticMolecule(
        name="uracil", smiles="O=c1ccn([H])c(=O)[nH]1", formula="C4H4N2O2",
        category="nucleobase",
        aliases=("ura",),
        relevance=(
            "RNA pyrimidine base; detected in Murchison meteorite and formed by "
            "UV photolysis of cytosine aqueous solutions."
        ),
    ),
    PrebioticMolecule(
        name="thymine", smiles="Cc1c[nH]c(=O)[nH]c1=O", formula="C5H6N2O2",
        category="nucleobase",
        aliases=("thy", "5-methyluracil"),
        relevance="DNA-specific pyrimidine; formed from uracil methylation.",
    ),
    PrebioticMolecule(
        name="hypoxanthine", smiles="O=c1[nH]cnc2[nH]cnc12", formula="C5H4N4O",
        category="nucleobase",
        aliases=("hyp",),
        relevance=(
            "Purine found in Murchison meteorite; may have preceded adenine in early "
            "RNA-like polymers (Szostak hypothesis)."
        ),
    ),

    # ── Sugars ────────────────────────────────────────────────────────────
    PrebioticMolecule(
        name="ribose", smiles="OC[C@H]1O[C@@H](O)[C@H](O)[C@@H]1O", formula="C5H10O5",
        category="sugar",
        aliases=("d-ribose", "ribo"),
        relevance=(
            "Backbone sugar of RNA. Synthesised by formose reaction; selectively "
            "enriched by mineral (borate) catalysis (Ricardo et al. 2004)."
        ),
    ),
    PrebioticMolecule(
        name="deoxyribose", smiles="OC[C@H]1O[C@@H](O)C[C@@H]1O", formula="C5H10O4",
        category="sugar",
        aliases=("2-deoxyribose", "dribo"),
        relevance="Backbone sugar of DNA; derived from ribose by enzymatic reduction in biology.",
    ),
    PrebioticMolecule(
        name="glycolaldehyde", smiles="OCC=O", formula="C2H4O2",
        category="sugar",
        aliases=("glycoaldehyde", "ga"),
        relevance=(
            "Simplest sugar aldehyde (C2); detected in interstellar medium. "
            "First product of formose reaction with self-condensation of formaldehyde."
        ),
    ),
    PrebioticMolecule(
        name="glyceraldehyde", smiles="OC[C@@H](O)C=O", formula="C3H6O3",
        category="sugar",
        aliases=("glyceral",),
        relevance=(
            "C3 sugar; central intermediate of formose reaction and glycolysis. "
            "Detected in carbonaceous meteorites."
        ),
    ),

    # ── Nucleotides / nucleosides ─────────────────────────────────────────
    PrebioticMolecule(
        name="adenosine_monophosphate", formula="C10H14N5O7P",
        smiles="Nc1ncnc2n(cnc12)[C@@H]1O[C@H](COP(=O)(O)O)[C@@H](O)[C@H]1O",
        category="nucleotide",
        aliases=("amp", "adenosine-5-monophosphate"),
        relevance=(
            "AMP; building block for ATP and RNA. Synthesised prebiotically "
            "by nucleobase + ribose + phosphate condensation."
        ),
    ),
    PrebioticMolecule(
        name="adenosine_triphosphate", formula="C10H16N5O13P3",
        smiles="Nc1ncnc2n(cnc12)[C@@H]1O[C@H](COP(=O)(O)OP(=O)(O)OP(=O)(O)O)[C@@H](O)[C@H]1O",
        category="nucleotide",
        aliases=("atp",),
        relevance=(
            "Universal energy currency; may have originated at hydrothermal vents "
            "via condensation of AMP with polyphosphate."
        ),
    ),

    # ── Lipid / membrane precursors ───────────────────────────────────────
    PrebioticMolecule(
        name="decanoic_acid", smiles="CCCCCCCCCC(=O)O", formula="C10H20O2",
        category="lipid",
        aliases=("capric_acid", "c10_fatty_acid"),
        relevance=(
            "C10 fatty acid; forms stable vesicles at mildly acidic pH. "
            "Detected in Murchison meteorite (Deamer 1985)."
        ),
    ),
    PrebioticMolecule(
        name="oleic_acid", smiles="CCCCCCCCC=CCCCCCCCC(=O)O", formula="C18H34O2",
        category="lipid",
        aliases=("c18:1", "9-octadecenoic acid"),
        relevance=(
            "Monounsaturated fatty acid; relevant for protocell membrane formation "
            "under hydrothermal conditions."
        ),
    ),

    # ── Cofactor / metabolite precursors ──────────────────────────────────
    PrebioticMolecule(
        name="nicotinamide", smiles="NC(=O)c1ccncc1", formula="C6H6N2O",
        category="cofactor",
        aliases=("niacinamide", "vitamin_b3_amide"),
        relevance=(
            "Amide form of nicotinic acid; precursor to NAD⁺. Synthesised from "
            "β-alanine and aspartate under prebiotic conditions."
        ),
    ),
    PrebioticMolecule(
        name="acetic_acid", smiles="CC(=O)O", formula="C2H4O2",
        category="metabolite",
        aliases=("acetate", "ethanoic_acid"),
        relevance=(
            "Simplest carboxylic acid thioester analogue; central to iron-sulfur-world "
            "acetyl-CoA pathway hypothesis (Wächtershäuser)."
        ),
    ),
    PrebioticMolecule(
        name="pyruvic_acid", smiles="CC(=O)C(=O)O", formula="C3H4O3",
        category="metabolite",
        aliases=("pyruvate",),
        relevance=(
            "Keto acid; hub of early metabolism. Formed from oxidative decarboxylation "
            "of malate in the reverse TCA cycle."
        ),
    ),

    # ── Mineral surface proxies ───────────────────────────────────────────
    PrebioticMolecule(
        name="silicic_acid", smiles="O[Si](O)(O)O", formula="H4SiO4",
        category="mineral_proxy",
        aliases=("orthosilicic acid",),
        relevance=(
            "Represents clay/silica surface chemistry. Adsorbs amino acids and "
            "nucleotides; promotes condensation reactions."
        ),
    ),
    PrebioticMolecule(
        name="boric_acid", smiles="OB(O)O", formula="H3BO3",
        category="mineral_proxy",
        aliases=("h3bo3", "boron"),
        relevance=(
            "Stabilises ribose over other pentoses in formose reaction "
            "(Ricardo et al. 2004, Science)."
        ),
    ),
]

# ---------------------------------------------------------------------------
# Lookup helpers
# ---------------------------------------------------------------------------
_NAME_INDEX: Dict[str, PrebioticMolecule] = {}
for _mol in _LIBRARY:
    _NAME_INDEX[_mol.name.lower()] = _mol
    for _alias in _mol.aliases:
        _NAME_INDEX[_alias.lower()] = _mol


def get_molecule(name: str) -> Optional[PrebioticMolecule]:
    """Return a :class:`PrebioticMolecule` by name or alias (case-insensitive)."""
    return _NAME_INDEX.get(name.strip().lower())


def list_molecules(category: Optional[str] = None) -> List[PrebioticMolecule]:
    """Return all molecules, optionally filtered by *category*."""
    if category is None:
        return list(_LIBRARY)
    cat = category.strip().lower()
    return [m for m in _LIBRARY if m.category.lower() == cat]


def list_categories() -> List[str]:
    """Return sorted list of unique category names."""
    return sorted({m.category for m in _LIBRARY})
