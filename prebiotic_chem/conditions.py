"""
prebiotic_chem.conditions
=========================
Prebiotic environment scenario definitions.

Each ``PrebioticScenario`` bundles the physical chemistry parameters that
characterise a canonical early-Earth environment (temperature, pH, ionic
strength, dominant solvent, atmospheric composition, etc.) and provides
helper methods to translate those conditions into xTB/CREST flags and ORCA
solvent keywords.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

from .constants import (
    PREBIOTIC_TEMP_RANGES,
    PREBIOTIC_PH_RANGES,
    PREBIOTIC_IONIC_STRENGTH,
    EARLY_EARTH_ATMOSPHERE,
    DEFAULT_XTB_METHOD,
    DEFAULT_CREST_METHOD,
)


# ---------------------------------------------------------------------------
# Scenario dataclass
# ---------------------------------------------------------------------------

@dataclass
class PrebioticScenario:
    """
    Bundles the physical conditions of a canonical prebiotic environment.

    Parameters
    ----------
    name:
        Identifier key (matches keys in ``constants.PREBIOTIC_TEMP_RANGES``).
    display_name:
        Human-readable name for reports.
    temperature_K:
        Representative temperature in Kelvin (mid-point of the range by default).
    pH:
        Representative pH.
    ionic_strength_M:
        Representative ionic strength in mol L⁻¹.
    solvent:
        Dominant solvent name (used to look up ORCA solvation keyword).
    atmosphere:
        Key into ``EARLY_EARTH_ATMOSPHERE`` dict; ``None`` means standard.
    description:
        Short scientific narrative.
    relevant_reactions:
        Comma-separated list of reaction family names relevant to this scenario.
    xtb_method:
        GFN-xTB level keyword (default GFN2-xTB).
    crest_method:
        CREST conformer-search method keyword.
    extra_xtb_flags:
        Additional command-line flags passed verbatim to xTB.
    """

    name: str
    display_name: str
    temperature_K: float
    pH: float
    ionic_strength_M: float
    solvent: str = "water"
    atmosphere: Optional[str] = None
    description: str = ""
    relevant_reactions: List[str] = field(default_factory=list)
    xtb_method: str = DEFAULT_XTB_METHOD
    crest_method: str = DEFAULT_CREST_METHOD
    extra_xtb_flags: List[str] = field(default_factory=list)

    # ------------------------------------------------------------------ #
    # Derived helpers                                                       #
    # ------------------------------------------------------------------ #

    @property
    def temperature_C(self) -> float:
        """Temperature in degrees Celsius."""
        return self.temperature_K - 273.15

    @property
    def temp_range_K(self) -> Tuple[float, float]:
        return PREBIOTIC_TEMP_RANGES.get(self.name, (self.temperature_K, self.temperature_K))

    @property
    def pH_range(self) -> Tuple[float, float]:
        return PREBIOTIC_PH_RANGES.get(self.name, (self.pH, self.pH))

    @property
    def ionic_strength_range_M(self) -> Tuple[float, float]:
        return PREBIOTIC_IONIC_STRENGTH.get(self.name, (self.ionic_strength_M, self.ionic_strength_M))

    def xtb_flags(self) -> List[str]:
        """
        Build the complete list of xTB command-line flags for this scenario.

        GFN2-xTB supports implicit ALPB solvation (--alpb <solvent>) and
        temperature via --etemp.  Ionic strength is approximated by the GBSA
        salt parameter when ALPB is unavailable.
        """
        flags: List[str] = [self.xtb_method]
        # Temperature
        flags += ["--etemp", f"{self.temperature_K:.2f}"]
        # Implicit solvation (ALPB — available for GFN2-xTB)
        if "gfn2" in self.xtb_method:
            flags += ["--alpb", self.solvent]
        elif "gfn1" in self.xtb_method:
            flags += ["--gbsa", self.solvent]
        # Any user-specified extras
        flags.extend(self.extra_xtb_flags)
        return flags

    def crest_flags(self) -> List[str]:
        """CREST flags for conformer sampling under this scenario."""
        flags: List[str] = [self.crest_method]
        flags += ["-T", str(int(self.temperature_K))]
        if "gfn2" in self.crest_method:
            flags += ["--alpb", self.solvent]
        return flags

    def orca_solvent_keyword(self) -> str:
        """
        Return the ORCA CPCM/SMD solvent keyword string.

        For example: ``CPCM(Water)`` or ``SMD(Water)``.
        """
        orca_name = _SOLVENT_MAP.get(self.solvent.lower(), self.solvent.capitalize())
        return f"CPCM({orca_name})"

    def summary(self) -> str:
        t_lo, t_hi = self.temp_range_K
        ph_lo, ph_hi = self.pH_range
        is_lo, is_hi = self.ionic_strength_range_M
        lines = [
            f"Scenario    : {self.display_name}",
            f"Temperature : {self.temperature_C:.1f} °C  "
            f"(range {t_lo - 273.15:.0f}–{t_hi - 273.15:.0f} °C)",
            f"pH          : {self.pH:.1f}  (range {ph_lo:.1f}–{ph_hi:.1f})",
            f"Ionic str.  : {self.ionic_strength_M:.3f} M  "
            f"(range {is_lo:.3f}–{is_hi:.3f} M)",
            f"Solvent     : {self.solvent}",
        ]
        if self.description:
            lines.append(f"Description : {self.description}")
        if self.relevant_reactions:
            lines.append(f"Reactions   : {', '.join(self.relevant_reactions)}")
        return "\n".join(lines)


# ---------------------------------------------------------------------------
# Solvent name map (IAK/ORCA name ← general name)
# ---------------------------------------------------------------------------

_SOLVENT_MAP: Dict[str, str] = {
    "water":    "Water",
    "h2o":      "Water",
    "methanol": "Methanol",
    "ethanol":  "Ethanol",
    "dmso":     "DMSO",
    "acetone":  "Acetone",
}


# ---------------------------------------------------------------------------
# Canonical scenario catalogue
# ---------------------------------------------------------------------------

_WARM_LITTLE_POND = PrebioticScenario(
    name="warm_little_pond",
    display_name="Warm Little Pond (Darwin's Pond)",
    temperature_K=333.15,   # ~60 °C representative
    pH=6.5,
    ionic_strength_M=0.15,
    solvent="water",
    atmosphere="modern_estimate",
    description=(
        "Shallow surface body of water subjected to wet–dry cycles. "
        "Concentrates organics by evaporation; allows condensation reactions "
        "during dry phase. Proposed by Darwin (1871) and championed by Deamer."
    ),
    relevant_reactions=["strecker", "peptide_condensation", "phosphorylation"],
)

_HYDROTHERMAL_VENT = PrebioticScenario(
    name="hydrothermal_vent",
    display_name="Submarine Black-Smoker Hydrothermal Vent",
    temperature_K=573.15,   # 300 °C representative
    pH=3.5,
    ionic_strength_M=0.5,
    solvent="water",
    atmosphere=None,
    description=(
        "High-temperature, acidic, sulfur-rich environment at mid-ocean ridges. "
        "Provides mineral catalysts (FeS, pyrite) and chemical gradients. "
        "Proposed by Wächtershäuser (1988) as the birthplace of metabolism."
    ),
    relevant_reactions=["iron_sulfur_world", "co2_reduction", "thioester_coupling"],
    extra_xtb_flags=["--chrg", "0"],
)

_ALKALINE_VENT = PrebioticScenario(
    name="alkaline_vent",
    display_name="Alkaline Hydrothermal Vent (White Smoker / Lost City)",
    temperature_K=373.15,   # ~100 °C representative
    pH=10.0,
    ionic_strength_M=0.2,
    solvent="water",
    description=(
        "Moderate-temperature, alkaline, hydrogen-rich environment. "
        "Natural proton gradients across thin inorganic membranes could drive "
        "chemiosmotic ATP synthesis (Russell & Martin 2004; Lane & Martin 2010)."
    ),
    relevant_reactions=["co2_reduction", "atp_synthesis", "amino_acid_formation"],
)

_ICE_EUTECTIC = PrebioticScenario(
    name="ice_eutectic",
    display_name="Ice-Eutectic / Frozen Lake",
    temperature_K=255.15,   # −18 °C representative
    pH=5.5,
    ionic_strength_M=2.0,
    solvent="water",
    atmosphere=None,
    description=(
        "Freeze concentration of solutes in eutectic brines within sea-ice. "
        "Greatly increases local reactant concentrations; stabilises RNA oligomers. "
        "Explored by Monnard, Trinks, and others."
    ),
    relevant_reactions=["hcn_polymerisation", "oligonucleotide_ligation"],
    extra_xtb_flags=["--etemp", "255.15"],
)

_TIDAL_POOL = PrebioticScenario(
    name="tidal_pool",
    display_name="Tidal Pool / Littoral Zone",
    temperature_K=318.15,   # ~45 °C representative
    pH=7.5,
    ionic_strength_M=0.4,
    solvent="water",
    description=(
        "Near-shore tidal environment with cycles of dilution and concentration. "
        "Mineral surfaces (clays, iron oxides) concentrate and catalyse reactions."
    ),
    relevant_reactions=["strecker", "phosphorylation", "nucleotide_condensation"],
)

_MILLER_UREY = PrebioticScenario(
    name="miller_urey",
    display_name="Miller–Urey Reducing Atmosphere",
    temperature_K=298.15,
    pH=7.0,
    ionic_strength_M=0.05,
    solvent="water",
    atmosphere="miller_urey_original",
    description=(
        "Classic 1953 experiment: reducing gas mixture (CH4, NH3, H2, H2O) "
        "subjected to electrical sparks. Produced amino acids, sugars, and "
        "nucleobase precursors. Remains the paradigm for abiogenesis experiments."
    ),
    relevant_reactions=["strecker", "formose", "amino_acid_formation"],
)

_VOLCANIC_SPRING = PrebioticScenario(
    name="volcanic_spring",
    display_name="Volcanic Hot Spring",
    temperature_K=363.15,   # 90 °C representative
    pH=4.0,
    ionic_strength_M=0.8,
    solvent="water",
    description=(
        "Geothermal springs rich in sulfur, metals, and CO2. "
        "Similar to modern Kamchatka, Yellowstone, or Taupo environments. "
        "Mineral-catalysed synthesis of amino acids and nucleobases possible."
    ),
    relevant_reactions=["strecker", "thioester_coupling", "hcn_polymerisation"],
)


# Public catalogue
SCENARIOS: Dict[str, PrebioticScenario] = {
    s.name: s
    for s in [
        _WARM_LITTLE_POND,
        _HYDROTHERMAL_VENT,
        _ALKALINE_VENT,
        _ICE_EUTECTIC,
        _TIDAL_POOL,
        _MILLER_UREY,
        _VOLCANIC_SPRING,
    ]
}


def get_scenario(name: str) -> Optional[PrebioticScenario]:
    """Return a scenario by name key (case-insensitive).  Returns ``None`` if not found."""
    return SCENARIOS.get(name.strip().lower())


def list_scenarios() -> List[PrebioticScenario]:
    """Return all predefined prebiotic scenarios."""
    return list(SCENARIOS.values())
