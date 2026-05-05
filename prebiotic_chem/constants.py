"""
prebiotic_chem.constants
========================
Physical constants, unit conversions, and early-Earth environmental parameters
used throughout the prebiotic chemistry computational workflow.
"""

from __future__ import annotations

# ---------------------------------------------------------------------------
# Fundamental physical constants (SI unless noted)
# ---------------------------------------------------------------------------
AVOGADRO    = 6.02214076e23          # mol⁻¹
BOLTZMANN_J = 1.380649e-23           # J K⁻¹
GAS_CONST_J = 8.314462618            # J mol⁻¹ K⁻¹
GAS_CONST_KCAL = GAS_CONST_J / 4184 # kcal mol⁻¹ K⁻¹
PLANCK      = 6.62607015e-34         # J s
SPEED_LIGHT = 2.99792458e10          # cm s⁻¹ (for wavenumber conversions)

# ---------------------------------------------------------------------------
# Unit conversions
# ---------------------------------------------------------------------------
HARTREE_TO_KCAL  = 627.509474       # kcal mol⁻¹ per Hartree
KCAL_TO_KJ       = 4.184            # kJ per kcal
HARTREE_TO_KJ    = HARTREE_TO_KCAL * KCAL_TO_KJ
EV_TO_HARTREE    = 0.036749405469   # Hartree per eV
BOHR_TO_ANGSTROM = 0.529177210903   # Å per Bohr
HARTREE_TO_EV    = 27.2113862       # eV per Hartree

# ---------------------------------------------------------------------------
# Standard conditions
# ---------------------------------------------------------------------------
STANDARD_TEMPERATURE_K = 298.15     # K (25 °C)
STANDARD_PRESSURE_BAR  = 1.01325    # bar

# ---------------------------------------------------------------------------
# Early-Earth / prebiotic environment ranges
# Reference: Deamer & Szostak (2010), Lazcano & Miller (1996),
#            Baaske et al. (2007), Martin et al. (2008)
# ---------------------------------------------------------------------------

# Temperature ranges for canonical prebiotic scenarios (Kelvin)
PREBIOTIC_TEMP_RANGES: dict[str, tuple[float, float]] = {
    "warm_little_pond":      (293.15, 363.15),   # 20–90 °C
    "hydrothermal_vent":     (323.15, 673.15),   # 50–400 °C
    "alkaline_vent":         (323.15, 423.15),   # 50–150 °C
    "ice_eutectic":          (243.15, 273.15),   # −30–0 °C
    "tidal_pool":            (288.15, 353.15),   # 15–80 °C
    "deep_ocean":            (274.15, 288.15),   # 1–15 °C
    "volcanic_spring":       (333.15, 423.15),   # 60–150 °C
    "miller_urey":           (298.15, 373.15),   # 25–100 °C
}

# pH ranges for canonical prebiotic scenarios
PREBIOTIC_PH_RANGES: dict[str, tuple[float, float]] = {
    "warm_little_pond":      (5.0, 8.0),
    "hydrothermal_vent":     (2.0, 5.0),
    "alkaline_vent":         (9.0, 11.0),
    "ice_eutectic":          (4.0, 7.0),
    "tidal_pool":            (6.0, 9.0),
    "deep_ocean":            (7.0, 8.5),
    "volcanic_spring":       (2.0, 6.0),
    "miller_urey":           (6.0, 8.0),
}

# Ionic-strength ranges (mol L⁻¹) — dominated by NaCl, Mg²⁺, Ca²⁺
PREBIOTIC_IONIC_STRENGTH: dict[str, tuple[float, float]] = {
    "warm_little_pond":      (0.01, 0.5),
    "hydrothermal_vent":     (0.05, 1.0),
    "alkaline_vent":         (0.05, 0.5),
    "ice_eutectic":          (0.5,  5.0),
    "tidal_pool":            (0.1,  1.0),
    "deep_ocean":            (0.5,  0.7),
    "volcanic_spring":       (0.1,  2.0),
    "miller_urey":           (0.0,  0.1),
}

# Representative atmospheric composition (mixing ratios, early Earth)
# Miller-Urey (1953) and modern estimates post-Kasting (1993)
EARLY_EARTH_ATMOSPHERE: dict[str, dict[str, float]] = {
    "miller_urey_original": {
        "H2":  0.20, "CH4": 0.20, "NH3": 0.10, "H2O": 0.50,
    },
    "modern_estimate": {
        "N2":  0.80, "CO2": 0.10, "H2O": 0.09, "H2":  0.01,
    },
    "co2_rich": {
        "CO2": 0.70, "N2":  0.20, "H2O": 0.08, "H2S": 0.02,
    },
}

# ---------------------------------------------------------------------------
# xTB / CREST method keywords for prebiotic conditions
# ---------------------------------------------------------------------------
# GFN2-xTB is the recommended semi-empirical method for prebiotic molecules
# (good balance of accuracy and cost for H, C, N, O, P, S, halogens)
DEFAULT_XTB_METHOD   = "--gfn2"
DEFAULT_CREST_METHOD = "--gfn2"
DEFAULT_ORCA_METHOD  = "B97-3c Opt Freq"   # fast and reliable composite DFT

# Energy window for conformer retention (kcal mol⁻¹)
DEFAULT_XTB_EWIN     = 5.0    # kcal mol⁻¹  (ChemRefine-inspired)
DEFAULT_CREST_EWIN   = 3.0    # kcal mol⁻¹
DEFAULT_TEMPERATURE  = 298.15 # K — standard for Boltzmann analysis

# Elements commonly found in prebiotic molecules
PREBIOTIC_ELEMENTS = frozenset({"H", "C", "N", "O", "P", "S", "Na", "Mg", "Ca", "K", "Cl", "Fe"})

# Periodic table (symbol → atomic number)
_SYMBOLS = [
    "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne",
    "Na", "Mg", "Al", "Si", "P",  "S",  "Cl", "Ar", "K",  "Ca",
    "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y",  "Zr",
    "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
    "Sb", "Te", "I",  "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
    "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
    "Lu", "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
    "Pa", "U",  "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm",
    "Md", "No", "Lr",
]
ATOMIC_NUMBERS: dict[str, int] = {sym: idx + 1 for idx, sym in enumerate(_SYMBOLS)}
