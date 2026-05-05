# prebiotic_chem

> **Professional computational chemistry app focused on prebiotic conditions**

`prebiotic_chem` is a self-contained Python package that brings a **professional
computational chemist's workflow** to the study of early-Earth prebiotic chemistry.
It is inspired by the [ChemRefine](https://github.com/sterling-group/ChemRefine)
framework (Sterling Group, University of Texas at Dallas) and integrates seamlessly
with the **xTB**, **CREST**, and **ORCA** quantum-chemistry codes used throughout
the IAK_ASIf repository.

---

## Features

| Feature | Description |
|---|---|
| **Molecule library** | 30+ curated prebiotic molecules (HCN, amino acids, nucleobases, sugars, nucleotides, lipids, cofactors) with SMILES strings and prebiotic relevance notes |
| **Scenario library** | 7 canonical early-Earth environments: warm little pond, hydrothermal vent, alkaline vent, ice eutectic, tidal pool, Miller–Urey, volcanic spring |
| **Reaction network** | 13 key prebiotic reaction steps organised into 6 pathways (Strecker, formose, HCN→adenine, nucleotide assembly, iron-sulfur world, protocell vesicle formation) |
| **Computational pipeline** | 5-step workflow: xTB pre-opt → CREST conformer sampling → energy-window + RMSD filtering → ORCA DFT opt+freq → ORCA high-level single point |
| **Thermodynamics** | Boltzmann population analysis, free-energy calculations, temperature sweeps for reaction feasibility |
| **Structure analysis** | Kabsch RMSD, H-bond detection, nucleophilic/electrophilic site identification, prebiotic interaction scoring |
| **Reports** | Plain-text ASCII reports, CSV population tables, H-bond tables |

---

## Installation

### Prerequisites

```bash
# Core scientific Python
pip install numpy

# Optional (required for SMILES → 3-D coordinate generation)
conda install -c conda-forge rdkit

# External computational chemistry tools (must be in PATH)
# xTB  — https://xtb-docs.readthedocs.io
# CREST — https://crest-lab.github.io/crest-docs/
# ORCA  — https://www.faccts.de/orca/ (free for academics)
```

### Install prebiotic_chem

```bash
# From the IAK_ASIf repository root:
pip install -e ./prebiotic_chem   # editable install
# or simply add the repo root to PYTHONPATH:
export PYTHONPATH=/path/to/IAK_ASIf:$PYTHONPATH
```

---

## Quick Start

### List available scenarios and molecules

```bash
python -m prebiotic_chem --list-scenarios
python -m prebiotic_chem --list-molecules
python -m prebiotic_chem --list-pathways
```

### Show detailed scenario information

```bash
python -m prebiotic_chem --scenario warm_little_pond --info
python -m prebiotic_chem --scenario hydrothermal_vent --info
```

### Convert SMILES to a 3-D XYZ structure

```bash
python -m prebiotic_chem --smiles "NCC(=O)O" --out-xyz glycine.xyz
```

### Run the full computational pipeline

```bash
python -m prebiotic_chem \
    --xyz glycine.xyz \
    --scenario warm_little_pond \
    --run-dir glycine_wlp \
    --cores 8 \
    --orca-opt-method "B97-3c Opt Freq" \
    --orca-sp-method "B3LYP D3BJ def2-TZVP TightSCF"
```

### Run only xTB + CREST (no ORCA)

```bash
python -m prebiotic_chem \
    --xyz adenine.xyz \
    --scenario miller_urey \
    --run-dir adenine_mu \
    --skip-orca-opt --skip-orca-sp
```

### Boltzmann population analysis from CSV

```bash
# CSV must have columns: id, energy_hartree
python -m prebiotic_chem \
    --boltzmann-csv conformer_energies.csv \
    --temperature 333
```

### Temperature sweep for a reaction

```bash
python -m prebiotic_chem \
    --temp-sweep \
    --delta-H -5.2 \
    --delta-S 12.0 \
    --scenario hydrothermal_vent
```

### Generate a comprehensive report

```bash
python -m prebiotic_chem \
    --scenario warm_little_pond \
    --molecule glycine \
    --pathway strecker_amino_acid \
    --run-dir glycine_wlp \
    --report glycine_report.txt
```

---

## Python API

```python
import prebiotic_chem as pc

# ── Scenarios ─────────────────────────────────────────────────────────
sc = pc.get_scenario("warm_little_pond")
print(sc.summary())
print(sc.xtb_flags())           # ['--gfn2', '--etemp', '333.15', '--alpb', 'water']
print(sc.orca_solvent_keyword()) # 'CPCM(Water)'

# ── Molecules ─────────────────────────────────────────────────────────
mol = pc.get_molecule("glycine")
print(mol.smiles)     # 'NCC(=O)O'
print(mol.relevance)

# All nucleobases
bases = pc.list_molecules("nucleobase")

# ── Reactions ─────────────────────────────────────────────────────────
pw = pc.get_pathway("strecker_amino_acid")
print(pw.summary())

relevant = pc.steps_for_scenario("warm_little_pond")

# ── Thermodynamics ────────────────────────────────────────────────────
import numpy as np
energies = np.array([-76.123, -76.118, -76.110])
pops = pc.boltzmann_populations(energies, ids=["conf1","conf2","conf3"], temperature_K=333)
for p in pops:
    print(f"{p['id']}: ΔE={p['rel_energy_kcal']:.2f} kcal/mol  pop={p['population_pct']:.1f}%")

dG, Keq = pc.reaction_free_energy(delta_H_kcal=-5.2, delta_S_cal_mol_K=12.0, temperature_K=333)
print(f"ΔG = {dG:.2f} kcal/mol, K_eq = {Keq:.2e}")

sweep = pc.temperature_sweep(-5.2, 12.0)

# ── Structure analysis ────────────────────────────────────────────────
syms, coords, _ = pc.read_xyz("glycine.xyz")
hbonds = pc.detect_hbonds(syms, coords)
score  = pc.score_prebiotic_geometry(syms, coords)
nuc    = pc.nucleophilic_sites(syms)

# ── SMILES → XYZ ─────────────────────────────────────────────────────
syms, coords = pc.smiles_to_xyz("NCC(=O)O", output_path="glycine_3d.xyz")

# ── Full pipeline ─────────────────────────────────────────────────────
pl = pc.PrebioticPipeline(
    run_dir="my_run",
    scenario=sc,
    input_xyz="glycine.xyz",
    n_cores=8,
    skip_orca_opt=True,   # xTB + CREST only for a quick test
    skip_orca_sp=True,
)
pl.run()

# ── Reports ───────────────────────────────────────────────────────────
report = pc.PrebioticReport(scenario=sc, molecule=mol, pathway=pw, run_dir="my_run")
print(report.generate())
report.write("glycine_report.txt")
```

---

## Pipeline Architecture

```
Input XYZ / SMILES
       │
       ▼
 ┌─────────────────┐
 │ Step 1: xTB     │  GFN2-xTB, ALPB implicit solvation,
 │ pre-optimisation│  temperature-dependent (--etemp T)
 └────────┬────────┘
          │
          ▼
 ┌─────────────────┐
 │ Step 2: CREST   │  GFN2-xTB conformer sampling,
 │ conformer search│  scenario temperature & solvent
 └────────┬────────┘
          │
          ▼
 ┌─────────────────┐
 │ Step 3: Filter  │  Energy window (≤ 5 kcal/mol default)
 │ & cluster       │  + Kabsch-RMSD clustering (0.5 Å default)
 └────────┬────────┘
          │
          ▼
 ┌─────────────────┐
 │ Step 4: ORCA    │  B97-3c Opt Freq (default),
 │ DFT opt + freq  │  CPCM(Water) implicit solvation
 └────────┬────────┘
          │  Boltzmann filter (99 % population)
          ▼
 ┌─────────────────┐
 │ Step 5: ORCA    │  B3LYP-D3BJ/def2-TZVP (default),
 │ high-level SP   │  CPCM solvation
 └────────┬────────┘
          │
          ▼
  Summary CSV + report
```

---

## Prebiotic Scenarios

| Key | Display name | T range (°C) | pH | Ionic strength |
|-----|-------------|-------------|-----|----------------|
| `warm_little_pond` | Warm Little Pond (Darwin's Pond) | 20–90 | 5–8 | 0.01–0.5 M |
| `hydrothermal_vent` | Submarine Black-Smoker Vent | 50–400 | 2–5 | 0.05–1.0 M |
| `alkaline_vent` | Alkaline Vent (Lost City) | 50–150 | 9–11 | 0.05–0.5 M |
| `ice_eutectic` | Ice-Eutectic / Frozen Lake | −30–0 | 4–7 | 0.5–5.0 M |
| `tidal_pool` | Tidal Pool / Littoral Zone | 15–80 | 6–9 | 0.1–1.0 M |
| `miller_urey` | Miller–Urey Reducing Atmosphere | 25–100 | 6–8 | 0–0.1 M |
| `volcanic_spring` | Volcanic Hot Spring | 60–150 | 2–6 | 0.1–2.0 M |

---

## Molecule Categories

- `feedstock` — HCN, formaldehyde, NH₃, CO₂, CH₄, H₂S, cyanamide, phosphoric acid
- `solvent` — water
- `amino_acid` — glycine, alanine, serine, aspartic acid, glutamic acid
- `nucleobase` — adenine, guanine, cytosine, uracil, thymine, hypoxanthine
- `sugar` — ribose, deoxyribose, glycolaldehyde, glyceraldehyde
- `nucleotide` — AMP, ATP
- `lipid` — decanoic acid, oleic acid
- `cofactor` — nicotinamide
- `metabolite` — acetic acid, pyruvic acid
- `mineral_proxy` — silicic acid, boric acid

---

## Reaction Pathways

| Key | Pathway | Steps |
|-----|---------|-------|
| `strecker_amino_acid` | Strecker Amino Acid Synthesis | Imine → aminonitrile → amino acid |
| `formose_sugars` | Formose Reaction | Glycolaldehyde → glyceraldehyde → ribose |
| `hcn_to_adenine` | HCN → Adenine | DAMN → pentamer → adenine |
| `nucleotide_synthesis` | Nucleotide Assembly | N-glycosidic bond + phosphorylation |
| `iron_sulfur_acetate` | Iron–Sulfur World | CO₂ → acetate (FeS catalysis) |
| `protocell_assembly` | Protocell Vesicle Formation | Fatty acid self-assembly |

---

## References

- Miller, S.L. (1953). *Science* **117**, 528.
- Oró, J. (1960). *Biochem. Biophys. Res. Commun.* **2**, 407.
- Butlerow, A. (1861). *Liebigs Ann.* **120**, 295.
- Wächtershäuser, G. (1988). *Microbiol. Rev.* **52**, 452.
- Ricardo, A. et al. (2004). *Science* **303**, 196.
- Deamer, D.W. (1985). *Nature* **317**, 792.
- Szostak, J.W. et al. (2001). *Nature* **409**, 387.
- Russell, M.J. & Martin, W. (2004). *Trends Biochem. Sci.* **29**, 358.
- Becker, S. et al. (2019). *Science* **366**, 76.
- Migliaro, I. et al. (2025). *ChemRxiv* 10.26434/chemrxiv-2025-cvg1x (**ChemRefine**).

---

## License

This package is part of the IAK_ASIf repository. See the repository root for license information.
