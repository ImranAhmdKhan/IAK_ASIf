# IAK — Intelligent Automated Kluster-generator Pipeline

IAK is a graphical desktop application for automated generation and refinement of molecular clusters used in non-covalent interaction studies and related computational chemistry workflows. It combines random cluster generation, geometry scoring, RMSD-based clustering, and optional high-level quantum-chemistry refinement using **xTB**, **CREST**, and **ORCA**.

---

## Table of Contents

1. [Features](#features)
2. [Requirements](#requirements)
3. [Installation](#installation)
4. [Quick Start](#quick-start)
5. [Usage Guide](#usage-guide)
   - [Launching the GUI](#launching-the-gui)
   - [Input Files](#input-files)
   - [Run Modes](#run-modes)
   - [Advanced Settings](#advanced-settings)
   - [Output Structure](#output-structure)
6. [External Engines](#external-engines)
7. [Programmatic (Headless) Use](#programmatic-headless-use)
8. [Troubleshooting](#troubleshooting)
9. [License](#license)

---

## Features

- **Cluster generation** — randomly place one or more guest molecules around a host (anchor) molecule with clash detection and user-controlled stoichiometry (e.g. 1 host : 3 guests).
- **Geometry scoring** — score candidates by inter-fragment distances and contact quality.
- **RMSD clustering** — keep only structurally unique clusters (Kabsch RMSD).
- **xTB pre-optimisation** — fast semi-empirical GFN2-xTB geometry optimisations.
- **CREST conformational search** — metadynamics-based conformer/rotamer ensemble generation.
- **ORCA refinement** — DFT (default: B97-3c) optimisation and frequency calculation.
- **Tkinter GUI** — no command-line experience required; progress bars, live log, energy plots, and 3-D structure viewer.
- **Cross-platform** — Linux, macOS, and Windows (xTB/CREST run under WSL on Windows).
- **Resumable jobs** — state is saved after each stage; restarting the pipeline picks up where it left off.

---

## Requirements

| Dependency | Version | Notes |
|---|---|---|
| Python | ≥ 3.9 | Tested on 3.10 and 3.11 |
| numpy | any recent | cluster geometry math |
| matplotlib | any recent | energy plots |
| tkinter | bundled with Python | GUI; may need `python3-tk` on Linux |

**Optional external engines** (see [External Engines](#external-engines)):

| Engine | Purpose |
|---|---|
| [xTB](https://github.com/grimme-lab/xtb) | Semi-empirical pre-optimisation |
| [CREST](https://github.com/grimme-lab/crest) | Conformational search |
| [ORCA](https://www.faccts.de/orca/) | DFT refinement |

IAK can download and unpack xTB and CREST automatically (Linux/WSL only). ORCA must be installed manually.

---

## Installation

### 1. Clone the repository

```bash
git clone https://github.com/ImranAhmdKhan/IAK_ASIf.git
cd IAK_ASIf
```

### 2. Create and activate a virtual environment (recommended)

```bash
python3 -m venv .venv
source .venv/bin/activate      # Linux / macOS
# .venv\Scripts\activate       # Windows PowerShell
```

### 3. Install Python dependencies

```bash
pip install numpy matplotlib
```

On Ubuntu/Debian, also install tkinter if it is missing:

```bash
sudo apt install python3-tk
```

---

## Quick Start

```bash
# Run the GUI
python IAK_ASIf.py

# Or run as a module
python -m iak
```

The main window opens. Load your host and guest XYZ files, set the stoichiometry, choose a run mode, and click **Run**.

---

## Usage Guide

### Launching the GUI

```bash
python IAK_ASIf.py          # from the repository root
# OR
python -m iak               # equivalent
```

### Input Files

All molecular inputs must be in **XYZ format** (standard `.xyz`):

```
<number of atoms>
<comment line (may be blank)>
<Element>  <x>  <y>  <z>
...
```

- **Anchor (host) molecule** — the central molecule that guest molecules are placed around.
- **Guest molecule(s)** — one or more molecules to be randomly distributed around the anchor. Multiple guest types are supported simultaneously.

Use the **Browse** buttons in the GUI to select the files.

### Stoichiometry

Set how many copies of each molecule appear in the cluster:

- **Anchor copies** — typically `1`; increase for host–host clusters.
- **Guest copies** — number of guest molecules per cluster (e.g. `3` for a 1:3 complex).

### Run Modes

| Mode | Clusters generated | xTB runs | CREST runs | Recommended for |
|---|---|---|---|---|
| **Fast** | 50 | 5 | 1 | Quick exploration |
| **Balanced** *(default)* | 200 | 20 | 5 | General use |
| **Thorough** | 1 000 | 50 | 10 | Publication-quality results |

### Advanced Settings

Click **Advanced / Settings** to adjust:

| Setting | Default | Description |
|---|---|---|
| `n_generate` | 200 | Raw clusters to generate |
| `n_keep_scored` | 50 | Top-scored clusters passed to RMSD filter |
| `n_keep_clustered` | 40 | Unique clusters after RMSD pruning |
| `n_run_xtb` | 20 | Clusters submitted to xTB |
| `n_run_crest` | 5 | xTB minima submitted to CREST |
| `rmsd_cutoff` | 0.5 Å | Minimum RMSD for two structures to be considered distinct |
| `clash_cutoff` | 1.2 Å | Minimum allowed inter-atomic distance during generation |
| `xtb_method` | `--gfn2` | xTB Hamiltonian |
| `crest_method` | `--gfn2` | CREST Hamiltonian |
| `orca_method` | `B97-3c Opt Freq` | ORCA keywords |
| `charge` | 0 | Total system charge |
| `multiplicity` | 1 | Spin multiplicity |
| `cores` | 4 | CPU threads for xTB/ORCA |
| `maxcore` | 2000 | Memory per core for ORCA (MB) |
| `random_seed` | 42 | Reproducibility seed |
| `reaction_type` | Non-covalent | Affects scoring weights |

### Reaction Types

Select from the drop-down:

- Non-covalent *(default)*
- Covalent
- Substitution
- Addition
- Nucleophilic
- Electrophilic

### Output Structure

All results are written to the **output directory** you specify (default: `IAK_output/` next to the input files):

```
IAK_output/
├── 01_Inputs_and_Clusters/   # Raw input copies + generated clusters (.xyz)
├── 02_xTB_Results/           # xTB-optimised geometries and energies
├── 03_CREST_Results/         # CREST conformer ensembles
├── 04_ORCA_Refinement/       # ORCA optimised structures, frequencies, energies
├── 05_Top_Models_Comparison/ # Final ranked structures + energy summary
├── logs/                     # Timestamped log files
├── state.json                # Pipeline checkpoint (enables resume)
└── provenance_manifest.json  # Full record of every calculation
```

The GUI displays an **energy ranking table** and **bar chart** of the top structures when the run finishes.

---

## External Engines

### xTB

IAK attempts to download a pre-compiled xTB binary automatically on first use (Linux/WSL). Alternatively, install manually:

```bash
# Ubuntu/Debian via conda-forge
conda install -c conda-forge xtb
# or download from https://github.com/grimme-lab/xtb/releases
```

Make sure `xtb` is on your `PATH`, or place the binary in `iak/iak_engine/xtb/bin/xtb`.

### CREST

IAK attempts to download CREST automatically (Linux/WSL). Manual installation:

```bash
conda install -c conda-forge crest
# or https://github.com/grimme-lab/crest/releases
```

Make sure `crest` is on your `PATH`, or place the binary in `iak/iak_engine/crest/crest`.

### ORCA

ORCA is **not** auto-downloaded (licence required). Download from:
<https://www.faccts.de/orca/>

After installation, either:
- Add the ORCA directory to your `PATH`, **or**
- Place the ORCA directory inside `iak/iak_engine/` so IAK detects it automatically.

On **Windows**, IAK can call an ORCA installation inside WSL.

---

## Programmatic (Headless) Use

You can run the pipeline without the GUI:

```python
from iak.models import Config, RunMode
from iak.pipeline import Pipeline

config = Config.from_mode(RunMode.BALANCED)
config.charge = 0
config.multiplicity = 1

pipe = Pipeline(
    fragA_path="host.xyz",
    fragB_paths=["guest.xyz"],
    n_guests_list=[3],          # 1 host + 3 guests
    config=config,
    out_dir="my_output",
    reaction_type="Non-covalent",
    n_anchor=1,
)

pipe.run()                       # blocking; prints progress to stdout
```

Pass a `progress_cb` callable `(percent: float, message: str) -> None` to receive live progress updates.

---

## Troubleshooting

| Problem | Solution |
|---|---|
| `ModuleNotFoundError: tkinter` | `sudo apt install python3-tk` (Linux) |
| xTB / CREST not found | Install binaries and add to `PATH`, or place in `iak/iak_engine/` |
| ORCA not found on Windows | Ensure ORCA is installed in WSL and accessible as `orca` |
| Job hangs at CREST stage | CREST may not be compatible with your kernel; update or skip CREST in settings |
| Output directory is full | Check disk space; each CREST run can produce large trajectory files |
| Resume does not work | Delete `state.json` in the output directory to force a clean restart |

---

## License

See [LICENSE](LICENSE) for details.
