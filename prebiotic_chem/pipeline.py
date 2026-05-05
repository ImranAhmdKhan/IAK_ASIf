"""
prebiotic_chem.pipeline
=======================
Multi-step computational workflow orchestrator for prebiotic chemistry.

Implements the canonical ChemRefine-inspired pipeline adapted for prebiotic
conditions:

    Step 1 — xTB pre-optimisation (fast, implicit solvation)
    Step 2 — CREST conformer sampling (GFN2-xTB, ALPB water)
    Step 3 — Energy-window / Boltzmann filtering
    Step 4 — ORCA geometry optimisation + frequency (DFT, with CPCM)
    Step 5 — ORCA single-point energy at higher level of theory

Each step writes its output to ``<run_dir>/<step_name>/`` and produces a
plain CSV summary.  The pipeline is resumable: if output already exists it
skips that step.

Usage example::

    from prebiotic_chem.conditions import get_scenario
    from prebiotic_chem.pipeline import PrebioticPipeline

    scenario = get_scenario("warm_little_pond")
    pl = PrebioticPipeline(
        run_dir="my_run",
        scenario=scenario,
        input_xyz="molecule.xyz",
    )
    pl.run()
"""

from __future__ import annotations

import csv
import logging
import os
import shutil
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from .analysis import rmsd_cluster, score_prebiotic_geometry, molecular_formula
from .conditions import PrebioticScenario
from .constants import DEFAULT_TEMPERATURE, HARTREE_TO_KCAL
from .io_utils import (
    read_xyz,
    read_xyz_trajectory,
    run_crest,
    run_xtb,
    write_orca_input,
    write_xyz,
    parse_orca_energy,
    parse_orca_thermochemistry,
)
from .thermodynamics import (
    boltzmann_populations,
    filter_energy_window,
    filter_by_boltzmann_cutoff,
)

logger = logging.getLogger(__name__)


class PrebioticPipeline:
    """
    Orchestrates the full prebiotic computational chemistry workflow.

    Parameters
    ----------
    run_dir:
        Root directory for all pipeline outputs.
    scenario:
        :class:`~prebiotic_chem.conditions.PrebioticScenario` defining the
        environmental conditions.
    input_xyz:
        Path to the starting XYZ structure.
    charge:
        Total molecular charge.
    multiplicity:
        Spin multiplicity.
    n_cores:
        CPU cores to use for xTB, CREST, and ORCA.
    maxcore_mb:
        Memory per ORCA core in MB.
    orca_opt_method:
        ORCA keyword line for geometry optimisation.
    orca_sp_method:
        ORCA keyword line for the high-level single-point step.
    energy_window_kcal:
        Energy window (kcal mol⁻¹) for conformer filtering after CREST.
    rmsd_cutoff:
        RMSD cutoff (Å) for clustering before ORCA.
    max_orca_structures:
        Maximum number of structures to submit to the ORCA steps.
    boltzmann_cutoff_pct:
        Cumulative Boltzmann population cutoff for final ensemble selection.
    skip_xtb:
        Skip the xTB pre-optimisation step.
    skip_crest:
        Skip the CREST conformer sampling step.
    skip_orca_opt:
        Skip the ORCA geometry optimisation step.
    skip_orca_sp:
        Skip the ORCA single-point step.
    """

    def __init__(
        self,
        run_dir: str | Path,
        scenario: PrebioticScenario,
        input_xyz: str | Path,
        charge: int = 0,
        multiplicity: int = 1,
        n_cores: int = 4,
        maxcore_mb: int = 2000,
        orca_opt_method: str = "B97-3c Opt Freq",
        orca_sp_method: str = "B3LYP D3BJ def2-TZVP TightSCF",
        energy_window_kcal: float = 5.0,
        rmsd_cutoff: float = 0.5,
        max_orca_structures: int = 10,
        boltzmann_cutoff_pct: float = 99.0,
        skip_xtb: bool = False,
        skip_crest: bool = False,
        skip_orca_opt: bool = False,
        skip_orca_sp: bool = False,
    ) -> None:
        self.run_dir = Path(run_dir)
        self.scenario = scenario
        self.input_xyz = Path(input_xyz)
        self.charge = charge
        self.multiplicity = multiplicity
        self.n_cores = n_cores
        self.maxcore_mb = maxcore_mb
        self.orca_opt_method = orca_opt_method
        self.orca_sp_method = orca_sp_method
        self.energy_window_kcal = energy_window_kcal
        self.rmsd_cutoff = rmsd_cutoff
        self.max_orca_structures = max_orca_structures
        self.boltzmann_cutoff_pct = boltzmann_cutoff_pct
        self.skip_xtb = skip_xtb
        self.skip_crest = skip_crest
        self.skip_orca_opt = skip_orca_opt
        self.skip_orca_sp = skip_orca_sp

        self.run_dir.mkdir(parents=True, exist_ok=True)
        self._setup_logging()

    # ------------------------------------------------------------------ #
    # Logging                                                               #
    # ------------------------------------------------------------------ #

    def _setup_logging(self) -> None:
        log_file = self.run_dir / "pipeline.log"
        fh = logging.FileHandler(log_file, mode="a")
        fh.setLevel(logging.DEBUG)
        fh.setFormatter(logging.Formatter("%(asctime)s %(levelname)-8s %(message)s"))
        logger.addHandler(fh)
        if not logger.handlers or not any(
            isinstance(h, logging.StreamHandler) for h in logger.handlers
        ):
            sh = logging.StreamHandler()
            sh.setLevel(logging.INFO)
            sh.setFormatter(logging.Formatter("%(levelname)-8s %(message)s"))
            logger.addHandler(sh)
        logger.setLevel(logging.DEBUG)

    # ------------------------------------------------------------------ #
    # Top-level run                                                         #
    # ------------------------------------------------------------------ #

    def run(self) -> Dict:
        """
        Execute the full pipeline and return a summary dict.
        """
        logger.info("=" * 65)
        logger.info("Prebiotic Chemistry Pipeline")
        logger.info("Scenario : %s", self.scenario.display_name)
        logger.info("T = %.1f °C | pH %.1f | I = %.3f M",
                    self.scenario.temperature_C,
                    self.scenario.pH,
                    self.scenario.ionic_strength_M)
        logger.info("Input XYZ: %s", self.input_xyz)
        logger.info("=" * 65)

        # --- Step 1: xTB pre-optimisation ---
        xtb_xyz = self._step_xtb_preopt()

        # --- Step 2: CREST conformer sampling ---
        conformers = self._step_crest(xtb_xyz)

        # --- Step 3: Filtering ---
        filtered = self._step_filter(conformers)

        # --- Step 4: ORCA geometry optimisation ---
        orca_opt_results = self._step_orca_opt(filtered)

        # --- Step 5: ORCA single-point ---
        orca_sp_results  = self._step_orca_sp(orca_opt_results)

        # --- Summary ---
        summary = self._write_summary(orca_sp_results or orca_opt_results)
        logger.info("Pipeline complete.  Results in: %s", self.run_dir)
        return summary

    # ------------------------------------------------------------------ #
    # Step 1: xTB pre-optimisation                                          #
    # ------------------------------------------------------------------ #

    def _step_xtb_preopt(self) -> Path:
        step_dir = self.run_dir / "01_xtb_preopt"
        step_dir.mkdir(exist_ok=True)
        result_xyz = step_dir / "xtbopt.xyz"

        if self.skip_xtb:
            logger.info("[Step 1] xTB pre-opt skipped (--skip-xtb).")
            return self.input_xyz

        if result_xyz.exists():
            logger.info("[Step 1] xTB output found — skipping.")
            return result_xyz

        # Copy input to step dir
        shutil.copy2(self.input_xyz, step_dir / self.input_xyz.name)
        flags = self.scenario.xtb_flags() + ["--opt"]

        logger.info("[Step 1] Running xTB pre-optimisation ...")
        res = run_xtb(
            xyz_path=step_dir / self.input_xyz.name,
            flags=flags,
            work_dir=step_dir,
            n_cores=self.n_cores,
        )
        if res["returncode"] != 0 or not res["converged"]:
            logger.warning("[Step 1] xTB did not converge — using input geometry.")
            return self.input_xyz

        energy = res["energy_hartree"]
        logger.info(
            "[Step 1] xTB optimisation converged.  E = %.8f Eh", energy or 0.0
        )
        return res["opt_xyz"] or self.input_xyz

    # ------------------------------------------------------------------ #
    # Step 2: CREST conformer sampling                                      #
    # ------------------------------------------------------------------ #

    def _step_crest(
        self, input_xyz: Path
    ) -> List[Tuple[List[str], np.ndarray, Optional[float]]]:
        step_dir = self.run_dir / "02_crest"
        step_dir.mkdir(exist_ok=True)
        ensemble_path = step_dir / "crest_conformers.xyz"

        if self.skip_crest:
            logger.info("[Step 2] CREST skipped (--skip-crest).")
            syms, coords, _ = read_xyz(input_xyz)
            return [(syms, coords, None)]

        if ensemble_path.exists():
            logger.info("[Step 2] CREST ensemble found — skipping.")
            frames = read_xyz_trajectory(ensemble_path)
            return [(s, c, self._energy_from_comment(comment))
                    for s, c, comment in frames]

        shutil.copy2(input_xyz, step_dir / input_xyz.name)
        flags = self.scenario.crest_flags()
        logger.info("[Step 2] Running CREST conformer sampling ...")
        res = run_crest(
            xyz_path=step_dir / input_xyz.name,
            flags=flags,
            work_dir=step_dir,
            n_cores=self.n_cores,
        )
        if res["returncode"] != 0:
            logger.warning("[Step 2] CREST failed — using single input geometry.")
            syms, coords, _ = read_xyz(input_xyz)
            return [(syms, coords, None)]

        logger.info("[Step 2] CREST found %d conformers.", res["n_conformers"])
        return res["conformers"] if res["conformers"] else [(read_xyz(input_xyz)[:2] + (None,))]

    # ------------------------------------------------------------------ #
    # Step 3: Filtering                                                     #
    # ------------------------------------------------------------------ #

    def _step_filter(
        self, conformers: List[Tuple[List[str], np.ndarray, Optional[float]]]
    ) -> List[Tuple[List[str], np.ndarray, Optional[float]]]:
        step_dir = self.run_dir / "03_filtered"
        step_dir.mkdir(exist_ok=True)

        logger.info("[Step 3] Filtering %d conformers ...", len(conformers))

        # Separate those with and without energies
        with_e  = [(s, c, e) for s, c, e in conformers if e is not None]
        no_e    = [(s, c, e) for s, c, e in conformers if e is None]

        if with_e:
            energies = [e for _, _, e in with_e]
            kept_idx = filter_energy_window(energies, window_kcal=self.energy_window_kcal)
            with_e = [with_e[i] for i in kept_idx]
            logger.info("[Step 3] Energy window (%.1f kcal/mol): %d kept",
                        self.energy_window_kcal, len(with_e))

            # RMSD clustering
            coords_list = [c for _, c, _ in with_e]
            kept_cluster = rmsd_cluster(coords_list, cutoff=self.rmsd_cutoff)
            with_e = [with_e[i] for i in kept_cluster]
            logger.info("[Step 3] After RMSD clustering (%.2f Å): %d unique",
                        self.rmsd_cutoff, len(with_e))

        filtered = (with_e + no_e)[: self.max_orca_structures]
        logger.info("[Step 3] Sending %d structures to ORCA.", len(filtered))

        # Write filtered ensemble
        for idx, (syms, coords, _e) in enumerate(filtered):
            write_xyz(
                step_dir / f"conf_{idx + 1:04d}.xyz",
                syms, coords,
                comment=f"Filtered conformer {idx + 1}",
            )
        return filtered

    # ------------------------------------------------------------------ #
    # Step 4: ORCA geometry optimisation                                    #
    # ------------------------------------------------------------------ #

    def _step_orca_opt(
        self, conformers: List[Tuple[List[str], np.ndarray, Optional[float]]]
    ) -> List[Dict]:
        step_dir = self.run_dir / "04_orca_opt"
        step_dir.mkdir(exist_ok=True)

        if self.skip_orca_opt:
            logger.info("[Step 4] ORCA opt skipped (--skip-orca-opt).")
            return []

        solvent_kw = self.scenario.orca_solvent_keyword()
        results = []
        for idx, (syms, coords, _e) in enumerate(conformers, 1):
            label = f"conf_{idx:04d}"
            inp_file = step_dir / f"{label}.inp"
            out_file = step_dir / f"{label}.out"

            if out_file.exists():
                logger.info("[Step 4] %s — output found, skipping.", label)
                e = parse_orca_energy(out_file)
                thermo = parse_orca_thermochemistry(out_file)
                results.append({
                    "label": label, "symbols": syms,
                    "energy_hartree": e, "thermochemistry": thermo,
                    "out_file": out_file,
                })
                continue

            write_orca_input(
                inp_file, syms, coords,
                method=self.orca_opt_method,
                charge=self.charge,
                multiplicity=self.multiplicity,
                n_cores=self.n_cores,
                maxcore_mb=self.maxcore_mb,
                solvent_keyword=solvent_kw,
                comment=f"Prebiotic scenario: {self.scenario.display_name}",
            )
            logger.info("[Step 4] ORCA opt: %s", label)
            ok, e = self._run_orca(inp_file, out_file)
            if not ok:
                logger.warning("[Step 4] %s — ORCA failed.", label)
                continue
            thermo = parse_orca_thermochemistry(out_file)
            results.append({
                "label": label, "symbols": syms,
                "energy_hartree": e, "thermochemistry": thermo,
                "out_file": out_file,
            })

        logger.info("[Step 4] ORCA opt complete: %d structures.", len(results))
        self._write_step_csv(step_dir / "results.csv", results)
        return results

    # ------------------------------------------------------------------ #
    # Step 5: ORCA single-point                                             #
    # ------------------------------------------------------------------ #

    def _step_orca_sp(self, opt_results: List[Dict]) -> List[Dict]:
        step_dir = self.run_dir / "05_orca_sp"
        step_dir.mkdir(exist_ok=True)

        if self.skip_orca_sp or not opt_results:
            logger.info("[Step 5] ORCA SP skipped.")
            return []

        # Apply Boltzmann filtering from opt energies
        energies = [r["energy_hartree"] for r in opt_results
                    if r["energy_hartree"] is not None]
        if not energies:
            logger.warning("[Step 5] No valid opt energies — skipping SP.")
            return []

        selected_idx = filter_by_boltzmann_cutoff(
            energies,
            temperature_K=self.scenario.temperature_K,
            cumulative_pct=self.boltzmann_cutoff_pct,
        )
        selected = [opt_results[i] for i in selected_idx]
        logger.info("[Step 5] Boltzmann filter (%.0f %%, T=%.1f K): %d → %d structures",
                    self.boltzmann_cutoff_pct, self.scenario.temperature_K,
                    len(opt_results), len(selected))

        solvent_kw = self.scenario.orca_solvent_keyword()
        sp_results = []
        for res in selected:
            label = f"{res['label']}_sp"
            inp_file = step_dir / f"{label}.inp"
            out_file = step_dir / f"{label}.out"

            # Read optimised geometry from step 4 (parse from .out or use stored coords)
            syms = res["symbols"]
            coords_from_opt = self._read_orca_opt_coords(res["out_file"]) or \
                              np.zeros((len(syms), 3))

            if out_file.exists():
                logger.info("[Step 5] %s — output found, skipping.", label)
                e = parse_orca_energy(out_file)
                sp_results.append({"label": label, "symbols": syms,
                                   "energy_hartree": e, "out_file": out_file})
                continue

            write_orca_input(
                inp_file, syms, coords_from_opt,
                method=self.orca_sp_method,
                charge=self.charge,
                multiplicity=self.multiplicity,
                n_cores=self.n_cores,
                maxcore_mb=self.maxcore_mb,
                solvent_keyword=solvent_kw,
                comment=f"High-level SP — {self.scenario.display_name}",
            )
            logger.info("[Step 5] ORCA SP: %s", label)
            ok, e = self._run_orca(inp_file, out_file)
            if not ok:
                logger.warning("[Step 5] %s — ORCA SP failed.", label)
                continue
            sp_results.append({"label": label, "symbols": syms,
                                "energy_hartree": e, "out_file": out_file})

        logger.info("[Step 5] ORCA SP complete: %d structures.", len(sp_results))
        self._write_step_csv(step_dir / "results.csv", sp_results)
        return sp_results

    # ------------------------------------------------------------------ #
    # Summary                                                               #
    # ------------------------------------------------------------------ #

    def _write_summary(self, final_results: List[Dict]) -> Dict:
        summary_path = self.run_dir / "summary.csv"

        if not final_results:
            return {"status": "no_results"}

        energies = [r["energy_hartree"] for r in final_results
                    if r.get("energy_hartree") is not None]
        ids      = [r["label"] for r in final_results
                    if r.get("energy_hartree") is not None]

        if not energies:
            return {"status": "no_energies"}

        pops = boltzmann_populations(
            energies, ids, temperature_K=self.scenario.temperature_K
        )

        fieldnames = ["rank", "id", "energy_hartree",
                      "rel_energy_kcal", "population_pct"]
        with open(summary_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            for rank, row in enumerate(pops, 1):
                writer.writerow({
                    "rank":             rank,
                    "id":               row["id"],
                    "energy_hartree":   f"{row['energy_hartree']:.8f}",
                    "rel_energy_kcal":  f"{row['rel_energy_kcal']:.3f}",
                    "population_pct":   f"{row['population_pct']:.2f}",
                })

        logger.info("[Summary] Written to %s", summary_path)
        for rank, row in enumerate(pops[:5], 1):
            logger.info(
                "  #%d  %-20s  ΔE = %6.2f kcal/mol  pop = %5.1f %%",
                rank, row["id"], row["rel_energy_kcal"], row["population_pct"],
            )
        return {"status": "ok", "populations": pops, "summary_csv": str(summary_path)}

    # ------------------------------------------------------------------ #
    # Private helpers                                                        #
    # ------------------------------------------------------------------ #

    @staticmethod
    def _run_orca(inp_file: Path, out_file: Path) -> Tuple[bool, Optional[float]]:
        """Call the ORCA binary and parse the output energy."""
        import subprocess
        result = subprocess.run(
            ["orca", str(inp_file)],
            capture_output=True, text=True,
            cwd=str(inp_file.parent),
        )
        out_text = result.stdout + result.stderr
        out_file.write_text(out_text)
        e = parse_orca_energy(out_file)
        ok = result.returncode == 0 and e is not None
        return ok, e

    @staticmethod
    def _read_orca_opt_coords(
        orca_out: Path,
    ) -> Optional[np.ndarray]:
        """Extract the final Cartesian coordinates from an ORCA .out file."""
        import re as _re
        text = orca_out.read_text(errors="replace")
        blocks = _re.findall(
            r"CARTESIAN COORDINATES \(ANGSTROEM\)\s*\n-+\n((?:.*?\n)+?)-+",
            text, _re.DOTALL,
        )
        if not blocks:
            return None
        last_block = blocks[-1]
        coords = []
        for line in last_block.splitlines():
            parts = line.split()
            if len(parts) == 4:
                try:
                    coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
                except ValueError:
                    pass
        return np.array(coords, dtype=float) if coords else None

    @staticmethod
    def _energy_from_comment(comment: str) -> Optional[float]:
        import re as _re
        m = _re.search(r"[-]?\d+\.\d+", comment)
        return float(m.group()) if m else None

    @staticmethod
    def _write_step_csv(path: Path, results: List[Dict]) -> None:
        if not results:
            return
        fieldnames = ["label", "energy_hartree"]
        with open(path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
            w.writeheader()
            for r in results:
                w.writerow({
                    "label": r.get("label", ""),
                    "energy_hartree": (
                        f"{r['energy_hartree']:.8f}"
                        if r.get("energy_hartree") is not None else "N/A"
                    ),
                })
