from __future__ import annotations

import json
import logging
import math
import os
import random
import shutil
import subprocess
import sys
import time
import traceback
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Tuple

from . import engines as _engines
from .engines import is_tool_available, get_wsl_path
from .constants import normalize_reaction_type, EH2KCAL, ORCA_CMD
from .models import Molecule, Config
from .models import read_multi_xyz
from .chemistry import (
    generate_cluster,
    score_geometry,
    kabsch_rmsd,
    atom_counts_for_molecule,
    atom_counts_to_formula,
    compare_atom_counts,
    expected_cluster_atom_counts,
)
from .utils import _now
from .job_summary import JobSummaryRecorder


class Pipeline:
    def __init__(self, fragA_path, fragB_paths, n_guests_list, config, out_dir, ratio_label=None, progress_cb=None, reaction_type="Non-covalent", n_anchor=1):
        self.fragA_path = fragA_path
        self.fragB_paths = fragB_paths if fragB_paths else []
        self.n_anchor = max(0, int(n_anchor))
        self.fragA = Molecule.from_xyz(fragA_path) if fragA_path and os.path.exists(fragA_path) else Molecule([])
        self.fragBs = [Molecule.from_xyz(p) if p and os.path.exists(p) else Molecule([]) for p in self.fragB_paths]
        self.n_guests_list = [max(0, int(x)) for x in (n_guests_list if isinstance(n_guests_list, list) else [n_guests_list])]
        self.reaction_type = normalize_reaction_type(reaction_type)
        self.config = config
        self.progress_cb = progress_cb
        self.progress_started = time.time()
        self.progress_percent = 0.0
        self.expected_atom_counts = expected_cluster_atom_counts(self.fragA, self.fragBs, self.n_anchor, self.n_guests_list)

        self.out_dir = os.path.abspath(out_dir)
        self.ratio_str = ratio_label if ratio_label else f"{self.n_anchor}_{'_'.join(map(str, self.n_guests_list))}"
        self.dirs = {
            "inputs": os.path.join(self.out_dir, "01_Inputs_and_Clusters"),
            "xtb": os.path.join(self.out_dir, "02_xTB_Results"),
            "crest": os.path.join(self.out_dir, "03_CREST_Results"),
            "orca": os.path.join(self.out_dir, "04_ORCA_Refinement"),
            "comparison": os.path.join(self.out_dir, "05_Top_Models_Comparison"),
            "logs": os.path.join(self.out_dir, "logs"),
        }
        for d in self.dirs.values():
            os.makedirs(d, exist_ok=True)
        if fragA_path and os.path.exists(fragA_path):
            shutil.copy2(fragA_path, os.path.join(self.dirs["inputs"], "raw_input_anchor.xyz"))
        for i, bp in enumerate(self.fragB_paths):
            if bp and os.path.exists(bp):
                shutil.copy2(bp, os.path.join(self.dirs["inputs"], f"raw_input_guest_{i}.xyz"))
        self.state_file = os.path.join(self.out_dir, "state.json")
        self.prov_file = os.path.join(self.out_dir, "provenance_manifest.json")
        self.state = json.load(open(self.state_file, "r", encoding="utf-8")) if os.path.exists(self.state_file) else {"gen": [], "filt": [], "xtb": {}, "crest": {}, "orca": {}}
        # Prune cached paths whose files have since been deleted so stages re-run cleanly.
        self.state["gen"] = [p for p in self.state.get("gen", []) if os.path.exists(p)]
        self.state["filt"] = [p for p in self.state.get("filt", []) if os.path.exists(p)]
        self.provenance = json.load(open(self.prov_file, "r", encoding="utf-8")) if os.path.exists(self.prov_file) else []
        self.job_summary = JobSummaryRecorder(self)

    def _validate_generated_cluster(self, mol: Molecule) -> Tuple[bool, str]:
        observed = atom_counts_for_molecule(mol)
        return compare_atom_counts(self.expected_atom_counts, observed)

    def save(self):
        with open(self.state_file, "w", encoding="utf-8") as f:
            json.dump(self.state, f, indent=4)
        with open(self.prov_file, "w", encoding="utf-8") as f:
            json.dump(self.provenance, f, indent=4)

    def _progress(self, stage, percent=None, status="running", message=""):
        if percent is not None:
            self.progress_percent = max(0.0, min(100.0, float(percent)))
        payload = {
            "job": self.ratio_str.replace("_", ":"),
            "stage": stage,
            "percent": self.progress_percent,
            "status": status,
            "message": message,
            "elapsed_seconds": time.time() - self.progress_started,
        }
        if self.progress_cb:
            try:
                self.progress_cb(payload)
            except Exception:
                pass

    def _record_engine_result(self, engine, source, flags, result, elapsed):
        self.job_summary.add_engine_result(engine, source, flags, result, elapsed)
        self.job_summary.write()

    def run(self, run_xtb=False, run_crest=False, run_orca=False, log_cb=None, status_cb=None):
        def _log(m):
            logging.getLogger("IAK").info(m)
            if log_cb:
                log_cb(m + "\n")

        self.progress_started = time.time()
        self.job_summary.capture_inputs(run_xtb, run_crest, run_orca)
        self.job_summary.log_event("Validation", "Captured job inputs, settings, and engine availability.")
        self.job_summary.write()
        self._progress("Validation", 3.0, "running", "Checking input correctness and engine readiness.")

        try:
            if run_orca and not is_tool_available("orca"):
                raise RuntimeError("Run ORCA DFT is checked, but ORCA is not installed or detected!\nPlease install it using the LOAD LOCAL ENGINE button.")

            if _engines.ORCA_DIR and " " in os.path.abspath(_engines.ORCA_DIR):
                _log("\n[CRITICAL WARNING] Your ORCA installation folder is inside a path with SPACES. IAK will self-heal and fall back if MPI fails.\n")

            _log("\nSCIENTIFIC NOTICE: Structures produced by this workflow are the lowest found under the current sampling settings. They are NOT guaranteed to be global minima.")
            _log(f"[INFO] Reaction mechanism: {self.reaction_type}")
            _log(f"[INFO] Expected atom composition: {atom_counts_to_formula(self.expected_atom_counts)}")

            if self.n_anchor > 0 and self.fragA.n_atoms() == 0:
                raise RuntimeError("Anchor stoichiometry is > 0 but Anchor A XYZ is missing or empty.")

            if sum(self.n_guests_list) == 0 and self.n_anchor <= 1:
                self._progress("Sampling/filtering", 12.0, "running", "Single Molecule Mode detected.")
                _log("\n[INFO] Single Molecule Mode detected. Bypassing cluster generation phase...")
                anchor_path = os.path.join(self.dirs["inputs"], "raw_input_anchor.xyz")
                if not self.state["gen"]:
                    self.state["gen"] = [anchor_path]
                if not self.state["filt"]:
                    self.state["filt"] = [anchor_path]
                    self.save()
                self.state["preopt_done"] = True
                self.save()
            else:
                if run_xtb and self.config.preopt_inputs and not self.state.get("preopt_done"):
                    self._progress("xTB", 14.0, "running", "Pre-optimizing input monomers with xTB.")
                    _log("Stabilizing input geometries using xTB before clustering...")
                    
                    if self.fragA.n_atoms() > 0:
                        p = os.path.join(self.dirs["inputs"], "preopt_anchor_start.xyz")
                        self.fragA.to_xyz(p)
                        res = self._run_engine_via_wrapper(p, "xtb", f"{self.config.xtb_method} --opt", log_cb, status_cb)
                        if res["status"] == "success":
                            self.fragA = Molecule.from_xyz(res["path"])
                            shutil.copy2(res["path"], os.path.join(self.dirs["inputs"], "stabilized_anchor.xyz"))
                            _log("-> Anchor geometry stabilized.")
                            
                    for i, mol in enumerate(self.fragBs):
                        if mol.n_atoms() == 0:
                            continue
                        p = os.path.join(self.dirs["inputs"], f"preopt_guest_{i}_start.xyz")
                        mol.to_xyz(p)
                        res = self._run_engine_via_wrapper(p, "xtb", f"{self.config.xtb_method} --opt", log_cb, status_cb)
                        if res["status"] == "success":
                            self.fragBs[i] = Molecule.from_xyz(res["path"])
                            shutil.copy2(res["path"], os.path.join(self.dirs["inputs"], f"stabilized_guest_{i}.xyz"))
                            _log(f"-> Guest {i} geometry stabilized.")
                            
                    self.state["preopt_done"] = True
                    self.save()

                if not self.state["gen"]:
                    self._progress("Sampling/filtering", 18.0, "running", f"Generating {self.config.n_generate} H-bonded seeds.")
                    _log(f"Generating {self.config.n_generate} H-bonded seeds...")
                    mols = []
                    for i in range(self.config.n_generate):
                        if i % max(1, self.config.n_generate // 20) == 0:
                            self._progress("Sampling/filtering", 18.0 + 7.0 * (i / max(self.config.n_generate, 1)), "running", f"Generated {i}/{self.config.n_generate} seeds.")
                        m = generate_cluster(
                            self.fragA,
                            self.fragBs,
                            self.n_guests_list,
                            random.Random(self.config.random_seed + i),
                            max_att=self.config.max_placement_attempts,
                            n_anchor=self.n_anchor,
                            reaction_type=self.reaction_type,
                        )
                        if m and m.n_atoms() > 0:
                            valid_atoms, atom_message = self._validate_generated_cluster(m)
                            if not valid_atoms:
                                _log(f"[ATOM-BALANCE WARNING] {atom_message}")
                            p = os.path.join(self.dirs["inputs"], f"cluster_{self.ratio_str}_raw_{i:03d}.xyz")
                            m.to_xyz(p)
                            mols.append(p)
                    self.state["gen"] = mols
                    self.save()
                    if not mols:
                        raise RuntimeError(f"Failed to generate ANY valid collision-free clusters for Ratio {self.ratio_str}.")

                if not self.state["filt"]:
                    self._progress("Sampling/filtering", 26.0, "running", "Scoring geometries and filtering RMSD duplicates.")
                    _log("Scoring geometries and filtering out RMSD duplicates...")
                    mols = [Molecule.from_xyz(p) for p in self.state["gen"]]
                    for m in mols:
                        score_geometry(m)
                    mols.sort(key=lambda x: x.score, reverse=True)
                    unique = []
                    for m in mols[: self.config.n_keep_scored]:
                        if not any(
                            u.n_atoms() == m.n_atoms() and
                            kabsch_rmsd(m.coords_array(), u.coords_array()) < self.config.rmsd_cutoff
                            for u in unique
                        ):
                            unique.append(m)
                    for i, m in enumerate(unique[: self.config.n_keep_clustered]):
                        p = os.path.join(self.dirs["inputs"], f"cluster_{self.ratio_str}_filt_{i:03d}.xyz")
                        m.to_xyz(p, f"Score: {m.score:.2f}")
                        self.state["filt"].append(p)
                    self.save()

            self._progress("Sampling/filtering", 29.0, "completed", "Cluster seeding and filtering complete.")

            if run_xtb:
                _log("\n=======================================================")
                _log(f"Running xTB Optimization on top {self.config.n_run_xtb} clustered geometries...")
                _log("=======================================================")
                targets = self.state["filt"][: self.config.n_run_xtb]
                total = max(len(targets), 1)
                for idx, p in enumerate(targets, start=1):
                    self._progress("xTB", 30.0 + 20.0 * ((idx - 1) / total), "running", f"xTB optimization {idx}/{total}.")
                    stem = Path(p).stem
                    if stem not in self.state["xtb"] or self.state["xtb"][stem].get("status") != "success":
                        res = self._run_engine_via_wrapper(p, "xtb", f"{self.config.xtb_method} --opt", log_cb, status_cb)
                        if res["status"] == "success":
                            final_name = os.path.join(self.dirs["xtb"], f"xtbopt_{stem}.xyz")
                            shutil.copy2(res["path"], final_name)
                            res["path"] = final_name
                        self.state["xtb"][stem] = res
                        self.save()
                self._progress("xTB", 51.0, "completed", "xTB optimizations complete.")

            if run_crest:
                self._progress("CREST", 52.0, "running", "Filtering xTB structures and promoting to CREST.")
                _log("\n=======================================================")
                _log("Filtering xTB results and promoting to CREST...")
                xtb_success = [(k, v) for k, v in self.state["xtb"].items() if v.get("status") == "success"]
                if not xtb_success:
                    _log("[ERROR] No successful xTB runs available for CREST. Did xTB fail?")
                else:
                    xtb_success.sort(key=lambda x: x[1].get("energy", 0))
                    min_e = xtb_success[0][1].get("energy", 0)
                    promoted_xtb = [(k, v) for k, v in xtb_success if (v["energy"] - min_e) * EH2KCAL <= self.config.xtb_ewin_kcal][: self.config.n_run_crest]
                    _log(f"Promoted {len(promoted_xtb)} xTB structures within {self.config.xtb_ewin_kcal} kcal/mol of minimum.")
                    total = max(len(promoted_xtb), 1)
                    for idx, (k, v) in enumerate(promoted_xtb, start=1):
                        self._progress("CREST", 55.0 + 15.0 * ((idx - 1) / total), "running", f"CREST search {idx}/{total}.")
                        if k not in self.state["crest"] or self.state["crest"][k].get("status") != "success":
                            crest_flags = f"{self.config.crest_method} -T {self.config.cores}"
                            res = self._run_engine_via_wrapper(v["path"], "crest", crest_flags, log_cb, status_cb)
                            if res["status"] != "success":
                                _log("[WARNING] CREST trial MD did not converge. Retrying with reduced time step (--tstep 1 --mdlen 0.5)...")
                                res = self._run_engine_via_wrapper(v["path"], "crest", f"{crest_flags} --tstep 1 --mdlen 0.5", log_cb, status_cb)
                            if res["status"] == "success":
                                final_conf = os.path.join(self.dirs["crest"], f"crest_conformers_{k}.xyz")
                                shutil.copy2(res["path"], final_conf)
                                res["path"] = final_conf
                                if "best_path" in res and os.path.exists(res["best_path"]):
                                    final_best = os.path.join(self.dirs["crest"], f"crest_best_{k}.xyz")
                                    shutil.copy2(res["best_path"], final_best)
                                    res["best_path"] = final_best
                            self.state["crest"][k] = res
                            self.save()
                self._progress("CREST", 71.0, "completed", "CREST conformer search complete.")

            if run_orca:
                self._progress("ORCA", 72.0, "running", "Aggregating CREST conformers and promoting to ORCA.")
                _log("\n=======================================================")
                _log("Aggregating CREST conformers and promoting to ORCA DFT...")
                all_confs = []
                for k, v in self.state["crest"].items():
                    if v.get("status") == "success" and "path" in v and os.path.exists(v["path"]):
                        mols = read_multi_xyz(v["path"])
                        for i, m in enumerate(mols):
                            m.lineage = [k, f"conf_{i}"]
                            all_confs.append(m)
                if not all_confs:
                    _log("[WARNING] No CREST conformers found. Falling back to best xTB-optimized structure(s) for ORCA.")
                    for k, v in self.state["xtb"].items():
                        if v.get("status") == "success" and "path" in v and os.path.exists(v["path"]):
                            m = Molecule.from_xyz(v["path"])
                            m.lineage = [k, "xtb_direct"]
                            m.energy_eh = v.get("energy", 0.0)
                            all_confs.append(m)
                if not all_confs:
                    _log("[ERROR] No structures available for ORCA. Ensure xTB or CREST ran successfully.")
                else:
                    all_confs.sort(key=lambda x: x.energy_eh)
                    min_crest_e = all_confs[0].energy_eh
                    unique_promoted = []
                    for m in all_confs:
                        rel_e = (m.energy_eh - min_crest_e) * EH2KCAL
                        if rel_e > self.config.crest_ewin_kcal:
                            continue
                        is_dup = any(
                            u.n_atoms() == m.n_atoms() and
                            kabsch_rmsd(m.coords_array(), u.coords_array()) < 0.25
                            for u in unique_promoted
                        )
                        if not is_dup and len(unique_promoted) < 10:
                            unique_promoted.append(m)
                    _log(f"Promoted top {len(unique_promoted)} highly unique conformers to ORCA.")
                    total = max(len(unique_promoted), 1)
                    for i, m in enumerate(unique_promoted):
                        self._progress("ORCA", 75.0 + 20.0 * (i / total), "running", f"ORCA refinement {i + 1}/{total}.")
                        stem = f"orca_target_{i:03d}"
                        if stem not in self.state["orca"] or self.state["orca"][stem].get("status") != "success":
                            inp_path = os.path.join(self.dirs["inputs"], f"{stem}.xyz")
                            m.to_xyz(inp_path)
                            res = self._run_orca_via_sandbox(inp_path, stem, log_cb, status_cb)
                            if res["status"] == "success":
                                res["lineage"] = " -> ".join(m.lineage)
                                if res["imag"] == 0:
                                    best_name = os.path.join(self.dirs["orca"], f"ORCA_MINIMUM_{stem}.xyz")
                                    shutil.copy2(res["path"], best_name)
                                    res["best_path"] = best_name
                                    _log("-> True Local Minimum Confirmed!")
                            self.state["orca"][stem] = res
                            self.save()

            self._progress("Reports", 96.0, "running", "Extracting top models and writing reports.")
            self._extract_top_models(log_cb)
            self._generate_markdown_report()
            self.job_summary.status = "SUCCESS"
            self.job_summary.finished_at = _now()
            self.job_summary.log_event("Complete", "Pipeline finished and final artifacts were written.")
            self.job_summary.write()
            self._progress("Complete", 100.0, "completed", "Job complete.")
            _log(f"\nPipeline Complete for Ratio {self.ratio_str.replace('_', ':')}. Check the '5. TOP 3 COMPARE' tab.")
        except Exception as exc:
            self.job_summary.status = "FAILED"
            self.job_summary.finished_at = _now()
            self.job_summary.error = str(exc)
            self.job_summary.traceback = traceback.format_exc()
            self.job_summary.log_event("Failure", str(exc), level="ERROR")
            self.job_summary.write()
            self._progress("Failed", max(self.progress_percent, 1.0), "failed", str(exc))
            raise

    def _run_engine_via_wrapper(self, xyz_path, engine, flags, log_cb, status_cb):
        fname = os.path.basename(xyz_path)
        stem = Path(xyz_path).stem
        wd = os.path.abspath(os.path.join(self.dirs.get(engine, self.dirs["inputs"]), f"job_{stem}"))
        os.makedirs(wd, exist_ok=True)
        shutil.copy2(xyz_path, os.path.join(wd, fname))
        # Add charge/multiplicity flags before platform branching so both Windows and Linux get them
        if engine in ("xtb", "crest"):
            if self.config.charge != 0:
                flags += f" --chrg {self.config.charge}"
            if self.config.multiplicity != 1:
                flags += f" --uhf {self.config.multiplicity - 1}"
        if sys.platform == "win32":
            sh_path = os.path.abspath(os.path.join(wd, f"run_{engine}.sh"))
            with open(sh_path, "w", newline="\n", encoding="utf-8") as f:
                f.write("#!/bin/bash\n")
                wsl_wd = get_wsl_path(wd).replace("'", "'\\''")
                f.write(f"SANDBOX=\"/tmp/iak_{engine}_{stem}_$RANDOM\"\n")
                f.write("mkdir -p \"$SANDBOX\"\n")
                f.write(f"cp '{wsl_wd}/{fname}' \"$SANDBOX/\"\n")
                f.write("cd \"$SANDBOX\"\n")
                if engine == "xtb" and _engines.XTB_DIR:
                    wsl_xtb = get_wsl_path(os.path.abspath(_engines.XTB_DIR))
                    f.write(f"export PATH='{wsl_xtb}:/usr/bin:/bin:/usr/local/bin:$PATH'\n")
                    f.write(f"export XTBPATH='{get_wsl_path(os.path.abspath(os.path.join(_engines.XTB_DIR, '..', 'share', 'xtb')))}'\n")
                    exec_cmd = f"'{wsl_xtb}/xtb'"
                elif engine == "crest" and _engines.CREST_DIR:
                    wsl_crest = get_wsl_path(os.path.abspath(_engines.CREST_DIR))
                    f.write(f"export PATH='{wsl_crest}:/usr/bin:/bin:/usr/local/bin:$PATH'\n")
                    exec_cmd = f"'{wsl_crest}/crest'"
                else:
                    f.write("source ~/.bashrc 2>/dev/null\n")
                    f.write("source ~/.profile 2>/dev/null\n")
                    f.write("export PATH='/usr/bin:/bin:/usr/local/bin:$PATH'\n")
                    exec_cmd = engine
                f.write(f"export OMP_STACKSIZE=1G\nexport OMP_NUM_THREADS={self.config.cores}\nulimit -s unlimited\n")
                if engine == "crest":
                    f.write("export OMPI_ALLOW_RUN_AS_ROOT=1\nexport OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1\n")
                    f.write("export OMPI_MCA_btl_vader_single_copy_mechanism=none\nexport OMPI_MCA_btl=\"^openib\"\n")
                    f.write("export OMPI_MCA_rmaps_base_oversubscribe=1\nexport OMPI_MCA_hwloc_base_binding_policy=none\n")
                f.write(f"{exec_cmd} '{fname}' {flags}\n")
                f.write(f"/bin/cp -r * '{wsl_wd}/' 2>/dev/null\n")
                f.write("cd /tmp\n")
                f.write("/bin/rm -rf \"$SANDBOX\"\n")
            cmd = f"wsl -e bash \"{get_wsl_path(sh_path)}\""
        else:
            # Build env-var prefix and select correct binary on Linux/macOS
            env_exports = (
                f"export OMP_STACKSIZE=1G; "
                f"export OMP_NUM_THREADS={self.config.cores}; "
                "ulimit -s unlimited 2>/dev/null; "
            )
            if engine == "crest":
                env_exports += (
                    "export OMPI_ALLOW_RUN_AS_ROOT=1; "
                    "export OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1; "
                    "export OMPI_MCA_btl_vader_single_copy_mechanism=none; "
                    "export OMPI_MCA_btl='^openib'; "
                    "export OMPI_MCA_rmaps_base_oversubscribe=1; "
                    "export OMPI_MCA_hwloc_base_binding_policy=none; "
                )
            if engine == "xtb" and _engines.XTB_DIR:
                xtb_abs = os.path.abspath(_engines.XTB_DIR)
                xtb_share = os.path.abspath(os.path.join(_engines.XTB_DIR, "..", "share", "xtb"))
                env_exports += (
                    f"export PATH='{xtb_abs}:/usr/bin:/bin:/usr/local/bin:$PATH'; "
                    f"export XTBPATH='{xtb_share}'; "
                )
                exec_cmd_linux = f"'{os.path.join(xtb_abs, 'xtb')}'"
            elif engine == "crest" and _engines.CREST_DIR:
                crest_abs = os.path.abspath(_engines.CREST_DIR)
                env_exports += f"export PATH='{crest_abs}:/usr/bin:/bin:/usr/local/bin:$PATH'; "
                exec_cmd_linux = f"'{os.path.join(crest_abs, 'crest')}'"
            else:
                exec_cmd_linux = engine
            cmd = f"cd '{wd}' && {env_exports}{exec_cmd_linux} '{fname}' {flags}"

        if status_cb:
            status_cb(engine, 1)
        if log_cb:
            log_cb(f"\n>>> Executing Sandbox Engine: {engine} {flags}\n")
        started = time.time()
        try:
            proc = subprocess.Popen(cmd, shell=True, cwd=wd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, encoding="utf-8", errors="replace", bufsize=1)
            output = []
            for line in iter(proc.stdout.readline, ""):
                output.append(line)
                if log_cb:
                    log_cb(line)
            proc.wait()
            if engine == "xtb":
                opt_file = os.path.join(wd, "xtbopt.xyz")
                if not os.path.exists(opt_file):
                    opt_file = os.path.join(wd, fname)
                energy, success = 0.0, False
                for l in reversed(output):
                    if "total energy" in l.lower():
                        try:
                            parts = l.split()
                            energy = float(parts[parts.index("Eh") - 1]) if "Eh" in parts else float(parts[-2])
                            success = True
                            break
                        except Exception:
                            pass
                res = {"status": "success", "energy": energy, "path": opt_file} if success or os.path.exists(os.path.join(wd, "xtbopt.xyz")) else {"status": "failed"}
                self._record_engine_result(engine, stem, flags, res, time.time() - started)
                return res
            if engine == "crest":
                conf_file = os.path.join(wd, "crest_conformers.xyz")
                best_file = os.path.join(wd, "crest_best.xyz")
                res = {"status": "success", "path": conf_file, "best_path": best_file if os.path.exists(best_file) else conf_file} if os.path.exists(conf_file) else {"status": "failed"}
                self._record_engine_result(engine, stem, flags, res, time.time() - started)
                return res
        except Exception as e:
            if log_cb:
                log_cb(f"\n[Python Error] {e}\n")
        finally:
            if status_cb:
                status_cb(engine, -1)
        res = {"status": "failed"}
        self._record_engine_result(engine, stem, flags, res, time.time() - started)
        return res

    def _run_orca_via_sandbox(self, xyz_path, stem, log_cb, status_cb):
        wd = os.path.abspath(os.path.join(self.dirs["orca"], f"job_{stem}"))
        os.makedirs(wd, exist_ok=True)
        wsl_wd = get_wsl_path(wd).replace("'", "'\\''")
        mol = Molecule.from_xyz(xyz_path)

        if sys.platform == "win32" and not _engines.ORCA_IS_WINDOWS:
            has_mpi = subprocess.run("wsl -e bash -c 'export PATH=\"/usr/bin:/bin:/usr/local/bin:$PATH\"; which mpirun'", shell=True, capture_output=True).returncode == 0
        else:
            has_mpi = shutil.which("mpirun") is not None

        def _execute_orca(use_mpi: bool, skip_freq: bool):
            inp_file = os.path.join(wd, f"{stem}.inp")
            with open(inp_file, "w", encoding="utf-8", newline="\n") as f:
                base_method = self.config.orca_method.replace("Opt Freq", "").replace("TightOpt", "").replace("Opt", "").strip()
                task = "Opt" if skip_freq else "Opt Freq"
                f.write(f"! {base_method} {task}\n")
                f.write(f"%maxcore {self.config.maxcore}\n")
                if use_mpi:
                    f.write(f"%pal nprocs {self.config.cores} end\n")
                f.write("%geom\n  MaxIter 350\nend\n")
                f.write("%scf\n  MaxIter 350\nend\n")
                f.write(f"* xyz {self.config.charge} {self.config.multiplicity}\n")
                for a in mol.atoms:
                    f.write(f"{a.symbol:<4} {a.x:15.6f} {a.y:15.6f} {a.z:15.6f}\n")
                f.write("*\n")

            if sys.platform == "win32" and not _engines.ORCA_IS_WINDOWS:
                sh_path = os.path.abspath(os.path.join(wd, "run_orca.sh"))
                with open(sh_path, "w", newline="\n", encoding="utf-8") as f:
                    f.write("#!/bin/bash\n")
                    f.write(f"SANDBOX=\"/tmp/iak_orca_{stem}_$RANDOM\"\n")
                    f.write("mkdir -p \"$SANDBOX\"\n")
                    f.write(f"cp '{wsl_wd}/{stem}.inp' \"$SANDBOX/\"\n")
                    f.write("cd \"$SANDBOX\"\n")
                    f.write("export PATH='/usr/bin:/bin:/usr/local/bin:$PATH'\n")
                    f.write("ulimit -s unlimited 2>/dev/null\n")
                    f.write("export OMP_NUM_THREADS=1\nexport OMP_STACKSIZE=1G\n")
                    if _engines.ORCA_DIR:
                        wsl_orca = get_wsl_path(os.path.abspath(_engines.ORCA_DIR))
                        f.write("SAFE_ORCA=\"/tmp/iak_orca_bin\"\n")
                        f.write("/bin/rm -f $SAFE_ORCA\n")
                        f.write(f"/bin/ln -s '{wsl_orca}' $SAFE_ORCA\n")
                        f.write("export PATH=$SAFE_ORCA:$PATH\n")
                        f.write("export LD_LIBRARY_PATH=$SAFE_ORCA:/usr/lib/x86_64-linux-gnu:/usr/lib:$LD_LIBRARY_PATH\n")
                        exec_cmd = "orca"
                    else:
                        f.write("source ~/.bashrc 2>/dev/null\nsource ~/.profile 2>/dev/null\n")
                        exec_cmd = "orca"
                    if use_mpi:
                        f.write("REAL_MPI=$(which mpirun)\n")
                        f.write("if [ -z \"$REAL_MPI\" ]; then REAL_MPI=\"/usr/bin/mpirun\"; fi\n")
                        f.write("echo '#!/bin/bash' > ./mpirun\n")
                        f.write("echo \"$REAL_MPI \\\"\\$@\\\"\" >> ./mpirun\n")
                        f.write("chmod +x ./mpirun\n")
                        f.write("export PATH=\".:$PATH\"\n")
                    f.write("export OMPI_ALLOW_RUN_AS_ROOT=1\nexport OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1\n")
                    f.write("export OMPI_MCA_btl_vader_single_copy_mechanism=none\nexport OMPI_MCA_btl=\"^openib\"\n")
                    f.write("export OMPI_MCA_rmaps_base_oversubscribe=1\nexport OMPI_MCA_hwloc_base_binding_policy=none\n")
                    f.write(f"{exec_cmd} {stem}.inp > '{wsl_wd}/{stem}.out' 2>&1\n")
                    f.write(f"/bin/cp {stem}_trj.xyz '{wsl_wd}/' 2>/dev/null\n")
                    f.write(f"/bin/cp *xyz '{wsl_wd}/' 2>/dev/null\n")
                    f.write("cd /tmp\n")
                    f.write("/bin/rm -rf \"$SANDBOX\"\n")
                cmd = f"wsl -e bash \"{get_wsl_path(sh_path)}\""
            elif sys.platform == "win32" and _engines.ORCA_IS_WINDOWS:
                cmd = f"cd /d \"{wd}\" && \"{os.path.join(os.path.abspath(_engines.ORCA_DIR), 'orca.exe')}\" {stem}.inp > {stem}.out 2>&1"
            else:
                cmd = f"cd '{wd}' && {ORCA_CMD} {stem}.inp > {stem}.out 2>&1"

            try:
                proc = subprocess.Popen(cmd, shell=True, cwd=wd)
                out_file = os.path.join(wd, f"{stem}.out")
                last_pos = 0
                while proc.poll() is None:
                    if os.path.exists(out_file):
                        try:
                            with open(out_file, "r", encoding="utf-8", errors="replace") as f:
                                f.seek(last_pos)
                                chunk = f.read()
                                if chunk and log_cb:
                                    for line in chunk.splitlines(True):
                                        log_cb(line)
                                last_pos = f.tell()
                        except Exception:
                            pass
                    time.sleep(0.5)
                if os.path.exists(out_file):
                    with open(out_file, "r", encoding="utf-8", errors="replace") as f:
                        f.seek(last_pos)
                        chunk = f.read()
                        if chunk and log_cb:
                            for line in chunk.splitlines(True):
                                log_cb(line)
                trj_file = os.path.join(wd, f"{stem}_trj.xyz")
                energy, gibbs, imag, success = 0.0, 0.0, 0, False
                mpirun_crashed, scf_crashed = False, False
                time.sleep(1.0)
                if os.path.exists(out_file):
                    with open(out_file, "r", encoding="utf-8", errors="replace") as f:
                        content = f.read()
                    if "FINAL SINGLE POINT ENERGY" in content:
                        for line in content.split("\n"):
                            if "FINAL SINGLE POINT ENERGY" in line:
                                energy = float(line.split()[-1])
                            if "Final Gibbs free energy" in line:
                                gibbs = float(line.split()[-2])
                            if "*** imaginary mode ***" in line:
                                imag += 1
                            if "ORCA TERMINATED NORMALLY" in line:
                                success = True
                    content_lower = content.lower()
                    if "mpirun: not found" in content_lower or "aborting the run" in content_lower or "mpi" in content_lower or "hwloc" in content_lower:
                        mpirun_crashed = True
                    if "scf not converged" in content_lower or "optimization did not converge" in content_lower or "std::bad_alloc" in content_lower or "command not found" in content_lower:
                        scf_crashed = True
                    if not success and not mpirun_crashed:
                        scf_crashed = True
                else:
                    scf_crashed = True
                if success and energy != 0.0:
                    opt_path = trj_file if os.path.exists(trj_file) else os.path.join(wd, f"{stem}.xyz")
                    return True, {"status": "success", "energy": energy, "gibbs": gibbs, "imag": imag, "path": opt_path}, False, False
                return False, {"status": "failed"}, mpirun_crashed, scf_crashed
            except Exception:
                return False, {"status": "failed"}, True, True

        if status_cb:
            status_cb("orca", 1)
        started = time.time()
        try:
            if log_cb:
                log_cb(f"\n>>> [Tier 1] Executing ORCA (Parallel {self.config.cores}-Core, Opt Freq): {stem}\n")
            success, res, mpi_crash, scf_crash = _execute_orca(use_mpi=has_mpi, skip_freq=False)
            if not success:
                if mpi_crash and has_mpi:
                    if log_cb:
                        log_cb(f"\n>>> [Tier 2] MPI Failed! Restarting in Serial Mode (1-Core, Opt Freq): {stem}\n")
                    success, res, mpi_crash, scf_crash = _execute_orca(use_mpi=False, skip_freq=False)
                if not success:
                    use_mpi_for_opt = has_mpi and not mpi_crash
                    cores = f"{self.config.cores}-Core" if use_mpi_for_opt else "1-Core"
                    if log_cb:
                        log_cb(f"\n>>> [Tier 3] Math/Memory Failed! Stripping Freq, forcing Opt-Only ({cores}): {stem}\n")
                    success, res, mpi_crash, scf_crash = _execute_orca(use_mpi=use_mpi_for_opt, skip_freq=True)
                    if not success and use_mpi_for_opt:
                        if log_cb:
                            log_cb(f"\n>>> [Tier 4] Ultimate Fallback: Serial Mode Opt-Only (1-Core): {stem}\n")
                        success, res, mpi_crash, scf_crash = _execute_orca(use_mpi=False, skip_freq=True)
            self._record_engine_result("orca", stem, self.config.orca_method, res, time.time() - started)
            return res
        finally:
            if status_cb:
                status_cb("orca", -1)

    def _extract_top_models(self, log_cb):
        comp_dir = self.dirs["comparison"]
        report = []
        best_xtb_k, best_xtb_e, best_xtb_path = None, float("inf"), None
        for k, v in self.state["xtb"].items():
            if v.get("status") == "success" and "energy" in v and v["energy"] < best_xtb_e:
                best_xtb_e, best_xtb_k, best_xtb_path = v["energy"], k, v["path"]
        if best_xtb_path and os.path.exists(best_xtb_path):
            shutil.copy2(best_xtb_path, os.path.join(comp_dir, "1_BEST_xTB.xyz"))
            report.append(f"xTB,{best_xtb_k},{best_xtb_e:.6f},N/A,N/A")

        best_crest_e, best_crest_path = float("inf"), None
        for k, v in self.state["crest"].items():
            if v.get("status") == "success" and "best_path" in v and os.path.exists(v["best_path"]):
                mols = read_multi_xyz(v["best_path"])
                if mols and mols[0].energy_eh < best_crest_e:
                    best_crest_e, best_crest_path = mols[0].energy_eh, v["best_path"]
        if best_crest_path:
            shutil.copy2(best_crest_path, os.path.join(comp_dir, "2_BEST_CREST.xyz"))
            report.append(f"CREST,-,{best_crest_e:.6f},N/A,N/A")

        best_orca_k, best_orca_e, best_orca_g, best_orca_path = None, float("inf"), 0.0, None
        for k, v in self.state["orca"].items():
            if v.get("status") == "success" and v.get("imag") == 0 and "energy" in v:
                rank_metric = v.get("gibbs") if v.get("gibbs", 0) != 0.0 else v["energy"]
                if rank_metric < best_orca_e:
                    best_orca_e, best_orca_g, best_orca_k, best_orca_path = rank_metric, v.get("gibbs", 0.0), k, v.get("best_path", v["path"])
        if best_orca_path is None:
            for k, v in self.state["orca"].items():
                if v.get("status") == "success" and "energy" in v:
                    rank_metric = v.get("gibbs") if v.get("gibbs", 0) != 0.0 else v["energy"]
                    if rank_metric < best_orca_e:
                        best_orca_e, best_orca_g, best_orca_k, best_orca_path = rank_metric, v.get("gibbs", 0.0), k, v.get("best_path", v["path"])
        if best_orca_path and os.path.exists(best_orca_path):
            shutil.copy2(best_orca_path, os.path.join(comp_dir, "3_BEST_ORCA.xyz"))
            binding_str = "N/A"
            if self.config.energy_a is not None and self.config.energy_b is not None:
                total_guests = sum(self.n_guests_list)
                delta_e_eh = best_orca_e - (self.config.energy_a + (total_guests * self.config.energy_b))
                binding_str = f"{delta_e_eh * EH2KCAL:.2f}"
            report.append(f"ORCA,{best_orca_k},{best_orca_e:.6f},{best_orca_g:.6f},{binding_str}")
        if report:
            with open(os.path.join(comp_dir, "Energy_Comparison.csv"), "w", encoding="utf-8") as f:
                f.write("Level_of_Theory,Source_ID,Total_Energy_Eh,Gibbs_Free_Energy_Eh,Binding_Energy_kcal_mol\n")
                for line in report:
                    f.write(line + "\n")
            if log_cb:
                log_cb("\n-> Extracted Top Models into 05_Top_Models_Comparison\n")

    def _generate_markdown_report(self):
        report_path = os.path.join(self.dirs["comparison"], "Post_CREST_Report.md")
        with open(report_path, "w", encoding="utf-8") as f:
            f.write(f"# Workflow Analysis Report: Ratio {self.ratio_str.replace('_', ':')}\n\n")
            f.write("## Scientific Disclaimer\n")
            f.write("> The structures reported herein represent local minima found within the defined sampling depth and selected level of theory. They are **not guaranteed global minima**.\n\n")
            all_orca = self.state.get("orca", {})
            valid_orca = [(k, v) for k, v in all_orca.items() if v.get("status") == "success"]
            true_minima = [(k, v) for k, v in valid_orca if v.get("imag") == 0]
            f.write("## 1. Summary Metrics\n")
            f.write(f"- **ORCA Method Used:** {self.config.orca_method}\n")
            f.write(f"- **System Charge:** {self.config.charge}  |  **Multiplicity:** {self.config.multiplicity}\n")
            f.write(f"- **ORCA Conformers Evaluated:** {len(all_orca)}\n")
            f.write(f"- **Successful ORCA Completions:** {len(valid_orca)}\n")
            f.write(f"- **True Minima (0 Imag Freqs):** {len(true_minima)}\n\n")
            if not true_minima and valid_orca:
                f.write("> **WARNING:** No true minima were found. The best structure has been extracted anyway for inspection.\n\n")
            elif not valid_orca and len(all_orca) > 0:
                f.write("> **CRITICAL FAILURE:** All ORCA jobs crashed or failed to converge. Check `.out` files in `04_ORCA_Refinement`.\n\n")
            f.write("## 2. Minima Ranking (Gibbs Free Energy)\n")
            if valid_orca:
                valid_orca.sort(key=lambda x: x[1].get("gibbs") if x[1].get("gibbs", 0) != 0 else x[1].get("energy"))
                best_g = valid_orca[0][1].get("gibbs") if valid_orca[0][1].get("gibbs", 0) != 0 else valid_orca[0][1].get("energy")
                f.write("| ORCA ID | Parent Lineage | ΔG (kcal/mol) | Imag Freqs | Status |\n")
                f.write("|---|---|---|---|---|\n")
                for k, v in valid_orca:
                    val = v.get("gibbs") if v.get("gibbs", 0) != 0 else v.get("energy")
                    rel_g = (val - best_g) * EH2KCAL
                    f.write(f"| {k} | {v.get('lineage','N/A')} | {rel_g:.2f} | {v.get('imag')} | Success |\n")
                for k, v in all_orca.items():
                    if v.get("status") != "success":
                        f.write(f"| {k} | {v.get('lineage','N/A')} | N/A | N/A | FAILED |\n")
            else:
                f.write("*No valid minima found to rank.*\n")

