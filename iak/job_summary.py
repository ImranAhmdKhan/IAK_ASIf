from __future__ import annotations

import csv
import json
import os
import sys
from typing import Any, Dict, List

from .constants import normalize_reaction_type, REACTION_TYPE_CHOICES
from .engines import is_tool_available
from .utils import _now, _safe_int, _validate_xyz
from .chemistry import validate_charge_multiplicity


class JobSummaryRecorder:
    def __init__(self, pipeline):
        self.pipeline = pipeline
        self.ratio = getattr(pipeline, "ratio_str", "unknown").replace("_", ":")
        self.started_at = _now()
        self.finished_at = ""
        self.status = "RUNNING"
        self.error = ""
        self.traceback = ""
        self.inputs: List[Dict[str, Any]] = []
        self.workflow: Dict[str, Any] = {}
        self.events: List[Dict[str, Any]] = []
        self.engine_results: List[Dict[str, Any]] = []

    @property
    def logs_dir(self):
        return self.pipeline.dirs.get("logs", os.path.join(self.pipeline.out_dir, "logs"))

    @property
    def comparison_dir(self):
        return self.pipeline.dirs.get("comparison", os.path.join(self.pipeline.out_dir, "05_Top_Models_Comparison"))

    @property
    def json_path(self):
        return os.path.join(self.logs_dir, f"Job_Summary_{self.ratio.replace(':', '_')}.json")

    @property
    def csv_path(self):
        return os.path.join(self.logs_dir, f"Job_Engine_Results_{self.ratio.replace(':', '_')}.csv")

    @property
    def markdown_path(self):
        return os.path.join(self.comparison_dir, "Job_Summary_Report.md")

    def capture_inputs(self, run_xtb, run_crest, run_orca):
        cfg = self.pipeline.config
        n_anchor = int(getattr(self.pipeline, "n_anchor", 1))
        n_guests_list = [int(x) for x in getattr(self.pipeline, "n_guests_list", [])]
        guest_paths = list(getattr(self.pipeline, "fragB_paths", []) or [])
        guest_required = sum(max(x, 0) for x in n_guests_list) > 0
        charge = _safe_int(cfg.charge)
        mult = _safe_int(cfg.multiplicity)
        cores = _safe_int(cfg.cores)
        maxcore = _safe_int(cfg.maxcore)
        reaction_type = normalize_reaction_type(getattr(self.pipeline, "reaction_type", "Non-covalent"))

        ratio_valid = n_anchor >= 0 and all(x >= 0 for x in n_guests_list)
        electron_check = {
            "valid": mult is not None and mult >= 1,
            "total_electrons": None,
            "unknown_symbols": [],
            "message": "Electron consistency check not executed.",
        }
        if charge is not None and mult is not None and mult >= 1:
            electron_check = validate_charge_multiplicity(
                self.pipeline.fragA,
                self.pipeline.fragBs,
                n_anchor,
                n_guests_list,
                charge,
                mult,
            )

        self.inputs = [
            _validate_xyz(getattr(self.pipeline, "fragA_path", ""), "Anchor A XYZ", n_anchor > 0),
            self._setting(
                "Stoichiometric ratio",
                self.ratio,
                True,
                ratio_valid,
                f"Anchor copies={n_anchor}; guest copies={n_guests_list or [0]}.",
            ),
            self._setting(
                "Charge",
                cfg.charge,
                True,
                charge is not None and (electron_check["total_electrons"] is None or electron_check["total_electrons"] > 0),
                "Integer charge accepted." if charge is not None else "Charge must be an integer.",
            ),
            self._setting(
                "Multiplicity",
                cfg.multiplicity,
                True,
                mult is not None and mult >= 1 and electron_check["valid"],
                electron_check["message"] if mult is not None and mult >= 1 else "Multiplicity must be >= 1.",
            ),
            self._setting(
                "Total Electrons",
                electron_check["total_electrons"] if electron_check["total_electrons"] is not None else "unknown",
                True,
                electron_check["valid"],
                electron_check["message"],
            ),
            self._setting("CPU cores", cfg.cores, True, cores is not None and cores >= 1, "Core count accepted." if cores is not None and cores >= 1 else "Cores must be >= 1."),
            self._setting("RAM per core", cfg.maxcore, True, maxcore is not None and maxcore >= 128, "MaxCore accepted." if maxcore is not None and maxcore >= 128 else "MaxCore is unusually low."),
            self._setting("ORCA method", cfg.orca_method, run_orca, bool(str(cfg.orca_method).strip()) or not run_orca, "Theory level recorded." if cfg.orca_method else "ORCA method is empty."),
            self._setting("Reaction Type", reaction_type, True, reaction_type in REACTION_TYPE_CHOICES, "Reaction mechanism selected."),
        ]

        if guest_paths:
            for idx, guest_path in enumerate(guest_paths):
                required = idx < len(n_guests_list) and n_guests_list[idx] > 0
                self.inputs.append(_validate_xyz(guest_path, f"Guest B{idx + 1} XYZ", required))
        else:
            self.inputs.append(
                self._setting(
                    "Guest queue",
                    "(empty)",
                    guest_required,
                    not guest_required,
                    "No guests queued." if not guest_required else "Guests are required by ratio but no guest files are queued.",
                )
            )

        if run_xtb:
            self.inputs.append(self._engine("xTB engine", "xtb"))
        if run_crest:
            self.inputs.append(self._engine("CREST engine", "crest"))
        if run_orca:
            self.inputs.append(self._engine("ORCA engine", "orca"))
        self.workflow = {
            "ratio": self.ratio,
            "output_directory": self.pipeline.out_dir,
            "run_xtb": bool(run_xtb),
            "run_crest": bool(run_crest),
            "run_orca": bool(run_orca),
            "reaction_type": reaction_type,
            "preopt_inputs": bool(cfg.preopt_inputs),
            "n_generate": cfg.n_generate,
            "n_keep_scored": cfg.n_keep_scored,
            "n_keep_clustered": cfg.n_keep_clustered,
            "n_run_xtb": cfg.n_run_xtb,
            "n_run_crest": cfg.n_run_crest,
            "n_anchor": n_anchor,
            "n_guests_list": n_guests_list,
            "charge": cfg.charge,
            "multiplicity": cfg.multiplicity,
            "orca_method": cfg.orca_method,
        }

    def _setting(self, label, value, required, valid, message):
        return {"label": label, "path": str(value), "required": required, "provided": value not in (None, ""), "valid": valid, "result": "CORRECT" if valid else "INCORRECT", "message": message}

    def _engine(self, label, engine):
        available = is_tool_available(engine)
        return self._setting(label, engine, True, available, "Engine detected." if available else "Engine missing or not on PATH/WSL/local engine folder.")

    def log_event(self, stage, message, level="INFO"):
        self.events.append({"time": _now(), "level": level, "stage": stage, "message": str(message)})

    def add_engine_result(self, engine, source, flags, result, elapsed):
        self.engine_results.append(
            {
                "time": _now(),
                "engine": engine.upper(),
                "source": source,
                "flags": flags,
                "status": result.get("status", "unknown"),
                "energy_eh": result.get("energy", ""),
                "gibbs_eh": result.get("gibbs", ""),
                "imaginary_frequencies": result.get("imag", ""),
                "path": result.get("path", ""),
                "best_path": result.get("best_path", ""),
                "elapsed_seconds": round(elapsed, 2),
            }
        )

    def to_dict(self):
        return {
            "schema": "IAK_JOB_SUMMARY_V11_2",
            "ratio": self.ratio,
            "started_at": self.started_at,
            "finished_at": self.finished_at,
            "status": self.status,
            "error": self.error,
            "traceback": self.traceback,
            "workflow": self.workflow,
            "inputs": self.inputs,
            "events": self.events,
            "engine_results": self.engine_results,
        }

    def write(self):
        os.makedirs(self.logs_dir, exist_ok=True)
        os.makedirs(self.comparison_dir, exist_ok=True)
        with open(self.json_path, "w", encoding="utf-8") as handle:
            json.dump(self.to_dict(), handle, indent=2)
        with open(self.csv_path, "w", newline="", encoding="utf-8") as handle:
            fieldnames = ["time", "engine", "source", "flags", "status", "energy_eh", "gibbs_eh", "imaginary_frequencies", "path", "best_path", "elapsed_seconds"]
            writer = csv.DictWriter(handle, fieldnames=fieldnames)
            writer.writeheader()
            for row in self.engine_results:
                writer.writerow(row)
        with open(self.markdown_path, "w", encoding="utf-8") as handle:
            handle.write(self._markdown())
        opts = getattr(self.pipeline.config, "report_options", {})
        if opts.get("ai_prompt", True):
            with open(os.path.join(self.comparison_dir, "AI_Summary_Prompt.txt"), "w", encoding="utf-8") as handle:
                handle.write(
                    "You are helping summarize an IAK computational chemistry job.\n"
                    "Use only the supplied JSON/CSV/report files from this folder.\n"
                    f"Ratio: {self.ratio}\n"
                    f"Status: {self.status}\n"
                    "Host-selected report sections: "
                    + ", ".join(k for k, v in opts.items() if v)
                    + "\n\n"
                    "Requested output: concise scientific summary, key thermodynamic values if present, failures/warnings, and a CSV-ready table.\n"
                )

    def _markdown(self):
        opts = getattr(self.pipeline.config, "report_options", {})
        enabled = lambda key: opts.get(key, True)
        correct = [x for x in self.inputs if x.get("valid")]
        incorrect = [x for x in self.inputs if not x.get("valid")]
        lines = [
            f"# IAK Job Summary Report: Ratio {self.ratio}",
            "",
            f"- **Status:** {self.status}",
            f"- **Started:** {self.started_at}",
            f"- **Finished:** {self.finished_at or 'Still running'}",
            f"- **Output Directory:** `{self.pipeline.out_dir}`",
            "",
            "## High-Level Input Assessment",
            "",
            f"- Correct inputs/settings: **{len(correct)}**",
            f"- Incorrect or missing inputs/settings: **{len(incorrect)}**",
        ]
        if self.error:
            lines.extend(["", f"> **Job Error:** {self.error}"])
        if enabled("inputs"):
            lines.extend(["", "## Input Correctness Table", "", "| Input / Setting | Given Value | Correct or Incorrect | Notes |", "|---|---:|---|---|"])
            for item in self.inputs:
                lines.append(f"| {item.get('label','')} | `{item.get('path','')}` | **{item.get('result','')}** | {item.get('message','')} |")
        if enabled("workflow"):
            lines.extend(["", "## Requested Workflow", "", "| Setting | Value |", "|---|---:|"])
            for key, value in self.workflow.items():
                lines.append(f"| {key} | `{value}` |")
        if enabled("engine_results"):
            lines.extend(["", "## Engine-Level Results", "", "| Engine | Source | Status | E(eh) | G(eh) | Imag | Elapsed |", "|---|---|---|---:|---:|---:|---:|"])
            if self.engine_results:
                for result in self.engine_results:
                    lines.append(f"| {result.get('engine','')} | `{result.get('source','')}` | {result.get('status','')} | {result.get('energy_eh','')} | {result.get('gibbs_eh','')} | {result.get('imaginary_frequencies','')} | {result.get('elapsed_seconds','')}s |")
            else:
                lines.append("| N/A | N/A | No engine jobs recorded yet | | | | |")
        if enabled("system_info"):
            lines.extend(["", "## System Info", "", f"- Platform: `{sys.platform}`", f"- Python: `{sys.version.split()[0]}`", f"- Working directory: `{os.getcwd()}`"])
        if enabled("timeline"):
            lines.extend(["", "## Timeline", "", "| Time | Level | Stage | Message |", "|---|---|---|---|"])
            for event in self.events[-100:]:
                lines.append(f"| {event.get('time','')} | {event.get('level','')} | {event.get('stage','')} | {event.get('message','')} |")
        return "\n".join(lines) + "\n"

