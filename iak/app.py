from __future__ import annotations

import argparse
import csv
import dataclasses
import json
import logging
import math
import os
import queue
import random
import re
import shutil
import ssl
import subprocess
import sys
import tarfile
import threading
import time
import traceback
import urllib.request
import webbrowser
from enum import Enum
from pathlib import Path
from tkinter import Toplevel, Text, colorchooser, filedialog, messagebox, simpledialog, ttk
from typing import Any, Callable, Dict, Iterable, List, Optional, Tuple

import numpy as np
import tkinter as tk

try:
    import matplotlib
    import matplotlib.patches as mpatches
    import matplotlib.path as mpath
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
    from matplotlib.figure import Figure
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 - registers 3D projection

    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False

from .constants import (EH2KCAL, KCAL2KJ, REACTION_TYPE_CHOICES, normalize_reaction_type)
from .colors import _C
from .plotting import (_draw_fancy_arrow, MATPLOTLIB_AVAILABLE)
from .logging_setup import setup_logging
from . import engines as _engines
from .engines import is_tool_available, inject_embedded_engines, get_wsl_path, ENGINE_DIR, XTB_URLS, CREST_URLS
from .models import Molecule, Config, RunMode
from .chemistry import (
    atom_counts_to_formula,
    validate_charge_multiplicity,
    atom_counts_for_molecule,
    expected_cluster_atom_counts,
    validate_reactant_product_atom_balance,
)
from .utils import _fmt_duration, _safe_float, _safe_int
from .pipeline import Pipeline


class IAKApp:
    def __init__(self):
        self.root = tk.Tk()
        self.root.title("IAK v11.2 The Job-Report and Visual Analytics Edition")
        self.root.geometry("1500x1040")
        self.root.configure(bg=_C["bg"])
        self._q = queue.Queue()
        self._vars = {k: tk.StringVar(value=v) for k, v in {"a": "", "b": "", "ratio": "1:0, 1:1, 1:2", "mode": "balanced", "out": "run", "reaction_type": "Non-covalent"}.items()}
        self._vars_pes = {
            "scan_mode": tk.StringVar(value="coord"),
            "reactant": tk.StringVar(value=""),
            "product": tk.StringVar(value=""),
            "engine": tk.StringVar(value="xtb"),
            "c1_a1": tk.StringVar(value="1"), "c1_a2": tk.StringVar(value="2"),
            "c1_start": tk.StringVar(value="1.0"), "c1_end": tk.StringVar(value="2.5"), "c1_steps": tk.StringVar(value="10"),
            "use_c2": tk.BooleanVar(value=False),
            "c2_a1": tk.StringVar(value="3"), "c2_a2": tk.StringVar(value="4"),
            "c2_start": tk.StringVar(value="1.5"), "c2_end": tk.StringVar(value="3.0"), "c2_steps": tk.StringVar(value="10"),
            "run_ts": tk.BooleanVar(value=True),
        }
        self._guest_list: list[str] = []  # Multi-guest B paths
        for k, v in {
            "cores": "4",
            "maxcore": "2000",
            "charge": "0",
            "mult": "1",
            "xtb_method": "--gfn2",
            "crest_method": "--gfn2",
            "orca_method": "B97-3c Opt Freq",
            "n_generate": "200",
            "n_keep_scored": "50",
            "n_keep_clustered": "40",
            "n_run_xtb": "20",
            "n_run_crest": "5",
            "rmsd_cutoff": "0.5",
            "xtb_ewin_kcal": "5.0",
            "crest_ewin_kcal": "3.0",
            "random_seed": "42",
            "max_placement_attempts": "50",
            "e_a": "",
            "e_b": "",
        }.items():
            self._vars[k] = tk.StringVar(value=v)
        self._vars_graph = {
            "dir1": tk.StringVar(value=""),
            "name1": tk.StringVar(value="Series 1"),
            "dir2": tk.StringVar(value=""),
            "name2": tk.StringVar(value="Series 2"),
            "method": tk.StringVar(value="ORCA"),
            "metric": tk.StringVar(value="Gibbs Free Energy (ΔG)"),
            "x_labels": tk.StringVar(value=""),
            "top_space": tk.BooleanVar(value=False),
            "conn_style": tk.StringVar(value="Smooth (Bezier)"),
            "delta_arrow": tk.StringVar(value="Straight"),
            "label_box": tk.StringVar(value="Rounded Box"),
            "series1_color": tk.StringVar(value="#17a2b8"),
            "series2_color": tk.StringVar(value="#ffc107"),
            "crest_color": tk.StringVar(value="#6c757d"),
            "orca_color": tk.StringVar(value="#388bfd"),
            "graph_bg": tk.StringVar(value="#ffffff"),
            "axis_color": tk.StringVar(value="#000000"),
            "grid_color": tk.StringVar(value="#b8b8b8"),
            "watermark_color": tk.StringVar(value="#6e6e6e"),
            "crop_xmin": tk.StringVar(value=""),
            "crop_xmax": tk.StringVar(value=""),
            "crop_ymin": tk.StringVar(value=""),
            "crop_ymax": tk.StringVar(value=""),
            "model_style": tk.StringVar(value="Ball-and-stick"),
            "model_bg": tk.StringVar(value="#ffffff"),
            "bond_color": tk.StringVar(value="#6b7280"),
            "atom_H": tk.StringVar(value="#f8fafc"),
            "atom_C": tk.StringVar(value="#4b5563"),
            "atom_N": tk.StringVar(value="#2563eb"),
            "atom_O": tk.StringVar(value="#dc2626"),
            "atom_F": tk.StringVar(value="#22c55e"),
            "atom_S": tk.StringVar(value="#facc15"),
            "atom_Cl": tk.StringVar(value="#16a34a"),
            "atom_other": tk.StringVar(value="#f97316"),
            "bond_mode": tk.StringVar(value="Covalent radii"),
            "bond_cutoff": tk.StringVar(value="1.85"),
            "bond_tolerance": tk.StringVar(value="0.45"),
            "show_bond_distances": tk.BooleanVar(value=False),
            "show_atom_labels": tk.BooleanVar(value=False),
            "show_3d_axes": tk.BooleanVar(value=True),
            "model_elev": tk.StringVar(value="18"),
            "model_azim": tk.StringVar(value="38"),
            "image_dpi": tk.StringVar(value="300"),
            "conformer_top_n": tk.StringVar(value="10"),
        }
        self._report_vars = {
            "inputs": tk.BooleanVar(value=True),
            "workflow": tk.BooleanVar(value=True),
            "engine_results": tk.BooleanVar(value=True),
            "thermodynamics": tk.BooleanVar(value=True),
            "images": tk.BooleanVar(value=True),
            "timeline": tk.BooleanVar(value=True),
            "system_info": tk.BooleanVar(value=True),
            "ai_prompt": tk.BooleanVar(value=True),
        }
        self.active_xtb, self.active_crest, self.active_orca = 0, 0, 0
        self.preview_tw = None
        self.preview_file = None
        self.is_running = False
        self.start_time = 0
        self._color_buttons = {}
        self._build_ui()
        setup_logging(gui_queue=self._q)
        self.root.after(100, self._poll_log)
        self.root.after(800, self._show_capabilities_popup)

    def mainloop(self):
        self.root.mainloop()

    def _show_capabilities_popup(self):
        self.cap_w = tk.Toplevel(self.root)
        self.cap_w.title("IAK Workflow Capabilities")
        self.root.update_idletasks()
        x = self.root.winfo_x() + (self.root.winfo_width() // 2) - 360
        y = self.root.winfo_y() + (self.root.winfo_height() // 2) - 270
        self.cap_w.geometry(f"720x520+{x}+{y}")
        self.cap_w.configure(bg=_C["panel"])
        self.cap_w.grab_set()
        self.cap_w.transient(self.root)
        tk.Label(self.cap_w, text="IAK v11.2 Workflow Capabilities", font=("Segoe UI", 16, "bold"), bg=_C["panel"], fg=_C["accent"]).pack(pady=(20, 10))
        caps = [
            "1. Intelligent H-bond logic docking for supramolecular assembly.",
            "2. Single Molecule Mode for isolated anchor runs.",
            "3. RMSD-based duplicate filtering and steric screening.",
            "4. xTB, CREST, and self-healing ORCA refinement pipeline.",
            "5. Per-job summary reports with input correctness checks.",
            "6. Live percent-complete and ETA progress monitor.",
            "7. JACS-grade graphing with custom colors, pan, zoom, and crop.",
            "8. Chemcraft-like molecular image styles and atom/bond color control.",
        ]
        for cap in caps:
            tk.Label(self.cap_w, text=cap, font=("Segoe UI", 11), bg=_C["panel"], fg=_C["text"], justify="left", anchor="w").pack(fill="x", padx=40, pady=8)

        def next_step():
            self.cap_w.destroy()
            self.root.after(100, self._show_guidelines_popup)

        tk.Button(self.cap_w, text="Next: Usage Guidelines", command=next_step, bg=_C["run"], fg="white", font=("Segoe UI", 11, "bold"), relief="flat", cursor="hand2", width=25).pack(pady=25)

    def _show_guidelines_popup(self):
        self.guide_w = tk.Toplevel(self.root)
        self.guide_w.title("IAK Usage Guidelines")
        self.root.update_idletasks()
        x = self.root.winfo_x() + (self.root.winfo_width() // 2) - 360
        y = self.root.winfo_y() + (self.root.winfo_height() // 2) - 300
        self.guide_w.geometry(f"720x600+{x}+{y}")
        self.guide_w.configure(bg=_C["panel"])
        self.guide_w.grab_set()
        self.guide_w.transient(self.root)
        tk.Label(self.guide_w, text="How to Use IAK", font=("Segoe UI", 16, "bold"), bg=_C["panel"], fg=_C["accent"]).pack(pady=(20, 5))
        guide_text = (
            "1. Select your Anchor (A) and Guest (B) .xyz files. Leave B empty for single molecule.\n"
            "2. Define Ratios such as 1:0, 1:1, 1:4. Use commas for batch queue.\n"
            "3. Set hardware resources and charge/multiplicity.\n"
            "4. Run the pipeline. Job reports are written into each ratio folder.\n"
            "5. Use Trend Analysis & Graphs for custom JACS-style plots and model images."
        )
        tk.Label(self.guide_w, text=guide_text, font=("Segoe UI", 11), bg=_C["panel"], fg=_C["text"], justify="left").pack(padx=40, pady=10)
        tk.Label(self.guide_w, text="Manual Engine Downloads:", font=("Segoe UI", 12, "bold"), bg=_C["panel"], fg=_C["yellow"]).pack(pady=(10, 5))

        def make_link(parent, text, url):
            lbl = tk.Label(parent, text=text, font=("Segoe UI", 10, "underline"), bg=_C["panel"], fg=_C["accent"], cursor="hand2")
            lbl.pack(pady=2)
            lbl.bind("<Button-1>", lambda e: webbrowser.open_new(url))

        make_link(self.guide_w, "xTB Releases (grimme-lab)", "https://github.com/grimme-lab/xtb/releases")
        make_link(self.guide_w, "CREST Releases (grimme-lab)", "https://github.com/grimme-lab/crest/releases")
        make_link(self.guide_w, "ORCA Forum", "https://orcaforum.kofo.mpg.de")
        tk.Label(self.guide_w, text="'Where molecules become meaning, and computation becomes discovery — inspired by Dr. Imran A. Khan,\nbrought to life by his team Asif Raza and Huma Basheer,\ncreating a platform where chemistry meets intelligence, innovation, and scientific excellence.'", font=("Segoe UI", 10, "bold", "italic"), bg=_C["panel"], fg="green", justify="center").pack(pady=30)
        tk.Button(self.guide_w, text="Enter Workspace", command=lambda: (self.guide_w.destroy(), self.root.after(100, self._show_startup_check)), bg=_C["run"], fg="white", font=("Segoe UI", 11, "bold"), relief="flat", cursor="hand2", width=20).pack(pady=10)

    def _show_startup_check(self):
        self.sw = tk.Toplevel(self.root)
        self.sw.title("System Verification")
        self.root.update_idletasks()
        x = self.root.winfo_x() + (self.root.winfo_width() // 2) - 225
        y = self.root.winfo_y() + (self.root.winfo_height() // 2) - 150
        self.sw.geometry(f"450x300+{x}+{y}")
        self.sw.configure(bg=_C["panel"])
        self.sw.grab_set()
        self.sw.transient(self.root)
        tk.Label(self.sw, text="Checking Computational Engines...", font=("Segoe UI", 14, "bold"), bg=_C["panel"], fg=_C["accent"]).pack(pady=(20, 10))
        self.check_vars = {}
        self.check_btns = {}
        frame = tk.Frame(self.sw, bg=_C["panel"])
        frame.pack(fill="both", expand=True, padx=30)
        for engine in ["xTB", "CREST", "ORCA"]:
            row = tk.Frame(frame, bg=_C["panel"])
            row.pack(fill="x", pady=8)
            tk.Label(row, text=f"{engine}:", font=("Segoe UI", 11, "bold"), width=8, anchor="w", bg=_C["panel"], fg=_C["text"]).pack(side="left")
            status_lbl = tk.Label(row, text="Checking...", font=("Segoe UI", 11, "italic"), width=15, anchor="w", bg=_C["panel"], fg=_C["muted"])
            status_lbl.pack(side="left")
            btn = tk.Button(row, text="Add Manually", command=lambda e=engine: self._load_local_from_startup(e), bg=_C["yellow"], fg="black", font=("Segoe UI", 9, "bold"), relief="flat", cursor="hand2")
            btn.pack(side="right")
            btn.pack_forget()
            self.check_vars[engine] = status_lbl
            self.check_btns[engine] = btn
        self.sw_continue = tk.Button(self.sw, text="Continue to IAK", command=self.sw.destroy, bg=_C["run"], fg="white", font=("Segoe UI", 11, "bold"), relief="flat", cursor="hand2", width=20)
        self.sw_continue.pack(pady=20)
        self.sw_continue.pack_forget()
        threading.Thread(target=self._perform_startup_checks, daemon=True).start()

    def _perform_startup_checks(self):
        time.sleep(0.5)
        all_ready = True
        for engine in ["xTB", "CREST", "ORCA"]:
            is_avail = is_tool_available(engine.lower())
            if not is_avail:
                all_ready = False
            self.root.after(0, lambda e=engine, a=is_avail: self._update_startup_ui(e, a))
            time.sleep(0.4)
        self.root.after(0, lambda: self._finalize_startup_ui(all_ready))

    def _update_startup_ui(self, engine, is_avail):
        if not hasattr(self, "sw") or not self.sw.winfo_exists():
            return
        lbl, btn = self.check_vars.get(engine), self.check_btns.get(engine)
        if is_avail:
            lbl.config(text="Running", fg=_C["green"], font=("Segoe UI", 11, "bold"))
            btn.pack_forget()
        else:
            lbl.config(text="Missing", fg=_C["red"], font=("Segoe UI", 11, "bold"))
            btn.pack(side="right")

    def _finalize_startup_ui(self, all_ready):
        if not hasattr(self, "sw") or not self.sw.winfo_exists():
            return
        self.sw_continue.pack(pady=20)
        if all_ready:
            self.sw_continue.config(text="All Engines Ready - Start", bg=_C["green"])

    def _load_local_from_startup(self, engine):
        fp = filedialog.askopenfilename(parent=self.sw, title=f"Select {engine} Archive", filetypes=[("Archives", "*.tar.xz *.tar.gz *.tgz *.zip")])
        if fp:
            self.check_vars[engine].config(text="Extracting...", fg=_C["yellow"], font=("Segoe UI", 11, "italic"))
            self.check_btns[engine].pack_forget()
            threading.Thread(target=self._extract_local_worker, args=(fp, engine, True), daemon=True).start()

    def _choose_output_base(self):
        current = self._vars["out"].get().strip(" \"'") or "run"
        current_path = Path(current)
        initial_dir = str(current_path.parent) if current_path.parent != Path(".") else os.getcwd()
        parent = filedialog.askdirectory(title="Choose parent folder for IAK output", initialdir=initial_dir)
        if not parent:
            return
        default_name = current_path.name if current_path.name else "run"
        name = simpledialog.askstring(
            "Output Name",
            "Enter the output base name you want:",
            initialvalue=default_name,
            parent=self.root,
        )
        if not name:
            return
        clean_name = re.sub(r'[<>:"/\\\\|?*]+', "_", name.strip())
        if not clean_name:
            clean_name = "run"
        self._vars["out"].set(os.path.join(parent, clean_name))

    def _apply_mode_preset(self):
        mode = self._vars["mode"].get().strip().lower()
        if mode == "custom":
            return
        try:
            cfg = Config.from_mode(RunMode[mode.upper()])
        except Exception:
            cfg = Config.from_mode(RunMode.BALANCED)
        for key in [
            "n_generate",
            "n_keep_scored",
            "n_keep_clustered",
            "n_run_xtb",
            "n_run_crest",
            "rmsd_cutoff",
            "xtb_ewin_kcal",
            "crest_ewin_kcal",
            "random_seed",
            "max_placement_attempts",
            "xtb_method",
            "crest_method",
            "orca_method",
        ]:
            if key in self._vars:
                self._vars[key].set(str(getattr(cfg, key)))

    def _poll_log(self):
        while not self._q.empty():
            self._append_text(self._q.get() + "\n")
        self.root.after(100, self._poll_log)

    def _append_text(self, t):
        self.term.config(state="normal")
        self.term.insert("end", t)
        self.term.see("end")
        self.term.config(state="disabled")

    def _status_cb(self, mode, delta):
        if mode == "xtb":
            self.active_xtb += delta
        elif mode == "crest":
            self.active_crest += delta
        elif mode == "orca":
            self.active_orca += delta
        self.root.after(0, self._update_status_ui)

    def _update_status_ui(self):
        total = self.active_xtb + self.active_crest + self.active_orca
        color = _C["green"] if total > 0 else _C["dim"]
        icon = "ACTIVE" if total > 0 else "IDLE"
        self.live_status_lbl.config(text=f"[{icon}] JOBS: {self.active_xtb} xTB | {self.active_crest} CREST | {self.active_orca} ORCA", fg=color)

    def _update_timer(self):
        if self.is_running:
            elapsed = int(time.time() - self.start_time)
            hrs, mins, secs = elapsed // 3600, (elapsed % 3600) // 60, elapsed % 60
            self.timer_lbl.config(text=f"{hrs:02d}:{mins:02d}:{secs:02d}")
        self.root.after(1000, self._update_timer)

    def _update_installation_labels(self):
        xtb_ok, crest_ok, orca_ok = is_tool_available("xtb"), is_tool_available("crest"), is_tool_available("orca")
        self._xtb_lbl.config(text="xTB: Ready" if xtb_ok else "xTB: Missing", fg=_C["green"] if xtb_ok else _C["red"])
        self._crest_lbl.config(text="CREST: Ready" if crest_ok else "CREST: Missing", fg=_C["green"] if crest_ok else _C["red"])
        self._orca_lbl.config(text="ORCA: Ready" if orca_ok else "ORCA: Missing", fg=_C["green"] if orca_ok else _C["red"])

    def _build_ui(self):
        style = ttk.Style()
        style.theme_use("clam")
        style.configure("TNotebook.Tab", padding=[15, 5], font=("Segoe UI", 10, "bold"))
        hdr = tk.Frame(self.root, bg=_C["hdr_bg"], pady=10)
        hdr.pack(fill="x")
        tk.Label(hdr, text="IAK PIPELINE v11.2", font=("Segoe UI", 20, "bold"), fg=_C["accent"], bg=_C["hdr_bg"]).pack(side="left", padx=20)
        stat_f = tk.Frame(hdr, bg=_C["hdr_bg"])
        stat_f.pack(side="left", padx=20)
        self._xtb_lbl = tk.Label(stat_f, font=("Segoe UI", 10, "bold"), bg=_C["hdr_bg"])
        self._xtb_lbl.pack(side="left", padx=10)
        self._crest_lbl = tk.Label(stat_f, font=("Segoe UI", 10, "bold"), bg=_C["hdr_bg"])
        self._crest_lbl.pack(side="left", padx=10)
        self._orca_lbl = tk.Label(stat_f, font=("Segoe UI", 10, "bold"), bg=_C["hdr_bg"])
        self._orca_lbl.pack(side="left", padx=10)
        self.timer_lbl = tk.Label(hdr, text="00:00:00", font=("Consolas", 12, "bold"), fg=_C["accent"], bg=_C["hdr_bg"])
        self.timer_lbl.pack(side="right", padx=(10, 20))
        self.live_status_lbl = tk.Label(hdr, text="[IDLE] JOBS: 0 xTB | 0 CREST | 0 ORCA", font=("Segoe UI", 10, "bold"), fg=_C["dim"], bg=_C["hdr_bg"])
        self.live_status_lbl.pack(side="right", padx=10)
        self._update_installation_labels()
        self.nb = ttk.Notebook(self.root)
        self.nb.pack(fill="both", expand=True, padx=15, pady=15)
        self.tab_main = tk.Frame(self.nb, bg=_C["bg"])
        self.tab_pes = tk.Frame(self.nb, bg=_C["bg"])
        self.tab_res = tk.Frame(self.nb, bg=_C["bg"])
        self.tab_graph = tk.Frame(self.nb, bg=_C["bg"])
        self.tab_table = tk.Frame(self.nb, bg=_C["bg"])
        self.nb.add(self.tab_main, text=" Workflow Pipeline ")
        self.nb.add(self.tab_pes, text=" PES & TS Search ")
        self.nb.add(self.tab_res, text=" Generated Results ")
        if MATPLOTLIB_AVAILABLE:
            self.nb.add(self.tab_graph, text=" Trend Analysis & Graphs ")
            self.nb.add(self.tab_table, text=" Thermodynamic Table ")
        self._build_pipeline_tab(self.tab_main)
        self._build_pes_tab(self.tab_pes)
        self._build_results_tab(self.tab_res)
        if MATPLOTLIB_AVAILABLE:
            self._build_graph_tab(self.tab_graph)
            self._build_table_tab(self.tab_table)
        self.root.after(1000, self._update_timer)

    def _build_pipeline_tab(self, parent):
        main = tk.PanedWindow(parent, orient=tk.HORIZONTAL, sashwidth=6, sashrelief=tk.RAISED, bg=_C["border"], bd=0)
        main.pack(fill="both", expand=True, padx=10, pady=10)
        left_wrapper = tk.Frame(main, bg=_C["bg"])
        left_canvas = tk.Canvas(left_wrapper, bg=_C["bg"], highlightthickness=0)
        left_scroll = tk.Scrollbar(left_wrapper, orient="vertical", command=left_canvas.yview)
        left = tk.Frame(left_canvas, bg=_C["bg"])
        left_canvas.configure(yscrollcommand=left_scroll.set)
        left_canvas.pack(side="left", fill="both", expand=True)
        left_scroll.pack(side="right", fill="y")
        canvas_window = left_canvas.create_window((0, 0), window=left, anchor="nw")
        left_canvas.bind("<Configure>", lambda e: left_canvas.itemconfig(canvas_window, width=e.width))
        left.bind("<Configure>", lambda e: left_canvas.configure(scrollregion=left_canvas.bbox("all")))

        wf = tk.LabelFrame(left, text=" 1. WORKFLOW SETUP ", fg=_C["accent"], bg=_C["panel"], font=("Segoe UI", 10, "bold"), padx=15, pady=10)
        wf.pack(fill="x", pady=(0, 6))
        tk.Label(wf, text="Tip: Leave Guest list empty and set Ratio to '1:0' for single-molecule runs.",
                 bg=_C["panel"], fg=_C["muted"], font=("Segoe UI", 9, "italic")).pack(anchor="w", pady=(0, 5))

        # ── Anchor (A) ──────────────────────────────────────────────────────────
        for lbl, var in [("Anchor (A) XYZ:", "a")]:
            f = tk.Frame(wf, bg=_C["panel"])
            f.pack(fill="x", pady=2)
            tk.Label(f, text=lbl, bg=_C["panel"], fg="white", width=26, anchor="w").pack(side="left")
            tk.Entry(f, textvariable=self._vars[var], bg=_C["entry"], fg="white", bd=0).pack(side="left", fill="x", expand=True, padx=5)
            tk.Button(f, text="...", command=lambda v=var: self._vars[v].set(
                filedialog.askopenfilename(filetypes=[("XYZ", "*.xyz"), ("All", "*")])), bg=_C["accent"]).pack(side="left")

        # ── Multi-Guest (B) panel ────────────────────────────────────────────────
        gf = tk.LabelFrame(wf, text=" Guest (B) XYZ Files — Multiple Guests Supported ",
                           fg=_C["green"], bg=_C["panel"], font=("Segoe UI", 9, "bold"), padx=8, pady=6)
        gf.pack(fill="x", pady=(4, 2))

        # Listbox showing queued guests
        lb_frame = tk.Frame(gf, bg=_C["panel"])
        lb_frame.pack(fill="x")
        self._guest_lb = tk.Listbox(lb_frame, bg=_C["entry"], fg="white", selectmode="extended",
                                    height=4, font=("Consolas", 9), bd=0, highlightthickness=0,
                                    activestyle="dotbox")
        self._guest_lb.pack(side="left", fill="x", expand=True)
        lb_sb = tk.Scrollbar(lb_frame, orient="vertical", command=self._guest_lb.yview)
        lb_sb.pack(side="right", fill="y")
        self._guest_lb.config(yscrollcommand=lb_sb.set)

        # Populate listbox from _guest_list on startup
        def _refresh_guest_lb():
            self._guest_lb.delete(0, "end")
            for p in self._guest_list:
                self._guest_lb.insert("end", os.path.basename(p))
        _refresh_guest_lb()

        # Buttons: Add / Remove / Clear
        btn_row = tk.Frame(gf, bg=_C["panel"])
        btn_row.pack(fill="x", pady=(4, 0))

        def _add_guests():
            paths = filedialog.askopenfilenames(
                title="Select Guest (B) XYZ files",
                filetypes=[("XYZ files", "*.xyz"), ("All files", "*")])
            for p in paths:
                p = os.path.abspath(p)
                if p not in self._guest_list:
                    self._guest_list.append(p)
            # also update legacy _vars["b"] to first guest for backward compat
            if self._guest_list:
                self._vars["b"].set(self._guest_list[0])
            _refresh_guest_lb()
            self._update_guest_count_label()

        def _remove_selected_guests():
            sel = list(self._guest_lb.curselection())
            for idx in reversed(sel):
                self._guest_list.pop(idx)
            if self._guest_list:
                self._vars["b"].set(self._guest_list[0])
            else:
                self._vars["b"].set("")
            _refresh_guest_lb()
            self._update_guest_count_label()

        def _clear_all_guests():
            self._guest_list.clear()
            self._vars["b"].set("")
            _refresh_guest_lb()
            self._update_guest_count_label()

        def _move_up():
            sel = list(self._guest_lb.curselection())
            if not sel or sel[0] == 0:
                return
            for idx in sel:
                self._guest_list[idx-1], self._guest_list[idx] = self._guest_list[idx], self._guest_list[idx-1]
            _refresh_guest_lb()
            for idx in sel:
                self._guest_lb.selection_set(idx-1)

        def _move_down():
            sel = list(self._guest_lb.curselection())
            if not sel or sel[-1] == len(self._guest_list)-1:
                return
            for idx in reversed(sel):
                self._guest_list[idx], self._guest_list[idx+1] = self._guest_list[idx+1], self._guest_list[idx]
            _refresh_guest_lb()
            for idx in sel:
                self._guest_lb.selection_set(idx+1)

        for txt, cmd, col in [
            ("+ Add Guest(s)",  _add_guests,            _C["green"]),
            ("- Remove Sel.",   _remove_selected_guests, _C["red"]),
            ("Clear All",       _clear_all_guests,       _C["dim"]),
            ("▲ Up",            _move_up,                _C["accent"]),
            ("▼ Down",          _move_down,              _C["accent"]),
        ]:
            tk.Button(btn_row, text=txt, command=cmd, bg=col, fg="white",
                      font=("Segoe UI", 8, "bold"), relief="flat", cursor="hand2",
                      padx=8, pady=2).pack(side="left", padx=2)

        self._guest_count_lbl = tk.Label(btn_row,
            text="0 guest(s) queued", bg=_C["panel"], fg=_C["muted"],
            font=("Segoe UI", 8, "italic"))
        self._guest_count_lbl.pack(side="left", padx=12)

        # ── Ratios / Output ──────────────────────────────────────────────────────
        for lbl, var in [("Reaction Type:", "reaction_type"), ("Ratios (e.g., 1:0, 1:1, 2:1):", "ratio"), ("Base Output Name:", "out")]:
            f = tk.Frame(wf, bg=_C["panel"])
            f.pack(fill="x", pady=2)
            tk.Label(f, text=lbl, bg=_C["panel"], fg="white", width=26, anchor="w").pack(side="left")
            if var == "reaction_type":
                cb = ttk.Combobox(f, textvariable=self._vars[var], values=REACTION_TYPE_CHOICES, state="readonly")
                cb.pack(side="left", fill="x", expand=True, padx=5)
            else:
                tk.Entry(f, textvariable=self._vars[var], bg=_C["entry"], fg="white", bd=0).pack(side="left", fill="x", expand=True, padx=5)
            if var == "out":
                tk.Button(f, text="...", command=self._choose_output_base, bg=_C["accent"]).pack(side="left")

        mode_row = tk.Frame(wf, bg=_C["panel"])
        mode_row.pack(fill="x", pady=2)
        tk.Label(mode_row, text="Preset Mode:", bg=_C["panel"], fg="white", width=26, anchor="w").pack(side="left")
        mode_cb = ttk.Combobox(mode_row, textvariable=self._vars["mode"], values=["fast", "balanced", "thorough", "custom"], width=16, state="readonly")
        mode_cb.pack(side="left", padx=5)
        mode_cb.bind("<<ComboboxSelected>>", lambda _e: self._apply_mode_preset())
        tk.Label(mode_row, text="Custom values below always override the preset.", bg=_C["panel"], fg=_C["muted"], font=("Segoe UI", 8, "italic")).pack(side="left", padx=8)

        cf = tk.LabelFrame(left, text=" 2. CHEMISTRY SETTINGS ", fg=_C["accent"], bg=_C["panel"], font=("Segoe UI", 10, "bold"), padx=15, pady=10)
        cf.pack(fill="x", pady=(0, 6))
        for lbl, var, width in [("System Charge:", "charge", 10), ("Multiplicity:", "mult", 10), ("xTB Flags:", "xtb_method", 18), ("CREST Flags:", "crest_method", 18), ("ORCA Method:", "orca_method", 40)]:
            row = tk.Frame(cf, bg=_C["panel"])
            row.pack(fill="x", pady=2)
            tk.Label(row, text=lbl, bg=_C["panel"], fg="white", width=22, anchor="w").pack(side="left")
            tk.Entry(row, textvariable=self._vars[var], bg=_C["entry"], fg="white", width=width, bd=0).pack(side="left", fill="x", expand=var in {"orca_method", "xtb_method", "crest_method"}, padx=5)

        tf = tk.LabelFrame(left, text=" 3. THERMODYNAMICS (Optional) ", fg=_C["accent"], bg=_C["panel"], font=("Segoe UI", 10, "bold"), padx=15, pady=10)
        tf.pack(fill="x", pady=(0, 6))
        tk.Label(tf, text="Input isolated ORCA energies (Eh) for Binding Affinity.", bg=_C["panel"], fg=_C["muted"], font=("Segoe UI", 9, "italic")).pack(anchor="w", pady=(0, 5))
        for lbl, var in [("Monomer A Energy:", "e_a"), ("Monomer B Energy:", "e_b")]:
            row = tk.Frame(tf, bg=_C["panel"])
            row.pack(fill="x", pady=2)
            tk.Label(row, text=lbl, bg=_C["panel"], fg="white", width=22, anchor="w").pack(side="left")
            tk.Entry(row, textvariable=self._vars[var], bg=_C["entry"], fg="white", bd=0).pack(side="left", fill="x", expand=True, padx=5)

        hf = tk.LabelFrame(left, text=" 4. HARDWARE RESOURCES ", fg=_C["accent"], bg=_C["panel"], font=("Segoe UI", 10, "bold"), padx=15, pady=10)
        hf.pack(fill="x", pady=(0, 6))
        for lbl, var, note in [("CPU Cores:", "cores", "(16, 32+ for large systems)"), ("RAM/Core (MB):", "maxcore", "(Total = Cores * RAM)")]:
            row = tk.Frame(hf, bg=_C["panel"])
            row.pack(fill="x", pady=2)
            tk.Label(row, text=lbl, bg=_C["panel"], fg="white", width=22, anchor="w").pack(side="left")
            tk.Entry(row, textvariable=self._vars[var], bg=_C["entry"], fg="white", width=10, bd=0).pack(side="left", padx=5)
            tk.Label(row, text=note, bg=_C["panel"], fg=_C["muted"]).pack(side="left", padx=10)

        af = tk.LabelFrame(left, text=" 5. HOST-CONTROLLED SAMPLING AND FILTERS ", fg=_C["accent"], bg=_C["panel"], font=("Segoe UI", 10, "bold"), padx=15, pady=10)
        af.pack(fill="x", pady=(0, 6))
        host_rows = [
            [("Generate Seeds", "n_generate"), ("Score Keep", "n_keep_scored"), ("Cluster Keep", "n_keep_clustered")],
            [("Run xTB", "n_run_xtb"), ("Run CREST", "n_run_crest"), ("Random Seed", "random_seed")],
            [("RMSD Cutoff", "rmsd_cutoff"), ("xTB Window kcal", "xtb_ewin_kcal"), ("CREST Window kcal", "crest_ewin_kcal")],
            [("Placement Attempts", "max_placement_attempts")],
        ]
        for row_items in host_rows:
            row = tk.Frame(af, bg=_C["panel"])
            row.pack(fill="x", pady=2)
            for lbl, var in row_items:
                tk.Label(row, text=f"{lbl}:", bg=_C["panel"], fg="white", width=17, anchor="w").pack(side="left", padx=(0, 2))
                tk.Entry(row, textvariable=self._vars[var], bg=_C["entry"], fg="white", width=9, bd=0).pack(side="left", padx=(0, 10), ipady=2)

        sf = tk.LabelFrame(left, text=" 6. HOST-CONTROLLED REPORT CONTENT ", fg=_C["accent"], bg=_C["panel"], font=("Segoe UI", 10, "bold"), padx=15, pady=10)
        sf.pack(fill="x", pady=(0, 6))
        report_items = [
            ("Inputs", "inputs"),
            ("Workflow", "workflow"),
            ("Engine Results", "engine_results"),
            ("Thermodynamics", "thermodynamics"),
            ("Images", "images"),
            ("Timeline", "timeline"),
            ("System Info", "system_info"),
            ("AI Prompt", "ai_prompt"),
        ]
        for label, key in report_items:
            tk.Checkbutton(sf, text=label, variable=self._report_vars[key], bg=_C["panel"], fg="white", selectcolor="black").pack(side="left", padx=4)

        rf = tk.Frame(left, bg=_C["panel"])
        rf.pack(fill="x", pady=15)
        self.run_preopt = tk.BooleanVar(value=True)
        self.run_xtb = tk.BooleanVar(value=True)
        self.run_crest = tk.BooleanVar(value=True)
        self.run_orca = tk.BooleanVar(value=True)
        for text, var, fg in [("Pre-Opt", self.run_preopt, _C["yellow"]), ("Run xTB", self.run_xtb, "white"), ("Run CREST", self.run_crest, "white"), ("Run ORCA DFT", self.run_orca, _C["accent"])]:
            tk.Checkbutton(rf, text=text, variable=var, bg=_C["panel"], fg=fg, selectcolor="black").pack(side="left", padx=5)
        btn_f = tk.Frame(left, bg=_C["panel"])
        btn_f.pack(fill="x", pady=5)
        tk.Button(btn_f, text="LOAD LOCAL ENGINE (.tar.xz / .zip)", command=self._load_local, bg="#d29922", fg="white", font=("Segoe UI", 9, "bold")).pack(side="left", fill="x", expand=True, padx=(2, 0))
        tk.Button(btn_f, text="ANALYZE FOLDER", command=self._analyze_existing_folder, bg="#6e40c9", fg="white", font=("Segoe UI", 9, "bold")).pack(side="left", fill="x", expand=True, padx=(6, 2))
        tk.Button(btn_f, text="ANALYZE FILES (.out/.xyz)", command=self._analyze_existing_files, bg="#8957e5", fg="white", font=("Segoe UI", 9, "bold")).pack(side="left", fill="x", expand=True, padx=(6, 2))
        tk.Button(btn_f, text="RESUME STOPPED JOB", command=self._resume_existing_job, bg="#238636", fg="white", font=("Segoe UI", 9, "bold")).pack(side="left", fill="x", expand=True, padx=(6, 2))
        self.go_btn = tk.Button(left, text="START BATCH PIPELINE", command=self._start, bg=_C["run"], fg="white", font=("Segoe UI", 12, "bold"), pady=10)
        self.go_btn.pack(fill="x", pady=10)

        right = tk.LabelFrame(main, text=" LIVE PIPELINE TERMINAL ", fg=_C["accent"], bg=_C["panel"], font=("Segoe UI", 10, "bold"))
        self.term = tk.Text(right, bg="black", fg=_C["green"], font=("Consolas", 10), state="disabled")
        self.term.pack(fill="both", expand=True, padx=5, pady=5)
        main.add(left_wrapper, minsize=560, stretch="always")
        main.add(right, minsize=460, stretch="always")
        self._build_progress_monitor(parent)

    def _build_progress_monitor(self, parent):
        frame = tk.LabelFrame(parent, text=" HIGH-LEVEL JOB PROGRESS, PERCENT COMPLETE, AND ETA ", fg=_C["accent"], bg=_C["panel"], font=("Segoe UI", 10, "bold"), padx=10, pady=8)
        frame.pack(fill="x", padx=10, pady=(0, 10))
        top = tk.Frame(frame, bg=_C["panel"])
        top.pack(fill="x")
        self.iak_progress_stage_lbl = tk.Label(top, text="Waiting for job", bg=_C["panel"], fg=_C["text"], font=("Segoe UI", 10, "bold"), anchor="w")
        self.iak_progress_stage_lbl.pack(side="left", fill="x", expand=True)
        self.iak_progress_eta_lbl = tk.Label(top, text="ETA: --:--:--", bg=_C["panel"], fg=_C["muted"], font=("Consolas", 10, "bold"), width=18)
        self.iak_progress_eta_lbl.pack(side="right", padx=6)
        self.iak_progress_percent_lbl = tk.Label(top, text="0.0%", bg=_C["panel"], fg=_C["green"], font=("Consolas", 11, "bold"), width=8)
        self.iak_progress_percent_lbl.pack(side="right", padx=6)
        self.iak_progress_bar = ttk.Progressbar(frame, orient="horizontal", mode="determinate", maximum=100)
        self.iak_progress_bar.pack(fill="x", pady=(6, 8))
        cols = ("stage", "status", "percent", "eta")
        self.iak_progress_tree = ttk.Treeview(frame, columns=cols, show="headings", height=4)
        for col, text, width in [("stage", "Pipeline Stage", 220), ("status", "Status", 110), ("percent", "Complete", 90), ("eta", "ETA", 120)]:
            self.iak_progress_tree.heading(col, text=text)
            self.iak_progress_tree.column(col, width=width, anchor="center")
        self.iak_progress_tree.pack(fill="x")
        self._progress_items = {}
        self._reset_progress_ui()

    def _reset_progress_ui(self):
        if not hasattr(self, "iak_progress_tree"):
            return
        for row in self.iak_progress_tree.get_children():
            self.iak_progress_tree.delete(row)
        self._progress_items = {}
        for stage in ["Validation", "Sampling/filtering", "xTB", "CREST", "ORCA", "Reports"]:
            item = self.iak_progress_tree.insert("", "end", values=(stage, "Waiting", "0.0%", "--:--:--"))
            self._progress_items[stage] = item
        self.iak_progress_bar["value"] = 0
        self.iak_progress_stage_lbl.config(text="Waiting for job")
        self.iak_progress_percent_lbl.config(text="0.0%")
        self.iak_progress_eta_lbl.config(text="ETA: --:--:--")

    def _progress_cb(self, payload):
        self.root.after(0, lambda p=dict(payload): self._apply_progress_payload(p))

    def _normalize_stage(self, stage):
        text = (stage or "").lower()
        if "validation" in text:
            return "Validation"
        if "sample" in text or "filter" in text or "generat" in text:
            return "Sampling/filtering"
        if "xtb" in text:
            return "xTB"
        if "crest" in text:
            return "CREST"
        if "orca" in text:
            return "ORCA"
        if "report" in text or "complete" in text:
            return "Reports"
        return stage or "Pipeline"

    def _apply_progress_payload(self, payload):
        percent = max(0.0, min(100.0, float(payload.get("percent", 0.0))))
        elapsed = float(payload.get("elapsed_seconds", 0.0))
        eta = "complete"
        if 0.0 < percent < 100.0:
            eta = _fmt_duration(elapsed * (100.0 - percent) / percent)
        elif percent <= 0.0:
            eta = "estimating"
        stage = self._normalize_stage(payload.get("stage", "Pipeline"))
        status = payload.get("status", "running").title()
        message = payload.get("message") or stage
        job = payload.get("job", "")
        self.iak_progress_bar["value"] = percent
        self.iak_progress_percent_lbl.config(text=f"{percent:5.1f}%")
        self.iak_progress_eta_lbl.config(text=f"ETA: {eta}")
        self.iak_progress_stage_lbl.config(text=f"Ratio {job} | {stage}: {message}" if job else f"{stage}: {message}")
        item = self._progress_items.get(stage)
        if not item:
            item = self.iak_progress_tree.insert("", "end", values=(stage, "Waiting", "0.0%", "--:--:--"))
            self._progress_items[stage] = item
        self.iak_progress_tree.item(item, values=(stage, status, f"{percent:5.1f}%", eta))
        self.iak_progress_tree.see(item)

    def _get_report_options(self):
        return {key: bool(var.get()) for key, var in self._report_vars.items()}

    def _build_pes_tab(self, parent):
        left = tk.Frame(parent, bg=_C["bg"])
        left.pack(side="left", fill="both", expand=True, padx=10, pady=10)
        right = tk.Frame(parent, bg=_C["panel"], width=600)
        right.pack(side="right", fill="both", expand=False, padx=10, pady=10)
        
        # --- Left Panel: Settings ---
        lbl_cfg = tk.LabelFrame(left, text=" 1. PES & TS CONFIGURATION ", fg=_C["accent"], bg=_C["panel"], font=("Segoe UI", 10, "bold"), padx=15, pady=10)
        lbl_cfg.pack(fill="x", pady=5)
        
        # Scan Mode
        f_mode = tk.Frame(lbl_cfg, bg=_C["panel"])
        f_mode.pack(fill="x", pady=4)
        tk.Label(f_mode, text="Scan Mode:", bg=_C["panel"], fg="white", width=20, anchor="w").pack(side="left")
        tk.Radiobutton(f_mode, text="Coordinate Grid (1D/2D)", variable=self._vars_pes["scan_mode"], value="coord", bg=_C["panel"], fg="white", selectcolor=_C["bg"]).pack(side="left")
        tk.Radiobutton(f_mode, text="Reactant -> Product Path", variable=self._vars_pes["scan_mode"], value="path", bg=_C["panel"], fg="white", selectcolor=_C["bg"]).pack(side="left")

        # Reactant XYZ
        f_react = tk.Frame(lbl_cfg, bg=_C["panel"])
        f_react.pack(fill="x", pady=4)
        tk.Label(f_react, text="Reactant XYZ:", bg=_C["panel"], fg="white", width=20, anchor="w").pack(side="left")
        tk.Entry(f_react, textvariable=self._vars_pes["reactant"], bg=_C["entry"], fg="white", bd=0).pack(side="left", fill="x", expand=True, padx=5)
        tk.Button(f_react, text="...", command=lambda: self._vars_pes["reactant"].set(filedialog.askopenfilename()), bg=_C["accent"]).pack(side="left")
        
        # Product XYZ (For Path Mode)
        f_prod = tk.Frame(lbl_cfg, bg=_C["panel"])
        f_prod.pack(fill="x", pady=4)
        tk.Label(f_prod, text="Product XYZ (Path Mode):", bg=_C["panel"], fg="white", width=20, anchor="w").pack(side="left")
        tk.Entry(f_prod, textvariable=self._vars_pes["product"], bg=_C["entry"], fg="white", bd=0).pack(side="left", fill="x", expand=True, padx=5)
        tk.Button(f_prod, text="...", command=lambda: self._vars_pes["product"].set(filedialog.askopenfilename()), bg=_C["accent"]).pack(side="left")
        
        # Engine
        f_eng = tk.Frame(lbl_cfg, bg=_C["panel"])
        f_eng.pack(fill="x", pady=4)
        tk.Label(f_eng, text="Engine:", bg=_C["panel"], fg="white", width=20, anchor="w").pack(side="left")
        ttk.Combobox(f_eng, textvariable=self._vars_pes["engine"], values=["xtb", "orca"], state="readonly", width=10).pack(side="left", padx=5)
        
        # Coord 1
        lbl_c1 = tk.LabelFrame(left, text=" Coordinate 1 (Required) ", fg=_C["green"], bg=_C["panel"], font=("Segoe UI", 9, "bold"), padx=10, pady=5)
        lbl_c1.pack(fill="x", pady=5)
        f_c1_a = tk.Frame(lbl_c1, bg=_C["panel"])
        f_c1_a.pack(fill="x", pady=2)
        tk.Label(f_c1_a, text="Atom 1 Index:", bg=_C["panel"], fg="white", width=15).pack(side="left")
        tk.Entry(f_c1_a, textvariable=self._vars_pes["c1_a1"], bg=_C["entry"], fg="white", width=5).pack(side="left")
        tk.Label(f_c1_a, text="Atom 2 Index:", bg=_C["panel"], fg="white", width=15).pack(side="left")
        tk.Entry(f_c1_a, textvariable=self._vars_pes["c1_a2"], bg=_C["entry"], fg="white", width=5).pack(side="left")
        
        f_c1_d = tk.Frame(lbl_c1, bg=_C["panel"])
        f_c1_d.pack(fill="x", pady=2)
        tk.Label(f_c1_d, text="Start Dist (Å):", bg=_C["panel"], fg="white", width=15).pack(side="left")
        tk.Entry(f_c1_d, textvariable=self._vars_pes["c1_start"], bg=_C["entry"], fg="white", width=5).pack(side="left")
        tk.Label(f_c1_d, text="End Dist (Å):", bg=_C["panel"], fg="white", width=15).pack(side="left")
        tk.Entry(f_c1_d, textvariable=self._vars_pes["c1_end"], bg=_C["entry"], fg="white", width=5).pack(side="left")
        tk.Label(f_c1_d, text="Steps:", bg=_C["panel"], fg="white", width=10).pack(side="left")
        tk.Entry(f_c1_d, textvariable=self._vars_pes["c1_steps"], bg=_C["entry"], fg="white", width=5).pack(side="left")
        
        # Coord 2
        lbl_c2 = tk.LabelFrame(left, text=" Coordinate 2 (For 3D Plot) ", fg=_C["yellow"], bg=_C["panel"], font=("Segoe UI", 9, "bold"), padx=10, pady=5)
        lbl_c2.pack(fill="x", pady=5)
        tk.Checkbutton(lbl_c2, text="Enable 2D Grid Scan", variable=self._vars_pes["use_c2"], bg=_C["panel"], fg="white", selectcolor=_C["bg"]).pack(anchor="w", pady=2)
        f_c2_a = tk.Frame(lbl_c2, bg=_C["panel"])
        f_c2_a.pack(fill="x", pady=2)
        tk.Label(f_c2_a, text="Atom 1 Index:", bg=_C["panel"], fg="white", width=15).pack(side="left")
        tk.Entry(f_c2_a, textvariable=self._vars_pes["c2_a1"], bg=_C["entry"], fg="white", width=5).pack(side="left")
        tk.Label(f_c2_a, text="Atom 2 Index:", bg=_C["panel"], fg="white", width=15).pack(side="left")
        tk.Entry(f_c2_a, textvariable=self._vars_pes["c2_a2"], bg=_C["entry"], fg="white", width=5).pack(side="left")
        
        f_c2_d = tk.Frame(lbl_c2, bg=_C["panel"])
        f_c2_d.pack(fill="x", pady=2)
        tk.Label(f_c2_d, text="Start Dist (Å):", bg=_C["panel"], fg="white", width=15).pack(side="left")
        tk.Entry(f_c2_d, textvariable=self._vars_pes["c2_start"], bg=_C["entry"], fg="white", width=5).pack(side="left")
        tk.Label(f_c2_d, text="End Dist (Å):", bg=_C["panel"], fg="white", width=15).pack(side="left")
        tk.Entry(f_c2_d, textvariable=self._vars_pes["c2_end"], bg=_C["entry"], fg="white", width=5).pack(side="left")
        tk.Label(f_c2_d, text="Steps:", bg=_C["panel"], fg="white", width=10).pack(side="left")
        tk.Entry(f_c2_d, textvariable=self._vars_pes["c2_steps"], bg=_C["entry"], fg="white", width=5).pack(side="left")
        
        # TS Options
        lbl_ts = tk.LabelFrame(left, text=" 2. TRANSITION STATE REFINEMENT ", fg=_C["red"], bg=_C["panel"], font=("Segoe UI", 10, "bold"), padx=15, pady=10)
        lbl_ts.pack(fill="x", pady=5)
        tk.Checkbutton(lbl_ts, text="Find Saddle Point & Run TS Optimization (OptTS + Freq)", variable=self._vars_pes["run_ts"], bg=_C["panel"], fg="white", selectcolor=_C["bg"]).pack(anchor="w", pady=2)
        
        # Run Button
        self.btn_run_pes = tk.Button(left, text="START PES & TS SEARCH", command=self._start_pes, bg=_C["run"], fg="white", font=("Segoe UI", 14, "bold"), pady=10)
        self.btn_run_pes.pack(fill="x", pady=15)
        
        # --- Right Panel: Terminal ---
        tk.Label(right, text=" PES/TS TERMINAL LOG ", bg=_C["panel"], fg="white", font=("Consolas", 10, "bold")).pack(pady=5)
        self.pes_term = tk.Text(right, bg=_C["entry"], fg="#00ff00", font=("Consolas", 9), state="disabled", wrap="word")
        self.pes_term.pack(fill="both", expand=True, padx=5, pady=5)

    def _build_results_tab(self, parent):
        top = tk.Frame(parent, bg=_C["panel"], pady=10)
        top.pack(fill="x")
        tk.Button(top, text=" REFRESH ", command=self._refresh_results, bg=_C["accent"], fg="white", font=("Segoe UI", 10, "bold"), padx=10).pack(side="left", padx=10)
        tk.Button(top, text=" ANALYZE FOLDER ", command=self._analyze_existing_folder, bg="#6e40c9", fg="white", font=("Segoe UI", 10, "bold"), padx=10).pack(side="left", padx=10)
        tk.Button(top, text=" ANALYZE FILES (.out/.xyz) ", command=self._analyze_existing_files, bg="#8957e5", fg="white", font=("Segoe UI", 10, "bold"), padx=10).pack(side="left", padx=10)
        tk.Button(top, text=" AVOGADRO ", command=self._open_in_avogadro, bg="#6e40c9", fg="white", font=("Segoe UI", 10, "bold"), padx=10).pack(side="left", padx=10)
        tk.Button(top, text=" EXPORT ", command=self._export_file, bg=_C["green"], fg="white", font=("Segoe UI", 10, "bold"), padx=10).pack(side="left", padx=10)
        tk.Label(top, text="Double-click XYZ for interactive model viewer; hover for quick preview.", bg=_C["panel"], fg=_C["muted"], font=("Segoe UI", 10, "italic")).pack(side="right", padx=20)
        container = tk.Frame(parent, bg=_C["bg"])
        container.pack(fill="both", expand=True, padx=10, pady=10)
        for i in range(6):
            container.columnconfigure(i, weight=1)
        container.rowconfigure(0, weight=1)
        self.listboxes = []
        self.folders = ["01_Inputs_and_Clusters", "02_xTB_Results", "03_CREST_Results", "04_ORCA_Refinement", "05_Top_Models_Comparison", "06_File_Analysis_Images"]

        def make_listbox(col, title, fg_color=_C["green"]):
            f = tk.LabelFrame(container, text=title, bg=_C["panel"], fg=_C["accent"], font=("Segoe UI", 9, "bold"))
            f.grid(row=0, column=col, sticky="nsew", padx=2, pady=5)
            lb = tk.Listbox(f, bg=_C["entry"], fg=fg_color, font=("Consolas", 9), selectbackground=_C["accent"], exportselection=False)
            lb.pack(side="left", fill="both", expand=True, padx=2, pady=5)
            sb = ttk.Scrollbar(f, orient="vertical", command=lb.yview)
            sb.pack(side="right", fill="y")
            lb.config(yscrollcommand=sb.set)
            lb.bind("<Double-Button-1>", lambda e, l=lb, c=col: self._open_file(l, c))
            lb.bind("<Motion>", lambda e, l=lb, c=col: self._on_hover(e, l, c))
            lb.bind("<Leave>", lambda e: self._hide_preview())
            self.listboxes.append(lb)

        make_listbox(0, " 1. Clusters ")
        make_listbox(1, " 2. xTB Opt ")
        make_listbox(2, " 3. CREST ")
        make_listbox(3, " 4. ORCA DFT ")
        make_listbox(4, " 5. TOP 3 COMPARE ", fg_color=_C["yellow"])
        make_listbox(5, " 6. 3D IMAGES ", fg_color=_C["accent"])

    def _build_graph_tab(self, parent):
        top = tk.Frame(parent, bg=_C["panel"], pady=10)
        top.pack(fill="x")
        r1 = tk.Frame(top, bg=_C["panel"])
        r1.pack(fill="x", pady=2)
        tk.Label(r1, text="Series 1 Dir:", bg=_C["panel"], fg="white", font=("Segoe UI", 9, "bold"), width=12, anchor="w").pack(side="left", padx=(10, 0))
        tk.Entry(r1, textvariable=self._vars_graph["dir1"], bg=_C["entry"], fg="white", width=25, bd=0).pack(side="left", padx=5, ipady=3)
        tk.Button(r1, text="...", command=lambda: self._vars_graph["dir1"].set(filedialog.askdirectory()), bg=_C["accent"], fg="white", bd=0, padx=5).pack(side="left")
        tk.Label(r1, text="Legend 1:", bg=_C["panel"], fg="white", font=("Segoe UI", 9, "bold")).pack(side="left", padx=(15, 5))
        tk.Entry(r1, textvariable=self._vars_graph["name1"], bg=_C["entry"], fg="white", width=15, bd=0).pack(side="left", padx=5, ipady=3)
        tk.Label(r1, text="Series 2 Dir:", bg=_C["panel"], fg="white", font=("Segoe UI", 9, "bold"), width=12, anchor="w").pack(side="left", padx=(25, 0))
        tk.Entry(r1, textvariable=self._vars_graph["dir2"], bg=_C["entry"], fg="white", width=25, bd=0).pack(side="left", padx=5, ipady=3)
        tk.Button(r1, text="...", command=lambda: self._vars_graph["dir2"].set(filedialog.askdirectory()), bg=_C["accent"], fg="white", bd=0, padx=5).pack(side="left")
        tk.Label(r1, text="Legend 2:", bg=_C["panel"], fg="white", font=("Segoe UI", 9, "bold")).pack(side="left", padx=(15, 5))
        tk.Entry(r1, textvariable=self._vars_graph["name2"], bg=_C["entry"], fg="white", width=15, bd=0).pack(side="left", padx=5, ipady=3)

        r2 = tk.Frame(top, bg=_C["panel"])
        r2.pack(fill="x", pady=(10, 2))
        tk.Label(r2, text="Custom X-Axis Labels:", bg=_C["panel"], fg="white", font=("Segoe UI", 9, "bold"), width=18, anchor="w").pack(side="left", padx=(10, 0))
        tk.Entry(r2, textvariable=self._vars_graph["x_labels"], bg=_C["entry"], fg="white", width=30, bd=0).pack(side="left", padx=5, ipady=3)
        tk.Checkbutton(r2, text="Add Blank Space for 3D Models", variable=self._vars_graph["top_space"], bg=_C["panel"], fg=_C["yellow"], selectcolor="black", font=("Segoe UI", 9, "bold")).pack(side="left", padx=10)

        r4 = tk.Frame(top, bg=_C["panel"])
        r4.pack(fill="x", pady=(10, 2))
        tk.Label(r4, text="Connection Line:", bg=_C["panel"], fg="white", font=("Segoe UI", 9, "bold")).pack(side="left", padx=(10, 5))
        ttk.Combobox(r4, textvariable=self._vars_graph["conn_style"], values=["Smooth (Bezier)", "Straight", "Wavy", "Zigzag", "Curled (Arc)"], width=15, state="readonly").pack(side="left", padx=5)
        tk.Label(r4, text="Delta Arrow Style:", bg=_C["panel"], fg="white", font=("Segoe UI", 9, "bold")).pack(side="left", padx=(15, 5))
        ttk.Combobox(r4, textvariable=self._vars_graph["delta_arrow"], values=["Straight", "Wavy", "Zigzag", "Curled (Arc)", "No Arrow"], width=12, state="readonly").pack(side="left", padx=5)
        tk.Label(r4, text="Label Box:", bg=_C["panel"], fg="white", font=("Segoe UI", 9, "bold")).pack(side="left", padx=(15, 5))
        ttk.Combobox(r4, textvariable=self._vars_graph["label_box"], values=["Rounded Box", "Square Box", "No Box"], width=12, state="readonly").pack(side="left", padx=5)

        r3 = tk.Frame(top, bg=_C["panel"])
        r3.pack(fill="x", pady=(15, 2))
        tk.Label(r3, text="Energy Source:", bg=_C["panel"], fg="white", font=("Segoe UI", 9, "bold")).pack(side="left", padx=(10, 5))
        ttk.Combobox(r3, textvariable=self._vars_graph["method"], values=["ORCA", "CREST", "xTB"], width=10, state="readonly").pack(side="left", padx=5)
        tk.Label(r3, text="Metric:", bg=_C["panel"], fg="white", font=("Segoe UI", 9, "bold")).pack(side="left", padx=(15, 5))
        ttk.Combobox(r3, textvariable=self._vars_graph["metric"], values=["Total Energy (ΔE)", "Gibbs Free Energy (ΔG)"], width=20, state="readonly").pack(side="left", padx=5)
        tk.Button(r3, text=" Plot Publication Energy Profile ", command=self._plot_reaction_profile, bg=_C["green"], fg="white", font=("Segoe UI", 10, "bold"), padx=15).pack(side="left", padx=(20, 10))
        tk.Button(r3, text=" Plot CREST vs ORCA Method Acc. ", command=self._plot_methods, bg=_C["accent"], fg="white", font=("Segoe UI", 10, "bold"), padx=15).pack(side="left", padx=10)
        tk.Button(r3, text=" Plot All Optimized ", command=self._plot_all_optimized, bg="#8957e5", fg="white", font=("Segoe UI", 10, "bold"), padx=12).pack(side="left", padx=8)
        tk.Button(r3, text=" Plot Top Conformers ", command=self._plot_top_conformers, bg="#6e40c9", fg="white", font=("Segoe UI", 10, "bold"), padx=12).pack(side="left", padx=8)
        tk.Label(r3, text="Top N:", bg=_C["panel"], fg="white", font=("Segoe UI", 9, "bold")).pack(side="left", padx=(8, 2))
        tk.Entry(r3, textvariable=self._vars_graph["conformer_top_n"], bg=_C["entry"], fg="white", width=5, bd=0).pack(side="left", ipady=2)
        tk.Button(r3, text=" Copy to Clipboard ", command=lambda: self._copy_to_clipboard(self.fig), bg="#6e40c9", fg="white", font=("Segoe UI", 10, "bold"), padx=10).pack(side="right", padx=20)
        tk.Button(r3, text=" Download High-Res Graph ", command=self._export_high_res_graph, bg=_C["yellow"], fg="black", font=("Segoe UI", 10, "bold"), padx=15).pack(side="right", padx=10)
        self._build_visual_controls(top)

        self.fig = Figure(figsize=(12, 7), dpi=100)
        self.ax = self.fig.add_subplot(111)
        self.ax2 = None
        self.canvas = FigureCanvasTkAgg(self.fig, master=parent)
        self.canvas.get_tk_widget().pack(fill="both", expand=True, padx=10, pady=10)
        self.graph_toolbar = NavigationToolbar2Tk(self.canvas, parent, pack_toolbar=False)
        self.graph_toolbar.update()
        self.graph_toolbar.pack(fill="x", padx=10, pady=(0, 6))

    def _build_visual_controls(self, top):
        row = tk.Frame(top, bg=_C["panel"])
        row.pack(fill="x", pady=(8, 2))
        tk.Label(row, text="JACS visual colors:", bg=_C["panel"], fg=_C["yellow"], font=("Segoe UI", 9, "bold")).pack(side="left", padx=(10, 6))
        for label, key in [("Series 1", "series1_color"), ("Series 2", "series2_color"), ("ORCA", "orca_color"), ("CREST", "crest_color"), ("BG", "graph_bg"), ("Axis", "axis_color"), ("Grid", "grid_color")]:
            self._color_button(row, label, key)
        row2 = tk.Frame(top, bg=_C["panel"])
        row2.pack(fill="x", pady=(4, 2))
        tk.Label(row2, text="Model:", bg=_C["panel"], fg=_C["text"], font=("Segoe UI", 9, "bold")).pack(side="left", padx=(10, 5))
        ttk.Combobox(row2, textvariable=self._vars_graph["model_style"], values=["Ball-and-stick", "CPK", "Wireframe", "Licorice"], width=14, state="readonly").pack(side="left", padx=3)
        for label, key in [("C", "atom_C"), ("H", "atom_H"), ("N", "atom_N"), ("O", "atom_O"), ("F", "atom_F"), ("S", "atom_S"), ("Bond", "bond_color"), ("Model BG", "model_bg")]:
            self._color_button(row2, label, key, width=5)
        row3 = tk.Frame(top, bg=_C["panel"])
        row3.pack(fill="x", pady=(4, 2))
        tk.Label(row3, text="3D bonds:", bg=_C["panel"], fg=_C["text"], font=("Segoe UI", 9, "bold")).pack(side="left", padx=(10, 5))
        ttk.Combobox(
            row3,
            textvariable=self._vars_graph["bond_mode"],
            values=["Covalent radii", "Distance cutoff", "Hydrogen bonds", "All close contacts"],
            width=16,
            state="readonly",
        ).pack(side="left", padx=4)
        tk.Label(row3, text="Cutoff A:", bg=_C["panel"], fg=_C["muted"]).pack(side="left", padx=(8, 2))
        tk.Entry(row3, textvariable=self._vars_graph["bond_cutoff"], bg=_C["entry"], fg="white", width=5, bd=0).pack(side="left", ipady=2)
        tk.Label(row3, text="Tol A:", bg=_C["panel"], fg=_C["muted"]).pack(side="left", padx=(8, 2))
        tk.Entry(row3, textvariable=self._vars_graph["bond_tolerance"], bg=_C["entry"], fg="white", width=5, bd=0).pack(side="left", ipady=2)
        tk.Checkbutton(row3, text="Distances", variable=self._vars_graph["show_bond_distances"], bg=_C["panel"], fg=_C["yellow"], selectcolor="black").pack(side="left", padx=8)
        tk.Checkbutton(row3, text="Atom labels", variable=self._vars_graph["show_atom_labels"], bg=_C["panel"], fg=_C["yellow"], selectcolor="black").pack(side="left", padx=4)
        tk.Checkbutton(row3, text="Axes", variable=self._vars_graph["show_3d_axes"], bg=_C["panel"], fg=_C["yellow"], selectcolor="black").pack(side="left", padx=4)
        for label, key in [("Elev", "model_elev"), ("Azim", "model_azim"), ("DPI", "image_dpi")]:
            tk.Label(row3, text=f"{label}:", bg=_C["panel"], fg=_C["muted"]).pack(side="left", padx=(8, 2))
            tk.Entry(row3, textvariable=self._vars_graph[key], bg=_C["entry"], fg="white", width=5, bd=0).pack(side="left", ipady=2)
        row4 = tk.Frame(top, bg=_C["panel"])
        row4.pack(fill="x", pady=(4, 2))
        tk.Label(row4, text="Crop axes:", bg=_C["panel"], fg=_C["text"], font=("Segoe UI", 9, "bold")).pack(side="left", padx=(10, 5))
        for label, key in [("X min", "crop_xmin"), ("X max", "crop_xmax"), ("Y min", "crop_ymin"), ("Y max", "crop_ymax")]:
            tk.Label(row4, text=label, bg=_C["panel"], fg=_C["muted"]).pack(side="left", padx=(5, 2))
            tk.Entry(row4, textvariable=self._vars_graph[key], bg=_C["entry"], fg="white", width=7, bd=0).pack(side="left", ipady=2)
        tk.Button(row4, text="Apply Crop", command=self._apply_graph_crop, bg=_C["accent"], fg="white", relief="flat", font=("Segoe UI", 9, "bold")).pack(side="left", padx=8)
        tk.Button(row4, text="Reset View", command=self._reset_graph_view, bg=_C["dim"], fg="white", relief="flat", font=("Segoe UI", 9, "bold")).pack(side="left", padx=4)
        tk.Button(row4, text="Pan", command=lambda: self.graph_toolbar.pan(), bg=_C["panel"], fg="white", relief="groove").pack(side="left", padx=4)
        tk.Button(row4, text="Zoom", command=lambda: self.graph_toolbar.zoom(), bg=_C["panel"], fg="white", relief="groove").pack(side="left", padx=4)

    def _color_button(self, parent, label, key, width=8):
        frame = tk.Frame(parent, bg=_C["panel"])
        frame.pack(side="left", padx=3)
        tk.Label(frame, text=label, bg=_C["panel"], fg=_C["text"], font=("Segoe UI", 8, "bold")).pack(side="left", padx=(0, 2))
        btn = tk.Button(frame, text="", width=width, bg=self._vars_graph[key].get(), relief="flat", command=lambda k=key: self._pick_color(k))
        btn.pack(side="left")
        self._color_buttons[key] = btn

    def _pick_color(self, key):
        picked = colorchooser.askcolor(color=self._vars_graph[key].get(), parent=self.root)[1]
        if picked:
            self._vars_graph[key].set(picked)
            if key in self._color_buttons:
                self._color_buttons[key].config(bg=picked)
            if hasattr(self, "canvas"):
                self._apply_custom_visual_theme()
                self.canvas.draw_idle()

    def _build_table_tab(self, parent):
        top = tk.Frame(parent, bg=_C["panel"], pady=10)
        top.pack(fill="x")
        tk.Button(top, text=" Refresh from Graph Dirs ", command=self._refresh_thermo_table, bg=_C["accent"], fg="white", font=("Segoe UI", 10, "bold"), padx=10).pack(side="left", padx=10)
        tk.Button(top, text=" Export CSV ", command=self._export_thermo_csv, bg=_C["green"], fg="white", font=("Segoe UI", 10, "bold"), padx=10).pack(side="left", padx=10)
        tk.Label(top, text="Scans the Series 1 & 2 directories from the Graph tab to tabulate all energies.", bg=_C["panel"], fg=_C["muted"], font=("Segoe UI", 10, "italic")).pack(side="left", padx=20)
        cols = ("Series", "Ratio", "Method", "Electronic (Eh)", "Gibbs (Eh)", "Rel. Energy (kcal/mol)")
        style = ttk.Style()
        style.configure("Custom.Treeview", font=("Consolas", 10), rowheight=25)
        style.configure("Custom.Treeview.Heading", font=("Segoe UI", 10, "bold"))
        self.thermo_tree = ttk.Treeview(parent, columns=cols, show="headings", style="Custom.Treeview")
        for col in cols:
            self.thermo_tree.heading(col, text=col)
            self.thermo_tree.column(col, width=150, anchor="center")
        scroll = ttk.Scrollbar(parent, orient="vertical", command=self.thermo_tree.yview)
        self.thermo_tree.configure(yscrollcommand=scroll.set)
        self.thermo_tree.pack(side="left", fill="both", expand=True, padx=10, pady=10)
        scroll.pack(side="right", fill="y", pady=10)

    def _extract_all_thermo(self, base_dir, method_col):
        data = []
        if not base_dir or not os.path.isdir(base_dir):
            return []
        for root, dirs, files in os.walk(base_dir):
            if "Energy_Comparison.csv" in files:
                csv_path = os.path.join(root, "Energy_Comparison.csv")
                parent_folder = os.path.basename(os.path.dirname(root))
                match = re.search(r"(\d+)[_:]+(\d+)", parent_folder)
                if match:
                    ratio_str = f"{match.group(1)}:{match.group(2)}"
                    sort_key = int(match.group(1)) + int(match.group(2))
                else:
                    ratio_str = parent_folder
                    sort_key = 999
                try:
                    with open(csv_path, "r", encoding="utf-8") as f:
                        for line in f.readlines()[1:]:
                            parts = line.strip().split(",")
                            if len(parts) >= 4 and parts[0] == method_col:
                                data.append({"ratio": ratio_str, "ele": parts[2].strip(), "gibbs": parts[3].strip(), "sort": sort_key})
                except Exception as e:
                    print(f"Error parsing {csv_path}: {e}")
        data.sort(key=lambda x: x["sort"])
        return data

    def _refresh_thermo_table(self):
        for row in self.thermo_tree.get_children():
            self.thermo_tree.delete(row)
        for d, name in [(self._vars_graph["dir1"].get(), self._vars_graph["name1"].get().strip() or "Series 1"), (self._vars_graph["dir2"].get(), self._vars_graph["name2"].get().strip() or "Series 2")]:
            if not d:
                continue
            data = self._extract_all_thermo(d, self._vars_graph["method"].get())
            if not data:
                continue
            base_e, is_gibbs = None, False
            for item in data:
                if item["gibbs"] not in ("N/A", "0.000000"):
                    base_e = float(item["gibbs"])
                    is_gibbs = True
                    break
                if item["ele"] not in ("N/A", "0.000000"):
                    base_e = float(item["ele"])
                    break
            for item in data:
                rel_e = "N/A"
                if base_e is not None:
                    val = float(item["gibbs"]) if is_gibbs and item["gibbs"] != "N/A" else float(item["ele"])
                    rel_e = f"{(val - base_e) * EH2KCAL:.2f}"
                self.thermo_tree.insert("", "end", values=(name, item["ratio"], self._vars_graph["method"].get(), item["ele"], item["gibbs"], rel_e))

    def _export_thermo_csv(self):
        dest = filedialog.asksaveasfilename(defaultextension=".csv", filetypes=[("CSV", "*.csv")], title="Export Thermodynamic Data")
        if dest:
            with open(dest, "w", newline="", encoding="utf-8") as f:
                writer = csv.writer(f)
                writer.writerow(["Series", "Ratio", "Method", "Electronic (Eh)", "Gibbs (Eh)", "Rel. Energy (kcal/mol)"])
                for row_id in self.thermo_tree.get_children():
                    writer.writerow(self.thermo_tree.item(row_id)["values"])
            messagebox.showinfo("Success", f"Data exported to {dest}")

    def _extract_dir_energies(self, base_dir, method_col, metric_type):
        data = {}
        if not base_dir or not os.path.isdir(base_dir):
            return []
        col_idx = 3 if "Gibbs" in metric_type else 2
        for root, dirs, files in os.walk(base_dir):
            if "Energy_Comparison.csv" in files:
                csv_path = os.path.join(root, "Energy_Comparison.csv")
                parent_folder = os.path.basename(os.path.dirname(root))
                match = re.search(r"(\d+)[_:]+(\d+)", parent_folder)
                if match:
                    ratio_str = f"{match.group(1)}:{match.group(2)}"
                    sort_key = int(match.group(1)) + int(match.group(2))
                else:
                    ratio_str = parent_folder
                    sort_key = 999
                try:
                    with open(csv_path, "r", encoding="utf-8") as f:
                        for line in f.readlines()[1:]:
                            parts = line.strip().split(",")
                            if len(parts) > col_idx and parts[0] == method_col:
                                val_str = parts[col_idx].strip()
                                if val_str not in ("N/A", "0.000000"):
                                    data[ratio_str] = {"energy": float(val_str), "sort": sort_key}
                except Exception as e:
                    print(f"Error parsing {csv_path}: {e}")
        return [(k, v["energy"]) for k, v in sorted(data.items(), key=lambda x: x[1]["sort"])]

    def _apply_thesis_theme(self):
        self.ax.clear()
        if self.ax2 is not None:
            self.ax2.remove()
            self.ax2 = None
        bg = self._vars_graph["graph_bg"].get()
        axis_c = self._vars_graph["axis_color"].get()
        self.ax.set_facecolor(bg)
        self.fig.patch.set_facecolor(bg)
        self.ax.tick_params(colors=axis_c, labelsize=11)
        self.ax.xaxis.label.set_color(axis_c)
        self.ax.yaxis.label.set_color(axis_c)
        self.ax.title.set_color(axis_c)
        for spine in self.ax.spines.values():
            spine.set_edgecolor(axis_c)
            spine.set_linewidth(1.5)

    def _make_bbox(self, style, edge_c="#333", alpha=0.85, lw=0.7):
        if style == "No Box":
            return None
        bs = "round,pad=0.25" if style == "Rounded Box" else "square,pad=0.25"
        return dict(boxstyle=bs, facecolor=self._vars_graph["graph_bg"].get(), edgecolor=edge_c, alpha=alpha, linewidth=lw)

    def _plot_reaction_profile(self):
        dir1, dir2 = self._vars_graph["dir1"].get(), self._vars_graph["dir2"].get()
        name1 = self._vars_graph["name1"].get().strip() or "Series 1"
        name2 = self._vars_graph["name2"].get().strip() or "Series 2"
        method = self._vars_graph["method"].get()
        metric = self._vars_graph["metric"].get()
        custom_x = [l.strip() for l in self._vars_graph["x_labels"].get().split(",") if l.strip()]
        conn_style = self._vars_graph["conn_style"].get()
        arrow_style = self._vars_graph["delta_arrow"].get()
        lbl_box_style = self._vars_graph["label_box"].get()
        if method in ["CREST", "xTB"] and "Gibbs" in metric:
            messagebox.showwarning("Metric Auto-Corrected", f"{method} does not compute Gibbs Free Energy. Switching to Total Electronic Energy.")
            self._vars_graph["metric"].set("Total Energy (ΔE)")
            metric = "Total Energy (ΔE)"
        self._apply_thesis_theme()
        data1 = self._extract_dir_energies(dir1, method, metric) if dir1 else []
        data2 = self._extract_dir_energies(dir2, method, metric) if dir2 else []
        if not data1 and not data2:
            messagebox.showerror("Error", "No valid data found.")
            return
        all_labels = []
        for d in [data1, data2]:
            for lbl, _ in d:
                if lbl not in all_labels:
                    all_labels.append(lbl)

        def sort_ratio(lbl):
            try:
                a, b = map(int, lbl.replace("_", ":").split(":"))
                return a + b
            except Exception:
                return 999

        all_labels.sort(key=sort_ratio)
        final_labels = custom_x if custom_x and len(custom_x) >= len(all_labels) else all_labels
        label_to_x = {lbl: i for i, lbl in enumerate(all_labels)}
        all_energies = []
        for d in [data1, data2]:
            if d:
                base_e = d[0][1]
                all_energies.extend([(e - base_e) * EH2KCAL for _, e in d])
        y_span = max(all_energies) - min(all_energies) if all_energies else 1.0
        if y_span == 0:
            y_span = 1.0
        y_pad = y_span * 0.05

        def plot_series(data, color, label, is_series_one):
            if not data:
                return
            base_e = data[0][1]
            rel_energies = [(e - base_e) * EH2KCAL for _, e in data]
            x_positions = [label_to_x[lbl] for lbl, _ in data]
            for i in range(len(data)):
                x_pos, y_val = x_positions[i], rel_energies[i]
                self.ax.hlines(y_val, x_pos - 0.35, x_pos + 0.35, color=color, lw=5, label=label if i == 0 else "")
                val_pad = y_pad if is_series_one else -y_pad
                va_align = "bottom" if is_series_one else "top"
                txt_kw = dict(ha="center", va=va_align, color=color, fontweight="bold", fontsize=11, zorder=9)
                bbox = self._make_bbox(lbl_box_style, edge_c=color, lw=1.5)
                if bbox:
                    txt_kw["bbox"] = bbox
                t = self.ax.text(x_pos, y_val + val_pad, f"{y_val:.1f}", **txt_kw)
                try:
                    t.set_draggable(True)
                except Exception:
                    pass
                if i > 0:
                    prev_x, prev_y = x_positions[i - 1], rel_energies[i - 1]
                    start_x, start_y = prev_x + 0.35, prev_y
                    end_x, end_y = x_pos - 0.35, y_val
                    kw_line = dict(color=color, alpha=0.60, lw=2.5, zorder=3)
                    if conn_style == "Straight":
                        self.ax.plot([start_x, end_x], [start_y, end_y], **kw_line)
                    elif conn_style in ("Smooth (Bezier)", "Smooth"):
                        xs = np.linspace(start_x, end_x, 60)
                        ts = (xs - start_x) / max(end_x - start_x, 1e-9)
                        self.ax.plot(xs, start_y + (end_y - start_y) * (3 * ts**2 - 2 * ts**3), **kw_line)
                    elif conn_style == "Wavy":
                        _draw_wavy(self.ax, start_x, start_y, end_x, end_y, n_waves=6, amp_frac=0.04, **kw_line)
                    elif conn_style == "Zigzag":
                        _draw_zigzag(self.ax, start_x, start_y, end_x, end_y, n_zigs=8, amp_frac=0.04, **kw_line)
                    elif conn_style == "Curled (Arc)":
                        _draw_curled(self.ax, start_x, start_y, end_x, end_y, rad=0.30, **kw_line)
                    delta = y_val - prev_y
                    mid_x, mid_y = (start_x + end_x) / 2, (start_y + end_y) / 2
                    if arrow_style != "No Arrow":
                        _draw_fancy_arrow(self.ax, mid_x, mid_y + (y_pad if is_series_one else -y_pad) * 0.8, mid_x, mid_y, style=arrow_style, color=color, lw=1.5, alpha=0.8, zorder=11)
                    d_kw = dict(ha="center", va="center", color=color, fontsize=10, fontweight="bold", zorder=12)
                    d_bbox = self._make_bbox(lbl_box_style, edge_c=color, alpha=0.95, lw=1.5)
                    if d_bbox:
                        d_kw["bbox"] = d_bbox
                    t_del = self.ax.text(mid_x, mid_y + (y_pad if is_series_one else -y_pad) * 0.8, f"$\\Delta$={delta:.1f}", **d_kw)
                    try:
                        t_del.set_draggable(True)
                    except Exception:
                        pass

        plot_series(data1, self._vars_graph["series1_color"].get(), name1, True)
        plot_series(data2, self._vars_graph["series2_color"].get(), name2, False)
        self.ax.set_xticks(range(len(all_labels)))
        self.ax.set_xticklabels(final_labels[: len(all_labels)], fontweight="bold", fontsize=12)
        self.ax.set_xlim(left=-0.75, right=len(all_labels) - 0.25)
        if self._vars_graph["top_space"].get():
            y_min, y_max = self.ax.get_ylim()
            self.ax.set_ylim(y_min - y_span * 0.1, y_max + y_span * 0.5)
        self.ax2 = self.ax.twinx()
        y1_min, y1_max = self.ax.get_ylim()
        self.ax2.set_ylim(y1_min * KCAL2KJ, y1_max * KCAL2KJ)
        self.ax2.set_ylabel("Relative Energy (kJ/mol)", fontweight="bold", fontsize=14, color=self._vars_graph["axis_color"].get())
        self.ax2.tick_params(colors=self._vars_graph["axis_color"].get(), labelsize=11)
        self.ax2.spines["right"].set_linewidth(1.5)
        symb = "ΔG" if "Gibbs" in metric else "ΔE"
        self.ax.set_ylabel(f"Relative Energy {symb} (kcal/mol)", fontweight="bold", fontsize=14)
        self.ax.set_title(f"Reaction Energy Profile Diagram ({symb})", fontweight="bold", fontsize=18, pad=20)
        self.ax.grid(axis="y", linestyle="-", alpha=0.25, color=self._vars_graph["grid_color"].get())
        self.ax.legend(loc="lower left" if self._vars_graph["top_space"].get() else "best", framealpha=0.95, fontsize=12, edgecolor=self._vars_graph["axis_color"].get())
        watermark = f"Level of Theory: {method}"
        if method == "ORCA":
            watermark += f" ({self._vars['orca_method'].get()})"
        self.ax.text(0.99, 0.02, watermark, ha="right", va="bottom", color=self._vars_graph["watermark_color"].get(), fontsize=10, transform=self.ax.transAxes, style="italic")
        self.fig.tight_layout()
        self.canvas.draw()

    def _plot_methods(self):
        dir1 = self._vars_graph["dir1"].get()
        name1 = self._vars_graph["name1"].get().strip() or "Series"
        metric = self._vars_graph["metric"].get()
        custom_x = [l.strip() for l in self._vars_graph["x_labels"].get().split(",") if l.strip()]
        if "Gibbs" in metric:
            messagebox.showwarning("Metric Auto-Corrected", "CREST vs ORCA comparison uses Total Energy.")
            self._vars_graph["metric"].set("Total Energy (ΔE)")
            metric = "Total Energy (ΔE)"
        if not dir1:
            messagebox.showerror("Error", "Please specify Series 1 Directory.")
            return
        crest_data = self._extract_dir_energies(dir1, "CREST", metric)
        orca_data = self._extract_dir_energies(dir1, "ORCA", metric)
        self._apply_thesis_theme()
        plotted_any = False
        all_labels = []
        if crest_data:
            plotted_any = True
            all_labels = [l for l, _ in crest_data]
            x_c = np.arange(len(all_labels))
            base_c = crest_data[0][1]
            y_c = [(e - base_c) * EH2KCAL for _, e in crest_data]
            self.ax.plot(x_c, y_c, marker="o", markersize=10, linestyle="--", linewidth=2.5, color=self._vars_graph["crest_color"].get(), label="CREST Trend")
        if orca_data:
            plotted_any = True
            all_labels = [l for l, _ in orca_data]
            x_o = np.arange(len(all_labels))
            base_o = orca_data[0][1]
            y_o = [(e - base_o) * EH2KCAL for _, e in orca_data]
            self.ax.plot(x_o, y_o, marker="s", markersize=10, linestyle="-", linewidth=3.5, color=self._vars_graph["orca_color"].get(), label="ORCA DFT Trend")
        if not plotted_any:
            messagebox.showerror("Error", "No valid data found.")
            return
        final_labels = custom_x if custom_x and len(custom_x) >= len(all_labels) else all_labels
        self.ax.set_xticks(np.arange(len(all_labels)))
        self.ax.set_xticklabels(final_labels[: len(all_labels)], fontweight="bold", fontsize=12)
        if self._vars_graph["top_space"].get():
            y_min, y_max = self.ax.get_ylim()
            self.ax.set_ylim(y_min, y_max + ((y_max - y_min) * 0.5))
        self.ax.set_ylabel("Relative Energy ΔE (kcal/mol)", fontsize=14, fontweight="bold", color=self._vars_graph["axis_color"].get())
        self.ax.set_title(f"Method Accuracy Comparison: {name1}", fontsize=18, fontweight="bold", color=self._vars_graph["axis_color"].get(), pad=20)
        self.ax.legend(loc="best", framealpha=0.95, fontsize=12, edgecolor=self._vars_graph["axis_color"].get())
        self.ax.grid(True, alpha=0.3, color=self._vars_graph["grid_color"].get(), linestyle="-")
        self.ax.text(0.99, 0.02, f"Level of Theory (DFT): {self._vars['orca_method'].get()}", ha="right", va="bottom", color=self._vars_graph["watermark_color"].get(), fontsize=10, transform=self.ax.transAxes, style="italic")
        self.fig.tight_layout()
        self.canvas.draw()

    def _collect_plot_records(self, base_dir, method_filter="ALL", metric="Total Energy (ΔE)", top_n=0):
        if not base_dir or not os.path.isdir(base_dir):
            return []
        job_dirs = self._find_iak_job_dirs(base_dir)
        if not job_dirs:
            job_dirs = [base_dir]
        records = []
        method_filter = method_filter.upper()
        for job_dir in job_dirs:
            job_label = os.path.basename(job_dir.rstrip("\\/")) or "job"
            job_records = self._collect_existing_job_records(job_dir)
            for method, method_records in job_records.items():
                if method_filter != "ALL" and method.upper() != method_filter:
                    continue
                for rec in method_records:
                    energy = rec.get("gibbs") if ("Gibbs" in metric and method == "ORCA" and rec.get("gibbs")) else rec.get("energy")
                    if energy is None:
                        continue
                    records.append({
                        "method": method,
                        "job": job_label,
                        "source": rec.get("source", ""),
                        "energy": float(energy),
                        "label": f"{job_label}:{rec.get('source','')}",
                    })
            if method_filter in ("ALL", "CREST"):
                crest_dir = os.path.join(job_dir, "03_CREST_Results")
                if os.path.isdir(crest_dir):
                    for conf_path in Path(crest_dir).glob("crest_conformers_*.xyz"):
                        for idx, mol in enumerate(read_multi_xyz(str(conf_path))):
                            if mol.energy_eh:
                                records.append({
                                    "method": "CREST",
                                    "job": job_label,
                                    "source": f"{conf_path.stem}_conf_{idx}",
                                    "energy": float(mol.energy_eh),
                                    "label": f"{job_label}:{conf_path.stem}_{idx}",
                                })
        if top_n and top_n > 0:
            limited = []
            for method in ["xTB", "CREST", "ORCA"]:
                group = [r for r in records if r["method"] == method]
                group.sort(key=lambda item: item["energy"])
                limited.extend(group[:top_n])
            records = limited
        return records

    def _plot_records_by_method(self, records, title):
        if not records:
            messagebox.showerror("No Data", "No optimized/conformer energy records were found for plotting.")
            return
        self._apply_thesis_theme()
        colors = {"xTB": self._vars_graph["series1_color"].get(), "CREST": self._vars_graph["crest_color"].get(), "ORCA": self._vars_graph["orca_color"].get()}
        x_labels = []
        x_pos = 0
        for method in ["xTB", "CREST", "ORCA"]:
            group = [r for r in records if r["method"] == method]
            if not group:
                continue
            group.sort(key=lambda item: item["energy"])
            base = group[0]["energy"]
            xs, ys = [], []
            for rec in group:
                xs.append(x_pos)
                ys.append((rec["energy"] - base) * EH2KCAL)
                x_labels.append(rec["label"])
                x_pos += 1
            self.ax.plot(xs, ys, marker="o", linewidth=2.2, markersize=7, color=colors.get(method, "#333333"), label=f"{method} relative")
        self.ax.set_xticks(range(len(x_labels)))
        self.ax.set_xticklabels(x_labels, rotation=60, ha="right", fontsize=8)
        self.ax.set_ylabel("Relative Energy (kcal/mol)", fontweight="bold", fontsize=13)
        self.ax.set_title(title, fontweight="bold", fontsize=16, pad=16)
        self.ax.grid(axis="y", linestyle="-", alpha=0.25, color=self._vars_graph["grid_color"].get())
        self.ax.legend(loc="best", framealpha=0.95, fontsize=10)
        self.fig.tight_layout()
        self.canvas.draw()

    def _plot_all_optimized(self):
        base_dir = self._vars_graph["dir1"].get() or self._vars["out"].get()
        records = self._collect_plot_records(base_dir, method_filter="ALL", metric=self._vars_graph["metric"].get(), top_n=0)
        self._plot_records_by_method(records, "All Optimized Structures: xTB vs CREST vs ORCA")

    def _plot_top_conformers(self):
        base_dir = self._vars_graph["dir1"].get() or self._vars["out"].get()
        try:
            top_n = max(1, int(self._vars_graph["conformer_top_n"].get()))
        except Exception:
            top_n = 10
        method = self._vars_graph["method"].get()
        records = self._collect_plot_records(base_dir, method_filter=method, metric=self._vars_graph["metric"].get(), top_n=top_n)
        self._plot_records_by_method(records, f"Top {top_n} {method} Conformers / Optimized Structures")

    def _apply_custom_visual_theme(self):
        if not hasattr(self, "ax"):
            return
        bg = self._vars_graph["graph_bg"].get()
        axis_c = self._vars_graph["axis_color"].get()
        grid_c = self._vars_graph["grid_color"].get()
        self.fig.patch.set_facecolor(bg)
        for axis in [self.ax] + ([self.ax2] if getattr(self, "ax2", None) else []):
            axis.set_facecolor(bg)
            axis.tick_params(colors=axis_c)
            axis.xaxis.label.set_color(axis_c)
            axis.yaxis.label.set_color(axis_c)
            axis.title.set_color(axis_c)
            for spine in axis.spines.values():
                spine.set_edgecolor(axis_c)
            for gridline in axis.get_xgridlines() + axis.get_ygridlines():
                gridline.set_color(grid_c)
                gridline.set_alpha(0.35)

    def _parse_crop(self, value):
        value = (value or "").strip()
        return None if not value else float(value)

    def _apply_graph_crop(self):
        try:
            xmin = self._parse_crop(self._vars_graph["crop_xmin"].get())
            xmax = self._parse_crop(self._vars_graph["crop_xmax"].get())
            ymin = self._parse_crop(self._vars_graph["crop_ymin"].get())
            ymax = self._parse_crop(self._vars_graph["crop_ymax"].get())
        except Exception:
            return messagebox.showerror("Crop Error", "Crop values must be numeric or blank.")
        if xmin is not None or xmax is not None:
            current = self.ax.get_xlim()
            self.ax.set_xlim(xmin if xmin is not None else current[0], xmax if xmax is not None else current[1])
        if ymin is not None or ymax is not None:
            current = self.ax.get_ylim()
            self.ax.set_ylim(ymin if ymin is not None else current[0], ymax if ymax is not None else current[1])
            if self.ax2 is not None:
                yl = self.ax.get_ylim()
                self.ax2.set_ylim(yl[0] * KCAL2KJ, yl[1] * KCAL2KJ)
        self._apply_custom_visual_theme()
        self.canvas.draw_idle()

    def _reset_graph_view(self):
        self.ax.relim()
        self.ax.autoscale_view()
        if self.ax2 is not None:
            yl = self.ax.get_ylim()
            self.ax2.set_ylim(yl[0] * KCAL2KJ, yl[1] * KCAL2KJ)
        for key in ["crop_xmin", "crop_xmax", "crop_ymin", "crop_ymax"]:
            self._vars_graph[key].set("")
        self._apply_custom_visual_theme()
        self.canvas.draw_idle()

    def _copy_to_clipboard(self, figure):
        import tempfile

        tmp_path = None
        try:
            with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as f:
                tmp_path = f.name
            figure.savefig(tmp_path, dpi=300, bbox_inches="tight")
            ps_script = (
                "Add-Type -AssemblyName System.Windows.Forms;"
                "Add-Type -AssemblyName System.Drawing;"
                f"$img = [System.Drawing.Image]::FromFile('{tmp_path}');"
                "[System.Windows.Forms.Clipboard]::SetImage($img);"
                "$img.Dispose();"
            )
            result = subprocess.run(["powershell", "-NoProfile", "-Command", ps_script], capture_output=True, timeout=10)
            if result.returncode == 0:
                messagebox.showinfo("Clipboard", "Plot copied to Windows Clipboard.")
            else:
                messagebox.showerror("Clipboard Error", result.stderr.decode(errors="replace"))
        except Exception as e:
            messagebox.showerror("Clipboard Error", str(e))
        finally:
            if tmp_path and os.path.exists(tmp_path):
                try:
                    os.unlink(tmp_path)
                except Exception:
                    pass

    def _export_high_res_graph(self):
        dest = filedialog.asksaveasfilename(defaultextension=".png", filetypes=[("PNG Image", "*.png"), ("PDF Document", "*.pdf"), ("JPEG Image", "*.jpg"), ("SVG Vector", "*.svg")], title="Export High-Res Graph")
        if dest:
            self.fig.savefig(dest, dpi=600, bbox_inches="tight")
            messagebox.showinfo("Export Successful", f"High-resolution publication graph saved to:\n{dest}")

    def _get_selected_filepath(self):
        for col, lb in enumerate(self.listboxes):
            sel = lb.curselection()
            if sel:
                fname = lb.get(sel[0])
                if fname.startswith("("):
                    return None
                return os.path.join(os.path.abspath(self._vars["out"].get().strip(" \"'")), self.folders[col], fname)
        return None

    def _open_file(self, listbox, col):
        sel = listbox.curselection()
        if not sel:
            return
        fname = listbox.get(sel[0])
        if fname.startswith("("):
            return
        path = os.path.join(os.path.abspath(self._vars["out"].get().strip(" \"'")), self.folders[col], fname)
        if not os.path.exists(path):
            return messagebox.showerror("Missing File", f"File not found:\n{path}")
        if path.lower().endswith(".xyz") and MATPLOTLIB_AVAILABLE:
            return self._open_model_viewer(path)
        try:
            if sys.platform == "win32":
                os.startfile(path)
            else:
                webbrowser.open(Path(path).as_uri())
        except Exception as exc:
            messagebox.showerror("Open Error", str(exc))

    def _open_in_avogadro(self):
        path = self._get_selected_filepath()
        if not path:
            return messagebox.showwarning("Select File", "Please select a file first.")
        cmd = None
        if shutil.which("avogadro"):
            cmd = ["avogadro", path]
        elif shutil.which("avogadro2"):
            cmd = ["avogadro2", path]
        elif sys.platform == "win32":
            cps = [
                r"C:\Program Files\Avogadro\bin\avogadro.exe",
                r"C:\Program Files (x86)\Avogadro\bin\avogadro.exe",
                r"C:\Program Files\Avogadro\avogadro.exe",
                r"C:\Program Files (x86)\Avogadro\avogadro.exe",
                r"C:\Program Files\Avogadro2\bin\avogadro2.exe",
                r"C:\Program Files\Avogadro2\avogadro2.exe",
            ]
            for cp in cps:
                if os.path.exists(cp):
                    cmd = [cp, path]
                    break
        if not cmd:
            return messagebox.showerror("Not Found", "Avogadro not found in standard paths.")
        try:
            subprocess.Popen(cmd)
        except Exception as e:
            messagebox.showerror("Error", f"Launch failed: {e}")

    def _export_file(self):
        path = self._get_selected_filepath()
        if not path:
            return messagebox.showwarning("Select File", "Please select a file first.")
        dest = filedialog.asksaveasfilename(defaultextension=Path(path).suffix, initialfile=os.path.basename(path), title="Export")
        if dest:
            shutil.copy2(path, dest)

    def _on_hover(self, event, listbox, col):
        idx = listbox.nearest(event.y)
        bbox = listbox.bbox(idx)
        if not bbox or not (bbox[1] <= event.y <= bbox[1] + bbox[3]):
            return self._hide_preview()
        fname = listbox.get(idx)
        if fname.startswith("("):
            return self._hide_preview()
        filepath = os.path.join(os.path.abspath(self._vars["out"].get().strip(" \"'")), self.folders[col], fname)
        if not os.path.exists(filepath) or not fname.endswith(".xyz"):
            return
        if self.preview_file == filepath and self.preview_tw is not None:
            return
        self._show_preview(event.x_root, event.y_root, filepath)

    def _hide_preview(self):
        if self.preview_tw:
            self.preview_tw.destroy()
            self.preview_tw = None
        self.preview_file = None

    def _read_xyz_atoms(self, filepath):
        with open(filepath, "r", encoding="utf-8", errors="replace") as f:
            lines = f.readlines()
        n = int(lines[0].strip())
        return [(p[0], float(p[1]), float(p[2]), float(p[3])) for l in lines[2 : 2 + n] for p in [l.split()] if len(p) >= 4]

    def _model_projection(self, atoms):
        coords = np.array([[a[1], a[2], a[3]] for a in atoms], dtype=float)
        coords -= np.mean(coords, axis=0)
        rx = np.array([[1, 0, 0], [0, math.cos(0.35), -math.sin(0.35)], [0, math.sin(0.35), math.cos(0.35)]])
        ry = np.array([[math.cos(0.61), 0, math.sin(0.61)], [0, 1, 0], [-math.sin(0.61), 0, math.cos(0.61)]])
        return coords @ rx @ ry

    def _atom_color(self, symbol):
        if symbol == "Cl":
            key = "atom_Cl"
        elif symbol in {"H", "C", "N", "O", "F", "S"}:
            key = f"atom_{symbol}"
        else:
            key = "atom_other"
        return self._vars_graph[key].get()

    def _covalent_radius(self, symbol):
        radii = {
            "H": 0.31,
            "B": 0.85,
            "C": 0.76,
            "N": 0.71,
            "O": 0.66,
            "F": 0.57,
            "P": 1.07,
            "S": 1.05,
            "Cl": 1.02,
            "Br": 1.20,
            "I": 1.39,
            "Si": 1.11,
            "Na": 1.66,
            "K": 2.03,
            "Mg": 1.41,
            "Ca": 1.76,
            "Fe": 1.24,
            "Cu": 1.32,
            "Zn": 1.22,
        }
        return radii.get(symbol, 0.80)

    def _vdw_radius(self, symbol):
        radii = {
            "H": 1.20,
            "C": 1.70,
            "N": 1.55,
            "O": 1.52,
            "F": 1.47,
            "P": 1.80,
            "S": 1.80,
            "Cl": 1.75,
            "Br": 1.85,
            "I": 1.98,
        }
        return radii.get(symbol, 1.65)

    def _float_var(self, key, fallback):
        try:
            return float(self._vars_graph[key].get())
        except Exception:
            return fallback

    def _detect_bonds(self, atoms, coords):
        mode = self._vars_graph["bond_mode"].get()
        cutoff = self._float_var("bond_cutoff", 1.85)
        tolerance = self._float_var("bond_tolerance", 0.45)
        bonds = []
        for i in range(len(atoms)):
            for j in range(i + 1, len(atoms)):
                dist = float(np.linalg.norm(coords[i] - coords[j]))
                sym_i, sym_j = atoms[i][0], atoms[j][0]
                include = False
                kind = "bond"
                if mode == "Distance cutoff":
                    include = dist <= cutoff
                elif mode == "Hydrogen bonds":
                    pair = {sym_i, sym_j}
                    include = "H" in pair and bool(pair.intersection({"O", "N", "F", "S"})) and 1.3 <= dist <= cutoff
                    kind = "hbond"
                elif mode == "All close contacts":
                    include = dist <= cutoff
                    kind = "contact"
                else:
                    include = dist <= self._covalent_radius(sym_i) + self._covalent_radius(sym_j) + tolerance
                if include and dist > 0.25:
                    bonds.append((i, j, dist, kind))
        return bonds

    def _set_3d_equal_axes(self, axis, coords):
        mins = coords.min(axis=0)
        maxs = coords.max(axis=0)
        centers = (mins + maxs) / 2.0
        radius = max(float(np.max(maxs - mins)) / 2.0, 1.0)
        pad = radius * 0.22
        radius += pad
        axis.set_xlim(centers[0] - radius, centers[0] + radius)
        axis.set_ylim(centers[1] - radius, centers[1] + radius)
        axis.set_zlim(centers[2] - radius, centers[2] + radius)
        try:
            axis.set_box_aspect((1, 1, 1))
        except Exception:
            pass

    def _render_xyz_model(self, axis, filepath):
        axis.clear()
        axis.set_facecolor(self._vars_graph["model_bg"].get())
        axis.figure.patch.set_facecolor(self._vars_graph["model_bg"].get())
        atoms = self._read_xyz_atoms(filepath)
        coords = np.array([[a[1], a[2], a[3]] for a in atoms], dtype=float)
        coords -= np.mean(coords, axis=0)
        style = self._vars_graph["model_style"].get()
        bond_color = self._vars_graph["bond_color"].get()
        bond_lw = {"Wireframe": 1.1, "Ball-and-stick": 2.6, "Licorice": 6.0, "CPK": 0.8}.get(style, 2.6)
        atom_scale = {"Wireframe": 18, "Ball-and-stick": 180, "Licorice": 90, "CPK": 520}.get(style, 180)
        bonds = self._detect_bonds(atoms, coords)
        if style != "CPK":
            for i, j, dist, kind in bonds:
                linestyle = "--" if kind in {"hbond", "contact"} else "-"
                alpha = 0.55 if kind in {"hbond", "contact"} else 0.95
                axis.plot(
                    [coords[i, 0], coords[j, 0]],
                    [coords[i, 1], coords[j, 1]],
                    [coords[i, 2], coords[j, 2]],
                    color=bond_color,
                    linewidth=bond_lw,
                    linestyle=linestyle,
                    solid_capstyle="round",
                    alpha=alpha,
                    zorder=1,
                )
                if self._vars_graph["show_bond_distances"].get():
                    mid = (coords[i] + coords[j]) / 2.0
                    axis.text(
                        mid[0],
                        mid[1],
                        mid[2],
                        f"{dist:.2f} A",
                        color=self._vars_graph["axis_color"].get(),
                        fontsize=7,
                        ha="center",
                        va="center",
                    )
        for idx in range(len(atoms)):
            symbol = atoms[idx][0]
            base = self._vdw_radius(symbol) if style == "CPK" else self._covalent_radius(symbol) + 0.30
            size = atom_scale * (base**2)
            if style == "Wireframe":
                size *= 0.35
            elif style == "CPK":
                size *= 0.75
            elif style == "Licorice":
                size *= 0.55
            axis.scatter(
                coords[idx, 0],
                coords[idx, 1],
                coords[idx, 2],
                s=size,
                color=self._atom_color(symbol),
                edgecolors="#111827",
                linewidths=0.65,
                depthshade=True,
                zorder=3,
            )
            if self._vars_graph["show_atom_labels"].get():
                axis.text(
                    coords[idx, 0],
                    coords[idx, 1],
                    coords[idx, 2],
                    f"{symbol}{idx + 1}",
                    color=self._vars_graph["axis_color"].get(),
                    fontsize=8,
                )
        self._set_3d_equal_axes(axis, coords)
        axis.view_init(elev=self._float_var("model_elev", 18.0), azim=self._float_var("model_azim", 38.0))
        if self._vars_graph["show_3d_axes"].get():
            axis.set_xlabel("X (A)", color=self._vars_graph["axis_color"].get(), fontsize=8)
            axis.set_ylabel("Y (A)", color=self._vars_graph["axis_color"].get(), fontsize=8)
            axis.set_zlabel("Z (A)", color=self._vars_graph["axis_color"].get(), fontsize=8)
            axis.tick_params(colors=self._vars_graph["axis_color"].get(), labelsize=7)
        else:
            axis.set_axis_off()
        for pane in [axis.xaxis.pane, axis.yaxis.pane, axis.zaxis.pane]:
            pane.set_facecolor(self._vars_graph["model_bg"].get())
            pane.set_edgecolor("#d0d7de")

    def _show_preview(self, x, y, filepath):
        self._hide_preview()
        self.preview_file = filepath
        tw = tk.Toplevel(self.root)
        self.preview_tw = tw
        tw.wm_overrideredirect(True)
        tw.geometry(f"+{x + 15}+{y + 15}")
        cv = tk.Canvas(tw, width=270, height=270, bg=self._vars_graph["model_bg"].get(), highlightthickness=2, highlightbackground=_C["accent"])
        cv.pack()
        try:
            atoms = self._read_xyz_atoms(filepath)
            coords = self._model_projection(atoms)
            xy, z = coords[:, :2], coords[:, 2]
            span = max(float(np.max(np.abs(xy))), 1.0)
            xy = xy * (105.0 / span) + 135.0
            distances = np.linalg.norm(coords[:, None, :] - coords[None, :, :], axis=2)
            style = self._vars_graph["model_style"].get()
            bond_width = {"Wireframe": 1, "Ball-and-stick": 3, "Licorice": 6, "CPK": 1}.get(style, 3)
            if style != "CPK":
                for i in range(len(atoms)):
                    for j in range(i + 1, len(atoms)):
                        if 0.4 < distances[i, j] < 1.85:
                            cv.create_line(xy[i, 0], xy[i, 1], xy[j, 0], xy[j, 1], fill=self._vars_graph["bond_color"].get(), width=bond_width)
            for idx in sorted(range(len(atoms)), key=lambda i: z[i]):
                symbol = atoms[idx][0]
                radius = 4 if symbol == "H" else 8
                if style == "CPK":
                    radius = 11 if symbol == "H" else 18
                elif style == "Wireframe":
                    radius = 3
                elif style == "Licorice":
                    radius = 5 if symbol == "H" else 7
                cv.create_oval(xy[idx, 0] - radius, xy[idx, 1] - radius, xy[idx, 0] + radius, xy[idx, 1] + radius, fill=self._atom_color(symbol), outline="#111827")
        except Exception:
            cv.create_text(135, 135, text="Preview unavailable", fill="white")

    def _open_model_viewer(self, filepath):
        win = tk.Toplevel(self.root)
        win.title(f"IAK Model Viewer - {os.path.basename(filepath)}")
        win.geometry("900x720")
        win.configure(bg=_C["bg"])
        controls = tk.Frame(win, bg=_C["panel"])
        controls.pack(fill="x", padx=8, pady=8)
        tk.Label(controls, text="Model style:", bg=_C["panel"], fg=_C["text"], font=("Segoe UI", 9, "bold")).pack(side="left", padx=(5, 4))
        ttk.Combobox(controls, textvariable=self._vars_graph["model_style"], values=["Ball-and-stick", "CPK", "Wireframe", "Licorice"], width=15, state="readonly").pack(side="left", padx=4)
        tk.Label(controls, text="Bonds:", bg=_C["panel"], fg=_C["text"], font=("Segoe UI", 9, "bold")).pack(side="left", padx=(10, 4))
        ttk.Combobox(
            controls,
            textvariable=self._vars_graph["bond_mode"],
            values=["Covalent radii", "Distance cutoff", "Hydrogen bonds", "All close contacts"],
            width=16,
            state="readonly",
        ).pack(side="left", padx=4)
        tk.Label(controls, text="Cutoff A:", bg=_C["panel"], fg=_C["muted"]).pack(side="left", padx=(8, 2))
        tk.Entry(controls, textvariable=self._vars_graph["bond_cutoff"], bg=_C["entry"], fg="white", width=5, bd=0).pack(side="left", ipady=2)
        tk.Label(controls, text="Tol A:", bg=_C["panel"], fg=_C["muted"]).pack(side="left", padx=(8, 2))
        tk.Entry(controls, textvariable=self._vars_graph["bond_tolerance"], bg=_C["entry"], fg="white", width=5, bd=0).pack(side="left", ipady=2)
        tk.Checkbutton(controls, text="Distances", variable=self._vars_graph["show_bond_distances"], bg=_C["panel"], fg=_C["yellow"], selectcolor="black").pack(side="left", padx=8)
        tk.Checkbutton(controls, text="Atom labels", variable=self._vars_graph["show_atom_labels"], bg=_C["panel"], fg=_C["yellow"], selectcolor="black").pack(side="left", padx=4)
        tk.Checkbutton(controls, text="Axes", variable=self._vars_graph["show_3d_axes"], bg=_C["panel"], fg=_C["yellow"], selectcolor="black").pack(side="left", padx=4)
        for label, key in [("Elev", "model_elev"), ("Azim", "model_azim")]:
            tk.Label(controls, text=f"{label}:", bg=_C["panel"], fg=_C["muted"]).pack(side="left", padx=(8, 2))
            tk.Entry(controls, textvariable=self._vars_graph[key], bg=_C["entry"], fg="white", width=5, bd=0).pack(side="left", ipady=2)
        fig = Figure(figsize=(7, 6), dpi=110)
        ax = fig.add_subplot(111, projection="3d")
        canvas = FigureCanvasTkAgg(fig, master=win)
        toolbar = NavigationToolbar2Tk(canvas, win, pack_toolbar=False)

        def redraw():
            self._render_xyz_model(ax, filepath)
            canvas.draw_idle()

        tk.Button(controls, text="Redraw", command=redraw, bg=_C["accent"], fg="white", relief="flat", font=("Segoe UI", 9, "bold")).pack(side="left", padx=8)
        tk.Button(controls, text="Export Image", command=lambda: self._export_model_figure(fig), bg=_C["green"], fg="white", relief="flat", font=("Segoe UI", 9, "bold")).pack(side="left", padx=4)
        tk.Label(controls, text="Drag to rotate; use toolbar pan/zoom to inspect and crop.", bg=_C["panel"], fg=_C["muted"], font=("Segoe UI", 9, "italic")).pack(side="left", padx=12)
        canvas.get_tk_widget().pack(fill="both", expand=True, padx=8, pady=(0, 8))
        toolbar.update()
        toolbar.pack(fill="x", padx=8, pady=(0, 8))
        redraw()

    def _export_model_figure(self, fig):
        dest = filedialog.asksaveasfilename(defaultextension=".png", filetypes=[("PNG Image", "*.png"), ("PDF Document", "*.pdf"), ("SVG Vector", "*.svg"), ("JPEG Image", "*.jpg")], title="Export Model Image")
        if dest:
            image_dpi = int(self._float_var("image_dpi", 600))
            image_dpi = max(72, min(image_dpi, 1200))
            fig.savefig(dest, dpi=image_dpi, bbox_inches="tight")
            messagebox.showinfo("Export Successful", f"Model image saved to:\n{dest}")

    def _refresh_results(self):
        out_dir = os.path.abspath(self._vars["out"].get().strip(" \"'"))
        for col, lb in enumerate(self.listboxes):
            lb.delete(0, tk.END)
            path = os.path.join(out_dir, self.folders[col])
            if os.path.exists(path):
                files = sorted([f for f in os.listdir(path) if f.lower().endswith((".xyz", ".json", ".csv", ".md", ".png", ".jpg", ".jpeg", ".svg", ".pdf"))])
                for f in files:
                    lb.insert(tk.END, f)
                if not files:
                    lb.insert(tk.END, "(Empty)")

    def _is_iak_job_dir(self, folder):
        if not os.path.isdir(folder):
            return False
        markers = [
            "01_Inputs_and_Clusters",
            "02_xTB_Results",
            "03_CREST_Results",
            "04_ORCA_Refinement",
            "05_Top_Models_Comparison",
            "state.json",
        ]
        return any(os.path.exists(os.path.join(folder, marker)) for marker in markers)

    def _find_iak_job_dirs(self, folder):
        folder = os.path.abspath(folder)
        if self._is_iak_job_dir(folder):
            return [folder]
        jobs = []
        for root, dirs, files in os.walk(folder):
            if self._is_iak_job_dir(root):
                jobs.append(root)
                dirs[:] = []
        return sorted(set(jobs))

    def _resolve_existing_path(self, job_dir, path_value):
        if not path_value:
            return None
        path_value = str(path_value)
        if os.path.exists(path_value):
            return path_value
        basename = os.path.basename(path_value)
        if not basename:
            return None
        for root, _, files in os.walk(job_dir):
            if basename in files:
                return os.path.join(root, basename)
        return None

    def _xyz_comment_energy(self, path):
        try:
            with open(path, "r", encoding="utf-8", errors="replace") as handle:
                lines = handle.readlines()
            if len(lines) < 2:
                return None
            return float(lines[1].split()[0])
        except Exception:
            return None

    def _parse_orca_out(self, out_path):
        energy, gibbs, imag, success = None, None, 0, False
        try:
            with open(out_path, "r", encoding="utf-8", errors="replace") as handle:
                content = handle.read()
            for line in content.splitlines():
                if "FINAL SINGLE POINT ENERGY" in line:
                    try:
                        energy = float(line.split()[-1])
                    except Exception:
                        pass
                elif "Final Gibbs free energy" in line:
                    try:
                        gibbs = float(line.split()[-2])
                    except Exception:
                        pass
                elif "*** imaginary mode ***" in line:
                    imag += 1
                elif "ORCA TERMINATED NORMALLY" in line:
                    success = True
            xyz_path = None
            folder = os.path.dirname(out_path)
            stem = Path(out_path).stem
            candidates = [
                os.path.join(folder, f"{stem}_trj.xyz"),
                os.path.join(folder, f"{stem}.xyz"),
            ]
            candidates.extend(str(p) for p in Path(folder).glob("*_trj.xyz"))
            candidates.extend(str(p) for p in Path(folder).glob("*.xyz"))
            for candidate in candidates:
                if os.path.exists(candidate):
                    xyz_path = candidate
                    break
            if xyz_path is None:
                extracted = os.path.join(folder, f"{stem}_final_from_out.xyz")
                if self._extract_orca_final_xyz(out_path, extracted):
                    xyz_path = extracted
            if energy is None:
                success = False
            return {
                "status": "success" if success and energy is not None else "failed",
                "energy": energy,
                "gibbs": gibbs or 0.0,
                "imag": imag,
                "path": xyz_path,
                "out_path": out_path,
            }
        except Exception:
            return {"status": "failed", "energy": None, "gibbs": 0.0, "imag": 0, "path": None, "out_path": out_path}

    def _extract_orca_final_xyz(self, out_path, dest_xyz):
        try:
            with open(out_path, "r", encoding="utf-8", errors="replace") as handle:
                lines = handle.readlines()
            blocks = []
            idx = 0
            while idx < len(lines):
                line = lines[idx]
                upper = line.upper()
                if "CARTESIAN COORDINATES" in upper and ("ANGSTROEM" in upper or "ANGSTROM" in upper):
                    idx += 1
                    block = []
                    while idx < len(lines):
                        raw = lines[idx].strip()
                        idx += 1
                        if not raw:
                            if block:
                                break
                            continue
                        if set(raw) <= {"-"}:
                            continue
                        parts = raw.split()
                        if len(parts) >= 4 and re.match(r"^[A-Z][a-z]?$", parts[0]):
                            try:
                                block.append((parts[0], float(parts[1]), float(parts[2]), float(parts[3])))
                                continue
                            except Exception:
                                pass
                        if block:
                            break
                    if block:
                        blocks.append(block)
                else:
                    idx += 1
            if not blocks:
                return None
            atoms = blocks[-1]
            os.makedirs(os.path.dirname(dest_xyz), exist_ok=True)
            with open(dest_xyz, "w", encoding="utf-8", newline="\n") as handle:
                handle.write(f"{len(atoms)}\nExtracted final coordinates from {os.path.basename(out_path)}\n")
                for sym, x, y, z in atoms:
                    handle.write(f"{sym:<4} {x:15.6f} {y:15.6f} {z:15.6f}\n")
            return dest_xyz
        except Exception:
            return None

    def _parse_orca_thermo(self, out_path):
        data = {
            "final_single_point_energy_eh": None,
            "zero_point_energy_eh": None,
            "total_thermal_energy_eh": None,
            "total_enthalpy_eh": None,
            "final_entropy_term_eh": None,
            "final_gibbs_free_energy_eh": None,
            "temperature_k": None,
            "pressure_atm": None,
            "imaginary_frequencies": 0,
            "normal_termination": False,
        }
        labels = {
            "FINAL SINGLE POINT ENERGY": "final_single_point_energy_eh",
            "Zero point energy": "zero_point_energy_eh",
            "Total thermal energy": "total_thermal_energy_eh",
            "Total Enthalpy": "total_enthalpy_eh",
            "Final entropy term": "final_entropy_term_eh",
            "Final Gibbs free energy": "final_gibbs_free_energy_eh",
        }
        try:
            with open(out_path, "r", encoding="utf-8", errors="replace") as handle:
                for line in handle:
                    stripped = line.strip()
                    for label, key in labels.items():
                        if label.lower() in stripped.lower():
                            floats = re.findall(r"[-+]?\d+\.\d+(?:[Ee][-+]?\d+)?", stripped)
                            if floats:
                                data[key] = float(floats[-1])
                    if "*** imaginary mode ***" in stripped:
                        data["imaginary_frequencies"] += 1
                    if "ORCA TERMINATED NORMALLY" in stripped:
                        data["normal_termination"] = True
                    if "Temperature" in stripped and "K" in stripped:
                        floats = re.findall(r"[-+]?\d+\.\d+(?:[Ee][-+]?\d+)?", stripped)
                        if floats:
                            data["temperature_k"] = float(floats[-1])
                    if "Pressure" in stripped and ("atm" in stripped.lower() or "ATM" in stripped):
                        floats = re.findall(r"[-+]?\d+\.\d+(?:[Ee][-+]?\d+)?", stripped)
                        if floats:
                            data["pressure_atm"] = float(floats[-1])
        except Exception:
            pass
        return data

    def _copy_unique_result(self, src, dest):
        if not src or not os.path.exists(src):
            return None
        os.makedirs(os.path.dirname(dest), exist_ok=True)
        if os.path.abspath(src) == os.path.abspath(dest):
            return dest
        shutil.copy2(src, dest)
        return dest

    def _collect_existing_job_records(self, job_dir):
        state_path = os.path.join(job_dir, "state.json")
        state = {"xtb": {}, "crest": {}, "orca": {}}
        if os.path.exists(state_path):
            try:
                with open(state_path, "r", encoding="utf-8") as handle:
                    state.update(json.load(handle))
            except Exception:
                pass

        records = {"xTB": [], "CREST": [], "ORCA": []}

        for key, value in state.get("xtb", {}).items():
            if value.get("status") == "success" and value.get("energy") is not None:
                path = self._resolve_existing_path(job_dir, value.get("path"))
                records["xTB"].append({"source": key, "energy": float(value["energy"]), "gibbs": None, "imag": None, "path": path})

        xtb_dir = os.path.join(job_dir, "02_xTB_Results")
        if not records["xTB"] and os.path.isdir(xtb_dir):
            for path in Path(xtb_dir).rglob("*.xyz"):
                energy = self._xyz_comment_energy(str(path))
                if energy is not None:
                    records["xTB"].append({"source": path.stem, "energy": energy, "gibbs": None, "imag": None, "path": str(path)})

        for key, value in state.get("crest", {}).items():
            if value.get("status") == "success":
                path = self._resolve_existing_path(job_dir, value.get("best_path") or value.get("path"))
                if path:
                    mols = read_multi_xyz(path)
                    energy = mols[0].energy_eh if mols else self._xyz_comment_energy(path)
                    if energy is not None:
                        records["CREST"].append({"source": key, "energy": float(energy), "gibbs": None, "imag": None, "path": path})

        crest_dir = os.path.join(job_dir, "03_CREST_Results")
        if os.path.isdir(crest_dir):
            for path in list(Path(crest_dir).glob("crest_best_*.xyz")) + list(Path(crest_dir).glob("crest_conformers_*.xyz")):
                mols = read_multi_xyz(str(path))
                energy = mols[0].energy_eh if mols else self._xyz_comment_energy(str(path))
                if energy is not None:
                    records["CREST"].append({"source": path.stem, "energy": float(energy), "gibbs": None, "imag": None, "path": str(path)})

        for key, value in state.get("orca", {}).items():
            if value.get("status") == "success" and value.get("energy") is not None:
                path = self._resolve_existing_path(job_dir, value.get("best_path") or value.get("path"))
                records["ORCA"].append({
                    "source": key,
                    "energy": float(value["energy"]),
                    "gibbs": float(value.get("gibbs", 0.0) or 0.0),
                    "imag": int(value.get("imag", 0) or 0),
                    "path": path,
                })

        orca_dir = os.path.join(job_dir, "04_ORCA_Refinement")
        if os.path.isdir(orca_dir):
            for out_path in Path(orca_dir).rglob("*.out"):
                parsed = self._parse_orca_out(str(out_path))
                if parsed["status"] == "success" and parsed["energy"] is not None:
                    records["ORCA"].append({
                        "source": out_path.stem,
                        "energy": float(parsed["energy"]),
                        "gibbs": float(parsed.get("gibbs", 0.0) or 0.0),
                        "imag": int(parsed.get("imag", 0) or 0),
                        "path": parsed.get("path"),
                    })
        return records

    def _write_existing_analysis_outputs(self, job_dir):
        comp_dir = os.path.join(job_dir, "05_Top_Models_Comparison")
        logs_dir = os.path.join(job_dir, "logs")
        os.makedirs(comp_dir, exist_ok=True)
        os.makedirs(logs_dir, exist_ok=True)
        records = self._collect_existing_job_records(job_dir)
        rows = []

        xtb_best = min(records["xTB"], key=lambda item: item["energy"], default=None)
        if xtb_best:
            self._copy_unique_result(xtb_best.get("path"), os.path.join(comp_dir, "1_BEST_xTB.xyz"))
            rows.append(["xTB", xtb_best["source"], f"{xtb_best['energy']:.6f}", "N/A", "N/A"])

        crest_best = min(records["CREST"], key=lambda item: item["energy"], default=None)
        if crest_best:
            self._copy_unique_result(crest_best.get("path"), os.path.join(comp_dir, "2_BEST_CREST.xyz"))
            rows.append(["CREST", crest_best["source"], f"{crest_best['energy']:.6f}", "N/A", "N/A"])

        orca_candidates = records["ORCA"]
        true_minima = [item for item in orca_candidates if item.get("imag", 0) == 0]
        rank_pool = true_minima or orca_candidates
        orca_best = min(rank_pool, key=lambda item: item["gibbs"] if item.get("gibbs", 0.0) else item["energy"], default=None)
        if orca_best:
            self._copy_unique_result(orca_best.get("path"), os.path.join(comp_dir, "3_BEST_ORCA.xyz"))
            rows.append(["ORCA", orca_best["source"], f"{orca_best['energy']:.6f}", f"{orca_best.get('gibbs', 0.0):.6f}", "N/A"])

        csv_path = os.path.join(comp_dir, "Energy_Comparison.csv")
        with open(csv_path, "w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(["Level_of_Theory", "Source_ID", "Total_Energy_Eh", "Gibbs_Free_Energy_Eh", "Binding_Energy_kcal_mol"])
            writer.writerows(rows)

        report_path = os.path.join(comp_dir, "Existing_Folder_Analysis_Report.md")
        with open(report_path, "w", encoding="utf-8") as handle:
            handle.write(f"# Existing Folder Analysis\n\n")
            handle.write(f"- **Analyzed folder:** `{job_dir}`\n")
            handle.write(f"- **Analysis time:** {_now()}\n")
            handle.write(f"- **xTB records found:** {len(records['xTB'])}\n")
            handle.write(f"- **CREST records found:** {len(records['CREST'])}\n")
            handle.write(f"- **ORCA records found:** {len(records['ORCA'])}\n")
            handle.write(f"- **ORCA true minima (0 imaginary frequencies):** {len(true_minima)}\n\n")
            handle.write("## Best Extracted Models\n\n")
            handle.write("| Method | Source | Energy (Eh) | Gibbs (Eh) | Notes |\n")
            handle.write("|---|---|---:|---:|---|\n")
            if xtb_best:
                handle.write(f"| xTB | {xtb_best['source']} | {xtb_best['energy']:.6f} | N/A | Best available xTB record |\n")
            if crest_best:
                handle.write(f"| CREST | {crest_best['source']} | {crest_best['energy']:.6f} | N/A | Best available CREST record |\n")
            if orca_best:
                handle.write(f"| ORCA | {orca_best['source']} | {orca_best['energy']:.6f} | {orca_best.get('gibbs', 0.0):.6f} | Ranked by Gibbs when available |\n")
            if not rows:
                handle.write("| N/A | N/A | N/A | N/A | No parseable energy records were found |\n")

        with open(os.path.join(logs_dir, "Existing_Folder_Analysis.json"), "w", encoding="utf-8") as handle:
            json.dump({"folder": job_dir, "time": _now(), "records": records}, handle, indent=2)
        file_analysis = self._generate_file_level_analysis(job_dir)
        return {"job_dir": job_dir, "rows": len(rows), "report": report_path, "csv": csv_path, "file_analysis": file_analysis}

    def _count_xyz_atoms(self, xyz_path):
        try:
            with open(xyz_path, "r", encoding="utf-8", errors="replace") as handle:
                return int(handle.readline().strip())
        except Exception:
            return None

    def _unique_analysis_path(self, folder, stem, suffix, used_names):
        clean_stem = re.sub(r"[^A-Za-z0-9_.-]+", "_", stem).strip("_") or "analysis"
        candidate = f"{clean_stem}{suffix}"
        idx = 1
        while candidate.lower() in used_names or os.path.exists(os.path.join(folder, candidate)):
            candidate = f"{clean_stem}_{idx}{suffix}"
            idx += 1
        used_names.add(candidate.lower())
        return os.path.join(folder, candidate)

    def _generate_structure_image(self, xyz_path, image_path):
        if not MATPLOTLIB_AVAILABLE or not xyz_path or not os.path.exists(xyz_path):
            return None
        try:
            image_dpi = int(self._float_var("image_dpi", 300))
            image_dpi = max(72, min(image_dpi, 1200))
            fig = Figure(figsize=(5.6, 5.2), dpi=min(image_dpi, 220))
            ax = fig.add_subplot(111, projection="3d")
            self._render_xyz_model(ax, xyz_path)
            ax.set_title(os.path.basename(xyz_path), fontsize=8, color=self._vars_graph["axis_color"].get())
            fig.savefig(image_path, dpi=image_dpi, bbox_inches="tight")
            return image_path
        except Exception:
            return None

    def _analysis_value(self, value):
        if value is None or value == "":
            return "N/A"
        if isinstance(value, float):
            return f"{value:.8f}"
        return str(value)

    def _relative_report_path(self, target, report_dir):
        try:
            return os.path.relpath(target, report_dir).replace("\\", "/")
        except Exception:
            return target.replace("\\", "/")

    def _analyze_xyz_file_record(self, xyz_path, images_dir, image_names):
        energy = self._xyz_comment_energy(xyz_path)
        image_path = self._unique_analysis_path(images_dir, Path(xyz_path).stem, ".png", image_names)
        generated_image = self._generate_structure_image(xyz_path, image_path)
        return {
            "source_file": xyz_path,
            "file_type": "XYZ",
            "status": "parsed",
            "atoms": self._count_xyz_atoms(xyz_path),
            "xyz_path": xyz_path,
            "image_path": generated_image,
            "final_single_point_energy_eh": energy,
            "zero_point_energy_eh": None,
            "total_thermal_energy_eh": None,
            "total_enthalpy_eh": None,
            "final_entropy_term_eh": None,
            "final_gibbs_free_energy_eh": None,
            "temperature_k": None,
            "pressure_atm": None,
            "imaginary_frequencies": None,
            "normal_termination": None,
            "notes": "Energy read from XYZ comment line." if energy is not None else "Structure image generated; no thermodynamic data present in XYZ.",
        }

    def _analyze_orca_out_record(self, out_path, images_dir, image_names):
        thermo = self._parse_orca_thermo(out_path)
        parsed = self._parse_orca_out(out_path)
        xyz_path = parsed.get("path")
        image_path = None
        if xyz_path and os.path.exists(xyz_path):
            image_path = self._unique_analysis_path(images_dir, Path(out_path).stem, ".png", image_names)
            image_path = self._generate_structure_image(xyz_path, image_path)
        status = "parsed" if thermo.get("final_single_point_energy_eh") is not None or xyz_path else "partial"
        notes = []
        if thermo.get("final_gibbs_free_energy_eh") is not None:
            notes.append("Thermodynamic Gibbs energy found.")
        if thermo.get("total_enthalpy_eh") is not None:
            notes.append("Thermodynamic enthalpy found.")
        if image_path:
            notes.append("Molecular image generated from final/output coordinates.")
        elif not xyz_path:
            notes.append("No final Cartesian coordinate block or matching XYZ was found.")
        return {
            "source_file": out_path,
            "file_type": "ORCA OUT",
            "status": status,
            "atoms": self._count_xyz_atoms(xyz_path) if xyz_path else None,
            "xyz_path": xyz_path,
            "image_path": image_path,
            "final_single_point_energy_eh": thermo.get("final_single_point_energy_eh") if thermo.get("final_single_point_energy_eh") is not None else parsed.get("energy"),
            "zero_point_energy_eh": thermo.get("zero_point_energy_eh"),
            "total_thermal_energy_eh": thermo.get("total_thermal_energy_eh"),
            "total_enthalpy_eh": thermo.get("total_enthalpy_eh"),
            "final_entropy_term_eh": thermo.get("final_entropy_term_eh"),
            "final_gibbs_free_energy_eh": thermo.get("final_gibbs_free_energy_eh") if thermo.get("final_gibbs_free_energy_eh") is not None else parsed.get("gibbs"),
            "temperature_k": thermo.get("temperature_k"),
            "pressure_atm": thermo.get("pressure_atm"),
            "imaginary_frequencies": thermo.get("imaginary_frequencies"),
            "normal_termination": thermo.get("normal_termination"),
            "notes": " ".join(notes) if notes else "ORCA output parsed.",
        }

    def _generate_file_level_analysis(self, analysis_dir):
        comp_dir = os.path.join(analysis_dir, "05_Top_Models_Comparison")
        logs_dir = os.path.join(analysis_dir, "logs")
        images_dir = os.path.join(analysis_dir, "06_File_Analysis_Images")
        opts = self._get_report_options() if hasattr(self, "_get_report_options") else {}
        include_thermo = opts.get("thermodynamics", True)
        include_images = opts.get("images", True)
        os.makedirs(comp_dir, exist_ok=True)
        os.makedirs(logs_dir, exist_ok=True)
        os.makedirs(images_dir, exist_ok=True)

        source_files = []
        for path in Path(analysis_dir).rglob("*"):
            if not path.is_file():
                continue
            lower_path = str(path).lower()
            if "06_file_analysis_images" in lower_path or "__pycache__" in lower_path:
                continue
            if path.name.lower().endswith("_final_from_out.xyz"):
                continue
            if path.suffix.lower() in (".out", ".xyz"):
                source_files.append(str(path))
        source_files = sorted(set(source_files))

        image_names = set()
        rows = []
        for file_path in source_files:
            if file_path.lower().endswith(".out"):
                rows.append(self._analyze_orca_out_record(file_path, images_dir, image_names))
            elif file_path.lower().endswith(".xyz"):
                rows.append(self._analyze_xyz_file_record(file_path, images_dir, image_names))

        csv_path = os.path.join(logs_dir, "File_Analysis_Summary.csv")
        fields = [
            "source_file",
            "file_type",
            "status",
            "atoms",
            "xyz_path",
            "notes",
        ]
        thermo_fields = [
            "final_single_point_energy_eh",
            "zero_point_energy_eh",
            "total_thermal_energy_eh",
            "total_enthalpy_eh",
            "final_entropy_term_eh",
            "final_gibbs_free_energy_eh",
            "temperature_k",
            "pressure_atm",
            "imaginary_frequencies",
            "normal_termination",
        ]
        if include_thermo:
            fields[4:4] = thermo_fields
        if include_images:
            fields.insert(-1, "image_path")
        else:
            for row in rows:
                row["image_path"] = ""
        with open(csv_path, "w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()
            for row in rows:
                writer.writerow({key: row.get(key, "") for key in fields})

        report_path = os.path.join(comp_dir, "File_Analysis_Report.md")
        with open(report_path, "w", encoding="utf-8") as handle:
            handle.write("# File-Level Structure and Thermodynamic Analysis\n\n")
            handle.write(f"- **Analyzed folder:** `{analysis_dir}`\n")
            handle.write(f"- **Analysis time:** {_now()}\n")
            handle.write(f"- **Files analyzed:** {len(rows)}\n")
            if include_images:
                handle.write(f"- **Images generated:** {sum(1 for row in rows if row.get('image_path'))}\n")
            if include_thermo:
                handle.write(f"- **Thermodynamic ORCA records:** {sum(1 for row in rows if row.get('final_gibbs_free_energy_eh') not in (None, 0.0) or row.get('total_enthalpy_eh') is not None)}\n")
            handle.write("\n")
            handle.write("## Summary Table\n\n")
            if include_thermo:
                handle.write("| File | Type | Atoms | E(eh) | H(eh) | G(eh) | ZPE(eh) | Imag | Normal | Notes |\n")
                handle.write("|---|---|---:|---:|---:|---:|---:|---:|---|---|\n")
            else:
                handle.write("| File | Type | Atoms | Notes |\n")
                handle.write("|---|---|---:|---|\n")
            for row in rows:
                if include_thermo:
                    handle.write(
                        f"| {os.path.basename(row['source_file'])} | {row['file_type']} | {self._analysis_value(row.get('atoms'))} | "
                        f"{self._analysis_value(row.get('final_single_point_energy_eh'))} | "
                        f"{self._analysis_value(row.get('total_enthalpy_eh'))} | "
                        f"{self._analysis_value(row.get('final_gibbs_free_energy_eh'))} | "
                        f"{self._analysis_value(row.get('zero_point_energy_eh'))} | "
                        f"{self._analysis_value(row.get('imaginary_frequencies'))} | "
                        f"{self._analysis_value(row.get('normal_termination'))} | {row.get('notes','')} |\n"
                    )
                else:
                    handle.write(f"| {os.path.basename(row['source_file'])} | {row['file_type']} | {self._analysis_value(row.get('atoms'))} | {row.get('notes','')} |\n")
            if rows and include_images:
                handle.write("\n## Generated Molecular Images\n\n")
                for row in rows:
                    if row.get("image_path"):
                        rel_image = self._relative_report_path(row["image_path"], comp_dir)
                        handle.write(f"### {os.path.basename(row['source_file'])}\n\n")
                        handle.write(f"![{os.path.basename(row['source_file'])}]({rel_image})\n\n")
            else:
                handle.write("\nNo `.out` or `.xyz` files were found for file-level analysis.\n")

        return {"csv": csv_path, "report": report_path, "images_dir": images_dir, "file_count": len(rows)}

    def _copy_loose_file(self, src, dest_dir, used_names):
        os.makedirs(dest_dir, exist_ok=True)
        name = os.path.basename(src)
        stem, suffix = os.path.splitext(name)
        candidate = name
        idx = 1
        while candidate.lower() in used_names or os.path.exists(os.path.join(dest_dir, candidate)):
            candidate = f"{stem}_{idx}{suffix}"
            idx += 1
        used_names.add(candidate.lower())
        dest = os.path.join(dest_dir, candidate)
        shutil.copy2(src, dest)
        return dest

    def _prepare_selected_files_analysis_dir(self, files, analysis_name="IAK_selected_files_analysis"):
        parents = [os.path.dirname(os.path.abspath(path)) for path in files]
        try:
            base_dir = os.path.commonpath(parents)
            if not os.path.isdir(base_dir):
                base_dir = parents[0]
        except Exception:
            base_dir = parents[0]
        analysis_dir = os.path.join(base_dir, analysis_name)
        inputs_dir = os.path.join(analysis_dir, "01_Inputs_and_Clusters")
        orca_dir = os.path.join(analysis_dir, "04_ORCA_Refinement")
        os.makedirs(inputs_dir, exist_ok=True)
        os.makedirs(os.path.join(analysis_dir, "02_xTB_Results"), exist_ok=True)
        os.makedirs(os.path.join(analysis_dir, "03_CREST_Results"), exist_ok=True)
        os.makedirs(orca_dir, exist_ok=True)
        os.makedirs(os.path.join(analysis_dir, "05_Top_Models_Comparison"), exist_ok=True)
        used = set()
        copied_xyz, copied_out = 0, 0
        for src in files:
            lower = src.lower()
            if lower.endswith(".xyz"):
                self._copy_loose_file(src, inputs_dir, used)
                copied_xyz += 1
            elif lower.endswith(".out"):
                self._copy_loose_file(src, orca_dir, used)
                copied_out += 1
        return analysis_dir, copied_xyz, copied_out

    def _create_loose_folder_analysis(self, folder):
        analysis_dir = os.path.join(folder, "IAK_existing_analysis")
        inputs_dir = os.path.join(analysis_dir, "01_Inputs_and_Clusters")
        orca_dir = os.path.join(analysis_dir, "04_ORCA_Refinement")
        os.makedirs(inputs_dir, exist_ok=True)
        os.makedirs(os.path.join(analysis_dir, "02_xTB_Results"), exist_ok=True)
        os.makedirs(os.path.join(analysis_dir, "03_CREST_Results"), exist_ok=True)
        os.makedirs(orca_dir, exist_ok=True)
        os.makedirs(os.path.join(analysis_dir, "05_Top_Models_Comparison"), exist_ok=True)
        used = set()
        copied_xyz, copied_out = 0, 0
        for root, dirs, files in os.walk(folder):
            if os.path.abspath(root).startswith(os.path.abspath(analysis_dir)):
                continue
            for name in files:
                src = os.path.join(root, name)
                lower = name.lower()
                if lower.endswith(".xyz") and copied_xyz < 300:
                    self._copy_loose_file(src, inputs_dir, used)
                    copied_xyz += 1
                elif lower.endswith(".out") and copied_out < 300:
                    self._copy_loose_file(src, orca_dir, used)
                    copied_out += 1
                elif lower == "energy_comparison.csv":
                    self._copy_loose_file(src, os.path.join(analysis_dir, "05_Top_Models_Comparison"), used)
        return analysis_dir

    def _analyze_existing_folder(self):
        folder = filedialog.askdirectory(title="Select existing IAK result folder or folder containing XYZ/ORCA files")
        if not folder:
            return
        folder = os.path.abspath(folder)
        self._reset_progress_ui()
        self._apply_progress_payload({"stage": "Validation", "percent": 5.0, "status": "running", "message": "Scanning existing folder.", "elapsed_seconds": 0.0})
        self._append_text(f"\n[Existing Analysis] Scanning folder: {folder}\n")
        try:
            job_dirs = self._find_iak_job_dirs(folder)
            if not job_dirs:
                self._append_text("[Existing Analysis] No IAK folder structure found. Creating loose-file analysis folder...\n")
                job_dirs = [self._create_loose_folder_analysis(folder)]

            results = []
            total = max(len(job_dirs), 1)
            started = time.time()
            for idx, job_dir in enumerate(job_dirs, start=1):
                percent = 10.0 + 80.0 * ((idx - 1) / total)
                self._apply_progress_payload({
                    "stage": "Reports",
                    "percent": percent,
                    "status": "running",
                    "message": f"Analyzing existing job {idx}/{total}.",
                    "elapsed_seconds": time.time() - started,
                })
                results.append(self._write_existing_analysis_outputs(job_dir))

            first_job = results[0]["job_dir"]
            self._vars["out"].set(first_job)
            if MATPLOTLIB_AVAILABLE:
                self._vars_graph["dir1"].set(folder if len(results) > 1 else first_job)
                self._vars_graph["name1"].set(os.path.basename(folder.rstrip("\\/")) or "Existing Data")
                if hasattr(self, "thermo_tree"):
                    self._refresh_thermo_table()
            self._refresh_results()
            self.nb.select(self.tab_res)
            self._apply_progress_payload({
                "stage": "Complete",
                "percent": 100.0,
                "status": "completed",
                "message": f"Analyzed {len(results)} existing folder(s).",
                "elapsed_seconds": time.time() - started,
            })
            self._append_text(f"[Existing Analysis] Complete. Loaded results from: {first_job}\n")
            messagebox.showinfo(
                "Existing Folder Analysis Complete",
                f"Analyzed {len(results)} folder(s).\n\nLoaded first result folder:\n{first_job}\n\nReports and Energy_Comparison.csv files were generated where possible.",
            )
        except Exception as exc:
            self._append_text(f"[Existing Analysis Error] {exc}\n")
            self._apply_progress_payload({"stage": "Failed", "percent": 100.0, "status": "failed", "message": str(exc), "elapsed_seconds": 0.0})
            messagebox.showerror("Existing Folder Analysis Error", str(exc))

    def _analyze_existing_files(self):
        files = filedialog.askopenfilenames(
            title="Select existing ORCA .out and/or XYZ files",
            filetypes=[
                ("ORCA output and XYZ", "*.out *.xyz"),
                ("ORCA output", "*.out"),
                ("XYZ structures", "*.xyz"),
                ("All files", "*.*"),
            ],
        )
        files = [os.path.abspath(path) for path in files if path.lower().endswith((".out", ".xyz"))]
        if not files:
            return
        self._reset_progress_ui()
        started = time.time()
        self._apply_progress_payload({
            "stage": "Validation",
            "percent": 5.0,
            "status": "running",
            "message": f"Preparing {len(files)} selected file(s).",
            "elapsed_seconds": 0.0,
        })
        try:
            analysis_dir, copied_xyz, copied_out = self._prepare_selected_files_analysis_dir(files)
            self._append_text(
                f"\n[Existing File Analysis] Selected {len(files)} file(s): "
                f"{copied_xyz} XYZ and {copied_out} ORCA OUT copied into {analysis_dir}\n"
            )
            self._apply_progress_payload({
                "stage": "Reports",
                "percent": 55.0,
                "status": "running",
                "message": "Parsing selected files and writing analysis outputs.",
                "elapsed_seconds": time.time() - started,
            })
            result = self._write_existing_analysis_outputs(analysis_dir)
            self._vars["out"].set(analysis_dir)
            if MATPLOTLIB_AVAILABLE:
                self._vars_graph["dir1"].set(analysis_dir)
                self._vars_graph["name1"].set("Selected Files")
                if hasattr(self, "thermo_tree"):
                    self._refresh_thermo_table()
            self._refresh_results()
            self.nb.select(self.tab_res)
            self._apply_progress_payload({
                "stage": "Complete",
                "percent": 100.0,
                "status": "completed",
                "message": "Selected file analysis complete.",
                "elapsed_seconds": time.time() - started,
            })
            messagebox.showinfo(
                "Existing File Analysis Complete",
                "Selected files were analyzed.\n\n"
                f"Analysis folder:\n{analysis_dir}\n\n"
                f"XYZ files: {copied_xyz}\nORCA .out files: {copied_out}\n"
                f"Energy rows parsed: {result['rows']}",
            )
        except Exception as exc:
            self._append_text(f"[Existing File Analysis Error] {exc}\n")
            self._apply_progress_payload({
                "stage": "Failed",
                "percent": 100.0,
                "status": "failed",
                "message": str(exc),
                "elapsed_seconds": time.time() - started,
            })
            messagebox.showerror("Existing File Analysis Error", str(exc))

    def _resume_existing_job(self):
        folder = filedialog.askdirectory(title="Select stopped/failed IAK job folder or parent folder")
        if not folder:
            return
        folder = os.path.abspath(folder)
        job_dirs = self._find_iak_job_dirs(folder)
        if not job_dirs:
            return messagebox.showerror("Resume Job", "No IAK job folder was found in the selected location.")
        job_dir = job_dirs[0]
        inputs_dir = os.path.join(job_dir, "01_Inputs_and_Clusters")
        anchor = os.path.join(inputs_dir, "raw_input_anchor.xyz")
        guest = os.path.join(inputs_dir, "raw_input_guest.xyz")
        if not os.path.exists(anchor):
            candidates = list(Path(inputs_dir).glob("*.xyz")) if os.path.isdir(inputs_dir) else []
            anchor = str(candidates[0]) if candidates else ""
        if anchor and os.path.exists(anchor):
            self._vars["a"].set(anchor)
        if os.path.exists(guest):
            self._vars["b"].set(guest)
        else:
            self._vars["b"].set("")
        folder_name = os.path.basename(job_dir.rstrip("\\/"))
        match = re.search(r"(\d+)[_:](\d+)$", folder_name)
        if match:
            self._vars["ratio"].set(f"{match.group(1)}:{match.group(2)}")
        self._vars["out"].set(job_dir)
        self._refresh_results()
        self.nb.select(self.tab_main)
        self._append_text(
            f"\n[Resume] Loaded existing job folder:\n{job_dir}\n"
            "Press START BATCH PIPELINE to continue. Successful stages in state.json will be skipped automatically.\n"
        )
        messagebox.showinfo(
            "Resume Ready",
            "Existing job loaded.\n\nPress START BATCH PIPELINE to continue from the saved state.\nAlready successful xTB/CREST/ORCA jobs will be skipped.",
        )

    def _update_guest_count_label(self):
        n = len(self._guest_list)
        color = _C["green"] if n > 0 else _C["muted"]
        self._guest_count_lbl.config(
            text=f"{n} guest(s) queued" if n != 1 else "1 guest queued",
            fg=color)

    def _start(self):
        v = {}
        for k, var in self._vars.items():
            val = var.get().strip(" \"'")
            if k in ["a", "b"] and val:
                val = os.path.abspath(val)
            v[k] = val
        inject_embedded_engines()
        self._update_installation_labels()
        if not os.path.isfile(v["a"]):
            self._append_text(f"\n[Error] Anchor XYZ file not found or invalid path: {v['a']}\n")
            return messagebox.showerror("Error", f"Anchor XYZ file not found:\n{v['a']}")

        # Resolve effective guest list: multi-guest panel takes priority over legacy field
        effective_guests = [p for p in self._guest_list if os.path.isfile(p)]
        if not effective_guests and v["b"] and os.path.isfile(v["b"]):
            effective_guests = [v["b"]]

        if not effective_guests:
            if "0" not in v["ratio"]:
                ans = messagebox.askyesno("Single Molecule Mode",
                    "No Guest (B) file found in queue or legacy field.\n"
                    "Run as single isolated molecule (Ratio 1:0)?")
                if ans:
                    v["ratio"] = "1:0"
                else:
                    return
        v["_guests"] = effective_guests  # pass resolved list to worker
        self._reset_progress_ui()
        self.go_btn.config(state="disabled", text="RUNNING BATCH QUEUE...")
        self.nb.select(0)
        g_count = len(effective_guests)
        self._append_text(
            f"\n[System] Initialization complete. "
            f"{g_count} guest molecule(s) queued. "
            f"Booting Batch Pipeline Thread...\n")
        self.start_time = time.time()
        self.is_running = True
        threading.Thread(target=self._worker, args=(v,), daemon=True).start()

    def _worker(self, v):
        try:
            mode_name = v.get("mode", "balanced").lower()
            base_mode = RunMode.BALANCED if mode_name == "custom" else RunMode[mode_name.upper()]
            config = Config.from_mode(base_mode)
            config.preopt_inputs = self.run_preopt.get()
            try:
                config.cores = int(v.get("cores", 4))
                config.maxcore = int(v.get("maxcore", 2000))
                config.charge = int(v.get("charge", 0))
                config.multiplicity = int(v.get("mult", 1))
                config.n_generate = int(v.get("n_generate", config.n_generate))
                config.n_keep_scored = int(v.get("n_keep_scored", config.n_keep_scored))
                config.n_keep_clustered = int(v.get("n_keep_clustered", config.n_keep_clustered))
                config.n_run_xtb = int(v.get("n_run_xtb", config.n_run_xtb))
                config.n_run_crest = int(v.get("n_run_crest", config.n_run_crest))
                config.random_seed = int(v.get("random_seed", config.random_seed))
                config.max_placement_attempts = int(v.get("max_placement_attempts", config.max_placement_attempts))
                config.rmsd_cutoff = float(v.get("rmsd_cutoff", config.rmsd_cutoff))
                config.xtb_ewin_kcal = float(v.get("xtb_ewin_kcal", config.xtb_ewin_kcal))
                config.crest_ewin_kcal = float(v.get("crest_ewin_kcal", config.crest_ewin_kcal))
            except ValueError:
                raise RuntimeError("Host-controlled numerical inputs must be valid numbers.")
            if config.n_generate < 1 or config.n_keep_scored < 1 or config.n_keep_clustered < 1:
                raise RuntimeError("Generate/keep counts must be at least 1.")
            if config.n_run_xtb < 0 or config.n_run_crest < 0:
                raise RuntimeError("xTB/CREST run counts cannot be negative.")
            if config.max_placement_attempts < 1:
                raise RuntimeError("Placement attempts must be at least 1.")
            config.xtb_method = v.get("xtb_method", config.xtb_method).strip() or "--gfn2"
            config.crest_method = v.get("crest_method", config.crest_method).strip() or "--gfn2"
            config.orca_method = v.get("orca_method", "B97-3c Opt Freq")
            config.report_options = self._get_report_options()
            reaction_type = normalize_reaction_type(v.get("reaction_type", "Non-covalent"))
            try:
                if v.get("e_a"):
                    config.energy_a = float(v.get("e_a"))
                if v.get("e_b"):
                    config.energy_b = float(v.get("e_b"))
            except ValueError:
                raise RuntimeError("Monomer Energies must be valid decimals. Leave empty if not calculating Binding Energy.")
            ratio_pairs = []
            for r in v["ratio"].split(","):
                r = r.strip()
                if not r:
                    continue
                if ":" in r:
                    parts = tuple(int(p) for p in r.split(":"))
                    ratio_pairs.append(parts)
                else:
                    ratio_pairs.append((1, int(r)))
            if not ratio_pairs:
                raise RuntimeError("Please provide at least one valid ratio (e.g., 1:0, 1:1).")

            # Resolve guest files: multi-guest list (may be empty → single-molecule mode)
            guest_paths = v.get("_guests") or ([] if not v.get("b") else [v["b"]])
            if not guest_paths:
                guest_paths = [None]  # sentinel → single-molecule

            base_out_dir = v["out"]
            prefix = re.sub(r"\d+[_:\-]*\d*$", "", base_out_dir)
            if prefix == base_out_dir and not prefix.endswith("_"):
                prefix += "_"

            total_jobs = len(ratio_pairs)
            job_idx = 0
            for r_tuple in ratio_pairs:
                job_idx += 1
                n_A = r_tuple[0]
                n_guests_list = list(r_tuple[1:]) if len(r_tuple) > 1 else [0]
                if n_A < 0 or any(x < 0 for x in n_guests_list):
                    raise RuntimeError(f"Negative stoichiometry is not allowed: {r_tuple}")
                if n_A == 0 and sum(n_guests_list) == 0:
                    raise RuntimeError(f"Invalid ratio {r_tuple}: at least one species count must be > 0.")
                
                # Check tuple length against guest paths
                if len(n_guests_list) > len(guest_paths) and guest_paths != [None]:
                    self.root.after(0, lambda r=r_tuple: self._append_text(
                        f"\n[WARNING] Ratio {r} has more guest components ({len(n_guests_list)}) "
                        f"than queued guests ({len(guest_paths)}). Extra components ignored.\n"))
                    n_guests_list = n_guests_list[:len(guest_paths)]
                elif len(n_guests_list) < len(guest_paths) and guest_paths != [None]:
                    n_guests_list.extend([0] * (len(guest_paths) - len(n_guests_list)))

                # If ratio signifies single molecule (either no guest or 0 copies of all guests)
                is_single = (sum(n_guests_list) == 0)

                try:
                    current_out_name = f"{prefix}{n_A}_" + "_".join(map(str, n_guests_list))
                    current_out = os.path.abspath(current_out_name)
                    ratio_label = f"{n_A}_" + "_".join(map(str, n_guests_list))
                    self.root.after(0, lambda o=current_out_name: self._vars["out"].set(o))
                    self.root.after(0, lambda r=r_tuple, o=current_out_name, ji=job_idx, jt=total_jobs:
                        self._append_text(
                            f"\n{'=' * 70}\n"
                            f"[BATCH QUEUE] Job {ji}/{jt} | Ratio {r} | Folder: {o}\n"
                            f"{'=' * 70}\n"))
                    
                    eff_a_path = v["a"] if n_A > 0 else None
                    eff_b_paths = guest_paths if not is_single else []

                    # Validate and auto-correct spin multiplicity based on total electrons.
                    anchor_mol = Molecule.from_xyz(eff_a_path) if eff_a_path and os.path.exists(eff_a_path) else Molecule([])
                    guest_mols = [Molecule.from_xyz(p) for p in eff_b_paths if p and os.path.exists(p)]
                    job_config = dataclasses.replace(config)
                    cm_check = validate_charge_multiplicity(
                        anchor_mol,
                        guest_mols,
                        n_A,
                        n_guests_list,
                        job_config.charge,
                        job_config.multiplicity,
                    )
                    if not cm_check["valid"]:
                        total_electrons = cm_check.get("total_electrons")
                        suggested = int(cm_check.get("suggested_multiplicity", job_config.multiplicity))
                        if total_electrons is not None and total_electrons > 0 and suggested >= 1:
                            old_mult = job_config.multiplicity
                            job_config.multiplicity = suggested
                            self.root.after(
                                0,
                                lambda r=r_tuple, m=cm_check["message"], om=old_mult, sm=suggested: self._append_text(
                                    f"\n[SPIN AUTO-CORRECTION] Ratio {r}: {m} "
                                    f"Using multiplicity {om} -> {sm}.\n"
                                ),
                            )
                        else:
                            raise RuntimeError(f"Ratio {r_tuple}: {cm_check['message']}")
                    else:
                        self.root.after(
                            0,
                            lambda r=r_tuple, te=cm_check.get("total_electrons"): self._append_text(
                                f"[Validation] Ratio {r}: total electrons = {te}; charge/multiplicity validated.\n"
                            ),
                        )
                    
                    pipe = Pipeline(
                        eff_a_path,
                        eff_b_paths,
                        n_guests_list,
                        job_config,
                        current_out,
                        ratio_label=ratio_label,
                        progress_cb=self._progress_cb,
                        reaction_type=reaction_type,
                        n_anchor=n_A,
                    )
                    pipe.run(run_xtb=self.run_xtb.get(), run_crest=self.run_crest.get(),
                             run_orca=self.run_orca.get(), log_cb=self._append_text,
                             status_cb=self._status_cb)
                    self.root.after(0, self._refresh_results)
                except Exception as iter_e:
                    self.root.after(0, lambda e=iter_e, r=r_tuple:
                        self._append_text(
                            f"\n[CRITICAL BATCH ERROR] Ratio {r} aborted: "
                            f"{str(e)}. Proceeding to next job...\n"))
            self.root.after(0, lambda jt=total_jobs: messagebox.showinfo(
                "Batch Complete",
                f"All {jt} job(s) completed!\n\n"
                f"Guests processed: {len(guest_paths) if guest_paths != [None] else 0}\n"
                f"Ratios processed: {len(ratio_pairs)}\n\n"
                "Check the respective folders for Top Models, Energy CSVs, and Job Summary Reports."))
        except Exception as e:
            self.root.after(0, lambda e=e: self._append_text(f"\n[CRITICAL ERROR] Pipeline aborted: {str(e)}\n"))
            self.root.after(0, lambda e=e: messagebox.showerror("Error", str(e)))
        finally:
            self.is_running = False
            self.root.after(0, lambda: self.go_btn.config(state="normal", text="START BATCH PIPELINE"))

    def _pes_log(self, text):
        self.root.after(0, lambda: self.pes_term.config(state="normal"))
        self.root.after(0, lambda: self.pes_term.insert("end", text + "\n"))
        self.root.after(0, lambda: self.pes_term.see("end"))
        self.root.after(0, lambda: self.pes_term.config(state="disabled"))

    def _start_pes(self):
        v = {k: var.get() for k, var in self._vars_pes.items()}
        if not v["reactant"] or not os.path.exists(v["reactant"]):
            return messagebox.showerror("Error", "Valid Reactant XYZ is required for PES Scan.")
        self.btn_run_pes.config(state="disabled", text="RUNNING PES/TS...")
        self._pes_log("="*60)
        self._pes_log(" INITIALIZING PES & TS SEARCH PIPELINE ")
        self._pes_log("="*60)
        threading.Thread(target=self._worker_pes, args=(v,), daemon=True).start()

    def _worker_pes(self, v):
        try:
            # 1. Parse Input Geometry
            react_path = os.path.abspath(v["reactant"])
            react_name = os.path.splitext(os.path.basename(react_path))[0]
            work_dir = os.path.abspath(f"PES_{react_name}_{int(time.time())}")
            os.makedirs(work_dir, exist_ok=True)
            shutil.copy2(react_path, os.path.join(work_dir, "reactant.xyz"))
            reaction_type = normalize_reaction_type(self._vars.get("reaction_type").get() if "reaction_type" in self._vars else "Non-covalent")

            charge = _safe_int(self._vars.get("charge").get() if "charge" in self._vars else 0)
            multiplicity = _safe_int(self._vars.get("mult").get() if "mult" in self._vars else 1)
            if charge is None:
                raise RuntimeError("System charge must be an integer for PES/TS calculations.")
            if multiplicity is None or multiplicity < 1:
                raise RuntimeError("Multiplicity must be a positive integer for PES/TS calculations.")

            reactant_mol = Molecule.from_xyz(react_path)
            pes_cm_check = validate_charge_multiplicity(reactant_mol, [], 1, [], charge, multiplicity)
            if not pes_cm_check["valid"]:
                total_electrons = pes_cm_check.get("total_electrons")
                suggested = int(pes_cm_check.get("suggested_multiplicity", multiplicity))
                if total_electrons is not None and total_electrons > 0 and suggested >= 1:
                    self._pes_log(
                        f"[SPIN AUTO-CORRECTION] {pes_cm_check['message']} "
                        f"Using multiplicity {multiplicity} -> {suggested}."
                    )
                    multiplicity = suggested
                else:
                    raise RuntimeError(pes_cm_check["message"])
            else:
                self._pes_log(f"[Validation] Total electrons = {pes_cm_check.get('total_electrons')}; charge/multiplicity validated.")
            
            scan_mode = v.get("scan_mode", "coord")
            engine = v["engine"]
            
            if scan_mode == "path":
                if not v["product"] or not os.path.exists(v["product"]):
                    raise RuntimeError("Product XYZ is required for Reactant -> Product Path Interpolation.")
                prod_path = os.path.abspath(v["product"])
                shutil.copy2(prod_path, os.path.join(work_dir, "product.xyz"))
                product_mol = Molecule.from_xyz(prod_path)
                balanced, balance_message = validate_reactant_product_atom_balance(reactant_mol, product_mol, reaction_type)
                if not balanced:
                    raise RuntimeError(balance_message)
                self._pes_log(f"[Validation] {balance_message}")
                
                self._pes_log(f"Working Directory: {work_dir}")
                self._pes_log(f"Engine: {engine.upper()}")
                self._pes_log(f"Mode: Reactant -> Product Path Interpolation ({reaction_type})")
                
                if engine == "xtb":
                    self._run_xtb_path(work_dir, "reactant.xyz", "product.xyz", charge, multiplicity)
                elif engine == "orca":
                    self._run_orca_neb(work_dir, "reactant.xyz", "product.xyz", charge, multiplicity)
                else:
                    raise RuntimeError(f"Engine {engine} not supported for Path Interpolation.")
            else:
                try:
                    c1_a1, c1_a2 = int(v["c1_a1"]), int(v["c1_a2"])
                    c1_s, c1_e, c1_steps = float(v["c1_start"]), float(v["c1_end"]), int(v["c1_steps"])
                except ValueError:
                    raise RuntimeError("Invalid values for Coordinate 1. Must be numerical.")
                
                use_c2 = int(v["use_c2"]) == 1
                if use_c2:
                    try:
                        c2_a1, c2_a2 = int(v["c2_a1"]), int(v["c2_a2"])
                        c2_s, c2_e, c2_steps = float(v["c2_start"]), float(v["c2_end"]), int(v["c2_steps"])
                    except ValueError:
                        raise RuntimeError("Invalid values for Coordinate 2. Must be numerical.")
                else:
                    c2_a1=c2_a2=c2_s=c2_e=c2_steps=None
                
                self._pes_log(f"Working Directory: {work_dir}")
                self._pes_log(f"Engine: {engine.upper()}")
                self._pes_log("Mode: Coordinate Grid Scan")
                
                if engine == "xtb":
                    self._run_xtb_pes(work_dir, "reactant.xyz", c1_a1, c1_a2, c1_s, c1_e, c1_steps, use_c2, c2_a1, c2_a2, c2_s, c2_e, c2_steps, charge, multiplicity)
                elif engine == "orca":
                    self._run_orca_pes(work_dir, "reactant.xyz", c1_a1, c1_a2, c1_s, c1_e, c1_steps, use_c2, c2_a1, c2_a2, c2_s, c2_e, c2_steps, charge, multiplicity)
                else:
                    raise RuntimeError(f"Engine {engine} not supported for PES.")
                
            # 4. Extract TS Guess and Run OptTS if checked
            run_ts = int(v["run_ts"]) == 1
            if run_ts:
                self._pes_log("\n--- Proceeding to Transition State Refinement ---")
                ts_guess_path = os.path.join(work_dir, "ts_guess.xyz")
                if not os.path.exists(ts_guess_path):
                    raise RuntimeError("Saddle point (ts_guess.xyz) not found from PES scan. Cannot run TS Opt.")
                if engine == "xtb":
                    self._run_xtb_ts(work_dir, ts_guess_path, charge, multiplicity)
                else:
                    self._run_orca_ts(work_dir, ts_guess_path, charge, multiplicity)
                    
            self._pes_log("\n[SUCCESS] PES/TS Workflow Completed!")
            self.root.after(0, lambda: messagebox.showinfo("PES/TS Complete", f"Workflow finished.\nOutputs saved in:\n{work_dir}"))
        except Exception as e:
            self._pes_log(f"\n[ERROR] PES Pipeline failed: {str(e)}")
            self.root.after(0, lambda: messagebox.showerror("PES Error", str(e)))
        finally:
            self.root.after(0, lambda: self.btn_run_pes.config(state="normal", text="START PES & TS SEARCH"))

    def _run_xtb_pes(self, work_dir, mol_file, c1_a1, c1_a2, c1_s, c1_e, c1_steps, use_c2, c2_a1, c2_a2, c2_s, c2_e, c2_steps, charge, multiplicity):
        scan_inp = os.path.join(work_dir, "scan.inp")
        with open(scan_inp, "w") as f:
            f.write("$scan\n")
            f.write(f"  {c1_a1} {c1_a2} : {c1_s}, {c1_e}, {c1_steps}\n")
            if use_c2:
                f.write(f"  {c2_a1} {c2_a2} : {c2_s}, {c2_e}, {c2_steps}\n")
            f.write("$end\n")
        self._pes_log("Running xTB relaxed scan... (This may take a while)")
        
        from platform import system as _os_sys
        xtb_exe = "xtb.exe" if _os_sys() == "Windows" else "xtb"
        if not is_tool_available("xtb"):
             raise RuntimeError("xTB not found on path.")
             
        cmd = f"\"{xtb_exe}\" {mol_file} --opt --input scan.inp"
        if charge != 0:
            cmd += f" --chrg {charge}"
        if multiplicity != 1:
            cmd += f" --uhf {multiplicity - 1}"
        p = subprocess.run(cmd, shell=True, cwd=work_dir, capture_output=True, text=True)
        with open(os.path.join(work_dir, "xtb_pes.log"), "w") as f:
            f.write(p.stdout)
            f.write(p.stderr)
            
        self._pes_log("xTB Scan finished. Parsing xtbscan.log...")
        scan_log = os.path.join(work_dir, "xtbscan.log")
        if not os.path.exists(scan_log):
             raise RuntimeError("xtbscan.log not found. The xTB scan failed.")
             
        energies = []
        geoms = []
        with open(scan_log, "r") as f:
             lines = f.readlines()
        
        i = 0
        while i < len(lines):
             line = lines[i].strip()
             if not line:
                 i += 1; continue
             try:
                 n_atoms = int(line)
                 if "energy:" in lines[i+1]:
                     en = float(re.search(r"energy:\s*([-0-9.]+)", lines[i+1]).group(1))
                 else:
                     en = 0.0
                 energies.append(en)
                 block = lines[i:i+2+n_atoms]
                 geoms.append("".join(block))
                 i += 2 + n_atoms
             except Exception:
                 i += 1
                 
        if not energies:
             raise RuntimeError("No step energies found in xtbscan.log.")
             
        self._pes_log(f"Extracted {len(energies)} points from scan.")
        max_idx = np.argmax(energies)
        self._pes_log(f"Highest energy point found at Step {max_idx+1} (E = {energies[max_idx]} Eh). Saving as ts_guess.xyz")
        
        with open(os.path.join(work_dir, "ts_guess.xyz"), "w") as f:
            f.write(geoms[max_idx])
            
        if use_c2 and MATPLOTLIB_AVAILABLE:
            self._generate_3d_pes_plot(work_dir, c1_steps, c2_steps, energies)

    def _run_orca_pes(self, work_dir, mol_file, c1_a1, c1_a2, c1_s, c1_e, c1_steps, use_c2, c2_a1, c2_a2, c2_s, c2_e, c2_steps, charge, multiplicity):
        self._pes_log("Running ORCA relaxed scan... (This may take a very long time)")
        mol = Molecule.from_xyz(os.path.join(work_dir, mol_file))
        inp_path = os.path.join(work_dir, "orca_pes.inp")
        
        # NOTE: ORCA atom indices are 0-based for scanning!
        with open(inp_path, "w") as f:
            f.write("! B97-3c Opt\n")
            f.write("%geom Scan\n")
            f.write(f"  B {c1_a1-1} {c1_a2-1} = {c1_s}, {c1_e}, {c1_steps}\n")
            if use_c2:
                f.write(f"  B {c2_a1-1} {c2_a2-1} = {c2_s}, {c2_e}, {c2_steps}\n")
            f.write("end end\n")
            f.write(f"* xyz {charge} {multiplicity}\n")
            for a in mol.atoms:
                f.write(f"{a.symbol} {a.x} {a.y} {a.z}\n")
            f.write("*\n")
            
        cmd = f"orca orca_pes.inp > orca_pes.out"
        subprocess.run(cmd, shell=True, cwd=work_dir)
        
        out_path = os.path.join(work_dir, "orca_pes.out")
        if not os.path.exists(out_path):
             raise RuntimeError("orca_pes.out not found. The ORCA scan failed.")
             
        trj_path = os.path.join(work_dir, "orca_pes.trj")
        if not os.path.exists(trj_path):
             self._pes_log("Warning: orca_pes.trj not found. Attempting to parse ORCA output for TS guess manually.")
             # Very basic fallback
             return

        # Parse orca.trj
        energies = []
        geoms = []
        with open(trj_path, "r") as f:
             lines = f.readlines()
        i = 0
        while i < len(lines):
             line = lines[i].strip()
             if not line:
                 i += 1; continue
             try:
                 n_atoms = int(line)
                 en_match = re.search(r"[-0-9.]+", lines[i+1])
                 en = float(en_match.group(0)) if en_match else 0.0
                 energies.append(en)
                 block = lines[i:i+2+n_atoms]
                 geoms.append("".join(block))
                 i += 2 + n_atoms
             except Exception:
                 i += 1
                 
        if not energies:
             raise RuntimeError("No step energies found in orca_pes.trj.")
             
        self._pes_log(f"Extracted {len(energies)} points from ORCA scan.")
        max_idx = np.argmax(energies)
        self._pes_log(f"Highest energy point found at Step {max_idx+1} (E = {energies[max_idx]} Eh). Saving as ts_guess.xyz")
        
        with open(os.path.join(work_dir, "ts_guess.xyz"), "w") as f:
            f.write(geoms[max_idx])
            
        if use_c2 and MATPLOTLIB_AVAILABLE:
            self._generate_3d_pes_plot(work_dir, c1_steps, c2_steps, energies)

    def _run_xtb_ts(self, work_dir, ts_guess_path, charge, multiplicity):
        self._pes_log("Running xTB Transition State Optimization (--opt ts)...")
        from platform import system as _os_sys
        xtb_exe = "xtb.exe" if _os_sys() == "Windows" else "xtb"
        ts_mol = os.path.basename(ts_guess_path)
        cmd = f"\"{xtb_exe}\" {ts_mol} --opt ts --hess"
        if charge != 0:
            cmd += f" --chrg {charge}"
        if multiplicity != 1:
            cmd += f" --uhf {multiplicity - 1}"
        p = subprocess.run(cmd, shell=True, cwd=work_dir, capture_output=True, text=True)
        with open(os.path.join(work_dir, "xtb_ts.log"), "w") as f:
            f.write(p.stdout)
            f.write(p.stderr)
        
        imag_freqs = []
        for line in p.stdout.splitlines():
            if re.match(r"^\s*\d+\s+-\d+\.\d+", line):
                parts = line.split()
                if len(parts) >= 2 and float(parts[1]) < 0:
                     imag_freqs.append(float(parts[1]))
                     
        self._pes_log(f"xTB TS Output: Found {len(imag_freqs)} imaginary frequencies.")
        if len(imag_freqs) == 1:
             self._pes_log(f"SUCCESS: True TS confirmed! Imaginary mode: {imag_freqs[0]} cm-1")
        else:
             self._pes_log(f"WARNING: Found {len(imag_freqs)} imaginary modes. This might not be a valid TS.")

    def _run_orca_ts(self, work_dir, ts_guess_path, charge, multiplicity):
        self._pes_log("Generating ORCA OptTS input...")
        mol = Molecule.from_xyz(ts_guess_path)
        inp_path = os.path.join(work_dir, "orca_ts.inp")
        with open(inp_path, "w") as f:
            f.write("! B97-3c OptTS Freq\n")
            f.write(f"* xyz {charge} {multiplicity}\n")
            for a in mol.atoms:
                f.write(f"{a.symbol} {a.x} {a.y} {a.z}\n")
            f.write("*\n")
            
        cmd = f"orca orca_ts.inp > orca_ts.out"
        self._pes_log("Running ORCA OptTS... (this may take a long time)")
        subprocess.run(cmd, shell=True, cwd=work_dir)
        
        out_path = os.path.join(work_dir, "orca_ts.out")
        if not os.path.exists(out_path):
            raise RuntimeError("ORCA TS output not found.")
            
        with open(out_path, "r") as f:
            content = f.read()
            
        imag_freqs = re.findall(r"^\s*\d+:\s+-\d+\.\d+\s+cm\*\*-1", content, re.MULTILINE)
        self._pes_log(f"ORCA TS Output: Found {len(imag_freqs)} imaginary frequencies.")
        if len(imag_freqs) == 1:
             self._pes_log(f"SUCCESS: True TS confirmed! Imaginary mode: {imag_freqs[0].strip()}")
        else:
             self._pes_log(f"WARNING: Found {len(imag_freqs)} imaginary modes. This might not be a valid TS.")

    def _run_xtb_path(self, work_dir, reactant_file, product_file, charge, multiplicity):
        self._pes_log("Running xTB Reactant -> Product Path Interpolation (--path)...")
        from platform import system as _os_sys
        xtb_exe = "xtb.exe" if _os_sys() == "Windows" else "xtb"
        if not is_tool_available("xtb"):
             raise RuntimeError("xTB not found on path.")
             
        cmd = f"\"{xtb_exe}\" {reactant_file} --path {product_file}"
        if charge != 0:
            cmd += f" --chrg {charge}"
        if multiplicity != 1:
            cmd += f" --uhf {multiplicity - 1}"
        p = subprocess.run(cmd, shell=True, cwd=work_dir, capture_output=True, text=True)
        with open(os.path.join(work_dir, "xtb_path.log"), "w") as f:
            f.write(p.stdout)
            f.write(p.stderr)
            
        self._pes_log("xTB Path Interpolation finished. Parsing xtbpath.xyz...")
        path_log = os.path.join(work_dir, "xtbpath.xyz")
        if not os.path.exists(path_log):
             raise RuntimeError("xtbpath.xyz not found. The xTB path interpolation failed.")
             
        # Extract highest energy point (or just take point #5 if arbitrary)
        # xTB NEB logs energy in the comment line of each frame in xtbpath.xyz
        energies = []
        geoms = []
        with open(path_log, "r") as f:
             lines = f.readlines()
        
        i = 0
        while i < len(lines):
             line = lines[i].strip()
             if not line:
                 i += 1; continue
             try:
                 n_atoms = int(line)
                 en_match = re.search(r"energy:\s*([-0-9.]+)", lines[i+1].lower())
                 en = float(en_match.group(1)) if en_match else 0.0
                 energies.append(en)
                 block = lines[i:i+2+n_atoms]
                 geoms.append("".join(block))
                 i += 2 + n_atoms
             except Exception:
                 i += 1
                 
        if not energies:
             raise RuntimeError("No structures parsed from xtbpath.xyz.")
             
        self._pes_log(f"Extracted {len(energies)} frames from path.")
        max_idx = np.argmax(energies)
        self._pes_log(f"Highest energy frame found at Step {max_idx+1} (E = {energies[max_idx]} Eh). Saving as ts_guess.xyz")
        
        with open(os.path.join(work_dir, "ts_guess.xyz"), "w") as f:
            f.write(geoms[max_idx])
            
        # Optional: Generate a simple 1D plot
        if MATPLOTLIB_AVAILABLE:
            self._generate_1d_pes_plot(work_dir, energies)

    def _run_orca_neb(self, work_dir, reactant_file, product_file, charge, multiplicity):
        self._pes_log("Running ORCA NEB-TS Interpolation... (This will take a very long time)")
        mol_react = Molecule.from_xyz(os.path.join(work_dir, reactant_file))
        inp_path = os.path.join(work_dir, "orca_neb.inp")
        
        with open(inp_path, "w") as f:
            f.write("! B97-3c OptTS NEB-TS\n")
            f.write(f"%neb NEB_End_XYZFile \"{product_file}\" end\n")
            f.write(f"* xyz {charge} {multiplicity}\n")
            for a in mol_react.atoms:
                f.write(f"{a.symbol} {a.x} {a.y} {a.z}\n")
            f.write("*\n")
            
        cmd = f"orca orca_neb.inp > orca_neb.out"
        subprocess.run(cmd, shell=True, cwd=work_dir)
        
        out_path = os.path.join(work_dir, "orca_neb.out")
        if not os.path.exists(out_path):
             raise RuntimeError("orca_neb.out not found. The ORCA NEB failed.")
             
        # ORCA NEB-TS automatically performs the TS optimization and outputs it
        # We can extract the final geometry as ts_guess.xyz
        # For simplicity, if ORCA finishes properly, we just look for orca_neb.xyz
        final_xyz = os.path.join(work_dir, "orca_neb.xyz")
        if os.path.exists(final_xyz):
             shutil.copy2(final_xyz, os.path.join(work_dir, "ts_guess.xyz"))
             self._pes_log("Found final ORCA NEB-TS geometry. Saved as ts_guess.xyz")
        else:
             self._pes_log("Warning: orca_neb.xyz not found. Cannot set ts_guess.xyz.")

    def _generate_1d_pes_plot(self, work_dir, energies):
        try:
            import matplotlib.pyplot as plt
            fig = plt.figure(figsize=(8,6))
            plt.plot(range(1, len(energies)+1), energies, marker='o', linestyle='-', color='b')
            plt.title("Reactant -> Product Interpolation Path")
            plt.xlabel("Interpolation Step")
            plt.ylabel("Energy (Eh)")
            plt.grid(True)
            plot_path = os.path.join(work_dir, "Path_1D_Plot.png")
            plt.savefig(plot_path, dpi=300)
            plt.close()
            self._pes_log(f"1D Path Plot saved to {plot_path}")
        except Exception as e:
            self._pes_log(f"Failed to generate 1D plot: {str(e)}")

    def _generate_3d_pes_plot(self, work_dir, c1_steps, c2_steps, energies):
        try:
            import matplotlib.pyplot as plt
            from mpl_toolkits.mplot3d import Axes3D
            
            N = len(energies)
            # Adjust dimensions depending on inclusive/exclusive steps. Usually step count = N intervals -> N+1 points.
            # But let's just attempt a safe reshape.
            # Try to infer dims if c1_steps and c2_steps don't match exactly
            import math
            possible_cols = [int(c1_steps), int(c1_steps)+1, int(c2_steps), int(c2_steps)+1]
            cols = None
            for c in possible_cols:
                if c > 0 and N % c == 0:
                    cols = c
                    break
                    
            if cols is None:
                self._pes_log(f"Warning: Cannot reshape energy array ({N} points) for a 3D plot smoothly. Plot skipped.")
                return
                
            rows = N // cols
            Z = np.array(energies).reshape(rows, cols)
            X, Y = np.meshgrid(range(cols), range(rows))
            
            fig = plt.figure(figsize=(8,6))
            ax = fig.add_subplot(111, projection='3d')
            surf = ax.plot_surface(X, Y, Z, cmap='viridis')
            fig.colorbar(surf, ax=ax, shrink=0.5, aspect=5)
            ax.set_title("2D Relaxed Scan PES")
            ax.set_xlabel("Coord 1 (Steps)")
            ax.set_ylabel("Coord 2 (Steps)")
            ax.set_zlabel("Energy (Eh)")
            
            plot_path = os.path.join(work_dir, "PES_3D_Plot.png")
            plt.savefig(plot_path, dpi=300)
            plt.close()
            self._pes_log(f"3D PES Plot saved to {plot_path}")
        except Exception as e:
            self._pes_log(f"Failed to generate 3D plot: {str(e)}")

    def _load_local(self):
        fp = filedialog.askopenfilename(title="Select Engine Archive", filetypes=[("Archives", "*.tar.xz *.tar.gz *.tgz *.zip")])
        if fp:
            fname = os.path.basename(fp).lower()
            engine_type = "xTB" if "xtb" in fname else "CREST" if "crest" in fname else "ORCA" if "orca" in fname else "Unknown"
            self._append_text(f"\n[Local Load] Selected {engine_type} archive: {fname}\n")
            threading.Thread(target=self._extract_local_worker, args=(fp, engine_type), daemon=True).start()

    def _extract_local_worker(self, fp, engine_type, from_startup=False):
        try:
            os.makedirs(ENGINE_DIR, exist_ok=True)
            self._append_text(f"  Extracting {fp} into {ENGINE_DIR}...\n")
            if fp.lower().endswith(".zip"):
                import zipfile

                with zipfile.ZipFile(fp, "r") as zip_ref:
                    zip_ref.extractall(ENGINE_DIR)
            else:
                success_wsl = False
                if sys.platform == "win32":
                    wsl_fp = get_wsl_path(os.path.abspath(fp))
                    wsl_dir = get_wsl_path(os.path.abspath(ENGINE_DIR))
                    self._append_text("  [WSL] Using native Linux extraction to safely handle symlinks...\n")
                    cmd = f"wsl -e bash -c \"cd '{wsl_dir}' && tar -xf '{wsl_fp}'\""
                    rc = subprocess.run(cmd, shell=True, capture_output=True, text=True)
                    if rc.returncode == 0:
                        success_wsl = True
                    else:
                        self._append_text(f"  [WSL Tar Warning] Extraction issue: {rc.stderr.strip()}\n  Falling back to pure Python extractor...\n")
                if not success_wsl:
                    with tarfile.open(fp) as f:
                        for member in f.getmembers():
                            try:
                                f.extract(member, ENGINE_DIR)
                            except Exception:
                                pass
            inject_embedded_engines()
            self._append_text(f"[Success] {engine_type} locally installed and linked!\n")
            self.root.after(0, self._update_installation_labels)
            if from_startup:
                is_avail = is_tool_available(engine_type.lower())
                self.root.after(0, lambda e=engine_type, a=is_avail: self._update_startup_ui(e, a))
            else:
                self.root.after(0, lambda: messagebox.showinfo("Success", f"{engine_type} successfully installed from local file!"))
        except Exception as e:
            self._append_text(f"[Error] Failed to extract local archive: {str(e)}\n")
            if from_startup:
                self.root.after(0, lambda e=engine_type: self._update_startup_ui(e, False))


if __name__ == "__main__":
    IAKApp().mainloop()
