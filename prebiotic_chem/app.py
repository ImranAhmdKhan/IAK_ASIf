"""
prebiotic_chem.app
==================
Command-line interface for the Prebiotic Chemistry Computational App.

Usage examples::

    # List available prebiotic scenarios
    python -m prebiotic_chem --list-scenarios

    # List available molecules
    python -m prebiotic_chem --list-molecules

    # List reaction pathways
    python -m prebiotic_chem --list-pathways

    # Show detail for a specific scenario
    python -m prebiotic_chem --scenario warm_little_pond --info

    # Generate a 3-D structure from SMILES (requires RDKit)
    python -m prebiotic_chem --smiles "NCC(=O)O" --out-xyz glycine.xyz

    # Run the full pipeline on an XYZ file
    python -m prebiotic_chem \\
        --xyz molecule.xyz \\
        --scenario warm_little_pond \\
        --run-dir my_calculation \\
        --cores 8

    # Run pipeline but skip ORCA (xTB + CREST only)
    python -m prebiotic_chem \\
        --xyz molecule.xyz \\
        --scenario miller_urey \\
        --run-dir my_run \\
        --skip-orca-opt --skip-orca-sp

    # Boltzmann analysis from a CSV / list of energies
    python -m prebiotic_chem --boltzmann-csv energies.csv --temperature 333

    # Temperature sweep for a reaction (ΔH and ΔS required)
    python -m prebiotic_chem \\
        --temp-sweep \\
        --delta-H -5.2 \\
        --delta-S 12.0 \\
        --scenario hydrothermal_vent

    # Generate a full report
    python -m prebiotic_chem \\
        --scenario warm_little_pond \\
        --molecule glycine \\
        --pathway strecker_amino_acid \\
        --report report.txt
"""

from __future__ import annotations

import argparse
import csv
import logging
import sys
from pathlib import Path
from typing import List, Optional

# ---------------------------------------------------------------------------
# Module-level logger (configured after argument parsing)
# ---------------------------------------------------------------------------
logger = logging.getLogger("prebiotic_chem")


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="prebiotic_chem",
        description=(
            "Prebiotic Chemistry Computational App — "
            "professional-grade computational chemistry focused on early-Earth "
            "prebiotic conditions.\n\n"
            "Inspired by ChemRefine (Sterling Group, University of Texas at Dallas)."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    # ── Query / listing ──────────────────────────────────────────────────
    query = p.add_argument_group("Query and information")
    query.add_argument(
        "--list-scenarios", action="store_true",
        help="List all predefined prebiotic environment scenarios and exit.",
    )
    query.add_argument(
        "--list-molecules", action="store_true",
        help="List the built-in prebiotic molecule library and exit.",
    )
    query.add_argument(
        "--list-pathways", action="store_true",
        help="List predefined prebiotic reaction pathways and exit.",
    )
    query.add_argument(
        "--info", action="store_true",
        help="Show detailed information for the selected scenario / molecule / pathway.",
    )

    # ── Context selection ────────────────────────────────────────────────
    ctx = p.add_argument_group("Context selection")
    ctx.add_argument(
        "--scenario", metavar="NAME",
        default="warm_little_pond",
        help=(
            "Prebiotic environment scenario to use. "
            "Run --list-scenarios to see all options. "
            "Default: warm_little_pond."
        ),
    )
    ctx.add_argument(
        "--molecule", metavar="NAME",
        help="Molecule name from the built-in library (for info / report).",
    )
    ctx.add_argument(
        "--pathway", metavar="NAME",
        help="Reaction pathway name (for info / report).",
    )

    # ── Structure input ──────────────────────────────────────────────────
    struct = p.add_argument_group("Structure input")
    struct.add_argument(
        "--xyz", metavar="FILE",
        help="Input XYZ geometry file.",
    )
    struct.add_argument(
        "--smiles", metavar="SMILES",
        help="Generate 3-D structure from SMILES string (requires RDKit).",
    )
    struct.add_argument(
        "--out-xyz", metavar="FILE",
        help="Output XYZ file for --smiles conversion.",
    )

    # ── Pipeline control ─────────────────────────────────────────────────
    pipe = p.add_argument_group("Pipeline control")
    pipe.add_argument(
        "--run-dir", metavar="DIR", default="prebiotic_run",
        help="Root directory for pipeline output. Default: prebiotic_run/",
    )
    pipe.add_argument(
        "--charge", type=int, default=0,
        help="Total molecular charge. Default: 0.",
    )
    pipe.add_argument(
        "--mult", type=int, default=1,
        help="Spin multiplicity. Default: 1.",
    )
    pipe.add_argument(
        "--cores", type=int, default=4,
        help="CPU cores for xTB / CREST / ORCA. Default: 4.",
    )
    pipe.add_argument(
        "--maxcore", type=int, default=2000,
        help="Memory per ORCA core in MB. Default: 2000.",
    )
    pipe.add_argument(
        "--orca-opt-method", metavar="KEYWORD",
        default="B97-3c Opt Freq",
        help="ORCA keyword line for geometry optimisation. Default: 'B97-3c Opt Freq'.",
    )
    pipe.add_argument(
        "--orca-sp-method", metavar="KEYWORD",
        default="B3LYP D3BJ def2-TZVP TightSCF",
        help="ORCA keyword line for single-point. Default: 'B3LYP D3BJ def2-TZVP TightSCF'.",
    )
    pipe.add_argument(
        "--energy-window", type=float, default=5.0, metavar="KCAL",
        help="Energy window (kcal/mol) for conformer filtering. Default: 5.0.",
    )
    pipe.add_argument(
        "--rmsd-cutoff", type=float, default=0.5, metavar="ANG",
        help="RMSD cutoff (Å) for conformer clustering. Default: 0.5.",
    )
    pipe.add_argument(
        "--max-orca", type=int, default=10,
        help="Maximum structures to submit to ORCA. Default: 10.",
    )
    pipe.add_argument(
        "--boltzmann-pct", type=float, default=99.0,
        help="Cumulative Boltzmann population cutoff (%%). Default: 99.0.",
    )
    pipe.add_argument(
        "--skip-xtb",      action="store_true", help="Skip xTB pre-optimisation.",
    )
    pipe.add_argument(
        "--skip-crest",    action="store_true", help="Skip CREST conformer sampling.",
    )
    pipe.add_argument(
        "--skip-orca-opt", action="store_true", help="Skip ORCA geometry optimisation.",
    )
    pipe.add_argument(
        "--skip-orca-sp",  action="store_true", help="Skip ORCA single-point calculation.",
    )

    # ── Thermodynamic analysis ───────────────────────────────────────────
    thermo = p.add_argument_group("Thermodynamic analysis")
    thermo.add_argument(
        "--boltzmann-csv", metavar="FILE",
        help=(
            "CSV file with columns 'id' and 'energy_hartree'. "
            "Performs Boltzmann population analysis and prints results."
        ),
    )
    thermo.add_argument(
        "--temperature", type=float, metavar="K",
        help="Temperature in Kelvin for Boltzmann analysis. Overrides scenario value.",
    )
    thermo.add_argument(
        "--temp-sweep", action="store_true",
        help="Run ΔG / K_eq temperature sweep (requires --delta-H and --delta-S).",
    )
    thermo.add_argument(
        "--delta-H", type=float, metavar="KCAL_MOL",
        help="Reaction enthalpy (kcal/mol) for --temp-sweep.",
    )
    thermo.add_argument(
        "--delta-S", type=float, metavar="CAL_MOL_K",
        help="Reaction entropy (cal/mol/K) for --temp-sweep.",
    )

    # ── Report ───────────────────────────────────────────────────────────
    rep = p.add_argument_group("Report generation")
    rep.add_argument(
        "--report", metavar="FILE",
        help="Write a comprehensive plain-text report to FILE.",
    )

    # ── Misc ─────────────────────────────────────────────────────────────
    p.add_argument(
        "--verbose", "-v", action="store_true",
        help="Enable verbose (DEBUG) logging.",
    )

    return p


# ---------------------------------------------------------------------------
# Action handlers
# ---------------------------------------------------------------------------

def _action_list_scenarios() -> None:
    from .conditions import list_scenarios
    print("\nAvailable prebiotic scenarios:")
    print("-" * 50)
    for s in list_scenarios():
        t_lo, t_hi = s.temp_range_K
        print(
            f"  {s.name:<25}  {s.display_name}  "
            f"[T: {t_lo - 273.15:.0f}–{t_hi - 273.15:.0f} °C]"
        )
    print()


def _action_list_molecules(category: Optional[str] = None) -> None:
    from .molecules import list_molecules, list_categories
    cats = list_categories()
    print(f"\nPrebiotic molecule library ({len(list_molecules())} molecules):")
    print("-" * 60)
    for cat in cats:
        mols = list_molecules(cat)
        print(f"\n  [{cat.upper()}]")
        for m in mols:
            print(f"    {m.name:<35}  {m.formula}")
    print()


def _action_list_pathways() -> None:
    from .reactions import PATHWAYS
    print("\nPrebiotic reaction pathways:")
    print("-" * 50)
    for name, pw in PATHWAYS.items():
        print(f"  {name:<30}  {pw.display_name}")
    print()


def _action_info(args: argparse.Namespace) -> None:
    from .conditions import get_scenario
    from .molecules import get_molecule
    from .reactions import get_pathway

    if args.scenario:
        sc = get_scenario(args.scenario)
        if sc:
            print()
            print(sc.summary())
        else:
            print(f"Scenario '{args.scenario}' not found.")

    if args.molecule:
        m = get_molecule(args.molecule)
        if m:
            print(f"\n{m.name}  ({m.formula})  — {m.category}")
            print(f"  SMILES      : {m.smiles}")
            if m.aliases:
                print(f"  Aliases     : {', '.join(m.aliases)}")
            if m.relevance:
                print(f"  Relevance   : {m.relevance}")
        else:
            print(f"Molecule '{args.molecule}' not found.")

    if args.pathway:
        pw = get_pathway(args.pathway)
        if pw:
            print()
            print(pw.summary())
        else:
            print(f"Pathway '{args.pathway}' not found.")


def _action_smiles_to_xyz(args: argparse.Namespace) -> None:
    from .io_utils import smiles_to_xyz
    out = args.out_xyz or "output.xyz"
    syms, coords = smiles_to_xyz(args.smiles, output_path=out)
    print(f"Generated 3-D structure: {out}  ({len(syms)} atoms)")


def _action_boltzmann_csv(args: argparse.Namespace) -> None:
    from .conditions import get_scenario
    from .thermodynamics import boltzmann_populations

    sc = get_scenario(args.scenario)
    T = args.temperature or (sc.temperature_K if sc else 298.15)

    rows = list(csv.DictReader(open(args.boltzmann_csv)))
    ids      = [r["id"] for r in rows]
    energies = [float(r["energy_hartree"]) for r in rows]

    pops = boltzmann_populations(energies, ids, temperature_K=T)
    print(f"\nBoltzmann populations at T = {T:.1f} K")
    print("-" * 60)
    print(f"  {'ID':<20}  {'ΔE (kcal/mol)':>14}  {'pop (%)':>8}")
    print("  " + "-" * 48)
    for row in pops:
        print(f"  {row['id']:<20}  {row['rel_energy_kcal']:>14.3f}  {row['population_pct']:>8.2f}")
    print()


def _action_temp_sweep(args: argparse.Namespace) -> None:
    from .report import temperature_sweep_report
    print(f"\nTemperature sweep  ΔH = {args.delta_H:.2f} kcal/mol,  "
          f"ΔS = {args.delta_S:.2f} cal/mol/K")
    print("-" * 60)
    print(temperature_sweep_report(args.delta_H, args.delta_S))
    print()


def _action_pipeline(args: argparse.Namespace) -> None:
    from .conditions import get_scenario
    from .pipeline import PrebioticPipeline

    sc = get_scenario(args.scenario)
    if sc is None:
        print(f"ERROR: scenario '{args.scenario}' not found.  "
              "Run --list-scenarios to see options.", file=sys.stderr)
        sys.exit(1)

    pl = PrebioticPipeline(
        run_dir=args.run_dir,
        scenario=sc,
        input_xyz=args.xyz,
        charge=args.charge,
        multiplicity=args.mult,
        n_cores=args.cores,
        maxcore_mb=args.maxcore,
        orca_opt_method=args.orca_opt_method,
        orca_sp_method=args.orca_sp_method,
        energy_window_kcal=args.energy_window,
        rmsd_cutoff=args.rmsd_cutoff,
        max_orca_structures=args.max_orca,
        boltzmann_cutoff_pct=args.boltzmann_pct,
        skip_xtb=args.skip_xtb,
        skip_crest=args.skip_crest,
        skip_orca_opt=args.skip_orca_opt,
        skip_orca_sp=args.skip_orca_sp,
    )
    pl.run()


def _action_report(args: argparse.Namespace) -> None:
    from .conditions import get_scenario
    from .molecules import get_molecule
    from .reactions import get_pathway
    from .report import PrebioticReport

    sc  = get_scenario(args.scenario)
    mol = get_molecule(args.molecule) if args.molecule else None
    pw  = get_pathway(args.pathway)  if args.pathway  else None
    rd  = Path(args.run_dir) if args.run_dir else None

    report = PrebioticReport(scenario=sc, molecule=mol, pathway=pw, run_dir=rd)
    out    = report.write(args.report)
    print(f"Report written to: {out}")


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main(argv: Optional[List[str]] = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)

    # Logging
    level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        format="%(levelname)-8s %(message)s",
        level=level,
        stream=sys.stdout,
    )

    # ── Query / listing ──────────────────────────────────────────────────
    if args.list_scenarios:
        _action_list_scenarios()
        return 0
    if args.list_molecules:
        _action_list_molecules()
        return 0
    if args.list_pathways:
        _action_list_pathways()
        return 0
    if args.info:
        _action_info(args)
        return 0

    # ── SMILES → XYZ ────────────────────────────────────────────────────
    if args.smiles:
        _action_smiles_to_xyz(args)
        return 0

    # ── Boltzmann CSV ────────────────────────────────────────────────────
    if args.boltzmann_csv:
        _action_boltzmann_csv(args)
        return 0

    # ── Temperature sweep ────────────────────────────────────────────────
    if args.temp_sweep:
        if args.delta_H is None or args.delta_S is None:
            print("ERROR: --temp-sweep requires --delta-H and --delta-S.", file=sys.stderr)
            return 1
        _action_temp_sweep(args)
        return 0

    # ── Report (without pipeline) ────────────────────────────────────────
    if args.report and not args.xyz:
        _action_report(args)
        return 0

    # ── Pipeline ─────────────────────────────────────────────────────────
    if args.xyz:
        _action_pipeline(args)
        if args.report:
            _action_report(args)
        return 0

    # ── Default: show help ───────────────────────────────────────────────
    parser.print_help()
    return 0
