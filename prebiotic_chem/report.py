"""
prebiotic_chem.report
=====================
Report generation for prebiotic chemistry computational results.

Produces:
- Plain-text ASCII report for terminal / logfile consumption
- CSV table of Boltzmann populations
- H-bond analysis table
- Reaction pathway overview
- Temperature-sweep thermodynamics table
"""

from __future__ import annotations

import csv
import io
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from .analysis import detect_hbonds, molecular_formula
from .conditions import PrebioticScenario
from .molecules import PrebioticMolecule
from .reactions import ReactionPathway
from .thermodynamics import boltzmann_populations, temperature_sweep

# Separator lines for ASCII reports
_SEP  = "=" * 70
_SEP2 = "-" * 70


def _centre(text: str, width: int = 70) -> str:
    return text.center(width)


# ---------------------------------------------------------------------------
# Main report class
# ---------------------------------------------------------------------------

class PrebioticReport:
    """
    Assembles a comprehensive plain-text report for a prebiotic calculation run.

    Parameters
    ----------
    scenario:
        The prebiotic environment scenario used.
    molecule:
        (Optional) molecule metadata record.
    pathway:
        (Optional) reaction pathway analysed.
    run_dir:
        The pipeline run directory (used to locate CSV files).
    """

    def __init__(
        self,
        scenario: PrebioticScenario,
        molecule: Optional[PrebioticMolecule] = None,
        pathway: Optional[ReactionPathway] = None,
        run_dir: Optional[str | Path] = None,
    ) -> None:
        self.scenario  = scenario
        self.molecule  = molecule
        self.pathway   = pathway
        self.run_dir   = Path(run_dir) if run_dir else None

    # ------------------------------------------------------------------ #
    # Public API                                                            #
    # ------------------------------------------------------------------ #

    def generate(self) -> str:
        """Return the full report as a single string."""
        sections = [
            self._header(),
            self._scenario_section(),
        ]
        if self.molecule:
            sections.append(self._molecule_section())
        if self.pathway:
            sections.append(self._pathway_section())
        if self.run_dir:
            sections.append(self._results_section())
        sections.append(self._footer())
        return "\n\n".join(sections)

    def write(self, path: str | Path) -> Path:
        """Write the report to a text file."""
        path = Path(path)
        path.write_text(self.generate())
        return path

    # ------------------------------------------------------------------ #
    # Sections                                                              #
    # ------------------------------------------------------------------ #

    def _header(self) -> str:
        lines = [
            _SEP,
            _centre("PREBIOTIC CHEMISTRY COMPUTATIONAL REPORT"),
            _centre("prebiotic_chem — Professional Computational Chemistry"),
            _SEP,
        ]
        return "\n".join(lines)

    def _scenario_section(self) -> str:
        s = self.scenario
        t_lo, t_hi = s.temp_range_K
        ph_lo, ph_hi = s.pH_range
        is_lo, is_hi = s.ionic_strength_range_M

        lines = [
            "PREBIOTIC ENVIRONMENT",
            _SEP2,
            f"  Scenario        : {s.display_name}",
            f"  Temperature     : {s.temperature_C:.1f} °C  "
            f"(range {t_lo - 273.15:.0f}–{t_hi - 273.15:.0f} °C)",
            f"  pH              : {s.pH:.1f}  (range {ph_lo:.1f}–{ph_hi:.1f})",
            f"  Ionic strength  : {s.ionic_strength_M:.3f} M  "
            f"(range {is_lo:.3f}–{is_hi:.3f} M)",
            f"  Solvent         : {s.solvent}",
            f"  Atmosphere      : {s.atmosphere or 'N/A'}",
            f"  xTB flags       : {' '.join(s.xtb_flags())}",
            f"  ORCA solvation  : {s.orca_solvent_keyword()}",
        ]
        if s.description:
            lines.append("")
            for chunk in _wrap(s.description, width=66, indent="  "):
                lines.append(chunk)
        if s.relevant_reactions:
            lines.append(f"  Reactions       : {', '.join(s.relevant_reactions)}")
        return "\n".join(lines)

    def _molecule_section(self) -> str:
        m = self.molecule
        lines = [
            "MOLECULE",
            _SEP2,
            f"  Name     : {m.name}",
            f"  Formula  : {m.formula}",
            f"  SMILES   : {m.smiles}",
            f"  Category : {m.category}",
        ]
        if m.aliases:
            lines.append(f"  Aliases  : {', '.join(m.aliases)}")
        if m.relevance:
            lines.append("")
            lines.append("  Prebiotic relevance:")
            for chunk in _wrap(m.relevance, width=64, indent="    "):
                lines.append(chunk)
        return "\n".join(lines)

    def _pathway_section(self) -> str:
        pw = self.pathway
        lines = [
            "REACTION PATHWAY",
            _SEP2,
            f"  Pathway  : {pw.display_name}",
        ]
        if pw.description:
            lines.append("")
            for chunk in _wrap(pw.description, width=66, indent="  "):
                lines.append(chunk)
        lines.append("")
        for i, step in enumerate(pw.steps, 1):
            lhs = " + ".join(step.reactants)
            rhs = " + ".join(step.products)
            dG_txt = (f"  ΔG ≈ {step.delta_G_kcal:.1f} kcal/mol"
                      if step.delta_G_kcal is not None else "")
            lines.append(f"  Step {i}: {step.display_name}")
            lines.append(f"    {lhs} → {rhs}{dG_txt}")
            if step.conditions:
                lines.append(f"    Conditions: {step.conditions}")
            if step.references:
                lines.append(f"    References: {'; '.join(step.references)}")
        return "\n".join(lines)

    def _results_section(self) -> str:
        lines = ["COMPUTATIONAL RESULTS", _SEP2]
        summary_csv = self.run_dir / "summary.csv"
        if summary_csv.exists():
            rows = list(csv.DictReader(open(summary_csv)))
            if rows:
                lines.append("")
                lines.append("  Conformer Boltzmann Populations")
                lines.append("  " + "-" * 60)
                header = f"  {'Rank':>4}  {'ID':<20}  {'ΔE (kcal/mol)':>14}  {'Pop (%)':>8}"
                lines.append(header)
                lines.append("  " + "-" * 60)
                for row in rows:
                    lines.append(
                        f"  {row['rank']:>4}  {row['id']:<20}  "
                        f"{row['rel_energy_kcal']:>14}  {row['population_pct']:>8}"
                    )
        else:
            lines.append("  (No results CSV found — run pipeline first.)")
        return "\n".join(lines)

    def _footer(self) -> str:
        return "\n".join([
            _SEP,
            _centre("Generated by prebiotic_chem"),
            _centre("Inspired by ChemRefine (Sterling Group, UT Dallas)"),
            _SEP,
        ])


# ---------------------------------------------------------------------------
# Standalone report helpers
# ---------------------------------------------------------------------------

def hbond_report(
    symbols: Sequence[str],
    coords: np.ndarray,
) -> str:
    """Return a plain-text H-bond analysis table for a geometry."""
    hbonds = detect_hbonds(symbols, coords)
    if not hbonds:
        return "No hydrogen bonds detected."
    lines = [
        f"{'H-idx':>6}  {'D-idx':>6}  {'A-idx':>6}  "
        f"{'H…A (Å)':>9}  {'D…A (Å)':>9}  {'∠D-H…A (°)':>11}",
        "-" * 58,
    ]
    for hb in hbonds:
        lines.append(
            f"{hb['donor_H']:>6}  {hb['donor_heavy']:>6}  {hb['acceptor']:>6}  "
            f"{hb['dist_HA_ang']:>9.3f}  {hb['dist_DA_ang']:>9.3f}  "
            f"{hb['angle_DHA_deg']:>11.1f}"
        )
    return "\n".join(lines)


def temperature_sweep_report(
    delta_H_kcal: float,
    delta_S_cal_mol_K: float,
    temperatures_K: Optional[Sequence[float]] = None,
) -> str:
    """Return a temperature-sweep thermodynamics table as a string."""
    rows = temperature_sweep(delta_H_kcal, delta_S_cal_mol_K, temperatures_K)
    lines = [
        f"{'T (K)':>7}  {'T (°C)':>7}  {'ΔG (kcal/mol)':>15}  "
        f"{'K_eq':>12}  {'Spontaneous':>12}",
        "-" * 60,
    ]
    for row in rows:
        lines.append(
            f"{row['temperature_K']:>7.1f}  {row['temperature_C']:>7.1f}  "
            f"{row['delta_G_kcal']:>15.3f}  {row['K_eq']:>12.3e}  "
            f"{'yes' if row['spontaneous'] else 'no':>12}"
        )
    return "\n".join(lines)


def write_population_csv(
    energies_hartree: Sequence[float],
    ids: Optional[Sequence[str]],
    temperature_K: float,
    path: str | Path,
) -> Path:
    """Write a Boltzmann population CSV and return the path."""
    pops = boltzmann_populations(energies_hartree, ids, temperature_K)
    path = Path(path)
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(
            f, fieldnames=["id", "energy_hartree", "rel_energy_kcal",
                           "rel_energy_kJ", "boltzmann_weight", "population_pct"]
        )
        writer.writeheader()
        for row in pops:
            writer.writerow({k: (f"{v:.6f}" if isinstance(v, float) else v)
                             for k, v in row.items()})
    return path


# ---------------------------------------------------------------------------
# Text utilities
# ---------------------------------------------------------------------------

def _wrap(text: str, width: int = 70, indent: str = "") -> List[str]:
    """Naive word-wrapper."""
    words = text.split()
    lines: List[str] = []
    current = indent
    for word in words:
        if current == indent:
            current += word
        elif len(current) + 1 + len(word) <= width:
            current += " " + word
        else:
            lines.append(current)
            current = indent + word
    if current != indent:
        lines.append(current)
    return lines
