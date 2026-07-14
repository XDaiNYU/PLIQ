#!/usr/bin/env python3
"""Run pliq ver5 BINANA detail for CASP15 T1124 model 0 (rxatm) and print TP/FP/FN."""
from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REPO = ROOT.parent
sys.path.insert(0, str(ROOT))

from pliq.edockq import binana_recall as BR
from pliq.edockq import helpers as H
from pliq.edockq import ligand_mapping as LM
from pliq.edockq.otmol_runtime import ensure_local_otmol_on_path

ensure_local_otmol_on_path()

CASE = "T1124"
MODEL = 0
LEVEL = "rxatm"
RECEPTOR_MODE = LM.RECEPTOR_SIG_RID_ATOM
SFX = "_rxatm"

REF_LIG = REPO / "posebench/prepared_pliq_af3/references/casp15/T1124/T1124_ligand_chainB.pdb"
REF_REC = REPO / "posebench/prepared_pliq_af3/references/casp15/T1124/T1124_protein_renum.pdb"
DOCK_LIG = (
    REPO
    / "posebench/prepared_pliq_af3/predictions/alphafold3/casp15/raw/posebusters_T1124"
    / f"posebusters_T1124_model_{MODEL}_aligned_ligand.pdb"
)
DOCK_REC = (
    REPO
    / "posebench/prepared_pliq_af3/predictions/alphafold3/casp15/raw/posebusters_T1124"
    / f"posebusters_T1124_model_{MODEL}_aligned_protein.pdb"
)
BINANA_PY = REPO / "tools/binana/python"
OUT_DIR = REPO / "pliq_ver5/output"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def main() -> None:
    mol_ref = H.pdb_to_rdkit_with_prody(str(REF_LIG), "B")
    mol_dock = H.pdb_to_rdkit_with_prody(str(DOCK_LIG), "B")
    if mol_ref is None or mol_dock is None:
        raise RuntimeError("Failed to load ligand mols via ProDy (chain B)")

    dock_to_ref, assignment, rmsd, alpha, bci, n_try, bci_ok = (
        LM.align_ligands_otmol_bci_constrained(mol_ref, mol_dock, max_attempts=20)
    )
    diag = LM.otmol_mapping_diagnostics(mol_ref, mol_dock, assignment)
    print(
        f"OTMol: mapped={len(dock_to_ref)} rmsd={rmsd:.4f} alpha={alpha} BCI={bci} "
        f"bci_zero_ok={bci_ok} attempts={n_try} atomicnum_mismatch={diag['atomicnum_mismatch']}"
    )

    row = BR.compute_binana_recall_row(
        str(REF_LIG),
        str(REF_REC),
        str(DOCK_LIG),
        str(DOCK_REC),
        BINANA_PY,
        mol_ref=mol_ref,
        mol_dock=mol_dock,
        otmol_assignment=assignment,
        dock_to_ref=dock_to_ref,
        otmol_BCI=bci,
        otmol_alpha=alpha,
        otmol_rmsd=rmsd,
        otmol_align_attempts=n_try,
        otmol_bci_zero_ok=bci_ok,
    )
    print(
        f"BINANA→RDKit coord map: {row['binana_ligand_binana_to_rdkit_n']} atoms | "
        f"CSV otmol_BCI={row['binana_otmol_BCI']} mismatch={row['binana_otmol_atomicnum_mismatch']}"
    )

    # One-row CSV with all ver5 columns
    csv_path = OUT_DIR / f"{CASE}_model{MODEL}_binana_ver5_row.csv"
    with open(csv_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(row.keys()))
        w.writeheader()
        w.writerow(row)
    print(f"Wrote row CSV: {csv_path}")

    # Human-readable report for six types + hydrophobic @ rxatm
    report_lines = [
        f"# BINANA detail report: {CASE} model {MODEL} @ {LEVEL}",
        f"OTMol map pairs: {len(dock_to_ref)} | BCI={bci} | bci_zero_ok={bci_ok} | "
        f"atomicnum_mismatch={diag['atomicnum_mismatch']}",
        f"BINANA PDB-index → RDKit map: {row['binana_ligand_binana_to_rdkit_n']} atoms",
        "",
    ]

    types = BR.BINANA_SIX_TYPE_KEYS + [("hydrophobicContacts", "hydrophobic")]

    detail_json_path = OUT_DIR / f"{CASE}_model{MODEL}_binana_details_rxatm.json"
    all_details: dict = {}

    for json_key, short in types:
        tp = json.loads(row[f"binana_{short}_tp_details{SFX}"])
        fn = json.loads(row[f"binana_{short}_fn_details{SFX}"])
        fp = json.loads(row[f"binana_{short}_fp_details{SFX}"])
        stats = row[f"binana_{short}_stats{SFX}"]
        all_details[short] = {"stats": stats, "tp": tp, "fn": fn, "fp": fp}
        report_lines.append(
            BR.format_interaction_summary_human(
                tp, fn, fp, interaction_label=f"{short} ({json_key}) stats={stats}"
            )
        )
        report_lines.append("")

    detail_json_path.write_text(
        json.dumps(all_details, ensure_ascii=True, indent=2), encoding="utf-8"
    )
    print(f"Wrote detail JSON: {detail_json_path}")

    report_path = OUT_DIR / f"{CASE}_model{MODEL}_binana_report_rxatm.txt"
    report_path.write_text("\n".join(report_lines), encoding="utf-8")
    print(f"Wrote report: {report_path}")
    print("\n" + "\n".join(report_lines[:80]))


if __name__ == "__main__":
    main()
