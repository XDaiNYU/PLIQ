#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Run OTMol FNAT validation on prepared PoseBench cases and print a short report.

Usage:
  python scripts/validate_fnat_otmol_cases.py
  python scripts/validate_fnat_otmol_cases.py --data-root /path/to/hpc_upload_..._singlelig
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from pliq.edockq import helpers as H
from pliq.edockq import ligand_mapping as LM

DEFAULT_CASES = [
    ("dockgen", "1hg0_1_SIN_1"),
    ("astex_diverse", "1HP0_AD3"),
    ("dockgen", "2hk9_1_SKM_0"),
]


def run_one(data_root: Path, dataset: str, pdb_id: str, model_k: int = 0) -> dict:
    ref_lig = data_root / "references" / dataset / pdb_id / f"{pdb_id}_ligand_chainB.pdb"
    ref_pro = data_root / "references" / dataset / pdb_id / f"{pdb_id}_protein_renum.pdb"
    dock_dir = data_root / "predictions" / "alphafold3" / dataset / "raw" / f"posebusters_{pdb_id}"
    pred_lig = dock_dir / f"posebusters_{pdb_id}_model_{model_k}_aligned_ligand.pdb"
    pred_pro = dock_dir / f"posebusters_{pdb_id}_model_{model_k}_aligned_protein.pdb"

    mol_ref = H.pdb_to_rdkit_with_prody(str(ref_lig), "B")
    mol_dock = H.pdb_to_rdkit_with_prody(str(pred_lig), "B")
    prot_ref = H.pdb_chain_to_rdkit_mol(str(ref_pro), "A")
    prot_dock = H.pdb_chain_to_rdkit_mol(str(pred_pro), "A")

    dock_to_ref, _a, rmsd, alpha, bci, n_try, bci0 = LM.align_ligands_otmol_bci_constrained(
        mol_ref, mol_dock, max_attempts=20
    )
    otm = H.calculate_fnat_backbone_sidechain_otmol(
        mol_ref, mol_dock, prot_ref, prot_dock, dock_to_ref, threshold=5.0
    )

    return {
        "pdb_id": pdb_id,
        "dataset": dataset,
        "ligand_rmsd": rmsd,
        "mapped_pairs": len(dock_to_ref),
        "ref_atoms": mol_ref.GetNumAtoms(),
        "dock_atoms": mol_dock.GetNumAtoms(),
        "fnat_all_otmol": otm["fnat_all"],
        "fnat_backbone_otmol": otm["fnat_backbone"],
        "fnat_sidechain_otmol": otm["fnat_sidechain"],
        "tp_otmol": otm["tp_all"],
        "ref_n_otmol": otm["total_ref_all"],
        "bci_zero": bci0,
        "otmol_attempts": n_try,
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--data-root",
        type=Path,
        default=ROOT.parent / "posebench" / "hpc_upload_pliq_processed_20260602_singlelig",
    )
    args = ap.parse_args()
    data_root = args.data_root.resolve()
    if not data_root.is_dir():
        print(f"ERROR: data root not found: {data_root}", file=sys.stderr)
        return 1

    LM.require_otmol_importable()
    print(f"PLIQ OTMol FNAT validation  data={data_root}\n")
    print(
        f"{'case':<22} {'map':>4} {'rmsd':>7} {'fnat_all':>8} {'fnat_bb':>8} "
        f"{'fnat_sc':>8} {'tp':>5} {'ref_n':>5}"
    )
    print("-" * 85)

    ok = 0
    for dataset, pdb_id in DEFAULT_CASES:
        try:
            r = run_one(data_root, dataset, pdb_id)
            ok += 1
            print(
                f"{pdb_id:<22} {r['mapped_pairs']:>4} {r['ligand_rmsd']:>7.3f} "
                f"{r['fnat_all_otmol']:>8.4f} {r['fnat_backbone_otmol']:>8.4f} "
                f"{r['fnat_sidechain_otmol']:>8.4f} {r['tp_otmol']:>5} {r['ref_n_otmol']:>5}"
            )
        except Exception as e:
            print(f"{pdb_id:<22} FAIL: {e}")

    print("-" * 85)
    print(f"OK: {ok}/{len(DEFAULT_CASES)}")
    return 0 if ok == len(DEFAULT_CASES) else 1


if __name__ == "__main__":
    raise SystemExit(main())
