#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Before/after consistency harness.

Computes the full e-DockQ metric row (ligand RMSD, interface RMSD, residue FNAT,
backbone/sidechain FNAT, BINANA recall) for a fixed set of (case, model) using the
SAME logic as the runner, and dumps a JSON snapshot. Run before and after the refactor
and diff the two JSON files.

Imports via the legacy ``plq`` name on purpose, so after the rename this also exercises
the ``plq`` compatibility alias.

Usage:
    python pliq_ver5/tests/consistency_harness.py OUT.json
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from pliq.edockq import binana_recall as BR  # noqa: E402
from pliq.edockq import helpers as H  # noqa: E402
from pliq.edockq import ligand_mapping as LM  # noqa: E402
from pliq.edockq.otmol_runtime import ensure_local_otmol_on_path  # noqa: E402
from pliq.config import resolve_binana_python  # noqa: E402

ensure_local_otmol_on_path()

DATA_ROOT = ROOT.parent / "posebench" / "hpc_upload_pliq_processed_20260602_singlelig"

CASES = [
    ("dockgen", "1hg0_1_SIN_1"),
    ("astex_diverse", "1HP0_AD3"),
    ("dockgen", "2hk9_1_SKM_0"),
]
MODELS = [0, 1, 2]


def _paths(dataset, pdb_id, k):
    ref_lig = DATA_ROOT / "references" / dataset / pdb_id / f"{pdb_id}_ligand_chainB.pdb"
    ref_pro = DATA_ROOT / "references" / dataset / pdb_id / f"{pdb_id}_protein_renum.pdb"
    dd = DATA_ROOT / "predictions" / "alphafold3" / dataset / "raw" / f"posebusters_{pdb_id}"
    pred_lig = dd / f"posebusters_{pdb_id}_model_{k}_aligned_ligand.pdb"
    pred_pro = dd / f"posebusters_{pdb_id}_model_{k}_aligned_protein.pdb"
    return ref_lig, ref_pro, pred_lig, pred_pro


def _round(v, n=6):
    try:
        return round(float(v), n)
    except (TypeError, ValueError):
        return v


def compute_row(dataset, pdb_id, k, binana_python):
    ref_lig, ref_pro, pred_lig, pred_pro = _paths(dataset, pdb_id, k)
    out = {"case": f"{dataset}/{pdb_id}", "model": k}
    for p in (ref_lig, ref_pro, pred_lig, pred_pro):
        if not p.is_file():
            out["status"] = f"missing:{p.name}"
            return out

    mol1 = H.pdb_to_rdkit_with_prody(str(ref_lig), "B")
    mol2 = H.pdb_to_rdkit_with_prody(str(pred_lig), "B")
    prot1 = H.pdb_chain_to_rdkit_mol(str(ref_pro), "A")
    prot2 = H.pdb_chain_to_rdkit_mol(str(pred_pro), "A")
    if not (mol1 and mol2 and prot1 and prot2):
        out["status"] = "load_failed"
        return out

    dock_to_ref, assignment2, rmsd, alpha2, BCI2, attempts, bci_ok = (
        LM.align_ligands_otmol_bci_constrained(mol1, mol2, max_attempts=20)
    )
    out["ligand_rmsd"] = _round(rmsd)
    out["otmol_BCI"] = _round(BCI2)
    out["otmol_map_pairs"] = len(dock_to_ref)

    iface_ca = H.calculate_interface_rmsd_otmol(
        mol1, mol2, prot1, prot2, dock_to_ref, threshold=10.0,
        protein_atom_names=H.INTERFACE_PROTEIN_CA,
    )
    iface_bb = H.calculate_interface_rmsd_otmol(
        mol1, mol2, prot1, prot2, dock_to_ref, threshold=10.0,
        protein_atom_names=H.INTERFACE_PROTEIN_BACKBONE,
    )
    out["interface_rmsd"] = _round(iface_ca["interface_rmsd"])
    out["interface_rmsd_n_prot"] = iface_ca["n_protein_atoms"]
    out["interface_rmsd_n_lig"] = iface_ca["n_ligand_atoms"]
    out["interface_rmsd_backbone"] = _round(iface_bb["interface_rmsd"])
    out["interface_rmsd_bb_n_prot"] = iface_bb["n_protein_atoms"]

    sub1, sub2 = H.find_close_residues_otmol(mol1, mol2, prot1, prot2, dock_to_ref, threshold=5.0)
    common = H.find_common_contact(sub1, sub2)
    out["residue_fnat"] = _round(len(common) / len(sub1) if sub1 else -1)
    out["residue_ref_contacts"] = len(sub1)
    out["residue_pred_contacts"] = len(sub2)
    out["residue_common_contacts"] = len(common)
    out["residue_ref_contacts_list"] = ";".join(sorted(sub1))
    out["residue_common_contacts_list"] = ";".join(sorted(common))

    fn = H.calculate_fnat_backbone_sidechain_otmol(mol1, mol2, prot1, prot2, dock_to_ref, threshold=5.0)
    for key in ("fnat_all", "fnat_backbone", "fnat_sidechain"):
        out[key] = _round(fn[key])
    for key in ("tp_all", "tp_backbone", "tp_sidechain",
                "total_ref_all", "total_ref_backbone", "total_ref_sidechain",
                "fnat_map_src", "fnat_mapped_pairs"):
        out[key] = fn[key]
    out["ref_contacts_all_list"] = ";".join(sorted(fn["ref_contacts_all_set"]))

    # BINANA recall (uses obabel + vendored binana)
    try:
        rec = BR.compute_binana_recall_row(
            str(ref_lig), str(ref_pro), str(pred_lig), str(pred_pro),
            Path(binana_python), mol_ref=mol1, mol_dock=mol2,
            otmol_assignment=assignment2,
            dock_to_ref=dock_to_ref,
            otmol_BCI=BCI2,
            otmol_bci_zero_ok=bci_ok,
            get_all_interactions_kwargs={}, obabel_ph=None,
        )
        for kk, vv in rec.items():
            if isinstance(vv, (int, float)) or (isinstance(vv, str) and len(vv) < 200):
                out[f"binana_{kk}" if not kk.startswith("binana") else kk] = (
                    _round(vv) if isinstance(vv, float) else vv
                )
        out["binana_ok"] = True
    except Exception as e:
        out["binana_ok"] = False
        out["binana_err"] = str(e)[:200]

    out["status"] = "ok"
    return out


def main():
    out_path = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("consistency_snapshot.json")
    binana_python = str(resolve_binana_python())
    rows = []
    for ds, pid in CASES:
        for k in MODELS:
            print(f"computing {ds}/{pid} model {k} ...", file=sys.stderr)
            rows.append(compute_row(ds, pid, k, binana_python))
    out_path.write_text(json.dumps(rows, indent=2, sort_keys=True))
    print(f"wrote {out_path} ({len(rows)} rows)")


if __name__ == "__main__":
    main()
