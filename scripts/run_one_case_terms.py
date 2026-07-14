#!/usr/bin/env python3
"""Run one case: PLIQ raw metrics + local eDockQ scoring (eDockQ1 + eDockQ6_1..4)."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pliq.edockq import binana_recall as BR
from pliq.edockq import helpers as H
from pliq.edockq import ligand_mapping as LM
from pliq.edockq.otmol_runtime import ensure_local_otmol_on_path
from pliq_edockq_scoring import PB_SOURCE, score_edockq_stepwise

ensure_local_otmol_on_path()

REPO = ROOT.parent
BINANA = ROOT / "pliq" / "_vendor" / "binana_python"
LEVELS = ["default", "rxatm", "rxidx", "rxres"]


def _leads_pep(case: str, model: int) -> dict:
    base = REPO / "LEADS_PEP" / "pliq_input_models0_49_unfixed_ref"
    pose = base / f"posebusters_{case}"
    return {
        "case": case,
        "model": model,
        "ref_lig": base / "reference_inputs" / case / "ligand.pdb",
        "ref_pro": base / "reference_inputs" / case / "protein.pdb",
        "dock_lig": pose / f"posebusters_{case}_model_{model}_aligned_ligand.pdb",
        "dock_pro": pose / f"posebusters_{case}_model_{model}_aligned_protein.pdb",
    }


def _casp15(case: str, model: int) -> dict:
    refd = REPO / "posebench/prepared_pliq_af3/references/casp15" / case
    pose = (
        REPO
        / "posebench/prepared_pliq_af3/predictions/alphafold3/casp15/raw"
        / f"posebusters_{case}"
    )
    return {
        "case": case,
        "model": model,
        "ref_lig": refd / f"{case}_ligand_chainB.pdb",
        "ref_pro": refd / f"{case}_protein_renum.pdb",
        "dock_lig": pose / f"posebusters_{case}_model_{model}_aligned_ligand.pdb",
        "dock_pro": pose / f"posebusters_{case}_model_{model}_aligned_protein.pdb",
    }


CONFIGS = {
    "1B9J": lambda: _leads_pep("1B9J", 0),
    "T1124": lambda: _casp15("T1124", 0),
}


def _pick(row: pd.Series, keys: list[str]) -> dict:
    return {k: row[k] for k in keys if k in row.index}


def main() -> None:
    key = sys.argv[1] if len(sys.argv) > 1 else "1B9J"
    if key not in CONFIGS:
        raise SystemExit(f"Unknown case '{key}'. Available: {', '.join(CONFIGS)}")
    cfg = CONFIGS[key]()
    CASE, MODEL = cfg["case"], cfg["model"]
    REF_LIG, REF_PRO = cfg["ref_lig"], cfg["ref_pro"]
    DOCK_LIG, DOCK_PRO = cfg["dock_lig"], cfg["dock_pro"]
    for p in (REF_LIG, REF_PRO, DOCK_LIG, DOCK_PRO):
        if not p.exists():
            raise FileNotFoundError(p)

    mol1 = H.pdb_to_rdkit_with_prody(str(REF_LIG), "B")
    mol2 = H.pdb_to_rdkit_with_prody(str(DOCK_LIG), "B")
    prot1 = H.pdb_chain_to_rdkit_mol(str(REF_PRO), "A")
    prot2 = H.pdb_chain_to_rdkit_mol(str(DOCK_PRO), "A")

    dock_to_ref, assignment2, rmsd, alpha2, BCI2, attempts, bci_ok = (
        LM.align_ligands_otmol_bci_constrained(mol1, mol2, max_attempts=20)
    )
    otmol_diag = LM.otmol_mapping_diagnostics(mol1, mol2, assignment2)

    iface_ca = H.calculate_interface_rmsd_otmol(
        mol1, mol2, prot1, prot2, dock_to_ref, protein_atom_names=H.INTERFACE_PROTEIN_CA
    )
    iface_bb = H.calculate_interface_rmsd_otmol(
        mol1, mol2, prot1, prot2, dock_to_ref, protein_atom_names=H.INTERFACE_PROTEIN_BACKBONE
    )

    sub1, sub2 = H.find_close_residues_otmol(mol1, mol2, prot1, prot2, dock_to_ref, threshold=5.0)
    common_res = H.find_common_contact(sub1, sub2)
    fnat = H.calculate_fnat_backbone_sidechain_otmol(mol1, mol2, prot1, prot2, dock_to_ref, threshold=5.0)

    rec: dict = {
        "pdb_id": CASE,
        "af_model_id": MODEL,
        "ligand_rmsd": rmsd,
        "interface_rmsd": iface_ca["interface_rmsd"],
        "interface_rmsd_backbone": iface_bb["interface_rmsd"],
        "otmol_BCI": BCI2,
        "otmol_map_pairs": otmol_diag["map_pairs"],
        "residue_fnat": len(common_res) / len(sub1) if sub1 else -1,
        "fnat_all": fnat["fnat_all"],
        "fnat_backbone": fnat["fnat_backbone"],
        "fnat_sidechain": fnat["fnat_sidechain"],
        "ref_contacts_all_list": "; ".join(sorted(fnat["ref_contacts_all_set"])),
        "pred_contacts_all_list": "; ".join(sorted(fnat["pred_contacts_all_set"])),
        "ref_contacts_backbone_list": "; ".join(sorted(fnat["ref_contacts_backbone_set"])),
        "pred_contacts_backbone_list": "; ".join(sorted(fnat["pred_contacts_backbone_set"])),
        "ref_contacts_sidechain_list": "; ".join(sorted(fnat["ref_contacts_sidechain_set"])),
        "pred_contacts_sidechain_list": "; ".join(sorted(fnat["pred_contacts_sidechain_set"])),
        "residue_ref_contacts_list": "; ".join(sorted(sub1)),
        "residue_pred_contacts_list": "; ".join(sorted(sub2)),
        "residue_common_contacts_list": "; ".join(sorted(common_res)),
        "tm_TMscore": 1.0,
    }
    for short, col in {
        "bond_lengths": True,
        "bond_angles": True,
        "aromatic_ring_flatness": True,
        "double_bond_flatness": True,
        "internal_steric_clash": True,
        "energy_ratio": 50.0,
        "minimum_distance_to_protein": True,
        "volume_overlap_with_protein": True,
    }.items():
        rec[f"pb_{short}"] = col

    import contextlib
    import io

    buf_err = io.StringIO()
    with contextlib.redirect_stderr(buf_err):
        rec.update(
            BR.compute_binana_recall_row(
            str(REF_LIG),
            str(REF_PRO),
            str(DOCK_LIG),
            str(DOCK_PRO),
            BINANA,
            mol_ref=mol1,
            mol_dock=mol2,
            otmol_assignment=assignment2,
            dock_to_ref=dock_to_ref,
            otmol_BCI=BCI2,
            otmol_bci_zero_ok=bci_ok,
            include_detail_columns=False,
            )
        )

    scored = score_edockq_stepwise(pd.DataFrame([rec]))
    row = scored.iloc[0]

    out = {
        "case": f"{CASE} model {MODEL}",
        "geometry": {
            "ligand_rmsd": round(float(row["ligand_rmsd"]), 4),
            "interface_rmsd_CA": round(float(row["interface_rmsd"]), 4),
            "fnat_BBSC": round(float(row["fnat"]), 4),
            "kernel_i": round(float(row["kernel_i"]), 4),
            "kernel_l": round(float(row["kernel_l"]), 4),
            "eDockQ1": round(float(row["eDockQ1"]), 4),
        },
        "term_4_posebusters": round(float(row["term_4"]), 4),
        "TMscore": round(float(row["TMscore"]), 4),
        "eDockQ6_by_level": {},
    }

    for i, lv in enumerate(LEVELS, start=1):
        out["eDockQ6_by_level"][f"eDockQ6_{i}_{lv}"] = {
            "B6_six_pooled_micro_f1": round(float(row[f"fnat_binana6_{lv}"]), 6),
            "H_hydrophobic_f1": round(float(row[f"term_binana_hydrophobic_f1_{lv}"]), 4),
            "B6xH": round(float(row[f"fnat_binana_{lv}"]), 6),
            "eDockQ6": round(float(row[f"eDockQ6_{i}"]), 4),
        }

    per_type = {}
    for t in ["hbond", "salt", "pipi", "tstack", "cationpi", "halogen", "hydrophobic"]:
        per_type[t] = {
            lv: round(float(row[f"binana_{t}_{lv}_f1"]), 4)
            if pd.notna(row.get(f"binana_{t}_{lv}_f1"))
            else None
            for lv in LEVELS
        }
    out["per_type_f1_by_level"] = per_type

    print(json.dumps(out, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
