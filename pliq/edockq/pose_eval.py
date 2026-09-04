# -*- coding: utf-8 -*-
"""Evaluate a single docking pose from four PDB files (e-DockQ + BINANA)."""
from __future__ import annotations

import os
from pathlib import Path
from typing import Any, Dict, Optional, Union

from . import binana_recall as BR
from . import helpers as H
from . import ligand_mapping as LM
from .binana_recall import DEFAULT_BINANA_OBABEL_PH, DEFAULT_BINANA_STRIP_H, resolve_binana_h_settings
from .otmol_runtime import ensure_local_otmol_on_path

ensure_local_otmol_on_path()


def _load_ligand_mol(pdb_path: Union[str, Path], chain_id: str):
    temp_lig = None
    try:
        fd, temp_lig = os.mkstemp(suffix=".pdb")
        os.close(fd)
        return H.pdb_to_rdkit_with_prody(str(pdb_path), chain_id, temp_lig)
    except Exception:
        return H.pdb_to_rdkit_with_prody(str(pdb_path), chain_id)
    finally:
        if temp_lig and os.path.exists(temp_lig):
            try:
                os.remove(temp_lig)
            except OSError:
                pass


def evaluate_pose_from_pdbs(
    ref_lig_pdb: Union[str, Path],
    ref_pro_pdb: Union[str, Path],
    dock_lig_pdb: Union[str, Path],
    dock_pro_pdb: Union[str, Path],
    *,
    pdb_id: str = "",
    af_model_id: int = 0,
    ref_lig_chain: str = "B",
    ref_pro_chain: str = "A",
    dock_lig_chain: str = "B",
    dock_pro_chain: str = "A",
    run_binana: bool = True,
    binana_python: Optional[Union[str, Path]] = None,
    binana_obabel_ph: Optional[float] = DEFAULT_BINANA_OBABEL_PH,
    binana_strip_h: bool = DEFAULT_BINANA_STRIP_H,
) -> Dict[str, Any]:
    """
    Compute e-DockQ raw metrics (+ BINANA with built-in default cutoffs) for one pose.

    Ligand mapping uses OTMol only; reflection is always disabled in alignment presets.

    BINANA hydrogen (Open Babel before PDB→PDBQT):
    - Default: ``obabel -d`` strips H from ref+dock protein+ligand.
    - ``binana_obabel_ph=pH``: use ``obabel -p`` at that pH instead (no -d).
    - ``binana_strip_h=False`` with ``binana_obabel_ph=None``: keep inputs as-is (legacy).
    """
    ref_lig_pdb = Path(ref_lig_pdb).resolve()
    ref_pro_pdb = Path(ref_pro_pdb).resolve()
    dock_lig_pdb = Path(dock_lig_pdb).resolve()
    dock_pro_pdb = Path(dock_pro_pdb).resolve()
    for p in (ref_lig_pdb, ref_pro_pdb, dock_lig_pdb, dock_pro_pdb):
        if not p.is_file():
            raise FileNotFoundError(p)

    pdb_file1 = str(ref_lig_pdb)
    prot_pdb_file1 = str(ref_pro_pdb)
    pdb_file2 = str(dock_lig_pdb)
    prot_pdb_file2 = str(dock_pro_pdb)

    mol1 = _load_ligand_mol(ref_lig_pdb, ref_lig_chain)
    mol2 = _load_ligand_mol(dock_lig_pdb, dock_lig_chain)

    ligand1_heavy_atoms = mol1.GetNumAtoms() if mol1 else -1
    ligand2_heavy_atoms = mol2.GetNumAtoms() if mol2 else -1
    alpha2, BCI2 = None, None
    rmsd = None
    dock_to_ref: dict = {}
    assignment2 = None
    otmol_bci_zero_ok = None
    otmol_align_attempts = None

    if mol1 and mol2:
        dock_to_ref, assignment2, rmsd, alpha2, BCI2, otmol_align_attempts, otmol_bci_zero_ok = (
            LM.align_ligands_otmol_bci_constrained(mol1, mol2, max_attempts=20)
        )
        otmol_diag = LM.otmol_mapping_diagnostics(mol1, mol2, assignment2)
    else:
        otmol_diag = LM.otmol_mapping_diagnostics(None, None, None)

    prot1 = H.pdb_chain_to_rdkit_mol(prot_pdb_file1, ref_pro_chain)
    prot2 = H.pdb_chain_to_rdkit_mol(prot_pdb_file2, dock_pro_chain)
    if not prot1 or not prot2 or not mol1 or not mol2:
        raise ValueError("Missing protein or ligand")

    iface_ca = H.calculate_interface_rmsd_otmol(
        mol1, mol2, prot1, prot2, dock_to_ref,
        threshold=10.0, protein_atom_names=H.INTERFACE_PROTEIN_CA,
    )
    iface_bb = H.calculate_interface_rmsd_otmol(
        mol1, mol2, prot1, prot2, dock_to_ref,
        threshold=10.0, protein_atom_names=H.INTERFACE_PROTEIN_BACKBONE,
    )

    sub_results1, sub_results2 = H.find_close_residues_otmol(
        mol1, mol2, prot1, prot2, dock_to_ref, threshold=5.0
    )
    common_contacts = H.find_common_contact(sub_results1, sub_results2)
    residue_fnat = len(common_contacts) / len(sub_results1) if sub_results1 else -1

    fnat_results = H.calculate_fnat_backbone_sidechain_otmol(
        mol1, mol2, prot1, prot2, dock_to_ref, threshold=5.0
    )
    common_all = fnat_results["common_all"]
    common_backbone = fnat_results["common_backbone"]
    common_sidechain = fnat_results["common_sidechain"]

    rec: Dict[str, Any] = {
        "pdb_id": pdb_id or ref_lig_pdb.stem,
        "af_model_id": af_model_id,
        "ligand_rmsd": rmsd,
        "otmol_alpha": alpha2,
        "otmol_BCI": BCI2,
        "otmol_bci_zero_ok": otmol_bci_zero_ok,
        "otmol_align_attempts": otmol_align_attempts,
        "otmol_map_ref_atoms": otmol_diag["ref_atoms"],
        "otmol_map_dock_atoms": otmol_diag["dock_atoms"],
        "otmol_map_pairs": otmol_diag["map_pairs"],
        "otmol_map_unmapped_dock": otmol_diag["unmapped_dock"],
        "otmol_map_atomicnum_mismatch": otmol_diag["atomicnum_mismatch"],
        "ligand_dock_to_ref_pairs": LM.format_ligand_dock_to_ref_pairs(dock_to_ref),
        "ligand_dock_to_ref_map": LM.format_ligand_dock_to_ref_map_json(mol1, mol2, dock_to_ref)
        if mol1 and mol2 and dock_to_ref
        else "[]",
        "ligand1_heavy_atoms": ligand1_heavy_atoms,
        "ligand2_heavy_atoms": ligand2_heavy_atoms,
        "interface_rmsd": iface_ca["interface_rmsd"],
        "interface_rmsd_backbone": iface_bb["interface_rmsd"],
        "interface_rmsd_src": "otmol",
        "residue_fnat": residue_fnat,
        "fnat_map_src": fnat_results.get("fnat_map_src", "otmol"),
        "fnat_mapped_pairs": fnat_results.get("fnat_mapped_pairs", len(dock_to_ref)),
        "fnat_all": fnat_results["fnat_all"],
        "fnat_backbone": fnat_results["fnat_backbone"],
        "fnat_sidechain": fnat_results["fnat_sidechain"],
        "tp_all": fnat_results["tp_all"],
        "tp_backbone": fnat_results["tp_backbone"],
        "tp_sidechain": fnat_results["tp_sidechain"],
        "total_ref_all": fnat_results["total_ref_all"],
        "total_ref_backbone": fnat_results["total_ref_backbone"],
        "total_ref_sidechain": fnat_results["total_ref_sidechain"],
        "ref_contacts_all": fnat_results["ref_contacts_all"],
        "pred_contacts_all": fnat_results["pred_contacts_all"],
        "ref_contacts_backbone": fnat_results["ref_contacts_backbone"],
        "pred_contacts_backbone": fnat_results["pred_contacts_backbone"],
        "ref_contacts_sidechain": fnat_results["ref_contacts_sidechain"],
        "pred_contacts_sidechain": fnat_results["pred_contacts_sidechain"],
        "residue_ref_contacts": len(sub_results1),
        "residue_pred_contacts": len(sub_results2),
        "residue_common_contacts": len(common_contacts),
        "ref_contacts_all_list": "; ".join(sorted(fnat_results["ref_contacts_all_set"])),
        "pred_contacts_all_list": "; ".join(sorted(fnat_results["pred_contacts_all_set"])),
        "common_contacts_all_list": "; ".join(sorted(common_all)),
        "ref_contacts_backbone_list": "; ".join(sorted(fnat_results["ref_contacts_backbone_set"])),
        "pred_contacts_backbone_list": "; ".join(sorted(fnat_results["pred_contacts_backbone_set"])),
        "common_contacts_backbone_list": "; ".join(sorted(common_backbone)),
        "ref_contacts_sidechain_list": "; ".join(sorted(fnat_results["ref_contacts_sidechain_set"])),
        "pred_contacts_sidechain_list": "; ".join(sorted(fnat_results["pred_contacts_sidechain_set"])),
        "common_contacts_sidechain_list": "; ".join(sorted(common_sidechain)),
        "residue_ref_contacts_list": "; ".join(sorted(sub_results1)),
        "residue_pred_contacts_list": "; ".join(sorted(sub_results2)),
        "residue_common_contacts_list": "; ".join(sorted(common_contacts)),
    }
    rec.update(
        LM.pliq_mapping_csv_columns(
            mol_ref=mol1,
            mol_dock=mol2,
            dock_to_ref=dock_to_ref,
            bci=BCI2,
            bci_zero_ok=otmol_bci_zero_ok,
        )
    )
    rec["binana_error"] = ""
    bpath = (
        Path(binana_python)
        if binana_python
        else Path(__file__).resolve().parent.parent / "_vendor" / "binana_python"
    )
    if run_binana:
        try:
            rec.update(
                BR.compute_binana_recall_row(
                    pdb_file1,
                    prot_pdb_file1,
                    pdb_file2,
                    prot_pdb_file2,
                    bpath,
                    mol_ref=mol1,
                    mol_dock=mol2,
                    otmol_assignment=assignment2,
                    dock_to_ref=dock_to_ref,
                    otmol_BCI=BCI2,
                    otmol_bci_zero_ok=otmol_bci_zero_ok,
                    get_all_interactions_kwargs={},
                    obabel_ph=binana_obabel_ph,
                    strip_h=binana_strip_h,
                )
            )
        except Exception as e_bin:
            rec.update(BR.empty_binana_row())
            rec["binana_error"] = str(e_bin)[:500]
            rec.update(
                LM.pliq_mapping_csv_columns(
                    mol_ref=mol1,
                    mol_dock=mol2,
                    dock_to_ref=dock_to_ref,
                    bci=BCI2,
                    bci_zero_ok=otmol_bci_zero_ok,
                )
            )
    else:
        rec.update(BR.empty_binana_row())
        rec["binana_error"] = "disabled"
    return rec
