# -*- coding: utf-8 -*-
"""Run e-DockQ pipeline over CSV cases and save results to CSV."""
from __future__ import annotations

import os
from typing import Optional

import pandas as pd
from tqdm import tqdm
from pathlib import Path

from . import binana_recall as BR
from .otmol_runtime import ensure_local_otmol_on_path

ensure_local_otmol_on_path()

from . import helpers as H
from . import ligand_mapping as LM


def run_edockq(
    csv_file,
    path_ref_lig,
    path_ref_pro,
    path_docked,
    out_csv,
    n_models=20,
    n_cases=0,
    csv_sep="\t",
    subfolder_col="Subfolder",
    run_binana=True,
    binana_python=None,
    interaction_criteria_xlsx=None,
    binana_criteria_tier="medium",
    binana_obabel_ph: Optional[float] = BR.DEFAULT_BINANA_OBABEL_PH,
    binana_strip_h: bool = BR.DEFAULT_BINANA_STRIP_H,
):
    """
    Run e-DockQ (Ligand RMSD with otmol, Interface RMSD, FNAT residue/backbone/sidechain) and save CSV.
    If run_binana is True, also runs BINANA on ref vs docked PDBs and adds interaction recall fractions
    (hydrogen bonds, salt bridges, pi-pi, T-stack, cation-pi, halogen, hydrophobic) using ligand/protein
    atom-pair signatures. Ligand RMSD and BINANA mapping use **OTMol only** (no MCS coordinate RMSD fallback).

    Cutoffs for BINANA are read from ``interaction_criteria.xlsx`` or ``interaction_criteria.json``
    (``binana_criteria_tier``: ``low`` = loose, ``high`` = strict, default ``medium``;
    legacy ``medium1`` = ``medium``) and passed into
    ``get_all_interactions`` without editing vendored BINANA scripts.

    BINANA hydrogen handling (no ``obabel -p`` by default):
    - ``binana_strip_h=False`` (scheme A): ref/dock PDBs as-is; dock keeps predictor H if present.
    - ``binana_strip_h=True`` (scheme B): ``obabel -d`` on ref+dock protein+ligand before BINANA.

    ``binana_obabel_ph``: optional Open Babel ``-p`` pH (usually leave ``None``).

    n_cases: 0 = all rows; else take first n_cases.
    """
    from pliq.binana_interaction_criteria import (
        DEFAULT_BINANA_TIER,
        load_binana_get_all_interaction_kwargs_from_path,
    )

    # None/{} means: keep BINANA built-in defaults (no external cutoff override).
    binana_gkwargs = {}
    if run_binana:
        if interaction_criteria_xlsx:
            crit_path = Path(interaction_criteria_xlsx).resolve()
            binana_gkwargs = load_binana_get_all_interaction_kwargs_from_path(
                crit_path, tier=str(binana_criteria_tier or DEFAULT_BINANA_TIER)
            )

    df1 = pd.read_csv(csv_file, sep=csv_sep)
    if n_cases and n_cases > 0:
        df1 = df1.iloc[: n_cases]
    df1 = df1.copy()
    if subfolder_col in df1.columns and "pdb_id" not in df1.columns:
        df1["pdb_id"] = df1[subfolder_col]

    results = []
    list_fail = []
    Path(out_csv).parent.mkdir(parents=True, exist_ok=True)

    for index, row in tqdm(df1.iterrows(), total=len(df1), desc="e-DockQ"):
        case_name = row["pdb_id"]

        for af_model_id in range(n_models):
            alpha2, BCI2 = None, None
            try:
                pdb_file1 = f"{path_ref_lig}/{case_name}/{case_name}_ligand_chainB.pdb"
                pdb_file2 = f"{path_docked}/posebusters_{case_name}/posebusters_{case_name}_model_{af_model_id}_aligned_ligand.pdb"
                prot_pdb_file1 = f"{path_ref_pro}/{case_name}/{case_name}_protein_renum.pdb"
                prot_pdb_file2 = f"{path_docked}/posebusters_{case_name}/posebusters_{case_name}_model_{af_model_id}_aligned_protein.pdb"

                temp_lig = None
                try:
                    fd, temp_lig = os.mkstemp(suffix=".pdb")
                    os.close(fd)
                    mol1 = H.pdb_to_rdkit_with_prody(pdb_file1, "B", temp_lig)
                except Exception:
                    mol1 = H.pdb_to_rdkit_with_prody(pdb_file1, "B")
                finally:
                    if temp_lig and os.path.exists(temp_lig):
                        try:
                            os.remove(temp_lig)
                        except Exception:
                            pass

                temp_lig2 = None
                try:
                    fd2, temp_lig2 = os.mkstemp(suffix=".pdb")
                    os.close(fd2)
                    mol2 = H.pdb_to_rdkit_with_prody(pdb_file2, "B", temp_lig2)
                except Exception:
                    mol2 = H.pdb_to_rdkit_with_prody(pdb_file2, "B")
                finally:
                    if temp_lig2 and os.path.exists(temp_lig2):
                        try:
                            os.remove(temp_lig2)
                        except Exception:
                            pass

                ligand1_heavy_atoms = mol1.GetNumAtoms() if mol1 else -1
                ligand2_heavy_atoms = mol2.GetNumAtoms() if mol2 else -1
                rmsd = None
                dock_to_ref: dict = {}
                assignment2 = None
                otmol_bci_zero_ok = None
                otmol_align_attempts = None

                if mol1 and mol2:
                    # All mappings (ligand RMSD, interface, FNAT, BINANA) are OTMol-only.
                    dock_to_ref, assignment2, rmsd, alpha2, BCI2, otmol_align_attempts, otmol_bci_zero_ok = (
                        LM.align_ligands_otmol_bci_constrained(mol1, mol2, max_attempts=20)
                    )
                    otmol_diag = LM.otmol_mapping_diagnostics(mol1, mol2, assignment2)
                else:
                    otmol_diag = LM.otmol_mapping_diagnostics(None, None, None)

                ligand_rmsd = rmsd

                prot1 = H.pdb_chain_to_rdkit_mol(prot_pdb_file1, "A")
                prot2 = H.pdb_chain_to_rdkit_mol(prot_pdb_file2, "A")
                if not prot1 or not prot2 or not mol1 or not mol2:
                    raise ValueError("Missing protein or ligand")
                # Interface RMSD via OTMol ligand mapping (consistent with ligand RMSD / FNAT).
                # Protein backbone atoms are paired by residue + atom name (no MCS); ligand atoms
                # are paired by the OTMol dock->ref assignment. CA-only reproduces the legacy CA
                # behaviour, while N/CA/C/O is now available without further changes.
                iface_ca = H.calculate_interface_rmsd_otmol(
                    mol1, mol2, prot1, prot2, dock_to_ref,
                    threshold=10.0, protein_atom_names=H.INTERFACE_PROTEIN_CA,
                )
                iface_bb = H.calculate_interface_rmsd_otmol(
                    mol1, mol2, prot1, prot2, dock_to_ref,
                    threshold=10.0, protein_atom_names=H.INTERFACE_PROTEIN_BACKBONE,
                )
                interface_rmsd_value = iface_ca["interface_rmsd"]
                interface_rmsd_backbone_value = iface_bb["interface_rmsd"]

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

                rec = {
                    "pdb_id": case_name,
                    "af_model_id": af_model_id,
                    "ligand_rmsd": ligand_rmsd,
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
                    "ligand_dock_to_ref_map": LM.format_ligand_dock_to_ref_map_json(
                        mol1, mol2, dock_to_ref
                    )
                    if mol1 and mol2 and dock_to_ref
                    else "[]",
                    "ligand1_heavy_atoms": ligand1_heavy_atoms,
                    "ligand2_heavy_atoms": ligand2_heavy_atoms,
                    "interface_rmsd": interface_rmsd_value,
                    "interface_rmsd_backbone": interface_rmsd_backbone_value,
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
                                get_all_interactions_kwargs=binana_gkwargs,
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
                results.append(rec)
                pd.DataFrame(results).to_csv(out_csv, index=False)
            except Exception as e:
                # skip_marked: never silently drop — record the case/model + reason.
                list_fail.append(
                    {"pdb_id": case_name, "af_model_id": af_model_id, "reason": str(e)[:500]}
                )

    pd.DataFrame(results).to_csv(out_csv, index=False)
    if list_fail:
        fail_csv = Path(out_csv).with_name(Path(out_csv).stem + "_failures.csv")
        pd.DataFrame(list_fail).to_csv(fail_csv, index=False)
        print(f"e-DockQ: {len(list_fail)} case/model failures recorded -> {fail_csv}")
    return results, list_fail
