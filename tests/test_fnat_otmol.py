#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Validate OTMol-based FNAT on real prepared PoseBench cases."""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from pliq.edockq import helpers as H
from pliq.edockq import ligand_mapping as LM

DATA_ROOT = (
    ROOT.parent
    / "posebench"
    / "hpc_upload_pliq_processed_20260602_singlelig"
)

CASES = [
    ("dockgen", "1hg0_1_SIN_1", 0),
    ("astex_diverse", "1HP0_AD3", 0),
    ("dockgen", "2hk9_1_SKM_0", 0),
]


def _paths(dataset: str, pdb_id: str, model_k: int):
    ref_lig = DATA_ROOT / "references" / dataset / pdb_id / f"{pdb_id}_ligand_chainB.pdb"
    ref_pro = DATA_ROOT / "references" / dataset / pdb_id / f"{pdb_id}_protein_renum.pdb"
    dock_dir = DATA_ROOT / "predictions" / "alphafold3" / dataset / "raw" / f"posebusters_{pdb_id}"
    pred_lig = dock_dir / f"posebusters_{pdb_id}_model_{model_k}_aligned_ligand.pdb"
    pred_pro = dock_dir / f"posebusters_{pdb_id}_model_{model_k}_aligned_protein.pdb"
    return ref_lig, ref_pro, pred_lig, pred_pro


def _load_pair(dataset: str, pdb_id: str, model_k: int = 0):
    ref_lig, ref_pro, pred_lig, pred_pro = _paths(dataset, pdb_id, model_k)
    for p in (ref_lig, ref_pro, pred_lig, pred_pro):
        if not p.is_file():
            pytest.skip(f"missing {p}")
    mol_ref = H.pdb_to_rdkit_with_prody(str(ref_lig), "B")
    mol_dock = H.pdb_to_rdkit_with_prody(str(pred_lig), "B")
    prot_ref = H.pdb_chain_to_rdkit_mol(str(ref_pro), "A")
    prot_dock = H.pdb_chain_to_rdkit_mol(str(pred_pro), "A")
    assert mol_ref and mol_dock and prot_ref and prot_dock
    return mol_ref, mol_dock, prot_ref, prot_dock


@pytest.mark.parametrize("dataset,pdb_id,model_k", CASES)
def test_otmol_fnat_runs_and_maps_ligand_only(dataset, pdb_id, model_k):
    mol_ref, mol_dock, prot_ref, prot_dock = _load_pair(dataset, pdb_id, model_k)
    dock_to_ref, _assign, rmsd, _a, _bci, _n, _ok = LM.align_ligands_otmol_bci_constrained(
        mol_ref, mol_dock, max_attempts=5
    )
    assert dock_to_ref, "OTMol mapping empty"
    assert rmsd is not None
    n_ref = mol_ref.GetNumAtoms()
    n_dock = mol_dock.GetNumAtoms()
    for d_idx, r_idx in dock_to_ref.items():
        assert 0 <= d_idx < n_dock, f"bad dock idx {d_idx}"
        assert 0 <= r_idx < n_ref, f"bad ref idx {r_idx}"

    fnat = H.calculate_fnat_backbone_sidechain_otmol(
        mol_ref, mol_dock, prot_ref, prot_dock, dock_to_ref, threshold=5.0
    )
    assert fnat["fnat_map_src"] == "otmol"
    assert fnat["fnat_mapped_pairs"] == len(dock_to_ref)
    assert fnat["tp_all"] <= fnat["total_ref_all"]
    assert fnat["common_all"].issubset(fnat["ref_contacts_all_set"])
    assert fnat["common_all"].issubset(fnat["pred_contacts_all_set"])
    # contact tags must be ref-atom keyed
    for tag in fnat["ref_contacts_all_set"]:
        assert tag.startswith("mol_ref"), tag


@pytest.mark.parametrize("dataset,pdb_id,model_k", CASES[:2])
def test_otmol_fnat_valid_range(dataset, pdb_id, model_k):
    """OTMol FNAT must produce valid (in-range) values; mapping source is OTMol-only."""
    mol_ref, mol_dock, prot_ref, prot_dock = _load_pair(dataset, pdb_id, model_k)
    dock_to_ref, *_ = LM.align_ligands_otmol_bci_constrained(mol_ref, mol_dock, max_attempts=5)
    otm = H.calculate_fnat_backbone_sidechain_otmol(
        mol_ref, mol_dock, prot_ref, prot_dock, dock_to_ref, threshold=5.0
    )
    assert otm["fnat_map_src"] == "otmol"
    assert otm["total_ref_all"] >= 0
    if otm["total_ref_all"] > 0:
        assert 0.0 <= otm["fnat_all"] <= 1.0
        assert otm["tp_all"] <= otm["total_ref_all"]
