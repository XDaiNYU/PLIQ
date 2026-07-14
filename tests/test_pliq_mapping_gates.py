"""Mapping quality flags for CSV (soft gate — mark, do not hard-fail the row)."""
from __future__ import annotations

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pliq.edockq import binana_recall as BR
from pliq.edockq import helpers as H
from pliq.edockq import ligand_mapping as LM
from pliq.edockq.otmol_runtime import ensure_local_otmol_on_path

ensure_local_otmol_on_path()

CASE_BASE = ROOT.parent / "data_posebusters" / "7ZU2_DHT"
BINANA = ROOT / "pliq" / "_vendor" / "binana_python"


def test_check_otmol_bci_zero_flags_nonzero():
    ok, msg = LM.check_otmol_bci_zero(0.12, bci_zero_ok=False)
    assert not ok
    assert "BCI" in msg


@pytest.mark.skipif(not CASE_BASE.is_dir(), reason="7ZU2_DHT fixture not present")
def test_posebusters_case_mapping_columns_all_ok():
    ref_lig = CASE_BASE / "7ZU2_DHT_ref" / "7ZU2_DHT_ligand_chainB.pdb"
    dock_lig = CASE_BASE / "7ZU2_DHT_dock" / "posebusters_7ZU2_DHT_model_15_aligned_ligand.pdb"
    mol_ref = H.pdb_to_rdkit_with_prody(str(ref_lig), "B")
    mol_dock = H.pdb_to_rdkit_with_prody(str(dock_lig), "B")
    d2r, asg, rmsd, alpha, bci, _attempts, bci_ok = LM.align_ligands_otmol_bci_constrained(
        mol_ref, mol_dock, max_attempts=20
    )
    assert bci_ok is True

    ref_table = BR.load_binana_ligand_atom_table(str(ref_lig), BINANA, obabel_ph=None)
    dock_table = BR.load_binana_ligand_atom_table(str(dock_lig), BINANA, obabel_ph=None)
    ref_map = LM.build_binana_pdbindex_to_rdkit_map(ref_table, mol_ref)
    dock_map = LM.build_binana_pdbindex_to_rdkit_map(dock_table, mol_dock)

    cols = LM.pliq_mapping_csv_columns(
        mol_ref=mol_ref,
        mol_dock=mol_dock,
        dock_to_ref=d2r,
        bci=bci,
        bci_zero_ok=bci_ok,
        ref_binana_atoms=ref_table,
        dock_binana_atoms=dock_table,
        binana_to_rdkit_ref=ref_map,
        binana_to_rdkit_dock=dock_map,
    )
    assert cols["pliq_mapping_trusted"] is True
    assert cols["otmol_dock_to_ref_complete_ok"] is True
    assert cols["binana_rdkit_ref_map_ok"] is True
    assert cols["binana_rdkit_dock_map_ok"] is True
    assert cols["binana_rdkit_ref_map_stats"] == f"{mol_ref.GetNumAtoms()}/{mol_ref.GetNumAtoms()}"


def test_incomplete_binana_map_not_ok():
    if not CASE_BASE.is_dir():
        pytest.skip("7ZU2_DHT fixture not present")
    mol = H.pdb_to_rdkit_with_prody(
        str(CASE_BASE / "7ZU2_DHT_ref" / "7ZU2_DHT_ligand_chainB.pdb"), "B"
    )
    table = BR.load_binana_ligand_atom_table(
        str(CASE_BASE / "7ZU2_DHT_ref" / "7ZU2_DHT_ligand_chainB.pdb"), BINANA, obabel_ph=None
    )
    partial = {k: v for k, v in list(LM.build_binana_pdbindex_to_rdkit_map(table, mol).items())[:-1]}
    st = LM.assess_binana_rdkit_map(table, mol, partial, side="ref")
    assert st["ok"] is False
    assert "unmapped" in st["message"]
