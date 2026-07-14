"""BINANA ligand signatures must use one OTMol dock_to_ref for both ref and dock poses."""
from __future__ import annotations

import contextlib
import io
import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pliq.edockq import binana_recall as BR
from pliq.edockq import helpers as H
from pliq.edockq import ligand_mapping as LM
from pliq.edockq.otmol_runtime import ensure_local_otmol_on_path
from pliq.binana_interaction_criteria import load_binana_get_all_interaction_kwargs_from_json

ensure_local_otmol_on_path()

REPO = ROOT.parent
CASE_BASE = REPO / "data_posebusters" / "7ZU2_DHT"
BINANA = ROOT / "pliq" / "_vendor" / "binana_python"
CRIT = ROOT / "interaction_criteria.json"


@pytest.mark.skipif(not CASE_BASE.is_dir(), reason="7ZU2_DHT fixture not present")
def test_binana_ref_and_dock_share_otmol_ligand_indices():
    """Ref/dock THR hbond must share the same L[ref_rdkit_idx] after unified OTMol mapping."""
    mol_ref = H.pdb_to_rdkit_with_prody(str(CASE_BASE / "7ZU2_DHT_ref" / "7ZU2_DHT_ligand_chainB.pdb"), "B")
    mol_dock = H.pdb_to_rdkit_with_prody(
        str(CASE_BASE / "7ZU2_DHT_dock" / "posebusters_7ZU2_DHT_model_15_aligned_ligand.pdb"), "B"
    )
    dock_to_ref, assignment, *_ = LM.align_ligands_otmol_bci_constrained(mol_ref, mol_dock, max_attempts=20)
    assert dock_to_ref

    gkw = load_binana_get_all_interaction_kwargs_from_json(CRIT)
    with contextlib.redirect_stderr(io.StringIO()):
        row = BR.compute_binana_recall_row(
            str(CASE_BASE / "7ZU2_DHT_ref" / "7ZU2_DHT_ligand_chainB.pdb"),
            str(CASE_BASE / "7ZU2_DHT_ref" / "7ZU2_DHT_protein_renum.pdb"),
            str(CASE_BASE / "7ZU2_DHT_dock" / "posebusters_7ZU2_DHT_model_15_aligned_ligand.pdb"),
            str(CASE_BASE / "7ZU2_DHT_dock" / "posebusters_7ZU2_DHT_model_15_aligned_protein.pdb"),
            BINANA,
            mol_ref=mol_ref,
            mol_dock=mol_dock,
            otmol_assignment=assignment,
            dock_to_ref=dock_to_ref,
            otmol_BCI=0.0,
            otmol_bci_zero_ok=True,
            get_all_interactions_kwargs=gkw,
            obabel_ph=None,
            include_detail_columns=True,
        )

    assert row["binana_ligand_map_src"] == "otmol"
    assert row["binana_hbond_tp_n"] >= 1

    tp_sigs = {x["signature"] for x in json.loads(row["binana_hbond_tp_details"])}
    fn_sigs = {x["signature"] for x in json.loads(row["binana_hbond_fn_details"])}
    fp_sigs = {x["signature"] for x in json.loads(row["binana_hbond_fp_details"])}

    # Matched pairs: ref and dock must use identical L[...] prefix.
    for sig in tp_sigs:
        assert sig.startswith("L[")
        lig_part = sig.split("]R[", 1)[0] + "]"

    # No orphan ref/dock hbond with different ligand index for same receptor atom class.
    for fs in fn_sigs | fp_sigs:
        lig_part, rec_part = fs.split("]R[", 1)
        lig_part += "]"
        rec_key = rec_part.rstrip("]")
        counterpart = lig_part + "R[" + rec_key + "]"
        if fs in fn_sigs:
            assert counterpart not in fp_sigs, f"split hbond: {fs} vs {counterpart}"

    thr_tp = [s for s in tp_sigs if ":THR:" in s and ":200:" in s]
    assert thr_tp, f"expected THR A200 hbond TP, got {tp_sigs}"
    assert all(s.startswith("L[18]") for s in thr_tp)


def test_align_ligands_not_called_inside_compute_binana_recall_row(monkeypatch):
    """compute_binana_recall_row must consume upstream dock_to_ref, not re-align."""
    if not CASE_BASE.is_dir():
        pytest.skip("7ZU2_DHT fixture not present")

    mol_ref = H.pdb_to_rdkit_with_prody(str(CASE_BASE / "7ZU2_DHT_ref" / "7ZU2_DHT_ligand_chainB.pdb"), "B")
    mol_dock = H.pdb_to_rdkit_with_prody(
        str(CASE_BASE / "7ZU2_DHT_dock" / "posebusters_7ZU2_DHT_model_15_aligned_ligand.pdb"), "B"
    )
    dock_to_ref, assignment, *_ = LM.align_ligands_otmol_bci_constrained(mol_ref, mol_dock, max_attempts=5)

    def _boom(*_a, **_k):
        raise AssertionError("OTMol must not run inside compute_binana_recall_row")

    monkeypatch.setattr(LM, "align_ligands_otmol_bci_constrained", _boom)

    with contextlib.redirect_stderr(io.StringIO()):
        BR.compute_binana_recall_row(
            str(CASE_BASE / "7ZU2_DHT_ref" / "7ZU2_DHT_ligand_chainB.pdb"),
            str(CASE_BASE / "7ZU2_DHT_ref" / "7ZU2_DHT_protein_renum.pdb"),
            str(CASE_BASE / "7ZU2_DHT_dock" / "posebusters_7ZU2_DHT_model_15_aligned_ligand.pdb"),
            str(CASE_BASE / "7ZU2_DHT_dock" / "posebusters_7ZU2_DHT_model_15_aligned_protein.pdb"),
            BINANA,
            mol_ref=mol_ref,
            mol_dock=mol_dock,
            dock_to_ref=dock_to_ref,
            otmol_assignment=assignment,
            otmol_BCI=0.0,
            otmol_bci_zero_ok=True,
            get_all_interactions_kwargs=load_binana_get_all_interaction_kwargs_from_json(CRIT),
            obabel_ph=None,
            include_detail_columns=False,
        )
