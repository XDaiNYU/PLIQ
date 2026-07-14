# -*- coding: utf-8 -*-
"""BINANA scheme A (as-is) vs scheme B (strip all H)."""
from __future__ import annotations

from pathlib import Path

import pytest

from pliq.config import resolve_binana_python
from pliq.edockq import binana_recall as BR
from pliq.edockq import helpers as H
from pliq.edockq import ligand_mapping as LM

ROOT = Path(__file__).resolve().parents[2]
LOCAL = ROOT / "data_posebusters"
BINANA = resolve_binana_python()


@pytest.fixture(scope="module")
def zu2_paths():
    base = LOCAL / "7ZU2_DHT"
    return dict(
        rl=base / "7ZU2_DHT_ref" / "7ZU2_DHT_ligand_chainB.pdb",
        rp=base / "7ZU2_DHT_ref" / "7ZU2_DHT_protein_renum.pdb",
        dl=base / "7ZU2_DHT_dock" / "posebusters_7ZU2_DHT_model_15_aligned_ligand.pdb",
        dp=base / "7ZU2_DHT_dock" / "posebusters_7ZU2_DHT_model_15_aligned_protein.pdb",
    )


def _binana_row(paths, *, strip_h: bool):
    mol1 = H.pdb_to_rdkit_with_prody(str(paths["rl"]), "B")
    mol2 = H.pdb_to_rdkit_with_prody(str(paths["dl"]), "B")
    d2r, asg, *_ = LM.align_ligands_otmol_bci_constrained(mol1, mol2, max_attempts=20)
    return BR.compute_binana_recall_row(
        str(paths["rl"]),
        str(paths["rp"]),
        str(paths["dl"]),
        str(paths["dp"]),
        BINANA,
        mol_ref=mol1,
        mol_dock=mol2,
        otmol_assignment=asg,
        dock_to_ref=d2r,
        otmol_BCI=0.0,
        otmol_bci_zero_ok=True,
        include_detail_columns=False,
        obabel_ph=None,
        strip_h=strip_h,
    )


@pytest.mark.skipif(not (LOCAL / "7ZU2_DHT").is_dir(), reason="7ZU2_DHT test structures missing")
def test_binana_h_mode_defaults(zu2_paths):
    row = _binana_row(zu2_paths, strip_h=False)
    assert row["binana_h_mode"] == "as_is"


@pytest.mark.skipif(not (LOCAL / "7ZU2_DHT").is_dir(), reason="7ZU2_DHT test structures missing")
def test_binana_strip_h_mode(zu2_paths):
    row = _binana_row(zu2_paths, strip_h=True)
    assert row["binana_h_mode"] == "strip_all"


@pytest.mark.skipif(not (LOCAL / "7ZU2_DHT").is_dir(), reason="7ZU2_DHT test structures missing")
def test_binana_strip_h_same_when_no_h(zu2_paths):
    a = _binana_row(zu2_paths, strip_h=False)
    b = _binana_row(zu2_paths, strip_h=True)
    for k in ("binana_hbond_tp_n", "binana_hydrophobic_tp_n"):
        assert a[k] == b[k]
