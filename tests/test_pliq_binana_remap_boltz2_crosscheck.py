# -*- coding: utf-8 -*-
"""Cross-check remap CSV + OTMol against PoseBuster (Boltz-2) stepwise reference."""
from __future__ import annotations

import warnings
from pathlib import Path

import pandas as pd
import pytest

from pliq.config import resolve_binana_python
from pliq.edockq import binana_recall as BR
from pliq.edockq import helpers as H
from pliq.edockq import ligand_mapping as LM

ROOT = Path(__file__).resolve().parents[1]
EDOCKQ = ROOT.parent
DATA = EDOCKQ / "data_posebusters"
BOLTZ_CSV = (
    EDOCKQ
    / "PLIQ_ver5_fixed_refined_collect"
    / "Posebusters_Boltz2"
    / "e-DockQ_stepwise_all.csv"
)
BINANA = resolve_binana_python()


@pytest.fixture(scope="module")
def boltz_df():
    if not BOLTZ_CSV.is_file():
        pytest.skip(f"missing {BOLTZ_CSV}")
    return pd.read_csv(BOLTZ_CSV, low_memory=False)


@pytest.fixture(autouse=True)
def _disable_pdbqt_verify(monkeypatch):
    """Historical Boltz-2 runs used obabel -p 7.4 without PDBQT key verification."""
    monkeypatch.setattr(BR, "verify_pdbqt_preserves_residue_keys", lambda *a, **k: None)


def _paths(case: str, k: int) -> tuple[Path, Path, Path, Path]:
    ref_lig = DATA / case / f"{case}_ligand_chainB.pdb"
    ref_pro = DATA / case / f"{case}_protein_renum.pdb"
    dock_dir = DATA / f"posebusters_{case}"
    dock_lig = dock_dir / f"posebusters_{case}_model_{k}_aligned_ligand.pdb"
    dock_pro = dock_dir / f"posebusters_{case}_model_{k}_aligned_protein.pdb"
    return ref_lig, ref_pro, dock_lig, dock_pro


@pytest.mark.skipif(not DATA.is_dir(), reason="data_posebusters missing")
def test_boltz2_5S8I_2LY_m0_otmol_matches_reference(boltz_df):
    case, k = "5S8I_2LY", 0
    ref_lig, ref_pro, dock_lig, dock_pro = _paths(case, k)
    for p in (ref_lig, ref_pro, dock_lig, dock_pro):
        if not p.is_file():
            pytest.skip(f"missing {p}")

    ref = boltz_df[(boltz_df.pdb_id == case) & (boltz_df.af_model_id == k)].iloc[0]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        mol_ref = H.pdb_to_rdkit_with_prody(str(ref_lig), "B")
        mol_dock = H.pdb_to_rdkit_with_prody(str(dock_lig), "B")
        d2r, asg, rmsd, _a, bci, _n, _ok = LM.align_ligands_otmol_bci_constrained(
            mol_ref, mol_dock, max_attempts=20
        )

    assert abs(float(rmsd) - float(ref.ligand_rmsd)) < 1e-6
    assert float(bci) == float(ref.otmol_BCI)


@pytest.mark.skipif(not DATA.is_dir(), reason="data_posebusters missing")
def test_boltz2_5S8I_2LY_m0_csv_remap_rebuilds_live_signatures(boltz_df):
    case, k = "5S8I_2LY", 0
    ref_lig, ref_pro, dock_lig, dock_pro = _paths(case, k)
    for p in (ref_lig, ref_pro, dock_lig, dock_pro):
        if not p.is_file():
            pytest.skip(f"missing {p}")

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        mol_ref = H.pdb_to_rdkit_with_prody(str(ref_lig), "B")
        mol_dock = H.pdb_to_rdkit_with_prody(str(dock_lig), "B")
        d2r, asg, *_ = LM.align_ligands_otmol_bci_constrained(mol_ref, mol_dock, max_attempts=20)
        row = BR.compute_binana_recall_row(
            str(ref_lig),
            str(ref_pro),
            str(dock_lig),
            str(dock_pro),
            BINANA,
            mol_ref=mol_ref,
            mol_dock=mol_dock,
            otmol_assignment=asg,
            dock_to_ref=d2r,
            obabel_ph=7.4,
            strip_h=False,
            include_detail_columns=True,
        )

    maps = LM.parse_pliq_binana_remap_from_csv_row(row)
    assert len(maps["dock_to_ref"]) == len(d2r)
    assert maps["binana_to_rdkit_ref"] and maps["binana_to_rdkit_dock"]

    mode = LM.RECEPTOR_SIG_RID_ATOM
    for side in ("ref", "dock"):
        for _json_key, short in BR.BINANA_INTERACTION_KEYS:
            raw = __import__("json").loads(row[f"binana_{short}_{side}_raw"])
            for inter in raw:
                offline = LM.rebuild_interaction_signature_from_csv_row(
                    inter, side=side, row=row, receptor_mode=mode
                )
                live = LM.interaction_signature_remapped(
                    inter,
                    mol_dock,
                    mol_ref,
                    d2r,
                    True,
                    receptor_mode=mode,
                    binana_pdbindex_to_rdkit_ref=maps["binana_to_rdkit_ref"],
                    binana_pdbindex_to_rdkit_dock=maps["binana_to_rdkit_dock"],
                    side=side,
                )
                assert offline == live

    LM.validate_offline_binana_remap_from_csv_row(row, receptor_mode=mode)
