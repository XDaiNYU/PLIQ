# -*- coding: utf-8 -*-
"""CSV-stored remap tables must rebuild BINANA ligand indices offline."""
from __future__ import annotations

import json
import warnings
from pathlib import Path

import pytest

from pliq.config import resolve_binana_python
from pliq.edockq import binana_recall as BR
from pliq.edockq import helpers as H
from pliq.edockq import ligand_mapping as LM

ROOT = Path(__file__).resolve().parents[1]
EXAMPLE = ROOT / "examples" / "7ZU2_DHT_model15"
BINANA = resolve_binana_python()

SRC = {
    "rl": EXAMPLE / "ref" / "7ZU2_DHT_ligand_chainB.pdb",
    "rp": EXAMPLE / "ref" / "7ZU2_DHT_protein_renum.pdb",
    "dl": EXAMPLE / "dock" / "posebusters_7ZU2_DHT_model_15_aligned_ligand.pdb",
    "dp": EXAMPLE / "dock" / "posebusters_7ZU2_DHT_model_15_aligned_protein.pdb",
}


@pytest.mark.skipif(not EXAMPLE.is_dir(), reason="7ZU2 example missing")
def test_pliq_binana_remap_csv_columns_roundtrip():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        mol_ref = H.pdb_to_rdkit_with_prody(str(SRC["rl"]), "B")
        mol_dock = H.pdb_to_rdkit_with_prody(str(SRC["dl"]), "B")
        d2r, asg, *_ = LM.align_ligands_otmol_bci_constrained(mol_ref, mol_dock, max_attempts=10)
        row = BR.compute_binana_recall_row(
            str(SRC["rl"]),
            str(SRC["rp"]),
            str(SRC["dl"]),
            str(SRC["dp"]),
            BINANA,
            mol_ref=mol_ref,
            mol_dock=mol_dock,
            otmol_assignment=asg,
            dock_to_ref=d2r,
            include_detail_columns=True,
            obabel_ph=None,
            strip_h=True,
        )

    for col in LM.PLIQ_BINANA_REMAP_CSV_COLUMNS:
        assert col in row, col
    assert row["pliq_binana_remap_schema"] == LM.PLIQ_BINANA_REMAP_SCHEMA
    assert row["ligand_dock_to_ref_pairs"]
    assert row["binana_pdbindex_to_rdkit_ref_pairs"]
    assert row["binana_pdbindex_to_rdkit_dock_pairs"]

    maps = LM.parse_pliq_binana_remap_from_csv_row(row)
    assert len(maps["dock_to_ref"]) == len(d2r)
    assert len(maps["binana_to_rdkit_ref"]) > 0
    assert len(maps["binana_to_rdkit_dock"]) > 0

    ref_raw = json.loads(row["binana_hbond_ref_raw"])
    assert ref_raw, "expected hbond ref_raw interactions"
    lig_atoms = ref_raw[0]["ligandAtoms"]
    offline = LM.remap_binana_ligand_atom_indices_from_csv_row(
        lig_atoms, side="ref", row=row
    )
    live = LM.ligand_ref_rdkit_indices_for_binana(
        lig_atoms,
        side="ref",
        mol_ref=mol_ref,
        mol_dock=mol_dock,
        dock_to_ref=d2r,
        binana_pdbindex_to_rdkit_ref=maps["binana_to_rdkit_ref"],
        binana_pdbindex_to_rdkit_dock=maps["binana_to_rdkit_dock"],
    )
    assert offline == sorted(set(live or []))
