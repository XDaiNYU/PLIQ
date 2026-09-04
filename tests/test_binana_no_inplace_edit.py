# -*- coding: utf-8 -*-
"""BINANA prep/PDBQT must preserve (resseq, resname, atomname) and not touch source files."""
from __future__ import annotations

import hashlib
import shutil
import warnings
from pathlib import Path

import pytest

from pliq.config import resolve_binana_python
from pliq.edockq import binana_recall as BR
from pliq.edockq import helpers as H
from pliq.edockq import ligand_mapping as LM

ROOT = Path(__file__).resolve().parents[1]
EXAMPLE = ROOT / "examples" / "7ZU2_DHT_model15"
NCAA = ROOT.parent / "NCAA_case" / "hpc_upload_pliq_ncaa_corrected" / "pliq_inputs_after_noligH" / "references" / "1nm6_B"
BINANA = resolve_binana_python()

SRC_7ZU2 = {
    "rl": EXAMPLE / "ref" / "7ZU2_DHT_ligand_chainB.pdb",
    "rp": EXAMPLE / "ref" / "7ZU2_DHT_protein_renum.pdb",
    "dl": EXAMPLE / "dock" / "posebusters_7ZU2_DHT_model_15_aligned_ligand.pdb",
    "dp": EXAMPLE / "dock" / "posebusters_7ZU2_DHT_model_15_aligned_protein.pdb",
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _run_binana_strip(paths: dict[str, Path]) -> None:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        mol1 = H.pdb_to_rdkit_with_prody(str(paths["rl"]), "B")
        mol2 = H.pdb_to_rdkit_with_prody(str(paths["dl"]), "B")
        d2r, asg, *_ = LM.align_ligands_otmol_bci_constrained(mol1, mol2, max_attempts=20)
        BR.compute_binana_recall_row(
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
            strip_h=True,
        )


@pytest.fixture
def zu2_copy(tmp_path):
    if not all(p.is_file() for p in SRC_7ZU2.values()):
        pytest.skip("7ZU2_DHT_model15 example PDBs missing")
    out = {}
    for key, src in SRC_7ZU2.items():
        sub = "ref" if key.startswith("r") else "dock"
        dst = tmp_path / sub / src.name
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)
        out[key] = dst
    out["hashes"] = {k: _sha256(p) for k, p in out.items()}
    return out


@pytest.mark.parametrize(
    "pdb_path,chain",
    [
        (EXAMPLE / "ref" / "7ZU2_DHT_protein_renum.pdb", "A"),
        (EXAMPLE / "ref" / "7ZU2_DHT_ligand_chainB.pdb", "B"),
        pytest.param(
            NCAA / "1nm6_B_protein_renum.pdb",
            "A",
            marks=pytest.mark.skipif(not (NCAA / "1nm6_B_protein_renum.pdb").is_file(), reason="no NCAA"),
        ),
        pytest.param(
            NCAA / "1nm6_B_ligand_chainB.pdb",
            "B",
            marks=pytest.mark.skipif(not (NCAA / "1nm6_B_ligand_chainB.pdb").is_file(), reason="no NCAA"),
        ),
    ],
)
def test_binana_prep_and_pdbqt_preserve_residue_keys(tmp_path, pdb_path, chain):
    import tempfile

    with tempfile.TemporaryDirectory() as td:
        td_path = Path(td)
        prep = BR._prepare_pdb_for_binana(str(pdb_path), chain, td_path, "entity", strip_h=True)
        assert BR.heavy_residue_atom_keys_from_pdb(str(pdb_path)) == BR.heavy_residue_atom_keys_from_pdb(
            prep
        ), "chain-forced temp copy must keep resseq/resname/atomname"
        pdbqt = td_path / "entity.pdbqt"
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            BR._run_obabel(prep, str(pdbqt), obabel_ph=None)
        BR.verify_pdbqt_preserves_residue_keys(prep, str(pdbqt), BINANA, label=pdb_path.name)


@pytest.mark.skipif(not EXAMPLE.is_dir(), reason="7ZU2 example missing")
def test_binana_strip_d_does_not_modify_source_pdbs(zu2_copy):
    before = zu2_copy["hashes"]
    paths = {k: zu2_copy[k] for k in SRC_7ZU2}
    _run_binana_strip(paths)
    for key, path in paths.items():
        assert _sha256(path) == before[key], f"{path.name} was modified on disk"
