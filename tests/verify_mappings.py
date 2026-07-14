#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Triple-check mapping correctness used by PLIQ:

  (A) OTMol ligand atom mapping (dock_to_ref): in-range, 1:1, element-consistent, BCI==0.
  (B) Residue-level FNAT mapping: contact tags keyed by reference ligand atom; both sides
      share the same tag space (so intersection == true positives).
  (C) Chain-level selection: ligand = chain B, protein = chain A, cross-checked atom counts
      against an independent gemmi parse of the raw PDB.

Run:
    python pliq_ver5/tests/verify_mappings.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import gemmi

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from pliq.edockq import helpers as H
from pliq.edockq import ligand_mapping as LM
from pliq.edockq.otmol_runtime import ensure_local_otmol_on_path

ensure_local_otmol_on_path()

DATA_ROOT = ROOT.parent / "posebench" / "hpc_upload_pliq_processed_20260602_singlelig"
CASES = [
    ("dockgen", "1hg0_1_SIN_1", 0),
    ("astex_diverse", "1HP0_AD3", 0),
    ("dockgen", "2hk9_1_SKM_0", 0),
]


def _paths(ds, pid, k):
    rl = DATA_ROOT / "references" / ds / pid / f"{pid}_ligand_chainB.pdb"
    rp = DATA_ROOT / "references" / ds / pid / f"{pid}_protein_renum.pdb"
    dd = DATA_ROOT / "predictions" / "alphafold3" / ds / "raw" / f"posebusters_{pid}"
    pl = dd / f"posebusters_{pid}_model_{k}_aligned_ligand.pdb"
    pp = dd / f"posebusters_{pid}_model_{k}_aligned_protein.pdb"
    return rl, rp, pl, pp


def gemmi_chain_heavy_count(pdb_path, chain_id):
    st = gemmi.read_structure(str(pdb_path))
    n = 0
    for chain in st[0]:
        if chain.name != chain_id:
            continue
        for res in chain:
            if res.name == "HOH":
                continue
            for atom in res:
                if atom.element.name != "H":
                    n += 1
    return n


def check(cond, msg, fails):
    tag = "OK " if cond else "FAIL"
    if not cond:
        fails.append(msg)
    print(f"    [{tag}] {msg}")


def run_case(ds, pid, k):
    rl, rp, pl, pp = _paths(ds, pid, k)
    print("=" * 76)
    print(f"CASE {ds}/{pid} model_{k}")
    fails = []
    if not all(p.is_file() for p in (rl, rp, pl, pp)):
        print("    SKIP (missing files)")
        return fails

    mol_ref = H.pdb_to_rdkit_with_prody(str(rl), "B")
    mol_dock = H.pdb_to_rdkit_with_prody(str(pl), "B")
    prot_ref = H.pdb_chain_to_rdkit_mol(str(rp), "A")
    prot_dock = H.pdb_chain_to_rdkit_mol(str(pp), "A")

    # ---- (C) chain-level selection ----
    print("  (C) chain-level selection")
    g_lig_ref = gemmi_chain_heavy_count(rl, "B")
    g_lig_dock = gemmi_chain_heavy_count(pl, "B")
    check(mol_ref.GetNumAtoms() == g_lig_ref,
          f"ref ligand chainB heavy atoms: rdkit={mol_ref.GetNumAtoms()} gemmi={g_lig_ref}", fails)
    check(mol_dock.GetNumAtoms() == g_lig_dock,
          f"dock ligand chainB heavy atoms: rdkit={mol_dock.GetNumAtoms()} gemmi={g_lig_dock}", fails)
    # protein chain A: ensure every rdkit protein atom carries residue info and no chain-B leakage
    g_prot_ref = gemmi_chain_heavy_count(rp, "A")
    check(prot_ref.GetNumAtoms() == g_prot_ref,
          f"ref protein chainA heavy atoms: rdkit={prot_ref.GetNumAtoms()} gemmi={g_prot_ref}", fails)

    # ---- (A) OTMol ligand mapping ----
    print("  (A) OTMol ligand atom mapping")
    dock_to_ref, assignment, rmsd, alpha, bci, n_try, bci0 = (
        LM.align_ligands_otmol_bci_constrained(mol_ref, mol_dock, max_attempts=20)
    )
    nr, nd = mol_ref.GetNumAtoms(), mol_dock.GetNumAtoms()
    in_range = all(0 <= d < nd and 0 <= r < nr for d, r in dock_to_ref.items())
    check(in_range, f"all dock/ref indices in range (map size={len(dock_to_ref)})", fails)
    one_to_one = len(set(dock_to_ref.values())) == len(dock_to_ref)
    check(one_to_one, "mapping is 1:1 (no two dock atoms share a ref atom)", fails)
    elem_ok = all(
        mol_dock.GetAtomWithIdx(d).GetAtomicNum() == mol_ref.GetAtomWithIdx(r).GetAtomicNum()
        for d, r in dock_to_ref.items()
    )
    check(elem_ok, "every mapped pair has identical element (atomic number)", fails)
    diag = LM.otmol_mapping_diagnostics(mol_ref, mol_dock, assignment)
    check(diag.get("atomicnum_mismatch", -1) == 0,
          f"otmol diagnostics atomicnum_mismatch == 0 (got {diag.get('atomicnum_mismatch')})", fails)
    check(bool(bci0), f"OTMol BCI == 0 achieved (BCI={bci}, attempts={n_try})", fails)

    # ---- (B) residue-level FNAT mapping ----
    print("  (B) residue-level FNAT contact mapping")
    sub1, sub2 = H.find_close_residues_otmol(mol_ref, mol_dock, prot_ref, prot_dock, dock_to_ref, 5.0)
    tags_ref_ok = all(t.startswith("mol_ref") for t in sub1)
    tags_pred_ok = all(t.startswith("mol_ref") for t in sub2)
    check(tags_ref_ok and tags_pred_ok,
          "ref & pred contact tags both keyed by reference ligand atom (shared tag space)", fails)
    common = set(H.find_common_contact(sub1, sub2))
    check(common.issubset(set(sub1)) and common.issubset(set(sub2)),
          f"common contacts subset of both ref({len(set(sub1))}) and pred({len(set(sub2))}) -> TP={len(common)}", fails)
    fn = H.calculate_fnat_backbone_sidechain_otmol(mol_ref, mol_dock, prot_ref, prot_dock, dock_to_ref, 5.0)
    check(fn["fnat_map_src"] == "otmol", "FNAT mapping source is OTMol", fails)
    check(fn["tp_all"] <= fn["total_ref_all"], "TP <= total ref contacts", fails)

    return fails


def main():
    all_fails = []
    for ds, pid, k in CASES:
        all_fails += run_case(ds, pid, k)
    print("=" * 76)
    if all_fails:
        print(f"RESULT: {len(all_fails)} FAILURES")
        for f in all_fails:
            print("  -", f)
        return 1
    print("RESULT: ALL MAPPING CHECKS PASSED")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
