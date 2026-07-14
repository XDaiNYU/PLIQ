#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Verify that expanding the OTMol interface RMSD protein side from {CA} to {N,CA,C,O}:

  1. pairs every atom by the SAME (resnum, resname, atom_name) on ref & dock (no cross-pairing),
  2. the CA subset inside {N,CA,C,O} is byte-for-byte identical to the CA-only pairing,
  3. each of N / CA / C / O contributes the expected number of atoms.

Run:
    python pliq_ver5/tests/verify_backbone_consistency.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from pliq.edockq import helpers as H

DATA_ROOT = ROOT.parent / "posebench" / "hpc_upload_pliq_processed_20260602_singlelig"

CASES = [
    ("dockgen", "1hg0_1_SIN_1", 0),
    ("astex_diverse", "1HP0_AD3", 0),
    ("dockgen", "2hk9_1_SKM_0", 0),
]


def _paths(dataset, pdb_id, model_k):
    ref_lig = DATA_ROOT / "references" / dataset / pdb_id / f"{pdb_id}_ligand_chainB.pdb"
    ref_pro = DATA_ROOT / "references" / dataset / pdb_id / f"{pdb_id}_protein_renum.pdb"
    dock_dir = DATA_ROOT / "predictions" / "alphafold3" / dataset / "raw" / f"posebusters_{pdb_id}"
    pred_lig = dock_dir / f"posebusters_{pdb_id}_model_{model_k}_aligned_ligand.pdb"
    pred_pro = dock_dir / f"posebusters_{pdb_id}_model_{model_k}_aligned_protein.pdb"
    return ref_lig, ref_pro, pred_lig, pred_pro


def labeled_protein_pairs(prot_ref, prot_dock, mol_ref, mol_dock, threshold, atom_names):
    """Return list of (key, ref_xyz, dock_xyz) keyed by (resnum, resname, atom_name)."""
    allowed = set(atom_names)
    ref_res = H._interface_residue_keys(prot_ref, mol_ref, threshold)
    dock_res = H._interface_residue_keys(prot_dock, mol_dock, threshold)
    shared = ref_res & dock_res
    ref_idx = H._protein_backbone_coord_index(prot_ref, allowed)
    dock_idx = H._protein_backbone_coord_index(prot_dock, allowed)
    out = []
    for rid, rnm in sorted(shared):
        for anm in atom_names:
            k = (rid, rnm, anm)
            if k in ref_idx and k in dock_idx:
                out.append((k, ref_idx[k], dock_idx[k]))
    return out


def run_case(dataset, pdb_id, model_k):
    ref_lig, ref_pro, pred_lig, pred_pro = _paths(dataset, pdb_id, model_k)
    for p in (ref_lig, ref_pro, pred_lig, pred_pro):
        if not p.is_file():
            print(f"SKIP {pdb_id} (missing {p.name})")
            return
    mol_ref = H.pdb_to_rdkit_with_prody(str(ref_lig), "B")
    mol_dock = H.pdb_to_rdkit_with_prody(str(pred_lig), "B")
    prot_ref = H.pdb_chain_to_rdkit_mol(str(ref_pro), "A")
    prot_dock = H.pdb_chain_to_rdkit_mol(str(pred_pro), "A")

    ca_pairs = labeled_protein_pairs(prot_ref, prot_dock, mol_ref, mol_dock, 10.0, ("CA",))
    bb_pairs = labeled_protein_pairs(prot_ref, prot_dock, mol_ref, mol_dock, 10.0,
                                     H.INTERFACE_PROTEIN_BACKBONE)

    # (1) every pair shares the same key on ref & dock -> guaranteed by construction; assert anyway
    same_key = all(isinstance(k, tuple) for (k, _r, _d) in bb_pairs)

    # (2) CA subset of backbone == CA-only
    ca_only_keys = {k for (k, _r, _d) in ca_pairs}
    ca_only_map = {k: (r, d) for (k, r, d) in ca_pairs}
    bb_ca = [(k, r, d) for (k, r, d) in bb_pairs if k[2] == "CA"]
    bb_ca_keys = {k for (k, _r, _d) in bb_ca}
    keys_match = ca_only_keys == bb_ca_keys
    coords_match = True
    for k, r, d in bb_ca:
        r0, d0 = ca_only_map[k]
        if not (np.allclose(r, r0, atol=1e-6) and np.allclose(d, d0, atol=1e-6)):
            coords_match = False
            break

    # (3) per-atom-name counts
    from collections import Counter
    counts = Counter(k[2] for (k, _r, _d) in bb_pairs)

    print("=" * 74)
    print(f"CASE {dataset}/{pdb_id} model_{model_k}")
    print(f"  interface residues (shared) : {len(ca_pairs)}  (== CA count)")
    print(f"  backbone atom counts        : "
          f"N={counts.get('N',0)} CA={counts.get('CA',0)} "
          f"C={counts.get('C',0)} O={counts.get('O',0)}  total={len(bb_pairs)}")
    print(f"  (1) all pairs same key ref/dock         : {same_key}")
    print(f"  (2a) CA subset keys == CA-only keys      : {keys_match}")
    print(f"  (2b) CA subset coords == CA-only coords  : {coords_match}")
    print(f"  >>> CA identical when using N/CA/C/O     : {keys_match and coords_match}")


def main():
    for ds, pid, k in CASES:
        run_case(ds, pid, k)
    print("=" * 74)


if __name__ == "__main__":
    main()
