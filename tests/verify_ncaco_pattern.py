#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Verify the interface RMSD protein side using N, CA, C, O *all four together* (one pattern).

Cross-checks the OTMol-path pairing against an INDEPENDENT gemmi parse of the raw PDBs:
  * every interface residue contributes exactly the four backbone atoms N, CA, C, O,
  * each paired atom maps ref<->dock on the SAME (resnum, resname, atom_name),
  * the coordinates extracted by the pipeline equal the raw PDB coordinates (ref & dock).

Run:
    python pliq_ver5/tests/verify_ncaco_pattern.py
"""
from __future__ import annotations

import sys
from collections import Counter, defaultdict
from pathlib import Path

import gemmi
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from pliq.edockq import helpers as H

DATA_ROOT = ROOT.parent / "posebench" / "hpc_upload_pliq_processed_20260602_singlelig"
PATTERN = ("N", "CA", "C", "O")

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


def gemmi_backbone_index(pdb_path, chain_id="A"):
    """Independent raw parse: (resnum, resname, atomname) -> (x, y, z) for N/CA/C/O of chain A."""
    st = gemmi.read_structure(str(pdb_path))
    idx = {}
    model = st[0]
    for chain in model:
        if chain.name != chain_id:
            continue
        for res in chain:
            for atom in res:
                anm = atom.name.strip()
                if anm not in PATTERN:
                    continue
                idx[(res.seqid.num, res.name.strip(), anm)] = (
                    atom.pos.x, atom.pos.y, atom.pos.z
                )
    return idx


def labeled_pairs(prot_ref, prot_dock, mol_ref, mol_dock, threshold, atom_names):
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

    pairs = labeled_pairs(prot_ref, prot_dock, mol_ref, mol_dock, 10.0, PATTERN)

    # Independent raw coordinates
    g_ref = gemmi_backbone_index(ref_pro, "A")
    g_dock = gemmi_backbone_index(pred_pro, "A")

    # group by residue to confirm the 4-together pattern per residue
    per_res = defaultdict(set)
    coord_ok = True
    key_ok = True
    for k, r_xyz, d_xyz in pairs:
        per_res[(k[0], k[1])].add(k[2])
        if k not in g_ref or k not in g_dock:
            key_ok = False
            continue
        if not (np.allclose(r_xyz, g_ref[k], atol=1e-3) and np.allclose(d_xyz, g_dock[k], atol=1e-3)):
            coord_ok = False

    counts = Counter(k[2] for (k, _r, _d) in pairs)
    full_quartets = sum(1 for s in per_res.values() if s == set(PATTERN))
    partial = {res: sorted(s) for res, s in per_res.items() if s != set(PATTERN)}

    rmsd = H.calculate_interface_rmsd_otmol(
        mol_ref, mol_dock, prot_ref, prot_dock,
        {},  # ligand mapping irrelevant to this protein-side check; supply empty
        threshold=10.0, protein_atom_names=PATTERN,
    )

    print("=" * 74)
    print(f"CASE {dataset}/{pdb_id} model_{model_k}")
    print(f"  interface residues                 : {len(per_res)}")
    print(f"  atom counts (pattern N,CA,C,O)     : "
          f"N={counts.get('N',0)} CA={counts.get('CA',0)} C={counts.get('C',0)} O={counts.get('O',0)}"
          f"  total={len(pairs)}")
    print(f"  residues with all 4 (N,CA,C,O)     : {full_quartets}/{len(per_res)}")
    if partial:
        print(f"  residues missing some backbone atom: {partial}")
    print(f"  ref/dock keyed by SAME residue+atom: {key_ok}")
    print(f"  coords == independent gemmi parse  : {coord_ok}")
    print(f"  protein backbone atoms used (RMSD) : {rmsd['n_protein_atoms']}")

    # show one sample residue with its 4 paired atoms
    if pairs:
        sample_res = next(iter(sorted(per_res)))
        print(f"  sample residue {sample_res}:")
        for k, r_xyz, d_xyz in pairs:
            if (k[0], k[1]) == sample_res:
                d = float(np.linalg.norm(np.array(r_xyz) - np.array(d_xyz)))
                print(f"      {k[2]:>2s}  ref={tuple(round(v,2) for v in r_xyz)}  "
                      f"dock={tuple(round(v,2) for v in d_xyz)}  |Δ|={d:.3f}")


def main():
    for ds, pid, k in CASES:
        run_case(ds, pid, k)
    print("=" * 74)


if __name__ == "__main__":
    main()
