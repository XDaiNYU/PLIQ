#!/usr/bin/env python3
"""Validate BINANA 6_3 (rxidx): ligand ref-RDKit idx + receptor amino identity (no protein idx)."""
from __future__ import annotations

import json
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from pliq.edockq import binana_recall as BR
from pliq.edockq import helpers as H
from pliq.edockq import ligand_mapping as LM

CASE = "1B9J"
MODEL = 0
BASE = Path(__file__).resolve().parents[2] / "LEADS_PEP" / "pliq_input_models0_49_unfixed_ref"
REF_LIG = BASE / "reference_inputs" / CASE / "ligand.pdb"
REF_PRO = BASE / "reference_inputs" / CASE / "protein.pdb"
DOCK_LIG = BASE / f"posebusters_{CASE}" / f"posebusters_{CASE}_model_{MODEL}_aligned_ligand.pdb"
DOCK_PRO = BASE / f"posebusters_{CASE}" / f"posebusters_{CASE}_model_{MODEL}_aligned_protein.pdb"
BINANA = ROOT / "pliq" / "_vendor" / "binana_python"


def _recall_for_mode(
    ref_js,
    dock_js,
    mol_ref,
    mol_dock,
    dock_to_ref,
    receptor_mode: str,
    binana_pdb_to_rdkit_dock,
    binana_pdb_to_rdkit_ref,
) -> dict[str, float | int]:
    out: dict[str, float | int] = {}
    for json_key, short in BR.BINANA_INTERACTION_KEYS:
        ref_c, dock_c = BR._counter_for_key_pair(
            ref_js,
            dock_js,
            json_key,
            mol_ref,
            mol_dock,
            dock_to_ref,
            bool(dock_to_ref),
            receptor_mode=receptor_mode,
            binana_pdbindex_to_rdkit=binana_pdb_to_rdkit_dock,
            binana_pdbindex_to_rdkit_ref=binana_pdb_to_rdkit_ref,
        )
        frac, ref_n, dock_n, tp = BR.multiset_recall(ref_c, dock_c)
        out[short] = {"frac": frac, "ref_n": ref_n, "dock_n": dock_n, "tp": tp}
    return out


def _sample_signatures(
    ref_js, dock_js, mol_ref, mol_dock, dock_to_ref, binana_pdb_to_rdkit_dock, binana_pdb_to_rdkit_ref, n=3
):
    key = "hydrogenBonds"
    print(f"\nSample 6_3 signatures ({key}, first {n} ref / dock):")
    for side, js in [("ref", ref_js), ("dock", dock_js)]:
        print(f"  -- {side} --")
        for inter in (js.get(key) or [])[:n]:
            if not isinstance(inter, dict):
                continue
            sig = LM.interaction_signature_remapped(
                inter,
                mol_dock,
                mol_ref,
                dock_to_ref,
                bool(dock_to_ref),
                receptor_mode=LM.RECEPTOR_SIG_RID_ATOM_IDX,
                binana_pdbindex_to_rdkit_ref=binana_pdb_to_rdkit_ref,
                binana_pdbindex_to_rdkit_dock=binana_pdb_to_rdkit_dock,
                side=side,
            )
            print(f"    {sig[:160]}")


def main() -> None:
    for p in (REF_LIG, REF_PRO, DOCK_LIG, DOCK_PRO):
        if not p.is_file():
            raise SystemExit(f"missing input: {p}")

    mol_ref = H.pdb_to_rdkit_with_prody(str(REF_LIG), "B")
    mol_dock = H.pdb_to_rdkit_with_prody(str(DOCK_LIG), "B")
    if not all([mol_ref, mol_dock]):
        raise SystemExit("failed to load RDKit mols")

    print(f"Loading {CASE} model {MODEL}...", flush=True)

    dock_to_ref, assignment, rmsd, alpha, bci, attempts, bci_ok = (
        LM.align_ligands_otmol_bci_constrained(mol_ref, mol_dock, max_attempts=5)
    )
    print(f"=== {CASE} model {MODEL} ===")
    print(f"ligand_rmsd={rmsd:.3f} Å  BCI={bci}  bci_zero_ok={bci_ok}  map_pairs={len(dock_to_ref)}")
    print(f"ligand_dock_to_ref_pairs (first 5): {LM.format_ligand_dock_to_ref_pairs(dock_to_ref)[:80]}...")

    print("Running BINANA on ref...", flush=True)
    ref_js = BR.run_binana_collect(str(REF_LIG), str(REF_PRO), BINANA)
    print("Running BINANA on dock...", flush=True)
    dock_js = BR.run_binana_collect(str(DOCK_LIG), str(DOCK_PRO), BINANA)
    dock_binana_atoms = BR.load_binana_ligand_atom_table(str(DOCK_LIG), BINANA)
    ref_binana_atoms = BR.load_binana_ligand_atom_table(str(REF_LIG), BINANA)
    binana_pdb_to_rdkit_dock = LM.build_binana_pdbindex_to_rdkit_map(dock_binana_atoms, mol_dock)
    binana_pdb_to_rdkit_ref = LM.build_binana_pdbindex_to_rdkit_map(ref_binana_atoms, mol_ref)

    old = _recall_for_mode(
        ref_js,
        dock_js,
        mol_ref,
        mol_dock,
        dock_to_ref,
        LM.RECEPTOR_SIG_RID_ATOM_IDX,
        binana_pdb_to_rdkit_dock,
        binana_pdb_to_rdkit_ref,
    )
    # Legacy rxidx used protein atomIndex on receptor; compare against 6_1 for sanity.
    edockq6_1 = _recall_for_mode(
        ref_js,
        dock_js,
        mol_ref,
        mol_dock,
        dock_to_ref,
        LM.RECEPTOR_SIG_CHAIN_ATOM,
        binana_pdb_to_rdkit_dock,
        binana_pdb_to_rdkit_ref,
    )

    print("\n6_3 (rxidx) per interaction type — ligand ref-RDKit idx + receptor amino identity:")
    print(f"{'type':12s} {'frac':>10s} {'ref_n':>6s} {'tp':>7s}  (6_1 frac for comparison)")
    pooled_tp = pooled_ref = 0
    for short in [s for _, s in BR.BINANA_INTERACTION_KEYS]:
        n = old[short]
        c1 = edockq6_1[short]
        print(
            f"{short:12s} {n['frac']:10.4f} {n['ref_n']:6d} {n['tp']:7d}  {c1['frac']:10.4f}"
        )
        if n["ref_n"] > 0:
            pooled_tp += n["tp"]
            pooled_ref += n["ref_n"]

    micro = pooled_tp / pooled_ref if pooled_ref else float("nan")
    print(f"\n6_3 micro recall (sum tp / sum ref): {micro:.4f}")

    _sample_signatures(
        ref_js, dock_js, mol_ref, mol_dock, dock_to_ref,
        binana_pdb_to_rdkit_dock, binana_pdb_to_rdkit_ref,
    )

    sample = json.loads(LM.format_ligand_dock_to_ref_map_json(mol_ref, mol_dock, dock_to_ref))[:2]
    print("\nSample ligand_dock_to_ref_map records:")
    print(json.dumps(sample, indent=2))


if __name__ == "__main__":
    main()
