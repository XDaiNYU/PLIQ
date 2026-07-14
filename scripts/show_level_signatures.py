#!/usr/bin/env python3
"""Show, for one real case, how the 4 BINANA levels label ref/dock signatures.

For each interaction type it prints the ref signature multiset and the dock
signature multiset under each receptor mode (default / rxatm / rxidx / rxres),
plus the TP (multiset intersection), so the level-3 difference is visible.
"""
from __future__ import annotations

import contextlib
import io
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pliq.edockq import binana_recall as BR
from pliq.edockq import helpers as H
from pliq.edockq import ligand_mapping as LM
from pliq.edockq.otmol_runtime import ensure_local_otmol_on_path

ensure_local_otmol_on_path()

CASE = sys.argv[1] if len(sys.argv) > 1 else "1B9J"
MODEL = int(sys.argv[2]) if len(sys.argv) > 2 else 0
BASE = ROOT.parent / "LEADS_PEP" / "pliq_input_models0_49_unfixed_ref"
REF_LIG = BASE / "reference_inputs" / CASE / "ligand.pdb"
REF_PRO = BASE / "reference_inputs" / CASE / "protein.pdb"
DOCK_LIG = BASE / f"posebusters_{CASE}" / f"posebusters_{CASE}_model_{MODEL}_aligned_ligand.pdb"
DOCK_PRO = BASE / f"posebusters_{CASE}" / f"posebusters_{CASE}_model_{MODEL}_aligned_protein.pdb"
BINANA = ROOT / "pliq" / "_vendor" / "binana_python"

# (label, receptor_mode)
LEVELS = [
    ("eDockQ6_1 default", LM.RECEPTOR_SIG_CHAIN_ATOM),
    ("eDockQ6_2 rxatm  ", LM.RECEPTOR_SIG_RID_ATOM),
    ("eDockQ6_3 rxidx  ", LM.RECEPTOR_SIG_RID_ATOM_IDX),
    ("eDockQ6_4 rxres  ", LM.RECEPTOR_SIG_RID_RES),
]

# Only show full signatures for these types; the summary table covers all.
SHOW_KEYS = [("hydrogenBonds", "hbond"), ("hydrophobicContacts", "hydrophobic")]
ALL_KEYS = [
    ("hydrogenBonds", "hbond"),
    ("saltBridges", "salt"),
    ("piPiStackingInteractions", "pipi"),
    ("tStackingInteractions", "tstack"),
    ("cationPiInteractions", "cationpi"),
    ("halogenBonds", "halogen"),
    ("hydrophobicContacts", "hydrophobic"),
]


def f1_of(rc, dc):
    tp = sum(min(rc[s], dc.get(s, 0)) for s in rc)
    ref_n, dock_n = sum(rc.values()), sum(dc.values())
    if ref_n == 0:
        return None, ref_n, dock_n, tp
    prec = tp / dock_n if dock_n else 0.0
    rec = tp / ref_n
    f1 = (2 * prec * rec / (prec + rec)) if (prec + rec) else 0.0
    return f1, ref_n, dock_n, tp


def sig_counter(
    records,
    *,
    side,
    mol_ref,
    mol_dock,
    dock_to_ref,
    remap,
    mode,
    b2r_dock,
    b2r_ref,
):
    c: Counter = Counter()
    for x in records or []:
        if not isinstance(x, dict):
            continue
        s = LM.interaction_signature_remapped(
            x,
            mol_dock,
            mol_ref,
            dock_to_ref,
            remap,
            receptor_mode=mode,
            binana_pdbindex_to_rdkit_ref=b2r_ref,
            binana_pdbindex_to_rdkit_dock=b2r_dock,
            side=side,
        )
        c[s] += 1
    return c


def main() -> None:
    mol1 = H.pdb_to_rdkit_with_prody(str(REF_LIG), "B")
    mol2 = H.pdb_to_rdkit_with_prody(str(DOCK_LIG), "B")

    dock_to_ref, assignment2, *_ = LM.align_ligands_otmol_bci_constrained(
        mol1, mol2, max_attempts=20
    )

    buf = io.StringIO()
    with contextlib.redirect_stderr(buf), contextlib.redirect_stdout(buf):
        ref_js = BR.run_binana_collect(str(REF_LIG), str(REF_PRO), BINANA)
        dock_js = BR.run_binana_collect(str(DOCK_LIG), str(DOCK_PRO), BINANA)
        ref_atoms = BR.load_binana_ligand_atom_table(str(REF_LIG), BINANA)
        dock_atoms = BR.load_binana_ligand_atom_table(str(DOCK_LIG), BINANA)
        b2r_ref = LM.build_binana_pdbindex_to_rdkit_map(ref_atoms, mol1)
        b2r_dock = LM.build_binana_pdbindex_to_rdkit_map(dock_atoms, mol2)

    print(f"=== {CASE} model {MODEL} ===")
    print(
        f"OTMol dock_to_ref pairs: {len(dock_to_ref)}   "
        f"binana->rdkit ref: {len(b2r_ref)} dock: {len(b2r_dock)}\n"
    )

    # ---- summary table: F1 (TP/ref_n/dock_n) per type per level ----
    print("SUMMARY  F1 [TP/ref/dock]   (levels differ ONLY by receptor granularity)")
    header = f"{'type':<12}" + "".join(f"{lbl.strip():<22}" for lbl, _ in LEVELS)
    print(header)
    for json_key, short in ALL_KEYS:
        cells = []
        for _, mode in LEVELS:
            rc = sig_counter(ref_js.get(json_key), side="ref", mol_ref=mol1, mol_dock=mol2,
                             dock_to_ref=dock_to_ref, remap=True, mode=mode,
                             b2r_ref=b2r_ref, b2r_dock=b2r_dock)
            dc = sig_counter(dock_js.get(json_key), side="dock", mol_ref=mol1, mol_dock=mol2,
                             dock_to_ref=dock_to_ref, remap=True, mode=mode,
                             b2r_ref=b2r_ref, b2r_dock=b2r_dock)
            f1, ref_n, dock_n, tp = f1_of(rc, dc)
            f1s = "  NA " if f1 is None else f"{f1:.3f}"
            cells.append(f"{f1s} [{tp}/{ref_n}/{dock_n}]")
        print(f"{short:<12}" + "".join(f"{c:<22}" for c in cells))
    print()

    for json_key, short in SHOW_KEYS:
        print(f"################ {short}  ({json_key}) ################")
        for label, mode in LEVELS:
            rc = sig_counter(
                ref_js.get(json_key), side="ref", mol_ref=mol1, mol_dock=mol2,
                dock_to_ref=dock_to_ref, remap=True, mode=mode,
                b2r_ref=b2r_ref, b2r_dock=b2r_dock,
            )
            dc = sig_counter(
                dock_js.get(json_key), side="dock", mol_ref=mol1, mol_dock=mol2,
                dock_to_ref=dock_to_ref, remap=True, mode=mode,
                b2r_ref=b2r_ref, b2r_dock=b2r_dock,
            )
            tp = sum(min(rc[s], dc.get(s, 0)) for s in rc)
            ref_n = sum(rc.values())
            dock_n = sum(dc.values())
            prec = tp / dock_n if dock_n else float("nan")
            rec = tp / ref_n if ref_n else float("nan")
            f1 = (2 * prec * rec / (prec + rec)) if (prec and rec and prec + rec) else 0.0
            print(f"\n--- {label} | ref_n={ref_n} dock_n={dock_n} TP={tp} F1={f1:.4f} ---")
            common = set(rc) & set(dc)
            print("  REF signatures:")
            for s in sorted(rc):
                mark = "  <= MATCH" if s in common else ""
                print(f"    [{rc[s]}] {s}{mark}")
            print("  DOCK signatures:")
            for s in sorted(dc):
                mark = "  <= MATCH" if s in common else ""
                print(f"    [{dc[s]}] {s}{mark}")
        print()


if __name__ == "__main__":
    main()
