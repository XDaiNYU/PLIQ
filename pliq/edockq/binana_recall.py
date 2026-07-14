# -*- coding: utf-8 -*-
"""
BINANA interaction recall: compare reference vs docked poses using atom-pair signatures.

pliq ver5 adds per-instance TP / FP / FN detail records (JSON in CSV) with ligand atoms
in dock frame, reference-ligand frame (OTMol-mapped), and receptor atom metadata.
"""
from __future__ import annotations

import io
import json
import math
import subprocess
import sys
from collections import Counter, defaultdict
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Any, Dict, List, Optional, Tuple

from . import helpers as H
from . import ligand_mapping as LM

# (JSON key from binana.output.dictionary.collect_all, short prefix for CSV columns)
BINANA_INTERACTION_KEYS: List[Tuple[str, str]] = [
    ("hydrogenBonds", "hbond"),
    ("saltBridges", "salt"),
    ("piPiStackingInteractions", "pipi"),
    ("tStackingInteractions", "tstack"),
    ("cationPiInteractions", "cationpi"),
    ("halogenBonds", "halogen"),
    ("hydrophobicContacts", "hydrophobic"),
]

# Six-type subset used in fnat_binana6 (excludes hydrophobic).
BINANA_SIX_TYPE_KEYS: List[Tuple[str, str]] = [
    t for t in BINANA_INTERACTION_KEYS if t[1] != "hydrophobic"
]

# (receptor_mode passed to ligand_mapping.interaction_signature_remapped, CSV column suffix)
RECEPTOR_SIGNATURE_MODES: List[Tuple[str, str]] = [
    (LM.RECEPTOR_SIG_CHAIN_ATOM, ""),
    (LM.RECEPTOR_SIG_RID_ATOM_IDX, "_rxidx"),
    (LM.RECEPTOR_SIG_RID_ATOM, "_rxatm"),
    (LM.RECEPTOR_SIG_RID_RES, "_rxres"),
]

# BINANA input: no Open Babel -p by default (scheme A: keep dock H if present).
# Optional obabel_ph adds polar H; binana_strip_h=True uses scheme B (obabel -d on all inputs).
DEFAULT_BINANA_OBABEL_PH = None
DEFAULT_BINANA_STRIP_H = False


def _strip_hydrogens_obabel(in_pdb: str, out_pdb: str) -> None:
    """Remove all hydrogens (scheme B). Preserves heavy-atom residue IDs and types."""
    subprocess.run(
        ["obabel", in_pdb, "-O", out_pdb, "-d"],
        check=True,
        capture_output=True,
        text=True,
    )


def _prepare_pdb_for_binana(
    pdb_path: str,
    chain_id: str,
    workdir: Path,
    basename: str,
    *,
    strip_h: bool,
) -> str:
    """Force chain ID; optionally strip all H with obabel -d before PDB→PDBQT."""
    chained = (
        H.write_forced_chain_pdb(pdb_path, chain_id, workdir / f"{basename}_chain.pdb")
        or pdb_path
    )
    if not strip_h:
        return chained
    stripped = workdir / f"{basename}_noh.pdb"
    _strip_hydrogens_obabel(chained, str(stripped))
    return str(stripped)


def _run_obabel(in_pdb: str, out_pdbqt: str, *, obabel_ph: Optional[float]) -> None:
    cmd = ["obabel", in_pdb, "-O", out_pdbqt]
    if obabel_ph is not None:
        cmd.extend(["-p", str(obabel_ph)])
    subprocess.run(
        cmd,
        check=True,
        capture_output=True,
        text=True,
    )


def run_binana_collect(
    ligand_pdb: str,
    receptor_pdb: str,
    binana_python: Path,
    get_all_interactions_kwargs: Optional[Dict[str, Any]] = None,
    *,
    obabel_ph: Optional[float] = DEFAULT_BINANA_OBABEL_PH,
    strip_h: bool = DEFAULT_BINANA_STRIP_H,
) -> Dict[str, Any]:
    """Convert PDB→PDBQT with obabel, run BINANA, return collect_all dict."""
    binana_python = Path(binana_python).resolve()
    if not binana_python.is_dir():
        raise FileNotFoundError(f"BINANA python path not found: {binana_python}")

    with TemporaryDirectory() as td:
        td_path = Path(td)
        # Canonicalize chains; scheme B optionally strips all H (obabel -d) first.
        lig_src = _prepare_pdb_for_binana(ligand_pdb, "B", td_path, "ligand", strip_h=strip_h)
        rec_src = _prepare_pdb_for_binana(receptor_pdb, "A", td_path, "receptor", strip_h=strip_h)
        lig_qt = td_path / "ligand.pdbqt"
        rec_qt = td_path / "receptor.pdbqt"
        _run_obabel(lig_src, str(lig_qt), obabel_ph=obabel_ph)
        _run_obabel(rec_src, str(rec_qt), obabel_ph=obabel_ph)

        old_path = sys.path[:]
        sys.path.insert(0, str(binana_python))
        try:
            from binana import interactions, load_ligand_receptor
            from binana.output import dictionary

            buf_out, buf_err = io.StringIO(), io.StringIO()
            with redirect_stdout(buf_out), redirect_stderr(buf_err):
                lig_mol, rec_mol = load_ligand_receptor.from_files(
                    str(lig_qt), str(rec_qt)
                )
                kwargs = dict(get_all_interactions_kwargs or {})
                all_inf = interactions.get_all_interactions(lig_mol, rec_mol, **kwargs)
                js = dictionary.collect_all(all_inf)
        finally:
            sys.path[:] = old_path

    return js


def load_binana_ligand_atom_table(
    ligand_pdb: str,
    binana_python: Path,
    *,
    obabel_ph: Optional[float] = DEFAULT_BINANA_OBABEL_PH,
    strip_h: bool = DEFAULT_BINANA_STRIP_H,
) -> Dict[int, Dict[str, Any]]:
    """
    Load BINANA ligand atoms from obabel PDBQT (same path as ``run_binana_collect``).

    Keys are PDB serial numbers (``atomIndex`` in BINANA JSON), values hold
    coordinates and element for mapping to the RDKit dock mol.
    """
    binana_python = Path(binana_python).resolve()
    if not binana_python.is_dir():
        raise FileNotFoundError(f"BINANA python path not found: {binana_python}")

    with TemporaryDirectory() as td:
        td_path = Path(td)
        lig_src = _prepare_pdb_for_binana(ligand_pdb, "B", td_path, "ligand", strip_h=strip_h)
        lig_qt = td_path / "ligand.pdbqt"
        _run_obabel(lig_src, str(lig_qt), obabel_ph=obabel_ph)

        old_path = sys.path[:]
        sys.path.insert(0, str(binana_python))
        try:
            from binana._structure.mol import Mol

            lig_mol = Mol()
            lig_mol.load_pdb_file(str(lig_qt))
        finally:
            sys.path[:] = old_path

    table: Dict[int, Dict[str, Any]] = {}
    for _idx, atom in lig_mol.all_atoms.items():
        pdb_idx = LM._safe_pdb_index(atom.pdb_index)
        if pdb_idx is None:
            continue
        table[pdb_idx] = {
            "x": float(atom.coordinates.x),
            "y": float(atom.coordinates.y),
            "z": float(atom.coordinates.z),
            "element": str(atom.element or "").upper(),
            "atomName": str(atom.atom_name or "").strip(),
            "resName": str(atom.residue or "").strip(),
            "resID": int(atom.resid) if atom.resid is not None else None,
            "chain": str(atom.chain or "").strip(),
        }
    return table


def _safe_int(x: Any) -> Optional[int]:
    if x is None:
        return None
    try:
        return int(x)
    except (TypeError, ValueError):
        return None


def _normalize_binana_atom(ba: Dict[str, Any], *, frame: str) -> Dict[str, Any]:
    """Normalize one BINANA atom record for JSON / CSV detail output."""
    return {
        "chain": str(ba.get("chain") or "").strip(),
        "resID": _safe_int(ba.get("resID")),
        "resName": str(ba.get("resName") or "").strip(),
        "atomName": str(ba.get("atomName") or "").strip(),
        "atomIndex": _safe_int(ba.get("atomIndex")),
        "frame": frame,
    }


def _rdkit_atom_record(mol: Any, idx: int, *, frame: str) -> Dict[str, Any]:
    ch, rid, rnm, anm = LM.rdkit_pdb_tuple(mol.GetAtomWithIdx(idx))
    return {
        "chain": ch,
        "resID": rid,
        "resName": rnm,
        "atomName": anm,
        "atomIndex": idx,
        "frame": frame,
    }


def _ligand_atoms_detail(
    lig_atoms: List[Any],
    *,
    side: str,
    mol_dock: Any,
    mol_ref: Any,
    dock_to_ref: Dict[int, int],
    remap_dock: bool,
    binana_pdbindex_to_rdkit: Optional[Dict[int, int]] = None,
    binana_pdbindex_to_rdkit_ref: Optional[Dict[int, int]] = None,
    binana_pdbindex_to_rdkit_dock: Optional[Dict[int, int]] = None,
) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
    """
    Return (ligand_atoms_structure_frame, ligand_atoms_ref_otmol_frame).

    For reference structure interactions, the second list uses OTMol reference-RDKit
    indices when ``dock_to_ref`` is available. For docked structure, the first list is
    dock BINANA labels; the second is OTMol-mapped reference-ligand labels.
    """
    if binana_pdbindex_to_rdkit_ref is None and side == "ref":
        binana_pdbindex_to_rdkit_ref = binana_pdbindex_to_rdkit
    if binana_pdbindex_to_rdkit_dock is None and side == "dock":
        binana_pdbindex_to_rdkit_dock = binana_pdbindex_to_rdkit

    struct_frame: List[Dict[str, Any]] = []
    ref_otmol_frame: List[Dict[str, Any]] = []

    for ba in lig_atoms:
        if not isinstance(ba, dict):
            continue
        struct_frame.append(_normalize_binana_atom(ba, frame=side))

        if remap_dock and mol_ref is not None and dock_to_ref:
            rj = LM.binana_ligand_atom_to_ref_rdkit_idx(
                ba,
                side=side,
                mol_ref=mol_ref,
                mol_dock=mol_dock,
                dock_to_ref=dock_to_ref,
                binana_pdbindex_to_rdkit_ref=binana_pdbindex_to_rdkit_ref,
                binana_pdbindex_to_rdkit_dock=binana_pdbindex_to_rdkit_dock,
            )
            if rj is not None and not LM._is_hydrogen_binana_atom(ba):
                rec = _rdkit_atom_record(mol_ref, rj, frame="ref_ligand_via_otmol")
                if side == "dock" and mol_dock is not None:
                    dj = LM.dock_rdkit_idx_for_binana_atom(
                        mol_dock,
                        ba,
                        binana_pdbindex_to_rdkit=binana_pdbindex_to_rdkit_dock,
                    )
                    if dj is not None:
                        rec["dock_rdkit_idx"] = dj
                        rec["dock_binana_atomIndex"] = _safe_int(ba.get("atomIndex"))
                        dock_atom = mol_dock.GetAtomWithIdx(dj)
                        ref_atom = mol_ref.GetAtomWithIdx(rj)
                        rec["dock_atomic_num"] = int(dock_atom.GetAtomicNum())
                        rec["ref_atomic_num"] = int(ref_atom.GetAtomicNum())
                ref_otmol_frame.append(rec)
                continue

        if side == "ref":
            ref_otmol_frame.append(_normalize_binana_atom(ba, frame="ref_ligand"))
        else:
            ref_otmol_frame.append(
                {
                    "chain": str(ba.get("chain") or "").strip(),
                    "resID": _safe_int(ba.get("resID")),
                    "resName": str(ba.get("resName") or "").strip(),
                    "atomName": str(ba.get("atomName") or "").strip(),
                    "atomIndex": _safe_int(ba.get("atomIndex")),
                    "frame": "dock_ligand_unmapped",
                }
            )

    return struct_frame, ref_otmol_frame


def _receptor_atoms_detail(rec_atoms: List[Any]) -> List[Dict[str, Any]]:
    return [
        _normalize_binana_atom(x, frame="receptor")
        for x in (rec_atoms or [])
        if isinstance(x, dict)
    ]


def _interaction_metrics(inter: Dict[str, Any]) -> Dict[str, Any]:
    other = inter.get("other")
    if isinstance(other, dict) and "metrics" in other:
        return dict(other.get("metrics") or {})
    if "metrics" in inter:
        return dict(inter.get("metrics") or {})
    return {}


def _interaction_instance_record(
    inter: Dict[str, Any],
    *,
    side: str,
    mol_dock: Any,
    mol_ref: Any,
    dock_to_ref: Dict[int, int],
    remap_dock: bool,
    receptor_mode: str,
    binana_pdbindex_to_rdkit: Optional[Dict[int, int]] = None,
    binana_pdbindex_to_rdkit_ref: Optional[Dict[int, int]] = None,
    binana_pdbindex_to_rdkit_dock: Optional[Dict[int, int]] = None,
) -> Dict[str, Any]:
    lig = inter.get("ligandAtoms") or []
    rec = inter.get("receptorAtoms") or []

    if binana_pdbindex_to_rdkit_ref is None and side == "ref":
        binana_pdbindex_to_rdkit_ref = binana_pdbindex_to_rdkit
    if binana_pdbindex_to_rdkit_dock is None and side == "dock":
        binana_pdbindex_to_rdkit_dock = binana_pdbindex_to_rdkit

    sig = LM.interaction_signature_remapped(
        inter,
        mol_dock,
        mol_ref,
        dock_to_ref,
        remap_dock,
        receptor_mode=receptor_mode,
        binana_pdbindex_to_rdkit_ref=binana_pdbindex_to_rdkit_ref,
        binana_pdbindex_to_rdkit_dock=binana_pdbindex_to_rdkit_dock,
        side=side,
    )
    lig_struct, lig_ref = _ligand_atoms_detail(
        lig,
        side=side,
        mol_dock=mol_dock,
        mol_ref=mol_ref,
        dock_to_ref=dock_to_ref,
        remap_dock=remap_dock,
        binana_pdbindex_to_rdkit_ref=binana_pdbindex_to_rdkit_ref,
        binana_pdbindex_to_rdkit_dock=binana_pdbindex_to_rdkit_dock,
    )

    rec_d = _receptor_atoms_detail(rec)
    rec_residues = sorted(
        {(a["resID"], a["resName"]) for a in rec_d if a.get("resID") is not None}
    )

    return {
        "signature": sig,
        "side": side,
        "ligand_atoms_structure": lig_struct,
        "ligand_atoms_ref_otmol": lig_ref,
        "receptor_atoms": rec_d,
        "receptor_residues": [
            {"resID": rid, "resName": rnm} for rid, rnm in rec_residues
        ],
        "metrics": _interaction_metrics(inter),
    }


def _instances_for_key(
    js: Dict[str, Any],
    key: str,
    *,
    side: str,
    mol_dock: Any,
    mol_ref: Any,
    dock_to_ref: Dict[int, int],
    remap_dock: bool,
    receptor_mode: str,
    binana_pdbindex_to_rdkit: Optional[Dict[int, int]] = None,
    binana_pdbindex_to_rdkit_ref: Optional[Dict[int, int]] = None,
    binana_pdbindex_to_rdkit_dock: Optional[Dict[int, int]] = None,
) -> List[Dict[str, Any]]:
    out: List[Dict[str, Any]] = []
    for x in js.get(key) or []:
        if isinstance(x, dict):
            out.append(
                _interaction_instance_record(
                    x,
                    side=side,
                    mol_dock=mol_dock,
                    mol_ref=mol_ref,
                    dock_to_ref=dock_to_ref,
                    remap_dock=remap_dock,
                    receptor_mode=receptor_mode,
                    binana_pdbindex_to_rdkit=binana_pdbindex_to_rdkit,
                    binana_pdbindex_to_rdkit_ref=binana_pdbindex_to_rdkit_ref,
                    binana_pdbindex_to_rdkit_dock=binana_pdbindex_to_rdkit_dock,
                )
            )
    return out


def classify_interaction_instances(
    ref_instances: List[Dict[str, Any]],
    dock_instances: List[Dict[str, Any]],
) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]], List[Dict[str, Any]]]:
    """
    Multiset TP / FN / FP classification by canonical signature.

    TP entries pair one ref and one dock instance with the same signature.
    FN = excess ref instances; FP = excess dock instances.
    """
    ref_by_sig: Dict[str, List[Dict[str, Any]]] = defaultdict(list)
    dock_by_sig: Dict[str, List[Dict[str, Any]]] = defaultdict(list)
    for r in ref_instances:
        ref_by_sig[r["signature"]].append(r)
    for d in dock_instances:
        dock_by_sig[d["signature"]].append(d)

    tp: List[Dict[str, Any]] = []
    fn: List[Dict[str, Any]] = []
    fp: List[Dict[str, Any]] = []

    all_sigs = set(ref_by_sig) | set(dock_by_sig)
    for sig in sorted(all_sigs):
        refs = ref_by_sig.get(sig, [])
        docks = dock_by_sig.get(sig, [])
        n_tp = min(len(refs), len(docks))
        for i in range(n_tp):
            tp.append(
                {
                    "signature": sig,
                    "ref": refs[i],
                    "dock": docks[i],
                }
            )
        for r in refs[n_tp:]:
            fn.append({"signature": sig, "ref": r})
        for d in docks[n_tp:]:
            fp.append({"signature": sig, "dock": d})

    return tp, fn, fp


def _counter_for_key_pair(
    ref_js: Dict[str, Any],
    dock_js: Dict[str, Any],
    key: str,
    mol_ref: Any,
    mol_dock: Any,
    dock_to_ref: Dict[int, int],
    remap_dock: bool,
    receptor_mode: str = LM.RECEPTOR_SIG_CHAIN_ATOM,
    binana_pdbindex_to_rdkit: Optional[Dict[int, int]] = None,
    binana_pdbindex_to_rdkit_ref: Optional[Dict[int, int]] = None,
) -> Tuple[Counter, Counter]:
    ref_c: Counter = Counter()
    dock_c: Counter = Counter()
    for x in ref_js.get(key) or []:
        if isinstance(x, dict):
            ref_c[
                LM.interaction_signature_remapped(
                    x,
                    mol_dock,
                    mol_ref,
                    dock_to_ref,
                    remap_dock,
                    receptor_mode=receptor_mode,
                    binana_pdbindex_to_rdkit_ref=binana_pdbindex_to_rdkit_ref,
                    binana_pdbindex_to_rdkit_dock=binana_pdbindex_to_rdkit,
                    side="ref",
                )
            ] += 1
    for x in dock_js.get(key) or []:
        if isinstance(x, dict):
            dock_c[
                LM.interaction_signature_remapped(
                    x,
                    mol_dock,
                    mol_ref,
                    dock_to_ref,
                    remap_dock,
                    receptor_mode=receptor_mode,
                    binana_pdbindex_to_rdkit_ref=binana_pdbindex_to_rdkit_ref,
                    binana_pdbindex_to_rdkit_dock=binana_pdbindex_to_rdkit,
                    side="dock",
                )
            ] += 1
    return ref_c, dock_c


def multiset_recall(ref_c: Counter, dock_c: Counter) -> Tuple[float, int, int, int]:
    ref_total = sum(ref_c.values())
    dock_total = sum(dock_c.values())
    if ref_total == 0:
        return -1.0, 0, dock_total, 0
    tp = sum(min(ref_c[s], dock_c.get(s, 0)) for s in ref_c)
    return tp / ref_total, ref_total, dock_total, tp


def format_binana_stats_quad(frac: float, ref_n: int, dock_n: int, tp: int) -> str:
    """One CSV cell: pair_frac/ref_n/dock_n/tp_n."""
    if isinstance(frac, float) and math.isnan(frac):
        fs = "nan"
    elif isinstance(frac, float) and abs(frac - round(frac)) < 1e-9:
        fs = str(int(round(frac)))
    elif isinstance(frac, float):
        fs = f"{frac:.6g}"
    else:
        fs = str(frac)
    return f"{fs}/{ref_n}/{dock_n}/{tp}"


def _json_compact(obj: Any) -> str:
    return json.dumps(obj, ensure_ascii=True, separators=(",", ":"))


def _sig_list_from_tp(tp: List[Dict[str, Any]]) -> str:
    return "; ".join(t["signature"] for t in tp)


def _detail_columns_for_type(
    ref_js: Dict[str, Any],
    dock_js: Dict[str, Any],
    json_key: str,
    short: str,
    sfx: str,
    mol_ref: Any,
    mol_dock: Any,
    dock_to_ref: Dict[int, int],
    remap_dock: bool,
    receptor_mode: str,
    binana_pdbindex_to_rdkit: Optional[Dict[int, int]] = None,
    binana_pdbindex_to_rdkit_ref: Optional[Dict[int, int]] = None,
) -> Dict[str, Any]:
    ref_inst = _instances_for_key(
        ref_js,
        json_key,
        side="ref",
        mol_dock=mol_dock,
        mol_ref=mol_ref,
        dock_to_ref=dock_to_ref,
        remap_dock=remap_dock,
        receptor_mode=receptor_mode,
        binana_pdbindex_to_rdkit_ref=binana_pdbindex_to_rdkit_ref,
        binana_pdbindex_to_rdkit_dock=binana_pdbindex_to_rdkit,
    )
    dock_inst = _instances_for_key(
        dock_js,
        json_key,
        side="dock",
        mol_dock=mol_dock,
        mol_ref=mol_ref,
        dock_to_ref=dock_to_ref,
        remap_dock=remap_dock,
        receptor_mode=receptor_mode,
        binana_pdbindex_to_rdkit_ref=binana_pdbindex_to_rdkit_ref,
        binana_pdbindex_to_rdkit_dock=binana_pdbindex_to_rdkit,
    )
    tp, fn, fp = classify_interaction_instances(ref_inst, dock_inst)
    return {
        f"binana_{short}_tp_details{sfx}": _json_compact(tp),
        f"binana_{short}_fn_details{sfx}": _json_compact(fn),
        f"binana_{short}_fp_details{sfx}": _json_compact(fp),
        f"binana_{short}_tp_sig_list{sfx}": _sig_list_from_tp(tp),
        f"binana_{short}_fn_n{sfx}": len(fn),
        f"binana_{short}_fp_n{sfx}": len(fp),
    }


def _raw_interaction_record(inter: Dict[str, Any]) -> Dict[str, Any]:
    """Compact raw record for one BINANA interaction: atoms + raw geometry metrics.

    ``metrics`` is BINANA's own per-interaction geometry, e.g. ``{"distance": ...}``
    for halogen / hydrophobic / salt bridges, ``{"distance": ..., "angle": ...}`` for
    hydrogen bonds, and ``{"angle": ...}`` for pi-pi / T-stacking / cation-pi.
    """
    lig = [
        _normalize_binana_atom(a, frame="ligand")
        for a in (inter.get("ligandAtoms") or [])
        if isinstance(a, dict)
    ]
    rec = [
        _normalize_binana_atom(a, frame="receptor")
        for a in (inter.get("receptorAtoms") or [])
        if isinstance(a, dict)
    ]
    return {
        "ligandAtoms": lig,
        "receptorAtoms": rec,
        "metrics": _interaction_metrics(inter),
    }


def _raw_interactions_columns(
    ref_js: Dict[str, Any],
    dock_js: Dict[str, Any],
    json_key: str,
    short: str,
) -> Dict[str, Any]:
    """Per-interaction raw BINANA data (distance/angle/…) for ref and dock poses.

    Emitted once per interaction type (independent of receptor signature mode), so
    every detected interaction's raw geometry is recorded exactly once.
    """
    ref_raw = [
        _raw_interaction_record(x) for x in (ref_js.get(json_key) or []) if isinstance(x, dict)
    ]
    dock_raw = [
        _raw_interaction_record(x) for x in (dock_js.get(json_key) or []) if isinstance(x, dict)
    ]
    return {
        f"binana_{short}_ref_raw": _json_compact(ref_raw),
        f"binana_{short}_dock_raw": _json_compact(dock_raw),
        f"binana_{short}_ref_raw_n": len(ref_raw),
        f"binana_{short}_dock_raw_n": len(dock_raw),
    }


def compute_binana_recall_row(
    ligand_ref_pdb: str,
    protein_ref_pdb: str,
    ligand_dock_pdb: str,
    protein_dock_pdb: str,
    binana_python: Path,
    mol_ref=None,
    mol_dock=None,
    otmol_assignment=None,
    otmol_BCI: Any = None,
    otmol_alpha: Any = None,
    otmol_rmsd: Any = None,
    otmol_align_attempts: int = -1,
    otmol_bci_zero_ok: Optional[bool] = None,
    get_all_interactions_kwargs: Optional[Dict[str, Any]] = None,
    dock_to_ref: Optional[Dict[int, int]] = None,
    *,
    obabel_ph: Optional[float] = DEFAULT_BINANA_OBABEL_PH,
    strip_h: bool = DEFAULT_BINANA_STRIP_H,
    include_detail_columns: bool = True,
) -> Dict[str, Any]:
    """Return flat dict of fractions, counts, and (ver5) TP/FP/FN detail JSON columns."""
    gkwargs = get_all_interactions_kwargs or {}
    ref_js = run_binana_collect(
        ligand_ref_pdb,
        protein_ref_pdb,
        binana_python,
        get_all_interactions_kwargs=gkwargs,
        obabel_ph=obabel_ph,
        strip_h=strip_h,
    )
    dock_js = run_binana_collect(
        ligand_dock_pdb,
        protein_dock_pdb,
        binana_python,
        get_all_interactions_kwargs=gkwargs,
        obabel_ph=obabel_ph,
        strip_h=strip_h,
    )

    # Reuse the single global ligand OTMol mapping if provided; only fall back to
    # re-deriving it from the raw assignment when no dock_to_ref was passed in.
    # OTMol alignment itself must happen once upstream (runner / HPC); this function
    # never calls align_ligands_otmol_bci_constrained.
    dock_to_ref = dict(dock_to_ref) if dock_to_ref else {}
    map_src = "none"
    binana_pdbindex_to_rdkit: Dict[int, int] = {}
    binana_pdbindex_to_rdkit_ref: Dict[int, int] = {}
    ref_binana_atoms: Dict[int, Dict[str, Any]] = {}
    dock_binana_atoms: Dict[int, Dict[str, Any]] = {}
    otmol_diag: Dict[str, int] = {}
    if mol_ref is not None and mol_dock is not None:
        if dock_to_ref:
            map_src = "otmol"
            otmol_diag = LM.otmol_mapping_diagnostics(mol_ref, mol_dock, otmol_assignment)
        elif otmol_assignment is not None:
            dock_to_ref = LM.dock_to_ref_from_otmol_assignment(
                mol_ref, mol_dock, otmol_assignment
            )
            map_src = "otmol" if dock_to_ref else "none"
            otmol_diag = LM.otmol_mapping_diagnostics(mol_ref, mol_dock, otmol_assignment)
        else:
            raise ValueError(
                "BINANA recall requires ligand dock_to_ref or an OTMol assignment."
            )
        try:
            dock_binana_atoms = load_binana_ligand_atom_table(
                ligand_dock_pdb, binana_python, obabel_ph=obabel_ph, strip_h=strip_h
            )
            binana_pdbindex_to_rdkit = LM.build_binana_pdbindex_to_rdkit_map(
                dock_binana_atoms, mol_dock
            )
            ref_binana_atoms = load_binana_ligand_atom_table(
                ligand_ref_pdb, binana_python, obabel_ph=obabel_ph, strip_h=strip_h
            )
            binana_pdbindex_to_rdkit_ref = LM.build_binana_pdbindex_to_rdkit_map(
                ref_binana_atoms, mol_ref
            )
        except Exception as e:
            raise RuntimeError(f"BINANA ligand atom table / coordinate map failed: {e}") from e
    remap_dock = bool(dock_to_ref)  # OTMol ligand map active for BOTH ref and dock BINANA sides

    out: Dict[str, Any] = {
        "binana_h_mode": "strip_all" if strip_h else "as_is",
        "binana_ligand_map_src": map_src,
        "binana_ligand_map_n": len(dock_to_ref),
        "binana_ligand_binana_to_rdkit_n": len(binana_pdbindex_to_rdkit),
        "binana_ligand_binana_to_rdkit_ref_n": len(binana_pdbindex_to_rdkit_ref),
        "binana_otmol_BCI": otmol_BCI,
        "binana_otmol_alpha": otmol_alpha,
        "binana_otmol_rmsd": otmol_rmsd,
        "binana_otmol_align_attempts": otmol_align_attempts,
        "binana_otmol_bci_zero_ok": otmol_bci_zero_ok,
        "binana_otmol_atomicnum_mismatch": otmol_diag.get("atomicnum_mismatch", -1),
    }
    for json_key, short in BINANA_INTERACTION_KEYS:
        if include_detail_columns:
            out.update(_raw_interactions_columns(ref_js, dock_js, json_key, short))
        for rmode, sfx in RECEPTOR_SIGNATURE_MODES:
            ref_c, dock_c = _counter_for_key_pair(
                ref_js,
                dock_js,
                json_key,
                mol_ref,
                mol_dock,
                dock_to_ref,
                remap_dock,
                receptor_mode=rmode,
                binana_pdbindex_to_rdkit=binana_pdbindex_to_rdkit,
                binana_pdbindex_to_rdkit_ref=binana_pdbindex_to_rdkit_ref,
            )
            frac, ref_n, dock_n, tp = multiset_recall(ref_c, dock_c)
            out[f"binana_{short}_pair_frac{sfx}"] = frac
            out[f"binana_{short}_ref_n{sfx}"] = ref_n
            out[f"binana_{short}_dock_n{sfx}"] = dock_n
            out[f"binana_{short}_tp_n{sfx}"] = tp
            out[f"binana_{short}_stats{sfx}"] = format_binana_stats_quad(
                frac, ref_n, dock_n, tp
            )
            if include_detail_columns:
                out.update(
                    _detail_columns_for_type(
                        ref_js,
                        dock_js,
                        json_key,
                        short,
                        sfx,
                        mol_ref,
                        mol_dock,
                        dock_to_ref,
                        remap_dock,
                        rmode,
                        binana_pdbindex_to_rdkit,
                        binana_pdbindex_to_rdkit_ref,
                    )
                )

    out["binana_ref_json_keys_ok"] = json.dumps(
        sorted(k for k, _ in BINANA_INTERACTION_KEYS if k in ref_js), ensure_ascii=True
    )
    out.update(
        LM.pliq_mapping_csv_columns(
            mol_ref=mol_ref,
            mol_dock=mol_dock,
            dock_to_ref=dock_to_ref,
            bci=otmol_BCI,
            bci_zero_ok=otmol_bci_zero_ok,
            ref_binana_atoms=ref_binana_atoms or None,
            dock_binana_atoms=dock_binana_atoms or None,
            binana_to_rdkit_ref=binana_pdbindex_to_rdkit_ref or None,
            binana_to_rdkit_dock=binana_pdbindex_to_rdkit or None,
        )
    )
    return out


def empty_binana_row(*, include_detail_columns: bool = True) -> Dict[str, Any]:
    out: Dict[str, Any] = {
        "binana_ligand_map_src": "",
        "binana_ligand_map_n": -1,
        "binana_ligand_binana_to_rdkit_n": -1,
        "binana_ligand_binana_to_rdkit_ref_n": -1,
        "binana_otmol_BCI": float("nan"),
        "binana_otmol_alpha": float("nan"),
        "binana_otmol_rmsd": float("nan"),
        "binana_otmol_align_attempts": -1,
        "binana_otmol_bci_zero_ok": "",
        "binana_otmol_atomicnum_mismatch": -1,
        "pliq_mapping_trusted": "",
        "otmol_dock_to_ref_complete_ok": "",
        "binana_rdkit_ref_map_ok": "",
        "binana_rdkit_dock_map_ok": "",
        "binana_rdkit_ref_map_stats": "",
        "binana_rdkit_dock_map_stats": "",
        "pliq_mapping_warnings": "",
    }
    for _, short in BINANA_INTERACTION_KEYS:
        if include_detail_columns:
            out[f"binana_{short}_ref_raw"] = "[]"
            out[f"binana_{short}_dock_raw"] = "[]"
            out[f"binana_{short}_ref_raw_n"] = -1
            out[f"binana_{short}_dock_raw_n"] = -1
        for _, sfx in RECEPTOR_SIGNATURE_MODES:
            out[f"binana_{short}_pair_frac{sfx}"] = float("nan")
            out[f"binana_{short}_ref_n{sfx}"] = -1
            out[f"binana_{short}_dock_n{sfx}"] = -1
            out[f"binana_{short}_tp_n{sfx}"] = -1
            out[f"binana_{short}_stats{sfx}"] = "nan/-1/-1/-1"
            if include_detail_columns:
                out[f"binana_{short}_tp_details{sfx}"] = "[]"
                out[f"binana_{short}_fn_details{sfx}"] = "[]"
                out[f"binana_{short}_fp_details{sfx}"] = "[]"
                out[f"binana_{short}_tp_sig_list{sfx}"] = ""
                out[f"binana_{short}_fn_n{sfx}"] = -1
                out[f"binana_{short}_fp_n{sfx}"] = -1
    out["binana_ref_json_keys_ok"] = "[]"
    return out


def format_interaction_summary_human(
    tp: List[Dict[str, Any]],
    fn: List[Dict[str, Any]],
    fp: List[Dict[str, Any]],
    *,
    interaction_label: str,
) -> str:
    """Readable text block for scripts / notebooks."""
    lines = [f"=== {interaction_label} ===", f"TP={len(tp)}  FN(missed ref)={len(fn)}  FP(false dock)={len(fp)}"]

    def _rec_atoms(atoms: List[Dict[str, Any]]) -> str:
        parts = []
        for a in atoms:
            parts.append(f"{a.get('resID')}:{a.get('resName')}:{a.get('atomName')}")
        return ", ".join(parts) if parts else "(none)"

    def _lig_atoms(atoms: List[Dict[str, Any]]) -> str:
        parts = []
        for a in atoms:
            parts.append(
                f"{a.get('resID')}:{a.get('resName')}:{a.get('atomName')}"
                f"[{a.get('frame', '?')}]"
            )
        return ", ".join(parts) if parts else "(none)"

    if tp:
        lines.append("\n-- TRUE POSITIVE --")
        for i, item in enumerate(tp, 1):
            ref = item["ref"]
            dock = item["dock"]
            lines.append(f"  [{i}] sig={item['signature'][:120]}...")
            lines.append(f"      ref  rec: {_rec_atoms(ref['receptor_atoms'])}")
            lines.append(f"      ref  lig: {_lig_atoms(ref['ligand_atoms_structure'])}")
            lines.append(f"      dock rec: {_rec_atoms(dock['receptor_atoms'])}")
            lines.append(
                f"      dock lig(struct): {_lig_atoms(dock['ligand_atoms_structure'])}"
            )
            lines.append(
                f"      dock lig(otmol→ref): {_lig_atoms(dock['ligand_atoms_ref_otmol'])}"
            )

    if fn:
        lines.append("\n-- FALSE NEGATIVE (ref only, missed) --")
        for i, item in enumerate(fn[:20], 1):
            ref = item["ref"]
            lines.append(f"  [{i}] sig={item['signature'][:120]}...")
            lines.append(f"      ref rec: {_rec_atoms(ref['receptor_atoms'])}")
            lines.append(f"      ref lig: {_lig_atoms(ref['ligand_atoms_structure'])}")
        if len(fn) > 20:
            lines.append(f"  ... +{len(fn) - 20} more FN")

    if fp:
        lines.append("\n-- FALSE POSITIVE (dock only) --")
        for i, item in enumerate(fp[:20], 1):
            dock = item["dock"]
            lines.append(f"  [{i}] sig={item['signature'][:120]}...")
            lines.append(f"      dock rec: {_rec_atoms(dock['receptor_atoms'])}")
            lines.append(
                f"      dock lig(struct): {_lig_atoms(dock['ligand_atoms_structure'])}"
            )
            lines.append(
                f"      dock lig(otmol→ref): {_lig_atoms(dock['ligand_atoms_ref_otmol'])}"
            )
        if len(fp) > 20:
            lines.append(f"  ... +{len(fp) - 20} more FP")

    return "\n".join(lines)
