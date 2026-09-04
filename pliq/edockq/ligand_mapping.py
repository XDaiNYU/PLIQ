# -*- coding: utf-8 -*-
"""Map docked-ligand RDKit atom indices to reference-ligand indices (OTMol-only in ver5)."""
from __future__ import annotations

import json
from typing import Any, Dict, List, Optional, Set, Tuple

from rdkit import Chem


def _chain_from_ri(ri) -> str:
    if ri is None:
        return "_"
    c = getattr(ri, "GetChainId", lambda: "")()
    if c is None or (isinstance(c, str) and c.strip() == ""):
        return "_"
    return str(c).strip()


def rdkit_pdb_tuple(atom: Chem.Atom) -> Tuple[str, int, str, str]:
    """(chain, resID, resName, atomName) for an RDKit atom with PDB residue info."""
    ri = atom.GetPDBResidueInfo()
    if not ri:
        return ("_", -1, "?", atom.GetSymbol() or "?")
    ch = _chain_from_ri(ri)
    rnm = str(ri.GetResidueName() or "").strip()
    try:
        rid = int(ri.GetResidueNumber())
    except (TypeError, ValueError):
        rid = -1
    anm = str(ri.GetName().strip() if hasattr(ri, "GetName") else "").strip()
    if not anm:
        anm = atom.GetSymbol() or "?"
    return (ch, rid, rnm, anm)


def _normalize_resname(rnm: str) -> str:
    """Uppercase residue name for case-insensitive amino-acid matching (6_3 receptor)."""
    return str(rnm or "").strip().upper()


# Receptor-side signature modes for BINANA recall (see binana_recall.RECEPTOR_SIGNATURE_MODES).
RECEPTOR_SIG_CHAIN_ATOM = "chain_atom"  # chain + resID + resName + atomName (legacy CSV columns)
RECEPTOR_SIG_RID_ATOM_IDX = "rid_atom_idx"  # 6_3: ligand ref-RDKit idx; receptor amino identity (no protein idx)
RECEPTOR_SIG_RID_ATOM = "rid_atom"  # resID + resName + atomName
RECEPTOR_SIG_RID_RES = "rid_res"  # resID + resName (unique residues in this interaction)


def _rec_rid_rnm_anm(ba: Dict[str, Any]) -> Tuple[int, str, str]:
    rid = ba.get("resID")
    try:
        rid_i = int(rid)
    except (TypeError, ValueError):
        rid_i = -1
    rnm = str(ba.get("resName", "") or "").strip()
    anm = str(ba.get("atomName", "") or "").strip()
    return rid_i, rnm, anm


def _rec_chain(ba: Dict[str, Any]) -> str:
    return str(ba.get("chain", "") or "").strip()


def binana_atom_tuple(ba: Dict[str, Any]) -> Tuple[str, int, str, str]:
    ch = str(ba.get("chain", "") or "").strip()
    rid = ba.get("resID")
    try:
        rid_i = int(rid)
    except (TypeError, ValueError):
        rid_i = -1
    rnm = str(ba.get("resName", "") or "").strip()
    anm = str(ba.get("atomName", "") or "").strip()
    return (ch, rid_i, rnm, anm)


def binana_receptor_amino_tuple(ba: Dict[str, Any]) -> Tuple[str, int, str, str]:
    """Receptor atom identity: chain, resID, resName (upper), atomName — no atomIndex."""
    ch, rid, rnm, anm = binana_atom_tuple(ba)
    return (ch, rid, _normalize_resname(rnm), anm)


def receptor_signature_part(rec: List[Any], mode: str) -> str:
    """Build canonical receptor substring for interaction signature (no L[...] wrapper).

    chain is included in every mode so that multi-chain receptors never collide on a
    shared resID. Mode ``rid_atom_idx`` (6_3) uses case-normalized amino-acid identity
    (chain, resID, resName, atomName) and never includes protein atomIndex.
    """
    items = [x for x in (rec or []) if isinstance(x, dict)]
    if mode == RECEPTOR_SIG_CHAIN_ATOM:
        tups = sorted({binana_atom_tuple(x) for x in items})
        return "|".join("%s:%d:%s:%s" % t for t in tups)
    if mode == RECEPTOR_SIG_RID_ATOM_IDX:
        tups = sorted({binana_receptor_amino_tuple(x) for x in items})
        return "|".join("%s:%d:%s:%s" % t for t in tups)
    if mode == RECEPTOR_SIG_RID_ATOM:
        tups = sorted({((_rec_chain(x),) + _rec_rid_rnm_anm(x)) for x in items})
        return "|".join("%s:%d:%s:%s" % t for t in tups)
    if mode == RECEPTOR_SIG_RID_RES:
        tups = sorted({((_rec_chain(x),) + _rec_rid_rnm_anm(x)[:2]) for x in items})
        return "|".join("%s:%d:%s" % t for t in tups)
    raise ValueError(f"unknown receptor signature mode: {mode!r}")


_ELEMENT_SYMBOL_TO_Z: Dict[str, int] = {
    "H": 1,
    "C": 6,
    "N": 7,
    "O": 8,
    "F": 9,
    "P": 15,
    "S": 16,
    "CL": 17,
    "BR": 35,
    "I": 53,
}


def _is_hydrogen_binana_atom(ba: Dict[str, Any]) -> bool:
    """True if a BINANA ligand atom record is a hydrogen.

    RDKit ligand mols are loaded with removeHs=True, so hydrogens have no RDKit/OTMol
    counterpart. BINANA interactions (esp. H-bonds/salt bridges) carry the donor H in
    their ligand-atom set; we key the ligand signature on the heavy donor/acceptor only.
    """
    el = str(ba.get("element", "") or "").strip().upper()
    if el:
        return el == "H"
    return str(ba.get("atomName", "") or "").strip().upper().startswith("H")


def _safe_pdb_index(x: Any) -> Optional[int]:
    if x is None:
        return None
    try:
        return int(x)
    except (TypeError, ValueError):
        return None


def build_binana_pdbindex_to_rdkit_map(
    binana_atoms: Dict[int, Dict[str, Any]],
    mol_dock: Chem.Mol,
    *,
    max_dist: float = 0.1,
) -> Dict[int, int]:
    """
  Map BINANA ligand ``atomIndex`` (PDB serial in PDBQT) → RDKit dock index.

  Uses 3D coordinates + matching element only (never assumes atomIndex == RDKit idx).
  Greedy one-to-one assignment by ascending distance. Same-atom pairs are ~0 Å;
  ``max_dist`` is an upper bound (0.1 Å) to reject spurious cross-atom matches.
  """
    if mol_dock is None or mol_dock.GetNumConformers() <= 0:
        return {}
    conf = mol_dock.GetConformer()
    candidates: List[Tuple[float, int, int]] = []
    for pdb_idx, meta in binana_atoms.items():
        z_bin = _ELEMENT_SYMBOL_TO_Z.get(str(meta.get("element") or "").upper())
        if z_bin is None:
            continue
        bx, by, bz = float(meta["x"]), float(meta["y"]), float(meta["z"])
        for j in range(mol_dock.GetNumAtoms()):
            if mol_dock.GetAtomWithIdx(j).GetAtomicNum() != z_bin:
                continue
            pos = conf.GetAtomPosition(j)
            d = ((pos.x - bx) ** 2 + (pos.y - by) ** 2 + (pos.z - bz) ** 2) ** 0.5
            if d <= max_dist:
                candidates.append((d, int(pdb_idx), j))
    candidates.sort(key=lambda t: t[0])
    out: Dict[int, int] = {}
    used_binana: set = set()
    used_rdkit: set = set()
    for d, pdb_idx, j in candidates:
        if pdb_idx in used_binana or j in used_rdkit:
            continue
        out[pdb_idx] = j
        used_binana.add(pdb_idx)
        used_rdkit.add(j)
    return out


DEFAULT_BINANA_RDKIT_MAX_DIST = 0.1


class PliqMappingError(ValueError):
    """Raised when OTMol or BINANA→RDKit correspondence is incomplete or invalid."""


def _binana_heavy_pdb_indices(binana_atoms: Dict[int, Dict[str, Any]]) -> Set[int]:
    out: Set[int] = set()
    for pdb_idx, meta in binana_atoms.items():
        sym = str(meta.get("element") or "").upper()
        if sym in ("", "H", "D"):
            continue
        z = _ELEMENT_SYMBOL_TO_Z.get(sym)
        if z is None or z == 1:
            continue
        out.add(int(pdb_idx))
    return out


def _rdkit_heavy_indices(mol: Chem.Mol) -> Set[int]:
    return {
        j
        for j in range(mol.GetNumAtoms())
        if mol.GetAtomWithIdx(j).GetAtomicNum() != 1
    }


def assess_binana_rdkit_map(
    binana_atoms: Dict[int, Dict[str, Any]],
    mol: Chem.Mol,
    binana_to_rdkit: Dict[int, int],
    *,
    side: str,
    max_dist: float = DEFAULT_BINANA_RDKIT_MAX_DIST,
) -> Dict[str, Any]:
    """Return BINANA→RDKit mapping quality (100% heavy-atom coverage required for ok=True)."""
    heavy_binana = _binana_heavy_pdb_indices(binana_atoms)
    heavy_rdkit = _rdkit_heavy_indices(mol)
    n_heavy = len(heavy_binana)
    n_rdkit = len(heavy_rdkit)
    mapped_keys = heavy_binana & set(binana_to_rdkit)
    n_mapped = len(mapped_keys)
    issues: List[str] = []

    if not heavy_binana:
        issues.append(f"{side}: BINANA table has no heavy atoms")
    if mol is None or n_rdkit == 0:
        issues.append(f"{side}: RDKit mol has no heavy atoms")

    unmapped_binana = sorted(heavy_binana - set(binana_to_rdkit))
    if unmapped_binana:
        issues.append(
            f"{side}: unmapped BINANA PDB idx {unmapped_binana[:8]}"
            f"{'...' if len(unmapped_binana) > 8 else ''}"
        )

    mapped_rdkit = {binana_to_rdkit[p] for p in mapped_keys}
    if len(mapped_rdkit) != n_mapped:
        issues.append(f"{side}: BINANA→RDKit map not one-to-one")

    missing_rdkit = sorted(heavy_rdkit - mapped_rdkit)
    if missing_rdkit:
        issues.append(
            f"{side}: RDKit idx not covered {missing_rdkit[:8]}"
            f"{'...' if len(missing_rdkit) > 8 else ''}"
        )

    if mol is not None and mol.GetNumConformers() > 0 and mapped_keys:
        conf = mol.GetConformer()
        for pdb_idx in sorted(mapped_keys):
            rj = binana_to_rdkit[pdb_idx]
            meta = binana_atoms[pdb_idx]
            pos = conf.GetAtomPosition(rj)
            d = (
                (pos.x - float(meta["x"])) ** 2
                + (pos.y - float(meta["y"])) ** 2
                + (pos.z - float(meta["z"])) ** 2
            ) ** 0.5
            if d > max_dist:
                issues.append(f"{side}: PDB#{pdb_idx}→RDKit#{rj} d={d:.4f}Å>{max_dist}Å")

    dist_bad = any("d=" in x for x in issues)
    ok = (
        n_heavy > 0
        and n_mapped == n_heavy
        and not missing_rdkit
        and len(mapped_rdkit) == n_mapped
        and not dist_bad
    )
    frac = (n_mapped / n_heavy) if n_heavy else 0.0
    return {
        "ok": ok,
        "side": side,
        "n_heavy_binana": n_heavy,
        "n_heavy_rdkit": n_rdkit,
        "n_mapped": n_mapped,
        "frac": frac,
        "stats": f"{n_mapped}/{n_heavy}" if n_heavy else "0/0",
        "message": "; ".join(issues),
    }


def validate_binana_rdkit_map_complete(
    binana_atoms: Dict[int, Dict[str, Any]],
    mol: Chem.Mol,
    binana_to_rdkit: Dict[int, int],
    *,
    side: str,
    max_dist: float = DEFAULT_BINANA_RDKIT_MAX_DIST,
) -> None:
    """Strict mode: raise if BINANA→RDKit heavy-atom map is incomplete."""
    st = assess_binana_rdkit_map(
        binana_atoms, mol, binana_to_rdkit, side=side, max_dist=max_dist
    )
    if not st["ok"]:
        raise PliqMappingError(st["message"] or f"BINANA→RDKit map incomplete on {side} ligand")


def build_and_validate_binana_pdbindex_to_rdkit_map(
    binana_atoms: Dict[int, Dict[str, Any]],
    mol: Chem.Mol,
    *,
    side: str,
    max_dist: float = DEFAULT_BINANA_RDKIT_MAX_DIST,
) -> Dict[int, int]:
    """Build BINANA PDB-index → RDKit map and require 100% heavy-atom coverage."""
    m = build_binana_pdbindex_to_rdkit_map(binana_atoms, mol, max_dist=max_dist)
    validate_binana_rdkit_map_complete(
        binana_atoms, mol, m, side=side, max_dist=max_dist
    )
    return m


def check_otmol_bci_zero(
    bci: Any,
    *,
    bci_zero_ok: Optional[bool] = None,
) -> Tuple[bool, str]:
    """Return whether OTMol BCI is exactly zero."""
    if bci_zero_ok is False:
        return False, f"OTMol BCI!=0 (bci_zero_ok=False, BCI={bci})"
    if bci is None:
        return False, "OTMol BCI missing"
    try:
        bci_f = float(bci)
    except (TypeError, ValueError):
        return False, f"OTMol BCI not numeric: {bci!r}"
    if abs(bci_f) >= 1e-15:
        return False, f"OTMol BCI={bci} (must be 0)"
    return True, ""


def check_otmol_dock_to_ref_complete(
    mol_ref: Chem.Mol,
    mol_dock: Chem.Mol,
    dock_to_ref: Dict[int, int],
) -> Tuple[bool, str]:
    """Return whether dock_to_ref is a full bijective element-consistent OTMol map."""
    if not dock_to_ref:
        return False, "dock_to_ref empty"
    try:
        verify_dock_to_ref_same_atom(mol_ref, mol_dock, dock_to_ref)
    except ValueError as e:
        return False, str(e)
    n_ref, n_dock = mol_ref.GetNumAtoms(), mol_dock.GetNumAtoms()
    if len(dock_to_ref) != n_dock:
        return False, f"dock_to_ref incomplete: {len(dock_to_ref)}/{n_dock} dock atoms"
    ref_vals = set(dock_to_ref.values())
    if len(ref_vals) != n_ref:
        return False, f"dock_to_ref incomplete: {len(ref_vals)}/{n_ref} ref atoms"
    return True, ""


def require_otmol_bci_zero(
    bci: Any,
    *,
    bci_zero_ok: Optional[bool] = None,
) -> None:
    """Strict mode: raise if OTMol BCI is not zero."""
    ok, msg = check_otmol_bci_zero(bci, bci_zero_ok=bci_zero_ok)
    if not ok:
        raise PliqMappingError(msg)


def require_otmol_dock_to_ref_complete(
    mol_ref: Chem.Mol,
    mol_dock: Chem.Mol,
    dock_to_ref: Dict[int, int],
) -> None:
    """Strict mode: raise if dock_to_ref is incomplete."""
    ok, msg = check_otmol_dock_to_ref_complete(mol_ref, mol_dock, dock_to_ref)
    if not ok:
        raise PliqMappingError(msg)


def pliq_mapping_csv_columns(
    *,
    mol_ref: Optional[Chem.Mol] = None,
    mol_dock: Optional[Chem.Mol] = None,
    dock_to_ref: Optional[Dict[int, int]] = None,
    bci: Any = None,
    bci_zero_ok: Optional[bool] = None,
    ref_binana_atoms: Optional[Dict[int, Dict[str, Any]]] = None,
    dock_binana_atoms: Optional[Dict[int, Dict[str, Any]]] = None,
    binana_to_rdkit_ref: Optional[Dict[int, int]] = None,
    binana_to_rdkit_dock: Optional[Dict[int, int]] = None,
) -> Dict[str, Any]:
    """
    Flat CSV fields flagging ligand-mapping quality.

    ``pliq_mapping_trusted`` is True only when OTMol BCI==0, dock_to_ref is complete,
    and BINANA→RDKit coordinate maps are 100% on both ref and dock ligands.
    """
    bci_ok, bci_msg = check_otmol_bci_zero(bci, bci_zero_ok=bci_zero_ok)
    d2r_ok, d2r_msg = (False, "mol_ref/mol_dock missing")
    if mol_ref is not None and mol_dock is not None and dock_to_ref is not None:
        d2r_ok, d2r_msg = check_otmol_dock_to_ref_complete(mol_ref, mol_dock, dock_to_ref)

    ref_st = {"ok": False, "stats": "", "message": "BINANA ref map not evaluated"}
    dock_st = {"ok": False, "stats": "", "message": "BINANA dock map not evaluated"}
    if (
        mol_ref is not None
        and ref_binana_atoms is not None
        and binana_to_rdkit_ref is not None
    ):
        ref_st = assess_binana_rdkit_map(
            ref_binana_atoms, mol_ref, binana_to_rdkit_ref, side="ref"
        )
    if (
        mol_dock is not None
        and dock_binana_atoms is not None
        and binana_to_rdkit_dock is not None
    ):
        dock_st = assess_binana_rdkit_map(
            dock_binana_atoms, mol_dock, binana_to_rdkit_dock, side="dock"
        )

    binana_evaluated = (
        ref_binana_atoms is not None
        and dock_binana_atoms is not None
        and binana_to_rdkit_ref is not None
        and binana_to_rdkit_dock is not None
    )
    trusted = bci_ok and d2r_ok
    if binana_evaluated:
        trusted = trusted and bool(ref_st["ok"]) and bool(dock_st["ok"])

    warnings = "; ".join(
        x
        for x in (
            [bci_msg, d2r_msg]
            + (
                [ref_st.get("message", ""), dock_st.get("message", "")]
                if binana_evaluated
                else []
            )
        )
        if x
    )
    return {
        "pliq_mapping_trusted": trusted,
        "otmol_dock_to_ref_complete_ok": d2r_ok,
        "binana_rdkit_ref_map_ok": bool(ref_st["ok"]) if binana_evaluated else "",
        "binana_rdkit_dock_map_ok": bool(dock_st["ok"]) if binana_evaluated else "",
        "binana_rdkit_ref_map_stats": ref_st.get("stats", "") if binana_evaluated else "",
        "binana_rdkit_dock_map_stats": dock_st.get("stats", "") if binana_evaluated else "",
        "pliq_mapping_warnings": warnings[:500],
    }


def _resolve_rdkit_idx_for_binana_atom(
    mol: Chem.Mol,
    ba: Dict[str, Any],
    binana_pdbindex_to_rdkit: Optional[Dict[int, int]] = None,
) -> Optional[int]:
    """Map one BINANA ligand atom to a unique RDKit index on ``mol``.

    Coordinate-derived BINANA PDB-index map takes priority (physical identity).
    PDB tuple (chain, resID, resName, atomName) is used only when it matches exactly
    one atom. Ambiguous tuple matches (e.g. multiple ``UNL`` oxygens named ``O``) never
    fall back to ``hits[0]``.
    """
    pdb_idx = _safe_pdb_index(ba.get("atomIndex"))
    if binana_pdbindex_to_rdkit is not None and pdb_idx is not None:
        j = binana_pdbindex_to_rdkit.get(pdb_idx)
        if j is not None:
            return j

    target = binana_atom_tuple(ba)
    hits: List[int] = []
    for atom in mol.GetAtoms():
        if rdkit_pdb_tuple(atom) == target:
            hits.append(atom.GetIdx())
    if len(hits) == 1:
        return hits[0]
    return None


def dock_rdkit_idx_for_binana_atom(
    mol_dock: Chem.Mol,
    ba: Dict[str, Any],
    binana_pdbindex_to_rdkit: Optional[Dict[int, int]] = None,
) -> Optional[int]:
    """
    Match BINANA ligand atom record to RDKit atom index on mol_dock.

    Order: coordinate-derived BINANA PDB-index map → unique PDB tuple on RDKit mol.
    Does **not** treat BINANA ``atomIndex`` as RDKit ``GetIdx()`` (different orderings).
    """
    return _resolve_rdkit_idx_for_binana_atom(
        mol_dock, ba, binana_pdbindex_to_rdkit=binana_pdbindex_to_rdkit
    )


def ref_rdkit_idx_for_binana_atom(
    mol_ref: Chem.Mol,
    ba: Dict[str, Any],
    binana_pdbindex_to_rdkit: Optional[Dict[int, int]] = None,
) -> Optional[int]:
    """Match BINANA ligand atom record to RDKit atom index on mol_ref (reference pose)."""
    return _resolve_rdkit_idx_for_binana_atom(
        mol_ref, ba, binana_pdbindex_to_rdkit=binana_pdbindex_to_rdkit
    )


def binana_ligand_atom_to_ref_rdkit_idx(
    ba: Dict[str, Any],
    *,
    side: str,
    mol_ref: Chem.Mol,
    mol_dock: Optional[Chem.Mol] = None,
    dock_to_ref: Optional[Dict[int, int]] = None,
    binana_pdbindex_to_rdkit_ref: Optional[Dict[int, int]] = None,
    binana_pdbindex_to_rdkit_dock: Optional[Dict[int, int]] = None,
) -> Optional[int]:
    """Map one BINANA ligand atom to reference-ligand RDKit index (OTMol frame).

    Both ref and dock poses resolve to the same ref-RDKit index space. Dock atoms
  project through ``dock_to_ref``; ref atoms resolve on ``mol_ref`` and must lie in
    the OTMol correspondence when ``dock_to_ref`` is provided.
    """
    if side == "ref":
        rj = ref_rdkit_idx_for_binana_atom(
            mol_ref, ba, binana_pdbindex_to_rdkit=binana_pdbindex_to_rdkit_ref
        )
    elif side == "dock":
        if mol_dock is None or not dock_to_ref:
            return None
        dj = dock_rdkit_idx_for_binana_atom(
            mol_dock, ba, binana_pdbindex_to_rdkit=binana_pdbindex_to_rdkit_dock
        )
        if dj is None:
            return None
        rj = dock_to_ref.get(dj)
    else:
        return None

    if rj is None or rj < 0 or rj >= mol_ref.GetNumAtoms():
        return None
    if dock_to_ref and int(rj) not in set(dock_to_ref.values()):
        return None
    return int(rj)


def dock_to_ref_from_otmol_assignment(
    mol_ref: Chem.Mol,
    mol_dock: Chem.Mol,
    assignment: Any,
) -> Dict[int, int]:
    """
    Convert OTMol `assignment` to dock_atom_idx -> ref_atom_idx.

    Per `otmol.tools._alignment.molecule_alignment`, the transport plan P has shape
    (n_ref, n_pose) and OTMol uses ``assignment = np.argmax(P, axis=1)``, so
    ``assignment`` is 1D of length n_ref and ``assignment[i_ref]`` is the RDKit atom
    index in the **pose** (docked) molecule.

    See: otmol_package/otmol/tools/_alignment.py (molecule_alignment, molecule_alignment_rdkit).
    """
    if assignment is None:
        return {}
    n_r = mol_ref.GetNumAtoms()
    n_d = mol_dock.GetNumAtoms()
    try:
        import numpy as np

        arr = np.asarray(assignment)
    except Exception:
        return {}

    if arr.ndim == 2 and arr.shape[0] == n_r and arr.shape[1] == n_d:
        arr = np.argmax(arr, axis=1)
    if arr.ndim != 1 or int(arr.size) != n_r:
        return {}

    dock_to_ref: Dict[int, int] = {}
    for i_ref in range(n_r):
        j_pose = int(arr[i_ref])
        if not (0 <= j_pose < n_d):
            continue
        if j_pose not in dock_to_ref:
            dock_to_ref[j_pose] = i_ref
        else:
            dock_to_ref[j_pose] = min(dock_to_ref[j_pose], i_ref)
    return dock_to_ref


def ligand_dock_to_ref_map_records(
    mol_ref: Chem.Mol,
    mol_dock: Chem.Mol,
    dock_to_ref: Dict[int, int],
) -> List[Dict[str, Any]]:
    """
    Human-readable records for the single global ligand OTMol mapping (dock RDKit idx -> ref).

    Each row: dock/ref RDKit indices, PDB labels, and atomic numbers for QA.
    """
    records: List[Dict[str, Any]] = []
    for dj in sorted(dock_to_ref):
        rj = dock_to_ref[dj]
        if not (0 <= dj < mol_dock.GetNumAtoms() and 0 <= rj < mol_ref.GetNumAtoms()):
            continue
        ad = mol_dock.GetAtomWithIdx(dj)
        ar = mol_ref.GetAtomWithIdx(rj)
        d_ch, d_rid, d_rnm, d_anm = rdkit_pdb_tuple(ad)
        r_ch, r_rid, r_rnm, r_anm = rdkit_pdb_tuple(ar)
        records.append(
            {
                "dock_rdkit_idx": int(dj),
                "ref_rdkit_idx": int(rj),
                "dock_chain": d_ch,
                "dock_resID": d_rid,
                "dock_resName": d_rnm,
                "dock_atomName": d_anm,
                "ref_chain": r_ch,
                "ref_resID": r_rid,
                "ref_resName": r_rnm,
                "ref_atomName": r_anm,
                "dock_atomic_num": int(ad.GetAtomicNum()),
                "ref_atomic_num": int(ar.GetAtomicNum()),
            }
        )
    return records


def format_ligand_dock_to_ref_map_json(
    mol_ref: Chem.Mol,
    mol_dock: Chem.Mol,
    dock_to_ref: Dict[int, int],
) -> str:
    """Compact JSON list of ``ligand_dock_to_ref_map_records`` for CSV storage."""
    return json.dumps(
        ligand_dock_to_ref_map_records(mol_ref, mol_dock, dock_to_ref),
        ensure_ascii=True,
        separators=(",", ":"),
    )


def format_ligand_dock_to_ref_pairs(dock_to_ref: Dict[int, int]) -> str:
    """Minimal ``dock_rdkit_idx:ref_rdkit_idx`` pairs, sorted by dock index."""
    return ";".join(f"{dj}:{rj}" for dj, rj in sorted(dock_to_ref.items()))


PLIQ_BINANA_REMAP_SCHEMA = "1"

PLIQ_BINANA_REMAP_CSV_COLUMNS: Tuple[str, ...] = (
    "pliq_binana_remap_schema",
    "ligand_dock_to_ref_pairs",
    "ligand_dock_to_ref_map",
    "binana_pdbindex_to_rdkit_ref_pairs",
    "binana_pdbindex_to_rdkit_dock_pairs",
    "binana_pdbindex_to_rdkit_ref_map",
    "binana_pdbindex_to_rdkit_dock_map",
)


def parse_index_pair_string(pairs: str) -> Dict[int, int]:
    """Parse ``left:right;left:right`` index maps from CSV (empty -> {})."""
    out: Dict[int, int] = {}
    text = (pairs or "").strip()
    if not text:
        return out
    for chunk in text.split(";"):
        chunk = chunk.strip()
        if not chunk or ":" not in chunk:
            continue
        left, right = chunk.split(":", 1)
        try:
            out[int(left.strip())] = int(right.strip())
        except ValueError:
            continue
    return out


def binana_pdbindex_to_rdkit_map_records(
    binana_atoms: Optional[Dict[int, Dict[str, Any]]],
    mol: Optional[Chem.Mol],
    pdb_to_rdkit: Optional[Dict[int, int]],
) -> List[Dict[str, Any]]:
    """Records for CSV: BINANA PDB serial -> RDKit idx (+ labels for offline QA)."""
    if not binana_atoms or not pdb_to_rdkit:
        return []
    records: List[Dict[str, Any]] = []
    for pdb_idx in sorted(pdb_to_rdkit):
        meta = binana_atoms.get(pdb_idx) or {}
        rj = pdb_to_rdkit[pdb_idx]
        rec: Dict[str, Any] = {
            "binana_pdb_index": int(pdb_idx),
            "rdkit_idx": int(rj),
            "element": str(meta.get("element") or "").upper(),
            "chain": str(meta.get("chain") or "").strip(),
            "resID": meta.get("resID"),
            "resName": str(meta.get("resName") or "").strip(),
            "atomName": str(meta.get("atomName") or "").strip(),
        }
        if mol is not None and 0 <= rj < mol.GetNumAtoms():
            ch, rid, rnm, anm = rdkit_pdb_tuple(mol.GetAtomWithIdx(rj))
            rec["rdkit_chain"] = ch
            rec["rdkit_resID"] = rid
            rec["rdkit_resName"] = rnm
            rec["rdkit_atomName"] = anm
            rec["rdkit_atomic_num"] = int(mol.GetAtomWithIdx(rj).GetAtomicNum())
        records.append(rec)
    return records


def format_binana_pdbindex_to_rdkit_map_json(
    binana_atoms: Optional[Dict[int, Dict[str, Any]]],
    mol: Optional[Chem.Mol],
    pdb_to_rdkit: Optional[Dict[int, int]],
) -> str:
    return json.dumps(
        binana_pdbindex_to_rdkit_map_records(binana_atoms, mol, pdb_to_rdkit),
        ensure_ascii=True,
        separators=(",", ":"),
    )


def format_binana_pdbindex_to_rdkit_pairs(pdb_to_rdkit: Optional[Dict[int, int]]) -> str:
    if not pdb_to_rdkit:
        return ""
    return ";".join(f"{p}:{j}" for p, j in sorted(pdb_to_rdkit.items()))


def empty_pliq_binana_remap_csv_columns() -> Dict[str, str]:
    return {
        "pliq_binana_remap_schema": PLIQ_BINANA_REMAP_SCHEMA,
        "ligand_dock_to_ref_pairs": "",
        "ligand_dock_to_ref_map": "[]",
        "binana_pdbindex_to_rdkit_ref_pairs": "",
        "binana_pdbindex_to_rdkit_dock_pairs": "",
        "binana_pdbindex_to_rdkit_ref_map": "[]",
        "binana_pdbindex_to_rdkit_dock_map": "[]",
    }


def pliq_binana_remap_csv_columns(
    mol_ref: Optional[Chem.Mol],
    mol_dock: Optional[Chem.Mol],
    dock_to_ref: Optional[Dict[int, int]],
    *,
    ref_binana_atoms: Optional[Dict[int, Dict[str, Any]]] = None,
    dock_binana_atoms: Optional[Dict[int, Dict[str, Any]]] = None,
    binana_to_rdkit_ref: Optional[Dict[int, int]] = None,
    binana_to_rdkit_dock: Optional[Dict[int, int]] = None,
) -> Dict[str, str]:
    """
    Flat CSV fields to rebuild PLIQ BINANA ligand remapping offline from one row.

    Chain for dock-side ligand atoms in raw BINANA JSON:
      ``binana_pdb_index`` -> (``binana_pdbindex_to_rdkit_dock``) -> dock RDKit idx
      -> (``ligand_dock_to_ref_pairs``) -> ref RDKit idx.
    Ref-side ligand atoms:
      ``binana_pdb_index`` -> (``binana_pdbindex_to_rdkit_ref``) -> ref RDKit idx.
    """
    out = empty_pliq_binana_remap_csv_columns()
    d2r = dict(dock_to_ref or {})
    if mol_ref is not None and mol_dock is not None and d2r:
        out["ligand_dock_to_ref_pairs"] = format_ligand_dock_to_ref_pairs(d2r)
        out["ligand_dock_to_ref_map"] = format_ligand_dock_to_ref_map_json(
            mol_ref, mol_dock, d2r
        )
    b_ref = dict(binana_to_rdkit_ref or {})
    b_dock = dict(binana_to_rdkit_dock or {})
    if b_ref:
        out["binana_pdbindex_to_rdkit_ref_pairs"] = format_binana_pdbindex_to_rdkit_pairs(b_ref)
        out["binana_pdbindex_to_rdkit_ref_map"] = format_binana_pdbindex_to_rdkit_map_json(
            ref_binana_atoms, mol_ref, b_ref
        )
    if b_dock:
        out["binana_pdbindex_to_rdkit_dock_pairs"] = format_binana_pdbindex_to_rdkit_pairs(b_dock)
        out["binana_pdbindex_to_rdkit_dock_map"] = format_binana_pdbindex_to_rdkit_map_json(
            dock_binana_atoms, mol_dock, b_dock
        )
    return out


def parse_pliq_binana_remap_from_csv_row(row: Dict[str, Any]) -> Dict[str, Any]:
    """Load remap tables from a PLIQ stepwise CSV row (dict-like)."""
    return {
        "schema": str(row.get("pliq_binana_remap_schema") or ""),
        "dock_to_ref": parse_index_pair_string(str(row.get("ligand_dock_to_ref_pairs") or "")),
        "dock_to_ref_records": json.loads(str(row.get("ligand_dock_to_ref_map") or "[]") or "[]"),
        "binana_to_rdkit_ref": parse_index_pair_string(
            str(row.get("binana_pdbindex_to_rdkit_ref_pairs") or "")
        ),
        "binana_to_rdkit_dock": parse_index_pair_string(
            str(row.get("binana_pdbindex_to_rdkit_dock_pairs") or "")
        ),
        "binana_to_rdkit_ref_records": json.loads(
            str(row.get("binana_pdbindex_to_rdkit_ref_map") or "[]") or "[]"
        ),
        "binana_to_rdkit_dock_records": json.loads(
            str(row.get("binana_pdbindex_to_rdkit_dock_map") or "[]") or "[]"
        ),
    }


def binana_ligand_atomindex_to_ref_rdkit_idx(
    binana_pdb_index: int,
    *,
    side: str,
    binana_to_rdkit_ref: Dict[int, int],
    binana_to_rdkit_dock: Dict[int, int],
    dock_to_ref: Dict[int, int],
) -> Optional[int]:
    """
    Map one BINANA ligand ``atomIndex`` (PDB serial) to reference-ligand RDKit idx.

    ``side`` is ``"ref"`` or ``"dock"`` (same convention as BINANA recall).
    """
    side = str(side).strip().lower()
    if side == "ref":
        return binana_to_rdkit_ref.get(int(binana_pdb_index))
    if side == "dock":
        dj = binana_to_rdkit_dock.get(int(binana_pdb_index))
        if dj is None:
            return None
        return dock_to_ref.get(int(dj))
    raise ValueError(f"side must be 'ref' or 'dock', got {side!r}")


def remap_binana_ligand_atom_indices_from_csv_row(
    ligand_atoms: List[Dict[str, Any]],
    *,
    side: str,
    row: Dict[str, Any],
) -> List[int]:
    """Return sorted unique ref-RDKit indices for ligand atoms in one BINANA interaction."""
    maps = parse_pliq_binana_remap_from_csv_row(row)
    idxs: List[int] = []
    for ba in ligand_atoms or []:
        if not isinstance(ba, dict):
            continue
        pdb_idx = _safe_pdb_index(ba.get("atomIndex"))
        if pdb_idx is None:
            continue
        if _is_hydrogen_binana_atom(ba):
            continue
        rj = binana_ligand_atomindex_to_ref_rdkit_idx(
            pdb_idx,
            side=side,
            binana_to_rdkit_ref=maps["binana_to_rdkit_ref"],
            binana_to_rdkit_dock=maps["binana_to_rdkit_dock"],
            dock_to_ref=maps["dock_to_ref"],
        )
        if rj is not None:
            idxs.append(int(rj))
    return sorted(set(idxs))


def rebuild_interaction_signature_from_csv_row(
    interaction: Dict[str, Any],
    *,
    side: str,
    row: Dict[str, Any],
    receptor_mode: str = RECEPTOR_SIG_RID_ATOM,
) -> str:
    """Rebuild PLIQ BINANA signature from CSV remap columns + raw interaction JSON."""
    lig = interaction.get("ligandAtoms") or []
    rec = interaction.get("receptorAtoms") or []
    lig_part = "|".join(
        str(i) for i in remap_binana_ligand_atom_indices_from_csv_row(lig, side=side, row=row)
    )
    rec_part = receptor_signature_part(rec, receptor_mode)
    return f"L[{lig_part}]R[{rec_part}]"


def validate_offline_binana_remap_from_csv_row(
    row: Dict[str, Any],
    *,
    receptor_mode: str = RECEPTOR_SIG_RID_ATOM,
    interaction_keys: Optional[List[Tuple[str, str]]] = None,
) -> None:
    """
    Raise ``ValueError`` if ``*_ref_raw`` / ``*_dock_raw`` signatures rebuilt from CSV remap
    columns disagree with stored signatures.

    When ``*_tp_details*`` / ``*_fn_details*`` / ``*_fp_details*`` are present, also checks
    those nested instance records (``ref`` / ``dock`` sub-objects from classify step).
    """
    keys = interaction_keys or []
    if not keys:
        from . import binana_recall as _br  # local import avoids cycle at module load

        keys = _br.BINANA_INTERACTION_KEYS
    suffix = ""
    if receptor_mode == RECEPTOR_SIG_RID_ATOM_IDX:
        suffix = "_rxidx"
    elif receptor_mode == RECEPTOR_SIG_RID_ATOM:
        suffix = "_rxatm"
    elif receptor_mode == RECEPTOR_SIG_RID_RES:
        suffix = "_rxres"

    def _check_instance(
        inst: Dict[str, Any],
        *,
        side: str,
        col: str,
        stored_sig: str,
    ) -> None:
        rebuilt = rebuild_interaction_signature_from_csv_row(
            {
                "ligandAtoms": inst.get("ligand_atoms_structure")
                or inst.get("ligandAtoms")
                or [],
                "receptorAtoms": inst.get("receptor_atoms")
                or inst.get("receptorAtoms")
                or [],
            },
            side=side,
            row=row,
            receptor_mode=receptor_mode,
        )
        if rebuilt != stored_sig:
            raise ValueError(
                f"{col}: signature mismatch side={side} stored={stored_sig!r} rebuilt={rebuilt!r}"
            )

    for _, short in keys:
        for bucket in ("tp", "fn", "fp"):
            col = f"binana_{short}_{bucket}_details{suffix}"
            raw = row.get(col)
            if not raw or raw == "[]":
                continue
            for item in json.loads(str(raw)):
                sig = str(item.get("signature") or "")
                if bucket == "tp":
                    if "ref" in item:
                        _check_instance(item["ref"], side="ref", col=col, stored_sig=sig)
                    if "dock" in item:
                        _check_instance(item["dock"], side="dock", col=col, stored_sig=sig)
                elif bucket == "fn" and "ref" in item:
                    _check_instance(item["ref"], side="ref", col=col, stored_sig=sig)
                elif bucket == "fp" and "dock" in item:
                    _check_instance(item["dock"], side="dock", col=col, stored_sig=sig)


def validate_rdkit_mol_for_otmol(mol: Optional[Chem.Mol], *, name: str) -> None:
    """Raise ValueError if an RDKit Mol is unusable for OTMol alignment."""
    if mol is None:
        raise ValueError(f"{name}: mol is None")
    if not isinstance(mol, Chem.Mol):
        raise ValueError(f"{name}: expected rdkit.Chem.Mol, got {type(mol)!r}")
    if mol.GetNumAtoms() <= 0:
        raise ValueError(f"{name}: mol has no atoms")
    if mol.GetNumConformers() <= 0:
        raise ValueError(f"{name}: mol has no conformer/3D coordinates")


def _atom_name_of(mol: Chem.Mol, idx: int) -> Optional[str]:
    ri = mol.GetAtomWithIdx(idx).GetPDBResidueInfo()
    if ri is not None and ri.GetName() is not None:
        return ri.GetName().strip()
    return None


def otmol_mapping_diagnostics(
    mol_ref: Optional[Chem.Mol],
    mol_dock: Optional[Chem.Mol],
    assignment: Any,
) -> Dict[str, int]:
    """Simple mapping diagnostics for debugging/QA.

    ``atomicnum_mismatch`` counts pairs whose elements differ (must be 0 — a non-zero
    value means OTMol mapped an atom to a chemically different atom). ``name_mismatch``
    counts pairs whose PDB atom names differ; this is informational only, because ref
    and dock files can use entirely different atom-naming schemes for the same molecule
    (so name equality is NOT a valid same-atom test — element + bond-consistency is).
    """
    if mol_ref is None or mol_dock is None or assignment is None:
        return {
            "ref_atoms": int(mol_ref.GetNumAtoms()) if isinstance(mol_ref, Chem.Mol) else -1,
            "dock_atoms": int(mol_dock.GetNumAtoms()) if isinstance(mol_dock, Chem.Mol) else -1,
            "map_pairs": 0,
            "unmapped_dock": -1 if not isinstance(mol_dock, Chem.Mol) else int(mol_dock.GetNumAtoms()),
            "atomicnum_mismatch": -1,
            "name_mismatch": -1,
            "duplicate_ref_targets": -1,
        }
    d2r = dock_to_ref_from_otmol_assignment(mol_ref, mol_dock, assignment)
    mismatch = 0
    name_mismatch = 0
    for d_idx, r_idx in d2r.items():
        if not (0 <= d_idx < mol_dock.GetNumAtoms() and 0 <= r_idx < mol_ref.GetNumAtoms()):
            mismatch += 1
            continue
        a_d = mol_dock.GetAtomWithIdx(d_idx).GetAtomicNum()
        a_r = mol_ref.GetAtomWithIdx(r_idx).GetAtomicNum()
        if a_d != a_r:
            mismatch += 1
        nd, nr = _atom_name_of(mol_dock, d_idx), _atom_name_of(mol_ref, r_idx)
        if nd is not None and nr is not None and nd != nr:
            name_mismatch += 1
    ref_targets = list(d2r.values())
    duplicate_ref_targets = len(ref_targets) - len(set(ref_targets))
    return {
        "ref_atoms": int(mol_ref.GetNumAtoms()),
        "dock_atoms": int(mol_dock.GetNumAtoms()),
        "map_pairs": int(len(d2r)),
        "unmapped_dock": int(max(0, mol_dock.GetNumAtoms() - len(d2r))),
        "atomicnum_mismatch": int(mismatch),
        "name_mismatch": int(name_mismatch),
        "duplicate_ref_targets": int(duplicate_ref_targets),
    }


def verify_dock_to_ref_same_atom(
    mol_ref: Chem.Mol,
    mol_dock: Chem.Mol,
    dock_to_ref: Dict[int, int],
) -> None:
    """Guarantee that each dock ligand atom maps to the SAME atom in the ref ligand.

    "Same atom" is enforced by (1) a one-to-one (bijective) index map and (2) identical
    element on both sides for every pair. A failure here means OTMol produced a wrong
    correspondence (mapped an atom to a chemically different / reused atom), which would
    silently corrupt ligand RMSD, interface RMSD, FNAT and BINANA — so we raise instead.

    Atom-name equality is intentionally NOT required: ref and dock can name the same
    molecule completely differently; element identity plus the BCI==0 bond-consistency
    already enforced by ``align_ligands_otmol_bci_constrained`` are the valid criteria.
    """
    if not dock_to_ref:
        raise ValueError("empty dock_to_ref: no ligand atom correspondence")
    n_ref, n_dock = mol_ref.GetNumAtoms(), mol_dock.GetNumAtoms()
    ref_targets = list(dock_to_ref.values())
    dup = len(ref_targets) - len(set(ref_targets))
    if dup:
        raise ValueError(
            f"dock_to_ref is not one-to-one: {dup} ref atom(s) reused by multiple dock atoms"
        )
    bad_idx, elem_mis = 0, 0
    for d_idx, r_idx in dock_to_ref.items():
        if not (0 <= d_idx < n_dock and 0 <= r_idx < n_ref):
            bad_idx += 1
            continue
        if mol_dock.GetAtomWithIdx(d_idx).GetAtomicNum() != mol_ref.GetAtomWithIdx(r_idx).GetAtomicNum():
            elem_mis += 1
    if bad_idx:
        raise ValueError(f"dock_to_ref has {bad_idx} out-of-range index pair(s)")
    if elem_mis:
        raise ValueError(
            f"dock_to_ref maps {elem_mis} atom(s) to a different element — not the same atom"
        )


def require_otmol_importable() -> None:
    """Raise ImportError if OTMol cannot be loaded (after vendored path shim)."""
    from .otmol_runtime import ensure_local_otmol_on_path

    ensure_local_otmol_on_path()
    try:
        import otmol  # noqa: F401
    except ImportError as e:
        raise ImportError(
            "PLIQ requires OTMol for ligand RMSD and BINANA ligand mapping. "
            "Install from https://github.com/weixiaoqimath/otmol "
            "(pip install git+https://github.com/weixiaoqimath/otmol.git) "
            "or add a local ``otmol_package`` checkout to PYTHONPATH."
        ) from e


def align_ligands_otmol_bci_constrained(
    mol_ref: Chem.Mol,
    mol_dock: Chem.Mol,
    *,
    max_attempts: int = 20,
    require_bci_zero: bool = False,
) -> Tuple[Dict[int, int], Any, float, Any, Any, int, bool]:
    """
    OTMol-only ligand alignment with RMSD on original coordinates (e-DockQ parity).

    Tries up to ``max_attempts`` hyperparameter presets (alpha grid, cst_D, method;
    reflection is always disabled) and stops early when bond-consistency index (BCI)
    is exactly zero.
    If no attempt yields BCI == 0, returns the best attempt (lowest BCI) after
    exhausting the budget.

    Returns:
        (dock_to_ref, assignment, rmsd, alpha, BCI, attempts_used, bci_is_zero)
    """
    validate_rdkit_mol_for_otmol(mol_ref, name="mol_ref")
    validate_rdkit_mol_for_otmol(mol_dock, name="mol_dock")
    require_otmol_importable()
    import numpy as np
    import otmol as otm

    presets: List[Dict[str, Any]] = [
        {"alpha_list": np.arange(0.01, 1.0, 0.01).tolist(), "cst_D": 0.5, "method": "fGW", "reflection": False},
        {"alpha_list": np.arange(0.01, 1.0, 0.02).tolist(), "cst_D": 0.5, "method": "fGW", "reflection": False},
        {"alpha_list": np.linspace(0.02, 0.98, 25).tolist(), "cst_D": 0.5, "method": "fGW", "reflection": False},
        {"alpha_list": np.arange(0.01, 1.0, 0.01).tolist(), "cst_D": 0.0, "method": "fGW", "reflection": False},
        {"alpha_list": np.arange(0.01, 1.0, 0.01).tolist(), "cst_D": 0.25, "method": "fGW", "reflection": False},
        {"alpha_list": np.arange(0.01, 1.0, 0.01).tolist(), "cst_D": 0.75, "method": "fGW", "reflection": False},
        {"alpha_list": np.arange(0.01, 1.0, 0.01).tolist(), "cst_D": 1.0, "method": "fGW", "reflection": False},
        {"alpha_list": [0.5], "cst_D": 0.5, "method": "fGW", "reflection": False},
        {"alpha_list": [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9], "cst_D": 0.0, "method": "fGW", "reflection": False},
        {"alpha_list": np.arange(0.01, 1.0, 0.03).tolist(), "cst_D": 0.5, "method": "fsGW", "reflection": False},
        {"alpha_list": np.arange(0.01, 1.0, 0.01).tolist(), "cst_D": 0.5, "method": "fsGW", "reflection": False},
        {"alpha_list": np.arange(0.01, 1.0, 0.01).tolist(), "cst_D": 0.0, "method": "fsGW", "reflection": False},
        {"alpha_list": np.arange(0.01, 1.0, 0.01).tolist(), "cst_D": 1.0, "method": "fsGW", "reflection": False},
        {"alpha_list": np.arange(0.005, 1.0, 0.005).tolist(), "cst_D": 0.5, "method": "fGW", "reflection": False},
        {"alpha_list": np.linspace(0.01, 0.99, 50).tolist(), "cst_D": 0.35, "method": "fGW", "reflection": False},
        {"alpha_list": np.linspace(0.01, 0.99, 50).tolist(), "cst_D": 0.65, "method": "fGW", "reflection": False},
        {"alpha_list": np.arange(0.01, 1.0, 0.01).tolist(), "cst_D": 0.5, "method": "fsGW", "reflection": False},
        {"alpha_list": [0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5], "cst_D": 0.5, "method": "fGW", "reflection": False},
        {"alpha_list": [0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.99], "cst_D": 0.5, "method": "fGW", "reflection": False},
        {"alpha_list": np.arange(0.02, 1.0, 0.04).tolist(), "cst_D": 0.5, "method": "fGW", "reflection": False},
    ]
    presets = presets[: int(max_attempts)]

    best: Optional[Tuple[Any, float, Any, Any, Dict[int, int]]] = None
    best_bci = float("inf")
    attempts_used = 0

    for i, preset in enumerate(presets):
        attempts_used = i + 1
        try:
            result = otm.tl.molecule_alignment_rdkit(
                mol_ref,
                mol_dock,
                minimize_mismatched_edges=True,
                use_original_coords_rmsd=True,
                alpha_list=preset["alpha_list"],
                cst_D=float(preset["cst_D"]),
                method=str(preset["method"]),
                reflection=False,
            )
        except Exception:
            continue
        if not result or len(result) < 4:
            continue
        assignment, rmsd, alpha, BCI = result[0], result[1], result[2], result[3]
        if assignment is None or rmsd is None:
            continue
        try:
            rmsd_f = float(rmsd)
            bci_f = float(BCI) if BCI is not None else float("inf")
        except (TypeError, ValueError):
            continue
        d2r = dock_to_ref_from_otmol_assignment(mol_ref, mol_dock, assignment)
        if not d2r:
            continue
        # Reject any correspondence that is not a valid same-atom map (bijective +
        # element-consistent); never let a chemically wrong mapping through.
        try:
            verify_dock_to_ref_same_atom(mol_ref, mol_dock, d2r)
        except ValueError:
            continue
        if abs(bci_f) < 1e-15:
            return d2r, assignment, rmsd_f, alpha, BCI, attempts_used, True
        if bci_f < best_bci:
            best_bci = bci_f
            best = (assignment, rmsd_f, alpha, BCI, d2r)

    if best is None:
        raise PliqMappingError(
            "OTMol molecule_alignment_rdkit did not return a usable assignment for this ligand pair "
            f"(tried {attempts_used} preset(s))."
        )
    assignment, rmsd_f, alpha, BCI, d2r = best
    if require_bci_zero:
        raise PliqMappingError(
            "OTMol BCI must be exactly 0 for PLIQ; "
            f"no zero-BCI preset succeeded in {attempts_used} attempt(s) (best BCI={BCI})."
        )
    return d2r, assignment, rmsd_f, alpha, BCI, attempts_used, False


def ligand_ref_rdkit_indices_for_binana(
    lig_atoms: List[Any],
    *,
    side: str,
    mol_ref: Chem.Mol,
    mol_dock: Optional[Chem.Mol] = None,
    dock_to_ref: Optional[Dict[int, int]] = None,
    binana_pdbindex_to_rdkit: Optional[Dict[int, int]] = None,
    binana_pdbindex_to_rdkit_ref: Optional[Dict[int, int]] = None,
    binana_pdbindex_to_rdkit_dock: Optional[Dict[int, int]] = None,
) -> Optional[List[int]]:
    """Map BINANA ligand atoms to reference-ligand RDKit indices (OTMol frame).

    Both ref and dock BINANA poses resolve into ref-RDKit index space. When
    ``dock_to_ref`` is set (normal PLIQ path), every heavy atom must map to an OTMol
    partner on ``mol_ref``. Legacy single ``binana_pdbindex_to_rdkit`` applies to the
    active pose (ref map on ref side, dock map on dock side).
    """
    if binana_pdbindex_to_rdkit_ref is None and side == "ref":
        binana_pdbindex_to_rdkit_ref = binana_pdbindex_to_rdkit
    if binana_pdbindex_to_rdkit_dock is None and side == "dock":
        binana_pdbindex_to_rdkit_dock = binana_pdbindex_to_rdkit

    out: List[int] = []
    for ba in lig_atoms:
        if not isinstance(ba, dict):
            return None
        if _is_hydrogen_binana_atom(ba):
            continue  # H has no RDKit/OTMol counterpart; key on heavy atoms only
        rj = binana_ligand_atom_to_ref_rdkit_idx(
            ba,
            side=side,
            mol_ref=mol_ref,
            mol_dock=mol_dock,
            dock_to_ref=dock_to_ref,
            binana_pdbindex_to_rdkit_ref=binana_pdbindex_to_rdkit_ref,
            binana_pdbindex_to_rdkit_dock=binana_pdbindex_to_rdkit_dock,
        )
        if rj is None:
            return None
        out.append(int(rj))
    return out


def interaction_signature_remapped(
    inter: Dict[str, Any],
    mol_dock: Chem.Mol,
    mol_ref: Chem.Mol,
    dock_to_ref: Dict[int, int],
    use_remapping: bool,
    receptor_mode: str = RECEPTOR_SIG_CHAIN_ATOM,
    binana_pdbindex_to_rdkit: Optional[Dict[int, int]] = None,
    binana_pdbindex_to_rdkit_ref: Optional[Dict[int, int]] = None,
    binana_pdbindex_to_rdkit_dock: Optional[Dict[int, int]] = None,
    *,
    side: str = "ref",
) -> str:
    """BINANA recall signature; ligand atoms in OTMol reference-RDKit index space.

    Both ref and dock poses use ``dock_to_ref`` when ``use_remapping`` is true so every
    signature ``L[...]`` index is a reference-ligand RDKit idx in the OTMol correspondence.
    ``receptor_mode`` only controls how the receptor substring is built.
    """
    lig = inter.get("ligandAtoms") or []
    rec = inter.get("receptorAtoms") or []

    if binana_pdbindex_to_rdkit_ref is None and side == "ref":
        binana_pdbindex_to_rdkit_ref = binana_pdbindex_to_rdkit
    if binana_pdbindex_to_rdkit_dock is None and side == "dock":
        binana_pdbindex_to_rdkit_dock = binana_pdbindex_to_rdkit

    lig_part: Optional[str] = None
    if mol_ref is not None:
        idxs = ligand_ref_rdkit_indices_for_binana(
            lig,
            side=side,
            mol_ref=mol_ref,
            mol_dock=mol_dock if use_remapping else None,
            dock_to_ref=dock_to_ref if use_remapping else None,
            binana_pdbindex_to_rdkit_ref=binana_pdbindex_to_rdkit_ref,
            binana_pdbindex_to_rdkit_dock=binana_pdbindex_to_rdkit_dock,
        )
        if idxs is not None:
            lig_part = "|".join(str(i) for i in sorted(set(idxs)))
    if lig_part is None:
        # Fallback only when index mapping is unavailable: raw BINANA PDB serials.
        lig_part = "|".join(
            str(_safe_pdb_index(x.get("atomIndex")) or -1)
            for x in sorted(
                (a for a in lig if isinstance(a, dict)),
                key=lambda d: _safe_pdb_index(d.get("atomIndex")) or -1,
            )
        )

    rec_part = receptor_signature_part(rec, receptor_mode)
    return f"L[{lig_part}]R[{rec_part}]"
