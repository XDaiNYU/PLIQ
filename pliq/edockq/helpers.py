# -*- coding: utf-8 -*-
"""Helper functions for e-DockQ: PDB loading, interface RMSD, FNAT (residue + backbone/sidechain).

Ligand atom correspondence is OTMol-only (dock->ref atom index map); there is no MCS anywhere.
Protein atoms are paired by residue id + atom name. No global superposition is applied.
"""
from __future__ import annotations

import os
import tempfile
from typing import Dict, List, Optional, Set, Tuple

import numpy as np
from rdkit import Chem
import prody


def _select_entity_force_chain(structure, chain_id):
    """Select the non-water entity and relabel it to ``chain_id``.

    The ligand is always read as chain 'B' and the protein as chain 'A', regardless
    of the chain written in the file (docking/BINANA tools sometimes emit HETATM
    ligands as chain 'X'). We prefer atoms already on ``chain_id``; if none exist we
    take all non-water atoms, then force every selected atom's chain to ``chain_id``
    so ref and dock share canonical chain labels.
    """
    if structure is None:
        return None
    sel = structure.select(f"chain {chain_id} and not resname HOH")
    if sel is None:
        sel = structure.select("not resname HOH")
    if sel is None:
        return None
    sel.setChids([chain_id] * sel.numAtoms())
    return sel


def write_forced_chain_pdb(file_path, chain_id, out_path):
    """Write a copy of ``file_path`` with the non-water entity relabeled to ``chain_id``.

    Used to canonicalize chains (ligand→'B', protein→'A') before feeding PDBs to
    obabel/BINANA so receptor signatures use canonical chains for both ref and dock.
    Returns ``out_path`` on success, or ``None`` if nothing could be selected.
    """
    structure = prody.parsePDB(file_path)
    sel = _select_entity_force_chain(structure, chain_id)
    if sel is None:
        return None
    prody.writePDB(str(out_path), sel)
    return str(out_path)


def pdb_to_rdkit_with_prody(file_path, chain_id, temp_pdb_path=None):
    """Load a ligand PDB and convert to RDKit, forcing the ligand onto ``chain_id`` (default 'B')."""
    structure = prody.parsePDB(file_path)
    chain = _select_entity_force_chain(structure, chain_id)
    if chain is None:
        return None
    if temp_pdb_path is None:
        fd, temp_pdb_path = tempfile.mkstemp(suffix='.pdb')
        try:
            prody.writePDB(temp_pdb_path, chain)
            rdkit_mol = Chem.MolFromPDBFile(temp_pdb_path, removeHs=True)
            return rdkit_mol
        finally:
            os.close(fd)
            if os.path.exists(temp_pdb_path):
                os.remove(temp_pdb_path)
    prody.writePDB(temp_pdb_path, chain)
    rdkit_mol = Chem.MolFromPDBFile(temp_pdb_path, removeHs=True)
    return rdkit_mol


def pdb_chain_to_rdkit_mol(file_path, chain_id="A"):
    """Load a protein PDB and convert to RDKit, forcing the protein onto ``chain_id`` (default 'A')."""
    structure = prody.parsePDB(file_path)
    chain = _select_entity_force_chain(structure, chain_id or "A")
    if chain is None:
        return None
    temp_fd, temp_path = tempfile.mkstemp(suffix='.pdb')
    try:
        prody.writePDB(temp_path, chain)
        rdkit_mol = Chem.MolFromPDBFile(temp_path, removeHs=False, sanitize=False)
        if rdkit_mol is None:
            return None
        prody_atoms = list(chain)
        rdkit_atoms = list(rdkit_mol.GetAtoms())
        if len(prody_atoms) != len(rdkit_atoms):
            return None
        for prody_atom, rdkit_atom in zip(prody_atoms, rdkit_atoms):
            rdkit_atom.SetProp("residue_name", prody_atom.getResname())
            rdkit_atom.SetProp("residue_number", str(prody_atom.getResnum()))
            rdkit_atom.SetProp("atom_name", prody_atom.getName())
        return rdkit_mol
    finally:
        os.close(temp_fd)
        os.remove(temp_path)


def combine_coords_and_calculate_rmsd(sub_iprot1_coords, sub_iprot2_coords, sub_imol1_coords1, sub_imol2_coords2):
    def ensure_2d(a):
        a = np.array(a)
        return a.reshape((-1, 3)) if a.size else np.zeros((0, 3))
    c1a = ensure_2d(sub_iprot1_coords)
    c1b = ensure_2d(sub_imol1_coords1)
    c2a = ensure_2d(sub_iprot2_coords)
    c2b = ensure_2d(sub_imol2_coords2)
    combined_coords1 = np.vstack((c1a, c1b))
    combined_coords2 = np.vstack((c2a, c2b))
    if combined_coords1.shape != combined_coords2.shape or combined_coords1.shape[0] == 0:
        return None
    return float(np.sqrt(np.mean(np.sum((combined_coords1 - combined_coords2) ** 2, axis=1))))


# Protein interface backbone atom names (CA-only reproduces the legacy MCS behaviour;
# the full backbone {N, CA, C, O} becomes available once OTMol drives ligand correspondence).
INTERFACE_PROTEIN_CA = ("CA",)
INTERFACE_PROTEIN_BACKBONE = ("N", "CA", "C", "O")


def _atom_name_from_pdb_info(atom) -> Optional[str]:
    ri = atom.GetPDBResidueInfo()
    if ri is not None:
        nm = ri.GetName()
        if nm is not None and nm.strip():
            return nm.strip()
    if atom.HasProp("atom_name"):
        nm = atom.GetProp("atom_name").strip()
        if nm:
            return nm
    return None


def _residue_key_from_pdb_info(atom) -> Optional[Tuple[int, str]]:
    ri = atom.GetPDBResidueInfo()
    if ri is not None:
        return (ri.GetResidueNumber(), ri.GetResidueName().strip())
    if atom.HasProp("residue_number") and atom.HasProp("residue_name"):
        try:
            return (int(atom.GetProp("residue_number")), atom.GetProp("residue_name").strip())
        except (TypeError, ValueError):
            return None
    return None


def _interface_residue_keys(protein, ligand, threshold: float) -> Set[Tuple[int, str]]:
    """Protein residues (resnum, resname) with any atom within ``threshold`` of any ligand atom."""
    prot_conf = protein.GetConformer()
    lig_conf = ligand.GetConformer()
    lig_positions = [lig_conf.GetAtomPosition(a.GetIdx()) for a in ligand.GetAtoms()]
    thr2 = threshold * threshold
    keys: Set[Tuple[int, str]] = set()
    for patom in protein.GetAtoms():
        rkey = _residue_key_from_pdb_info(patom)
        if rkey is None:
            continue
        if rkey in keys:
            continue
        ppos = prot_conf.GetAtomPosition(patom.GetIdx())
        for lpos in lig_positions:
            dx = ppos.x - lpos.x
            dy = ppos.y - lpos.y
            dz = ppos.z - lpos.z
            if dx * dx + dy * dy + dz * dz < thr2:
                keys.add(rkey)
                break
    return keys


def _protein_backbone_coord_index(
    protein, allowed_atom_names: Set[str]
) -> Dict[Tuple[int, str, str], Tuple[float, float, float]]:
    """(resnum, resname, atomName) -> (x, y, z) for the requested backbone atoms."""
    conf = protein.GetConformer()
    out: Dict[Tuple[int, str, str], Tuple[float, float, float]] = {}
    for atom in protein.GetAtoms():
        anm = _atom_name_from_pdb_info(atom)
        if anm is None or anm not in allowed_atom_names:
            continue
        rkey = _residue_key_from_pdb_info(atom)
        if rkey is None:
            continue
        pos = conf.GetAtomPosition(atom.GetIdx())
        out[(rkey[0], rkey[1], anm)] = (pos.x, pos.y, pos.z)
    return out


def interface_protein_backbone_pairs(
    prot_ref,
    prot_dock,
    mol_ref,
    mol_dock,
    threshold: float = 10.0,
    atom_names=INTERFACE_PROTEIN_CA,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Paired interface protein backbone coordinates (ref, dock).

    Residues are paired by (resnum, resname) on residues that are interface on BOTH sides;
    atoms are paired by atom name (so no MCS is needed). ``atom_names`` can be CA-only or the
    full backbone. No superposition is applied (inputs are pre-aligned).
    """
    allowed = set(atom_names)
    ref_res = _interface_residue_keys(prot_ref, mol_ref, threshold)
    dock_res = _interface_residue_keys(prot_dock, mol_dock, threshold)
    shared = ref_res & dock_res
    ref_idx = _protein_backbone_coord_index(prot_ref, allowed)
    dock_idx = _protein_backbone_coord_index(prot_dock, allowed)
    coords_ref: List[Tuple[float, float, float]] = []
    coords_dock: List[Tuple[float, float, float]] = []
    for rid, rnm in sorted(shared):
        for anm in atom_names:
            k = (rid, rnm, anm)
            if k in ref_idx and k in dock_idx:
                coords_ref.append(ref_idx[k])
                coords_dock.append(dock_idx[k])
    a = np.array(coords_ref) if coords_ref else np.zeros((0, 3))
    b = np.array(coords_dock) if coords_dock else np.zeros((0, 3))
    return a, b


def _ligand_interface_ref_atom_indices(mol_ref, prot_ref, threshold: float) -> Set[int]:
    """Reference-ligand atom indices within ``threshold`` of any reference protein atom."""
    lconf = mol_ref.GetConformer()
    pconf = prot_ref.GetConformer()
    prot_positions = [pconf.GetAtomPosition(a.GetIdx()) for a in prot_ref.GetAtoms()]
    thr2 = threshold * threshold
    iface: Set[int] = set()
    for la in mol_ref.GetAtoms():
        lpos = lconf.GetAtomPosition(la.GetIdx())
        for ppos in prot_positions:
            dx = lpos.x - ppos.x
            dy = lpos.y - ppos.y
            dz = lpos.z - ppos.z
            if dx * dx + dy * dy + dz * dz < thr2:
                iface.add(la.GetIdx())
                break
    return iface


def interface_ligand_pairs_otmol(
    mol_ref,
    mol_dock,
    prot_ref,
    dock_to_ref: Dict[int, int],
    threshold: float = 10.0,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Paired interface ligand coordinates (ref, dock) using the OTMol dock->ref atom mapping.

    Replaces the legacy MCS ligand matching so the ligand correspondence is consistent with
    ligand RMSD / FNAT. Interface ligand atoms are defined on the reference side.
    """
    ref_to_dock: Dict[int, int] = {}
    for d_idx, r_idx in dock_to_ref.items():
        if r_idx not in ref_to_dock:
            ref_to_dock[r_idx] = d_idx
    iface_ref = _ligand_interface_ref_atom_indices(mol_ref, prot_ref, threshold)
    rconf = mol_ref.GetConformer()
    dconf = mol_dock.GetConformer()
    coords_ref: List[Tuple[float, float, float]] = []
    coords_dock: List[Tuple[float, float, float]] = []
    for r_idx in sorted(iface_ref):
        d_idx = ref_to_dock.get(r_idx)
        if d_idx is None:
            continue
        if d_idx < 0 or d_idx >= mol_dock.GetNumAtoms():
            continue
        rp = rconf.GetAtomPosition(r_idx)
        dp = dconf.GetAtomPosition(d_idx)
        coords_ref.append((rp.x, rp.y, rp.z))
        coords_dock.append((dp.x, dp.y, dp.z))
    a = np.array(coords_ref) if coords_ref else np.zeros((0, 3))
    b = np.array(coords_dock) if coords_dock else np.zeros((0, 3))
    return a, b


def calculate_interface_rmsd_otmol(
    mol_ref,
    mol_dock,
    prot_ref,
    prot_dock,
    dock_to_ref: Dict[int, int],
    threshold: float = 10.0,
    protein_atom_names=INTERFACE_PROTEIN_CA,
) -> dict:
    """
    Interface RMSD using OTMol ligand mapping (+ residue/atom-name protein pairing).

    ``protein_atom_names`` defaults to CA (legacy parity); pass
    :data:`INTERFACE_PROTEIN_BACKBONE` for N/CA/C/O. No superposition is applied.
    """
    pref, pdock = interface_protein_backbone_pairs(
        prot_ref, prot_dock, mol_ref, mol_dock, threshold=threshold, atom_names=protein_atom_names
    )
    lref, ldock = interface_ligand_pairs_otmol(
        mol_ref, mol_dock, prot_ref, dock_to_ref, threshold=threshold
    )
    value = combine_coords_and_calculate_rmsd(pref, pdock, lref, ldock)
    return {
        "interface_rmsd": value,
        "n_protein_atoms": int(pref.shape[0]),
        "n_ligand_atoms": int(lref.shape[0]),
        "protein_atom_names": tuple(protein_atom_names),
    }


def find_common_contact(sub_results1, sub_results2):
    s2 = set(sub_results2)
    return [c for c in sub_results1 if c in s2]


def is_backbone_atom(atom):
    backbone_atoms = {'N', 'CA', 'C', 'O'}
    return atom.GetProp("atom_name") in backbone_atoms if atom.HasProp("atom_name") else False


def _fnat_contact_tag(ref_atom_idx: int, residue_id: int, residue_name: str, atom_type: str = "all") -> str:
    """Canonical contact id keyed by reference ligand atom index (OTMol mapping)."""
    tag = f"mol_ref{ref_atom_idx}_res{residue_id}_{residue_name}"
    if atom_type != "all":
        tag += f"_{atom_type}"
    return tag


def _collect_fnat_contacts_otmol(
    mol_ref: Chem.Mol,
    mol_dock: Chem.Mol,
    prot_ref: Chem.Mol,
    prot_dock: Chem.Mol,
    dock_to_ref: Dict[int, int],
    threshold: float = 5.0,
    atom_type: str = "all",
) -> Tuple[List[str], List[str]]:
    """
    Build parallel ref/pred contact lists using OTMol dock_idx -> ref_idx mapping.
    Both sides use the same ref-atom-based tag so intersection = true positives.
    """
    contacts_ref: List[str] = []
    contacts_pred: List[str] = []
    conf_ref = mol_ref.GetConformer()
    conf_dock = mol_dock.GetConformer()
    prot_conf_ref = prot_ref.GetConformer()
    prot_conf_dock = prot_dock.GetConformer()

    for dock_idx, ref_idx in dock_to_ref.items():
        if ref_idx < 0 or ref_idx >= mol_ref.GetNumAtoms():
            continue
        if dock_idx < 0 or dock_idx >= mol_dock.GetNumAtoms():
            continue
        ref_pos = conf_ref.GetAtomPosition(ref_idx)
        dock_pos = conf_dock.GetAtomPosition(dock_idx)

        close_ref: Set[str] = set()
        for prot_atom in prot_ref.GetAtoms():
            if (ref_pos - prot_conf_ref.GetAtomPosition(prot_atom.GetIdx())).Length() >= threshold:
                continue
            ri = prot_atom.GetPDBResidueInfo()
            if not ri:
                continue
            if atom_type == "backbone" and not is_backbone_atom(prot_atom):
                continue
            if atom_type == "sidechain" and is_backbone_atom(prot_atom):
                continue
            close_ref.add(
                _fnat_contact_tag(
                    ref_idx,
                    ri.GetResidueNumber(),
                    ri.GetResidueName().strip(),
                    atom_type,
                )
            )
        contacts_ref.extend(close_ref)

        close_pred: Set[str] = set()
        for prot_atom in prot_dock.GetAtoms():
            if (dock_pos - prot_conf_dock.GetAtomPosition(prot_atom.GetIdx())).Length() >= threshold:
                continue
            ri = prot_atom.GetPDBResidueInfo()
            if not ri:
                continue
            if atom_type == "backbone" and not is_backbone_atom(prot_atom):
                continue
            if atom_type == "sidechain" and is_backbone_atom(prot_atom):
                continue
            close_pred.add(
                _fnat_contact_tag(
                    ref_idx,
                    ri.GetResidueNumber(),
                    ri.GetResidueName().strip(),
                    atom_type,
                )
            )
        contacts_pred.extend(close_pred)

    return contacts_ref, contacts_pred


def find_close_residues_otmol(
    mol_ref: Chem.Mol,
    mol_dock: Chem.Mol,
    prot_ref: Chem.Mol,
    prot_dock: Chem.Mol,
    dock_to_ref: Dict[int, int],
    threshold: float = 5.0,
) -> Tuple[List[str], List[str]]:
    """Residue-level FNAT contacts (OTMol ligand atom mapping only)."""
    return _collect_fnat_contacts_otmol(
        mol_ref, mol_dock, prot_ref, prot_dock, dock_to_ref, threshold, "all"
    )


def calculate_fnat_backbone_sidechain_otmol(
    mol_ref: Chem.Mol,
    mol_dock: Chem.Mol,
    prot_ref: Chem.Mol,
    prot_dock: Chem.Mol,
    dock_to_ref: Dict[int, int],
    threshold: float = 5.0,
) -> dict:
    """BBSC FNAT using OTMol dock_to_ref (no MCS ligand matching)."""
    if not dock_to_ref:
        empty: Set[str] = set()
        return {
            "fnat_all": -1,
            "fnat_backbone": -1,
            "fnat_sidechain": -1,
            "tp_all": 0,
            "tp_backbone": 0,
            "tp_sidechain": 0,
            "total_ref_all": 0,
            "total_ref_backbone": 0,
            "total_ref_sidechain": 0,
            "ref_contacts_all": 0,
            "pred_contacts_all": 0,
            "ref_contacts_backbone": 0,
            "pred_contacts_backbone": 0,
            "ref_contacts_sidechain": 0,
            "pred_contacts_sidechain": 0,
            "ref_contacts_all_set": empty,
            "pred_contacts_all_set": empty,
            "ref_contacts_backbone_set": empty,
            "pred_contacts_backbone_set": empty,
            "ref_contacts_sidechain_set": empty,
            "pred_contacts_sidechain_set": empty,
            "common_all": empty,
            "common_backbone": empty,
            "common_sidechain": empty,
            "fnat_map_src": "otmol",
            "fnat_mapped_pairs": 0,
        }

    ref_all, pred_all = _collect_fnat_contacts_otmol(
        mol_ref, mol_dock, prot_ref, prot_dock, dock_to_ref, threshold, "all"
    )
    ref_bb, pred_bb = _collect_fnat_contacts_otmol(
        mol_ref, mol_dock, prot_ref, prot_dock, dock_to_ref, threshold, "backbone"
    )
    ref_sc, pred_sc = _collect_fnat_contacts_otmol(
        mol_ref, mol_dock, prot_ref, prot_dock, dock_to_ref, threshold, "sidechain"
    )
    ref_all_s = set(ref_all)
    pred_all_s = set(pred_all)
    ref_bb_s = set(ref_bb)
    pred_bb_s = set(pred_bb)
    ref_sc_s = set(ref_sc)
    pred_sc_s = set(pred_sc)
    common_all = ref_all_s & pred_all_s
    common_bb = ref_bb_s & pred_bb_s
    common_sc = ref_sc_s & pred_sc_s
    fnat_all = len(common_all) / len(ref_all_s) if ref_all_s else -1
    fnat_backbone = len(common_bb) / len(ref_bb_s) if ref_bb_s else -1
    fnat_sidechain = len(common_sc) / len(ref_sc_s) if ref_sc_s else -1
    return {
        "fnat_all": fnat_all,
        "fnat_backbone": fnat_backbone,
        "fnat_sidechain": fnat_sidechain,
        "tp_all": len(common_all),
        "tp_backbone": len(common_bb),
        "tp_sidechain": len(common_sc),
        "total_ref_all": len(ref_all_s),
        "total_ref_backbone": len(ref_bb_s),
        "total_ref_sidechain": len(ref_sc_s),
        "ref_contacts_all": len(ref_all_s),
        "pred_contacts_all": len(pred_all_s),
        "ref_contacts_backbone": len(ref_bb_s),
        "pred_contacts_backbone": len(pred_bb_s),
        "ref_contacts_sidechain": len(ref_sc_s),
        "pred_contacts_sidechain": len(pred_sc_s),
        "ref_contacts_all_set": ref_all_s,
        "pred_contacts_all_set": pred_all_s,
        "ref_contacts_backbone_set": ref_bb_s,
        "pred_contacts_backbone_set": pred_bb_s,
        "ref_contacts_sidechain_set": ref_sc_s,
        "pred_contacts_sidechain_set": pred_sc_s,
        "common_all": common_all,
        "common_backbone": common_bb,
        "common_sidechain": common_sc,
        "fnat_map_src": "otmol",
        "fnat_mapped_pairs": len(dock_to_ref),
    }
