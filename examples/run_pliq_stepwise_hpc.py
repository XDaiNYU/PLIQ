#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Stepwise PLIQ ver5 driver for HPC.

This is the ver5 successor of ``run_plq_stepwise_hpc_ver*.py``. Behaviour and CLI are
the SAME as before, only the package was renamed ``plq`` -> ``pliq``. It computes the
RAW per-(case, model) metrics (ligand RMSD, interface RMSD, FNAT, BINANA counts,
PoseBusters, TM-score) and writes them to CSV. Turn those raw rows into eDockQ1 /
eDockQ6_1..4 afterwards with::

    python -c "import pandas as pd, pliq_edockq_scoring as S; \
               df=pd.read_csv('e-DockQ_stepwise_01.csv'); \
               S.score_edockq_stepwise(df).to_csv('e-DockQ_scored_01.csv', index=False)"

All ligand atom correspondence is OTMol-only (a single dock->ref atom-index map computed
once per pose and reused by ligand RMSD, interface RMSD, FNAT and every BINANA signature).
There is no MCS anywhere. The map is verified bijective + element-consistent, so the same
physical atom in the docked ligand maps to the same physical atom in the reference ligand.

HPC usage
---------
Either:
1. ``pip install .`` the ``pliq_ver6`` package, then run this script from anywhere
   (``import pliq`` already on path); or
2. Run from the distribution root (the directory that contains ``pliq/config.py``); this
   script also prepends that root automatically (see ``_ensure_pliq_importable``).

You still need ``import otmol`` (pip: ``git+https://github.com/weixiaoqimath/otmol.git``;
sibling ``otmol_package/`` is still auto-added as fallback) and
``obabel`` on PATH for BINANA / PoseBusters conversions.

BINANA uses built-in defaults. Override cutoffs with
``--interaction-criteria-xlsx /path/to/criteria.(xlsx|json)`` (needs ``openpyxl`` for xlsx).

Paths per case (model k):
  {path_ref_lig}/{pdb_id}/{pdb_id}_ligand_chainB.pdb
  {path_ref_pro}/{pdb_id}/{pdb_id}_protein_renum.pdb
  {path_docked}/posebusters_{pdb_id}/posebusters_{pdb_id}_model_{k}_aligned_ligand.pdb
  {path_docked}/posebusters_{pdb_id}/posebusters_{pdb_id}_model_{k}_aligned_protein.pdb

Outputs:
  --saved-csv: growing aggregate table
  --per-model-dir: one CSV per (pdb_id, af_model_id)
"""
from __future__ import annotations

import sys
from pathlib import Path


def _ensure_pliq_importable() -> None:
    """
    Prepend the PLIQ distribution root (directory where ``pliq/config.py`` exists),
    so this script works even when the package was not ``pip install``ed.

    Tries: script dir; parents up to great-grandparent; and ``<parent>/pliq`` so a copy
    living next to a sibling folder named ``pliq`` still resolves imports.
    """
    here = Path(__file__).resolve().parent
    roots: list[Path] = []
    p = here
    for _ in range(6):
        roots.append(p)
        p = p.parent
    roots.append(here.parent / "pliq")
    seen: set[str] = set()
    for root in roots:
        try:
            root = root.resolve()
        except OSError:
            continue
        key = str(root)
        if key in seen:
            continue
        seen.add(key)
        if (root / "pliq" / "config.py").is_file():
            if key not in sys.path:
                sys.path.insert(0, key)
            return


_ensure_pliq_importable()

import argparse
import os
import re
from typing import Dict, List, Optional

import pandas as pd
from tqdm import tqdm

from pliq.binana_interaction_criteria import (
    load_binana_get_all_interaction_kwargs_from_path,
)
from pliq.config import resolve_binana_python, resolve_tmscore_paths
from pliq.edockq import binana_recall as BR
from pliq.edockq import helpers as H
from pliq.edockq import ligand_mapping as LM
from pliq.edockq.otmol_runtime import ensure_local_otmol_on_path

# Add sibling otmol_package / cwd otmol_package to sys.path before import, so
# ``import otmol`` works without a pip-installed wheel.
ensure_local_otmol_on_path()

try:
    from posebusters import PoseBusters
except ImportError:
    PoseBusters = None

try:
    from rdkit.Chem.rdmolfiles import MolFromPDBFile
except ImportError:
    MolFromPDBFile = None

from pliq.tmscore_run.runner import _compile_tmscore, _parse_tmscore_output, _run_tmscore


def _safe_stem(pdb_id: str) -> str:
    s = str(pdb_id).strip()
    s = re.sub(r"[^\w.\-]+", "_", s)
    return s or "unknown"


def run_posebusters_one(
    *,
    case_name: str,
    af_model_id: int,
    path_ref_lig: Path,
    path_docked: Path,
    use_mol2: bool = True,
    fallback_pdb: bool = True,
) -> Optional[Dict]:
    """Single (case, model) PoseBusters redock; keys prefixed with pb_."""
    if PoseBusters is None:
        return None
    true_pdb_path = path_ref_lig / case_name / f"{case_name}_ligand_chainB.pdb"
    docked_folder = path_docked / f"posebusters_{case_name}"
    pred_pdb_path = docked_folder / f"posebusters_{case_name}_model_{af_model_id}_aligned_ligand.pdb"
    cond_file = docked_folder / f"posebusters_{case_name}_model_{af_model_id}_aligned_protein.pdb"
    pred_mol2_path = pred_pdb_path.with_suffix(".mol2")
    true_mol2_path = true_pdb_path.with_suffix(".mol2")

    if not pred_pdb_path.is_file() or not true_pdb_path.is_file() or not cond_file.is_file():
        return None

    buster = PoseBusters(config="redock")
    ok = False
    df_result = None

    if use_mol2:
        try:
            import subprocess

            subprocess.run(
                ["obabel", str(pred_pdb_path), "-O", str(pred_mol2_path), "--ph", "7.4"],
                check=True,
                capture_output=True,
                text=True,
                timeout=30,
            )
            subprocess.run(
                ["obabel", str(true_pdb_path), "-O", str(true_mol2_path), "--ph", "7.4"],
                check=True,
                capture_output=True,
                text=True,
                timeout=30,
            )
            df_result = buster.bust([str(pred_mol2_path)], str(true_mol2_path), str(cond_file), full_report=True)
            ok = True
        except Exception:
            pass
        finally:
            for p in (pred_mol2_path, true_mol2_path):
                if p.exists():
                    try:
                        p.unlink()
                    except OSError:
                        pass

    if not ok and fallback_pdb and MolFromPDBFile is not None:
        try:
            mol_pred = MolFromPDBFile(str(pred_pdb_path), sanitize=False, removeHs=True)
            mol_true = MolFromPDBFile(str(true_pdb_path), sanitize=False, removeHs=True)
            mol_cond = MolFromPDBFile(str(cond_file), sanitize=False, removeHs=True)
            if mol_pred is not None and mol_true is not None and mol_cond is not None:
                df_result = buster.bust([mol_pred], mol_true, mol_cond, full_report=True)
                ok = True
        except Exception:
            pass

    if not ok or df_result is None or df_result.empty:
        return None

    row = df_result.iloc[0].to_dict()
    return {f"pb_{k}": v for k, v in row.items()}


def run_tmscore_one(tm_exe: Path, ref_prot: Path, pred_prot: Path) -> dict:
    out = {
        "tm_Length1": "",
        "tm_Length2": "",
        "tm_CommonRes": "",
        "tm_RMSD_common": "",
        "tm_TMscore": "",
        "tm_MaxSub": "",
        "tm_GDT_TS": "",
        "tm_GDT_HA": "",
        "tm_error": "",
    }
    if not pred_prot.is_file():
        out["tm_error"] = "Missing pred protein"
        return out
    if not ref_prot.is_file():
        out["tm_error"] = "Missing ref protein"
        return out
    try:
        raw = _run_tmscore(tm_exe, pred_prot, ref_prot)
        if "GLIBC_" in raw or "GLIBCXX_" in raw:
            out["tm_error"] = "ABI error (glibc/libstdc++ mismatch)"
            return out
        p = _parse_tmscore_output(raw)
        out["tm_Length1"] = p.get("Length1", "")
        out["tm_Length2"] = p.get("Length2", "")
        out["tm_CommonRes"] = p.get("CommonRes", "")
        out["tm_RMSD_common"] = p.get("RMSD_common", "")
        out["tm_TMscore"] = p.get("TMscore", "")
        out["tm_MaxSub"] = p.get("MaxSub", "")
        out["tm_GDT_TS"] = p.get("GDT_TS", "")
        out["tm_GDT_HA"] = p.get("GDT_HA", "")
        if not out["tm_TMscore"]:
            out["tm_error"] = "Parse failed"
    except Exception as e:
        out["tm_error"] = str(e)[:200]
    return out


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(
        description="Stepwise PLIQ ver5: CSV + iloc slice + explicit ref/docked paths "
        "(no data_posebusters required).",
    )
    ap.add_argument("--csv-file", type=Path, required=True)
    ap.add_argument("--csv-sep", default="\t")
    ap.add_argument("--subfolder-col", default="Subfolder", help="Copied to pdb_id if pdb_id missing")
    ap.add_argument("--slice-start", type=int, default=0)
    ap.add_argument("--slice-stop", type=int, default=None, help="Exclusive iloc end; omit = all rows from slice-start")
    ap.add_argument("--path-ref-lig", type=Path, required=True)
    ap.add_argument("--path-ref-pro", type=Path, default=None)
    ap.add_argument("--path-docked", type=Path, required=True)
    ap.add_argument("--saved-csv", type=Path, default=Path("./e-DockQ_stepwise.csv"))
    ap.add_argument("--per-model-dir", type=Path, default=Path("./per_model"))
    ap.add_argument(
        "--tm-root",
        type=Path,
        default=None,
        help="Directory for resolve_tmscore_paths local TMscore lookup (default: --path-ref-lig). "
        "No data_posebusters required; bundled TMscore.cpp + cache used if no binary here.",
    )
    ap.add_argument("--n-models", type=int, default=20)
    ap.add_argument("--no-binana", action="store_true")
    ap.add_argument("--no-posebusters", action="store_true")
    ap.add_argument("--no-tmscore", action="store_true")
    ap.add_argument(
        "--interaction-criteria-xlsx",
        type=Path,
        default=None,
        help="Optional BINANA cutoffs override .xlsx/.json. If omitted, use BINANA built-in defaults.",
    )
    ap.add_argument(
        "--binana-criteria-tier",
        default="medium",
        choices=("low", "medium", "high", "medium1"),
        help="Tier: low=loose, medium=default (medium1=legacy alias), high=strict; default medium.",
    )
    ap.add_argument(
        "--binana-strip-h",
        action="store_true",
        help="Scheme B: obabel -d strip H from ref+dock protein+ligand. Default: scheme A (as-is).",
    )
    ap.add_argument(
        "--binana-obabel-ph",
        type=float,
        default=None,
        help="Optional Open Babel -p pH before BINANA. Default: omit -p.",
    )
    ap.add_argument(
        "--no-binana-obabel-h",
        action="store_true",
        help="Alias for default (no Open Babel -p).",
    )
    args = ap.parse_args(argv)

    LM.require_otmol_importable()
    print("run_pliq_stepwise: OTMol import OK (ligand RMSD + BINANA mapping are OTMol-only).", file=sys.stderr)

    path_ref_lig = Path(args.path_ref_lig).resolve()
    path_ref_pro = Path(args.path_ref_pro).resolve() if args.path_ref_pro else path_ref_lig
    path_docked = Path(args.path_docked).resolve()

    if not args.csv_file.is_file():
        print(f"Cases file not found: {args.csv_file}", file=sys.stderr)
        return 2

    # ----- STEP 1: load CSV -----
    df1 = pd.read_csv(args.csv_file, sep=args.csv_sep)

    # ----- STEP 2: pdb_id column -----
    if args.subfolder_col in df1.columns and "pdb_id" not in df1.columns:
        df1 = df1.copy()
        df1["pdb_id"] = df1[args.subfolder_col]

    # ----- STEP 3: slice (iloc[start:stop]) -----
    if args.slice_stop is None:
        df_subset = df1.iloc[args.slice_start :]
    else:
        df_subset = df1.iloc[args.slice_start : args.slice_stop]

    # TM-score resolution without PLIQ_EDOCKQ_ROOT / data_posebusters
    tm_root = Path(args.tm_root).resolve() if args.tm_root else path_ref_lig
    tm_dir, tm_cpp, tm_exe, compile_tm = resolve_tmscore_paths(tm_root)
    if compile_tm and not args.no_tmscore:
        _compile_tmscore(tm_dir, tm_cpp, tm_exe)

    binana_python = resolve_binana_python()
    binana_gkwargs = {}
    if not args.no_binana:
        if args.interaction_criteria_xlsx:
            crit = Path(args.interaction_criteria_xlsx).resolve()
            binana_gkwargs = load_binana_get_all_interaction_kwargs_from_path(
                crit, tier=args.binana_criteria_tier
            )
    binana_obabel_ph = None if args.no_binana_obabel_h else args.binana_obabel_ph
    binana_strip_h = args.binana_strip_h
    args.per_model_dir.mkdir(parents=True, exist_ok=True)
    args.saved_csv.parent.mkdir(parents=True, exist_ok=True)

    prl = str(path_ref_lig)
    prp = str(path_ref_pro)
    pdk = str(path_docked)

    results: List[dict] = []
    list_fail: List[str] = []

    # ----- STEP 4: nested loops -----
    for _index, row in tqdm(df_subset.iterrows(), total=len(df_subset), desc="cases"):
        case_name = row["pdb_id"]

        for af_model_id in range(0, args.n_models):
            alpha2, BCI2 = None, None
            assignment2 = None
            dock_to_ref: dict = {}
            otmol_bci_zero_ok = None
            otmol_align_attempts = None

            try:
                pdb_file1 = f"{prl}/{case_name}/{case_name}_ligand_chainB.pdb"
                chain_id1 = "B"
                pdb_file2 = f"{pdk}/posebusters_{case_name}/posebusters_{case_name}_model_{af_model_id}_aligned_ligand.pdb"
                chain_id2 = "B"

                prot_pdb_file1 = f"{prp}/{case_name}/{case_name}_protein_renum.pdb"
                prot_chain_id1 = "A"
                prot_pdb_file2 = f"{pdk}/posebusters_{case_name}/posebusters_{case_name}_model_{af_model_id}_aligned_protein.pdb"
                prot_chain_id2 = "A"

                # ----- STEP 4a: ligand RMSD (OTMol mapping) -----
                temp_lig = None
                try:
                    fd, temp_lig = os.mkstemp(suffix=".pdb")
                    os.close(fd)
                    mol1 = H.pdb_to_rdkit_with_prody(pdb_file1, chain_id1, temp_lig)
                except Exception:
                    mol1 = H.pdb_to_rdkit_with_prody(pdb_file1, chain_id1)
                finally:
                    if temp_lig and os.path.exists(temp_lig):
                        try:
                            os.remove(temp_lig)
                        except OSError:
                            pass

                temp_lig2 = None
                try:
                    fd2, temp_lig2 = os.mkstemp(suffix=".pdb")
                    os.close(fd2)
                    mol2 = H.pdb_to_rdkit_with_prody(pdb_file2, chain_id2, temp_lig2)
                except Exception:
                    mol2 = H.pdb_to_rdkit_with_prody(pdb_file2, chain_id2)
                finally:
                    if temp_lig2 and os.path.exists(temp_lig2):
                        try:
                            os.remove(temp_lig2)
                        except OSError:
                            pass

                ligand1_heavy_atoms = mol1.GetNumAtoms() if mol1 else -1
                ligand2_heavy_atoms = mol2.GetNumAtoms() if mol2 else -1
                rmsd = None

                if mol1 and mol2:
                    # All mappings (ligand RMSD, interface, FNAT, BINANA) are OTMol-only.
                    dock_to_ref, assignment2, rmsd, alpha2, BCI2, otmol_align_attempts, otmol_bci_zero_ok = (
                        LM.align_ligands_otmol_bci_constrained(mol1, mol2, max_attempts=20)
                    )

                ligand_rmsd = rmsd

                # ----- STEP 4b: interface RMSD -----
                prot1 = H.pdb_chain_to_rdkit_mol(prot_pdb_file1, prot_chain_id1)
                prot2 = H.pdb_chain_to_rdkit_mol(prot_pdb_file2, prot_chain_id2)
                if not prot1 or not prot2 or not mol1 or not mol2:
                    raise ValueError("Missing protein or ligand")

                interface_rmsd_value = H.calculate_interface_rmsd_otmol(
                    mol1, mol2, prot1, prot2, dock_to_ref,
                    threshold=10.0, protein_atom_names=H.INTERFACE_PROTEIN_CA,
                )["interface_rmsd"]

                # ----- STEP 4c–4d: FNAT (OTMol ligand atom mapping only) -----
                sub_results1, sub_results2 = H.find_close_residues_otmol(
                    mol1, mol2, prot1, prot2, dock_to_ref, threshold=5.0
                )
                common_contacts = H.find_common_contact(sub_results1, sub_results2)
                residue_fnat = len(common_contacts) / len(sub_results1) if sub_results1 else -1

                fnat_results = H.calculate_fnat_backbone_sidechain_otmol(
                    mol1, mol2, prot1, prot2, dock_to_ref, threshold=5.0
                )
                common_all = fnat_results["common_all"]
                common_backbone = fnat_results["common_backbone"]
                common_sidechain = fnat_results["common_sidechain"]

                rec: dict = {
                    "pdb_id": case_name,
                    "af_model_id": af_model_id,
                    "ligand_rmsd": ligand_rmsd,
                    "otmol_alpha": alpha2,
                    "otmol_BCI": BCI2,
                    "otmol_bci_zero_ok": otmol_bci_zero_ok,
                    "otmol_align_attempts": otmol_align_attempts,
                    "ligand1_heavy_atoms": ligand1_heavy_atoms,
                    "ligand2_heavy_atoms": ligand2_heavy_atoms,
                    "interface_rmsd": interface_rmsd_value,
                    "residue_fnat": residue_fnat,
                    "fnat_map_src": fnat_results.get("fnat_map_src", "otmol"),
                    "fnat_mapped_pairs": fnat_results.get("fnat_mapped_pairs", len(dock_to_ref)),
                    "fnat_all": fnat_results["fnat_all"],
                    "fnat_backbone": fnat_results["fnat_backbone"],
                    "fnat_sidechain": fnat_results["fnat_sidechain"],
                    "tp_all": fnat_results["tp_all"],
                    "tp_backbone": fnat_results["tp_backbone"],
                    "tp_sidechain": fnat_results["tp_sidechain"],
                    "total_ref_all": fnat_results["total_ref_all"],
                    "total_ref_backbone": fnat_results["total_ref_backbone"],
                    "total_ref_sidechain": fnat_results["total_ref_sidechain"],
                    "ref_contacts_all": fnat_results["ref_contacts_all"],
                    "pred_contacts_all": fnat_results["pred_contacts_all"],
                    "ref_contacts_backbone": fnat_results["ref_contacts_backbone"],
                    "pred_contacts_backbone": fnat_results["pred_contacts_backbone"],
                    "ref_contacts_sidechain": fnat_results["ref_contacts_sidechain"],
                    "pred_contacts_sidechain": fnat_results["pred_contacts_sidechain"],
                    "residue_ref_contacts": len(sub_results1),
                    "residue_pred_contacts": len(sub_results2),
                    "residue_common_contacts": len(common_contacts),
                    "ref_contacts_all_list": "; ".join(sorted(fnat_results["ref_contacts_all_set"])),
                    "pred_contacts_all_list": "; ".join(sorted(fnat_results["pred_contacts_all_set"])),
                    "common_contacts_all_list": "; ".join(sorted(common_all)),
                    "ref_contacts_backbone_list": "; ".join(sorted(fnat_results["ref_contacts_backbone_set"])),
                    "pred_contacts_backbone_list": "; ".join(sorted(fnat_results["pred_contacts_backbone_set"])),
                    "common_contacts_backbone_list": "; ".join(sorted(common_backbone)),
                    "ref_contacts_sidechain_list": "; ".join(sorted(fnat_results["ref_contacts_sidechain_set"])),
                    "pred_contacts_sidechain_list": "; ".join(sorted(fnat_results["pred_contacts_sidechain_set"])),
                    "common_contacts_sidechain_list": "; ".join(sorted(common_sidechain)),
                    "residue_ref_contacts_list": "; ".join(sorted(sub_results1)),
                    "residue_pred_contacts_list": "; ".join(sorted(sub_results2)),
                    "residue_common_contacts_list": "; ".join(sorted(common_contacts)),
                }
                rec.update(
                    LM.pliq_mapping_csv_columns(
                        mol_ref=mol1,
                        mol_dock=mol2,
                        dock_to_ref=dock_to_ref,
                        bci=BCI2,
                        bci_zero_ok=otmol_bci_zero_ok,
                    )
                )

                # ----- STEP 4e: BINANA -----
                rec["binana_error"] = ""
                bpath = Path(binana_python)
                if not args.no_binana:
                    try:
                        rec.update(
                            BR.compute_binana_recall_row(
                                pdb_file1,
                                prot_pdb_file1,
                                pdb_file2,
                                prot_pdb_file2,
                                bpath,
                                mol_ref=mol1,
                                mol_dock=mol2,
                                otmol_assignment=assignment2,
                                dock_to_ref=dock_to_ref,
                                otmol_BCI=BCI2,
                                otmol_bci_zero_ok=otmol_bci_zero_ok,
                                get_all_interactions_kwargs=binana_gkwargs,
                                obabel_ph=binana_obabel_ph,
                                strip_h=binana_strip_h,
                            )
                        )
                    except Exception as e_bin:
                        rec.update(BR.empty_binana_row())
                        rec["binana_error"] = str(e_bin)[:500]
                        rec.update(
                            LM.pliq_mapping_csv_columns(
                                mol_ref=mol1,
                                mol_dock=mol2,
                                dock_to_ref=dock_to_ref,
                                bci=BCI2,
                                bci_zero_ok=otmol_bci_zero_ok,
                            )
                        )
                else:
                    rec.update(BR.empty_binana_row())
                    rec["binana_error"] = "disabled"

                # ----- STEP 4f: PoseBusters -----
                if not args.no_posebusters:
                    pb = run_posebusters_one(
                        case_name=str(case_name),
                        af_model_id=af_model_id,
                        path_ref_lig=path_ref_lig,
                        path_docked=path_docked,
                    )
                    if pb:
                        rec.update(pb)

                # ----- STEP 4g: TM-score -----
                if not args.no_tmscore:
                    ref_prot = path_ref_pro / case_name / f"{case_name}_protein_renum.pdb"
                    pred_prot = path_docked / f"posebusters_{case_name}" / f"posebusters_{case_name}_model_{af_model_id}_aligned_protein.pdb"
                    rec.update(run_tmscore_one(tm_exe, ref_prot, pred_prot))

                results.append(rec)
                stem = _safe_stem(str(case_name))
                pm = args.per_model_dir / f"{stem}_model_{af_model_id}.csv"
                pd.DataFrame([rec]).to_csv(pm, index=False)
                pd.DataFrame(results).to_csv(args.saved_csv, index=False)

            except Exception as e:
                print(f"Error processing {case_name}_{af_model_id}: {e}")
                list_fail.append(f"{case_name}_{af_model_id}")

    pd.DataFrame(results).to_csv(args.saved_csv, index=False)
    print(f"Done. Wrote {len(results)} rows to {args.saved_csv}. Failed: {len(list_fail)}")
    if list_fail:
        print(list_fail[:30], "..." if len(list_fail) > 30 else "")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
