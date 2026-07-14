# -*- coding: utf-8 -*-
"""End-to-end PLIQ: four PDB files in, scored CSV out."""
from __future__ import annotations

import os
import subprocess
from pathlib import Path
from typing import Optional, Union

import pandas as pd

from pliq.config import resolve_binana_python, resolve_edockq_root, resolve_tmscore_paths
from pliq.edockq.binana_recall import DEFAULT_BINANA_OBABEL_PH, DEFAULT_BINANA_STRIP_H
from pliq.edockq.pose_eval import evaluate_pose_from_pdbs


def _run_posebusters_row(
    ref_lig_pdb: Path,
    dock_lig_pdb: Path,
    dock_pro_pdb: Path,
    *,
    use_mol2: bool = True,
    fallback_pdb: bool = True,
) -> dict:
    try:
        from posebusters import PoseBusters
    except ImportError as e:
        raise ImportError("posebusters is required. Install with: pip install posebusters") from e

    try:
        from rdkit.Chem.rdmolfiles import MolFromPDBFile
    except ImportError:
        MolFromPDBFile = None

    buster = PoseBusters(config="redock")
    pred_mol2 = dock_lig_pdb.with_suffix(".mol2")
    true_mol2 = ref_lig_pdb.with_suffix(".mol2")
    ok = False
    last_err = ""

    if use_mol2:
        try:
            subprocess.run(
                ["obabel", str(dock_lig_pdb), "-O", str(pred_mol2), "--ph", "7.4"],
                check=True, capture_output=True, text=True, timeout=30,
            )
            subprocess.run(
                ["obabel", str(ref_lig_pdb), "-O", str(true_mol2), "--ph", "7.4"],
                check=True, capture_output=True, text=True, timeout=30,
            )
            df_result = buster.bust([str(pred_mol2)], str(true_mol2), str(dock_pro_pdb), full_report=True)
            ok = True
            row = df_result.iloc[0].to_dict()
        except Exception as e:
            last_err = f"mol2: {str(e)[:200]}"
        finally:
            for p in (pred_mol2, true_mol2):
                if p.exists():
                    try:
                        os.remove(p)
                    except OSError:
                        pass
    else:
        row = {}

    if not ok and fallback_pdb and MolFromPDBFile is not None:
        try:
            mol_pred = MolFromPDBFile(str(dock_lig_pdb), sanitize=False, removeHs=True)
            mol_true = MolFromPDBFile(str(ref_lig_pdb), sanitize=False, removeHs=True)
            mol_cond = MolFromPDBFile(str(dock_pro_pdb), sanitize=False, removeHs=True)
            if mol_pred is not None and mol_true is not None and mol_cond is not None:
                df_result = buster.bust([mol_pred], mol_true, mol_cond, full_report=True)
                row = df_result.iloc[0].to_dict()
                ok = True
        except Exception as e:
            last_err = f"{last_err}; pdb: {str(e)[:200]}" if last_err else f"pdb: {str(e)[:200]}"

    if not ok:
        raise RuntimeError(last_err or "posebusters_failed")
    return {f"pb_{k}": v for k, v in row.items()}


def _run_tmscore_row(
    ref_pro_pdb: Path,
    dock_pro_pdb: Path,
    *,
    edockq_root: Optional[Path] = None,
) -> dict:
    from pliq.tmscore_run.runner import _compile_tmscore, _parse_tmscore_output, _run_tmscore

    root = resolve_edockq_root(edockq_root)
    tm_dir, tm_cpp, tm_exe, compile_tm = resolve_tmscore_paths(root)
    if compile_tm:
        _compile_tmscore(tm_dir, tm_cpp, tm_exe)
    out = _run_tmscore(tm_exe, dock_pro_pdb, ref_pro_pdb)
    parsed = _parse_tmscore_output(out)
    return {
        "Length1": parsed.get("Length1", ""),
        "Length2": parsed.get("Length2", ""),
        "CommonRes": parsed.get("CommonRes", ""),
        "RMSD_common": parsed.get("RMSD_common", ""),
        "tm_TMscore": parsed.get("TMscore", ""),
        "MaxSub": parsed.get("MaxSub", ""),
        "GDT_TS": parsed.get("GDT_TS", ""),
        "GDT_HA": parsed.get("GDT_HA", ""),
    }


def run_pliq_from_pdbs(
    ref_lig_pdb: Union[str, Path],
    ref_pro_pdb: Union[str, Path],
    dock_lig_pdb: Union[str, Path],
    dock_pro_pdb: Union[str, Path],
    out_csv: Union[str, Path],
    *,
    pdb_id: str = "",
    af_model_id: int = 0,
    ref_lig_chain: str = "B",
    ref_pro_chain: str = "A",
    dock_lig_chain: str = "B",
    dock_pro_chain: str = "A",
    run_binana: bool = True,
    run_posebusters: bool = True,
    run_tmscore: bool = True,
    score_pliq: bool = True,
    binana_python: Optional[Union[str, Path]] = None,
    binana_obabel_ph: Optional[float] = DEFAULT_BINANA_OBABEL_PH,
    binana_strip_h: bool = DEFAULT_BINANA_STRIP_H,
    edockq_root: Optional[Path] = None,
) -> pd.DataFrame:
    """
    **Ultimate end-to-end PLIQ** for a single pose: four PDB files → one scored CSV row.

    Pipeline: OTMol e-DockQ + BINANA → PoseBusters → TM-score → eDockQ1 / eDockQ6 scoring.

    Uses OTMol from https://github.com/weixiaoqimath/otmol (reflection always off),
    BINANA built-in default cutoffs, optional PoseBusters + TM-score.

    BINANA hydrogen (no ``obabel -p``):
    - ``binana_strip_h=False`` (default, scheme A): ref as-is (no crystal H), dock keeps
      predictor protein/ligand H if present.
    - ``binana_strip_h=True`` (scheme B): ``obabel -d`` strips H from ref+dock protein+ligand.

    Parameters
    ----------
    ref_lig_pdb, ref_pro_pdb, dock_lig_pdb, dock_pro_pdb
        Paths to reference ligand/protein and docked ligand/protein PDB files.
    out_csv
        Output CSV path (one row).
    """
    ref_lig_pdb = Path(ref_lig_pdb).resolve()
    ref_pro_pdb = Path(ref_pro_pdb).resolve()
    dock_lig_pdb = Path(dock_lig_pdb).resolve()
    dock_pro_pdb = Path(dock_pro_pdb).resolve()
    out_csv = Path(out_csv).resolve()
    out_csv.parent.mkdir(parents=True, exist_ok=True)

    if binana_python is None:
        binana_python = resolve_binana_python()

    rec = evaluate_pose_from_pdbs(
        ref_lig_pdb,
        ref_pro_pdb,
        dock_lig_pdb,
        dock_pro_pdb,
        pdb_id=pdb_id,
        af_model_id=af_model_id,
        ref_lig_chain=ref_lig_chain,
        ref_pro_chain=ref_pro_chain,
        dock_lig_chain=dock_lig_chain,
        dock_pro_chain=dock_pro_chain,
        run_binana=run_binana,
        binana_python=binana_python,
        binana_obabel_ph=binana_obabel_ph,
        binana_strip_h=binana_strip_h,
    )

    if run_posebusters:
        try:
            rec.update(_run_posebusters_row(ref_lig_pdb, dock_lig_pdb, dock_pro_pdb))
        except Exception as e:
            rec["pb_error"] = str(e)[:500]
    else:
        rec["pb_error"] = "disabled"

    if run_tmscore:
        try:
            rec.update(_run_tmscore_row(ref_pro_pdb, dock_pro_pdb, edockq_root=edockq_root))
        except Exception as e:
            rec["tm_error"] = str(e)[:500]
            rec["tm_TMscore"] = ""
    else:
        rec["tm_error"] = "disabled"
        rec["tm_TMscore"] = ""

    df = pd.DataFrame([rec])
    can_score = score_pliq
    if can_score:
        from pliq_edockq_scoring import PB_SOURCE

        missing_pb = any(col not in df.columns for col in PB_SOURCE.values())
        if missing_pb or rec.get("pb_error") not in (None, "", "disabled"):
            can_score = False
        if "tm_TMscore" not in df.columns or rec.get("tm_error") not in (None, "", "disabled"):
            can_score = False

    if can_score:
        import pliq_edockq_scoring as S

        df = S.score_edockq_stepwise(df)

    df.to_csv(out_csv, index=False)
    return df
