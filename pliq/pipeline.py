# -*- coding: utf-8 -*-
"""End-to-end PLQ pipeline (same outputs as legacy run_all.py)."""
from __future__ import annotations

from pathlib import Path
from typing import Optional, Tuple, Union

import pandas as pd

from pliq.config import resolve_binana_python, resolve_edockq_root, resolve_tmscore_paths
from pliq.edockq.binana_recall import DEFAULT_BINANA_OBABEL_PH, DEFAULT_BINANA_STRIP_H


def run_full_pipeline(
    *,
    csv_file: Union[str, Path],
    edockq_root: Optional[Path] = None,
    output_dir: Union[str, Path],
    n_models: int = 20,
    n_cases: int = 0,
    csv_sep: str = "\t",
    subfolder_col: str = "Subfolder",
    run_binana: bool = True,
    interaction_criteria_xlsx: Optional[Union[str, Path]] = None,
    binana_criteria_tier: str = "medium",
    binana_obabel_ph: Optional[float] = DEFAULT_BINANA_OBABEL_PH,
    binana_strip_h: bool = DEFAULT_BINANA_STRIP_H,
) -> Tuple[Path, Path, Path, Path]:
    """
    Run e-DockQ (+ BINANA), PoseBusters, TM-score, then merge.

    Returns paths: (e_dockq_csv, posebusters_csv, tmscore_csv, merged_csv).
    """
    root = resolve_edockq_root(edockq_root)
    out_dir = Path(output_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    path_data = str(root / "data_posebusters")
    tm_dir, tm_cpp, tm_exe, compile_tm = resolve_tmscore_paths(root)
    binana_python = resolve_binana_python()

    out_edockq = out_dir / "e-DockQ_results.csv"
    out_pb = out_dir / "posebusters_results.csv"
    out_tm = out_dir / "TMscore_results.csv"
    out_merged = out_dir / "merged_results.csv"

    from pliq.edockq.runner import run_edockq
    from pliq.posebusters_run.runner import run_posebusters
    from pliq.tmscore_run.runner import run_tmscore

    csv_file = str(Path(csv_file).resolve())

    print("[1/3] Running e-DockQ (+ BINANA)...")
    res_edockq, fail_edockq = run_edockq(
        csv_file=csv_file,
        path_ref_lig=path_data,
        path_ref_pro=path_data,
        path_docked=path_data,
        out_csv=str(out_edockq),
        n_models=n_models,
        n_cases=n_cases if n_cases else 0,
        csv_sep=csv_sep,
        subfolder_col=subfolder_col,
        run_binana=run_binana,
        binana_python=str(binana_python),
        interaction_criteria_xlsx=str(interaction_criteria_xlsx)
        if interaction_criteria_xlsx
        else None,
        binana_criteria_tier=binana_criteria_tier,
        binana_obabel_ph=binana_obabel_ph,
        binana_strip_h=binana_strip_h,
    )
    print(f"  -> {out_edockq} ({len(res_edockq)} rows, {len(fail_edockq)} failed)")

    print("[2/3] Running PoseBusters...")
    res_pb, fail_pb = run_posebusters(
        csv_file=csv_file,
        path_ref_lig=path_data,
        path_docked=path_data,
        out_csv=str(out_pb),
        n_models=n_models,
        n_cases=n_cases if n_cases else 0,
        csv_sep=csv_sep,
        subfolder_col=subfolder_col,
    )
    print(f"  -> {out_pb} ({len(res_pb)} rows, {len(fail_pb)} failed)")

    print("[3/3] Running TM-score...")
    if compile_tm:
        print(f"  -> Building TM-score binary: {tm_exe} (g++ + bundled TMscore.cpp)")
    df_tm = run_tmscore(
        csv_file=csv_file,
        path_ref_pro=path_data,
        path_docked=path_data,
        out_csv=str(out_tm),
        tm_dir=str(tm_dir),
        tm_cpp=tm_cpp,
        tm_exe=tm_exe,
        n_models=n_models,
        n_cases=n_cases if n_cases else 0,
        csv_sep=csv_sep,
        subfolder_col=subfolder_col,
        compile_on_run=compile_tm,
    )
    print(f"  -> {out_tm} ({len(df_tm)} rows)")

    print("[Merge] Writing merged CSV...")
    merge_cols = ["pdb_id", "af_model_id"]
    df_edockq = pd.DataFrame(res_edockq) if res_edockq else pd.DataFrame()
    df_pb = pd.concat(res_pb, ignore_index=True) if res_pb else pd.DataFrame()
    merged = df_tm.copy()
    if not df_edockq.empty:
        edockq_cols = [c for c in df_edockq.columns if c not in merge_cols]
        merged = merged.merge(df_edockq[merge_cols + edockq_cols], on=merge_cols, how="outer")
    if not df_pb.empty and "pdb_id" in df_pb.columns and "af_model_id" in df_pb.columns:
        pb_cols = [c for c in df_pb.columns if c not in merge_cols]
        df_pb_renamed = df_pb[merge_cols + pb_cols].rename(columns={c: f"pb_{c}" for c in pb_cols})
        merged = merged.merge(df_pb_renamed, on=merge_cols, how="outer")
    merged.to_csv(out_merged, index=False)
    print(f"  -> {out_merged} ({len(merged)} rows)")

    # Per-case merged copies (same naming as before)
    for pid in merged["pdb_id"].dropna().unique():
        sub = merged[merged["pdb_id"] == pid]
        end2end = out_dir / f"end2end_{pid}_merged.csv"
        sub.to_csv(end2end, index=False)
        print(f"  -> {end2end} ({len(sub)} rows)")

    print("Done.")
    return out_edockq, out_pb, out_tm, out_merged
