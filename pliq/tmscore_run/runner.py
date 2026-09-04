# -*- coding: utf-8 -*-
"""Run TM-score (pred vs ref protein) over CSV cases and save results to CSV."""
import re
import subprocess
from pathlib import Path

import pandas as pd
from tqdm import tqdm


def _compile_tmscore(tm_dir, tm_cpp, tm_exe):
    """Compile TMscore.cpp; OS-specific command from :mod:`plq.platinfo`."""
    from pliq.platinfo import tmscore_compile_argv

    tm_cpp = Path(tm_cpp)
    tm_exe = Path(tm_exe)
    if not tm_cpp.is_file():
        raise FileNotFoundError(f"TMscore.cpp not found: {tm_cpp}")
    cmd = tmscore_compile_argv(tm_cpp, tm_exe)
    subprocess.run(cmd, check=True, timeout=120)


def _run_tmscore(tm_exe, pred_pdb, ref_pdb, timeout=60):
    """Run TMscore and return raw stdout+stderr."""
    result = subprocess.run(
        [str(tm_exe), str(pred_pdb), str(ref_pdb)],
        capture_output=True,
        text=True,
        timeout=timeout,
    )
    return result.stdout + result.stderr


def resolve_tmscore_protein_paths(ref_pro_pdb, pred_pro_pdb):
    """Prefer *_full.pdb (ver5-style whole protein) when present; else trimmed interface files."""
    ref_pro_pdb = Path(ref_pro_pdb)
    pred_pro_pdb = Path(pred_pro_pdb)
    ref_full = ref_pro_pdb.with_name(
        ref_pro_pdb.name.replace("_protein_renum.pdb", "_protein_renum_full.pdb")
    )
    pred_full = pred_pro_pdb.with_name(
        pred_pro_pdb.name.replace("_aligned_protein.pdb", "_aligned_protein_full.pdb")
    )
    return (
        ref_full if ref_full.is_file() else ref_pro_pdb,
        pred_full if pred_full.is_file() else pred_pro_pdb,
    )


def _parse_tmscore_output(out):
    """Parse TMscore stdout for Length1, Length2, CommonRes, RMSD, TM-score, MaxSub, GDT-TS, GDT-HA."""
    data = {}
    # Structure1: Length= xxx
    m1 = re.search(r"Structure1:\s*.*?Length=\s*(\d+)", out)
    data["Length1"] = m1.group(1) if m1 else ""
    m2 = re.search(r"Structure2:\s*.*?Length=\s*(\d+)", out)
    data["Length2"] = m2.group(1) if m2 else ""
    # Number of residues in common = xxx
    m3 = re.search(r"Number of residues in common\s*=\s*(\d+)", out)
    data["CommonRes"] = m3.group(1) if m3 else ""
    # RMSD of  the common residues = xxx
    m4 = re.search(r"RMSD of\s+.*?=\s*([\d.]+)", out)
    data["RMSD_common"] = m4.group(1) if m4 else ""
    # TM-score = 0.xxxx (lines that start with optional spaces then TM-score)
    m5 = re.search(r"^\s*TM-score\s*=\s*([01]\.[\d]+)", out, re.MULTILINE)
    data["TMscore"] = m5.group(1) if m5 else ""
    m6 = re.search(r"^\s*MaxSub-score\s*=\s*([01]\.[\d]+)", out, re.MULTILINE)
    data["MaxSub"] = m6.group(1) if m6 else ""
    m7 = re.search(r"^\s*GDT-TS-score\s*=\s*([01]\.[\d]+)", out, re.MULTILINE)
    data["GDT_TS"] = m7.group(1) if m7 else ""
    m8 = re.search(r"^\s*GDT-HA-score\s*=\s*([01]\.[\d]+)", out, re.MULTILINE)
    data["GDT_HA"] = m8.group(1) if m8 else ""
    return data


def run_tmscore(
    csv_file,
    path_ref_pro,
    path_docked,
    out_csv,
    tm_dir,
    tm_cpp,
    tm_exe,
    n_models=20,
    n_cases=0,
    csv_sep="\t",
    subfolder_col="Subfolder",
    compile_on_run=True,
):
    """
    Compile TMscore (if compile_on_run), loop over cases and models, run TMscore(pred_protein, ref_protein),
    parse output, append to CSV.
    """
    tm_dir = Path(tm_dir)
    tm_cpp = Path(tm_cpp)
    tm_exe = Path(tm_exe)
    if compile_on_run:
        _compile_tmscore(tm_dir, tm_cpp, tm_exe)

    df1 = pd.read_csv(csv_file, sep=csv_sep)
    if n_cases and n_cases > 0:
        df1 = df1.iloc[: n_cases]
    if subfolder_col in df1.columns and "pdb_id" not in df1.columns:
        df1 = df1.copy()
        df1["pdb_id"] = df1[subfolder_col]

    Path(out_csv).parent.mkdir(parents=True, exist_ok=True)
    path_ref_pro = Path(path_ref_pro)
    path_docked = Path(path_docked)

    rows = []
    for _, row in tqdm(df1.iterrows(), total=len(df1), desc="TM-score"):
        case_name = row["pdb_id"]
        ref_prot, _ = resolve_tmscore_protein_paths(
            path_ref_pro / case_name / f"{case_name}_protein_renum.pdb",
            path_docked / f"posebusters_{case_name}" / f"posebusters_{case_name}_model_0_aligned_protein.pdb",
        )
        for mid in range(n_models):
            _, pred_prot = resolve_tmscore_protein_paths(
                path_ref_pro / case_name / f"{case_name}_protein_renum.pdb",
                path_docked / f"posebusters_{case_name}" / f"posebusters_{case_name}_model_{mid}_aligned_protein.pdb",
            )
            Length1 = Length2 = CommonRes = RMSD_common = TMscore = MaxSub = GDT_TS = GDT_HA = ""
            error = ""
            if not pred_prot.exists():
                error = "Missing pred"
            elif not ref_prot.exists():
                error = "Missing ref"
            else:
                try:
                    out = _run_tmscore(tm_exe, pred_prot, ref_prot)
                    if "GLIBC_" in out or "GLIBCXX_" in out:
                        error = "ABI error (glibc/libstdc++ mismatch). Compile on compute node."
                    else:
                        parsed = _parse_tmscore_output(out)
                        Length1 = parsed.get("Length1", "")
                        Length2 = parsed.get("Length2", "")
                        CommonRes = parsed.get("CommonRes", "")
                        RMSD_common = parsed.get("RMSD_common", "")
                        TMscore = parsed.get("TMscore", "")
                        MaxSub = parsed.get("MaxSub", "")
                        GDT_TS = parsed.get("GDT_TS", "")
                        GDT_HA = parsed.get("GDT_HA", "")
                        if not TMscore:
                            error = "Parse failed"
                except subprocess.TimeoutExpired:
                    error = "Timeout"
                except Exception as e:
                    error = str(e)[:200]
            rows.append({
                "pdb_id": case_name,
                "af_model_id": mid,
                "pred_protein_pdb": str(pred_prot),
                "ref_protein_pdb": str(ref_prot),
                "Length1": Length1,
                "Length2": Length2,
                "CommonRes": CommonRes,
                "RMSD_common": RMSD_common,
                "TMscore": TMscore,
                "MaxSub": MaxSub,
                "GDT_TS": GDT_TS,
                "GDT_HA": GDT_HA,
                "error": error,
            })

    df_out = pd.DataFrame(rows)
    df_out.to_csv(out_csv, index=False)
    return df_out
