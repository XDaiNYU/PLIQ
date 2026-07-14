"""e-DockQ scoring: produces eDockQ1 and eDockQ6_1..4 only.

eDockQ1   = (fnat_BBSC + k_i + k_l) / 3
eDockQ6_k = term_4 × (B6_l × H_l + k_i + k_l)/3 × TM   (k = level 1..4)
  B6_l   = ONE micro-pooled F1 over the six non-hydrophobic BINANA types
           (sum tp/ref_n/dock_n across hbond…halogen, then precision/recall/F1)
  H_l    = F1(hydrophobic) at the same BINANA level (missing → 0)
  term_4 = pb_valid_fraction ** 2  (all 8 PoseBusters checks)
  k_i, k_l = interface / ligand RMSD kernels  1/(1+(r/1.5)^2), 1/(1+(r/8.5)^2)
  fnat   = BBSC F1 (backbone ∪ sidechain FNAT contacts)

BINANA levels (receptor signature granularity): default→6_1, rxatm→6_2, rxidx→6_3, rxres→6_4.
The ligand side of every BINANA signature uses the single OTMol dock→ref atom-index map.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

BINANA_TYPES = ["hbond", "salt", "pipi", "tstack", "cationpi", "halogen"]
BINANA_TYPES_PRF = BINANA_TYPES + ["hydrophobic"]
BINANA_LEVELS = [
    ("default", ""),
    ("rxidx", "_rxidx"),
    ("rxatm", "_rxatm"),
    ("rxres", "_rxres"),
]
# Notebook §3: default->_1, rxatm->_2, rxidx->_3, rxres->_4
EDOCK_LEVEL_MAP = [
    ("default", "1"),
    ("rxatm", "2"),
    ("rxidx", "3"),
    ("rxres", "4"),
]

PB_SOURCE = {
    "bond_lengths": "pb_bond_lengths",
    "bond_angles": "pb_bond_angles",
    "aromatic_ring_flatness": "pb_aromatic_ring_flatness",
    "double_bond_flatness": "pb_double_bond_flatness",
    "internal_steric_clash": "pb_internal_steric_clash",
    "energy_ratio": "pb_energy_ratio",
    "minimum_distance_to_protein": "pb_minimum_distance_to_protein",
    "volume_overlap_with_protein": "pb_volume_overlap_with_protein",
}


def parse_contact_list(x) -> set:
    if isinstance(x, set):
        return x
    if isinstance(x, float) and np.isnan(x):
        return set()
    if not isinstance(x, str):
        return set()
    if x.strip() == "":
        return set()
    return {s.strip() for s in x.split(";") if s.strip() != ""}


def force_set(x) -> set:
    if isinstance(x, set):
        return x
    if isinstance(x, str):
        return parse_contact_list(x)
    return set()


def compute_prf(df: pd.DataFrame, prefix: str) -> pd.DataFrame:
    out = df.copy()
    out[f"{prefix}_ref_set"] = out[f"{prefix}_ref_set"].apply(force_set)
    out[f"{prefix}_pred_set"] = out[f"{prefix}_pred_set"].apply(force_set)
    out[f"{prefix}_tp"] = out.apply(
        lambda r: len(r[f"{prefix}_ref_set"] & r[f"{prefix}_pred_set"]), axis=1
    )
    ref_n = out[f"{prefix}_ref_set"].apply(len)
    pred_n = out[f"{prefix}_pred_set"].apply(len)
    tp = out[f"{prefix}_tp"]
    out[f"{prefix}_precision"] = tp / pred_n.replace(0, np.nan)
    out[f"{prefix}_recall"] = tp / ref_n.replace(0, np.nan)
    out[f"{prefix}_f1"] = (
        2 * out[f"{prefix}_precision"] * out[f"{prefix}_recall"]
        / (out[f"{prefix}_precision"] + out[f"{prefix}_recall"])
    )
    return out


def compute_all_metrics(df1: pd.DataFrame) -> pd.DataFrame:
    """Notebook §1 — contact lists -> BBSC F1 (fnat)."""
    df1 = df1.copy()
    df1["ALL_ref_set"] = df1["ref_contacts_all_list"].apply(parse_contact_list)
    df1["ALL_pred_set"] = df1["pred_contacts_all_list"].apply(parse_contact_list)
    df1["BACKBONE_ref_set"] = df1["ref_contacts_backbone_list"].apply(parse_contact_list)
    df1["BACKBONE_pred_set"] = df1["pred_contacts_backbone_list"].apply(parse_contact_list)
    df1["SIDECHAIN_ref_set"] = df1["ref_contacts_sidechain_list"].apply(parse_contact_list)
    df1["SIDECHAIN_pred_set"] = df1["pred_contacts_sidechain_list"].apply(parse_contact_list)
    df1["RESATOM_ref_set"] = df1["residue_ref_contacts_list"].apply(parse_contact_list)
    df1["RESATOM_pred_set"] = df1["residue_pred_contacts_list"].apply(parse_contact_list)
    df1["BBSC_ref_set"] = df1.apply(
        lambda r: force_set(r["BACKBONE_ref_set"]) | force_set(r["SIDECHAIN_ref_set"]), axis=1
    )
    df1["BBSC_pred_set"] = df1.apply(
        lambda r: force_set(r["BACKBONE_pred_set"]) | force_set(r["SIDECHAIN_pred_set"]), axis=1
    )
    for mode in ["ALL", "BACKBONE", "SIDECHAIN", "BBSC", "RESATOM"]:
        df1 = compute_prf(df1, mode)
    return df1


def add_binana_prf_columns(df_in: pd.DataFrame) -> pd.DataFrame:
    """Notebook §4 — per-type per-level BINANA P/R/F1."""
    out = df_in.copy()
    for itype in BINANA_TYPES_PRF:
        for level_name, suffix in BINANA_LEVELS:
            ref_c = f"binana_{itype}_ref_n{suffix}"
            dock_c = f"binana_{itype}_dock_n{suffix}"
            tp_c = f"binana_{itype}_tp_n{suffix}"
            if ref_c not in out.columns or dock_c not in out.columns or tp_c not in out.columns:
                continue
            ref_n = pd.to_numeric(out[ref_c], errors="coerce")
            dock_n = pd.to_numeric(out[dock_c], errors="coerce")
            tp_n = pd.to_numeric(out[tp_c], errors="coerce")
            invalid = (ref_n > 0) & (dock_n.isna() | (dock_n == 0))
            dock_eff = dock_n.copy()
            tp_eff = tp_n.copy()
            dock_eff[invalid] = 1.0
            tp_eff[invalid] = 0.0
            prec = tp_eff / dock_eff.replace(0, np.nan)
            rec = tp_eff / ref_n.replace(0, np.nan)
            f1 = 2 * prec * rec / (prec + rec)
            # ref_n==0 → NaN (skipped in six-type product); ref_n>0 miss → F1=0
            f1 = f1.where(ref_n > 0, np.nan)
            f1 = f1.where((ref_n <= 0) | f1.notna(), 0.0)
            tag = f"binana_{itype}_{level_name}"
            out[f"{tag}_precision"] = prec
            out[f"{tag}_recall"] = rec
            out[f"{tag}_f1"] = f1
    return out


def _binana_safe_int(series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce").fillna(0).clip(lower=0)


def _row_sum_binana_counts(frame: pd.DataFrame, col_names: list[str]) -> pd.Series:
    acc = pd.Series(0.0, index=frame.index)
    for col in col_names:
        if col in frame.columns:
            acc = acc + _binana_safe_int(frame[col])
    return acc


def add_binana_six_pool_columns(df_in: pd.DataFrame) -> pd.DataFrame:
    """Pool the six non-hydrophobic BINANA types into ONE term per level (micro P/R/F1).

    The six types (hbond, salt, pipi, tstack, cationpi, halogen) are summed
    (Σtp, Σref_n, Σdock_n) and a single precision/recall/F1 is computed:
        precision = Σtp / Σdock_n,  recall = Σtp / Σref_n,  F1 = 2PR/(P+R)
    Edge cases (same convention as pliq_ver4):
        Σref_n == 0                        → P=R=F1 = 1.0  (nothing to reproduce)
        Σref_n > 0 and Σdock_n == 0        → 0
        Σref_n > 0, Σdock_n > 0, Σtp == 0  → 0
    Produces binana_six_{level}_micro_{precision,recall,f1}.
    """
    out = df_in.copy()
    for level_name, suffix in BINANA_LEVELS:
        cols_ref = [f"binana_{t}_ref_n{suffix}" for t in BINANA_TYPES]
        cols_dock = [f"binana_{t}_dock_n{suffix}" for t in BINANA_TYPES]
        cols_tp = [f"binana_{t}_tp_n{suffix}" for t in BINANA_TYPES]
        if not any(c in out.columns for c in cols_ref):
            continue

        sr = _row_sum_binana_counts(out, cols_ref)
        sd = _row_sum_binana_counts(out, cols_dock)
        st = _row_sum_binana_counts(out, cols_tp)

        prec = st / sd.replace(0, np.nan)
        rec = st / sr.replace(0, np.nan)
        f1 = 2 * prec * rec / (prec + rec)

        m_ref0 = sr == 0
        m_zero = ((sr > 0) & (sd == 0)) | ((sr > 0) & (sd > 0) & (st == 0))
        for s in (prec, rec, f1):
            s.loc[m_ref0] = 1.0
            s.loc[m_zero] = 0.0

        tag = f"binana_six_{level_name}_micro"
        out[f"{tag}_sum_ref_n"] = sr
        out[f"{tag}_sum_dock_n"] = sd
        out[f"{tag}_sum_tp_n"] = st
        out[f"{tag}_precision"] = prec
        out[f"{tag}_recall"] = rec
        out[f"{tag}_f1"] = f1
    return out


def _single_term_nan_to_zero(series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce").fillna(0.0)


def _pb_8_checks(row) -> list[bool]:
    return [
        bool(row["bond_lengths"]) if pd.notna(row["bond_lengths"]) else False,
        bool(row["bond_angles"]) if pd.notna(row["bond_angles"]) else False,
        bool(row["aromatic_ring_flatness"]) if pd.notna(row["aromatic_ring_flatness"]) else False,
        bool(row["double_bond_flatness"]) if pd.notna(row["double_bond_flatness"]) else False,
        bool(row["internal_steric_clash"]) if pd.notna(row["internal_steric_clash"]) else False,
        (float(row["energy_ratio"]) <= 100) if pd.notna(row["energy_ratio"]) else False,
        (float(row["minimum_distance_to_protein"]) > 0.75)
        if pd.notna(row["minimum_distance_to_protein"])
        else False,
        bool(row["volume_overlap_with_protein"]) if pd.notna(row["volume_overlap_with_protein"]) else False,
    ]


def compute_edockq1(df: pd.DataFrame) -> pd.Series:
    lrmsd = pd.to_numeric(df["ligand_rmsd"], errors="coerce")
    irmsd = pd.to_numeric(df["interface_rmsd_value"], errors="coerce")
    fnat = pd.to_numeric(df["fnat"], errors="coerce")
    return (fnat + 1 / (1 + (irmsd / 1.5) ** 2) + 1 / (1 + (lrmsd / 8.5) ** 2)) / 3


def score_edockq_stepwise(df_raw: pd.DataFrame) -> pd.DataFrame:
    """Compute only eDockQ1 and eDockQ6_1..4 (BINANA levels default/rxatm/rxidx/rxres).

    eDockQ1   = (fnat_BBSC + k_i + k_l) / 3
    eDockQ6_k = term_4 × (B6_l × H_l + k_i + k_l)/3 × TM
      B6_l = ONE micro-pooled F1 over the six non-hydrophobic BINANA types
             (Σtp, Σref_n, Σdock_n across hbond…halogen, then P/R/F1)
      H_l  = F1(hydrophobic) at the same level   (missing → 0)
      term_4 = (PoseBusters 8-check pass fraction)²
    """
    df = df_raw.copy()
    df = compute_all_metrics(df)
    df["fnat"] = df["BBSC_f1"]

    df["interface_rmsd_value"] = pd.to_numeric(df["interface_rmsd"], errors="coerce")
    df["TMscore"] = pd.to_numeric(df.get("tm_TMscore", df.get("TMscore")), errors="coerce")

    for short, col in PB_SOURCE.items():
        if col not in df.columns:
            raise KeyError(f"Missing PoseBusters column: {col}")
        df[short] = df[col]
    df = df.dropna(subset=list(PB_SOURCE.keys()), how="all")

    df = add_binana_prf_columns(df)   # per-type P/R/F1 (gives the hydrophobic term)
    df = add_binana_six_pool_columns(df)  # six non-hydrophobic types pooled into one F1

    lrmsd_b = pd.to_numeric(df["ligand_rmsd"], errors="coerce")
    irmsd_b = pd.to_numeric(df["interface_rmsd_value"], errors="coerce")
    ki_b = 1 / (1 + (irmsd_b / 1.5) ** 2)
    kl_b = 1 / (1 + (lrmsd_b / 8.5) ** 2)
    tm_b = pd.to_numeric(df["TMscore"], errors="coerce")
    df["kernel_i"] = ki_b
    df["kernel_l"] = kl_b
    df["eDockQ1"] = compute_edockq1(df)
    df["pb_valid_fraction"] = df.apply(lambda r: sum(_pb_8_checks(r)) / 8, axis=1)
    df["term_4"] = df["pb_valid_fraction"] ** 2
    t4_b = pd.to_numeric(df["term_4"], errors="coerce")

    # --- eDockQ6_1..4: (six-type pooled micro-F1) × hydrophobic F1, per BINANA level ---
    for level, num in EDOCK_LEVEL_MAP:
        b6_col = f"binana_six_{level}_micro_f1"
        hyd_col = f"binana_hydrophobic_{level}_f1"
        if b6_col not in df.columns:
            continue
        b6 = _single_term_nan_to_zero(df[b6_col])
        hterm = _single_term_nan_to_zero(df[hyd_col]) if hyd_col in df.columns else 0.0
        b6xh = b6 * hterm

        df[f"term_binana_six_micro_f1_{level}"] = b6
        df[f"term_binana_hydrophobic_f1_{level}"] = hterm
        df[f"fnat_binana6_{level}"] = b6
        df[f"fnat_binana_{level}"] = b6xh

        eq1b_h = (b6xh + ki_b + kl_b) / 3
        df[f"eDockQ1_binana_hydrophobic_{level}"] = eq1b_h
        df[f"eDockQ6_{num}"] = t4_b * eq1b_h * tm_b

    df = add_pliq_summary_columns(df)
    return df


# Default headline PLIQ score uses eDockQ6_1 (BINANA level: default / chain+resID+resName+atomName).
PLIQ_DEFAULT_SCORE_COL = "eDockQ6_1"


def add_pliq_summary_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Add ``pliq`` (default = eDockQ6_1) and ``pliq_term_*`` decomposition columns."""
    out = df.copy()
    if PLIQ_DEFAULT_SCORE_COL in out.columns:
        out["pliq"] = pd.to_numeric(out[PLIQ_DEFAULT_SCORE_COL], errors="coerce")
    if "eDockQ1" in out.columns:
        out["pliq_eDockQ1"] = pd.to_numeric(out["eDockQ1"], errors="coerce")

    term_map = {
        "pliq_term_fnat": "fnat",
        "pliq_term_kernel_i": "kernel_i",
        "pliq_term_kernel_l": "kernel_l",
        "pliq_term_B6": "term_binana_six_micro_f1_default",
        "pliq_term_H": "term_binana_hydrophobic_f1_default",
        "pliq_term_B6xH": "fnat_binana_default",
        "pliq_term_posebusters": "term_4",
        "pliq_term_TM": "TMscore",
    }
    for dst, src in term_map.items():
        if src in out.columns:
            out[dst] = pd.to_numeric(out[src], errors="coerce")
    return out


PLIQ_SUMMARY_COLUMNS = [
    "pliq",
    "pliq_eDockQ1",
    "eDockQ6_1",
    "eDockQ6_2",
    "eDockQ6_3",
    "eDockQ6_4",
    "pliq_term_fnat",
    "pliq_term_kernel_i",
    "pliq_term_kernel_l",
    "pliq_term_B6",
    "pliq_term_H",
    "pliq_term_B6xH",
    "pliq_term_posebusters",
    "pliq_term_TM",
]


PLIQ_DOCKQ_CORE = ["eDockQ1"]
PLIQ_DOCKQ_LEVEL_COLS = [f"eDockQ6_{n}" for n in "1234"]


def pliq_dockq_columns(df: pd.DataFrame) -> list[str]:
    """The final score columns produced by score_edockq_stepwise (eDockQ1 + eDockQ6_1..4)."""
    cols = [c for c in PLIQ_DOCKQ_CORE if c in df.columns]
    cols += [c for c in PLIQ_DOCKQ_LEVEL_COLS if c in df.columns]
    return cols


__all__ = [
    "BINANA_LEVELS",
    "BINANA_TYPES",
    "EDOCK_LEVEL_MAP",
    "PLIQ_DEFAULT_SCORE_COL",
    "PLIQ_SUMMARY_COLUMNS",
    "add_binana_prf_columns",
    "add_binana_six_pool_columns",
    "add_pliq_summary_columns",
    "compute_all_metrics",
    "compute_edockq1",
    "pliq_dockq_columns",
    "score_edockq_stepwise",
]
