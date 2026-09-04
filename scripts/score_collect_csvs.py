#!/usr/bin/env python3
"""Score stepwise PLIQ CSVs with score_edockq_stepwise (adds pliq = eDockQ6_1)."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
PLIQ_PKG = SCRIPT_DIR.parent
if str(PLIQ_PKG) not in sys.path:
    sys.path.insert(0, str(PLIQ_PKG))

from pliq_edockq_scoring import PLIQ_SUMMARY_COLUMNS, score_edockq_stepwise


def score_one(src: Path, dst: Path, *, overwrite: bool) -> None:
    if dst.is_file() and not overwrite:
        print(f"skip (exists): {dst}")
        return
    df = score_edockq_stepwise(pd.read_csv(src, low_memory=False))
    dst.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(dst, index=False)
    pliq_ok = df["pliq"].notna().sum() if "pliq" in df.columns else 0
    print(f"wrote {dst}  rows={len(df)}  pliq non-null={pliq_ok}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--collect-dir",
        type=Path,
        default=PLIQ_PKG.parent / "pliq_ver6_analysis" / "PLIQ_ver6_collect",
        help="Root folder with stepwise CSVs",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=None,
        help="Output directory (default: <collect-dir>/scored)",
    )
    parser.add_argument("--overwrite", action="store_true", help="Replace existing scored CSVs")
    parser.add_argument("--glob", default="*.csv", help="Only score files matching this glob")
    args = parser.parse_args()

    collect = args.collect_dir.resolve()
    out_root = (args.out_dir or collect / "scored").resolve()
    if not collect.is_dir():
        raise SystemExit(f"collect dir not found: {collect}")

    paths = sorted(p for p in collect.rglob(args.glob) if "scored" not in p.parts)
    if not paths:
        raise SystemExit(f"no CSV files under {collect}")

    print(f"collect: {collect}")
    print(f"out:     {out_root}")
    print(f"pliq columns added: {PLIQ_SUMMARY_COLUMNS[:6]} ...")
    for src in paths:
        rel = src.relative_to(collect)
        if rel.parts[0] == "scored":
            continue
        dst = out_root / rel.with_name(f"{rel.stem}_scored.csv")
        score_one(src, dst, overwrite=args.overwrite)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
