# -*- coding: utf-8 -*-
"""Command-line entry for PLIQ."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import List, Optional

from pliq import __version__
from pliq.pipeline import run_full_pipeline


def _default_cases_file() -> Path:
    return Path(__file__).resolve().parent.parent / "examples" / "cases_two.tsv"


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        prog="pliq",
        description="PL-Quality (PLIQ): e-DockQ + PoseBusters + TM-score + BINANA merge",
    )
    parser.add_argument("--version", action="version", version=f"pliq {__version__}")
    sub = parser.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("run-all", help="Run full pipeline and write CSVs + end2end_*_merged.csv per case")
    p.add_argument(
        "--cases",
        type=Path,
        default=None,
        help=f"TSV with Subfolder column (default: bundled examples/cases_two.tsv)",
    )
    p.add_argument(
        "--edockq-root",
        type=Path,
        default=None,
        help="Data root: must contain data_posebusters/. TM-score binary optional (see README). "
        "Env PLIQ_EDOCKQ_ROOT (legacy PLQ_EDOCKQ_ROOT); if unset, parents are searched.",
    )
    p.add_argument(
        "--output-dir",
        type=Path,
        default=Path("pliq_output"),
        help="Directory for e-DockQ_results.csv, posebusters_results.csv, TMscore_results.csv, merged_results.csv",
    )
    p.add_argument("--n-models", type=int, default=20)
    p.add_argument("--n-cases", type=int, default=0, help="0 = all rows in cases TSV")
    p.add_argument("--csv-sep", default="\t")
    p.add_argument("--subfolder-col", default="Subfolder")
    p.add_argument("--no-binana", action="store_true", help="Skip BINANA (faster)")
    p.add_argument(
        "--interaction-criteria-xlsx",
        type=Path,
        default=None,
        help="Optional BINANA cutoffs override: .xlsx (parameter/low/medium/high) or .json "
        "(flat or tiered). If omitted, use BINANA built-in defaults.",
    )
    p.add_argument(
        "--binana-criteria-tier",
        default="medium",
        choices=("low", "medium", "high", "medium1"),
        help="Tier: low=loose, medium=default (legacy: medium1=same), high=strict; default medium.",
    )
    p.add_argument(
        "--binana-strip-h",
        action="store_true",
        help="Scheme B: obabel -d strip H from ref+dock protein+ligand before BINANA. "
        "Default: scheme A (as-is, keep dock H if present). No obabel -p in either mode.",
    )
    p.add_argument(
        "--binana-obabel-ph",
        type=float,
        default=None,
        help="Optional Open Babel -p pH before BINANA (advanced; default: omit -p).",
    )
    p.add_argument(
        "--no-binana-obabel-h",
        action="store_true",
        help="Alias for default (no Open Babel -p). Kept for backward compatibility.",
    )

    sub.add_parser(
        "compile-tmscore",
        help="Compile bundled TMscore.cpp into the PLIQ cache (g++ on Linux, c++ on macOS; optional pre-step)",
    )

    sub.add_parser(
        "info",
        help="Print host OS (Linux/macOS/Windows), cache dir, data root, and TM-score resolution",
    )

    args = parser.parse_args(argv)

    if args.cmd == "info":
        from pliq.config import pliq_cache_dir, resolve_edockq_root, resolve_tmscore_paths
        from pliq.deps import OTMOL_GITHUB, OTMOL_INSTALL_CMD, TMSCORE_UPSTREAM, check_all
        from pliq.platinfo import describe_platform_line, tmscore_compile_argv

        print(describe_platform_line())
        print(f"PLIQ cache: {pliq_cache_dir()}")
        tm_exe = None
        try:
            root = resolve_edockq_root()
            print(f"PLIQ data root: {root}")
            tm_dir, tm_cpp, tm_exe, compile_on_run = resolve_tmscore_paths(root)
            print(f"TM-score: dir={tm_dir} cpp={tm_cpp} exe={tm_exe} compile_on_run={compile_on_run}")
            if compile_on_run:
                argv = tmscore_compile_argv(tm_cpp, tm_exe)
                print(f"TM-score compile (preview): {' '.join(argv)}")
        except FileNotFoundError as e:
            print(f"PLIQ data root: (not found) {e}")
        print()
        print("External dependencies:")
        print(f"  OTMol upstream:  {OTMOL_GITHUB}")
        print(f"  OTMol install:   {OTMOL_INSTALL_CMD}")
        print(f"  TM-score upstream: {TMSCORE_UPSTREAM} (no pip package; bundled TMscore.cpp)")
        for st in check_all(tm_exe):
            mark = "OK" if st.ok else "MISSING"
            print(f"  [{mark}] {st.name}: {st.detail}")
            if not st.ok and st.install_hint:
                print(f"         -> {st.install_hint}")
        return 0

    if args.cmd == "compile-tmscore":
        from pliq.config import compile_bundled_tmscore_to_cache

        try:
            exe = compile_bundled_tmscore_to_cache()
        except Exception as e:
            print(f"compile-tmscore failed: {e}", file=sys.stderr)
            return 1
        print(f"TM-score binary: {exe}")
        return 0

    if args.cmd == "run-all":
        cases = args.cases
        if cases is None:
            cases = _default_cases_file()
        if not cases.is_file():
            print(f"Cases file not found: {cases}", file=sys.stderr)
            return 2
        run_full_pipeline(
            csv_file=cases,
            edockq_root=args.edockq_root,
            output_dir=args.output_dir,
            n_models=args.n_models,
            n_cases=args.n_cases,
            csv_sep=args.csv_sep,
            subfolder_col=args.subfolder_col,
            run_binana=not args.no_binana,
            interaction_criteria_xlsx=args.interaction_criteria_xlsx,
            binana_criteria_tier=args.binana_criteria_tier,
            binana_obabel_ph=None if args.no_binana_obabel_h else args.binana_obabel_ph,
            binana_strip_h=args.binana_strip_h,
        )
        return 0

    return 1


if __name__ == "__main__":
    raise SystemExit(main())
