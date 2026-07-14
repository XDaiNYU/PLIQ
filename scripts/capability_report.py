#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PLIQ capability report.

For every metric *term* the pipeline emits, probe whether the required package /
binary is actually runnable on this host. PLIQ never fabricates or randomly fills a
term: if a term's engine is missing, it raises / records an explicit error instead of
silently skipping or inventing numbers (per-case failures are a separate concern).

Run:
    python pliq_ver5/scripts/capability_report.py
"""
from __future__ import annotations

import importlib
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from pliq.edockq.otmol_runtime import ensure_local_otmol_on_path


def _pkg(name):
    try:
        m = importlib.import_module(name)
        return True, getattr(m, "__version__", "")
    except Exception as e:  # noqa: BLE001
        return False, str(e)[:60]


def _otmol():
    ensure_local_otmol_on_path()
    return _pkg("otmol")


def _binary(name):
    p = shutil.which(name)
    return (p is not None), (p or "not found")


def _binana_python():
    try:
        from pliq.config import resolve_binana_python

        p = resolve_binana_python()
        return p.is_dir(), str(p)
    except Exception as e:  # noqa: BLE001
        return False, str(e)[:60]


def _tmscore():
    try:
        from pliq.config import pliq_cache_dir, resolve_tmscore_paths

        # Resolve against bundled cpp + cache (no data root needed for capability).
        tm_dir, tm_cpp, tm_exe, compile_on_run = resolve_tmscore_paths(pliq_cache_dir())
        cxx = shutil.which("g++") or shutil.which("c++") or shutil.which("clang++")
        runnable = tm_exe.is_file() or (tm_cpp.is_file() and cxx is not None)
        detail = f"exe={'present' if tm_exe.is_file() else 'compile-on-run'} cpp={tm_cpp.is_file()} cxx={cxx or 'MISSING'}"
        return runnable, detail
    except Exception as e:  # noqa: BLE001
        return False, str(e)[:80]


def main():
    rdkit_ok, rdkit_v = _pkg("rdkit")
    prody_ok, prody_v = _pkg("prody")
    numpy_ok, numpy_v = _pkg("numpy")
    ot_ok, ot_v = _pkg("ot")
    otmol_ok, otmol_v = _otmol()
    pb_ok, pb_v = _pkg("posebusters")
    obabel_ok, obabel_p = _binary("obabel")
    binana_ok, binana_p = _binana_python()
    tm_ok, tm_detail = _tmscore()

    # term -> (engine summary, runnable?, mapping engine)
    terms = [
        ("ligand_rmsd", "OTMol + POT + RDKit + ProDy + numpy",
         all([otmol_ok, ot_ok, rdkit_ok, prody_ok, numpy_ok]), "OTMol"),
        ("interface_rmsd / interface_rmsd_backbone", "RDKit + ProDy + numpy + OTMol map",
         all([rdkit_ok, prody_ok, numpy_ok, otmol_ok]), "OTMol (ligand) + residue+atom-name (protein)"),
        ("residue_fnat", "RDKit + OTMol map",
         all([rdkit_ok, otmol_ok]), "OTMol"),
        ("fnat_all / fnat_backbone / fnat_sidechain", "RDKit + OTMol map",
         all([rdkit_ok, otmol_ok]), "OTMol"),
        ("binana_* (recall, all interaction types)", "obabel + vendored BINANA + OTMol map",
         all([obabel_ok, binana_ok, otmol_ok]), "OTMol"),
        ("pb_* (PoseBusters redock)", "posebusters + obabel",
         all([pb_ok, obabel_ok]), "n/a"),
        ("tm_* (TM-score / GDT)", "TMscore binary (bundled cpp + g++/c++)",
         tm_ok, "n/a"),
    ]

    print("=" * 92)
    print("PLIQ DEPENDENCIES")
    print("=" * 92)
    for label, (ok, v) in [
        ("rdkit", (rdkit_ok, rdkit_v)), ("prody", (prody_ok, prody_v)),
        ("numpy", (numpy_ok, numpy_v)), ("POT (ot)", (ot_ok, ot_v)),
        ("otmol", (otmol_ok, otmol_v)), ("posebusters", (pb_ok, pb_v)),
        ("obabel (bin)", (obabel_ok, obabel_p)), ("BINANA (vendored)", (binana_ok, binana_p)),
        ("TM-score", (tm_ok, tm_detail)),
    ]:
        print(f"  [{'OK ' if ok else 'MISSING'}] {label:20s} {v}")

    print("\n" + "=" * 92)
    print(f"{'TERM':44s} {'RUNNABLE':9s} {'MAPPING':28s}")
    print("=" * 92)
    all_ok = True
    for term, engine, runnable, mapping in terms:
        all_ok = all_ok and runnable
        print(f"{term:44s} {'YES' if runnable else 'NO':9s} {mapping:28s}")
        print(f"    engine: {engine}")
    print("=" * 92)
    print("Every term above is computed by a real engine; on engine failure PLIQ raises or")
    print("records an explicit *_error / NaN sentinel \u2014 it never randomly fills a value.")
    print("RESULT:", "ALL TERMS RUNNABLE ON THIS HOST" if all_ok else "SOME TERMS NOT RUNNABLE (see MISSING above)")
    return 0 if all_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
