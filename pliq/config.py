# -*- coding: utf-8 -*-
"""Path resolution for PLQ."""
from __future__ import annotations

import os
from pathlib import Path
from typing import Optional, Tuple

from pliq.platinfo import getenv_compat, is_mac, is_windows

# Package directory: .../plq/plq/
PACKAGE_DIR = Path(__file__).resolve().parent
VENDOR_BINANA = PACKAGE_DIR / "_vendor" / "binana_python"
VENDOR_TMSCORE_CPP = PACKAGE_DIR / "_vendor" / "tmscore" / "TMscore.cpp"


def pliq_cache_dir() -> Path:
    """User cache directory for auto-built TMscore binary."""
    env = getenv_compat("CACHE")
    if env:
        return Path(env).resolve()
    if is_mac():
        base = Path.home() / "Library" / "Caches" / "pliq"
    elif is_windows():
        local = os.environ.get("LOCALAPPDATA") or str(Path.home() / "AppData" / "Local")
        base = Path(local) / "pliq" / "cache"
    else:
        xdg = os.environ.get("XDG_CACHE_HOME", "").strip()
        base = Path(xdg) / "pliq" if xdg else Path.home() / ".cache" / "pliq"
    base.mkdir(parents=True, exist_ok=True)
    return base


def resolve_edockq_root(explicit: Optional[Path] = None) -> Path:
    """
    eDockQ (or data) repository root: must contain ``data_posebusters/``.
    TM-score binary may live here but is optional if the bundled ``TMscore.cpp`` can be compiled.

    Resolution order:
    1. ``explicit`` argument
    2. ``PLIQ_EDOCKQ_ROOT`` (legacy ``PLQ_EDOCKQ_ROOT``)
    3. Walk parents until ``data_posebusters`` is found.
    """
    if explicit is not None:
        root = Path(explicit).resolve()
        if not (root / "data_posebusters").is_dir():
            raise FileNotFoundError(f"Not a valid PLIQ data root (no data_posebusters): {root}")
        return root

    env = getenv_compat("EDOCKQ_ROOT")
    if env:
        root = Path(env).resolve()
        if not (root / "data_posebusters").is_dir():
            raise FileNotFoundError(f"PLIQ_EDOCKQ_ROOT has no data_posebusters: {root}")
        return root

    here = Path(__file__).resolve()
    for p in here.parents:
        if (p / "data_posebusters").is_dir():
            return p

    raise FileNotFoundError(
        "Could not locate data root (need data_posebusters/). "
        "Set PLIQ_EDOCKQ_ROOT or use: pliq run-all --edockq-root /path/to/eDockQ"
    )


def _executable_candidate(root: Path) -> Path:
    if is_windows():
        return root / "TMscore.exe"
    return root / "TMscore"


def _is_runnable_binary(p: Path) -> bool:
    if not p.is_file():
        return False
    if is_windows():
        return True
    return os.access(p, os.X_OK)


def resolve_tmscore_paths(edockq_root: Path) -> Tuple[Path, Path, Path, bool]:
    """
    Resolve TM-score ``tm_dir``, ``tm_cpp``, ``tm_exe``, and whether to compile before run.

    Priority:
    1. ``PLIQ_TM_EXE`` (legacy ``PLQ_TM_EXE``) — use this binary; ``PLIQ_TM_CPP`` optional
       (defaults: edockq root, else bundled).
    2. Executable ``TMscore`` (or ``TMscore.exe``) under ``edockq_root``.
    3. Bundled ``TMscore.cpp`` → compile into ``pliq_cache_dir()/TMscore`` (needs ``g++``).

    Returns:
        (tm_dir, tm_cpp, tm_exe, compile_on_run)
    """
    bundled_cpp = VENDOR_TMSCORE_CPP
    if not bundled_cpp.is_file():
        raise FileNotFoundError(f"Bundled TMscore.cpp missing: {bundled_cpp}")

    env_exe = getenv_compat("TM_EXE")
    if env_exe:
        tm_exe = Path(env_exe).resolve()
        env_cpp = getenv_compat("TM_CPP")
        if env_cpp:
            tm_cpp = Path(env_cpp).resolve()
        elif (edockq_root / "TMscore.cpp").is_file():
            tm_cpp = edockq_root / "TMscore.cpp"
        else:
            tm_cpp = bundled_cpp
        tm_dir = tm_exe.parent
        return tm_dir, tm_cpp, tm_exe, False

    local_exe = _executable_candidate(edockq_root)
    local_cpp = edockq_root / "TMscore.cpp"
    if _is_runnable_binary(local_exe):
        cpp = local_cpp if local_cpp.is_file() else bundled_cpp
        return edockq_root, cpp, local_exe, False

    cache = pliq_cache_dir()
    tm_exe = cache / ("TMscore.exe" if is_windows() else "TMscore")
    compile_on_run = True
    if tm_exe.is_file():
        try:
            if tm_exe.stat().st_mtime >= bundled_cpp.stat().st_mtime:
                compile_on_run = False
        except OSError:
            pass
    return cache, bundled_cpp, tm_exe, compile_on_run


def resolve_binana_python() -> Path:
    if VENDOR_BINANA.is_dir():
        return VENDOR_BINANA
    raise FileNotFoundError(f"Bundled BINANA python not found: {VENDOR_BINANA}")


def compile_bundled_tmscore_to_cache() -> Path:
    """Compile vendored ``TMscore.cpp`` into :func:`pliq_cache_dir` (see :func:`pliq.platinfo.tmscore_compile_argv`)."""
    if not VENDOR_TMSCORE_CPP.is_file():
        raise FileNotFoundError(f"Bundled TMscore.cpp missing: {VENDOR_TMSCORE_CPP}")
    cache = pliq_cache_dir()
    tm_exe = cache / ("TMscore.exe" if is_windows() else "TMscore")
    from pliq.tmscore_run.runner import _compile_tmscore

    _compile_tmscore(cache, VENDOR_TMSCORE_CPP, tm_exe)
    return tm_exe
