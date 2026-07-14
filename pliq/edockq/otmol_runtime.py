# -*- coding: utf-8 -*-
"""Ensure ``import otmol`` works (pip install or local ``otmol_package``)."""
from __future__ import annotations

import sys
from pathlib import Path
from typing import Optional


def _find_otmol_repo_root() -> Optional[Path]:
    from pliq.platinfo import getenv_compat

    env = getenv_compat("EDOCKQ_ROOT")
    if env:
        p = Path(env).resolve()
        if (p / "otmol_package" / "otmol").is_dir():
            return p
    here = Path(__file__).resolve()
    for p in here.parents:
        if (p / "otmol_package" / "otmol").is_dir():
            return p
    return None


def ensure_local_otmol_on_path(repo_root: Optional[Path] = None) -> Path:
    """
    Prefer a pip-installed ``otmol`` (https://github.com/weixiaoqimath/otmol).
    If unavailable, prepend ``{repo_root}/otmol_package`` to ``sys.path`` when present.

    Returns the repo root used for the local fallback, or cwd if nothing found.
    """
    try:
        import otmol  # noqa: F401
    except ImportError:
        pass
    else:
        return Path(repo_root).resolve() if repo_root is not None else Path.cwd()

    if repo_root is not None:
        repo_root = Path(repo_root).resolve()
    else:
        found = _find_otmol_repo_root()
        repo_root = found if found is not None else Path.cwd()

    pkg = repo_root / "otmol_package"
    if pkg.is_dir():
        s = str(pkg)
        if s not in sys.path:
            sys.path.insert(0, s)
    return repo_root
