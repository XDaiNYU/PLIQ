# -*- coding: utf-8 -*-
"""External dependency URLs and install checks for PLIQ."""
from __future__ import annotations

import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional

# --- Canonical upstream sources (documented in README / pyproject.toml) ---

OTMOL_GITHUB = "https://github.com/weixiaoqimath/otmol"
OTMOL_PIP_SPEC = "otmol @ git+https://github.com/weixiaoqimath/otmol.git"
OTMOL_INSTALL_CMD = "pip install git+https://github.com/weixiaoqimath/otmol.git"

# TM-score is NOT on PyPI. PLIQ ships Zhang-lab TMscore.cpp and builds a local binary.
TMSCORE_UPSTREAM = "https://zhanggroup.org/TM-score/"
TMSCORE_BUNDLED_CPP = "pliq/_vendor/tmscore/TMscore.cpp"
TMSCORE_BUILD_CMD = "pliq compile-tmscore"


@dataclass
class DepStatus:
    name: str
    ok: bool
    detail: str
    install_hint: str = ""


def check_otmol() -> DepStatus:
    try:
        import otmol  # noqa: F401

        return DepStatus("OTMol", True, "import otmol OK")
    except ImportError as e:
        return DepStatus(
            "OTMol",
            False,
            str(e)[:200],
            f"{OTMOL_INSTALL_CMD}  # upstream: {OTMOL_GITHUB}",
        )


def check_obabel() -> DepStatus:
    exe = shutil.which("obabel")
    if exe:
        return DepStatus("Open Babel (obabel)", True, exe)
    return DepStatus(
        "Open Babel (obabel)",
        False,
        "not found on PATH",
        "conda install -c conda-forge openbabel  # or brew install open-babel",
    )


def check_tmscore(tm_exe: Optional[Path] = None) -> DepStatus:
    if tm_exe is not None and tm_exe.is_file():
        try:
            subprocess.run(
                [str(tm_exe), "-h"],
                capture_output=True,
                timeout=5,
            )
            return DepStatus("TM-score binary", True, str(tm_exe))
        except Exception as e:
            return DepStatus(
                "TM-score binary",
                False,
                f"{tm_exe} not runnable: {e}",
                TMSCORE_BUILD_CMD,
            )
    return DepStatus(
        "TM-score binary",
        False,
        "not built yet",
        f"{TMSCORE_BUILD_CMD}  # needs g++/c++ in PATH; source bundled in {TMSCORE_BUNDLED_CPP}",
    )


def check_all(tm_exe: Optional[Path] = None) -> List[DepStatus]:
    return [check_otmol(), check_obabel(), check_tmscore(tm_exe)]
