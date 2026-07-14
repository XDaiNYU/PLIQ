# -*- coding: utf-8 -*-
"""Host OS detection for PLIQ (Linux vs macOS vs Windows). Used for TM-score compile, cache paths, etc."""
from __future__ import annotations

import os
import sys
from pathlib import Path
from typing import List


def getenv_compat(suffix: str, default: str = "") -> str:
    """Read ``PLIQ_<suffix>`` (preferred) then legacy ``PLQ_<suffix>``; else ``default``."""
    v = os.environ.get("PLIQ_" + suffix, "").strip()
    if v:
        return v
    v = os.environ.get("PLQ_" + suffix, "").strip()
    return v if v else default


def pliq_os_kind() -> str:
    """
    Return ``"darwin"``, ``"linux"``, ``"windows"``, or ``"other"``.

    Uses ``sys.platform`` (not uname) so behavior matches CPython on the running machine.
    """
    if sys.platform == "darwin":
        return "darwin"
    if sys.platform.startswith("linux"):
        return "linux"
    if sys.platform == "win32":
        return "windows"
    return "other"


def is_mac() -> bool:
    return pliq_os_kind() == "darwin"


def is_windows() -> bool:
    return pliq_os_kind() == "windows"


def tmscore_compile_argv(tm_cpp: Path, tm_exe: Path) -> List[str]:
    """
    Compiler argv for bundled ``TMscore.cpp`` → ``tm_exe``.

    - **Linux**: ``g++ -O3 -ffast-math -o exe cpp -lm``
    - **macOS**: ``c++ -O3 -ffast-math -stdlib=libc++ -o exe cpp``
    - **Windows**: ``g++ -O3 -ffast-math -o exe cpp`` (MinGW / MSYS2)

    **Overrides**
    - ``PLIQ_CXX`` (legacy ``PLQ_CXX``): compiler executable (e.g. ``g++-12``, ``clang++``).
      PLIQ then appends ``-O3 -ffast-math`` and ``-o <exe> <cpp>`` only (no auto ``-lm``).
    - ``PLIQ_CXXFLAGS`` (legacy ``PLQ_CXXFLAGS``): extra flags (space-separated), inserted before ``-o``.
    """
    tm_cpp = Path(tm_cpp)
    tm_exe = Path(tm_exe)
    cxx = getenv_compat("CXX")
    xflags = [f for f in getenv_compat("CXXFLAGS").split() if f]

    mid = ["-O3", "-ffast-math"] + xflags
    end = ["-o", str(tm_exe), str(tm_cpp)]

    if cxx:
        return cxx.split() + mid + end

    kind = pliq_os_kind()
    if kind == "darwin":
        return ["c++"] + mid + ["-stdlib=libc++"] + end
    if kind == "linux":
        return ["g++"] + mid + end + ["-lm"]
    if kind == "windows":
        return ["g++"] + mid + end
    return ["g++"] + mid + end + ["-lm"]


def describe_platform_line() -> str:
    return f"PLIQ host: {pliq_os_kind()} (sys.platform={sys.platform!r})"
