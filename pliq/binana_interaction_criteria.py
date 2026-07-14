# -*- coding: utf-8 -*-
"""Load BINANA distance/angle cutoffs from JSON or Excel (PLQ hyperparameters).

**Tier columns** ``low`` / ``medium`` / ``high`` (Excel or nested JSON) follow:

- **``low``** — *low-strength* / **loose** detection: **larger** distance caps and
  **larger** angle slack → more atom pairs qualify.
- **``medium``** — default balance (default tier name and spreadsheet column).
- **``high``** — *high-strength* / **strict** detection: **smaller** distance caps and
  **smaller** angle slack → fewer, tighter interactions count.

Legacy: ``medium1`` is accepted as an **alias** for ``medium`` (CLI flag, JSON ``"tier"``,
nested object keys, or a spreadsheet column named ``medium1`` when ``medium`` is absent).

So for distance maxima (Å), typically ``high`` < ``medium`` < ``low``. For angle
tolerances (°), stricter = smaller allowed deviation → same ordering where applicable.

**Hydrogen / halogen** share BINANA's ``hydrogen_halogen_bond_angle_cutoff`` (donor–
(H|halogen)–acceptor); **larger** ° is looser, **smaller** ° is stricter.

Stricter than default: ``--binana-criteria-tier high``. Looser: ``low``.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, Mapping, Optional, Union

# Keys must match ``binana.interactions.get_all_interactions`` keyword names.
KNOWN_PARAMS = frozenset(
    {
        "closest_dist_cutoff",
        "close_dist_cutoff",
        "electrostatic_dist_cutoff",
        "active_site_flexibility_dist_cutoff",
        "hydrophobic_dist_cutoff",
        "hydrogen_bond_dist_cutoff",
        "hydrogen_halogen_bond_angle_cutoff",
        "halogen_bond_dist_cutoff",
        "pi_pi_general_dist_cutoff",
        "pi_stacking_angle_tol",
        "t_stacking_angle_tol",
        "t_stacking_closest_dist_cutoff",
        "cation_pi_dist_cutoff",
        "salt_bridge_dist_cutoff",
        "metal_coordination_dist_cutoff",
        "pi_padding",
    }
)

DEFAULT_BINANA_TIER = "medium"

_JSON_META_KEYS = frozenset({"tier", "description", "comment", "name", "version"})
_JSON_CONTAINER_KEYS = frozenset({"cutoffs", "parameters", "params"})


def _normalize_tier(tier: Optional[str]) -> str:
    t = (tier or DEFAULT_BINANA_TIER).strip().lower()
    if t == "medium1":
        t = "medium"
    if t not in ("low", "medium", "high"):
        raise ValueError(
            f"tier must be low|medium|high (or legacy alias 'medium1'), got {tier!r}"
        )
    return t


def _tier_json_lookup_keys(tier: str) -> tuple[str, ...]:
    """Keys to try on nested tier dicts (Excel-style capitalization + legacy medium1)."""
    t = tier.lower()
    if t == "medium":
        return (
            "medium",
            "Medium",
            "MEDIUM",
            "medium1",
            "Medium1",
            "MEDIUM1",
        )
    return (tier, tier.capitalize(), tier.upper())


def _resolve_excel_tier_column(df: Any, tier: str) -> Any:
    """Return the actual DataFrame column label for *tier* (case-insensitive; legacy *medium1*)."""
    tier = _normalize_tier(tier)
    lower_map = {str(c).strip().lower(): c for c in df.columns}
    tl = tier.lower()
    if tl in lower_map:
        return lower_map[tl]
    if tl == "medium" and "medium1" in lower_map:
        return lower_map["medium1"]
    raise ValueError(f"{tier!r}: no matching tier column; columns={list(df.columns)}")


def _float_if_number(x: Any) -> Optional[float]:
    if isinstance(x, bool):
        return None
    if isinstance(x, (int, float)):
        return float(x)
    return None


def _coerce_cutoff_value(raw: Any, *, tier: str) -> Optional[float]:
    """Scalar number, or dict with ``low`` / ``medium`` / ``high`` (``medium1`` = legacy alias)."""
    v = _float_if_number(raw)
    if v is not None:
        return v
    if isinstance(raw, Mapping):
        for tk in _tier_json_lookup_keys(tier):
            if tk in raw:
                return _float_if_number(raw[tk])
    return None


def _merge_criteria_mapping(obj: Mapping[str, Any]) -> Dict[str, Any]:
    """Flatten optional ``cutoffs`` / ``parameters`` nesting; skip metadata keys."""
    flat: Dict[str, Any] = {}
    for key, val in obj.items():
        if not isinstance(key, str):
            continue
        lk = key.strip()
        if not lk or lk.startswith("_"):
            continue
        if lk.lower() in _JSON_META_KEYS:
            continue
        if lk in _JSON_CONTAINER_KEYS and isinstance(val, Mapping):
            flat.update(_merge_criteria_mapping(val))
        else:
            flat[lk] = val
    return flat


def load_binana_get_all_interaction_kwargs_from_json(
    json_path: Path,
    *,
    tier: str = DEFAULT_BINANA_TIER,
) -> Dict[str, Any]:
    """
    Read cutoffs from JSON.

    Supported shapes:

    1. **Flat** — one number per BINANA kwarg (``tier`` is ignored)::

           {"hydrogen_bond_dist_cutoff": 3.6, "pi_pi_general_dist_cutoff": 7.0, ...}

    2. **Nested tiers** — each value is a number or an object with ``low`` / ``medium`` / ``high``
       (``low`` = loose / larger caps, ``high`` = strict / smaller caps). Legacy key ``medium1``::

           {"hydrogen_bond_dist_cutoff": {"low": 4.0, "medium": 3.6, "high": 3.2}}

    3. **Wrapper** — optional ``cutoffs`` / ``parameters`` object with the same keys inside.

    Unknown keys are ignored; values that are not numbers (after coercion) are skipped.
    """
    tier = _normalize_tier(tier)
    json_path = Path(json_path).resolve()
    if not json_path.is_file():
        raise FileNotFoundError(f"interaction criteria JSON not found: {json_path}")

    with open(json_path, encoding="utf-8") as f:
        root = json.load(f)

    if not isinstance(root, dict):
        raise ValueError(f"{json_path}: JSON root must be an object, got {type(root).__name__}")

    file_tier = root.get("tier") or root.get("binana_criteria_tier")
    if isinstance(file_tier, str) and file_tier.strip():
        tier = _normalize_tier(file_tier)

    merged = _merge_criteria_mapping(root)
    out: Dict[str, Any] = {}
    for key, raw in merged.items():
        if key not in KNOWN_PARAMS:
            continue
        val = _coerce_cutoff_value(raw, tier=tier)
        if val is None:
            continue
        out[key] = val
    return out


def load_binana_get_all_interaction_kwargs(
    xlsx_path: Path,
    *,
    tier: str = DEFAULT_BINANA_TIER,
) -> Dict[str, Any]:
    """
    Read one sheet with columns: ``parameter``, ``low``, ``medium``, ``high`` (numeric).

    Legacy workbooks may use ``medium1`` instead of ``medium`` for the middle column; that
    column is used when ``tier`` is ``medium`` and ``medium`` is missing.

    Tier semantics: ``low`` = loose (larger thresholds), ``high`` = strict (smaller).
    Returns kwargs for ``get_all_interactions`` using the chosen tier column.
    Unknown ``parameter`` rows are skipped; empty cells are skipped (BINANA defaults apply).
    """
    tier = _normalize_tier(tier)
    xlsx_path = Path(xlsx_path).resolve()
    if not xlsx_path.is_file():
        raise FileNotFoundError(f"interaction_criteria.xlsx not found: {xlsx_path}")

    try:
        import pandas as pd
    except ImportError as e:
        raise ImportError("Reading interaction_criteria.xlsx requires pandas") from e

    try:
        df = pd.read_excel(xlsx_path, engine="openpyxl")
    except ImportError as e:
        raise ImportError(
            "Reading .xlsx requires openpyxl; install with: pip install openpyxl"
        ) from e

    col_param = None
    for c in df.columns:
        if str(c).strip().lower() in ("parameter", "param", "binana_kwarg", "binana_parameter"):
            col_param = c
            break
    if col_param is None:
        raise ValueError(
            f"{xlsx_path}: expected a column named 'parameter' (or param/binana_kwarg); got {list(df.columns)}"
        )

    tier_col = _resolve_excel_tier_column(df, tier)

    out: Dict[str, Any] = {}
    for _, row in df.iterrows():
        key = str(row[col_param]).strip()
        if not key or key == "nan":
            continue
        if key not in KNOWN_PARAMS:
            continue
        val = row[tier_col]
        if val is None or (isinstance(val, float) and val != val):  # NaN
            continue
        try:
            out[key] = float(val)
        except (TypeError, ValueError):
            continue
    return out


def load_binana_get_all_interaction_kwargs_from_path(
    path: Union[str, Path],
    *,
    tier: str = DEFAULT_BINANA_TIER,
) -> Dict[str, Any]:
    """
    Load BINANA ``get_all_interactions`` kwargs from ``.json`` or Excel (``.xlsx`` / ``.xls``).

    JSON is parsed with :func:`load_binana_get_all_interaction_kwargs_from_json`;
    Excel with :func:`load_binana_get_all_interaction_kwargs`.
    """
    path = Path(path).resolve()
    if not path.is_file():
        raise FileNotFoundError(f"interaction criteria file not found: {path}")
    suf = path.suffix.lower()
    if suf == ".json":
        return load_binana_get_all_interaction_kwargs_from_json(path, tier=tier)
    if suf in (".xlsx", ".xls"):
        return load_binana_get_all_interaction_kwargs(path, tier=tier)
    raise ValueError(
        f"Unsupported interaction criteria file type {path.suffix!r}; use .json, .xlsx, or .xls"
    )
