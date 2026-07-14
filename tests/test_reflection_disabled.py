# -*- coding: utf-8 -*-
"""PLIQ v6: reflection must never be enabled in OTMol alignment presets."""
from __future__ import annotations

import inspect

from pliq.edockq import ligand_mapping as LM


def test_align_presets_never_use_reflection():
    src = inspect.getsource(LM.align_ligands_otmol_bci_constrained)
    assert '"reflection": True' not in src
    assert "reflection=True" not in src
    assert "reflection=False" in src
