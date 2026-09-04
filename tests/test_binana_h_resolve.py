# -*- coding: utf-8 -*-
from pliq.edockq.binana_recall import resolve_binana_h_settings


def test_default_is_strip_d():
    ph, strip = resolve_binana_h_settings(None)
    assert ph is None
    assert strip is True


def test_no_binana_obabel_h_alias():
    ph, strip = resolve_binana_h_settings(7.4, force_no_obabel_ph=True)
    assert ph is None
    assert strip is True


def test_obabel_ph_disables_strip():
    ph, strip = resolve_binana_h_settings(7.4)
    assert ph == 7.4
    assert strip is False


def test_keep_h_as_is():
    ph, strip = resolve_binana_h_settings(None, keep_h=True)
    assert ph is None
    assert strip is False
