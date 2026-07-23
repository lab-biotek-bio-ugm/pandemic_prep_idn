"""Tests for M11 spectral dye-channel separation and shared IUPAC utilities.

Guards the dye-channel fix (previously a hardcoded "PASSED" string) and the
export active-target filter that stopped excluded viruses leaking into reports.
"""
import pytest

from conftest import load_module  # noqa: E402

iupac = load_module("scripts/tools/iupac_utils.py", "iupac_utils")


def test_distinct_channels_pass():
    status, ok = iupac.check_dye_separation("FAM", "HEX")
    assert ok and "PASSED" in status


def test_same_dye_collides():
    status, ok = iupac.check_dye_separation("FAM", "FAM")
    assert not ok and "FAIL" in status


def test_different_dyes_same_channel_collide():
    # HEX and JOE both live in the 'yellow' channel -> indistinguishable.
    status, ok = iupac.check_dye_separation("FAM", "HEX", "JOE")
    assert not ok


def test_triplex_distinct_channels_pass():
    status, ok = iupac.check_dye_separation("FAM", "HEX", "Cy5")
    assert ok


def test_triplex_repeated_channel_collides():
    status, ok = iupac.check_dye_separation("FAM", "HEX", "HEX")
    assert not ok


def test_unknown_dye_is_not_silently_passed():
    status, ok = iupac.check_dye_separation("FAM", "NotARealDye")
    assert not ok and "UNKNOWN" in status


def test_dye_names_are_case_and_separator_insensitive():
    # "cy-5", "CY5", "Cy 5" all normalise to the same far-red channel.
    _, ok1 = iupac.check_dye_separation("FAM", "cy-5")
    _, ok2 = iupac.check_dye_separation("FAM", "CY5")
    assert ok1 and ok2


def test_export_multiplex_active_filter():
    es = load_module("scripts/export_supplementary.py", "export_supplementary")
    es.ACTIVE_TARGETS = ["influenza-a", "influenza-b"]
    assert es.is_multiplex_active("multiplex-influenza-a_x_influenza-b") is True
    # A triplex containing Measles must be excluded when Measles is not active.
    assert es.is_multiplex_active("multiplex-triplex-influenza-a_influenza-b_measles") is False
    assert es.is_multiplex_active("multiplex-influenza-a_x_measles") is False
