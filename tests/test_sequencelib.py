"""Tests for sequencelib complement and reverse_complement functions."""

from HCRProbeDesign import sequencelib


def test_complement_uppercase():
    assert "".join(sequencelib.complement("ACGT")) == "TGCA"


def test_complement_lowercase():
    assert "".join(sequencelib.complement("acgt")) == "tgca"


def test_complement_mixed_case():
    assert "".join(sequencelib.complement("AcGt")) == "TgCa"


def test_complement_with_n():
    assert "".join(sequencelib.complement("ANn")) == "TNn"


def test_reverse_complement_uppercase():
    assert sequencelib.reverse_complement("ACGT") == "ACGT"  # palindrome
    assert sequencelib.reverse_complement("AACC") == "GGTT"


def test_reverse_complement_lowercase():
    assert sequencelib.reverse_complement("acgt") == "acgt"  # palindrome
    assert sequencelib.reverse_complement("aacc") == "ggtt"


def test_reverse_complement_mixed_case():
    assert sequencelib.reverse_complement("AaCc") == "gGtT"


def test_complement_c_to_g_not_t():
    """Regression test: lowercase 'c' must complement to 'g', not 't'."""
    result = "".join(sequencelib.complement("c"))
    assert result == "g", f"complement('c') returned '{result}', expected 'g'"
    result = "".join(sequencelib.complement("C"))
    assert result == "G", f"complement('C') returned '{result}', expected 'G'"
