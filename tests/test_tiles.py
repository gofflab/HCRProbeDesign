"""Tests for Tile class overlap detection and __len__."""

from HCRProbeDesign.tiles import Tile


def test_len_equals_sequence_length():
    tile = Tile(sequence="ACGT" * 13, seqName="test", startPos=1)  # 52 bases
    assert len(tile) == 52


def test_overlapping_tiles():
    t1 = Tile(sequence="A" * 52, seqName="test", startPos=1)
    t2 = Tile(sequence="A" * 52, seqName="test", startPos=26)  # overlaps with t1
    assert t1.overlaps(t2)
    assert t2.overlaps(t1)


def test_adjacent_tiles_do_not_overlap():
    """Regression test: tiles at adjacent non-overlapping positions must not overlap."""
    t1 = Tile(sequence="A" * 52, seqName="test", startPos=1)   # covers 1..52, end=53
    t2 = Tile(sequence="A" * 52, seqName="test", startPos=53)  # covers 53..104
    assert not t1.overlaps(t2), "Adjacent tiles incorrectly flagged as overlapping"
    assert not t2.overlaps(t1), "Adjacent tiles incorrectly flagged as overlapping"


def test_distant_tiles_do_not_overlap():
    t1 = Tile(sequence="A" * 52, seqName="test", startPos=1)
    t2 = Tile(sequence="A" * 52, seqName="test", startPos=100)
    assert not t1.overlaps(t2)
    assert not t2.overlaps(t1)


def test_same_start_overlaps():
    t1 = Tile(sequence="A" * 52, seqName="test", startPos=10)
    t2 = Tile(sequence="A" * 52, seqName="test", startPos=10)
    assert t1.overlaps(t2)


def test_one_base_overlap():
    t1 = Tile(sequence="A" * 52, seqName="test", startPos=1)   # end=53
    t2 = Tile(sequence="A" * 52, seqName="test", startPos=52)  # starts at 52, within t1
    assert t1.overlaps(t2)
    assert t2.overlaps(t1)
