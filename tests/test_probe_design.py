import argparse
import sys

import pytest

from HCRProbeDesign import HCR
from HCRProbeDesign import probeDesign
from HCRProbeDesign import tiles


class _FakeHairpin:
    def __init__(self, tm=0.0, structure_found=False):
        self.tm = tm
        self.structure_found = structure_found


def test_fasta_input_creates_hcr_probes(monkeypatch, tmp_path, capsys):
    fasta_path = tmp_path / "input.fa"
    fasta_path.write_text(">target1\nACGTACGTAC\n")

    def fake_calc_hairpin(_seq):
        return _FakeHairpin()

    def fake_calc_tm(_seq):
        return 50.0

    monkeypatch.setattr(probeDesign.primer3, "calc_hairpin", fake_calc_hairpin)
    monkeypatch.setattr(probeDesign.primer3, "calc_tm", fake_calc_tm)
    monkeypatch.setattr(tiles.primer3, "calc_tm", fake_calc_tm)
    monkeypatch.setattr(probeDesign, "outputRunParams", lambda _args: None)

    argv = [
        "probeDesign",
        str(fasta_path),
        "-g",
        "--tileSize",
        "10",
        "--minGC",
        "0",
        "--maxGC",
        "100",
        "--minGibbs",
        "-1000",
        "--maxGibbs",
        "1000",
        "--targetGibbs",
        "0",
        "--maxRunLength",
        "999",
        "--maxProbes",
        "1",
        "--channel",
        "B1",
    ]
    monkeypatch.setattr(sys, "argv", argv)

    probeDesign.main()

    out = capsys.readouterr().out.strip().splitlines()
    assert out[0].startswith("name\tprobe\tstart\tlength\tP1\tP2\tchannel")
    assert len(out) == 2

    row = out[1].split("\t")
    assert row[0].startswith("target1:")
    assert row[1] == "gtacgtacgt"
    assert row[4].startswith(HCR.initiators["B1"]["odd"])
    assert row[5].endswith(HCR.initiators["B1"]["even"])
    assert row[6] == "B1"


def test_batch_fasta_input_creates_hcr_probes(monkeypatch, tmp_path, capsys):
    fasta_path = tmp_path / "input.fa"
    fasta_path.write_text(">target1\nACGTACGTAC\n>target2\nTGCATGCATG\n")

    def fake_calc_hairpin(_seq):
        return _FakeHairpin()

    def fake_calc_tm(_seq):
        return 50.0

    monkeypatch.setattr(probeDesign.primer3, "calc_hairpin", fake_calc_hairpin)
    monkeypatch.setattr(probeDesign.primer3, "calc_tm", fake_calc_tm)
    monkeypatch.setattr(tiles.primer3, "calc_tm", fake_calc_tm)
    monkeypatch.setattr(probeDesign, "outputRunParams", lambda _args: None)

    argv = [
        "probeDesignBatch",
        str(fasta_path),
        "-g",
        "--tileSize",
        "10",
        "--minGC",
        "0",
        "--maxGC",
        "100",
        "--minGibbs",
        "-1000",
        "--maxGibbs",
        "1000",
        "--targetGibbs",
        "0",
        "--maxRunLength",
        "999",
        "--maxProbes",
        "1",
        "--channel",
        "B1",
    ]
    monkeypatch.setattr(sys, "argv", argv)

    probeDesign.main_batch()

    out = capsys.readouterr().out.strip().splitlines()
    assert out[0].startswith("name\tprobe\tstart\tlength\tP1\tP2\tchannel")
    assert len(out) == 3

    names = [row.split("\t")[0] for row in out[1:]]
    assert any(name.startswith("target1:") for name in names)
    assert any(name.startswith("target2:") for name in names)


def test_batch_fasta_header_channel_overrides(monkeypatch, tmp_path, capsys):
    fasta_path = tmp_path / "input.fa"
    fasta_path.write_text(
        ">target1 channel=B1\nACGTACGTAC\n>target2 channel=B2\nTGCATGCATG\n"
    )

    def fake_calc_hairpin(_seq):
        return _FakeHairpin()

    def fake_calc_tm(_seq):
        return 50.0

    monkeypatch.setattr(probeDesign.primer3, "calc_hairpin", fake_calc_hairpin)
    monkeypatch.setattr(probeDesign.primer3, "calc_tm", fake_calc_tm)
    monkeypatch.setattr(tiles.primer3, "calc_tm", fake_calc_tm)
    monkeypatch.setattr(probeDesign, "outputRunParams", lambda _args: None)

    argv = [
        "probeDesignBatch",
        str(fasta_path),
        "-g",
        "--tileSize",
        "10",
        "--minGC",
        "0",
        "--maxGC",
        "100",
        "--minGibbs",
        "-1000",
        "--maxGibbs",
        "1000",
        "--targetGibbs",
        "0",
        "--maxRunLength",
        "999",
        "--maxProbes",
        "1",
        "--channel",
        "B3",
    ]
    monkeypatch.setattr(sys, "argv", argv)

    probeDesign.main_batch()

    out = capsys.readouterr().out.strip().splitlines()
    rows = [row.split("\t") for row in out[1:]]
    for row in rows:
        if row[0].startswith("target1:"):
            assert row[4].startswith(HCR.initiators["B1"]["odd"])
            assert row[5].endswith(HCR.initiators["B1"]["even"])
            assert row[6] == "B1"
        if row[0].startswith("target2:"):
            assert row[4].startswith(HCR.initiators["B2"]["odd"])
            assert row[5].endswith(HCR.initiators["B2"]["even"])
            assert row[6] == "B2"


def test_batch_fasta_header_channel_invalid(monkeypatch, tmp_path):
    fasta_path = tmp_path / "input.fa"
    fasta_path.write_text(">target1 channel=Z9\nACGTACGTAC\n")

    def fake_calc_hairpin(_seq):
        return _FakeHairpin()

    monkeypatch.setattr(probeDesign.primer3, "calc_hairpin", fake_calc_hairpin)
    monkeypatch.setattr(probeDesign, "outputRunParams", lambda _args: None)

    argv = [
        "probeDesignBatch",
        str(fasta_path),
        "-g",
        "--tileSize",
        "10",
        "--minGC",
        "0",
        "--maxGC",
        "100",
        "--minGibbs",
        "-1000",
        "--maxGibbs",
        "1000",
        "--targetGibbs",
        "0",
        "--maxRunLength",
        "999",
        "--maxProbes",
        "1",
    ]
    monkeypatch.setattr(sys, "argv", argv)

    with pytest.raises(ValueError, match="Channel 'Z9' is not defined"):
        probeDesign.main_batch()


def _make_assert_args(species="mouse", no_genomemask=True, index=None):
    return argparse.Namespace(species=species, no_genomemask=no_genomemask, index=index)


def test_assert_species_config_reads_user_data_dir(tmp_path, monkeypatch):
    data_dir = tmp_path / ".hcrprobedesign"
    monkeypatch.setenv("HCRPROBEDESIGN_DATA_DIR", str(data_dir))
    data_dir.mkdir()
    config_path = data_dir / "HCRconfig.yaml"
    config_path.write_text("species:\n  mouse:\n    bowtie2_index: indices/mm10/mm10\n")

    probeDesign._assert_species_config(_make_assert_args(species="mouse"))


def test_assert_species_config_missing_species_raises(tmp_path, monkeypatch):
    data_dir = tmp_path / ".hcrprobedesign"
    monkeypatch.setenv("HCRPROBEDESIGN_DATA_DIR", str(data_dir))
    data_dir.mkdir()
    (data_dir / "HCRconfig.yaml").write_text("species: {}\n")

    with pytest.raises(SystemExit) as excinfo:
        probeDesign._assert_species_config(_make_assert_args(species="mouse"))
    assert "is not registered" in str(excinfo.value)


def test_assert_species_config_skips_when_index_supplied(tmp_path, monkeypatch):
    data_dir = tmp_path / ".hcrprobedesign"
    monkeypatch.setenv("HCRPROBEDESIGN_DATA_DIR", str(data_dir))
    data_dir.mkdir()
    (data_dir / "HCRconfig.yaml").write_text("species: {}\n")

    probeDesign._assert_species_config(
        _make_assert_args(species="mouse", index="/some/index/prefix")
    )


def test_assert_species_config_ignores_stale_package_config(tmp_path, monkeypatch):
    """Regression test: assertion must not read the package-relative config."""
    data_dir = tmp_path / ".hcrprobedesign"
    monkeypatch.setenv("HCRPROBEDESIGN_DATA_DIR", str(data_dir))
    data_dir.mkdir()
    (data_dir / "HCRconfig.yaml").write_text(
        "species:\n  mouse:\n    bowtie2_index: indices/mm10/mm10\n"
    )

    # Simulate a stale package-level config that lacks the registered species.
    stale = tmp_path / "stale_HCRconfig.yaml"
    stale.write_text("species: {}\n")
    monkeypatch.setattr(probeDesign, "package_directory", str(tmp_path))

    probeDesign._assert_species_config(_make_assert_args(species="mouse"))


def test_batch_fasta_header_channel_separator(monkeypatch, tmp_path, capsys):
    fasta_path = tmp_path / "input.fa"
    fasta_path.write_text(">target1|channel=B2\nACGTACGTAC\n")

    def fake_calc_hairpin(_seq):
        return _FakeHairpin()

    def fake_calc_tm(_seq):
        return 50.0

    monkeypatch.setattr(probeDesign.primer3, "calc_hairpin", fake_calc_hairpin)
    monkeypatch.setattr(probeDesign.primer3, "calc_tm", fake_calc_tm)
    monkeypatch.setattr(tiles.primer3, "calc_tm", fake_calc_tm)
    monkeypatch.setattr(probeDesign, "outputRunParams", lambda _args: None)

    argv = [
        "probeDesignBatch",
        str(fasta_path),
        "-g",
        "--tileSize",
        "10",
        "--minGC",
        "0",
        "--maxGC",
        "100",
        "--minGibbs",
        "-1000",
        "--maxGibbs",
        "1000",
        "--targetGibbs",
        "0",
        "--maxRunLength",
        "999",
        "--maxProbes",
        "1",
    ]
    monkeypatch.setattr(sys, "argv", argv)

    probeDesign.main_batch()

    out = capsys.readouterr().out.strip().splitlines()
    row = out[1].split("\t")
    assert row[0].startswith("target1:")
    assert row[4].startswith(HCR.initiators["B2"]["odd"])
    assert row[5].endswith(HCR.initiators["B2"]["even"])
    assert row[6] == "B2"
