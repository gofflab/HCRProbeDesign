import os

import pytest

from HCRProbeDesign import referenceGenome as rg


def _write_fasta(path, name="chr1", seq="ACGT"):
    path.write_text(f">{name}\n{seq}\n")


def test_collect_fasta_inputs_dir_and_file(tmp_path):
    fasta_dir = tmp_path / "genomes"
    fasta_dir.mkdir()
    fa = fasta_dir / "a.fa"
    fasta = fasta_dir / "b.fasta"
    extra = tmp_path / "extra.fa"
    _write_fasta(fa)
    _write_fasta(fasta)
    _write_fasta(extra)

    result = rg.collect_fasta_inputs([str(fasta_dir), str(extra), str(fa)])

    assert result == [str(fa), str(fasta), str(extra)]


def test_collect_fasta_inputs_missing_path(tmp_path):
    missing = tmp_path / "missing.fa"
    with pytest.raises(FileNotFoundError):
        rg.collect_fasta_inputs([str(missing)])


def test_collect_fasta_inputs_empty_dir(tmp_path):
    empty_dir = tmp_path / "empty"
    empty_dir.mkdir()
    with pytest.raises(FileNotFoundError):
        rg.collect_fasta_inputs([str(empty_dir)])


def test_format_index_path_relative(monkeypatch, tmp_path):
    monkeypatch.setattr(rg, "PACKAGE_DIRECTORY", str(tmp_path))
    index_prefix = tmp_path / "indices" / "zfish" / "zfish"
    expected = os.path.relpath(str(index_prefix), str(tmp_path))
    assert rg.format_index_path(str(index_prefix)) == expected


def test_format_index_path_absolute(monkeypatch, tmp_path):
    monkeypatch.setattr(rg, "PACKAGE_DIRECTORY", str(tmp_path))
    outside = tmp_path.parent / "other" / "idx"
    outside.parent.mkdir(parents=True, exist_ok=True)
    assert rg.format_index_path(str(outside)) == str(outside)


def test_build_bowtie2_index_requires_build_tool(monkeypatch):
    monkeypatch.setattr(rg.shutil, "which", lambda _: None)
    with pytest.raises(RuntimeError):
        rg.build_bowtie2_index(["genome.fa"], "mouse")


def test_build_bowtie2_index_calls_bowtie2_build(monkeypatch, tmp_path):
    fa = tmp_path / "genome.fa"
    _write_fasta(fa)

    monkeypatch.setattr(rg.shutil, "which", lambda _: "/usr/bin/bowtie2-build")
    captured = {}

    def fake_check_call(cmd):
        captured["cmd"] = cmd
        return 0

    monkeypatch.setattr(rg.subprocess, "check_call", fake_check_call)

    prefix = rg.build_bowtie2_index(
        [str(fa)],
        "zfish",
        index_name="danio",
        indices_dir=str(tmp_path),
        threads=4,
    )

    expected_prefix = tmp_path / "zfish" / "danio"
    assert prefix == str(expected_prefix)
    assert captured["cmd"] == [
        "/usr/bin/bowtie2-build",
        "--threads",
        "4",
        str(fa),
        str(expected_prefix),
    ]


def test_build_bowtie2_index_force_overwrite(monkeypatch, tmp_path):
    fa = tmp_path / "genome.fa"
    _write_fasta(fa)
    species_dir = tmp_path / "mouse"
    species_dir.mkdir()
    index_prefix = species_dir / "mouse"
    existing = index_prefix.with_suffix(".1.bt2")
    existing.write_text("old index")

    monkeypatch.setattr(rg.shutil, "which", lambda _: "/usr/bin/bowtie2-build")
    monkeypatch.setattr(rg.subprocess, "check_call", lambda cmd: 0)

    with pytest.raises(FileExistsError):
        rg.build_bowtie2_index([str(fa)], "mouse", indices_dir=str(tmp_path))

    rg.build_bowtie2_index(
        [str(fa)],
        "mouse",
        indices_dir=str(tmp_path),
        force=True,
    )
    assert not existing.exists()


def test_register_species_adds_entry(tmp_path):
    config_path = tmp_path / "HCRconfig.yaml"
    config_path.write_text(
        "species:\n  mouse:\n    bowtie2_index: indices/mm10/mm10\n"
        "default_params:\n  tileSize: 52\n"
    )

    rg.register_species(str(config_path), "zfish", "/abs/index")
    config = rg.load_config(str(config_path))
    assert config["species"]["zfish"]["bowtie2_index"] == "/abs/index"
    assert config["default_params"]["tileSize"] == 52


def test_register_species_requires_force(tmp_path):
    config_path = tmp_path / "HCRconfig.yaml"
    config_path.write_text("species:\n  mouse:\n    bowtie2_index: indices/mm10/mm10\n")

    with pytest.raises(ValueError):
        rg.register_species(str(config_path), "mouse", "indices/mm10/mm10")
