import os

import pytest

from HCRProbeDesign import listReferences as lr
from HCRProbeDesign import _datadir


def _write_config(path, content):
    path.write_text(content)


def test_load_config_missing_file(tmp_path):
    missing = tmp_path / "missing.yaml"
    assert lr.load_config(str(missing)) == {}


def test_load_config_valid(tmp_path):
    config_path = tmp_path / "HCRconfig.yaml"
    _write_config(config_path, "species:\n  mouse:\n    bowtie2_index: indices/mm10/mm10\n")
    config = lr.load_config(str(config_path))
    assert "mouse" in config["species"]
    assert config["species"]["mouse"]["bowtie2_index"] == "indices/mm10/mm10"


def test_list_references_no_species(tmp_path):
    config_path = tmp_path / "HCRconfig.yaml"
    _write_config(config_path, "species: {}\ndefault_params:\n  tileSize: 52\n")
    info = lr.list_references(str(config_path))
    assert info["species"] == []
    assert info["default_params"]["tileSize"] == 52


def test_list_references_with_species(tmp_path):
    # Create a fake index file so we can verify the installed check
    species_dir = tmp_path / "indices" / "zfish"
    species_dir.mkdir(parents=True)
    index_prefix = species_dir / "zfish"
    (species_dir / "zfish.1.bt2").write_text("fake")

    config_path = tmp_path / "HCRconfig.yaml"
    _write_config(
        config_path,
        f"species:\n  zebrafish:\n    bowtie2_index: {index_prefix}\n"
        f"default_params:\n  tileSize: 52\n",
    )

    info = lr.list_references(str(config_path))
    assert len(info["species"]) == 1
    sp = info["species"][0]
    assert sp["name"] == "zebrafish"
    assert sp["installed"] is True
    assert len(sp["index_files"]) == 1


def test_list_references_missing_index(tmp_path):
    config_path = tmp_path / "HCRconfig.yaml"
    _write_config(
        config_path,
        "species:\n  human:\n    bowtie2_index: /nonexistent/path/hg38\n"
        "default_params: {}\n",
    )
    info = lr.list_references(str(config_path))
    sp = info["species"][0]
    assert sp["name"] == "human"
    assert sp["installed"] is False
    assert sp["index_files"] == []


def test_format_references_no_species(tmp_path):
    config_path = tmp_path / "HCRconfig.yaml"
    _write_config(config_path, "species: {}\ndefault_params:\n  tileSize: 52\n")
    info = lr.list_references(str(config_path))
    output = lr.format_references(info)
    assert "No reference genomes registered" in output
    assert "buildGenomeIndex" in output
    assert "tileSize" in output


def test_format_references_with_species(tmp_path):
    species_dir = tmp_path / "indices" / "mm10"
    species_dir.mkdir(parents=True)
    index_prefix = species_dir / "mm10"
    for suffix in ["1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2"]:
        (species_dir / f"mm10.{suffix}").write_text("fake")

    config_path = tmp_path / "HCRconfig.yaml"
    _write_config(
        config_path,
        f"species:\n  mouse:\n    bowtie2_index: {index_prefix}\n"
        f"default_params:\n  tileSize: 52\n  minGC: 45.0\n",
    )
    info = lr.list_references(str(config_path))
    output = lr.format_references(info)
    assert "mouse" in output
    assert "INSTALLED" in output
    assert "6 files" in output
    assert "designProbes --species mouse" in output
    assert "tileSize" in output
    assert "--tileSize" in output
    assert "minGC" in output


def test_format_references_missing_index_status(tmp_path):
    config_path = tmp_path / "HCRconfig.yaml"
    _write_config(
        config_path,
        "species:\n  human:\n    bowtie2_index: /nonexistent/hg38\ndefault_params: {}\n",
    )
    info = lr.list_references(str(config_path))
    output = lr.format_references(info)
    assert "MISSING" in output
    assert "none found" in output


def test_resolve_absolute_index_path_absolute():
    result = lr._resolve_absolute_index_path("/abs/path/index")
    assert result == "/abs/path/index"


def test_resolve_absolute_index_path_relative(monkeypatch, tmp_path):
    data_dir = str(tmp_path / "datadir")
    monkeypatch.setenv("HCRPROBEDESIGN_DATA_DIR", data_dir)
    result = lr._resolve_absolute_index_path("indices/mm10/mm10")
    expected = os.path.join(data_dir, "indices/mm10/mm10")
    assert result == expected


def test_check_index_files_empty(tmp_path):
    prefix = tmp_path / "noindex" / "prefix"
    assert lr._check_index_files(str(prefix)) == []


def test_check_index_files_found(tmp_path):
    prefix = tmp_path / "idx"
    (tmp_path / "idx.1.bt2").write_text("data")
    (tmp_path / "idx.2.bt2").write_text("data")
    result = lr._check_index_files(str(prefix))
    assert len(result) == 2
