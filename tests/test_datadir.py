import os
import glob

import pytest
import yaml

from HCRProbeDesign import _datadir


@pytest.fixture(autouse=True)
def isolate_data_dir(tmp_path, monkeypatch):
    """Point all tests at a temporary data directory."""
    monkeypatch.setenv("HCRPROBEDESIGN_DATA_DIR", str(tmp_path / "data"))


def test_get_data_dir_from_env(monkeypatch):
    monkeypatch.setenv("HCRPROBEDESIGN_DATA_DIR", "/custom/path")
    assert _datadir.get_data_dir() == "/custom/path"


def test_get_data_dir_default(monkeypatch):
    monkeypatch.delenv("HCRPROBEDESIGN_DATA_DIR", raising=False)
    expected = os.path.join(os.path.expanduser("~"), ".hcrprobedesign")
    assert _datadir.get_data_dir() == expected


def test_get_config_path():
    data_dir = _datadir.get_data_dir()
    assert _datadir.get_config_path() == os.path.join(data_dir, "HCRconfig.yaml")


def test_get_indices_dir():
    data_dir = _datadir.get_data_dir()
    assert _datadir.get_indices_dir() == os.path.join(data_dir, "indices")


def test_ensure_data_dir_creates_directory(tmp_path, monkeypatch):
    data_dir = tmp_path / "data"
    monkeypatch.setenv("HCRPROBEDESIGN_DATA_DIR", str(data_dir))

    assert not data_dir.exists()
    _datadir.ensure_data_dir()

    assert data_dir.exists()
    assert (data_dir / "HCRconfig.yaml").exists()
    assert (data_dir / "indices").is_dir()


def test_ensure_data_dir_seeds_config(tmp_path, monkeypatch):
    data_dir = tmp_path / "data"
    monkeypatch.setenv("HCRPROBEDESIGN_DATA_DIR", str(data_dir))

    _datadir.ensure_data_dir()

    config = yaml.safe_load((data_dir / "HCRconfig.yaml").read_text())
    assert "default_params" in config
    assert config["default_params"]["tileSize"] == 52


def test_ensure_data_dir_idempotent(tmp_path, monkeypatch):
    data_dir = tmp_path / "data"
    monkeypatch.setenv("HCRPROBEDESIGN_DATA_DIR", str(data_dir))

    _datadir.ensure_data_dir()
    # Modify config to verify it doesn't get overwritten
    config_path = data_dir / "HCRconfig.yaml"
    config = yaml.safe_load(config_path.read_text())
    config["species"]["test"] = {"bowtie2_index": "/test/index"}
    with open(config_path, "w") as fh:
        yaml.safe_dump(config, fh)

    _datadir.ensure_data_dir()

    config2 = yaml.safe_load(config_path.read_text())
    assert "test" in config2["species"]


def test_migrate_old_data(tmp_path, monkeypatch):
    data_dir = tmp_path / "data"
    monkeypatch.setenv("HCRPROBEDESIGN_DATA_DIR", str(data_dir))

    # Set up old package directory with species and index
    old_indices = tmp_path / "old_pkg" / "indices" / "mm10"
    old_indices.mkdir(parents=True)
    for suffix in ["1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2"]:
        (old_indices / f"mm10.{suffix}").write_text("fake")

    old_config_path = tmp_path / "old_pkg" / "HCRconfig.yaml"
    old_config = {
        "species": {"mouse": {"bowtie2_index": "indices/mm10/mm10"}},
        "default_params": {"tileSize": 52},
    }
    with open(old_config_path, "w") as fh:
        yaml.safe_dump(old_config, fh)

    # Point _datadir at the fake old package dir
    monkeypatch.setattr(_datadir, "_PACKAGE_DIRECTORY", str(tmp_path / "old_pkg"))
    monkeypatch.setattr(
        _datadir, "_SEED_CONFIG", str(old_config_path)
    )

    _datadir.ensure_data_dir()

    # Verify migration
    new_config = yaml.safe_load((data_dir / "HCRconfig.yaml").read_text())
    assert "mouse" in new_config["species"]

    new_index_dir = data_dir / "indices" / "mm10"
    assert new_index_dir.exists()
    bt2_files = glob.glob(str(new_index_dir / "*.bt2*"))
    assert len(bt2_files) == 6


def test_migrate_auto_registers_species_when_config_empty(tmp_path, monkeypatch):
    """Index files survive pip upgrade but old config was replaced with species: {}.
    Migration should auto-register species based on index directory names."""
    data_dir = tmp_path / "data"
    monkeypatch.setenv("HCRPROBEDESIGN_DATA_DIR", str(data_dir))

    # Set up old package dir with index files but EMPTY species config
    # (simulates pip replacing HCRconfig.yaml with the clean default)
    old_pkg = tmp_path / "old_pkg"
    old_indices = old_pkg / "indices" / "mm10"
    old_indices.mkdir(parents=True)
    for suffix in ["1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2"]:
        (old_indices / f"mm10.{suffix}").write_text("fake")

    old_config_path = old_pkg / "HCRconfig.yaml"
    old_config = {
        "species": {},
        "default_params": {"tileSize": 52},
    }
    with open(old_config_path, "w") as fh:
        yaml.safe_dump(old_config, fh)

    monkeypatch.setattr(_datadir, "_PACKAGE_DIRECTORY", str(old_pkg))
    monkeypatch.setattr(_datadir, "_SEED_CONFIG", str(old_config_path))

    _datadir.ensure_data_dir()

    # Index files should be migrated
    new_index_dir = data_dir / "indices" / "mm10"
    assert new_index_dir.exists()
    bt2_files = glob.glob(str(new_index_dir / "*.bt2*"))
    assert len(bt2_files) == 6

    # Species should be auto-registered from the directory name
    new_config = yaml.safe_load((data_dir / "HCRconfig.yaml").read_text())
    assert "mm10" in new_config["species"]
    assert new_config["species"]["mm10"]["bowtie2_index"] == os.path.join(
        str(data_dir), "indices", "mm10", "mm10"
    )


def test_has_old_data_empty(tmp_path, monkeypatch):
    monkeypatch.setattr(_datadir, "_PACKAGE_DIRECTORY", str(tmp_path))
    # Write an empty config
    config_path = tmp_path / "HCRconfig.yaml"
    config_path.write_text("species: {}\n")
    indices_dir = tmp_path / "indices"
    indices_dir.mkdir()

    has_species, index_dirs = _datadir._has_old_data()
    assert has_species is False
    assert index_dirs == []


def test_has_old_data_with_species(tmp_path, monkeypatch):
    monkeypatch.setattr(_datadir, "_PACKAGE_DIRECTORY", str(tmp_path))
    config_path = tmp_path / "HCRconfig.yaml"
    config_path.write_text("species:\n  mouse:\n    bowtie2_index: indices/mm10/mm10\n")

    species_dir = tmp_path / "indices" / "mm10"
    species_dir.mkdir(parents=True)
    (species_dir / "mm10.1.bt2").write_text("fake")

    has_species, index_dirs = _datadir._has_old_data()
    assert has_species is True
    assert "mm10" in index_dirs
