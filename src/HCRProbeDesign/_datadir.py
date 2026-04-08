"""Manage the user data directory for HCRProbeDesign.

Configuration and Bowtie2 indices are stored in a persistent user directory
(``~/.hcrprobedesign/`` by default) so they survive package upgrades.

Override the location by setting the ``HCRPROBEDESIGN_DATA_DIR`` environment
variable.
"""

import glob
import os
import shutil
import sys

import yaml

_PACKAGE_DIRECTORY = os.path.dirname(os.path.abspath(__file__))
_DEFAULT_DATA_DIR = os.path.join(os.path.expanduser("~"), ".hcrprobedesign")
_SEED_CONFIG = os.path.join(_PACKAGE_DIRECTORY, "HCRconfig.yaml")


def get_data_dir():
    """
    Return the path to the user data directory.

    Uses ``HCRPROBEDESIGN_DATA_DIR`` if set, otherwise ``~/.hcrprobedesign``.

    :return: Absolute path to the data directory.
    """
    return os.environ.get("HCRPROBEDESIGN_DATA_DIR", _DEFAULT_DATA_DIR)


def get_config_path():
    """
    Return the path to the user's ``HCRconfig.yaml``.

    :return: Absolute path to the config file.
    """
    return os.path.join(get_data_dir(), "HCRconfig.yaml")


def get_indices_dir():
    """
    Return the path to the user's indices directory.

    :return: Absolute path to the indices directory.
    """
    return os.path.join(get_data_dir(), "indices")


def _load_yaml(path):
    """Load a YAML file, returning an empty dict if missing."""
    if not os.path.exists(path):
        return {}
    with open(path, "r") as fh:
        return yaml.safe_load(fh) or {}


def _old_package_config_path():
    """Return the path to the old package-relative config file."""
    return os.path.join(_PACKAGE_DIRECTORY, "HCRconfig.yaml")


def _old_package_indices_dir():
    """Return the path to the old package-relative indices directory."""
    return os.path.join(_PACKAGE_DIRECTORY, "indices")


def _has_old_data():
    """
    Check whether there is user data in the old package-relative locations.

    :return: Tuple of (has_species_entries, list_of_index_dirs).
    """
    old_config = _load_yaml(_old_package_config_path())
    species = old_config.get("species", {}) or {}
    has_species = len(species) > 0

    old_indices = _old_package_indices_dir()
    index_dirs = []
    if os.path.isdir(old_indices):
        for entry in os.listdir(old_indices):
            entry_path = os.path.join(old_indices, entry)
            if os.path.isdir(entry_path):
                bt2_files = glob.glob(os.path.join(entry_path, "*.bt2*"))
                if bt2_files:
                    index_dirs.append(entry)

    return has_species, index_dirs


def _migrate_old_data():
    """
    Migrate config entries and index files from the package directory to the
    user data directory.

    Prints progress messages to stderr.
    """
    data_dir = get_data_dir()
    new_config_path = get_config_path()
    new_indices_dir = get_indices_dir()

    old_config_path = _old_package_config_path()
    old_indices_dir = _old_package_indices_dir()

    old_config = _load_yaml(old_config_path)
    new_config = _load_yaml(new_config_path)

    old_species = old_config.get("species", {}) or {}
    new_species = new_config.setdefault("species", {})

    migrated_any = False
    migrated_index_dirs = []

    # Migrate index directories
    if os.path.isdir(old_indices_dir):
        for entry in os.listdir(old_indices_dir):
            old_entry_path = os.path.join(old_indices_dir, entry)
            if not os.path.isdir(old_entry_path):
                continue
            bt2_files = glob.glob(os.path.join(old_entry_path, "*.bt2*"))
            if not bt2_files:
                continue

            new_entry_path = os.path.join(new_indices_dir, entry)
            if os.path.exists(new_entry_path):
                print(
                    f"  Skipping index '{entry}' (already exists in {new_indices_dir})",
                    file=sys.stderr,
                )
            else:
                print(
                    f"  Moving index '{entry}' -> {new_entry_path}",
                    file=sys.stderr,
                )
                shutil.copytree(old_entry_path, new_entry_path)
                migrated_any = True

            migrated_index_dirs.append(entry)

    # Migrate species config entries from old config
    for name, entry in old_species.items():
        if name in new_species:
            continue
        old_index = entry.get("bowtie2_index", "")
        # Rewrite relative paths to point to the new indices dir
        if not os.path.isabs(old_index):
            # e.g. "indices/mm10/mm10" -> new absolute path
            basename = old_index
            if basename.startswith("indices/") or basename.startswith("indices\\"):
                basename = basename[len("indices/"):]
            new_index = os.path.join(new_indices_dir, basename)
        else:
            new_index = old_index
        new_species[name] = {"bowtie2_index": new_index}
        migrated_any = True
        print(
            f"  Migrated species '{name}' -> {new_index}",
            file=sys.stderr,
        )

    # Auto-register species for migrated index dirs that have no config entry.
    # This handles the case where pip replaced the old HCRconfig.yaml with a
    # fresh copy (species: {}) during upgrade, so the config entries were lost
    # even though the index files survived.
    for dirname in migrated_index_dirs:
        # Check if any existing species entry already points to this index dir
        index_prefix = os.path.join(new_indices_dir, dirname, dirname)
        already_registered = any(
            sp.get("bowtie2_index", "") == index_prefix
            or os.path.basename(sp.get("bowtie2_index", "")).startswith(dirname)
            for sp in new_species.values()
        )
        if dirname not in new_species and not already_registered:
            new_species[dirname] = {"bowtie2_index": index_prefix}
            migrated_any = True
            print(
                f"  Auto-registered species '{dirname}' -> {index_prefix}",
                file=sys.stderr,
            )

    # Preserve default_params from old config if not present in new
    if "default_params" not in new_config and "default_params" in old_config:
        new_config["default_params"] = old_config["default_params"]
        migrated_any = True

    if migrated_any:
        with open(new_config_path, "w") as fh:
            yaml.safe_dump(new_config, fh, sort_keys=False)

    return migrated_any


def ensure_data_dir():
    """
    Ensure the user data directory exists and is seeded.

    On first run, creates ``~/.hcrprobedesign/`` with a seeded
    ``HCRconfig.yaml`` and ``indices/`` directory.  If data exists in the
    old package-relative location, it is migrated automatically.
    """
    data_dir = get_data_dir()
    config_path = get_config_path()
    indices_dir = get_indices_dir()

    if os.path.isdir(data_dir) and os.path.exists(config_path):
        return  # Already initialized

    first_init = not os.path.isdir(data_dir)
    os.makedirs(indices_dir, exist_ok=True)

    # Seed config from package default if no config exists yet
    if not os.path.exists(config_path):
        if os.path.exists(_SEED_CONFIG):
            shutil.copy2(_SEED_CONFIG, config_path)
        else:
            # Minimal fallback
            with open(config_path, "w") as fh:
                yaml.safe_dump({"species": {}, "default_params": {}}, fh)

    if first_init:
        print(
            f"Initialized HCRProbeDesign data directory: {data_dir}",
            file=sys.stderr,
        )

    # Check for data in old package location and migrate
    has_species, index_dirs = _has_old_data()
    if has_species or index_dirs:
        print(
            "Found reference data in old package directory. Migrating...",
            file=sys.stderr,
        )
        _migrate_old_data()
        print("Migration complete.", file=sys.stderr)
