"""List installed reference genomes and their configuration details."""

import argparse
import glob
import os
import sys

import yaml

PACKAGE_DIRECTORY = os.path.dirname(os.path.abspath(__file__))
DEFAULT_CONFIG_PATH = os.path.join(PACKAGE_DIRECTORY, "HCRconfig.yaml")


def load_config(config_path=DEFAULT_CONFIG_PATH):
    """
    Load the HCRconfig.yaml file.

    :param config_path: Path to the YAML configuration file.
    :return: Parsed config dictionary (empty if missing).
    """
    if not os.path.exists(config_path):
        return {}
    with open(config_path, "r") as handle:
        return yaml.safe_load(handle) or {}


def _resolve_absolute_index_path(index_prefix):
    """
    Resolve an index prefix to an absolute path.

    Package-relative paths are resolved against the package directory.

    :param index_prefix: Bowtie2 index prefix (relative or absolute).
    :return: Absolute path to the index prefix.
    """
    if os.path.isabs(index_prefix):
        return index_prefix
    return os.path.join(PACKAGE_DIRECTORY, index_prefix)


def _check_index_files(index_prefix):
    """
    Check whether Bowtie2 index files exist for the given prefix.

    :param index_prefix: Absolute path to the Bowtie2 index prefix.
    :return: List of existing index files.
    """
    return sorted(glob.glob(f"{index_prefix}*.bt2*"))


def list_references(config_path=DEFAULT_CONFIG_PATH):
    """
    Gather information about installed reference genomes.

    :param config_path: Path to HCRconfig.yaml.
    :return: Dictionary with 'species' list and 'default_params' dict.
    """
    config = load_config(config_path)
    species_config = config.get("species", {}) or {}
    default_params = config.get("default_params", {}) or {}

    species_info = []
    for name, entry in species_config.items():
        index_prefix = entry.get("bowtie2_index", "")
        abs_path = _resolve_absolute_index_path(index_prefix)
        index_files = _check_index_files(abs_path)
        species_info.append({
            "name": name,
            "bowtie2_index": index_prefix,
            "absolute_path": abs_path,
            "index_files": index_files,
            "installed": len(index_files) > 0,
        })

    return {"species": species_info, "default_params": default_params}


def format_references(info):
    """
    Format reference genome information as a human-readable string.

    :param info: Dictionary returned by :func:`list_references`.
    :return: Formatted string for display.
    """
    lines = []
    lines.append("=" * 60)
    lines.append("HCRProbeDesign - Installed Reference Genomes")
    lines.append("=" * 60)

    species_list = info.get("species", [])
    if not species_list:
        lines.append("")
        lines.append("No reference genomes registered.")
        lines.append("")
        lines.append("To add a reference genome, use one of:")
        lines.append("  buildGenomeIndex --species <name> --fasta <file>")
        lines.append("  fetchMouseIndex")
    else:
        for i, sp in enumerate(species_list):
            lines.append("")
            status = "INSTALLED" if sp["installed"] else "MISSING"
            lines.append(f"  [{i + 1}] {sp['name']}  ({status})")
            lines.append(f"      Index prefix : {sp['bowtie2_index']}")
            lines.append(f"      Absolute path: {sp['absolute_path']}")
            if sp["installed"]:
                lines.append(f"      Index files   : {len(sp['index_files'])} files")
            else:
                lines.append("      Index files   : none found")
            lines.append(f"      CLI usage     : designProbes --species {sp['name']}")

    lines.append("")
    lines.append("-" * 60)
    lines.append("Default Parameters (HCRconfig.yaml)")
    lines.append("-" * 60)

    default_params = info.get("default_params", {})
    if default_params:
        # Map config keys to their CLI flag equivalents
        param_flags = {
            "tileSize": "--tileSize",
            "minGC": "--minGC",
            "maxGC": "--maxGC",
            "targetGC": "--targetGC",
            "dTmMax": "--dTmMax",
            "minGibbs": "--minGibbs",
            "maxGibbs": "--maxGibbs",
            "targetGibbs": "--targetGibbs",
            "maxRunLength": "--maxRunLength",
            "maxProbes": "--maxProbes",
            "maxRunMismatches": "--maxRunMismatches",
            "num_hits_allowed": "--num_hits_allowed",
        }
        for key, value in default_params.items():
            flag = param_flags.get(key, f"--{key}")
            lines.append(f"  {key:<20s} = {value:<10}  (flag: {flag})")
    else:
        lines.append("  No default parameters defined.")

    lines.append("")
    lines.append("=" * 60)
    return "\n".join(lines)


def main():
    """CLI entry point for listing installed reference genomes."""
    parser = argparse.ArgumentParser(
        description="List installed reference genomes and default parameters.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--config",
        default=DEFAULT_CONFIG_PATH,
        help="Path to HCRconfig.yaml (default: package config)",
    )
    args = parser.parse_args()

    info = list_references(config_path=args.config)
    print(format_references(info))
