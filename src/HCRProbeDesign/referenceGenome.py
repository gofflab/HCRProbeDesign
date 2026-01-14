import argparse
import glob
import os
import shutil
import subprocess
import yaml

from . import index_path

PACKAGE_DIRECTORY = os.path.dirname(os.path.abspath(__file__))
DEFAULT_CONFIG_PATH = os.path.join(PACKAGE_DIRECTORY, "HCRconfig.yaml")
FASTA_EXTENSIONS = (".fa", ".fasta", ".fna", ".fa.gz", ".fasta.gz", ".fna.gz")


def load_config(config_path=DEFAULT_CONFIG_PATH):
    if not os.path.exists(config_path):
        return {}
    with open(config_path, "r") as handle:
        return yaml.safe_load(handle) or {}


def save_config(config, config_path=DEFAULT_CONFIG_PATH):
    with open(config_path, "w") as handle:
        yaml.safe_dump(config, handle, sort_keys=False)


def collect_fasta_inputs(paths):
    files = []
    for path in paths:
        if os.path.isdir(path):
            matches = []
            for ext in FASTA_EXTENSIONS:
                matches.extend(sorted(glob.glob(os.path.join(path, f"*{ext}"))))
            if not matches:
                raise FileNotFoundError(f"No FASTA files found in directory: {path}")
            files.extend(matches)
        else:
            if not os.path.exists(path):
                raise FileNotFoundError(f"FASTA path not found: {path}")
            files.append(path)

    seen = set()
    deduped = []
    for path in files:
        if path not in seen:
            seen.add(path)
            deduped.append(path)
    return deduped


def format_index_path(index_prefix):
    abs_prefix = os.path.abspath(index_prefix)
    abs_root = os.path.abspath(PACKAGE_DIRECTORY)
    try:
        common = os.path.commonpath([abs_prefix, abs_root])
    except ValueError:
        common = None
    if common == abs_root:
        return os.path.relpath(abs_prefix, abs_root)
    return abs_prefix


def build_bowtie2_index(fasta_paths, species, index_name=None, indices_dir=None, threads=1, force=False):
    bowtie2_build = shutil.which("bowtie2-build")
    if not bowtie2_build:
        raise RuntimeError("bowtie2-build not found in PATH")

    indices_dir = indices_dir or index_path()
    species_dir = os.path.join(indices_dir, species)
    os.makedirs(species_dir, exist_ok=True)

    index_name = index_name or species
    index_prefix = os.path.join(species_dir, index_name)
    existing = glob.glob(f"{index_prefix}*.bt2*")
    if existing and not force:
        raise FileExistsError(f"Index already exists at {index_prefix}. Use --force to overwrite.")
    if existing and force:
        for fname in existing:
            os.remove(fname)

    cmd = [bowtie2_build]
    if threads and threads > 1:
        cmd.extend(["--threads", str(threads)])
    cmd.extend([",".join(fasta_paths), index_prefix])
    subprocess.check_call(cmd)
    return index_prefix


def register_species(config_path, species, index_prefix, force=False):
    config = load_config(config_path)
    species_config = config.setdefault("species", {})
    if species in species_config and not force:
        raise ValueError(f"Species '{species}' already exists in config. Use --force to replace.")
    species_config[species] = {"bowtie2_index": format_index_path(index_prefix)}
    save_config(config, config_path)


def main():
    parser = argparse.ArgumentParser(
        description="Build a Bowtie2 index from a reference genome and register it."
    )
    parser.add_argument("--species", required=True, help="Species key to register in HCRconfig.yaml")
    parser.add_argument(
        "--fasta",
        required=True,
        action="append",
        help="FASTA file or directory (repeatable for multiple inputs)",
    )
    parser.add_argument("--index-name", help="Index basename (default: species)")
    parser.add_argument("--indices-dir", help="Output directory for indices (default: package indices)")
    parser.add_argument("--threads", type=int, default=1, help="Threads for bowtie2-build")
    parser.add_argument("--config", default=DEFAULT_CONFIG_PATH, help="Path to HCRconfig.yaml")
    parser.add_argument("--force", action="store_true", help="Overwrite existing index/config entry")
    args = parser.parse_args()

    fasta_paths = collect_fasta_inputs(args.fasta)
    index_prefix = build_bowtie2_index(
        fasta_paths,
        args.species,
        index_name=args.index_name,
        indices_dir=args.indices_dir,
        threads=args.threads,
        force=args.force,
    )
    register_species(args.config, args.species, index_prefix, force=args.force)
    print(f"Registered {args.species} with index {format_index_path(index_prefix)}")
