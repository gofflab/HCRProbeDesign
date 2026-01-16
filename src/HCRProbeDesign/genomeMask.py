"""Genome masking utilities based on Bowtie2 alignments."""

# Debating whether or not to attempt Bowtie2 mapping for speed.
# Local would suck for install.
# Remote would suck for maintenance.

import subprocess
import tempfile
import pysam
from collections import defaultdict
import os
import urllib.request
import shutil
from zipfile import ZipFile
import argparse
import yaml

package_directory = os.path.dirname(os.path.abspath(__file__))
indices_directory = os.path.join(package_directory, "indices")
# indexLookup = {
#     'mouse': os.path.join(indices_directory,'mm10/mm10')
# }

def _load_config(config_path=None):
    """
    Load the HCRconfig.yaml file if present.

    :return: Parsed config dict (empty if missing).
    """
    config_path = config_path or os.path.join(package_directory, "HCRconfig.yaml")
    if not os.path.exists(config_path):
        return {}
    with open(config_path, "r") as file:
        return yaml.safe_load(file) or {}

def _save_config(config, config_path=None):
    """
    Save the HCRconfig.yaml file.

    :param config: Config dictionary to write.
    :param config_path: Optional config path override.
    :return: None.
    """
    config_path = config_path or os.path.join(package_directory, "HCRconfig.yaml")
    with open(config_path, "w") as file:
        yaml.safe_dump(config, file, sort_keys=False)

def _format_index_path(index_prefix):
    """
    Format an index prefix relative to the package when possible.

    :param index_prefix: Bowtie2 index prefix path.
    :return: Relative path if within package, otherwise absolute path.
    """
    abs_prefix = os.path.abspath(index_prefix)
    abs_root = os.path.abspath(package_directory)
    try:
        common = os.path.commonpath([abs_prefix, abs_root])
    except ValueError:
        common = None
    if common == abs_root:
        return os.path.relpath(abs_prefix, abs_root)
    return abs_prefix

def _register_species(species, index_prefix, force=False, config_path=None):
    """
    Register a species and its Bowtie2 index prefix in the config file.

    :param species: Species key to register.
    :param index_prefix: Bowtie2 index prefix path.
    :param force: Overwrite an existing species entry if True.
    :param config_path: Optional config path override.
    :return: None.
    :raises ValueError: If the species exists and force is False.
    """
    config = _load_config(config_path=config_path)
    species_config = config.setdefault("species", {})
    if species in species_config and not force:
        raise ValueError(f"Species '{species}' already exists in config. Use --force to replace.")
    species_config[species] = {"bowtie2_index": _format_index_path(index_prefix)}
    _save_config(config, config_path=config_path)


def _resolve_index(species, index):
    """
    Resolve a Bowtie2 index prefix for the requested species.

    :param species: Species key in HCRconfig.yaml.
    :param index: Optional explicit Bowtie2 index prefix.
    :return: Index prefix path (relative or absolute).
    :raises ValueError: If species is not registered and index is not provided.
    """
    if index:
        return index
    config = _load_config()
    species_config = config.get("species", {}) or {}
    if species not in species_config:
        raise ValueError(
            f"Species '{species}' is not registered in HCRconfig.yaml. "
            "Run buildGenomeIndex --species <name> --fasta <file_or_dir> "
            "or supply --index /path/to/bowtie2/index/prefix."
        )
    return species_config[species]["bowtie2_index"]

#TODO: make genomemask() take transient index argment if not default in species
def genomemask(fasta_string,handleName="tmp",species="mouse",nAlignments = 3, index=None):
    """
    Run Bowtie2 to align probe tiles and write a SAM file to disk.

    :param fasta_string: FASTA formatted string with probe sequences.
    :param handleName: Prefix for FASTA/SAM output files.
    :param species: Species key in HCRconfig.yaml.
    :param nAlignments: Number of alignments to report per read.
    :param index: Optional Bowtie2 index prefix override.
    :return: Bowtie2 subprocess return code.
    """
    fasta_file = f'{handleName}_reads.fa'
    tmpFasta = open(fasta_file,mode="w")
    #print(tmpFasta.name)
    tmpFasta.write(fasta_string)
    tmpFasta.close()
    sam_file = f'{handleName}.sam'
    index = _resolve_index(species, index)
    print(index)
    if os.path.isabs(index):
        index_path = index
    else:
        index_path = package_directory + "/" + index
    res = subprocess.call(["bowtie2", f"-k{nAlignments}", "-x", index_path, "-f", fasta_file, "-S", sam_file])
    return res

def countHitsFromSam(samFile):
    '''
    For each read in the sam file, add 1 to the count of hits for that read
    
    :param samFile: the name of the sam file
    :return: A dictionary with the read name as the key and the number of hits as the value.
    '''
    hitCounts = defaultdict(int)
    sam = pysam.AlignmentFile(samFile,"r")
    for read in sam.fetch():
        #print(read.get_tag(tag="TM:i")) # This doesn't work because pysam doesn't read the XS:i: tag
        if read.is_unmapped:
            hitCounts[read.query_name] += 0
        else:
            hitCounts[read.query_name] += 1
    return hitCounts

def test():
    """Quick manual test for Bowtie2 masking and SAM parsing."""
    fasta_string=">read1\ntacgagcttactggacgagcgtgactctgac\n>read2\nctgagctgatgcgacgtatctgatgctgtacgtgacg\n>read3\ncgagctatcgtactgagcggagcgcgggcgatat"
    handleName = 'test'
    proc_info = genomemask(fasta_string,handleName=handleName,species="mouse")
    res = countHitsFromSam(f'{handleName}.sam')
    print(res)

def install_index(url='https://genome-idx.s3.amazonaws.com/bt/mm10.zip', genome="mm10", species="mouse"):
    """
    Download and extract a prebuilt Bowtie2 index into the package indices.

    :param url: URL to a zipped Bowtie2 index archive.
    :param genome: Genome name used to name the extraction directory.
    :return: None.
    """
    # Parse arguments if any
    #Argument handling
    parser = argparse.ArgumentParser(description="Bowtie2 index retrieval and installation",formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--url", default=url, help="URL of a Bowtie2 index zip archive")
    parser.add_argument("--genome", default=genome, help="Genome name for the extracted folder")
    parser.add_argument("--species", default=species, help="Species key to register in HCRconfig.yaml")
    parser.add_argument("--indices-dir", default=indices_directory, help="Directory to install indices")
    parser.add_argument("--config", default=os.path.join(package_directory, "HCRconfig.yaml"), help="Path to HCRconfig.yaml")
    parser.add_argument("--force", action="store_true", help="Overwrite existing species entry")
    args = parser.parse_args()
    print(f'Downloading Bowtie2 index from {args.url} ...')
    index_folder = os.path.abspath(args.indices_dir)
    os.makedirs(index_folder, exist_ok=True)
    fname = os.path.join(index_folder, f"{args.genome}.zip")
    print(fname)
    if not os.path.exists(fname):
        urllib.request.urlretrieve(args.url, fname)
    with ZipFile(fname, 'r') as zip:
        # printing all the contents of the zip file
        #zip.printdir()
        # extracting all the files
        extract_dir = os.path.join(index_folder, args.genome)
        print(f'Extracting index from {fname} to {extract_dir}...')
        zip.extractall(extract_dir)
    index_prefix = os.path.join(index_folder, args.genome, args.genome)
    _register_species(args.species, index_prefix, force=args.force, config_path=args.config)
    print(f"Registered {args.species} with index {_format_index_path(index_prefix)}")
    print('Done!')
