# HCRProbeDesign

## Overview
HCRProbeDesign is a command-line and Python package for designing HCR v3.0 split-initiator probe pairs
from a target FASTA sequence. It tiles target sequences, filters by GC content, melting temperature,
homopolymers, and hairpins, and can optionally enforce genome uniqueness with Bowtie2.

Key tools:
- `designProbes`: primary probe design CLI.
- `fetchMouseIndex`: download and install the mm10 Bowtie2 index.
- `buildGenomeIndex`: build and register a new reference genome index.

## Installation
### Prerequisites
- Python >= 3.6
- Bowtie2 (`bowtie2` and `bowtie2-build`) for genome masking and index building
- Python dependencies are installed via `pip` (primer3-py, biopython, beautifulsoup4, pysam, zipfile36, pyyaml)

### Install from source
```bash
git clone https://github.com/gofflab/HCRProbeDesign.git
cd HCRProbeDesign

# Optional: use the provided conda environment
conda env create -f environment.yaml
conda activate HCRProbeDesign

# Install the package
pip install -e .
```

## Adding a new reference genome
HCRProbeDesign keeps Bowtie2 index paths in `HCRconfig.yaml`. The `buildGenomeIndex` utility builds
an index and registers it automatically.

```bash
buildGenomeIndex --species zebrafish --fasta /path/to/genome.fa --threads 8
```

Notes:
- `--fasta` can be repeated and can point to a directory; all `.fa`, `.fasta`, or `.fna` files inside
  will be used.
- By default, indices are written under the package `indices/` directory and the config is updated.
- Use `--indices-dir` to write indices elsewhere and `--config` to update a specific config file.
- Use `--force` to overwrite an existing index or config entry.

If you are working with mouse (mm10), you can use the prebuilt index:
```bash
fetchMouseIndex
```

## Usage instructions
### 1) Prepare a FASTA file
`designProbes` reads the first FASTA record in the input file.

```text
>MyTarget
ACGTACGTACGTACGTACGTACGTACGTACGT
```

### 2) Run probe design
```bash
designProbes targets.fa --species mouse --channel B1 --output probes.tsv --idt probes.idt
```

Common flags:
- `--no-genomemask`: skip Bowtie2 genome masking (faster, but no uniqueness check).
- `--index /path/to/index`: override the Bowtie2 index path.
- `--tileSize 52`, `--minGC`, `--maxGC`, `--maxProbes`: tune probe selection.

Outputs:
- Probe designs are written to `--output` (stdout by default).
- `--idt` writes an IDT-friendly TSV for ordering.
- Bowtie2 genome masking writes a SAM file named `{targetName}.sam` in the working directory.
