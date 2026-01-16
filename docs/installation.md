# Installation

## Prerequisites
- Python >= 3.6
- Bowtie2 (`bowtie2` and `bowtie2-build`) for genome masking and index building

If you plan to use genome masking (enabled by default), register a reference
species with `fetchMouseIndex` or `buildGenomeIndex` after installing.

## Install from PyPI
```bash
pip install hcrprobedesign
```

## Install from source
```bash
git clone https://github.com/gofflab/HCRProbeDesign.git
cd HCRProbeDesign

# Optional: use the provided conda environment
conda env create -f environment.yaml
conda activate HCRProbeDesign

# Install the package
pip install -e .
```

## Documentation build dependencies
If you want to build the documentation locally:
```bash
pip install -r docs/requirements.txt
mkdocs serve
```
