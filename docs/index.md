# HCRProbeDesign

HCRProbeDesign is a command-line and Python package for designing HCR v3.0 split-initiator probe pairs
from a target FASTA sequence. It tiles target sequences, filters by GC content, melting temperature,
homopolymers, and hairpins, and can optionally enforce genome uniqueness with Bowtie2.

## Key tools
- `designProbes`: primary probe design CLI for single-record FASTA inputs
- `designProbesBatch`: batch probe design for multi-record FASTA inputs
- `fetchMouseIndex`: download, install, and register the mm10 Bowtie2 index
- `buildGenomeIndex`: build and register a new reference genome index

## Documentation layout
- Installation and prerequisites
- Usage walkthroughs and recommended workflows
- CLI reference with common flags
- Configuration details for `HCRconfig.yaml`
- API reference generated from inline docstrings
