# CLI Reference

## designProbes
Design probes for the first record in a FASTA file.

```bash
designProbes targets.fa --species mouse --channel B1 --output probes.tsv --idt probes.idt
```

Common flags:
- `--no-genomemask`: skip Bowtie2 uniqueness checks
- `--index /path/to/index`: override Bowtie2 index prefix
- `--tileSize 52`: tile size (probe length before splitting)
- `--minGC`, `--maxGC`: GC content bounds
- `--maxProbes`: maximum number of probes to emit
- `--calcPrice`: estimate oligo synthesis cost

Note: genome masking is enabled by default and requires a registered species.
Use `fetchMouseIndex` or `buildGenomeIndex` first, or pass `--index` to point
directly to a Bowtie2 index.

## designProbesBatch
Design probes for every record in a FASTA file.

```bash
designProbesBatch targets.fa --species mouse --channel B1 --output probes.tsv --idt probes.idt
```

## fetchMouseIndex
Download a prebuilt mm10 Bowtie2 index and register it under the package indices directory.

```bash
fetchMouseIndex
```

## listReferences
List all installed reference genomes and their configuration details, including default
parameters and the CLI flags needed to use them.

```bash
listReferences
```

Use `--config` to inspect an alternate configuration file:
```bash
listReferences --config /path/to/HCRconfig.yaml
```

The output shows, for each registered species:
- Species name and installation status
- Bowtie2 index prefix path (relative and absolute)
- Number of index files found on disk
- The `--species` flag to pass to `designProbes`

It also displays all default parameters from `HCRconfig.yaml` with their corresponding
CLI flags.

## buildGenomeIndex
Build a Bowtie2 index for a new reference genome and register it in `HCRconfig.yaml`.

```bash
buildGenomeIndex --species zebrafish --fasta /path/to/genome.fa --threads 8
```

For large genomes (> 4 billion bases), use the `--large-index` flag:
```bash
buildGenomeIndex --species eberryi --fasta /path/to/genome.fa --threads 8 --large-index
```
