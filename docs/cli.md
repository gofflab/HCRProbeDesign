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

## buildGenomeIndex
Build a Bowtie2 index for a new reference genome and register it in `HCRconfig.yaml`.

```bash
buildGenomeIndex --species zebrafish --fasta /path/to/genome.fa --threads 8
```
