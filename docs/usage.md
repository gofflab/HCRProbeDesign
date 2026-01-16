# Usage

## Basic workflow
1. Prepare a FASTA file with your target sequence.
2. Register a reference genome with `fetchMouseIndex` or `buildGenomeIndex`.
3. Run `designProbes` (single record) or `designProbesBatch` (multi-record).
4. Review the TSV output and optional IDT ordering sheet.
5. Optionally inspect the Bowtie2 SAM file produced during genome masking.

## Example FASTA
```text
>MyTarget
ACGTACGTACGTACGTACGTACGTACGTACGT
```

## Single-record probe design
Before running, register a species or provide `--index` (or skip genome masking with `--no-genomemask`).
```bash
designProbes targets.fa --species mouse --channel B1 --output probes.tsv --idt probes.idt
```

## Batch probe design
```bash
designProbesBatch targets.fa --species mouse --channel B1 --output probes.tsv --idt probes.idt
```

### Override channel per record
Add `channel=` to the FASTA header:
```text
>MyTarget channel=B2
ACGTACGTACGTACGTACGTACGT
```
