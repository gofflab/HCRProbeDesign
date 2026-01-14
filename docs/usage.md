# Usage

## Basic workflow
1. Prepare a FASTA file with your target sequence.
2. Run `designProbes` (single record) or `designProbesBatch` (multi-record).
3. Review the TSV output and optional IDT ordering sheet.
4. Optionally inspect the Bowtie2 SAM file produced during genome masking.

## Example FASTA
```text
>MyTarget
ACGTACGTACGTACGTACGTACGTACGTACGT
```

## Single-record probe design
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
