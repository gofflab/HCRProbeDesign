# Reference Genomes

## Why add a reference genome
HCRProbeDesign can check whether candidate probes map uniquely in a genome using Bowtie2.
To enable this for species or assemblies beyond the default mouse index, you need to add a
Bowtie2 index and register it in `HCRconfig.yaml`. This lets you:
- Run genome masking for non-mouse species or custom assemblies.
- Keep indices organized and reusable across runs.
- Avoid turning off uniqueness checks with `--no-genomemask`.

If you do not need uniqueness checks, you can skip this entirely and use `--no-genomemask`.

## How to add a new reference genome
### 1) Prepare a FASTA
Download or generate the reference FASTA for your species or assembly. Gzipped FASTA files
are supported.

### 2) Build and register the Bowtie2 index
Use `buildGenomeIndex` to build the index and register it in `HCRconfig.yaml`.

```bash
buildGenomeIndex --species zebrafish --fasta /path/to/genome.fa --threads 8
```

Notes:
- `--fasta` is repeatable and can point to a directory; all `.fa`, `.fasta`, or `.fna` files
  (including `.gz`) will be included.
- Indices are written under the package `indices/` directory by default.
- Use `--indices-dir` to write indices elsewhere and `--config` to update a specific config file.
- Use `--force` to overwrite an existing index or config entry.

### 3) Verify the configuration
`buildGenomeIndex` adds an entry like this to `HCRconfig.yaml`:

```yaml
species:
  zebrafish:
    bowtie2_index: indices/zebrafish/zebrafish
```

You can now run probe design with `--species zebrafish`. If needed, override the index
location per run with `--index /path/to/index_prefix`.

## Prebuilt mouse index
For mouse (mm10), you can install a prebuilt index instead of building one:

```bash
fetchMouseIndex
```
