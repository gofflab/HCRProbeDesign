# Probe Design Workflow

This page describes the exact probe design pipeline used by `designProbes` and
`designProbesBatch`, including the order of operations, default parameters, and
output formats. The steps below are listed in the same order the code applies
them.

## Pipeline overview (ordered)
1. **Read FASTA record(s).**
   - `designProbes` uses only the first FASTA record.
   - `designProbesBatch` processes every record.
   - A FASTA header can override the channel using `channel=...`.
   - Genome masking is enabled by default and requires a registered species
     (use `fetchMouseIndex` or `buildGenomeIndex`) or an explicit `--index`.
2. **Optional repeat masking (disabled by default).**
   - If enabled, the input sequence is masked using RepeatMasker and masked bases
     are converted to `N`. Tiles containing masked bases are discarded later.
3. **Tile the target sequence (reverse-complemented).**
   - The sequence is scanned with a step size of 1 to generate 52-nt tiles by
     default.
   - Each tile is **reverse-complemented** so probes are antisense to the target.
   - Any tile containing `N` is discarded immediately.
4. **Filter homopolymer runs (C and G).**
   - Tiles are removed if they contain long runs of C's or G's beyond the allowed
     thresholds.
5. **Filter hairpins.**
   - Each tile is screened with `primer3` for self-hairpin structures.
   - Tiles with a predicted hairpin melting temperature >= 45 C are removed.
6. **Genome uniqueness filtering (Bowtie2).**
   - Remaining tiles are aligned to the reference genome index.
   - Tiles with more than the allowed number of hits are removed.
7. **GC-content filtering.**
   - Tiles must fall within the allowed GC% range.
8. **Gibbs free energy filtering.**
   - Tiles with binding free energies outside the allowed range are removed.
9. **Split tiles into probe halves.**
   - Each tile is split into two 25-mers (for a 52-mer tile), dropping the two
     middle bases to create a short gap.
10. **Optional dTm filtering.**
    - If enabled, tiles are removed when the temperature difference between
      the two halves is too large.
11. **Select top N non-overlapping tiles.**
    - Tiles are ranked by how close their Gibbs free energy is to the target.
    - Overlapping tiles are skipped to keep probes spread out.
12. **Add HCR initiators and spacers.**
    - Channel-specific initiator sequences are appended to the probe halves.

## Default parameters
These are the defaults used if you do not override them on the command line.

- `--tileSize`: 52
- `--channel`: B1
- `--species`: mouse
- `--minGC`: 45.0
- `--maxGC`: 55.0
- `--targetGC`: 50.0 (not currently used in ranking or filtering)
- `--minGibbs`: -70.0
- `--maxGibbs`: -50.0
- `--targetGibbs`: -60.0
- `--dTmMax`: 5.0
- `--dTmFilter`: off
- `--maxRunLength`: 7
- `--maxRunMismatches`: 2
- `--maxProbes`: 20
- `--num-hits-allowed`: 1
- `--no-genomemask`: off (genome masking is on by default)
- `--no-repeatmask`: on (repeat masking is disabled by default)

Other defaults:
- Hairpin filter threshold: 45 C (not currently configurable).
- Tile step size: 1 nt (not currently configurable).

## Thermodynamic parameters (what they mean and why they are used)
- **GC%**: Fraction of G/C bases in the tile. GC-rich probes bind more strongly
  because G-C base pairs have three hydrogen bonds. Filtering keeps probes within
  a moderate binding-strength window and improves uniformity.
- **Tm (melting temperature)**: Predicted temperature at which half of probe-
  target duplexes would melt. Reported for each full tile using `primer3`.
  This provides a quick proxy for binding strength across candidates.
- **dTm**: The absolute difference in Tm between the two split probe halves.
  Large dTm values imply one half binds much more strongly than the other, which
  can reduce uniformity. Filtering is optional and off by default.
- **Gibbs free energy (Gibbs FE)**: Predicted free energy of RNA/DNA binding,
  computed at 37 C with a salt correction (0.33 M). More negative values indicate
  stronger binding. Filtering removes probes that are too weak or too strong, and
  ranking prefers probes closest to `--targetGibbs`.
- **Hairpin Tm**: Predicted melting temperature of a self-hairpin in the probe
  sequence. Strong hairpins can compete with target binding, so tiles with a
  hairpin Tm >= 45 C are removed.

## Output formats
### Primary TSV output
`designProbes` and `designProbesBatch` write a tab-delimited table with these
columns:

- `name`: Tile identifier formatted as `record:start-end`.
- `probe`: The reverse-complemented tile sequence.
- `start`: 1-based start position in the original target sequence.
- `length`: Tile length (default 52).
- `P1`: Probe half with the odd initiator appended.
- `P2`: Probe half with the even initiator appended.
- `channel`: HCR channel used (e.g., B1).
- `GC`: GC% of the tile.
- `Tm`: `primer3` melting temperature of the full tile.
- `dTm`: Absolute Tm difference between the two halves.
- `GibbsFE`: Calculated Gibbs free energy of binding (kcal/mol).

Example (single row):
```text
MyTarget:1-52\tacg...\t1\t52\tP1SEQ...\tP2SEQ...\tB1\t50.00\t70.12\t1.83\t-60.45
```

### IDT ordering output (optional)
If `--idt` is provided, an additional file is written with two columns:
`Name` and `Sequence`. Each tile produces two entries:
- `{name}:{channel}:odd` for P1
- `{name}:{channel}:even` for P2

### Genome masking artifacts
When genome masking is enabled, Bowtie2 writes a SAM file in the working
directory named `{targetName}.sam`. This file lists all alignments used to
compute hit counts.
