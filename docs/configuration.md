# Configuration

HCRProbeDesign stores reference genome index paths in `HCRconfig.yaml`.
The file is located in the package directory and is read at runtime.

## Example
```yaml
species:
  mouse:
    bowtie2_index: indices/mm10/mm10
  human:
    bowtie2_index: /path/to/hg38/index_prefix
```

## Viewing the current configuration
Run `listReferences` to display all registered species and default parameters:

```bash
listReferences
```

## Notes
- Paths can be absolute or package-relative.
- `buildGenomeIndex` updates this file automatically.
- `fetchMouseIndex` also registers the mouse index for you.
- Use `--config` to write to a different config file when building indices.
- Use `listReferences --config /path/to/HCRconfig.yaml` to inspect a specific config file.
