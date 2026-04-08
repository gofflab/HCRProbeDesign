# Configuration

HCRProbeDesign stores reference genome index paths and default parameters in
`HCRconfig.yaml`, located in the user data directory (`~/.hcrprobedesign/` by default).

## Data directory

All user data (configuration and Bowtie2 indices) lives in `~/.hcrprobedesign/`.
This directory is created automatically on first use and **persists across package
upgrades**, so you won't lose your indices when running `pip install -U hcrprobedesign`.

Override the location by setting the `HCRPROBEDESIGN_DATA_DIR` environment variable:
```bash
export HCRPROBEDESIGN_DATA_DIR=/path/to/custom/dir
```

### Migration from older versions

If you are upgrading from v0.3.0 or earlier (where data was stored inside the
package directory), HCRProbeDesign will automatically detect and migrate your
existing indices and species registrations on first run.

## Example
```yaml
species:
  mouse:
    bowtie2_index: indices/mm10/mm10
  human:
    bowtie2_index: /path/to/hg38/index_prefix
```

## Viewing the current configuration
Run `listReferences` to display all registered species, default parameters, and the
data directory location:

```bash
listReferences
```

## Notes
- Paths can be absolute or relative to the data directory.
- `buildGenomeIndex` updates this file automatically.
- `fetchMouseIndex` also registers the mouse index for you.
- Use `--config` to write to a different config file when building indices.
- Use `listReferences --config /path/to/HCRconfig.yaml` to inspect a specific config file.
