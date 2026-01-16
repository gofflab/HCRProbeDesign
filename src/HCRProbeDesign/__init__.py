"""Top-level package for HCRProbeDesign utilities and probe design workflows."""

from . import probeDesign
from . import thermo
from . import sequencelib

import os
try:
    from importlib import metadata as importlib_metadata
except ImportError:  # Python < 3.8
    importlib_metadata = None

_ROOT = os.path.abspath(os.path.dirname(__file__))

def index_path():
    """Return the default package-relative path for Bowtie2 indices."""
    return os.path.join(_ROOT, 'indices')

def _resolve_version():
    if importlib_metadata is not None:
        try:
            return importlib_metadata.version("hcrprobedesign")
        except Exception:
            pass
    try:
        import pkg_resources
        return pkg_resources.get_distribution("hcrprobedesign").version
    except Exception:
        return "unknown"

__version__ = _resolve_version()
del _resolve_version
