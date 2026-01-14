"""Top-level package for HCRProbeDesign utilities and probe design workflows."""

from . import probeDesign
from . import thermo
from . import sequencelib

import os

_ROOT = os.path.abspath(os.path.dirname(__file__))

def index_path():
    """Return the default package-relative path for Bowtie2 indices."""
    return os.path.join(_ROOT, 'indices')
