import os

from . import probeDesign, sequencelib, thermo

_ROOT = os.path.abspath(os.path.dirname(__file__))


def index_path():
    return os.path.join(_ROOT, "indices")
