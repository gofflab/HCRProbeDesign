from . import probeDesign
from . import thermo
from . import sequencelib

import os

_ROOT = os.path.abspath(os.path.dirname(__file__))

def index_path():
    return os.path.join(_ROOT, 'indices')
