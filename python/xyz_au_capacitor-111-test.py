from __future__ import division, print_function
import pandas as pd
import numpy as np
import glob
from general import load_coordinates
from general import print_xyz
import matplotlib.pyplot as plt
import math
import ase
from ase.build import surface
from ase.visualize import view
from ase import Atoms
from ase.build import bulk
from ase.io import read, write
from ase.build import make_supercell
from ase.build import fcc111
from pathlib import Path

"""
    Create Cu-BDT junction
"""

# Coordinates
#####################
### build surface ###
#####################
indices = (1, 1, 1)  # Miller index
layers = 6  # number of layers
vacuum = 10.0  # vacuum thickness

### ase.build.surface(lattice, indices, layers, vacuum=None(Ang), tol=1e-10, periodic=False)
Mobulk = bulk('Au', 'fcc', a=3.16, cubic=True)
atoms = surface(Mobulk, indices, layers, vacuum)
view(atoms)
### ###