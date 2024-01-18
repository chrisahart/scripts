import numpy as np
from matplotlib import pyplot as plt
from ase.io.cube import read_cube_data
import optparse
import os
from ase.build import fcc111
from ase.build import surface
from ase.visualize import view
from ase import Atoms
from ase.build import bulk
from ase.io import read, write
from ase.build import make_supercell
from ase.build import fcc111, fcc100, add_adsorbate
from ase.constraints import FixAtoms

"""

"""

folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-wire-ismael/jad/structures/jad_structures_replace-a'
filename_bulk = 'structure_bulk.xyz'
filename_em = 'structure_em.xyz'
file_jad = 'Au_Junc_wrapped_sorted.xyz'

a = 4.1983  # Jad

# Create the bulk Au(111) structure
slab = fcc100('Au', size=(4, 4, 7), a=a)
z_max = slab.positions[-1][2]

view(slab)
