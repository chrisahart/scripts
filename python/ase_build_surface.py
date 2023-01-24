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

"""

"""

folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Masters/2021-2022/Yike/Pt3Ni/structures'

a = 3.85
Pt3Rh = Atoms('NiPt3',
              scaled_positions=[(0, 0, 0),
                                (0.5, 0.5, 0),
                                (0.5, 0, 0.5),
                                (0, 0.5, 0.5)],
              cell=[a, a, a],
              pbc=True)


supercell = make_supercell(Pt3Rh, np.diag([2, 2, 2]))

s3 = surface(Pt3Rh, (1, 1, 1), layers=6)
s3.center(vacuum=10, axis=2)

write('{}/ase_struct_222_001_supercell.xyz'.format(folder), supercell)
write('{}/ase_struct_111_111_slab.xyz'.format(folder), s3)

# view(Pt3Rh)
# view(supercell)
view(s3)
