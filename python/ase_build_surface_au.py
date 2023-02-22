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

folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/structures/111'
file = 'au_111_layers-7.xyz'

a = 4.168
au_au_z = a/2

surface = fcc111('Au', a=a, size=(4,4,7))

write('{}/{}'.format(folder, file), surface)

view(surface)
