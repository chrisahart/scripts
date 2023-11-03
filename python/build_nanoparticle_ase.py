import ase
import numpy as np
import matplotlib.pyplot as plt
import math
from ase.cluster.cubic import FaceCenteredCubic, SimpleCubic, BodyCenteredCubic
from ase.io import read, write
from ase import visualize
from ase.visualize import view
from ase.cluster.octahedron import Octahedron
from ase.cluster.decahedron import Decahedron
from ase.io import write
from ase.cluster.wulff import wulff_construction

# https://wiki.fysik.dtu.dk/ase/ase/cluster/cluster.html

# Create nanoparticle with only (100) and (111) surfaces, inspired by Clotilde nanoparticle
folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/theory_support/hannah/structures/nanoparticle-ase'

# lc = 3.61000  # Cu
lc = 3.94  # Pt
surfaces = [(1, 0, 0),  (1, 1, 1)]

# Build by specifying number of layers
# layers = [8, 6]  # Pt 1289
layers = [4, 3]  # Pt 201
atoms = ase.cluster.cubic.FaceCenteredCubic('Pt', surfaces, layers,  latticeconstant=lc)

# Build by specifying surface energies
# energies = [1, 0.866]
# atoms = wulff_construction('Pt', surfaces, energies, size=201, structure='fcc', rounding='below')

num_atoms = atoms.get_global_number_of_atoms()
print('num_atoms', num_atoms)

# Save
# write('{}/struct-pt_atoms-{}.xyz'.format(folder, num_atoms), atoms)
write('{}/pyramid.xyz'.format(folder), atoms)

if __name__ == "__main__":
    view(atoms)
    # plt.show()
    print('Finished.')
