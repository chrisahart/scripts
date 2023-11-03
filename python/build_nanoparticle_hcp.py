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
from wulffpack import (SingleCrystal,
                       Decahedron,
                       Icosahedron)
from ase.build import bulk
from ase.io import write

# https://wiki.fysik.dtu.dk/ase/ase/cluster/cluster.html
# https://wulffpack.materialsmodeling.org/moduleref/single_crystal.html
# This uses WulffPack, which uses ASE to generates Wulff constructions
# Create nanoparticle with only (100) and (111) surfaces, inspired by Clotilde nanoparticle
folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/theory_support/hannah/structures/nanoparticle-ase'

# prim = bulk('Co', crystalstructure='hcp')
# surface_energies = {(1, 0, -1, 0): 1.1,
#                     (0, 0, 0, 1): 1.0,
#                     (1, 1, -2, 0): 1.0}
# particle = SingleCrystal(surface_energies,
#                          primitive_structure=prim,
#                          natoms=5000)
# particle.view()
# write('{}/test.xyz'.format(folder), particle.atoms)

surface_energies = {(1, 0, 0): 1.0,
                    (1, 1, 1): 0.866}
# prim = bulk('Pt', a=3.94, crystalstructure='fcc')
prim = bulk('Pt', a=2.76, b=2.76, c=4.79, crystalstructure='hcp')
particle = SingleCrystal(surface_energies,
                      primitive_structure=prim,
                      natoms=201)
print(particle.atoms)
write('{}/test.xyz'.format(folder), particle.atoms)
particle.view()

# lc = 3.61000  # Cu
# lc = 3.94  # Pt
# surfaces = [(1, 0, 0),  (1, 1, 1)]

# Build by specifying surface energies
# energies = [1, 0.866]
# atoms = wulff_construction('Pd', surfaces, energies, size=201, structure='hcp', rounding='above',
#                            latticeconstant=[2.88, 2.76])

# Hexagonal testing
# layers = [10, 11]  # Pt 201
# atoms = ase.cluster.HexagonalClosedPacked('Pt', surfaces, layers,  latticeconstant=[4, 3])

# num_atoms = atoms.get_global_number_of_atoms()
# print('num_atoms', num_atoms)

# Save
# write('{}/struct-pt_atoms-{}.xyz'.format(folder, num_atoms), atoms)
# write('{}/pyramid.xyz'.format(folder), atoms)

# if __name__ == "__main__":
#     view(atoms)
    # plt.show()
    # print('Finished.')
