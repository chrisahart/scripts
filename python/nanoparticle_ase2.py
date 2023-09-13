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

# https://wiki.fysik.dtu.dk/ase/ase/cluster/cluster.html

# Create nanoparticle with only (100) and (111) surfaces, inspired by Clotilde nanoparticle
folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/theory_support/hannah/structures/nanoparticle-ase'

# lc = 3.61000  # Cu
lc = 3.94  # Pt

surfaces = [(1, 0, 0),  (1, 1, 1)]
layers = [8, 6]
atoms = ase.cluster.cubic.FaceCenteredCubic('Pt', surfaces, layers,  latticeconstant=lc)
print(atoms.get_global_number_of_atoms())

# atoms_align_1 = np.array([598, 1283]) - 1
atoms_align_1 = np.array([492, 1181]) - 1

coord_1 = atoms.positions[atoms_align_1[0]]
coord_2 = atoms.positions[atoms_align_1[1]]

atoms.symbols[atoms_align_1[0]] = 'Au'
atoms.symbols[atoms_align_1[1]] = 'Au'

# Save
write('{}/testing.xyz'.format(folder), atoms)
view(atoms)

# Orientation can be fixed using align tool in Avogadro to align along y and z avis
# Unfortunately Avogadro crashes if too many atoms

print(coord_1, coord_2)
diff = coord_1 - coord_2
print(diff)

# diff = np.array([-2, 2, 1])
xyLength = np.sqrt(diff[0] * diff[0] + diff[1] * diff[1])
print(xyLength)
zAngle = np.degrees(np.arccos(diff[1] / xyLength))
print(zAngle)
vecLength = np.sqrt(diff[0] * diff[0] + diff[1] * diff[1] + diff[2] * diff[2])
print(vecLength)
xAngle = np.degrees(np.arccos(xyLength / vecLength))
print(xAngle)

atoms.rotate(xAngle, 'x')
atoms.rotate(zAngle, 'z')

view(atoms)

# Save
write('{}/testing-rotate.xyz'.format(folder), atoms)

if __name__ == "__main__":
    # view(atoms)
    # plt.show()
    print('Finished.')

