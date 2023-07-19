import ase
from ase.cluster.cubic import FaceCenteredCubic, SimpleCubic, BodyCenteredCubic
from ase.io import read, write
from ase import visualize
from ase.visualize import view
from ase.cluster.octahedron import Octahedron

from ase.cluster.decahedron import Decahedron
from ase.io import write
from ase.cluster import wulff_construction

lc = 3.61000

surfaces = [(1, 1, 1)]
layers = [3]
#
# surfaces = [(1, 0, 0), (1, 1, 0), (1, 1, 1)]
# layers = [6, 9, 5]

atoms = FaceCenteredCubic('Cu', surfaces, layers, latticeconstant=lc)


surfaces = [(1, 0, 0),  (1, 1, 1)]
layers = [4, 3]
lc = 3.61000
atoms = ase.cluster.cubic.FaceCenteredCubic('Cu', surfaces, layers,  latticeconstant=lc)

num_atoms = atoms.get_number_of_atoms()
print(num_atoms)

# surfaces = [(1, 0, 0), (1, 1, 1)]
# for a in range(0, 10):
#     for b in range(0, 10):
#         layers = [a, b]
#         atoms = FaceCenteredCubic('Cu', surfaces, layers, latticeconstant=lc)
#         num_atoms = atoms.get_number_of_atoms()
#         print(num_atoms)
#
#         if num_atoms == 116:
#             print('Done:', a, b)
#             exit()

# surfaces = [(1, 0, 0), (1, 1, 0), (1, 1, 1)]
# for a in range(0, 10):
#     for b in range(0, 10):
#         for c in range(0, 10):
#             layers = [a, b, c]
#             atoms = FaceCenteredCubic('Cu', surfaces, layers, latticeconstant=lc)
#             num_atoms = atoms.get_number_of_atoms()
#             print(num_atoms)
#
#             if num_atoms == 116:
#                 print('Done:', a, b, c)
#                 exit()

#Export surface
folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/theory_support/hannah/structures/nanoparticle-ase'
write('{}/struct-cu.xyz'.format(folder),atoms)
view(atoms)