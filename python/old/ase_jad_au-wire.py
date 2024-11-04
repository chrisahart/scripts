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
import ase
from ase.io import read, write

"""
Wrap .xyz using ASE
"""

folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-wire/jad/structures'
folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-wire-ismael/jad/structures/chris/test'
filename_input = "{}/Au_Junc_wrapped_sorted.xyz".format(folder)
filename_output_lead = "{}/Au_Junc_wrapped_sorted_lead.xyz".format(folder)
filename_output_all = "{}/Au_Junc_wrapped_sorted_all.xyz".format(folder)
atoms = ase.io.read(filename_input, format='xyz')

# Wrap cell
cell = [[1.7811878396732908E+01, 0.0000000000000000E+00, 0.0000000000000000E+00],
        [0.0000000000000000E+00, 1.5425539180689931E+01, 0.0000000000000000E+00],
        [0.0000000000000000E+00, 0.0000000000000000E+00, 3.0226752151153697E+01]]
atoms.cell = cell
atoms.pbc = True

lead_left_1 = atoms[[atom.index for atom in atoms if atom.symbol == 'Au' and atom.position[2] < 2]]
lead_left_2 = atoms[[atom.index for atom in atoms if atom.symbol == 'Au' and atom.position[2] < 2]]
lead_left_2.set_positions(lead_left_2.positions + [0, 0, 2])
lead_left = lead_left_1 + lead_left_2
lead_right = lead_left_1 + lead_left_2
lead_right.set_positions(lead_right.positions + [0, 0, 32])

# Select atoms with z > 2
selected_indices = [atom.index for atom in atoms if atom.position[2] > 2]
em = atoms[selected_indices]
em.positions[:, 2] += 2

all = lead_left + em + lead_right

z = 3.0226752151153697E+01 + 2*3
cell = [[1.7811878396732908E+01, 0.0000000000000000E+00, 0.0000000000000000E+00],
        [0.0000000000000000E+00, 1.5425539180689931E+01, 0.0000000000000000E+00],
        [0.0000000000000000E+00, 0.0000000000000000E+00, z]]
all.cell = cell
all.pbc = True

# Copy first 36 atoms
# print(file_coord.shape[0])
# for i in range(file_coord.shape[0]):
#     if file_coord['Z'][i] > 20:
#         file_coord['Z'][i] = file_coord['Z'][i]+20

# Save file
write(filename_output_lead, lead_left)
write(filename_output_all, all)

if __name__ == "__main__":
    print('Finished.')
    view(all)
    # plt.show()

