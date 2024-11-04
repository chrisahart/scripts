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
filename_input = "{}/Au_Junc.xyz".format(folder)
filename_output = "{}/Au_Junc_wrapped.xyz".format(folder)
atoms = ase.io.read(filename_input, format='xyz')

# Wrap cell
cell = [[1.7811878396732908E+01, 0.0000000000000000E+00, 0.0000000000000000E+00],
        [0.0000000000000000E+00, 1.5425539180689931E+01, 0.0000000000000000E+00],
        [0.0000000000000000E+00, 0.0000000000000000E+00, 3.0226752151153697E+01]]
atoms.cell = cell
atoms.pbc = True
atoms.wrap()

# Save file
write(filename_output, atoms)

if __name__ == "__main__":
    print('Finished.')
    view(atoms)
    # plt.show()

