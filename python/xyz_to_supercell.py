from __future__ import division, print_function
import pandas as pd
import numpy as np
import glob
from general import load_coordinates
from general import print_xyz
import matplotlib.pyplot as plt
import math

"""
    .xyz unit cell to supercell
"""

# Au-BDT junction
# origin_atom = 0
# cell = [4.1712875, 4.1712875, 4.1712875]
# supercell = [3, 3, 3]
# coordinate_max = [math.inf, math.inf, 9]
# folder = '/Volumes/Storage/Data/Work/Postdoc/Work/calculations/testing/cp2k-smeagol/bdt/structures'
# input_filename = '/unit_cell.xyz'
# output_filename = '/hollow-site/5-4_3x3/supercell_3x3-5.xyz'

# Au capacitor
# origin_atom = 0
# cell = [4.1712875, 4.1712875, 4.1712875]
# supercell = [3, 3, 3]
# coordinate_max = [math.inf, math.inf, math.inf]
# folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/au-capacitor/structures'
# input_filename = '/unit_cell.xyz'
# output_filename = 'supercell_3x3-3.xyz'

# Au capacitor test
origin_atom = 0
cell = [4.1712875, 4.1712875, 4.1712875]
supercell = [2, 2, 1]
coordinate_max = [math.inf, math.inf, math.inf]
folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/au-capacitor-100/structures'
input_filename = '/unit_cell.xyz'
output_filename = 'layer.xyz'

# Read number of atoms and labels from .xyz file
cols = ['Species', 'X', 'Y', 'Z']
file_coord, num_atoms, species = load_coordinates.load_file_coord(folder, input_filename, cols)
file_coord = file_coord.reset_index(drop=True)

# Center around atom
translate = -1 * np.array([file_coord['X'][origin_atom], file_coord['Y'][origin_atom], file_coord['Z'][origin_atom]])
for i in range(0, num_atoms):
    file_coord['X'][i] = file_coord['X'][i] + translate[0]
    file_coord['Y'][i] = file_coord['Y'][i] + translate[1]
    file_coord['Z'][i] = file_coord['Z'][i] + translate[2]

# Construct supercell
atoms = 0
for z in range(0, supercell[2]):
    for y in range(0, supercell[1]):
        for x in range(0, supercell[0]):
            for i in range(0, num_atoms):
                new_row = pd.Series({'X': cell[0] * x + file_coord['X'][i],
                                     'Y': cell[0] * y + file_coord['Y'][i],
                                     'Z': cell[0] * z + file_coord['Z'][i]})
                file_coord = pd.concat([file_coord, new_row.to_frame().T], ignore_index=True)

# Remove initial rows
file_coord = file_coord.drop(np.arange(start=0, stop=num_atoms))
file_coord = file_coord.reset_index(drop=True)

# Sort by Z
file_coord = file_coord.sort_values('Z')

# Delete columns greater than coordinate
num_atoms = file_coord.shape[0]
file_coord = file_coord[(file_coord['X'] <= coordinate_max[0])]
file_coord = file_coord[(file_coord['Z'] <= coordinate_max[1])]
file_coord = file_coord[(file_coord['Z'] <= coordinate_max[2])]
file_coord = file_coord.reset_index(drop=True)

# Print to file
num_atoms = file_coord.shape[0]
file_coord.insert(loc=0, column='A', value=['Au'] * num_atoms)
print_xyz.print_from_pandas(file_coord, num_atoms, '{}/{}'.format(folder, output_filename))
