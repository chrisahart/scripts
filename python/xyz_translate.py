from __future__ import division, print_function
import pandas as pd
import numpy as np
from general import load_coordinates
from general import print_xyz
import matplotlib.pyplot as plt

"""
    .xyz translate
    Reads .xyz file and translates atomic positions
"""

folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/theory_support/hannah/structures/carbide'
input_filename = 'hexagonal_supercell-2x2x2.xyz'
output_filename = 'hexagonal_supercell-2x2x2_translated2.xyz'

# Read number of atoms and labels from .xyz file
cols = ['Species', 'X', 'Y', 'Z']
file_coord, num_atoms, species = load_coordinates.load_file_coord(folder, input_filename, cols)
print(file_coord)

translate = np.array([11.03428, 18.96513, 17.05978])

# Translate atoms
for i in range(0, num_atoms):
    file_coord['X'][i] = file_coord['X'][i] + translate[0]
    file_coord['Y'][i] = file_coord['Y'][i] + translate[1]
    file_coord['Z'][i] = file_coord['Z'][i] + translate[2]

# Print to file
file_coord.insert(loc=0, column='A', value=pd.Series(species).values)
print_xyz.print_from_pandas(file_coord, num_atoms, '{}/{}'.format(folder, output_filename))

