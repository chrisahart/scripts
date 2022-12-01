from __future__ import division, print_function
import pandas as pd
import numpy as np
import glob
from scripts.formatting import load_coordinates
from scripts.formatting import print_xyz
import matplotlib.pyplot as plt

"""
    .xyz translate
    Reads .xyz file and translates atomic positions
"""

folder = '/Volumes/Storage/Data/Work/Postdoc/Work/testing/cp2k-smeagol/bdt/structures/hollow-site/7-6'
# input_filename = 'C6H4S2.xyz'
# output_filename = 'C6H4S2-translated.xyz'
input_filename = '442-6.xyz'
output_filename = '442-6-translated.xyz'

# Read number of atoms and labels from .xyz file
cols = ['Species', 'X', 'Y', 'Z']
file_coord, num_atoms, species = load_coordinates.load_file_coord(folder, input_filename)

# translate = np.array([8.34258,  10.42822, 13.7994])
translate = np.array([0, 0, 21.4317])

# Translate atoms
for i in range(1, num_atoms+1):
    file_coord['X'][i] = file_coord['X'][i] + translate[0]
    file_coord['Y'][i] = file_coord['Y'][i] + translate[1]
    file_coord['Z'][i] = file_coord['Z'][i] + translate[2]

# Print to file
file_coord.insert(loc=0, column='A', value=pd.Series(species).values)
print_xyz.print_from_pandas(file_coord, num_atoms, '{}/{}'.format(folder, output_filename))

