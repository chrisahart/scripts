from __future__ import division, print_function
import pandas as pd
import numpy as np
import glob
from general import load_coordinates
from general import print_xyz
import matplotlib.pyplot as plt

"""
    Sort .xyz sequentially along desired axes (useful for transport calculations)
"""

folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/structures/001_Au-S-2.75'
# input_filename = 'bulk.xyz'
# output_filename = 'bulk_sorted.xyz'
input_filename = 'em.xyz'
output_filename = 'em_sorted.xyz'
filename_output = '{}/{}'.format(folder, output_filename)

# Read number of atoms and labels from .xyz file
cols = ['Species', 'X', 'Y', 'Z']
file_coord, num_atoms, species = load_coordinates.load_file_coord(folder, input_filename, cols)
file_coord = file_coord.reset_index(drop=True)
print(file_coord)

file_coord.insert(loc=0, column='A', value=pd.Series(species).values)
# file_coord.sort_values(by=['Z'], inplace=True, ascending=[True])
file_coord.sort_values(by=['Z', 'X', 'Y'], inplace=True, ascending=[True, True, True])
file_coord = file_coord.reset_index(drop=True)
print(file_coord)

# Print to file
print_xyz.print_from_pandas(file_coord, num_atoms, filename_output)
