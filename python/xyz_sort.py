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

# Files
# folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/au-capacitor-100/struct_9-8/geom_9-8/vesta'
# input_filename = 'em_all.xyz'
# output_filename = 'em_all_sorted.xyz'

# folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Masters/2021-2022/Yike/Pt3Ni/yike_struct'
# input_filename = 'input.xyz'
# output_filename = 'input_sorted.xyz'

# folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Masters/2021-2022/Yike/Pt3Ni/chris_struct/interface_001/brute_force/scf-1e-2'
folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Masters/2021-2022/Yike/Pt3Ni/delete/yike_struct/yike_struct_fm'
input_filename = 'input.xyz'
output_filename = 'input_sorted.xyz'
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
