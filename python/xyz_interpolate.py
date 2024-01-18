from __future__ import division, print_function
import pandas as pd
import numpy as np
import glob
from general import load_coordinates
from general import parameters as param
from general import print_xyz
import matplotlib.pyplot as plt
import csv
import xyz_siesta_to_cp2k

"""
    Interpolate between given .xyz files
"""

folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/melamine/cp2k-smeagol/cu/structures'
folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/melamine/cp2k-smeagol/cu/structures/5A-tip/geo_opt'
filename_input_1 = 'gs.xyz'
filename_input_2 = 'c1.xyz'
# filename_input_1 = 'c1.xyz'
# filename_input_2 = 'c2.xyz'
filename_output = 'ts1-linear.xyz'
# filename_output = 'ts2.xyz'
cols = ['Species', 'X', 'Y', 'Z']
file_coord_1, num_atoms_1, species_1 = load_coordinates.load_file_coord(folder, filename_input_1, cols)
file_coord_2, num_atoms_2, species_2 = load_coordinates.load_file_coord(folder, filename_input_2, cols)
print(file_coord_1)

interpolated = 0.5 * file_coord_1 + 0.5 * file_coord_2
interpolated.insert(loc=0, column='Species', value=species_1)
print_xyz.print_from_pandas(interpolated, num_atoms_1, '{}/{}'.format(folder, filename_output))

# interpolated = file_coord1
# interpolated = 0.5 * file_coord1 + 0.5 * file_coord2
# interpolated.insert(loc=0, column='Species', value=species1)
# print_xyz.print_from_pandas(interpolated, num_atoms1, '{}/{}'.format(folder1, output_filename))