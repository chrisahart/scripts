from __future__ import division, print_function
import pandas as pd
import numpy as np
import glob
from general import load_coordinates
from general import print_xyz
import matplotlib.pyplot as plt
from ase.io import read, write

"""
    Sort .xyz sequentially along desired axes (useful for transport calculations)
"""

# folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-capacitor/layers-1-2-3-4/cp2k/weightrho/testing/kpoints-4-4-20_hlb-auto_NEnergReal-64/input'
# folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-wire/jad/structures'
# input_filename = 'Au_Junc_wrapped.xyz'
# output_filename = 'Au_Junc_wrapped_sorted.xyz'

# folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-wire-ismael/jad/structures/jad_structures_replace-a'
# input_filename = 'structure_a.xyz'
# output_filename = 'structure_a_sorted.xyz'

folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-wire-ismael/jad/structures/jad_structures_replace-a'
input_filename = 'structure_em_label.xyz'
output_filename = 'structure_em_label_sorted_wire-only.xyz'

filename_output = '{}/{}'.format(folder, output_filename)

# Load Jad structure
jad_structure = read('{}/{}'.format(folder, input_filename), format='xyz')
print(jad_structure)

# Read number of atoms and labels from .xyz file
# cols = ['Species', 'X', 'Y', 'Z']
# file_coord, num_atoms, species = load_coordinates.load_file_coord(folder, input_filename, cols)
# file_coord = file_coord.reset_index(drop=True)
# print(file_coord)

# print(file_coord.shape[0])
# for i in range(file_coord.shape[0]):
#     if file_coord['Z'][i] > 20:
#         file_coord['Z'][i] = file_coord['Z'][i]+20

# Delete species
# file_coord= file_coord[file_coord['Species'] != 'Au_al']
#
# file_coord.insert(loc=0, column='A', value=pd.Series(species).values)
# # file_coord.sort_values(by=['Z'], inplace=True, ascending=[True])
# file_coord.sort_values(by=['Z', 'X', 'Y'], inplace=True, ascending=[True, True, True])
# # file_coord.sort_values(by=['Y', 'X', 'Z'], inplace=True, ascending=[True, True, True])
# file_coord = file_coord.reset_index(drop=True)
# print(file_coord)
#
# # Print to file
# print_xyz.print_from_pandas(file_coord, num_atoms, filename_output, save_dp='%.4f')
