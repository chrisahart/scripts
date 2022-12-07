from __future__ import division, print_function
import pandas as pd
import numpy as np
import glob
from general import load_coordinates
from general import print_xyz
import matplotlib.pyplot as plt
import math

"""
    
"""
dp = 4

# Au capacitor
folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/au-capacitor/structures'
input_filename_1 = 'supercell_3x3-6.xyz'
output_filename_1 = 'au-capacitor_cp2k.xyz'
output_filename_2 = 'au-capacitor_cp2k_dummy.xyz'
output_filename_3 = 'au-capacitor_cp2k_negf_dummy.xyz'
output_filename_4 = 'au-capacitor_siesta_dummy.xyz'

cols = ['Species', 'X', 'Y', 'Z']
file_coord_1, num_atoms_1, species_1 = load_coordinates.load_file_coord(folder, input_filename_1, cols)
file_coord_1 = file_coord_1.reset_index(drop=True)
species_1 = species_1.reset_index(drop=True)
print(file_coord_1.shape[0])

# Left
max_z_left = 19
file_coord_1 = file_coord_1.drop(file_coord_1[file_coord_1['Z'] > max_z_left].index)

# Right
max_z_end = 3
file_coord_2 = file_coord_1.copy()
file_coord_2 = file_coord_2.drop(file_coord_2[file_coord_2['Z'] > max_z_end].index)
translate_right = np.array([0, 0, 2.08560*10])
print(file_coord_2.shape[0])
for i in range(0, file_coord_2.shape[0]):
    file_coord_2['X'][i] = file_coord_2['X'][i] + translate_right[0]
    file_coord_2['Y'][i] = file_coord_2['Y'][i] + translate_right[1]
    file_coord_2['Z'][i] = file_coord_2['Z'][i] + translate_right[2]

# Join all
system = pd.concat([file_coord_1, file_coord_2], ignore_index=True, sort=False)
system = system.reset_index(drop=True)
system.X = system.X.round(dp)
system.Y = system.Y.round(dp)
system.Z = system.Z.round(dp)
num_atoms_1 = system.shape[0]
print(system.shape[0])

# Print to file CP2K normal
cp2k = system.copy()
species = 5 * 18 * ['Au'] + 3 * 18 * ['Au'] + 4 * 18 * ['Au']
cp2k.insert(loc=0, column='A', value=species)
print_xyz.print_from_pandas(cp2k, num_atoms_1, '{}/{}'.format(folder, output_filename_1))

# Print to file CP2K with dummy atoms
species = 5 * 18 * ['Au'] + 3 * 18 * ['X'] + 4 * 18 * ['Au']
cp2k = system.copy()
cp2k.insert(loc=0, column='A', value=species)
print(cp2k)
print_xyz.print_from_pandas(cp2k, num_atoms_1, '{}/{}'.format(folder, output_filename_2))

# Print to file CP2K-NEGF
cp2k_negf = system.copy()
species = 5 * 18 * ['Au'] + 3 * 18 * ['X'] + 4 * 18 * ['Au']
labels = 1 * 18 * ['L2'] + 1 * 18 * ['L1'] + 1 * 18 * ['L0'] + 6 * 18 * ['S'] + 1 * 18 * ['R0'] + 1 * 18 * ['R1'] + 1 * 18 * ['R2']
cp2k_negf.insert(loc=0, column='A', value=species)
cp2k_negf = cp2k_negf.reset_index(drop=True)
cp2k_negf.insert(loc=4, column='B', value=labels)
print_xyz.print_from_pandas2(cp2k_negf, num_atoms_1, '{}/{}'.format(folder, output_filename_3))

# Print to file SIESTA
siesta_species = 5 * 18 * ['1'] + 3 * 18 * ['2'] + 4 * 18 * ['1']
species = 5 * 18 * ['Au'] + 3 * 18 * ['X'] + 4 * 18 * ['Au']
nums = np.linspace(start=1, stop=num_atoms_1, num=num_atoms_1, dtype=int)
siesta = system.copy()
siesta.insert(loc=3, column='A', value=siesta_species)
siesta.insert(loc=4, column='B', value=species)
siesta.insert(loc=5, column='C', value=nums)
print_xyz.print_from_pandas3(siesta, num_atoms_1, '{}/{}'.format(folder, output_filename_4))
