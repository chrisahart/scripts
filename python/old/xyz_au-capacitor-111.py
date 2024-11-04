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

# Au capacitor 111
folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/au-capacitor-111/structures/chris-struct'
input_filename_1 = 'sergey_left.xyz'
input_filename_2 = 'sergey_right.xyz'
repeats_left = 4
repeats_right = 4
au_au_z = 2.38417
layers = 3
layers_electrode = 4
dp = 6

# Output
output_filename_cp2k_vacuum = 'au-capacitor_cp2k.xyz'
output_filename_cp2k_dummy = 'au-capacitor_cp2k_dummy.xyz'
output_filename_cp2k_negf_vacuum = 'au-capacitor_cp2k_negf.xyz'
output_filename_cp2k_negf_dummy = 'au-capacitor_cp2k_negf_dummy.xyz'
output_filename_siesta = 'au-capacitor_siesta.xyz'
output_filename_siesta_dummy = 'au-capacitor_siesta_dummy.xyz'

cols = ['Species', 'X', 'Y', 'Z']
unit_1, num_atoms_1, species_1 = load_coordinates.load_file_coord(folder, input_filename_1, cols)
unit_1 = unit_1.reset_index(drop=True)
species_1 = species_1.reset_index(drop=True)
unit_2, num_atoms_2, species_2 = load_coordinates.load_file_coord(folder, input_filename_2, cols)
unit_2 = unit_2.reset_index(drop=True)
species_2 = species_2.reset_index(drop=True)
print(num_atoms_1, num_atoms_2)
num_atoms_slice = num_atoms_1

# System left including dummy atoms
system_left = unit_1.copy()
for n in range(1, repeats_left + 1):
    translate = [0, 0, au_au_z * layers * n]
    unit_1_repeat = unit_1.copy()
    for i in range(0, unit_1_repeat.shape[0]):
        unit_1_repeat['X'][i] = unit_1_repeat['X'][i] + translate[0]
        unit_1_repeat['Y'][i] = unit_1_repeat['Y'][i] + translate[1]
        unit_1_repeat['Z'][i] = unit_1_repeat['Z'][i] + translate[2]
        unit_1_repeat['X'][i] = np.format_float_positional(unit_1_repeat['X'][i], unique=False, precision=dp)
        unit_1_repeat['Y'][i] = np.format_float_positional(unit_1_repeat['Y'][i], unique=False, precision=dp)
        unit_1_repeat['Z'][i] = np.format_float_positional(unit_1_repeat['Z'][i], unique=False, precision=dp)
    system_left = pd.concat([system_left, unit_1_repeat], ignore_index=True, sort=False)
    system_left = system_left.reset_index(drop=True)
num_atoms_system_left = system_left.shape[0]

# Reset right position to z=0
unit_2_reset = unit_2.copy()
translate_reset = [0, 0, np.min(unit_2_reset['Z'])]
# print('before translate', unit_2_reset)
for i in range(0, num_atoms_2):
    unit_2_reset['X'][i] = unit_2_reset['X'][i] - translate_reset[0]
    unit_2_reset['Y'][i] = unit_2_reset['Y'][i] - translate_reset[1]
    unit_2_reset['Z'][i] = unit_2_reset['Z'][i] - translate_reset[2]

# System right including dummy atoms
system_right = pd.DataFrame(columns=['X', 'Y', 'Z'])
for n in range(repeats_left+1, repeats_left + repeats_right + 1):
    translate = [0, 0, au_au_z * layers * n]
    unit_2_repeat = unit_2_reset.copy()
    for i in range(0, unit_2_repeat.shape[0]):
        print(translate[2])
        unit_2_repeat['X'][i] = unit_2_repeat['X'][i] + translate[0]
        unit_2_repeat['Y'][i] = unit_2_repeat['Y'][i] + translate[1]
        unit_2_repeat['Z'][i] = unit_2_repeat['Z'][i] + translate[2]
        unit_2_repeat['X'][i] = np.format_float_positional(unit_2_repeat['X'][i], unique=False, precision=dp)
        unit_2_repeat['Y'][i] = np.format_float_positional(unit_2_repeat['Y'][i], unique=False, precision=dp)
        unit_2_repeat['Z'][i] = np.format_float_positional(unit_2_repeat['Z'][i], unique=False, precision=dp)
    system_right = pd.concat([system_right, unit_2_repeat], ignore_index=True, sort=False)
    system_right = system_right.reset_index(drop=True)
num_atoms_system_right = system_right.shape[0]
system = pd.concat([system_left, system_right], ignore_index=True, sort=False)
system = system.reset_index(drop=True)
num_atoms_1 = system.shape[0]

# Print to file CP2K vacuum
cp2k = system.copy()
species = 4 * num_atoms_slice * ['Au'] + 1 * num_atoms_slice * ['NaN'] + 4 * num_atoms_slice * ['Au']
cp2k.insert(loc=0, column='A', value=species)
cp2k = cp2k.drop(cp2k[cp2k.A == 'NaN'].index)
print_xyz.print_from_pandas(cp2k, cp2k.shape[0], '{}/{}'.format(folder, output_filename_cp2k_vacuum), save_dp='%.6f')

# Print to file CP2K with dummy atoms
species = 4 * num_atoms_slice * ['Au'] + 1 * num_atoms_slice * ['X'] + 4 * num_atoms_slice * ['Au']
cp2k = system.copy()
cp2k.insert(loc=0, column='A', value=species)
print_xyz.print_from_pandas(cp2k, cp2k.shape[0], '{}/{}'.format(folder, output_filename_cp2k_dummy), save_dp='%.6f')

# Print to file CP2K-NEGF with vacuum
cp2k_negf = system.copy()
species = 4 * num_atoms_slice * ['Au'] + 1 * num_atoms_slice * ['NaN'] + 4 * num_atoms_slice * ['Au']
labels = 1 * num_atoms_slice * ['L3'] + 1 * num_atoms_slice * ['L2'] + 1 * num_atoms_slice * ['L1'] + 1 * num_atoms_slice * ['L0'] + \
         1 * num_atoms_slice * ['NaN'] + \
         1 * num_atoms_slice * ['R0'] + 1 * num_atoms_slice * ['R1'] + 1 * num_atoms_slice * ['R2'] + 1 * num_atoms_slice * ['R3']
cp2k_negf.insert(loc=0, column='A', value=species)
cp2k_negf = cp2k_negf.reset_index(drop=True)
cp2k_negf.insert(loc=4, column='B', value=labels)
cp2k_negf = cp2k_negf.drop(cp2k_negf[cp2k_negf.A == 'NaN'].index)
cp2k_negf = cp2k_negf.reset_index(drop=True)
print_xyz.print_from_pandas2(cp2k_negf, cp2k_negf.shape[0], '{}/{}'.format(folder, output_filename_cp2k_negf_vacuum), save_dp='%.6f')

# Print to file CP2K-NEGF with dummy atoms
cp2k_negf = system.copy()
species = 4 * num_atoms_slice * ['Au'] + 1 * num_atoms_slice * ['X'] + 4 * num_atoms_slice * ['Au']
labels = 1 * num_atoms_slice * ['L3'] + 1 * num_atoms_slice * ['L2'] + 1 * num_atoms_slice * ['L1'] + 1 * num_atoms_slice * ['L0'] + \
         1 * num_atoms_slice * ['S'] + \
         1 * num_atoms_slice * ['R0'] + 1 * num_atoms_slice * ['R1'] + 1 * num_atoms_slice * ['R2'] + 1 * num_atoms_slice * ['R3']
cp2k_negf.insert(loc=0, column='A', value=species)
cp2k_negf = cp2k_negf.reset_index(drop=True)
cp2k_negf.insert(loc=4, column='B', value=labels)
print_xyz.print_from_pandas2(cp2k_negf, cp2k_negf.shape[0], '{}/{}'.format(folder, output_filename_cp2k_negf_dummy), save_dp='%.6f')

# Print to file SIESTA with vacuum
siesta = system.copy()
siesta_species = 4 * num_atoms_slice * ['1'] + 1 * num_atoms_slice * ['2'] + 4 * num_atoms_slice * ['1']
species = 4 * num_atoms_slice * ['Au'] + 1 * num_atoms_slice * ['X'] + 4 * num_atoms_slice * ['Au']
nums = np.linspace(start=1, stop=siesta.shape[0], num=siesta.shape[0], dtype=int)
siesta.insert(loc=3, column='A', value=siesta_species)
siesta.insert(loc=4, column='B', value=species)
siesta.insert(loc=5, column='C', value=nums)
siesta = siesta.drop(siesta[siesta.B == 'X'].index)
siesta = siesta.reset_index(drop=True)
print(siesta)
print_xyz.print_from_pandas3(siesta, siesta.shape[0], '{}/{}'.format(folder, output_filename_siesta), save_dp='%.6f')

# Print to file SIESTA with dummy atoms
siesta = system.copy()
siesta_species = 4 * num_atoms_slice * ['1'] + 1 * num_atoms_slice * ['2'] + 4 * num_atoms_slice * ['1']
species = 4 * num_atoms_slice * ['Au'] + 1 * num_atoms_slice * ['X'] + 4 * num_atoms_slice * ['Au']
nums = np.linspace(start=1, stop=siesta.shape[0], num=siesta.shape[0], dtype=int)
siesta.insert(loc=3, column='A', value=siesta_species)
siesta.insert(loc=4, column='B', value=species)
siesta.insert(loc=5, column='C', value=nums)
print_xyz.print_from_pandas3(siesta, siesta.shape[0], '{}/{}'.format(folder, output_filename_siesta_dummy), save_dp='%.6f')
