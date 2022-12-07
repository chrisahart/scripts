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
folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/au-capacitor-111/structures'
input_filename_1 = 'sergey_unit.xyz'
repeats = 8
au_au_z = 2.38417
translate_right = 1.68586
layers = 3
layers_electrode = 4
dp = 6

# Output
output_filename_cp2k_vacuum = 'au-capacitor_cp2k.xyz'
output_filename_cp2k_dummy = 'au-capacitor_cp2k_dummy.xyz'
output_filename_cp2k_negf_vacuum = 'au-capacitor_cp2k_negf.xyz'
output_filename_cp2k_negf_dummy = 'au-capacitor_cp2k_negf_dummy.xyz'
output_filename_siesta = 'au-capacitor_siesta_dummy.xyz'

cols = ['Species', 'X', 'Y', 'Z']
unit_1, num_atoms_1, species_1 = load_coordinates.load_file_coord(folder, input_filename_1, cols)
unit_1 = unit_1.reset_index(drop=True)
species_1 = species_1.reset_index(drop=True)

# System
system = unit_1.copy()
for n in range(1, repeats + 1):
    translate = [0, 0, au_au_z * layers * n]
    unit_1_repeat = unit_1.copy()
    for i in range(0, unit_1_repeat.shape[0]):
        unit_1_repeat['X'][i] = unit_1_repeat['X'][i] + translate[0]
        unit_1_repeat['Y'][i] = unit_1_repeat['Y'][i] + translate[1]
        unit_1_repeat['Z'][i] = unit_1_repeat['Z'][i] + translate[2]
        unit_1_repeat['X'][i] = np.format_float_positional(unit_1_repeat['X'][i], unique=False, precision=dp)
        unit_1_repeat['Y'][i] = np.format_float_positional(unit_1_repeat['Y'][i], unique=False, precision=dp)
        unit_1_repeat['Z'][i] = np.format_float_positional(unit_1_repeat['Z'][i], unique=False, precision=dp)
    system = pd.concat([system, unit_1_repeat], ignore_index=True, sort=False)
    system = system.reset_index(drop=True)
num_atoms_system = system.shape[0]

# Translate right
translate = [0, translate_right, 0]
right = num_atoms_1 * layers_electrode
for i in range(right, system.shape[0]):
    system['X'][i] = system['X'][i] + translate[0]
    system['Y'][i] = system['Y'][i] + translate[1]
    system['Z'][i] = system['Z'][i] + translate[2]

# # Print to file CP2K vacuum
cp2k = system.copy()
species = 4 * num_atoms_1 * ['Au'] + 1 * num_atoms_1 * ['NaN'] + 4 * num_atoms_1 * ['Au']
cp2k.insert(loc=0, column='A', value=species)
cp2k = cp2k.drop(cp2k[cp2k.A == 'NaN'].index)
print_xyz.print_from_pandas(cp2k, cp2k.shape[0], '{}/{}'.format(folder, output_filename_cp2k_vacuum), save_dp='%.6f')

# Print to file CP2K with dummy atoms
species = 4 * num_atoms_1 * ['Au'] + 1 * num_atoms_1 * ['X'] + 4 * num_atoms_1 * ['Au']
cp2k = system.copy()
cp2k.insert(loc=0, column='A', value=species)
print_xyz.print_from_pandas(cp2k, cp2k.shape[0], '{}/{}'.format(folder, output_filename_cp2k_dummy), save_dp='%.6f')

# Print to file CP2K-NEGF with vacuum
cp2k_negf = system.copy()
species = 4 * num_atoms_1 * ['Au'] + 1 * num_atoms_1 * ['NaN'] + 4 * num_atoms_1 * ['Au']
labels = 1 * num_atoms_1 * ['L3'] + 1 * num_atoms_1 * ['L2'] + 1 * num_atoms_1 * ['L1'] + 1 * num_atoms_1 * ['L0'] + \
         1 * num_atoms_1 * ['NaN'] + \
         1 * num_atoms_1 * ['R0'] + 1 * num_atoms_1 * ['R1'] + 1 * num_atoms_1 * ['R2'] + 1 * num_atoms_1 * ['R3']
cp2k_negf.insert(loc=0, column='A', value=species)
cp2k_negf = cp2k_negf.reset_index(drop=True)
cp2k_negf.insert(loc=4, column='B', value=labels)
cp2k_negf = cp2k_negf.drop(cp2k_negf[cp2k_negf.A == 'NaN'].index)
cp2k_negf = cp2k_negf.reset_index(drop=True)
print_xyz.print_from_pandas2(cp2k_negf, cp2k_negf.shape[0], '{}/{}'.format(folder, output_filename_cp2k_negf_vacuum), save_dp='%.6f')

# Print to file CP2K-NEGF with dummy atoms
cp2k_negf = system.copy()
species = 4 * num_atoms_1 * ['Au'] + 1 * num_atoms_1 * ['X'] + 4 * num_atoms_1 * ['Au']
labels = 1 * num_atoms_1 * ['L3'] + 1 * num_atoms_1 * ['L2'] + 1 * num_atoms_1 * ['L1'] + 1 * num_atoms_1 * ['L0'] + \
         1 * num_atoms_1 * ['S'] + \
         1 * num_atoms_1 * ['R0'] + 1 * num_atoms_1 * ['R1'] + 1 * num_atoms_1 * ['R2'] + 1 * num_atoms_1 * ['R3']
cp2k_negf.insert(loc=0, column='A', value=species)
cp2k_negf = cp2k_negf.reset_index(drop=True)
cp2k_negf.insert(loc=4, column='B', value=labels)
print_xyz.print_from_pandas2(cp2k_negf, cp2k_negf.shape[0], '{}/{}'.format(folder, output_filename_cp2k_negf_dummy), save_dp='%.6f')

#
# # Print to file SIESTA
# siesta_species = 5 * 18 * ['1'] + 3 * 18 * ['2'] + 4 * 18 * ['1']
# species = 5 * 18 * ['Au'] + 3 * 18 * ['X'] + 4 * 18 * ['Au']
# nums = np.linspace(start=1, stop=num_atoms_1, num=num_atoms_1, dtype=int)
# siesta = system.copy()
# siesta.insert(loc=3, column='A', value=siesta_species)
# siesta.insert(loc=4, column='B', value=species)
# siesta.insert(loc=5, column='C', value=nums)
# print_xyz.print_from_pandas3(siesta, num_atoms_1, '{}/{}'.format(folder, output_filename_4))
