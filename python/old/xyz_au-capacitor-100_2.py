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
folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/cp2k-smeagol/sz/transmission/exp/capacitor/delete/geom'
input_filename_1 = 'layer.xyz'
repeats = [4.5, 1.5, 4]
repeats_left = int(np.sum(repeats)-1)
au_au_z = 2.08600
au_au_xy = 4.1710
layers = 2
dp = 6
dp_save = '%.6f'

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
num_atoms_slice = num_atoms_1

species = int(repeats[0] * num_atoms_slice) * ['Au'] + \
          int(repeats[1] * num_atoms_slice) * ['X'] + \
          int(repeats[2] * num_atoms_slice) * ['Au']

species_species = int(repeats[0] * num_atoms_slice) * ['1'] + \
                  int(repeats[1] * num_atoms_slice) * ['2'] + \
                  int(repeats[2] * num_atoms_slice) * ['1']

negf_species = 1 * num_atoms_slice * ['L3'] + 1 * num_atoms_slice * ['L2'] + 1 * num_atoms_slice * ['L1'] + 1 * num_atoms_slice * ['L0'] + \
         1 * num_atoms_slice * ['S'] + \
         1 * num_atoms_slice * ['R0'] + 1 * num_atoms_slice * ['R1'] + 1 * num_atoms_slice * ['R2'] + 1 * num_atoms_slice * ['R3']

# System left including dummy atoms
system = unit_1.copy()
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
    system = pd.concat([system, unit_1_repeat], ignore_index=True, sort=False)
    system = system.reset_index(drop=True)
num_atoms_system = system.shape[0]

# Print to file CP2K vacuum
cp2k = system.copy()
cp2k.insert(loc=0, column='A', value=species)
cp2k = cp2k.reset_index(drop=True)
cp2k = cp2k.drop(cp2k[cp2k.A == 'X'].index)
print_xyz.print_from_pandas(cp2k, cp2k.shape[0], '{}/{}'.format(folder, output_filename_cp2k_vacuum), save_dp=dp_save)

# Print to file CP2K with dummy atoms
cp2k = system.copy()
cp2k.insert(loc=0, column='A', value=species)
cp2k = cp2k.reset_index(drop=True)
print_xyz.print_from_pandas(cp2k, cp2k.shape[0], '{}/{}'.format(folder, output_filename_cp2k_dummy), save_dp=dp_save)

# Print to file CP2K-NEGF with vacuum
# cp2k_negf = system.copy()
# labels = 1 * num_atoms_slice * ['L3'] + 1 * num_atoms_slice * ['L2'] + 1 * num_atoms_slice * ['L1'] + 1 * num_atoms_slice * ['L0'] + \
#          1 * num_atoms_slice * ['NaN'] + \
#          1 * num_atoms_slice * ['R0'] + 1 * num_atoms_slice * ['R1'] + 1 * num_atoms_slice * ['R2'] + 1 * num_atoms_slice * ['R3']
# cp2k_negf.insert(loc=0, column='A', value=species)
# cp2k_negf = cp2k_negf.reset_index(drop=True)
# cp2k_negf.insert(loc=4, column='B', value=labels)
# cp2k_negf = cp2k_negf.drop(cp2k_negf[cp2k_negf.A == 'NaN'].index)
# cp2k_negf = cp2k_negf.reset_index(drop=True)
# print_xyz.print_from_pandas2(cp2k_negf, cp2k_negf.shape[0], '{}/{}'.format(folder, output_filename_cp2k_negf_vacuum), save_dp=dp_save)

# Print to file CP2K-NEGF with dummy atoms
# cp2k_negf = system.copy()
# cp2k_negf.insert(loc=0, column='A', value=species)
# cp2k_negf = cp2k_negf.reset_index(drop=True)
# cp2k_negf.insert(loc=4, column='B', value=labels)
# print_xyz.print_from_pandas2(cp2k_negf, cp2k_negf.shape[0], '{}/{}'.format(folder, output_filename_cp2k_negf_dummy), save_dp=dp_save)

# Print to file SIESTA with vacuum
siesta = system.copy()
siesta.insert(loc=3, column='A', value=species_species)
siesta.insert(loc=4, column='B', value=species)
siesta = siesta.drop(siesta[siesta.B == 'X'].index)
# nums = np.linspace(start=1, stop=siesta.shape[0], num=siesta.shape[0], dtype=int)
# siesta.insert(loc=5, column='C', value=nums)
print(siesta)
print_xyz.print_from_pandas2(siesta, siesta.shape[0], '{}/{}'.format(folder, output_filename_siesta), save_dp=dp_save)

# Print to file SIESTA with dummy atoms
siesta = system.copy()
siesta.insert(loc=3, column='A', value=species_species)
siesta.insert(loc=4, column='B', value=species)
# nums = np.linspace(start=1, stop=siesta.shape[0], num=siesta.shape[0], dtype=int)
# siesta.insert(loc=5, column='C', value=nums)
print_xyz.print_from_pandas2(siesta, siesta.shape[0], '{}/{}'.format(folder, output_filename_siesta_dummy), save_dp=dp_save)

print('cell:', au_au_xy*2, au_au_xy*2, au_au_z*(repeats_left+1)*layers)
