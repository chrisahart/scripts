from __future__ import division, print_function
import pandas as pd
import numpy as np
import glob
from general import load_coordinates
from general import print_xyz
import matplotlib.pyplot as plt
import math

"""
Convert .xyz to all NEGF (CP2K-SMEAGOL, CP2K-NEGF, SIESTA-SMEAGOL)    
"""

folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/au-capacitor-100/struct_9-8/geom_9-8/vesta'
input_filename = 'em_all_sorted.xyz'

filename_input = '{}/{}'.format(folder, input_filename)
cols = ['Species', 'X', 'Y', 'Z', 'Label']
file_coord, num_atoms, species = load_coordinates.load_file_coord(folder, input_filename, cols)
file_coord = file_coord.reset_index(drop=True)
print(file_coord)

repeats = [4.5, 1.5, 4]
repeats_negf = 4 * [1] + [0.5, 1.5] + 4 * [1]
num_atoms_slice = 24
dp = 6
dp_save = '%.6f'

# Output
output_filename_cp2k_vacuum = 'au-capacitor_cp2k.xyz'
output_filename_cp2k_dummy = 'au-capacitor_cp2k_dummy.xyz'
output_filename_cp2k_negf_vacuum = 'au-capacitor_cp2k_negf.xyz'
output_filename_cp2k_negf_dummy = 'au-capacitor_cp2k_negf_dummy.xyz'
output_filename_siesta = 'au-capacitor_siesta.xyz'
output_filename_siesta_dummy = 'au-capacitor_siesta_dummy.xyz'

system = file_coord.copy()
for i in range(0, system.shape[0]):
    system['X'][i] = np.format_float_positional(system['X'][i], unique=False, precision=dp)
    system['Y'][i] = np.format_float_positional(system['Y'][i], unique=False, precision=dp)
    system['Z'][i] = np.format_float_positional(system['Z'][i], unique=False, precision=dp)

species = int(repeats[0] * num_atoms_slice) * ['Au'] + \
          int(repeats[1] * num_atoms_slice) * ['X'] + \
          int(repeats[2] * num_atoms_slice) * ['Au']

species_siesta = int(repeats[0] * num_atoms_slice) * ['1'] + \
                  int(repeats[1] * num_atoms_slice) * ['2'] + \
                  int(repeats[2] * num_atoms_slice) * ['1']

species_negf = int(repeats_negf[0] * num_atoms_slice) * ['L3'] + int(repeats_negf[1] * num_atoms_slice) * ['L2'] + \
               int(repeats_negf[2] * num_atoms_slice) * ['L1'] + int(repeats_negf[3] * num_atoms_slice) * ['L0'] + \
               int(repeats_negf[4] * num_atoms_slice) * ['S'] + int(repeats_negf[5] * num_atoms_slice) * ['NaN'] + \
               int(repeats_negf[6] * num_atoms_slice) * ['R0'] + int(repeats_negf[7] * num_atoms_slice) * ['R1'] + \
               int(repeats_negf[8] * num_atoms_slice) * ['R2'] + int(repeats_negf[9] * num_atoms_slice) * ['R3']
species_negf_dummy = [w.replace('NaN', 'S') for w in species_negf]

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
cp2k_negf = system.copy()
cp2k_negf.insert(loc=0, column='A', value=species)
cp2k_negf = cp2k_negf.reset_index(drop=True)
cp2k_negf.insert(loc=4, column='B', value=species_negf)
cp2k_negf = cp2k_negf.drop(cp2k_negf[cp2k_negf.B == 'NaN'].index)
cp2k_negf = cp2k_negf.reset_index(drop=True)
print_xyz.print_from_pandas2(cp2k_negf, cp2k_negf.shape[0], '{}/{}'.format(folder, output_filename_cp2k_negf_vacuum), save_dp=dp_save)

# Print to file CP2K-NEGF with dummy atoms
cp2k_negf = system.copy()
cp2k_negf.insert(loc=0, column='A', value=species)
cp2k_negf = cp2k_negf.reset_index(drop=True)
cp2k_negf.insert(loc=4, column='B', value=species_negf_dummy)
print_xyz.print_from_pandas2(cp2k_negf, cp2k_negf.shape[0], '{}/{}'.format(folder, output_filename_cp2k_negf_dummy), save_dp=dp_save)

# Print to file SIESTA with vacuum
siesta = system.copy()
siesta.insert(loc=3, column='A', value=species_siesta)
siesta = siesta.reset_index(drop=True)
siesta.insert(loc=4, column='B', value=species)
print_xyz.print_from_pandas_siesta(siesta, siesta.shape[0], '{}/{}'.format(folder, output_filename_siesta), save_dp=dp_save)

# Print to file SIESTA with dummy atoms
siesta = system.copy()
siesta.insert(loc=3, column='A', value=species_siesta)
siesta.insert(loc=4, column='B', value=species)
siesta = siesta.reset_index(drop=True)
print_xyz.print_from_pandas_siesta(siesta, siesta.shape[0], '{}/{}'.format(folder, output_filename_siesta_dummy), save_dp=dp_save)
