from __future__ import division, print_function
import pandas as pd
import numpy as np
import glob
from general import load_coordinates
from general import print_xyz
import matplotlib.pyplot as plt
import math

"""
    create 1D chain .xyz 
"""

folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/iv/li/cp2k-negf/lda-q3_atoms-28/single-points'
output_filename_1 = 'input_cp2k.xyz'
output_filename_2 = 'input_cp2k_negf.xyz'
output_filename_3 = 'input_siesta.xyz'

num_atoms = 28
species = ['Au']
bond_length = 2.8

# Create dataframe
file_coord = pd.DataFrame(data={'X': [], 'Y': [], 'Z': []})

# Construct chain
for i in range(0, num_atoms):
    new_row = pd.Series({'X': 0,  'Y': 0, 'Z': 0 + i * bond_length})
    file_coord = pd.concat([file_coord, new_row.to_frame().T], ignore_index=True)

# Round
file_coord.X = file_coord.X.round(6)
file_coord.Y = file_coord.Y.round(6)
file_coord.Z = file_coord.Z.round(6)

# Print to file CP2K
cp2k = file_coord.copy()
cp2k.insert(loc=0, column='A', value=species*num_atoms)
print_xyz.print_from_pandas(cp2k, num_atoms, '{}/{}'.format(folder, output_filename_1))

# Print to file CP2K-NEGF
cp2k_negf = file_coord.copy()
labels = 4 * ['L2'] + 4 * ['L1'] + 4 * ['L0'] + 4 * ['S'] + 4 * ['R0'] + 4 * ['R1'] + 4 * ['R2']
cp2k_negf.insert(loc=0, column='A', value=species*num_atoms)
cp2k_negf = cp2k_negf.reset_index(drop=True)
cp2k_negf.insert(loc=4, column='B', value=labels)
print_xyz.print_from_pandas2(cp2k_negf, num_atoms, '{}/{}'.format(folder, output_filename_2))

# Print to file SIESTA
nums = np.linspace(start=1, stop=num_atoms, num=num_atoms, dtype=int)
siesta_species = num_atoms*[1]
siesta = file_coord.copy()
siesta.insert(loc=3, column='A', value=siesta_species)
siesta.insert(loc=4, column='B', value=species*num_atoms)
siesta.insert(loc=5, column='C', value=nums)
print(siesta)
print_xyz.print_from_pandas3(siesta, num_atoms, '{}/{}'.format(folder, output_filename_3))
