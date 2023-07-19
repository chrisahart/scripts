from __future__ import division, print_function
import pandas as pd
import numpy as np
import glob
from general import load_coordinates
from general import print_xyz
import matplotlib.pyplot as plt
import math

"""
    Create junction
"""

folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/mengxuan/structures/ben_ant/structures/chris'
left = 'left.xyz'
molecule = 'molecule.xyz'
print_em = 'em.xyz'
max_z = np.NaN
x1 = 4.17129
y1 = 2.08564
y2 = 4.17129
Au_S = 2.45
au_au_z = 2.08564
dp_print = 3

# left
cols = ['Species', 'X', 'Y', 'Z']
file_coord_1, num_atoms_1, species_1 = load_coordinates.load_file_coord(folder, left, cols)
file_coord_1 = file_coord_1.reset_index(drop=True)

# molecule
file_coord_2, num_atoms_2, species_2 = load_coordinates.load_file_coord(folder, molecule, cols)
file_coord_2 = file_coord_2.reset_index(drop=True)

# right
file_coord_3 = file_coord_1.copy()
file_coord_3 = file_coord_3.drop(file_coord_3[file_coord_3['Z'] > max_z].index)
file_coord_3 = file_coord_3.reset_index(drop=True)
num_atoms_3 = len(file_coord_3)
species_3 = species_1[0:num_atoms_3]

# Reset molecule
translate_reset = -1 * np.array([0, 0, np.min(file_coord_2['Z'])])
# translate_reset = np.array([file_coord_2['X'][0], file_coord_2['Y'][0], file_coord_2['Z'][0]])
for i in range(0, num_atoms_2):
    file_coord_2['X'][i] = file_coord_2['X'][i] + translate_reset[0]
    file_coord_2['Y'][i] = file_coord_2['Y'][i] + translate_reset[1]
    file_coord_2['Z'][i] = file_coord_2['Z'][i] + translate_reset[2]

# Translate molecule
z2 = np.max(file_coord_1['Z'])
z1 = z2 + 2.76
# z1 = np.sqrt((Au_S ** 2) - (y2 - y1) ** 2) + z2
translate_molecule = np.array([0, 0, z1])
for i in range(0, num_atoms_2):
    file_coord_2['X'][i] = file_coord_2['X'][i] + translate_molecule[0]
    file_coord_2['Y'][i] = file_coord_2['Y'][i] + translate_molecule[1]
    file_coord_2['Z'][i] = file_coord_2['Z'][i] + translate_molecule[2]

# Translate right
z2 = np.max(file_coord_2['Z'])
z1 = z2 + 2.76
# z1 = np.sqrt((Au_S ** 2) - (y2 - y1) ** 2) + z2
translate_right = np.array([0, 0, z1])
for i in range(0, num_atoms_3):
    file_coord_3['X'][i] = file_coord_3['X'][i] + translate_right[0]
    file_coord_3['Y'][i] = file_coord_3['Y'][i] + translate_right[1]
    file_coord_3['Z'][i] = file_coord_3['Z'][i] + translate_right[2]

# Join coordinates
left_molecule = pd.concat([file_coord_1, file_coord_2], ignore_index=True, sort=False)
left_molecule_right = pd.concat([left_molecule, file_coord_3], ignore_index=True, sort=False)
left_molecule_right = left_molecule_right.reset_index(drop=True)

# Join coordinates
species_12 = pd.concat([species_1, species_2], ignore_index=True, sort=False)
species_123 = pd.concat([species_12, species_3], ignore_index=True, sort=False)
species_123 = species_123.reset_index(drop=True)

# Print em to file
num_atoms = left_molecule_right.shape[0]
left_molecule_right = left_molecule_right.round(dp_print)
left_molecule_right.insert(loc=0, column='A', value=pd.Series(species_123).values)
print_xyz.print_from_pandas(left_molecule_right, num_atoms, '{}/{}'.format(folder, print_em))

# Print cell information
print('num_atoms', num_atoms)
# print('bulk cell size', '{0:.3f}'.format(au_au_z*(bulk_size[0]+1)),
#       '{0:.3f}'.format(au_au_z*(bulk_size[1]+1)), '{0:.3f}'.format(au_au_z*(bulk_size[2]+1)))
# print('em num_atoms', num_atoms)
# print('em cell size', '{0:.3f}'.format(au_au_z*(em_size[0]+1)),
#       '{0:.3f}'.format(au_au_z*(em_size[1]+1)), '{0:.3f}'.format(au_au_z+(em_size[2]+1)))
