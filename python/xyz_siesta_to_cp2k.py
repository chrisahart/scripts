from __future__ import division, print_function
import pandas as pd
import numpy as np
import glob
from general import load_coordinates
from general import print_xyz
import matplotlib.pyplot as plt
import csv

"""
    Convert XYZ to SIESTA %block AtomicCoordinatesAndAtomicSpecies
"""


def siesta_nolabel_to_cp2k(filename_input, filename_output, species_label):

    cols = ['X', 'Y', 'Z', 'Species_num']
    file_coord = pd.read_csv(filename_input, names=cols, delim_whitespace=True)
    num_atoms = file_coord.shape[0]
    species_number = file_coord['Species_num']
    file_coord = file_coord.apply(pd.to_numeric, errors='coerce')
    file_coord = file_coord.dropna(axis='rows', thresh=2)
    file_coord = file_coord.dropna(axis='columns', thresh=1)
    file_coord = file_coord.drop('Species_num', axis=1)

    species = []
    for i in range(num_atoms):
        species.append(species_label[species_number[i]])

    # Print to file
    file_coord.insert(loc=0, column='A', value=pd.Series(species).values)
    print(file_coord)
    file_coord = file_coord[file_coord['A'].notna()]
    file_coord = file_coord.dropna(axis='columns', thresh=1)
    num_atoms = file_coord.shape[0]
    print(file_coord)
    print_xyz.print_from_pandas(file_coord, num_atoms, filename_output)


def siesta_label_to_cp2k(filename_input, filename_output):

    cols = ['X', 'Y', 'Z', 'Species_num', 'Species', 'Num']
    file_coord = pd.read_csv(filename_input, names=cols, delim_whitespace=True)
    num_atoms = file_coord.shape[0]
    species = file_coord['Species']
    file_coord = file_coord.apply(pd.to_numeric, errors='coerce')
    file_coord = file_coord.dropna(axis='rows', thresh=2)
    file_coord = file_coord.dropna(axis='columns', thresh=1)
    file_coord = file_coord.drop('Num', axis=1)
    file_coord = file_coord.drop('Species_num', axis=1)

    # Print to file
    file_coord.insert(loc=0, column='A', value=pd.Series(species).values)
    print_xyz.print_from_pandas(file_coord, num_atoms, filename_output)


# folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/stm/clotilde/BIAS_NOSIC_FULL-Done'
# species_label = {1: 'Au_s', 2: 'Au', 3: 'C', 4: 'O', 5: 'H', 6: 'Si', 7: 'C_mol', 8: 'C_gra'}

# folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/stm/clotilde/BIAS_NOSIC-VACUUM'
# species_label = {1: 'Au', 2: 'Au', 3: np.NaN, 4: np.NaN, 5: np.NaN, 6: np.NaN, 7: 'C', 8: 'C'}

folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/ag-chain/clotilde/structures'

input_filename = 'em.siesta'
output_filename = 'em.xyz'
# input_filename = 'bulk.siesta'
# output_filename = 'bulk.xyz'

filename_input = '{}/{}'.format(folder, input_filename)
filename_output = '{}/{}'.format(folder, output_filename)

if __name__ == "__main__":
    siesta_label_to_cp2k(filename_input, filename_output)
    # siesta_nolabel_to_cp2k(filename_input, filename_output, species_label)
    print('Finished.')



