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


def cp2k_to_siesta(input_filename, filename_output, siesta_species_to_atomic):

    cols = ['Species', 'X', 'Y', 'Z', 'Label']
    file_coord, num_atoms, species = load_coordinates.load_file_coord('', input_filename, cols)

    # Replace species with value
    species_val = species.copy()
    for i in range(0, num_atoms):
        species_val[i] = siesta_species_to_atomic[species[i]]

    # Print to file
    file_coord.insert(loc=3, column='A', value=pd.Series(species_val).values)
    file_coord.insert(loc=4, column='B', value=pd.Series(species).values)
    file_coord.to_csv(filename_output, index=False, header=False, quoting=csv.QUOTE_NONE, sep=" ", float_format='%.8f')


def siesta_nolabel_to_cp2k(filename_input, filename_output, cols, species_label, num_atoms, image):

    file_coord = pd.read_csv(filename_input, names=cols, delim_whitespace=True)
    file_coord = file_coord[num_atoms * (image - 1) + 1 * image:num_atoms * image + 1 * image]

    file_coord = file_coord.reset_index(drop=True)
    species_number = file_coord['Species']
    file_coord = file_coord.apply(pd.to_numeric, errors='coerce')
    file_coord = file_coord.dropna(axis='rows', thresh=2)
    file_coord = file_coord.dropna(axis='columns', thresh=1)
    file_coord = file_coord[['X', 'Y', 'Z']]

    species = []
    for i in range(num_atoms):
        species.append(species_label[int(species_number[i])])

    # Print to file
    file_coord.insert(loc=0, column='A', value=pd.Series(species).values)
    file_coord = file_coord[file_coord['A'].notna()]
    file_coord = file_coord.dropna(axis='columns', thresh=1)
    num_atoms = file_coord.shape[0]
    print_xyz.print_from_pandas(file_coord, num_atoms, filename_output, save_dp='%.8f')


def siesta_label_to_cp2k(filename_input, filename_output, cols, num_atoms, image):

    file_coord = pd.read_csv(filename_input, names=cols, delim_whitespace=True)
    file_coord = file_coord[num_atoms * (image - 1) + 1 * image:num_atoms * image + 1 * image]

    species = file_coord['Species']
    file_coord = file_coord.apply(pd.to_numeric, errors='coerce')
    file_coord = file_coord.dropna(axis='rows', thresh=2)
    file_coord = file_coord.dropna(axis='columns', thresh=1)
    file_coord = file_coord[['X', 'Y', 'Z']]

    # Print to file
    file_coord.insert(loc=0, column='A', value=pd.Series(species).values)
    print_xyz.print_from_pandas(file_coord, num_atoms, filename_output, save_dp='%.8f')


folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/neb/siesta1-smeagol/archer/neb-tutorial/SMEAGOL_NEB/smeagol_example/r1/input'
cols = ['Species', 'X', 'Y', 'Z', 'A', 'B', 'C']
siesta_atomic_to_species = {79: 'Au', 7: 'N', 6: 'C', 1: 'H'}
# siesta_species_to_atomic = {'Au': 79, 'N': 6, 'C': 1, 'H': 1}
siesta_species_to_atomic = {'Au': 4, 'N': 3, 'C': 2, 'H': 1}
image = 11
num_atoms = 183
input_filename = 'melamine.NEB'
output_filename_cp2k = 'melamine-{}.xyz'.format(image)
output_filename_siesta = 'melamine-{}_au-same.siesta'.format(image)

filename_input = '{}/{}'.format(folder, input_filename)
filename_output_cp2k = '{}/{}'.format(folder, output_filename_cp2k)
filename_output_siesta = '{}/{}'.format(folder, output_filename_siesta)

if __name__ == "__main__":
    siesta_nolabel_to_cp2k(filename_input, filename_output_cp2k, cols, siesta_atomic_to_species, num_atoms, image)
    cp2k_to_siesta(filename_output_cp2k, filename_output_siesta, siesta_species_to_atomic)
    print('Finished.')



