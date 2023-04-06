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


def siesta_to_cp2k(filename_input, filename_output):

    # Calculate number of atoms
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


# folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/AuH2/structures'
folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/AuH2/transport/cp2k-smeagol/young/transmission/input'

input_filename = 'em.siesta'
output_filename = 'em.xyz'

# input_filename = 'bulk.siesta'
# output_filename = 'bulk.xyz'

filename_input = '{}/{}'.format(folder, input_filename)
filename_output = '{}/{}'.format(folder, output_filename)

if __name__ == "__main__":
    print('Finished.')
    siesta_to_cp2k(filename_input, filename_output)



