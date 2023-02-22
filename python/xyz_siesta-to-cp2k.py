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

folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/Pt-H2-Pt/structures'
input_filename = 'em.siesta'
output_filename = 'em.xyz'

filename_input = '{}/{}'.format(folder, input_filename)
filename_output = '{}/{}'.format(folder, output_filename)

# If first row has number of atoms
# cols = ['X', 'Y', 'Z', 'Num', 'Species']
# file_coord = pd.read_csv(filename_input, names=cols, delim_whitespace=True)
# num_atoms = int(float(file_coord['X'][0]))
# species = file_coord['Species'][1:num_atoms + 1]
# file_coord = file_coord.apply(pd.to_numeric, errors='coerce')
# file_coord = file_coord.dropna(axis='rows', thresh=2)
# file_coord = file_coord.dropna(axis='columns', thresh=1)
# file_coord = file_coord.drop('Num', axis=1)
# print(file_coord)

# # Calculate number of atoms
# cols = ['X', 'Y', 'Z', 'Num', 'Species']
# file_coord = pd.read_csv(filename_input, names=cols, delim_whitespace=True)
# num_atoms = file_coord.shape[0]
# print(file_coord)
# species = file_coord['Species'][num_atoms]
# file_coord = file_coord.apply(pd.to_numeric, errors='coerce')
# file_coord = file_coord.dropna(axis='rows', thresh=2)
# file_coord = file_coord.dropna(axis='columns', thresh=1)
# file_coord = file_coord.drop('Num', axis=1)
# print(file_coord)

# Calculate number of atoms
cols = ['X', 'Y', 'Z', 'Species_num', 'Species', 'Num']
file_coord = pd.read_csv(filename_input, names=cols, delim_whitespace=True)
num_atoms = file_coord.shape[0]
print(file_coord)
species = file_coord['Species']
file_coord = file_coord.apply(pd.to_numeric, errors='coerce')
file_coord = file_coord.dropna(axis='rows', thresh=2)
file_coord = file_coord.dropna(axis='columns', thresh=1)
file_coord = file_coord.drop('Num', axis=1)
file_coord = file_coord.drop('Species_num', axis=1)
print(file_coord)

# Print to file
file_coord.insert(loc=0, column='A', value=pd.Series(species).values)
# file_coord.to_csv(filename_output, index=False, header=False, quoting=csv.QUOTE_NONE, sep=" ")
print_xyz.print_from_pandas(file_coord, num_atoms, filename_output)

