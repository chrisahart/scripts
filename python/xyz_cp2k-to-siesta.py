from __future__ import division, print_function
import pandas as pd
import numpy as np
import glob
from scripts.formatting import load_coordinates
from scripts.formatting import print_xyz
import matplotlib.pyplot as plt
import csv

"""
    Convert XYZ to SIESTA %block AtomicCoordinatesAndAtomicSpecies
"""

# folder = '/Volumes/Storage/Data/Work/Postdoc/Work/testing/cp2k-smeagol/bdt/structures/hollow-site/7-6'
# input_filename = 'em-slice-leads-lr.xyz'
# output_filename = 'em-slice-leads-lr.siesta'

# folder = '/Volumes/Storage/Data/Work/Postdoc/Work/calculations/transport/sergey'
# input_filename = 'sergey-bulk-right.xyz'
# output_filename = 'sergey-bulk-right.siesta'

folder = '/Volumes/Storage/Data/Work/Postdoc/Work/calculations/testing/cp2k-smeagol/bdt/structures/hollow-site/5-4_3x3-order'
input_filename = 'em.xyz'
output_filename = 'em.siesta'
# input_filename = 'supercell_3x3-4.xyz'
# output_filename = 'supercell_3x3-4.siesta'

filename_output = '{}/{}'.format(folder, output_filename)

# Read number of atoms and labels from .xyz file
# cols = ['Species', 'X', 'Y', 'Z']
# file_coord, num_atoms, species = load_coordinates.load_file_coord(folder, input_filename)

cols = ['Species', 'X', 'Y', 'Z', 'Label']
file_coord, num_atoms, species = load_coordinates.load_file_coord(folder, input_filename, cols)

# Replace species with value
species_val = species.copy()
for i in range(1, num_atoms+1):
    if species[i] == "Au":
        species_val[i] = 1
    elif species[i] == "S":
        species_val[i] = 2
    elif species[i] == "C":
        species_val[i] = 3
    elif species[i] == "H":
        species_val[i] = 4
    else:
        print('undefined element')

# Print to file
file_coord.insert(loc=3, column='A', value=pd.Series(species_val).values)
file_coord.insert(loc=4, column='B', value=pd.Series(species).values)
file_coord.to_csv(filename_output, index=False, header=False, quoting=csv.QUOTE_NONE, sep=" ")
