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


def cp2k_to_siesta(folder, input_filename, filename_output):

    cols = ['Species', 'X', 'Y', 'Z', 'Label']
    file_coord, num_atoms, species = load_coordinates.load_file_coord(folder, input_filename, cols)

    # Replace species with value
    species_val = species.copy()
    for i in range(0, num_atoms):
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


folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/structures/001_opt-cg'
input_filename = 'em.xyz'
output_filename = 'em.siesta'
filename_output = '{}/{}'.format(folder, output_filename)


if __name__ == "__main__":
    print('Finished.')
    cp2k_to_siesta(folder, input_filename, filename_output)
