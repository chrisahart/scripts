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
    # species_val = species.copy()
    # for i in range(0, num_atoms):
    #     if species[i] == "Cu_bulk":
    #         species_val[i] = 1
    #     elif species[i] == "Cu":
    #         species_val[i] = 1
    #     elif species[i] == "N":
    #         species_val[i] = 2
    #     elif species[i] == "C":
    #         species_val[i] = 3
    #     elif species[i] == "H":
    #         species_val[i] = 4
        # elif species[i] == "Au":
        #     species_val[i] = 5
        # elif species[i] == "Au_surf":
        #     species_val[i] = 6
        # elif species[i] == "Au_tip":
        #     species_val[i] = 7
        # elif species[i] == "N":
        #     species_val[i] = 8
        # else:
        #     print('undefined element', species[i])

    # Replace species with value
    species_val = species.copy()
    for i in range(0, num_atoms):
        if species[i] == "Cu_bulk":
            species_val[i] = 1
        elif species[i] == "Cu_1":
            species_val[i] = 1
        elif species[i] == "Cu_2":
            species_val[i] = 1
        elif species[i] == "Hf":
            species_val[i] = 2
        elif species[i] == "O":
            species_val[i] = 3
        else:
            print('undefined element', species[i])

    # Print to file
    file_coord.insert(loc=3, column='A', value=pd.Series(species_val).values)
    file_coord.insert(loc=4, column='B', value=pd.Series(species).values)
    file_coord.to_csv(filename_output, index=False, header=False, quoting=csv.QUOTE_NONE, sep=" ")


# folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/forces/au-chain/paper_ruoxing/siesta-smeagol/kpoints-2-2-20/input'
folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/structures/001/atoms-948'
# input_filename = 'em.xyz'
# output_filename = 'em.siesta'
input_filename = 'bulk.xyz'
output_filename = 'bulk.siesta'

folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/melamine/cp2k-smeagol/structures/5A-tip/from-h2/geo_opt/ts2/from-cu'
input_filename = 'ts2.xyz'
output_filename = 'ts2.siesta'

folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/melamine/cp2k-smeagol/structures/5A-tip/from-h2/geo_opt-3A/geo_opt'
input_filename = 'gs-3A.xyz'
output_filename = 'gs-3A.siesta'

folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/structures/pymatgen/cu/supercell-1-1-3-cu-1.86'
input_filename = 't-hfo2_junction.xyz'
output_filename = 't-hfo2_capacitor.siesta'

filename_output = '{}/{}'.format(folder, output_filename)


if __name__ == "__main__":
    print('Finished.')
    cp2k_to_siesta(folder, input_filename, filename_output)
