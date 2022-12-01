from __future__ import division, print_function
import pandas as pd
import numpy as np
import glob
from scripts.formatting import load_coordinates
from scripts.formatting import print_xyz
import matplotlib.pyplot as plt
import math

"""
    create 1D chain .xyz 
"""

# folder = '/Volumes/Storage/Data/Work/Postdoc/Work/testing/cp2k-smeagol/examples/SMEAGOL_TEST/Au_chain/dft-8-pdos-50atoms'
folder = '/Volumes/Storage/Data/Work/Postdoc/Work/testing/smeagol/cx1/Smeagol_Tutorial_Files/DFT_NEGF_Transport/Day_1/Inputs/testing/au_chain'
output_filename = 'input.xyz'
num_atoms = 20
species = ['Au']
bond_length = 2.8

# Create dataframe
file_coord = pd.DataFrame(data={'X': [], 'Y': [], 'Z': []})

# Construct chain
for i in range(0, num_atoms):
    new_row = pd.Series({'X': 0,  'Y': 0, 'Z': 0 + i * bond_length})
    file_coord = pd.concat([file_coord, new_row.to_frame().T], ignore_index=True)

print('cell z', bond_length*num_atoms)

# Print to file
file_coord.insert(loc=0, column='A', value=species*num_atoms)
print_xyz.print_from_pandas(file_coord, num_atoms, '{}/{}'.format(folder, output_filename))
