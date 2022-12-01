from __future__ import division, print_function
import pandas as pd
import numpy as np
import glob
from scripts.formatting import load_coordinates
from scripts.formatting import print_xyz
import matplotlib.pyplot as plt

"""
    Slice .xyz along desired axis (useful for SMEAGOL transport calculations)
"""

# Files
folder = '/Volumes/Storage/Data/Work/Postdoc/Work/testing/cp2k-smeagol/bdt/cp2k/bulk/cell_opt/DIRECT_CELL_OPT'
input_filename = 'supercell_443.xyz'
output_filename = 'supercell_443-sortz.xyz'

# Read number of atoms and labels from .xyz file
cols = ['Species', 'X', 'Y', 'Z']
file_coord, num_atoms, species = load_coordinates.load_file_coord(folder, input_filename)
file_coord = file_coord.reset_index(drop=True)

file_coord = file_coord.sort_values('Z')
