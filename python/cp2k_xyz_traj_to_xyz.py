from __future__ import division, print_function
import pandas as pd
import numpy as np
import glob
from general import load_coordinates
from general import print_xyz
import matplotlib.pyplot as plt
import csv
from ase.io import read, write
from ase import Atoms

"""
    Use python ASE to import .xyz trajectory and print each coordinates of each timestep to a file
"""

folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/supercell-1-1-2-bulk-6-cu-1.86/junction/geo_opt/free-hf-o/kpoints-2-2-USE_GUESS-V-0.4'
input_filename = 'output/trajectory-clean/0.4V-pos-1.xyz'

output_folder = 'output/geometries'
trajectory = read('{}/{}'.format(folder, input_filename), index=':')
print(f"Number of frames in the trajectory: {len(trajectory)}")

# Loop through each frame and save it as an individual .xyz file
for i, atoms in enumerate(trajectory):
    output_file = f'{folder}/{output_folder}/geom-{i}.xyz'
    print(f"Writing frame {i + 1} to {output_file}")

    # Write the current frame to the .xyz file
    atoms.write(output_file, format='xyz')
