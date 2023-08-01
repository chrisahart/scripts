from __future__ import division, print_function
import pandas as pd
import numpy as np
from general import load_coordinates
from general import print_xyz
import matplotlib.pyplot as plt

"""
    .xyz translate
    Reads .xyz file and translates atomic positions
"""

folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/theory_support/hannah/structures/nanoparticle_carbide/hexagonal/units-1'
# folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/theory_support/hannah/structures/nanoparticle_carbide/cubic/atoms-27'
input_filename = 'carbide_opt.xyz'
# output_filename = 'carbide_Ptlayer_translated_hollow_z2.xyz'
output_filename = 'carbide_opt_ring-centered2.xyz'

# Read number of atoms and labels from .xyz file
cols = ['Species', 'X', 'Y', 'Z']
file_coord, num_atoms, species = load_coordinates.load_file_coord(folder, input_filename, cols)
# translate = np.array([7.23824, -1.39300, 13.64856 + 2])
# translate = np.array([4.82549, -2.78600, 13.64856 + 2])
# translate = np.array([4.82549, -0.00000 , 13.64856 + 2])
# translate = np.array([2.41275,  -1.39300,  13.64856 + 4])
translate = np.array([15.89836,  16.56547,  23.57468 ])


# Set atom to position (0, 0, 0)
# atom = 13 - 1
# atom = 19 - 1
# atom = 1 - 1
atom = 9 - 1
# atom = 26 - 1
pos = [file_coord['X'][atom], file_coord['Y'][atom], file_coord['Z'][atom]]
print(atom, pos)
for i in range(0, num_atoms):
    file_coord['X'][i] = file_coord['X'][i] - pos[0]
    file_coord['Y'][i] = file_coord['Y'][i] - pos[1]
    file_coord['Z'][i] = file_coord['Z'][i] - pos[2]

# Calculate midpoint of 3 atoms for hollow site
# a = np.array((4.82549,  -2.78600,  13.64856)) # Top
# b = np.array((4.82549,  -0.00000,  13.64856)) # Top
a = np.array((4.82549,  -2.78600,  0))  # Bottom
b = np.array((4.82549,  -0.00000,  0))  # Bottom
midpoint = np.array(((a[0]+b[0])/2, (a[1]+b[1])/2, (a[2]+b[2])/2))
translate = midpoint
# translate[0] = translate[0] - 0.5
translate[2] = translate[2] - 2

# Calculate midpoint of 3 atoms for hollow site
# a = np.array((4.82549,  -2.78600, 13.64856))
# b = np.array((4.82549,  -0.00000, 13.64856))
# c = np.array((7.23824 , -1.39300, 13.64856))
# a = np.array((4.82549,  -2.78600, 13.64856))
# b = np.array((7.23824, -4.17900, 13.64856))
# c = np.array((7.23824, -1.39300, 13.64856))
# a = np.array((2.41275 , -1.39300,  13.64856))
# b = np.array((2.41275  , 1.39300  ,13.64856 ))
# c = np.array(( 4.82549 ,  0.00000,  13.64856))
# midpoint = np.array(((a[0]+b[0]+c[0])/3, (a[1]+b[1]+c[1])/3, (a[2]+b[2]+c[2])/3))
# translate = midpoint
# translate[2] = translate[2] + 4

# Translate atoms
for i in range(0, num_atoms):
    file_coord['X'][i] = file_coord['X'][i] + translate[0]
    file_coord['Y'][i] = file_coord['Y'][i] + translate[1]
    file_coord['Z'][i] = file_coord['Z'][i] + translate[2]

# Change atom species
print(species)
for i in range(0, len(species)):
    print(species[i])
    if species[i] == 'Pt': species[i] ='Pt_PtC'
    if species[i] == 'C': species[i] = 'C_PtC'

# Print to file
file_coord.insert(loc=0, column='A', value=pd.Series(species).values)
print_xyz.print_from_pandas(file_coord, num_atoms, '{}/{}'.format(folder, output_filename))

