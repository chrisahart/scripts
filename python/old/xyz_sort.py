from __future__ import division, print_function
import pandas as pd
import numpy as np
import glob
from general import load_coordinates
from general import print_xyz
import matplotlib.pyplot as plt

"""
    Sort .xyz sequentially along desired axes (useful for transport calculations)
"""

# folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-capacitor/layers-1-2-3-4/cp2k/weightrho/testing/kpoints-4-4-20_hlb-auto_NEnergReal-64/input'
# folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-wire/jad/structures'
# input_filename = 'Au_Junc_wrapped.xyz'
# output_filename = 'Au_Junc_wrapped_sorted.xyz'

# folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-wire-ismael/jad/structures/jad_structures_replace-a'
# input_filename = 'structure_a.xyz'
# output_filename = 'structure_a_sorted.xyz'

# folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/melamine/cp2k-smeagol/structures/5A-tip/from-h2/geo_opt/ts2/from-cu-neb/9'
# input_filename = 'molecule.xyz'
# output_filename = 'molecule_translated.xyz'

folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/melamine/cp2k-smeagol/structures/5A-tip/from-h2/geo_opt-3A/geo_opt'
input_filename = 'gs.xyz'
output_filename = 'gs-3A.xyz'

filename_output = '{}/{}'.format(folder, output_filename)

# Read number of atoms and labels from .xyz file
cols = ['Species', 'X', 'Y', 'Z']
file_coord, num_atoms, species = load_coordinates.load_file_coord(folder, input_filename, cols)
file_coord = file_coord.reset_index(drop=True)
print(file_coord)
file_coord = file_coord.round(decimals=4)
print(file_coord)

print(file_coord.shape[0])
# shift = 8.67000
shift_start = 20
shift = np.array([0,0,-2])
# shift = -np.array([file_coord['X'][0], file_coord['Y'][0], file_coord['Z'][0]])
# shift = shift + np.array([7.19200,   5.93400 , 21.04300- 9.00300])
# shift = shift + np.array([4.43400 + (7.19200-4.43400)/2 ,   4.4330 ,28.2700-(30.36000-28.2700)])
# shift = np.array([-1.78700,   -1.28100,   0.00000 ])
for i in range(file_coord.shape[0]):
    if file_coord['Z'][i] > shift_start:
        file_coord['X'][i] = file_coord['X'][i] + shift[0]
        file_coord['Y'][i] = file_coord['Y'][i] + shift[1]
        file_coord['Z'][i] = file_coord['Z'][i] + shift[2]

# Delete species
# Au_lead_index = [i for i, e in enumerate(species) if e == 'Au_lead']
# Au_wire_index = [i for i, e in enumerate(species) if e == 'Au_wire']
# Au_al_index = [i for i, e in enumerate(species) if e == 'Au_al']
# Au_bl_index = [i for i, e in enumerate(species) if e == 'Au_bl']
# Au_cl_index = [i for i, e in enumerate(species) if e == 'Au_cl']
# Au_br_index = [i for i, e in enumerate(species) if e == 'Au_br']
# Au_cr_index = [i for i, e in enumerate(species) if e == 'Au_cr']
# H_index = [i for i, e in enumerate(species) if e == 'H']
# O_index = [i for i, e in enumerate(species) if e == 'O']
# delete_index = H_index + O_index
# delete_index = Au_wire_index + H_index + O_index + Au_al_index + Au_bl_index + Au_cl_index + Au_br_index + Au_cr_index
# delete_index = Au_wire_index + H_index + O_index
# delete_index = O_index

file_coord.insert(loc=0, column='A', value=pd.Series(species).values)
# file_coord = file_coord.drop(file_coord.index[delete_index])
# file_coord.sort_values(by=['Z'], inplace=True, ascending=[True])
file_coord.sort_values(by=['Z', 'X', 'Y'], inplace=True, ascending=[True, True, True])
# file_coord.sort_values(by=['Y', 'X', 'Z'], inplace=True, ascending=[True, True, True])
file_coord = file_coord.reset_index(drop=True)
num_atoms = file_coord.shape[0]
print(file_coord)
print(num_atoms)

print(file_coord[35:40])
print(np.max(file_coord['X']))
print(np.max(file_coord['Y']))
print(np.max(file_coord['Z']))

# Print to file
print_xyz.print_from_pandas(file_coord, num_atoms, filename_output, save_dp='%.4f')
