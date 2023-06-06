from __future__ import division, print_function
import pandas as pd
import numpy as np
import glob
from general import load_coordinates
from general import print_xyz
import matplotlib.pyplot as plt
import math
import ase
from ase.build import surface
from ase.visualize import view
from ase import Atoms
from ase.build import bulk
from ase.io import read, write
from ase.build import make_supercell
from ase.build import fcc111
from pathlib import Path

"""
    Create Au-BDT junction
"""

# Coordinates
a = 4.146  # CP2K optimised with TYPE DIRECT_CELL_OPT, KEEP_SYMMETRY TRUE
au_au_z_001 = a/2
au_au_z_111 = 2.406
Au_S = 2.65

# Files
left = 'left.xyz'
right = 'right.xyz'
bulk_ase = 'bulk_ase.xyz'
bulk = 'bulk.xyz'
molecule = 'molecule.xyz'
print_em = 'em.xyz'
dp_print = 3

# Make 001
# folder_ref = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/structures/001'
# folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/structures/001_Au-S-{}'.format(Au_S)
# folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/structures/001/atoms-116'
folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/structures/001/atoms-246'
folder_ref = folder
Path("{}".format(folder)).mkdir(parents=True, exist_ok=True)
au_au_z = au_au_z_001
x1 = a
y1 = a/2
y2 = a
Au_cell = Atoms('Au4',
              scaled_positions=[(0, 0, 0),
                                (0.5, 0.5, 0),
                                (0.5, 0, 0.5),
                                (0, 0.5, 0.5)],
              cell=[a, a, a],
              pbc=True)
supercell = make_supercell(Au_cell, np.diag([4, 4, 1]))
surface = surface(supercell, (0, 0, 1), layers=10)
surface_bulk = surface[[atom.index for atom in surface if atom.symbol == 'Au' and atom.position[2] < 8]]
surface_l = surface[[atom.index for atom in surface if atom.symbol == 'Au' and atom.position[2] < 13]]
surface_r = surface[[atom.index for atom in surface if atom.symbol == 'Au' and atom.position[2] < 11]]
write('{}/{}'.format(folder, left), surface_l)
write('{}/{}'.format(folder, right), surface_r)
write('{}/{}'.format(folder, bulk_ase), surface_bulk)
del_rows = [0, 1]

# Make 111
# folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/structures/111'
# au_au_z = au_au_z_111
# x1 = a
# y1 = a/2
# y2 = a
# x1 = 5.89400
# surface = fcc111('Au', a=a, size=(4, 4, 7))
# surface_bulk = surface[[atom.index for atom in surface if atom.symbol == 'Au' and atom.position[2] < 5]]
# surface_l = surface[[atom.index for atom in surface if atom.symbol == 'Au' and atom.position[2] < 13]]
# surface_r = surface[[atom.index for atom in surface if atom.symbol == 'Au' and atom.position[2] < 13]]
# del_rows = [0, 1]
# ase.io.write('{}/{}'.format(folder, left), surface_l)
# ase.io.write('{}/{}'.format(folder, right), surface_r)
# ase.io.write('{}/{}'.format(folder, bulk_ase), surface_bulk)

cols = ['Species', 'X', 'Y', 'Z']
file_coord_1, num_atoms_1, species_1 = load_coordinates.load_file_coord(folder, left, cols, del_rows=del_rows)
file_coord_1 = file_coord_1.reset_index(drop=True)

# molecule
file_coord_2, num_atoms_2, species_2 = load_coordinates.load_file_coord(folder_ref, molecule, cols)
file_coord_2 = file_coord_2.reset_index(drop=True)

# right
file_coord_3, num_atoms_3, species_3 = load_coordinates.load_file_coord(folder, right, cols, del_rows=del_rows)
file_coord_3 = file_coord_3.reset_index(drop=True)

# bulk
file_coord_4, num_atoms_4, species_4 = load_coordinates.load_file_coord(folder, bulk_ase, cols, del_rows=del_rows)
file_coord_4 = file_coord_4.reset_index(drop=True)

# Print bulk to file
file_coord_4_print = file_coord_4.round(dp_print)
file_coord_4_print.insert(loc=0, column='A', value=pd.Series(species_4).values)
file_coord_4_print.sort_values(by=['Z', 'X', 'Y'], inplace=True, ascending=[True, True, True])
print_xyz.print_from_pandas(file_coord_4_print, num_atoms_4, '{}/{}'.format(folder, bulk))
bulk_size = np.array([np.round(np.max(file_coord_4['X'])/au_au_z),
                      np.round(np.max(file_coord_4['Y'])/au_au_z), np.round(np.max(file_coord_4['Z'])/au_au_z)])

# Reset molecule
translate_reset = -1 * np.array([file_coord_2['X'][0], file_coord_2['Y'][0], file_coord_2['Z'][0]])
for i in range(0, num_atoms_2):
    file_coord_2['X'][i] = file_coord_2['X'][i] + translate_reset[0]
    file_coord_2['Y'][i] = file_coord_2['Y'][i] + translate_reset[1]
    file_coord_2['Z'][i] = file_coord_2['Z'][i] + translate_reset[2]

# Translate molecule
z2 = np.max(file_coord_1['Z'])
z1 = np.sqrt((Au_S ** 2) - (y2 - y1) ** 2) + z2
translate_molecule = np.array([x1, y1, z1])
for i in range(0, num_atoms_2):
    file_coord_2['X'][i] = file_coord_2['X'][i] + translate_molecule[0]
    file_coord_2['Y'][i] = file_coord_2['Y'][i] + translate_molecule[1]
    file_coord_2['Z'][i] = file_coord_2['Z'][i] + translate_molecule[2]

# Translate right
z2 = np.max(file_coord_2['Z'])
z1 = np.sqrt((Au_S ** 2) - (y2 - y1) ** 2) + z2
translate_right = np.array([0, 0, z1])
for i in range(0, num_atoms_3):
    file_coord_3['X'][i] = file_coord_3['X'][i] + translate_right[0]
    file_coord_3['Y'][i] = file_coord_3['Y'][i] + translate_right[1]
    file_coord_3['Z'][i] = file_coord_3['Z'][i] + translate_right[2]

# Join coordinates
left_molecule = pd.concat([file_coord_1, file_coord_2], ignore_index=True, sort=False)
left_molecule_right = pd.concat([left_molecule, file_coord_3], ignore_index=True, sort=False)
left_molecule_right = left_molecule_right.reset_index(drop=True)

# Join species
species_12 = pd.concat([species_1, species_2], ignore_index=True, sort=False)
species_123 = pd.concat([species_12, species_3], ignore_index=True, sort=False)
species_123 = species_123.reset_index(drop=True)

# Join coordinates and species and sort by z coordinate
num_atoms = left_molecule_right.shape[0]
left_molecule_right = left_molecule_right.round(dp_print)
em_size = np.array([np.round(np.max(left_molecule_right['X'])/au_au_z),
                    np.round(np.max(left_molecule_right['Y'])/au_au_z),
                    np.max(left_molecule_right['Z'])])
left_molecule_right.insert(loc=0, column='A', value=pd.Series(species_123).values)
left_molecule_right.sort_values(by=['Z', 'X', 'Y'], inplace=True, ascending=[True, True, True])
left_molecule_right = left_molecule_right.reset_index(drop=True)

# Print to file CP2K
print_xyz.print_from_pandas(left_molecule_right, num_atoms, '{}/{}'.format(folder, print_em))

# Print cell information
print('bulk right num_atoms', num_atoms_3)
print('bulk cell size', '{0:.3f}'.format(au_au_z*(bulk_size[0]+1)),
      '{0:.3f}'.format(au_au_z*(bulk_size[1]+1)), '{0:.3f}'.format(au_au_z*(bulk_size[2]+1)))
print('em num_atoms', num_atoms)
print('em cell size', '{0:.3f}'.format(au_au_z*(em_size[0]+1)),
      '{0:.3f}'.format(au_au_z*(em_size[1]+1)), '{0:.3f}'.format(au_au_z+em_size[2]))

