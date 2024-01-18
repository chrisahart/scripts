import numpy as np
from matplotlib import pyplot as plt
from ase.io.cube import read_cube_data
import optparse
import os
from ase.build import fcc111
from ase.build import surface
from ase.visualize import view
from ase import Atoms
from ase.build import bulk
from ase.io import read, write
from ase.build import make_supercell
from ase.build import fcc111, add_adsorbate
from ase.constraints import FixAtoms

"""

"""

# folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-wire-ismael/jad/structures/jad_structures_replace-a'
folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/au-capacitor/cp2k-smeagol/au-111'
filename_bulk = 'structure_bulk.xyz'
filename_em = 'structure_em.xyz'
file_jad = 'Au_Junc_wrapped_sorted.xyz'

a = 4.1983  # Jad

# Create the bulk Au(111) structure
slab = fcc111('Au', size=(12, 12, 50), a=a)
z_max = slab.positions[-1][2]
view(slab)
# Define the desired cell size
cell_size = np.array([1.7811878396732908E+01, 1.5425539180689931E+01, 100])
slab.center(about=(0, 0, z_max/2))
# view(slab)

# Remove atoms from outside cell
atoms_within_cell = Atoms(cell=cell_size)
for atom in slab:
    if all(0 <= atom.position[i] < cell_size[i] for i in range(3)):
        atoms_within_cell.append(atom)
atoms_within_cell.set_cell(cell_size)
# view(atoms_within_cell)

# Cut cell
z_start = 2
z_end = 64
jad_structure_cut = Atoms(cell=cell_size)
for atom in atoms_within_cell:
    if atom.position[2] > z_start and atom.position[2] < z_end:
        jad_structure_cut.append(atom)
jad_structure_cut.set_cell(cell_size)
# view(jad_structure_cut)

# atoms_within_cell = atoms_within_cell[[atom.index for atom in surface if atom.position[2] < 3]]
# view(atoms_within_cell)
#
# Sort the atoms by z, then x, then y coordinates
sorted_indices = np.lexsort((jad_structure_cut.positions[:, 1], jad_structure_cut.positions[:, 0], jad_structure_cut.positions[:, 2]))
atoms_within_cell = jad_structure_cut[sorted_indices]

# Shift Au atoms
atom_0 = jad_structure_cut.positions[0]
jad_atom_0 = [0, 0, 0]
# jad_atom_0 = [0, 3.42790, 0]
jad_structure_cut.positions[:, 0] -= atom_0[0] - jad_atom_0[0]
jad_structure_cut.positions[:, 1] -= atom_0[1] - jad_atom_0[1]
jad_structure_cut.positions[:, 2] -= atom_0[2] - jad_atom_0[2]
view(jad_structure_cut)

# # Bulk left a only
# z_cut = 2
# bulk_left_a = Atoms(cell=cell_size)
# for atom in atoms_within_cell:
#     if atom.position[2] < z_cut:
#         bulk_left_a.append(atom)
# bulk_left_a.positions[:, 2] =+ 2.42389*3
#
# # Bulk right
# bulk_right = atoms_within_cell.copy()
# bulk_right.positions[:, 2] += 3.0226752151153697E+01 + 2.42389*3
# write('{}/{}'.format(folder, filename_bulk), atoms_within_cell, format='xyz')

# # Load Jad structure
# jad_structure = read('{}/{}'.format(folder, file_jad), format='xyz')
#
# # Remove Jad structure z=0
# z_cut = 2
# jad_structure_cut = Atoms(cell=cell_size)
# for atom in jad_structure:
#     if atom.position[2] > z_cut:
#         jad_structure_cut.append(atom)
# jad_structure_cut.set_cell(cell_size)
#
# # Translate Jad structure z=7.27167
# jad_structure_cut.positions[:, 2] += 2.42389*3
#
# # Combine Jad structure with leads
# all = atoms_within_cell + bulk_left_a + jad_structure_cut + bulk_right
write('{}/{}'.format(folder, filename_em), jad_structure_cut, format='xyz')
#
# view(atoms_within_cell)
