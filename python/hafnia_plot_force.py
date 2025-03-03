from __future__ import division, print_function
import pandas as pd
import numpy as np
import glob
from general import load_coordinates
from general import load_forces
from general import parameters as param
from general import print_xyz
import matplotlib.pyplot as plt
import csv

"""
    Analysis script for au h2 as in clotilde phonon paper, to confirm SMEAGOL forces
"""

# folder_dft = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/supercell-1-1-5-bulk-6-cu-1.86/junction/bias/force/kpoints-2-2-V-0'
folder_dft = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/supercell-1-1-5-bulk-6-cu-1.86/junction/bias/force/kpoints-2-2-V-0-mpi-128-rs-dft-alpha-0.1'
# folder_dft = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/supercell-1-1-5-bulk-6-cu-1.86/junction/bias/force/kpoints-3-3-V-0'
folder_V0 = folder_dft
file_dft = 'dft_wfn-frc-1.xyz'
file_V0 = '0V-frc-1.xyz'

folders = [folder_dft, folder_V0]

# Read position
cols = ['Species', 'X', 'Y', 'Z']
coord_dft, _, _, _, species, num_atoms, _ = load_coordinates.load_values_coord(folder_dft, file_dft, cols)
coord_V0, _, _, _, _, _, _ = load_coordinates.load_values_coord(folder_V0, file_V0, cols)
index = np.arange(start=1, stop=num_atoms + 1)
dft_time = np.arange(start=1, stop=coord_dft.shape[0] + 1)
V0_time = np.arange(start=1, stop=coord_V0.shape[0] + 1)
diff = 0.2
step = 0
dimensions = ['x', 'y', 'z']
dimension = 0
dimension_string = dimensions[dimension]
save = False
xlim = [48 - diff, 122 + diff]
marker_size = 12


atom_cu_bulk = np.array([])
atom_cu_bulk_index = np.array([])
atom_cu_1 = np.array([])
atom_cu_1_index = np.array([])
atom_cu_2 = np.array([])
atom_cu_2_index = np.array([])
atom_hf = np.array([])
atom_hf_index = np.array([])
atom_o = np.array([])
atom_o_index = np.array([])

for i in range(len(species)):
    if species[i] == 'Cu_bulk':
        atom_cu_bulk = np.append(atom_cu_bulk, [coord_dft[step, dimension, i], coord_V0[step, dimension, i]])
        atom_cu_bulk_index = np.append(atom_cu_bulk_index, i)
    if species[i] == 'Cu_1':
        atom_cu_1 = np.append(atom_cu_1, [coord_dft[step, dimension, i], coord_V0[step, dimension, i]])
        atom_cu_1_index = np.append(atom_cu_1_index, i)
    if species[i] == 'Cu_2':
        atom_cu_2 = np.append(atom_cu_2, [coord_dft[step, dimension, i], coord_V0[step, dimension, i]])
        atom_cu_2_index = np.append(atom_cu_2_index, i)
    if species[i] == 'Hf':
        atom_hf = np.append(atom_hf, [coord_dft[step, dimension, i], coord_V0[step, dimension, i]])
        atom_hf_index = np.append(atom_hf_index, i)
    if species[i] == 'O':
        atom_o = np.append(atom_o, [coord_dft[step, dimension, i], coord_V0[step, dimension, i]])
        atom_o_index = np.append(atom_o_index, i)

# Plot x force for single timestep all atoms
print(atom_cu_2.shape)
fig_plot_1, ax_plot_1 = plt.subplots(figsize=(10, 4))
# ax_plot_1.plot(index, (coord_dft[step, dimension, :]) * param.hartree_per_bohr_to_ev_per_angstrom, 'k--', fillstyle='none')
# ax_plot_1.plot(1+atom_cu_2_index, atom_cu_2[::2] * param.hartree_per_bohr_to_ev_per_angstrom, 'bo-', fillstyle='none', label='Cu', markersize=marker_size)
ax_plot_1.plot(1+atom_cu_2_index[0:8], atom_cu_2[0:16:2] * param.hartree_per_bohr_to_ev_per_angstrom, 'bo--', fillstyle='none', label='Cu', markersize=marker_size)
ax_plot_1.plot(1+atom_cu_2_index[8:16], atom_cu_2[16:32:2] * param.hartree_per_bohr_to_ev_per_angstrom, 'bo--', fillstyle='none', markersize=marker_size)
ax_plot_1.plot(1+atom_hf_index, atom_hf[::2] * param.hartree_per_bohr_to_ev_per_angstrom, 'o--', color='orange', fillstyle='none', label='Hf', markersize=marker_size)
ax_plot_1.plot(1+atom_o_index, atom_o[::2] * param.hartree_per_bohr_to_ev_per_angstrom, 'ro--', fillstyle='none', label='O', markersize=marker_size)

# ax_plot_1.plot(index, (coord_V0[step, dimension, :]) * param.hartree_per_bohr_to_ev_per_angstrom, 'k-', fillstyle='none')
# ax_plot_1.plot(1+atom_cu_2_index, atom_cu_2[::2] * param.hartree_per_bohr_to_ev_per_angstrom, 'bo-', fillstyle='none', label='Cu', markersize=marker_size)
ax_plot_1.plot(1+atom_cu_2_index[0:8], atom_cu_2[1:17:2] * param.hartree_per_bohr_to_ev_per_angstrom, 'bx-', fillstyle='none', label='Cu V=0', markersize=marker_size)
ax_plot_1.plot(1+atom_cu_2_index[8:16], atom_cu_2[17:64:2] * param.hartree_per_bohr_to_ev_per_angstrom, 'bx-', fillstyle='none', markersize=marker_size)
ax_plot_1.plot(1+atom_hf_index, atom_hf[1::2] * param.hartree_per_bohr_to_ev_per_angstrom, 'x-', color='orange', fillstyle='none', label='Hf V=0', markersize=marker_size)
ax_plot_1.plot(1+atom_o_index, atom_o[1::2] * param.hartree_per_bohr_to_ev_per_angstrom, 'rx-', fillstyle='none', label='O V=0', markersize=marker_size)
ax_plot_1.set_xlabel('Atom index')
ax_plot_1.set_ylabel('Force {} (eV / Å)'.format(dimension_string))
ax_plot_1.set_xlim([xlim[0], xlim[1]])
ax_plot_1.legend(frameon=True)
fig_plot_1.tight_layout()
if save:
    for i in range(len(folders)):
        fig_plot_1.savefig('{}/force_{}.png'.format(folders[i], dimension_string), dpi=param.save_dpi)

# Analysis
cu_left_error = (atom_cu_2[0:16:2] - atom_cu_2[1:17:2]) * param.hartree_per_bohr_to_ev_per_angstrom
cu_right_error = (atom_cu_2[16:32:2] - atom_cu_2[17:64:2]) * param.hartree_per_bohr_to_ev_per_angstrom
o_error = (atom_o[::2] - atom_o[1::2]) * param.hartree_per_bohr_to_ev_per_angstrom
hf_error = (atom_hf[1::2] - atom_hf[::2]) * param.hartree_per_bohr_to_ev_per_angstrom

# Plot error
fig_plot_2, ax_plot_2 = plt.subplots(figsize=(10, 4))
ax_plot_2.plot(1+atom_cu_2_index[0:8], cu_left_error, 'bo--', fillstyle='none', label='Cu', markersize=marker_size)
ax_plot_2.plot(1+atom_cu_2_index[8:16], cu_right_error, 'bo--', fillstyle='none', markersize=marker_size)
ax_plot_2.plot(1+atom_hf_index, hf_error, 'o--', color='orange', fillstyle='none', label='Hf', markersize=marker_size)
ax_plot_2.plot(1+atom_o_index, o_error, 'ro--', fillstyle='none', label='O', markersize=marker_size)
ax_plot_2.set_xlabel('Atom index')
ax_plot_2.set_ylabel('Force {} (eV / Å)'.format(dimension_string))
ax_plot_2.set_xlim([xlim[0], xlim[1]])
ax_plot_2.legend(frameon=True)
fig_plot_2.tight_layout()

# Plot both
rows, cols = 2, 1
fig_plot_3, ax_plot_3 = plt.subplots(rows, cols, sharex='col', sharey='row', figsize=(10, 8))
ax_plot_3[0].plot(1+atom_cu_2_index[0:8], atom_cu_2[0:16:2] * param.hartree_per_bohr_to_ev_per_angstrom, 'bo--', fillstyle='none', label='Cu', markersize=marker_size)
ax_plot_3[0].plot(1+atom_cu_2_index[8:16], atom_cu_2[16:32:2] * param.hartree_per_bohr_to_ev_per_angstrom, 'bo--', fillstyle='none', markersize=marker_size)
ax_plot_3[0].plot(1+atom_hf_index, atom_hf[::2] * param.hartree_per_bohr_to_ev_per_angstrom, 'o--', color='orange', fillstyle='none', label='Hf', markersize=marker_size)
ax_plot_3[0].plot(1+atom_o_index, atom_o[::2] * param.hartree_per_bohr_to_ev_per_angstrom, 'ro--', fillstyle='none', label='O', markersize=marker_size)
ax_plot_3[0].plot(1+atom_cu_2_index[0:8], atom_cu_2[1:17:2] * param.hartree_per_bohr_to_ev_per_angstrom, 'bx-', fillstyle='none', label='Cu V=0', markersize=marker_size)
ax_plot_3[0].plot(1+atom_cu_2_index[8:16], atom_cu_2[17:64:2] * param.hartree_per_bohr_to_ev_per_angstrom, 'bx-', fillstyle='none', markersize=marker_size)
ax_plot_3[0].plot(1+atom_hf_index, atom_hf[1::2] * param.hartree_per_bohr_to_ev_per_angstrom, 'x-', color='orange', fillstyle='none', label='Hf V=0', markersize=marker_size)
ax_plot_3[0].plot(1+atom_o_index, atom_o[1::2] * param.hartree_per_bohr_to_ev_per_angstrom, 'rx-', fillstyle='none', label='O V=0', markersize=marker_size)
# ax_plot_3[0].set_xlabel('Atom index')
ax_plot_3[0].set_ylabel('Force {} (eV / Å)'.format(dimension_string))
ax_plot_3[0].set_xlim([xlim[0], xlim[1]])
ax_plot_3[0].legend(frameon=True)
ax_plot_3[1].plot(1+atom_cu_2_index[0:8], cu_left_error, 'bo-', fillstyle='none', label='Cu', markersize=marker_size)
ax_plot_3[1].plot(1+atom_cu_2_index[8:16], cu_right_error, 'bo-', fillstyle='none', markersize=marker_size)
ax_plot_3[1].plot(1+atom_hf_index, hf_error, 'o-', color='orange', fillstyle='none', label='Hf', markersize=marker_size)
ax_plot_3[1].plot(1+atom_o_index, o_error, 'ro-', fillstyle='none', label='O', markersize=marker_size)
ax_plot_3[1].set_xlabel('Atom index')
ax_plot_3[1].set_ylabel('Force error {} (eV / Å)'.format(dimension_string))
ax_plot_3[1].set_xlim([xlim[0], xlim[1]])
ax_plot_3[1].legend(frameon=True)
fig_plot_3.subplots_adjust(hspace=0)
fig_plot_3.tight_layout()
if save:
    for i in range(len(folders)):
        fig_plot_3.savefig('{}/force_error_{}.png'.format(folders[i], dimension_string), dpi=param.save_dpi)

if __name__ == "__main__":
    print('Finished.')
    plt.show()
