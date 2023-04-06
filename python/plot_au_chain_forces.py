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
import xyz_siesta_to_cp2k

"""
    Analysis script for au chain as in Rouxing paper
"""

folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/forces/au-chain/paper_ruoxing/kpoints-2-2-20/'
file_dft_1 = 'au-chain_szv_dftenergy-frc-1.xyz'
file_V0_1 = 'au-chain_szv_zerobias-frc-1.xyz'
file_dft_2 = 'au-chain_szv_shiftx_dftenergy-frc-1.xyz'
file_V0_2 = 'au-chain_szv_shiftx_zerobias-frc-1.xyz'
dimension = 'X'

# Read position
cols = ['Species', 'X', 'Y', 'Z']
dft_1, num_atoms, _ = load_coordinates.load_file_coord(folder, file_dft_1, cols)
V0_1, num_atoms, _ = load_coordinates.load_file_coord(folder, file_V0_1, cols)
dft_2, num_atoms, _ = load_coordinates.load_file_coord(folder, file_dft_2, cols)
V0_2, num_atoms, _ = load_coordinates.load_file_coord(folder, file_V0_2, cols)
index = np.arange(start=1, stop=num_atoms+1)

# Plot z force
xlim = [0, 2]
fig_plot_1, ax_plot_1 = plt.subplots()
ax_plot_1.plot(index, dft_2['X'] - dft_1['X'], 'ko-', fillstyle='none', label='CP2K')
ax_plot_1.plot(index, V0_2['X'] - V0_1['X'], 'ro-', fillstyle='none', label='CP2K-SMEAGOL')
ax_plot_1.set_xlabel('Atom index')
ax_plot_1.set_ylabel('Force / au')
ax_plot_1.legend(frameon=False)
fig_plot_1.tight_layout()
fig_plot_1.savefig('{}/force_x.png'.format(folder), dpi=param.save_dpi)

# Plot all force
# rows, cols = 3, 1
# fig_plot_all, ax_plot_all = plt.subplots(rows, cols, sharex='col', sharey='row', figsize=(6, 8))
# ax_plot_all[0].plot(index, (dft_2['X']) * param.hartree_per_bohr_to_ev_per_angstrom, 'ko-', fillstyle='none', label='CP2K')
# ax_plot_all[0].plot(index, (V0_2['X']) * param.hartree_per_bohr_to_ev_per_angstrom, 'ro-', fillstyle='none', label='CP2K-SMEAGOL')
# ax_plot_all[1].plot(index, (dft_2['Y']) * param.hartree_per_bohr_to_ev_per_angstrom, 'ko-',  fillstyle='none', label='CP2K')
# ax_plot_all[1].plot(index, (V0_2['Y']) * param.hartree_per_bohr_to_ev_per_angstrom, 'ro-', fillstyle='none',   label='CP2K-SMEAGOL')
# ax_plot_all[2].plot(index, (dft_2['Z']) * param.hartree_per_bohr_to_ev_per_angstrom, 'ko-',  fillstyle='none', label='CP2K')
# ax_plot_all[2].plot(index, (V0_2['Z']) * param.hartree_per_bohr_to_ev_per_angstrom, 'ro-', fillstyle='none',  label='CP2K-SMEAGOL')
# ax_plot_all[0].set_ylabel('Force x (eV / Å)')
# ax_plot_all[1].set_ylabel('Force y (eV / Å)')
# ax_plot_all[2].set_ylabel('Force z (eV / Å)')
# ax_plot_all[2].set_xlabel('Atom index')
# ax_plot_all[0].legend(frameon=False)
# fig_plot_all.tight_layout()
# fig_plot_all.savefig('{}/force_all.png'.format(folder), dpi=param.save_dpi)

# Plot all force differences
rows, cols = 3, 1
fig_plot_all, ax_plot_all = plt.subplots(rows, cols, sharex='col', sharey='row', figsize=(6, 8))
ax_plot_all[0].plot(index, (dft_2['X'] - dft_1['X'])*param.hartree_per_bohr_to_ev_per_angstrom, 'ko-', fillstyle='none', label='CP2K')
ax_plot_all[0].plot(index, (V0_2['X'] - V0_1['X'])*param.hartree_per_bohr_to_ev_per_angstrom, 'ro-', fillstyle='none', label='CP2K-SMEAGOL')
ax_plot_all[1].plot(index, (dft_2['Y'] - dft_1['Y'])*param.hartree_per_bohr_to_ev_per_angstrom, 'ko-', fillstyle='none', label='CP2K')
ax_plot_all[1].plot(index, (V0_2['Y'] - V0_1['Y'])*param.hartree_per_bohr_to_ev_per_angstrom, 'ro-', fillstyle='none', label='CP2K-SMEAGOL')
ax_plot_all[2].plot(index, (dft_2['Z'] - dft_1['Z'])*param.hartree_per_bohr_to_ev_per_angstrom, 'ko-', fillstyle='none', label='CP2K')
ax_plot_all[2].plot(index, (V0_2['Z'] - V0_1['Z'])*param.hartree_per_bohr_to_ev_per_angstrom, 'ro-', fillstyle='none', label='CP2K-SMEAGOL')
ax_plot_all[0].set_ylabel('Force x (eV / Å)')
ax_plot_all[1].set_ylabel('Force y (eV / Å)')
ax_plot_all[2].set_ylabel('Force z (eV / Å)')
ax_plot_all[2].set_xlabel('Atom index')
ax_plot_all[0].legend(frameon=False)
fig_plot_all.tight_layout()
fig_plot_all.savefig('{}/force_all.png'.format(folder), dpi=param.save_dpi)

if __name__ == "__main__":
    print('Finished.')
    plt.show()
