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
    Analysis script for au h2 as in clotilde phonon paper, to confirm SMEAGOL forces
"""

folder_dft = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/forces/AuH2/archer/geo_opt/kpoints-3-3-20_dft_MAX_FORCE-1.0E-03_CUTOFF-600_scf-1e-8'
file_dft = 'dft_wfn-2-frc-1.xyz'
folder_V0 = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/forces/AuH2/archer/geo_opt/kpoints-3-3-20_hlb-auto_chris_MAX_FORCE-1.0E-03_CUTOFF-600_scf-1e-8'
file_V0 = '0V-2-frc-1.xyz'
folders = [folder_dft, folder_V0]

# Read position
cols = ['Species', 'X', 'Y', 'Z']
coord_dft, _, _, _, _, num_atoms, _ = load_coordinates.load_values_coord(folder_dft, file_dft, cols)
coord_V0, _, _, _, _, _, _ = load_coordinates.load_values_coord(folder_V0, file_V0, cols)
index = np.arange(start=1, stop=num_atoms+1)
dft_time = np.arange(start=1, stop=coord_dft.shape[0]+1)
V0_time = np.arange(start=1, stop=coord_V0.shape[0]+1)
diff = 0.2

# Plot z force for single timestep all atoms
xlim = [58-diff, 65+diff]
dimension = 0
fig_plot_1, ax_plot_1 = plt.subplots()
ax_plot_1.plot(index, (coord_dft[dimension, 0, :]) * param.hartree_per_bohr_to_ev_per_angstrom, 'ko-', fillstyle='none', label='CP2K')
ax_plot_1.plot(index, (coord_V0[dimension, 0, :]) * param.hartree_per_bohr_to_ev_per_angstrom, 'ro-', fillstyle='none', label='CP2K-SMEAGOL')
ax_plot_1.set_xlabel('Atom index')
ax_plot_1.set_ylabel('Force x (eV / Å)')
ax_plot_1.set_xlim([xlim[0], xlim[1]])
ax_plot_1.legend(frameon=False)
fig_plot_1.tight_layout()
for i in range(len(folders)):
    fig_plot_1.savefig('{}/force_x.png'.format(folders[i]), dpi=param.save_dpi)

# Plot z force for single atom all timesteps
atom = np.array([61]) - 1
atoms = np.array([61, 62, 63]) - 1
xlim = [1-diff, 18+diff]
fig_plot_2, ax_plot_2 = plt.subplots()
ax_plot_2.plot(dft_time, (coord_dft[:, 0, atom]) * param.hartree_per_bohr_to_ev_per_angstrom, 'ko-', fillstyle='none', label='CP2K')
ax_plot_2.plot(V0_time, (coord_V0[:, 0, atom]) * param.hartree_per_bohr_to_ev_per_angstrom, 'ro-', fillstyle='none', label='CP2K-SMEAGOL')
ax_plot_2.plot(dft_time, (coord_dft[:, 0, atom+1]) * param.hartree_per_bohr_to_ev_per_angstrom, 'ko--', fillstyle='none')
ax_plot_2.plot(V0_time, (coord_V0[:, 0, atom+1]) * param.hartree_per_bohr_to_ev_per_angstrom, 'ro--', fillstyle='none')
ax_plot_2.set_xlabel('GEO_OPT step')
ax_plot_2.set_ylabel('Force x (eV / Å)')
ax_plot_2.set_xlim([xlim[0], xlim[1]])
ax_plot_2.legend(frameon=False)
fig_plot_2.tight_layout()
for i in range(len(folders)):
    fig_plot_2.savefig('{}/force_x_time.png'.format(folders[i]), dpi=param.save_dpi)

# Plot all forces on atom
# rows, cols = 3, 1
# fig_plot_forces_all, ax_plot_forces_all = plt.subplots(rows, cols, sharex='col', sharey='row', figsize=(6, 8))
# ax_plot_forces_all[0].plot(dft_time, (coord_dft[:, 0, atom]) * param.hartree_per_bohr_to_ev_per_angstrom, 'ko-', fillstyle='none', label='CP2K')
# ax_plot_forces_all[0].plot(V0_time, (coord_V0[:, 0, atom]) * param.hartree_per_bohr_to_ev_per_angstrom, 'ro-', fillstyle='none', label='CP2K-SMEAGOL')
# ax_plot_forces_all[1].plot(dft_time, (coord_dft[:, 1, atom]) * param.hartree_per_bohr_to_ev_per_angstrom, 'ko-', fillstyle='none', label='CP2K')
# ax_plot_forces_all[1].plot(V0_time, (coord_V0[:, 1, atom]) * param.hartree_per_bohr_to_ev_per_angstrom, 'ro-', fillstyle='none', label='CP2K-SMEAGOL')
# ax_plot_forces_all[2].plot(dft_time, (coord_dft[:, 2, atom]) * param.hartree_per_bohr_to_ev_per_angstrom, 'ko-', fillstyle='none', label='CP2K')
# ax_plot_forces_all[2].plot(V0_time, (coord_V0[:, 2, atom]) * param.hartree_per_bohr_to_ev_per_angstrom, 'ro-', fillstyle='none', label='CP2K-SMEAGOL')
# ax_plot_forces_all[0].set_ylabel('Force x (eV / Å)')
# ax_plot_forces_all[1].set_ylabel('Force y (eV / Å)')
# ax_plot_forces_all[2].set_ylabel('Force z (eV / Å)')
# ax_plot_forces_all[2].set_xlabel('GEO_OPT step')
# ax_plot_forces_all[0].set_xlim([xlim[0], xlim[1]])
# ax_plot_forces_all[1].set_xlim([xlim[0], xlim[1]])
# ax_plot_forces_all[2].set_xlim([xlim[0], xlim[1]])
# ax_plot_forces_all[0].legend(frameon=False)
# fig_plot_forces_all.tight_layout()
# for i in range(len(folders)):
#     fig_plot_forces_all.savefig('{}/force_all_time.png'.format(folders[i]), dpi=param.save_dpi)

# Plot all x force on atoms
# rows, cols = 3, 1
# fig_plot_atoms_all, ax_plot_atoms_all = plt.subplots(rows, cols, sharex='col', sharey='row', figsize=(6, 8))
# ax_plot_atoms_all[0].plot(dft_time, (coord_dft[:, 0, atom-1]) * param.hartree_per_bohr_to_ev_per_angstrom, 'ko-', fillstyle='none', label='CP2K')
# ax_plot_atoms_all[0].plot(V0_time, (coord_V0[:, 0, atom-1]) * param.hartree_per_bohr_to_ev_per_angstrom, 'ro-', fillstyle='none', label='CP2K-SMEAGOL')
# ax_plot_atoms_all[1].plot(dft_time, (coord_dft[:, 0, atom]) * param.hartree_per_bohr_to_ev_per_angstrom, 'ko-', fillstyle='none', label='CP2K')
# ax_plot_atoms_all[1].plot(V0_time, (coord_V0[:, 0, atom]) * param.hartree_per_bohr_to_ev_per_angstrom, 'ro-', fillstyle='none', label='CP2K-SMEAGOL')
# ax_plot_atoms_all[2].plot(dft_time, (coord_dft[:, 0, atom-2]) * param.hartree_per_bohr_to_ev_per_angstrom, 'ko-', fillstyle='none', label='CP2K')
# ax_plot_atoms_all[2].plot(V0_time, (coord_V0[:, 0, atom-2]) * param.hartree_per_bohr_to_ev_per_angstrom, 'ro-', fillstyle='none', label='CP2K-SMEAGOL')
# ax_plot_atoms_all[0].set_ylabel('Force x atom 4 (eV / Å)')
# ax_plot_atoms_all[1].set_ylabel('Force x atom 5 (eV / Å)')
# ax_plot_atoms_all[2].set_ylabel('Force x atom 6 (eV / Å)')
# ax_plot_atoms_all[2].set_xlabel('GEO_OPT step')
# ax_plot_atoms_all[0].set_xlim([xlim[0], xlim[1]])
# ax_plot_atoms_all[1].set_xlim([xlim[0], xlim[1]])
# ax_plot_atoms_all[2].set_xlim([xlim[0], xlim[1]])
# ax_plot_atoms_all[0].legend(frameon=False)
# fig_plot_atoms_all.tight_layout()
# for i in range(len(folders)):
#     fig_plot_atoms_all.savefig('{}/force_all_time.png'.format(folders[i]), dpi=param.save_dpi)
    
if __name__ == "__main__":
    print('Finished.')
    plt.show()
