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
    Analysis script for au chain as in Rouxing paper
"""

folder_cp2k_1 = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/forces/au-chain/paper_ruoxing/cp2k-smeagol/energy_force/kpoints-2-2-20_constrained_Au-2.8/'
folder_cp2k_2 = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/forces/au-chain/paper_ruoxing/cp2k-smeagol/energy_force/kpoints-2-2-20_constrained_Au-2.8_omega-F/'
folder_siesta_1 = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/forces/au-chain/paper_ruoxing/siesta-smeagol/kpoints-2-2-20_constrained_positive_Au-2.8/'
folder_siesta_2 = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/forces/au-chain/paper_ruoxing/siesta-smeagol/kpoints-2-2-20_constrained_positive_Au-2.8_omega-F/'
# folder_cp2k_1 = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/forces/au-chain/paper_ruoxing/cp2k-smeagol/energy_force/kpoints-2-2-20_constrained/'
# folder_cp2k_2 = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/forces/au-chain/paper_ruoxing/cp2k-smeagol/energy_force/kpoints-2-2-20_constrained_omega-F/'
# folder_siesta_1 = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/forces/au-chain/paper_ruoxing/siesta-smeagol/kpoints-2-2-20_constrained_positive/'
# folder_siesta_2 = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/forces/au-chain/paper_ruoxing/siesta-smeagol/kpoints-2-2-20_constrained_positive_omega-F/'
labels_cp2k_1 = ['CP2K', 'CP2K-SMEAGOL V=0']
labels_cp2k_2 = ['CP2K', 'CP2K-SMEAGOL V=0 Ω=0']
labels_siesta_1 = ['SIESTA', 'SIESTA-SMEAGOL V=0']
labels_siesta_2 = ['SIESTA', 'SIESTA-SMEAGOL V=0 Ω=0']
plot_cp2k_1 = True
plot_cp2k_2 = True
plot_siesta_1 = False
plot_siesta_2 = False
use_ylim = [False, True, False]
ylim_xyz = np.array([(-2, 2), (-1, 1), (-1, 1)])

# Read CP2K
cols = ['Species', 'X', 'Y', 'Z']
cp2k_1_dft_1, num_atoms, _ = load_coordinates.load_file_coord(folder_cp2k_1, 'au-chain_szv_dftenergy-frc-1.xyz', cols)
cp2k_1_V0_1, _, _ = load_coordinates.load_file_coord(folder_cp2k_1, 'au-chain_szv_zerobias-frc-1.xyz', cols)
cp2k_1_dft_2, _, _ = load_coordinates.load_file_coord(folder_cp2k_1, 'au-chain_szv_shiftx_dftenergy-frc-1.xyz', cols)
cp2k_1_V0_2, _, _ = load_coordinates.load_file_coord(folder_cp2k_1, 'au-chain_szv_shiftx_zerobias-frc-1.xyz', cols)
cp2k_2_dft_1, _, _ = load_coordinates.load_file_coord(folder_cp2k_2, 'au-chain_szv_dftenergy-frc-1.xyz', cols)
cp2k_2_V0_1, _, _ = load_coordinates.load_file_coord(folder_cp2k_2, 'au-chain_szv_zerobias-frc-1.xyz', cols)
cp2k_2_dft_2, _, _ = load_coordinates.load_file_coord(folder_cp2k_2, 'au-chain_szv_shiftx_dftenergy-frc-1.xyz', cols)
cp2k_2_V0_2, _, _ = load_coordinates.load_file_coord(folder_cp2k_2, 'au-chain_szv_shiftx_zerobias-frc-1.xyz', cols)
index = np.arange(start=1, stop=num_atoms+1)

# Read SIESTA
cols_siesta = ['A', 'B', 'X', 'Y', 'Z']
siesta_1_dft_1 = pd.read_csv('{}/data_2-1.out'.format(folder_siesta_1), names=cols_siesta, delim_whitespace=True)
siesta_1_dft_2 = pd.read_csv('{}/data_2-2.out'.format(folder_siesta_1), names=cols_siesta, delim_whitespace=True)
siesta_1_V_1 = pd.read_csv('{}/data_3-1.out'.format(folder_siesta_1), names=cols_siesta, delim_whitespace=True)
siesta_1_V_2 = pd.read_csv('{}/data_3-2.out'.format(folder_siesta_1), names=cols_siesta, delim_whitespace=True)
siesta_2_dft_1 = pd.read_csv('{}/data_2-1.out'.format(folder_siesta_2), names=cols_siesta, delim_whitespace=True)
siesta_2_dft_2 = pd.read_csv('{}/data_2-2.out'.format(folder_siesta_2), names=cols_siesta, delim_whitespace=True)
siesta_2_V_1 = pd.read_csv('{}/data_3-1.out'.format(folder_siesta_2), names=cols_siesta, delim_whitespace=True)
siesta_2_V_2 = pd.read_csv('{}/data_3-2.out'.format(folder_siesta_2), names=cols_siesta, delim_whitespace=True)

# Plot z force
xlim = [0, 2]
fig_plot_1, ax_plot_1 = plt.subplots()
if plot_cp2k_1: ax_plot_1.plot(index, (cp2k_1_dft_2['X'] - cp2k_1_dft_1['X']) * param.hartree_per_bohr_to_ev_per_angstrom, 'ko-', fillstyle='none', label=labels_cp2k_1[0])
if plot_cp2k_1: ax_plot_1.plot(index, (cp2k_1_V0_2['X'] - cp2k_1_V0_1['X']) * param.hartree_per_bohr_to_ev_per_angstrom, 'ro-', fillstyle='none', label=labels_cp2k_1[1])
# if plot_cp2k_2: ax_plot_1.plot(index, (cp2k_2_dft_2['X'] - cp2k_2_dft_1['X']) * param.hartree_per_bohr_to_ev_per_angstrom, 'go-', fillstyle='none', label=labels_cp2k_2[0])
if plot_cp2k_2: ax_plot_1.plot(index, (cp2k_2_V0_2['X'] - cp2k_2_V0_1['X']) * param.hartree_per_bohr_to_ev_per_angstrom, 'bo-', fillstyle='none', label=labels_cp2k_2[1])
if plot_siesta_1: ax_plot_1.plot(index, (siesta_1_dft_2['X'] - siesta_1_dft_1['X']), 'ko--', fillstyle='none', label=labels_siesta_1[0])
if plot_siesta_1: ax_plot_1.plot(index, (siesta_1_V_2['X'] - siesta_1_V_1['X']), 'ro--', fillstyle='none', label=labels_siesta_1[1])
# if plot_siesta_2: ax_plot_1.plot(index, (siesta_2_dft_2['X'] - siesta_2_dft_1['X']), 'ko--', fillstyle='none', label=labels_siesta_2[0])
if plot_siesta_2: ax_plot_1.plot(index, (siesta_2_V_2['X'] - siesta_2_V_1['X']), 'bo--', fillstyle='none', label=labels_siesta_2[1])
ax_plot_1.set_xlabel('Atom index')
ax_plot_1.set_ylabel('Force x (eV / Å)')
if use_ylim[0]: ax_plot_1.set_ylim(ylim_xyz[0, :])
ax_plot_1.set_xlim([1, 9])
ax_plot_1.legend(frameon=False)
fig_plot_1.tight_layout()
if plot_cp2k_1: fig_plot_1.savefig('{}/force_x.png'.format(folder_cp2k_1), dpi=param.save_dpi)
if plot_siesta_1: fig_plot_1.savefig('{}/force_x.png'.format(folder_siesta_1), dpi=param.save_dpi)
if plot_cp2k_2: fig_plot_1.savefig('{}/force_x.png'.format(folder_cp2k_2), dpi=param.save_dpi)
if plot_siesta_2: fig_plot_1.savefig('{}/force_x.png'.format(folder_siesta_2), dpi=param.save_dpi)

# Plot all force
rows, cols = 3, 1
fig_plot_all, ax_plot_all = plt.subplots(rows, cols, sharex='col', sharey='row', figsize=(6, 8))
if plot_cp2k_1: ax_plot_all[0].plot(index, (cp2k_1_dft_2['X']) * param.hartree_per_bohr_to_ev_per_angstrom, 'ko-', fillstyle='none', label=labels_cp2k_1[0])
if plot_cp2k_1: ax_plot_all[0].plot(index, (cp2k_1_V0_2['X']) * param.hartree_per_bohr_to_ev_per_angstrom, 'ro-', fillstyle='none', label=labels_cp2k_1[1])
if plot_cp2k_1: ax_plot_all[1].plot(index, (cp2k_1_dft_2['Y']) * param.hartree_per_bohr_to_ev_per_angstrom, 'ko-',  fillstyle='none', label=labels_cp2k_1[0])
if plot_cp2k_1: ax_plot_all[1].plot(index, (cp2k_1_V0_2['Y']) * param.hartree_per_bohr_to_ev_per_angstrom, 'ro-', fillstyle='none',   label=labels_cp2k_1[1])
if plot_cp2k_1: ax_plot_all[2].plot(index, (cp2k_1_dft_2['Z']) * param.hartree_per_bohr_to_ev_per_angstrom, 'ko-',  fillstyle='none',  label=labels_cp2k_1[0])
if plot_cp2k_1: ax_plot_all[2].plot(index, (cp2k_1_V0_2['Z']) * param.hartree_per_bohr_to_ev_per_angstrom, 'ro-', fillstyle='none',  label=labels_cp2k_1[1])
# if plot_cp2k_2: ax_plot_all[0].plot(index, (cp2k_2_dft_2['X']) * param.hartree_per_bohr_to_ev_per_angstrom, 'ko-', fillstyle='none', label=labels_cp2k_2[0])
if plot_cp2k_2: ax_plot_all[0].plot(index, (cp2k_2_V0_2['X']) * param.hartree_per_bohr_to_ev_per_angstrom, 'bo-', fillstyle='none', label=labels_cp2k_2[1])
# if plot_cp2k_2: ax_plot_all[1].plot(index, (cp2k_2_dft_2['Y']) * param.hartree_per_bohr_to_ev_per_angstrom, 'ko-',  fillstyle='none', label=labels_cp2k_2[0])
if plot_cp2k_2: ax_plot_all[1].plot(index, (cp2k_2_V0_2['Y']) * param.hartree_per_bohr_to_ev_per_angstrom, 'bo-', fillstyle='none',   label=labels_cp2k_2[1])
# if plot_cp2k_2: ax_plot_all[2].plot(index, (cp2k_2_dft_2['Z']) * param.hartree_per_bohr_to_ev_per_angstrom, 'ko-',  fillstyle='none',  label=labels_cp2k_2[0])
if plot_cp2k_2: ax_plot_all[2].plot(index, (cp2k_2_V0_2['Z']) * param.hartree_per_bohr_to_ev_per_angstrom, 'bo-', fillstyle='none',  label=labels_cp2k_2[1])
if plot_siesta_1: ax_plot_all[0].plot(index, (siesta_1_dft_2['X']), 'ko--', fillstyle='none', label=labels_siesta_1[0])
if plot_siesta_1: ax_plot_all[0].plot(index, (siesta_1_V_2['X']), 'ro--', fillstyle='none', label=labels_siesta_1[1])
if plot_siesta_1: ax_plot_all[1].plot(index, (siesta_1_dft_2['Y'] ), 'ko--', fillstyle='none', label=labels_siesta_1[0])
if plot_siesta_1: ax_plot_all[1].plot(index, (siesta_1_V_2['Y']), 'ro--', fillstyle='none', label=labels_siesta_1[1])
if plot_siesta_1: ax_plot_all[2].plot(index, (siesta_1_dft_2['Z'] ), 'ko--', fillstyle='none', label=labels_siesta_1[0])
if plot_siesta_1: ax_plot_all[2].plot(index, (siesta_1_V_2['Z']) , 'ro--', fillstyle='none', label=labels_siesta_1[1])
# if plot_siesta_2: ax_plot_all[0].plot(index, (siesta_2_dft_2['X']), 'ko--', fillstyle='none', label=labels_siesta_2[0])
if plot_siesta_2: ax_plot_all[0].plot(index, (siesta_2_V_2['X']), 'bo--', fillstyle='none', label=labels_siesta_2[1])
# if plot_siesta_2: ax_plot_all[1].plot(index, (siesta_2_dft_2['Y'] ), 'ko--', fillstyle='none', label=labels_siesta_2[0])
if plot_siesta_2: ax_plot_all[1].plot(index, (siesta_2_V_2['Y']), 'bo--', fillstyle='none', label=labels_siesta_2[1])
# if plot_siesta_2: ax_plot_all[2].plot(index, (siesta_2_dft_2['Z'] ), 'ko--', fillstyle='none', label=labels_siesta_2[0])
if plot_siesta_2: ax_plot_all[2].plot(index, (siesta_2_V_2['Z']) , 'bo--', fillstyle='none', label=labels_siesta_2[1])
ax_plot_all[0].set_ylabel('Force x (eV / Å)')
ax_plot_all[1].set_ylabel('Force y (eV / Å)')
ax_plot_all[2].set_ylabel('Force z (eV / Å)')
ax_plot_all[2].set_xlabel('Atom index')
if use_ylim[0]: ax_plot_all[0].set_ylim(ylim_xyz[0, :])
if use_ylim[1]: ax_plot_all[1].set_ylim(ylim_xyz[1, :])
if use_ylim[2]: ax_plot_all[2].set_ylim(ylim_xyz[2, :])
ax_plot_all[0].set_xlim([1, 9])
ax_plot_all[1].set_xlim([1, 9])
ax_plot_all[2].set_xlim([1, 9])
ax_plot_all[0].legend(frameon=False)
fig_plot_all.tight_layout()
if plot_cp2k_1: fig_plot_all.savefig('{}/force_all.png'.format(folder_cp2k_1), dpi=param.save_dpi)
if plot_siesta_1: fig_plot_all.savefig('{}/force_all.png'.format(folder_siesta_1), dpi=param.save_dpi)
if plot_cp2k_2: fig_plot_all.savefig('{}/force_all.png'.format(folder_cp2k_2), dpi=param.save_dpi)
if plot_siesta_2: fig_plot_all.savefig('{}/force_all.png'.format(folder_siesta_2), dpi=param.save_dpi)

# Plot all force differences
rows, cols = 3, 1
fig_plot_all_diff, ax_plot_all_diff = plt.subplots(rows, cols, sharex='col', sharey='row', figsize=(6, 8))
if plot_cp2k_1: ax_plot_all_diff[0].plot(index, (cp2k_1_dft_2['X'] - cp2k_1_dft_1['X'])*param.hartree_per_bohr_to_ev_per_angstrom, 'ko-', fillstyle='none', label=labels_cp2k_1[0])
if plot_cp2k_1: ax_plot_all_diff[0].plot(index, (cp2k_1_V0_2['X'] - cp2k_1_V0_1['X'])*param.hartree_per_bohr_to_ev_per_angstrom, 'ro-', fillstyle='none', label=labels_cp2k_1[1])
if plot_cp2k_1: ax_plot_all_diff[1].plot(index, (cp2k_1_dft_2['Y'] - cp2k_1_dft_1['Y'])*param.hartree_per_bohr_to_ev_per_angstrom, 'ko-', fillstyle='none', label=labels_cp2k_1[0])
if plot_cp2k_1: ax_plot_all_diff[1].plot(index, (cp2k_1_V0_2['Y'] - cp2k_1_V0_1['Y'])*param.hartree_per_bohr_to_ev_per_angstrom, 'ro-', fillstyle='none', label=labels_cp2k_1[1])
if plot_cp2k_1: ax_plot_all_diff[2].plot(index, (cp2k_1_dft_2['Z'] - cp2k_1_dft_1['Z'])*param.hartree_per_bohr_to_ev_per_angstrom, 'ko-', fillstyle='none', label=labels_cp2k_1[0])
if plot_cp2k_1: ax_plot_all_diff[2].plot(index, (cp2k_1_V0_2['Z'] - cp2k_1_V0_1['Z'])*param.hartree_per_bohr_to_ev_per_angstrom, 'ro-', fillstyle='none', label=labels_cp2k_1[1])
# if plot_cp2k_2: ax_plot_all_diff[0].plot(index, (cp2k_2_dft_2['X'] - cp2k_2_dft_1['X'])*param.hartree_per_bohr_to_ev_per_angstrom, 'ko-', fillstyle='none', label=labels_cp2k_2[0])
if plot_cp2k_2: ax_plot_all_diff[0].plot(index, (cp2k_2_V0_2['X'] - cp2k_2_V0_1['X'])*param.hartree_per_bohr_to_ev_per_angstrom, 'bo-', fillstyle='none', label=labels_cp2k_2[1])
# if plot_cp2k_2: ax_plot_all_diff[1].plot(index, (cp2k_2_dft_2['Y'] - cp2k_2_dft_1['Y'])*param.hartree_per_bohr_to_ev_per_angstrom, 'ko-', fillstyle='none', label=labels_cp2k_2[0])
if plot_cp2k_2: ax_plot_all_diff[1].plot(index, (cp2k_2_V0_2['Y'] - cp2k_2_V0_1['Y'])*param.hartree_per_bohr_to_ev_per_angstrom, 'bo-', fillstyle='none', label=labels_cp2k_2[1])
# if plot_cp2k_2: ax_plot_all_diff[2].plot(index, (cp2k_2_dft_2['Z'] - cp2k_2_dft_1['Z'])*param.hartree_per_bohr_to_ev_per_angstrom, 'ko-', fillstyle='none', label=labels_cp2k_2[0])
if plot_cp2k_2: ax_plot_all_diff[2].plot(index, (cp2k_2_V0_2['Z'] - cp2k_2_V0_1['Z'])*param.hartree_per_bohr_to_ev_per_angstrom, 'bo-', fillstyle='none', label=labels_cp2k_2[1])
if plot_siesta_1: ax_plot_all_diff[0].plot(index, (siesta_1_dft_2['X'] - siesta_1_dft_1['X']), 'ko--', fillstyle='none', label=labels_siesta_1[0])
if plot_siesta_1: ax_plot_all_diff[0].plot(index, (siesta_1_V_2['X'] - siesta_1_V_1['X']), 'r--', fillstyle='none', label=labels_siesta_1[1])
if plot_siesta_1: ax_plot_all_diff[1].plot(index, (siesta_1_dft_2['Y'] - siesta_1_dft_1['Y']), 'ko--', fillstyle='none', label=labels_siesta_1[0])
if plot_siesta_1: ax_plot_all_diff[1].plot(index, (siesta_1_V_2['Y'] - siesta_1_V_1['Y']), 'r--', fillstyle='none', label=labels_siesta_1[1])
if plot_siesta_1: ax_plot_all_diff[2].plot(index, (siesta_1_dft_2['Z'] - siesta_1_dft_1['Z']), 'ko--', fillstyle='none', label=labels_siesta_1[0])
if plot_siesta_1: ax_plot_all_diff[2].plot(index, (siesta_1_V_2['Z'] - siesta_1_V_1['Z']), 'r--', fillstyle='none', label=labels_siesta_1[1])
# if plot_siesta_2: ax_plot_all_diff[0].plot(index, (siesta_2_dft_2['X'] - siesta_2_dft_1['X']), 'ko--', fillstyle='none', label=labels_siesta_2[0])
if plot_siesta_2: ax_plot_all_diff[0].plot(index, (siesta_2_V_2['X'] - siesta_2_V_1['X']), 'b--', fillstyle='none', label=labels_siesta_2[1])
# if plot_siesta_2: ax_plot_all_diff[1].plot(index, (siesta_2_dft_2['Y'] - siesta_2_dft_1['Y']), 'ko--', fillstyle='none', label=labels_siesta_2[0])
if plot_siesta_2: ax_plot_all_diff[1].plot(index, (siesta_2_V_2['Y'] - siesta_2_V_1['Y']), 'b--', fillstyle='none', label=labels_siesta_2[1])
# if plot_siesta_2: ax_plot_all_diff[2].plot(index, (siesta_2_dft_2['Z'] - siesta_2_dft_1['Z']), 'ko--', fillstyle='none', label=labels_siesta_2[0])
if plot_siesta_2: ax_plot_all_diff[2].plot(index, (siesta_2_V_2['Z'] - siesta_2_V_1['Z']), 'b--', fillstyle='none', label=labels_siesta_2[1])
ax_plot_all_diff[0].set_ylabel('Force x (eV / Å)')
ax_plot_all_diff[1].set_ylabel('Force y (eV / Å)')
ax_plot_all_diff[2].set_ylabel('Force z (eV / Å)')
ax_plot_all_diff[2].set_xlabel('Atom index')
if use_ylim[0]: ax_plot_all_diff[0].set_ylim(ylim_xyz[0, :])
if use_ylim[1]: ax_plot_all_diff[1].set_ylim(ylim_xyz[1, :])
if use_ylim[2]: ax_plot_all_diff[2].set_ylim(ylim_xyz[2, :])
ax_plot_all_diff[0].set_xlim([1, 9])
ax_plot_all_diff[1].set_xlim([1, 9])
ax_plot_all_diff[2].set_xlim([1, 9])
ax_plot_all_diff[0].legend(frameon=False)
fig_plot_all_diff.tight_layout()
if plot_cp2k_1: fig_plot_all_diff.savefig('{}/force_all_diff.png'.format(folder_cp2k_1), dpi=param.save_dpi)
if plot_siesta_1: fig_plot_all_diff.savefig('{}/force_all_diff.png'.format(folder_siesta_1), dpi=param.save_dpi)
if plot_cp2k_2: fig_plot_all_diff.savefig('{}/force_all_diff.png'.format(folder_cp2k_2), dpi=param.save_dpi)
if plot_siesta_2: fig_plot_all_diff.savefig('{}/force_all_diff.png'.format(folder_siesta_2), dpi=param.save_dpi)

if __name__ == "__main__":
    print('Finished.')
    plt.show()
