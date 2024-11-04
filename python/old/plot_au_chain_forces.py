from __future__ import division, print_function
import pandas as pd
import numpy as np
from general import load_coordinates
from general import parameters as param
import matplotlib.pyplot as plt

"""
    Analysis script for au chain as in Rouxing paper
"""

folder_cp2k_1 = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/au-chain/cp2k-smeagol/kpoints-2-2-20_constrained_energy-density-T_contour-double'
folder_cp2k_2 = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/au-chain/cp2k-smeagol/kpoints-2-2-20_constrained_energy-density-F_contour-double'
folder_siesta_1 = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/au-chain/siesta1-smeagol/kpoints-2-2-20_energy-density-T'
folder_siesta_2 = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/au-chain/siesta1-smeagol/kpoints-2-2-20_energy-density-F'
# folder_siesta_1 = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/au-chain/siesta3-smeagol/kpoints-2-2-20_constrained_positive_Au-2.8_alex_2'

labels_cp2k_1 = ['CP2K', 'CP2K+SMEAGOL']
labels_cp2k_2 = ['CP2K', 'CP2K+SMEAGOL Ω=0']
labels_siesta_1 = ['SIESTA', 'SIESTA+SMEAGOL']
labels_siesta_2 = ['SIESTA', 'SIESTA+SMEAGOL Ω=0']
colors = ['k', 'r', 'b']
# labels_siesta_1 = ['SIESTA3', 'SIESTA3-SMEAGOL V=0 Ω=0']
# colors = ['k', 'b', 'r']
plot_cp2k_1 = True
plot_cp2k_2 = True
plot_siesta_1 = False
plot_siesta_2 = False
# plot_cp2k_1 = False
# plot_cp2k_2 = False
# plot_siesta_1 = True
# plot_siesta_2 = True
use_ylim = [False, True, False]
legend_loc='lower left'
ylim_xyz = np.array([(-2, 2), (-1, 1), (-1, 1)])

# Read CP2K
cols = ['Species', 'X', 'Y', 'Z']
cp2k_1_dft_1, num_atoms, _ = load_coordinates.load_file_coord(folder_cp2k_1, 'au-chain_szv_dftenergy-frc-1.xyz', cols)
if plot_cp2k_1: cp2k_1_V0_1, _, _ = load_coordinates.load_file_coord(folder_cp2k_1, 'au-chain_szv_zerobias-frc-1.xyz', cols)
if plot_cp2k_1: cp2k_1_dft_2, _, _ = load_coordinates.load_file_coord(folder_cp2k_1, 'au-chain_szv_shiftx_dftenergy-frc-1.xyz', cols)
if plot_cp2k_1: cp2k_1_V0_2, _, _ = load_coordinates.load_file_coord(folder_cp2k_1, 'au-chain_szv_shiftx_zerobias-frc-1.xyz', cols)
if plot_cp2k_2: cp2k_2_dft_1, _, _ = load_coordinates.load_file_coord(folder_cp2k_2, 'au-chain_szv_dftenergy-frc-1.xyz', cols)
if plot_cp2k_2: cp2k_2_V0_1, _, _ = load_coordinates.load_file_coord(folder_cp2k_2, 'au-chain_szv_zerobias-frc-1.xyz', cols)
if plot_cp2k_2: cp2k_2_dft_2, _, _ = load_coordinates.load_file_coord(folder_cp2k_2, 'au-chain_szv_shiftx_dftenergy-frc-1.xyz', cols)
if plot_cp2k_2: cp2k_2_V0_2, _, _ = load_coordinates.load_file_coord(folder_cp2k_2, 'au-chain_szv_shiftx_zerobias-frc-1.xyz', cols)
index = np.arange(start=1, stop=num_atoms+1)

# Read SIESTA
cols_siesta_force = ['A', 'X', 'Y', 'Z']
cols_siesta_force_constrain = ['A', 'B', 'X', 'Y', 'Z']
if plot_siesta_1: siesta_1_dft_1 = pd.read_csv('{}/data_2-1.out'.format(folder_siesta_1), names=cols_siesta_force, delim_whitespace=True)
if plot_siesta_1: siesta_1_dft_2 = pd.read_csv('{}/data_2-2.out'.format(folder_siesta_1), names=cols_siesta_force_constrain, delim_whitespace=True)
if plot_siesta_1: siesta_1_V_1 = pd.read_csv('{}/data_3-1.out'.format(folder_siesta_1), names=cols_siesta_force, delim_whitespace=True)
if plot_siesta_1: siesta_1_V_2 = pd.read_csv('{}/data_3-2.out'.format(folder_siesta_1), names=cols_siesta_force_constrain, delim_whitespace=True)
if plot_siesta_2: siesta_2_dft_1 = pd.read_csv('{}/data_2-1.out'.format(folder_siesta_2), names=cols_siesta_force, delim_whitespace=True)
if plot_siesta_2: siesta_2_dft_2 = pd.read_csv('{}/data_2-2.out'.format(folder_siesta_2), names=cols_siesta_force_constrain, delim_whitespace=True)
if plot_siesta_2: siesta_2_V_1 = pd.read_csv('{}/data_3-1.out'.format(folder_siesta_2), names=cols_siesta_force, delim_whitespace=True)
if plot_siesta_2: siesta_2_V_2 = pd.read_csv('{}/data_3-2.out'.format(folder_siesta_2), names=cols_siesta_force_constrain, delim_whitespace=True)

# Plot component of force
xlim = [0, 2]
component = 'X'
fig_plot_1, ax_plot_1 = plt.subplots(figsize=(6, 4))
if plot_cp2k_1: ax_plot_1.plot(index, (cp2k_1_dft_2[component] - cp2k_1_dft_1[component]) * param.hartree_per_bohr_to_ev_per_angstrom, 'o-', color=colors[0], fillstyle='none', label=labels_cp2k_1[0])
if plot_cp2k_1: ax_plot_1.plot(index, (cp2k_1_V0_2[component] - cp2k_1_V0_1[component]) * param.hartree_per_bohr_to_ev_per_angstrom, 'o--', color=colors[1], fillstyle='none', label=labels_cp2k_1[1])
if plot_cp2k_2: ax_plot_1.plot(index, (cp2k_2_V0_2[component] - cp2k_2_V0_1[component]) * param.hartree_per_bohr_to_ev_per_angstrom, 'o--', color=colors[2], fillstyle='none', label=labels_cp2k_2[1])
if plot_siesta_1: ax_plot_1.plot(index, (siesta_1_dft_2[component] - siesta_1_dft_1[component]), 'o-', color=colors[0], fillstyle='none', label=labels_siesta_1[0])
if plot_siesta_1: ax_plot_1.plot(index, (siesta_1_V_2[component] - siesta_1_V_1[component]), 'o--', color=colors[1], fillstyle='none', label=labels_siesta_1[1])
if plot_siesta_2: ax_plot_1.plot(index, (siesta_2_V_2[component] - siesta_2_V_1[component]), 'o--', color=colors[2], fillstyle='none', label=labels_siesta_2[1])
ax_plot_1.set_xlabel('Atom index')
ax_plot_1.set_ylabel('Force x (eV / Å)')
if use_ylim[0]: ax_plot_1.set_ylim(ylim_xyz[0, :])
ax_plot_1.set_xlim([1, 9])
ax_plot_1.legend(frameon=False, loc=legend_loc)
fig_plot_1.tight_layout()
if plot_cp2k_1: fig_plot_1.savefig('{}/force_{}.png'.format(folder_cp2k_1, component), dpi=param.save_dpi)
if plot_siesta_1: fig_plot_1.savefig('{}/force_{}.png'.format(folder_siesta_1, component), dpi=param.save_dpi)
if plot_cp2k_2: fig_plot_1.savefig('{}/force_{}.png'.format(folder_cp2k_2, component), dpi=param.save_dpi)
if plot_siesta_2: fig_plot_1.savefig('{}/force_{}.png'.format(folder_siesta_2, component), dpi=param.save_dpi)

# Plot all force
rows, cols = 3, 1
fig_plot_all, ax_plot_all = plt.subplots(rows, cols, sharex='col', sharey='row', figsize=(6, 8))
if plot_cp2k_1: ax_plot_all[0].plot(index, (cp2k_1_dft_2['X']) * param.hartree_per_bohr_to_ev_per_angstrom, 'o-', color=colors[0], fillstyle='none', label=labels_cp2k_1[0])
if plot_cp2k_1: ax_plot_all[0].plot(index, (cp2k_1_V0_2['X']) * param.hartree_per_bohr_to_ev_per_angstrom, 'o--', color=colors[1], fillstyle='none', label=labels_cp2k_1[1])
if plot_cp2k_1: ax_plot_all[1].plot(index, (cp2k_1_dft_2['Y']) * param.hartree_per_bohr_to_ev_per_angstrom, 'o-', color=colors[0],  fillstyle='none', label=labels_cp2k_1[0])
if plot_cp2k_1: ax_plot_all[1].plot(index, (cp2k_1_V0_2['Y']) * param.hartree_per_bohr_to_ev_per_angstrom, 'o--', color=colors[1], fillstyle='none',   label=labels_cp2k_1[1])
if plot_cp2k_1: ax_plot_all[2].plot(index, (cp2k_1_dft_2['Z']) * param.hartree_per_bohr_to_ev_per_angstrom, 'o-', color=colors[0],  fillstyle='none',  label=labels_cp2k_1[0])
if plot_cp2k_1: ax_plot_all[2].plot(index, (cp2k_1_V0_2['Z']) * param.hartree_per_bohr_to_ev_per_angstrom, 'o--', color=colors[1], fillstyle='none',  label=labels_cp2k_1[1])
if plot_cp2k_2: ax_plot_all[0].plot(index, (cp2k_2_V0_2['X']) * param.hartree_per_bohr_to_ev_per_angstrom, 'o--', color=colors[2], fillstyle='none', label=labels_cp2k_2[1])
if plot_cp2k_2: ax_plot_all[1].plot(index, (cp2k_2_V0_2['Y']) * param.hartree_per_bohr_to_ev_per_angstrom, 'o--', color=colors[2], fillstyle='none',   label=labels_cp2k_2[1])
if plot_cp2k_2: ax_plot_all[2].plot(index, (cp2k_2_V0_2['Z']) * param.hartree_per_bohr_to_ev_per_angstrom, 'o--', color=colors[2], fillstyle='none',  label=labels_cp2k_2[1])
if plot_siesta_1: ax_plot_all[0].plot(index, (siesta_1_dft_2['X']), 'o-', color=colors[0], fillstyle='none', label=labels_siesta_1[0])
if plot_siesta_1: ax_plot_all[0].plot(index, (siesta_1_V_2['X']), 'o--', color=colors[1], fillstyle='none', label=labels_siesta_1[1])
if plot_siesta_1: ax_plot_all[1].plot(index, (siesta_1_dft_2['Y']), 'o-', color=colors[0], fillstyle='none', label=labels_siesta_1[0])
if plot_siesta_1: ax_plot_all[1].plot(index, (siesta_1_V_2['Y']), 'o--', color=colors[1], fillstyle='none', label=labels_siesta_1[1])
if plot_siesta_1: ax_plot_all[2].plot(index, (siesta_1_dft_2['Z']), 'o-', color=colors[0], fillstyle='none', label=labels_siesta_1[0])
if plot_siesta_1: ax_plot_all[2].plot(index, (siesta_1_V_2['Z']), 'o--', color=colors[1], fillstyle='none', label=labels_siesta_1[1])
if plot_siesta_2: ax_plot_all[0].plot(index, (siesta_2_V_2['X']), 'o--', color=colors[2], fillstyle='none', label=labels_siesta_2[1])
if plot_siesta_2: ax_plot_all[1].plot(index, (siesta_2_V_2['Y']), 'o--', color=colors[2], fillstyle='none', label=labels_siesta_2[1])
if plot_siesta_2: ax_plot_all[2].plot(index, (siesta_2_V_2['Z']), 'o--', color=colors[2], fillstyle='none', label=labels_siesta_2[1])
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
ax_plot_all[0].legend(frameon=False, loc=legend_loc)
fig_plot_all.tight_layout()
if plot_cp2k_1: fig_plot_all.savefig('{}/force_all.png'.format(folder_cp2k_1), dpi=param.save_dpi)
if plot_siesta_1: fig_plot_all.savefig('{}/force_all.png'.format(folder_siesta_1), dpi=param.save_dpi)
if plot_cp2k_2: fig_plot_all.savefig('{}/force_all.png'.format(folder_cp2k_2), dpi=param.save_dpi)
if plot_siesta_2: fig_plot_all.savefig('{}/force_all.png'.format(folder_siesta_2), dpi=param.save_dpi)

if __name__ == "__main__":
    print('Finished.')
    plt.show()
