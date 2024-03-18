from __future__ import division, print_function
import pandas as pd
import numpy as np
import glob
from general import load_coordinates
from general import parameters as param
from general import print_xyz
import matplotlib.pyplot as plt
import csv
import xyz_siesta_to_cp2k

"""
    Analysis script for Au-H2-Au system
"""

params = {'axes.formatter.limits': [-4, 4],
          'axes.labelsize': 'x-large',
          'axes.titlesize': 'x-large',
          'legend.fontsize': 'medium',
          'lines.markersize': '8',
          }
plt.rcParams.update(params)



def calc_distance(x1, y1, z1, x2, y2, z2):
    distance = np.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2 + (z1 - z2) ** 2)
    return distance


# SIESTA+SMEAGOL
# folder_siesta_smeagol = ['/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/au-h2/siesta1-smeagol-reference']
# file_siesta = ['/opt/final.siesta']
# file_cp2k = ['/opt/final.xyz']
folder_siesta_smeagol = ['/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/au-h2/siesta1-smeagol/geo_opt/NEnergReal-64_kpoints-3-3-20_omp-1_ParallelOverKNum-1']
file_siesta = ['/final.siesta']
file_cp2k = ['/final.xyz']

# CP2k+SMEAGOL
# folder_cp2k_smeagol = ['/Volumes/ELEMENTS/Storage/Postdoc/data/Work/Postdoc/Work/calculations/transport/2023/AuH2/transport/cp2k-smeagol/archer/archer_pbe/geo_opt/NEnergReal-64_kpoints-3-3-20_omp-1_ParallelOverKNum-1']
# folder_cp2k_smeagol = ['/Volumes/ELEMENTS/Storage/Postdoc/data/Work/Postdoc/Work/calculations/transport/2023/AuH2/transport/cp2k-smeagol/archer/archer_pbe/geo_opt/NEnergReal-64_kpoints-3-3-20_omp-2_ParallelOverKNum-3_scf-1e-6_EXTRAPOLATION-USE_PREV_P']
# folder_cp2k_smeagol = ['/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/au-h2/cp2k-smeagol/geo_opt/kpoints-3-3-20_omp-2_ParallelOverKNum-3_contour-single_dynamic-14',
#                        '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/au-h2/cp2k-smeagol/geo_opt/kpoints-3-3-20_omp-2_ParallelOverKNum-3_contour-single_dynamic-14_leads-6s',
#                        ]
folder_cp2k_smeagol = ['/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/au-h2/cp2k-smeagol/geo_opt/kpoints-3-3-20_omp-2_ParallelOverKNum-3_contour-single_dynamic-14']
file_cp2k_smeagol = ['V-pos-1.xyz'] * 5


# General
# labels = ['SIESTA1+SMEAGOL EM.WeightRho 0.5', 'CP2K+SMEAGOL EM.WeightRho 0.5']
# labels = ['SIESTA+SMEAGOL', 'CP2K+SMEAGOL']
# labels = ['Bai et al.', 'SIESTA1+SMEAGOL EM.WeightRho 0.5', 'CP2K+SMEAGOL EM.WeightRho 0.5', 'CP2K+SMEAGOL weighted double contour']
# labels = ['SIESTA1+SMEAGOL Bai et al.', 'CP2K+SMEAGOL EM.WeightRho 0.5', 'CP2K+SMEAGOL weighted double contour']
# labels = ['SIESTA1+SMEAGOL Bai et al.', 'CP2K+SMEAGOL EM.WeightRho 0.5', 'CP2K+SMEAGOL weighted double contour']
# labels = ['SIESTA1+SMEAGOL EM.WeightRho 0.5', 'CP2K+SMEAGOL EM.WeightRho 0.5', 'CP2K+SMEAGOL weighted double contour']
labels = ['SIESTA+SMEAGOL', 'CP2K+SMEAGOL']
data_cp2k_smeagol = [0.0, 0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9]
file_prefix_cp2k_smeagol = data_cp2k_smeagol
data_siesta_smeagol = data_cp2k_smeagol
plotting_colors = ['b', 'r', 'orange', 'm']
ylim_collective = [-0.103, 0.004]
ylim_collective = [-0.20, 0.004]
ylim_collective = [-0.18, 0.004]
# ylim_h2 = [0.865, 0.92]
ylim_h2 = [0.865, 0.91]
ylim_h2_change = [-0.0001, 0.035]
xlim = [-0.05, 2]

# Convert from final.siesta to final.xyz
for j in range(0, len(folder_siesta_smeagol)):
    for i in range(len(data_siesta_smeagol)):
        print(j, i)
        print(folder_siesta_smeagol[j])
        print(data_siesta_smeagol[i])
        xyz_siesta_to_cp2k.siesta_label_to_cp2k('{}/{}/{}'.format(folder_siesta_smeagol[j], data_siesta_smeagol[i], file_siesta[j]),
                                                '{}/{}/{}'.format(folder_siesta_smeagol[j], data_siesta_smeagol[i], file_cp2k[j]),
                                                cols=['X', 'Y', 'Z', 'Species_num', 'Species', 'B', 'C'])

# Calculate bond lengths and collective variable
file_2 = []
cols = ['Species', 'X', 'Y', 'Z']
atoms = np.array([59, 60, 61, 62, 63, 64]) - 1
bond_length_siesta = np.zeros((len(folder_siesta_smeagol), len(data_siesta_smeagol), atoms.shape[0]-1))
z_coordinate_siesta = np.zeros((len(folder_siesta_smeagol), len(data_siesta_smeagol), atoms.shape[0]-1))
z_bond_change_siesta = np.zeros((len(folder_siesta_smeagol), len(data_siesta_smeagol), atoms.shape[0]-1))
collective_variable_siesta = np.zeros((len(folder_siesta_smeagol), len(data_siesta_smeagol)))

# SIESTA
for k in range(0, len(folder_siesta_smeagol)):
    for i in range(len(data_siesta_smeagol)):
        file_coord_siesta, _, _ = load_coordinates.load_file_coord(folder_siesta_smeagol[k], '{}/{}'.format(data_siesta_smeagol[i], file_cp2k[k]), cols)
        for j in range(atoms.shape[0]-1):
            z_coordinate_siesta[k, i, j] = file_coord_siesta['Z'][atoms[j]]
            bond_length_siesta[k, i, j] = calc_distance(file_coord_siesta['X'][atoms[j]], file_coord_siesta['Y'][atoms[j]], file_coord_siesta['Z'][atoms[j]],
                                              file_coord_siesta['X'][atoms[j+1]], file_coord_siesta['Y'][atoms[j+1]], file_coord_siesta['Z'][atoms[j+1]])
        # print(i, bond_length_siesta[k, i, :])
    z_bond_change_siesta[k] = z_coordinate_siesta[k] - z_coordinate_siesta[k, 0]
    collective_variable_siesta[k] = np.mean(z_bond_change_siesta[k], axis=1)

# CP2K
bond_length_cp2k = np.zeros((len(folder_cp2k_smeagol), len(data_cp2k_smeagol), atoms.shape[0]-1))
z_coordinate_cp2k = np.zeros((len(folder_cp2k_smeagol), len(data_cp2k_smeagol), atoms.shape[0]-1))
z_bond_change_cp2k = np.zeros((len(folder_cp2k_smeagol), len(data_cp2k_smeagol), atoms.shape[0]-1))
collective_variable_cp2k = np.zeros((len(folder_cp2k_smeagol), len(data_cp2k_smeagol)))
for k in range(0, len(folder_cp2k_smeagol)):
    for i in range(len(data_cp2k_smeagol)):
        coord_dft, _, _, _, _, num_atoms, _ = load_coordinates.load_values_coord(folder_cp2k_smeagol[k], '{}/{}{}'.format(data_cp2k_smeagol[i], file_prefix_cp2k_smeagol[i], file_cp2k_smeagol[k]), cols)
        # file_coord_cp2k, _, _ = load_coordinates.load_file_coord(folder_cp2k_smeagol[k], '{}/{}'.format(data_cp2k_smeagol[i], file_cp2k_smeagol), cols)
        print(coord_dft.shape)
        for j in range(atoms.shape[0]-1):
            z_coordinate_cp2k[k, i, j] = coord_dft[-1, 2, atoms[j]]
            bond_length_cp2k[k, i, j] = calc_distance(coord_dft[-1, 0, atoms[j]], coord_dft[-1, 1, atoms[j]], coord_dft[-1, 2, atoms[j]],
                                              coord_dft[-1, 0, atoms[j+1]], coord_dft[-1, 1, atoms[j+1]], coord_dft[-1, 2, atoms[j+1]])
        # print(i, bond_length_cp2k[k, i, :])
    z_bond_change_cp2k[k] = z_coordinate_cp2k[k] - z_coordinate_cp2k[k, 0]
    collective_variable_cp2k[k] = np.mean(z_bond_change_cp2k[k], axis=1)

# print(bond_length_siesta.shape)
# print(bond_length_cp2k.shape)
print(bond_length_siesta)
print(bond_length_cp2k)

# Plot H-H bond length with bias
atom_pair = 2
fig_plot_1, ax_plot_1 = plt.subplots()
for k in range(0, len(folder_siesta_smeagol)):
    ax_plot_1.plot(data_siesta_smeagol, bond_length_siesta[k, :, atom_pair], 'o-', fillstyle='none', color=plotting_colors[k], label=labels[k])
for m in range(0, len(folder_cp2k_smeagol)):
    ax_plot_1.plot(data_cp2k_smeagol, bond_length_cp2k[m, :, atom_pair], 'o-', fillstyle='none', color=plotting_colors[k+k], label=labels[k+k])
ax_plot_1.set_xlim([xlim[0], xlim[1]])
ax_plot_1.set_ylim([ylim_h2[0], ylim_h2[1]])
ax_plot_1.set_xlabel('Bias / V')
ax_plot_1.set_ylabel('H-H bond length / Å')
ax_plot_1.legend(frameon=False)
fig_plot_1.tight_layout()
for k in range(0, len(folder_siesta_smeagol)):
    fig_plot_1.savefig('{}/h2_bond_length_all.png'.format(folder_siesta_smeagol[k]), dpi=param.save_dpi)
fig_plot_1.savefig('{}/h2_bond_length_all.png'.format(folder_cp2k_smeagol[0]), dpi=param.save_dpi)

# Plot H-H bond length with bias
fig_plot_3, ax_plot_3 = plt.subplots()
for k in range(0, len(folder_siesta_smeagol)):
    ax_plot_3.plot(data_siesta_smeagol, bond_length_siesta[k, :, 2] - bond_length_siesta[k, 0, 2], 'o-', fillstyle='none', color=plotting_colors[k], label=labels[k])
for m in range(0, len(folder_cp2k_smeagol)):
    ax_plot_3.plot(data_cp2k_smeagol, bond_length_cp2k[m, :, 2] - bond_length_cp2k[m, 0, 2], 'o-', fillstyle='none', color=plotting_colors[m+k], label=labels[m+k])
ax_plot_3.set_xlim([xlim[0], xlim[1]])
ax_plot_3.set_ylim([ylim_h2_change[0], ylim_h2_change[1]])
ax_plot_3.set_xlabel('Bias / V')
ax_plot_3.set_ylabel('H-H bond length / Å')
ax_plot_3.legend(frameon=False)
fig_plot_3.tight_layout()
for k in range(0, len(folder_siesta_smeagol)):
    fig_plot_3.savefig('{}/h2_bond_length_change_all.png'.format(folder_siesta_smeagol[k]), dpi=param.save_dpi)
fig_plot_3.savefig('{}/h2_bond_length_change_all.png'.format(folder_cp2k_smeagol[0]), dpi=param.save_dpi)

# Plot change of the collective variable with bias
fig_plot_2, ax_plot_2 = plt.subplots()
for m in range(0, len(folder_cp2k_smeagol)):
    ax_plot_2.plot(data_cp2k_smeagol, collective_variable_cp2k[m, :], 'o-', fillstyle='none', color=plotting_colors[m+len(folder_siesta_smeagol)], label=labels[m+len(folder_siesta_smeagol)])
for k in range(0, len(folder_siesta_smeagol)):
    ax_plot_2.plot(data_siesta_smeagol, collective_variable_siesta[k], 'o-', fillstyle='none', color=plotting_colors[k], label=labels[k])
ax_plot_2.set_xlim([xlim[0], xlim[1]])
ax_plot_2.set_ylim([ylim_collective[0], ylim_collective[1]])
ax_plot_2.set_xlabel('Bias / V')
ax_plot_2.set_ylabel('Mean displacement / Å')
ax_plot_2.legend(frameon=False)
fig_plot_2.tight_layout()
for k in range(0, len(folder_siesta_smeagol)):
    fig_plot_2.savefig('{}/collective_variable_all.png'.format(folder_siesta_smeagol[k]), dpi=param.save_dpi)
fig_plot_2.savefig('{}/collective_variable_all.png'.format(folder_cp2k_smeagol[0]), dpi=param.save_dpi)

# Plot both
rows, cols = 2, 1
fig_cube_both, ax_cube_z = plt.subplots(rows, cols, sharex='col', sharey='row', figsize=(6, 8))
# fig_cube_both, ax_cube_z = plt.subplots(rows, cols, sharex='col', sharey='row', figsize=(7, 8))
for m in range(0, len(folder_cp2k_smeagol)):
    ax_cube_z[1].plot(data_cp2k_smeagol, bond_length_cp2k[m, :, atom_pair], 'o-', fillstyle='none', color=plotting_colors[m+len(folder_siesta_smeagol)], label=labels[m+len(folder_siesta_smeagol)])
for k in range(0, len(folder_siesta_smeagol)):
    ax_cube_z[1].plot(data_siesta_smeagol, bond_length_siesta[k, :, atom_pair], 'o-', fillstyle='none', color=plotting_colors[k], label=labels[k])
ax_cube_z[1].set_xlim([xlim[0], xlim[1]])
ax_cube_z[1].set_ylim([ylim_h2[0], ylim_h2[1]])
ax_cube_z[1].set_ylabel('H-H bond length / Å')
ax_cube_z[1].set_xlabel('Bias / V')
# ax_cube_z[1].legend(frameon=False)
for m in range(0, len(folder_cp2k_smeagol)):
    ax_cube_z[0].plot(data_cp2k_smeagol, collective_variable_cp2k[m, :], 'o-', fillstyle='none', color=plotting_colors[m+len(folder_siesta_smeagol)], label=labels[m+len(folder_siesta_smeagol)])
for k in range(0, len(folder_siesta_smeagol)):
    ax_cube_z[0].plot(data_siesta_smeagol, collective_variable_siesta[k], 'o-', fillstyle='none', color=plotting_colors[k], label=labels[k])
ax_cube_z[0].set_xlim([xlim[0], xlim[1]])
ax_cube_z[0].set_ylim([ylim_collective[0], ylim_collective[1]])
ax_cube_z[0].set_ylabel('Mean displacement / Å')
ax_cube_z[0].legend(frameon=False)
fig_cube_both.tight_layout()
for k in range(0, len(folder_siesta_smeagol)):
    fig_cube_both.savefig('{}/both.png'.format(folder_siesta_smeagol[k]), dpi=param.save_dpi)
for m in range(0, len(folder_cp2k_smeagol)):
    fig_cube_both.savefig('{}/both.png'.format(folder_cp2k_smeagol[m]), dpi=param.save_dpi)

if __name__ == "__main__":
    print('Finished.')
    plt.show()
