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


def calc_distance(x1, y1, z1, x2, y2, z2):
    distance = np.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2 + (z1 - z2) ** 2)
    return distance


# Clotilde SIESTA, SIESTA
labels = ['[1]', 'SIESTA-SMEAGOL', 'CP2K-SMEAGOL']
folder = ['/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/AuH2/clotilde/positive',
          '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/AuH2/transport/siesta-smeagol/geo_opt/bias']
# folder_cp2k_smeagol = ['/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/AuH2/transport/cp2k-smeagol/archer/geo_opt/bias_nodes-4']
folder_cp2k_smeagol = ['/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/forces/AuH2/archer_lda/geo_opt/bias']
file_cp2k_smeagol = 'V-pos-1.xyz'
file_siesta = ['/opt/final.siesta',
               '/final.siesta']
file_cp2k = ['/opt/final.xyz',
             '/final.xyz']
data = [0.0, 0.1, 0.3, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.45, 1.5, 1.55, 1.6, 1.7, 1.8, 1.9]
# data_cp2k = [0.0, 0.1, 0.3, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,]
data_cp2k = [0.0, 0.1, 0.3, 0.5, 0.7, 0.9, 1.1]

plotting_colors = ['r', 'b', 'g', 'm', 'grey', 'orange', 'y']

# Convert from final.siesta to final.xyz
for j in range(0, len(folder)):
    for i in range(len(data)):
        print(j, i)
        xyz_siesta_to_cp2k.siesta_to_cp2k('{}/{}/{}'.format(folder[j], data[i], file_siesta[j]), '{}/{}/{}'.format(folder[j], data[i], file_cp2k[j]))

# Calculate bond lengths and collective variable
file_2 = []
cols = ['Species', 'X', 'Y', 'Z']
atoms = np.array([59, 60, 61, 62, 63, 64]) - 1
bond_length_siesta = np.zeros((len(folder), len(data), atoms.shape[0]-1))
z_coordinate_siesta = np.zeros((len(folder), len(data), atoms.shape[0]-1))
z_bond_change_siesta = np.zeros((len(folder), len(data), atoms.shape[0]-1))
collective_variable_siesta = np.zeros((len(folder), len(data)))

# SIESTA
for k in range(0, len(folder)):
    for i in range(len(data)):
        file_coord_siesta, _, _ = load_coordinates.load_file_coord(folder[k], '{}/{}'.format(data[i], file_cp2k[k]), cols)
        for j in range(atoms.shape[0]-1):
            z_coordinate_siesta[k, i, j] = file_coord_siesta['Z'][atoms[j]]
            bond_length_siesta[k, i, j] = calc_distance(file_coord_siesta['X'][atoms[j]], file_coord_siesta['Y'][atoms[j]], file_coord_siesta['Z'][atoms[j]],
                                              file_coord_siesta['X'][atoms[j+1]], file_coord_siesta['Y'][atoms[j+1]], file_coord_siesta['Z'][atoms[j+1]])
        # print(i, bond_length_siesta[k, i, :])
    z_bond_change_siesta[k] = z_coordinate_siesta[k] - z_coordinate_siesta[k, 0]
    collective_variable_siesta[k] = np.mean(z_bond_change_siesta[k], axis=1)

# CP2K
bond_length_cp2k = np.zeros((len(folder_cp2k_smeagol), len(data_cp2k), atoms.shape[0]-1))
z_coordinate_cp2k = np.zeros((len(folder_cp2k_smeagol), len(data_cp2k), atoms.shape[0]-1))
z_bond_change_cp2k = np.zeros((len(folder_cp2k_smeagol), len(data_cp2k), atoms.shape[0]-1))
collective_variable_cp2k = np.zeros((len(folder_cp2k_smeagol), len(data_cp2k)))
for k in range(0, len(folder_cp2k_smeagol)):
    for i in range(len(data_cp2k)):
        coord_dft, _, _, _, _, num_atoms, _ = load_coordinates.load_values_coord(folder_cp2k_smeagol[k], '{}/{}'.format(data_cp2k[i], file_cp2k_smeagol), cols)
        # file_coord_cp2k, _, _ = load_coordinates.load_file_coord(folder_cp2k_smeagol[k], '{}/{}'.format(data_cp2k[i], file_cp2k_smeagol), cols)
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
xlim = [0, 2]
ylim = [0.865, 0.91]
atom_pair = 2
fig_plot_1, ax_plot_1 = plt.subplots()
for k in range(0, len(folder)):
    ax_plot_1.plot(data, bond_length_siesta[k, :, atom_pair], 'o-', fillstyle='none', color=plotting_colors[k], label=labels[k])
ax_plot_1.plot(data_cp2k, bond_length_cp2k[0, :, atom_pair], 'o-', fillstyle='none', color=plotting_colors[2], label=labels[2])
ax_plot_1.set_xlim([xlim[0], xlim[1]])
ax_plot_1.set_ylim([ylim[0], ylim[1]])
ax_plot_1.set_xlabel('Bias / V')
ax_plot_1.set_ylabel('H-H bond length / Å')
ax_plot_1.legend(frameon=False)
fig_plot_1.tight_layout()
for k in range(0, len(folder)):
    fig_plot_1.savefig('{}/h2_bond_length_all.png'.format(folder[k]), dpi=param.save_dpi)
fig_plot_1.savefig('{}/h2_bond_length_all.png'.format(folder_cp2k_smeagol[0]), dpi=param.save_dpi)

# Plot H-H bond length with bias
xlim = [0, 2]
ylim = [-0.0001, 0.035]
fig_plot_3, ax_plot_3 = plt.subplots()
for k in range(0, len(folder)):
    ax_plot_3.plot(data, bond_length_siesta[k, :, 2] - bond_length_siesta[k, 0, 2], 'o-', fillstyle='none', color=plotting_colors[k], label=labels[k])
ax_plot_3.plot(data_cp2k, bond_length_cp2k[0, :, 2] - bond_length_cp2k[0, 0, 2], 'o-', fillstyle='none', color=plotting_colors[2], label=labels[2])
ax_plot_3.set_xlim([xlim[0], xlim[1]])
ax_plot_3.set_ylim([ylim[0], ylim[1]])
ax_plot_3.set_xlabel('Bias / V')
ax_plot_3.set_ylabel('H-H bond length / Å')
ax_plot_3.legend(frameon=False)
fig_plot_3.tight_layout()
for k in range(0, len(folder)):
    fig_plot_3.savefig('{}/h2_bond_length_change_all.png'.format(folder[k]), dpi=param.save_dpi)
fig_plot_3.savefig('{}/h2_bond_length_change_all.png'.format(folder_cp2k_smeagol[0]), dpi=param.save_dpi)

# Plot change of the collective variable with bias
ylim = [-0.1, 0.001]
fig_plot_2, ax_plot_2 = plt.subplots()
for k in range(0, len(folder)):
    ax_plot_2.plot(data, collective_variable_siesta[k], 'o-', fillstyle='none', color=plotting_colors[k], label=labels[k])
ax_plot_2.plot(data_cp2k, collective_variable_cp2k[0, :], 'o-', fillstyle='none', color=plotting_colors[2], label=labels[2])
ax_plot_2.set_xlim([xlim[0], xlim[1]])
ax_plot_2.set_ylim([ylim[0], ylim[1]])
ax_plot_2.set_xlabel('Bias / V')
ax_plot_2.set_ylabel('Mean displacement / Å')
ax_plot_2.legend(frameon=False)
fig_plot_2.tight_layout()
for k in range(0, len(folder)):
    fig_plot_2.savefig('{}/collective_variable_all.png'.format(folder[k]), dpi=param.save_dpi)
fig_plot_2.savefig('{}/collective_variable_all.png'.format(folder_cp2k_smeagol[0]), dpi=param.save_dpi)

if __name__ == "__main__":
    print('Finished.')
    plt.show()
