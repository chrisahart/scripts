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
labels = ['[1]', 'SIESTA-SMEAGOL']
folder = ['/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/AuH2/clotilde/positive',
          '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/AuH2/transport/siesta-smeagol/geo_opt/bias']
file_siesta = ['/opt/final.siesta',
               '/final.siesta']
file_cp2k = ['/opt/final.xyz',
             '/final.xyz']
data = [0.0, 0.1, 0.3, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.45, 1.5, 1.55, 1.6, 1.7, 1.8, 1.9]
# data = [0.0, 0.1, 0.3, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.45, 1.5, 1.55, 1.6, 1.9]
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
bond_length = np.zeros((len(folder), len(data), atoms.shape[0]-1))
z_coordinate = np.zeros((len(folder), len(data), atoms.shape[0]-1))
z_bond_change = np.zeros((len(folder), len(data), atoms.shape[0]-1))
collective_variable = np.zeros((len(folder), len(data)))
for k in range(0, len(folder)):
    for i in range(len(data)):
        file_coord, num_atoms, species = load_coordinates.load_file_coord(folder[k], '{}/{}'.format(data[i], file_cp2k[k]), cols)
        for j in range(atoms.shape[0]-1):
            z_coordinate[k, i, j] = file_coord['Z'][atoms[j]]
            bond_length[k, i, j] = calc_distance(file_coord['X'][atoms[j]], file_coord['Y'][atoms[j]], file_coord['Z'][atoms[j]],
                                              file_coord['X'][atoms[j+1]], file_coord['Y'][atoms[j+1]], file_coord['Z'][atoms[j+1]])
        print(i, bond_length[k, i, :])
    z_bond_change[k] = z_coordinate[k] - z_coordinate[k, 0]
    collective_variable[k] = np.mean(z_bond_change[k], axis=1)
print(z_bond_change)
print(collective_variable)

# Plot change of the H-H bond length with bias
xlim = [0, 2]
fig_plot_1, ax_plot_1 = plt.subplots()
for k in range(0, len(folder)):
    ax_plot_1.plot(data, bond_length[k, :, 2], 'o-', fillstyle='none', color=plotting_colors[k], label=labels[k])
ax_plot_1.set_xlim([xlim[0], xlim[1]])
ax_plot_1.set_xlabel('Bias / V')
ax_plot_1.set_ylabel('H-H bond length / Å')
ax_plot_1.legend(frameon=False)
fig_plot_1.tight_layout()
for k in range(0, len(folder)):
    fig_plot_1.savefig('{}/h2_bond_lengths.png'.format(folder[k]), dpi=param.save_dpi)

# Plot change of the H-H bond length with bias
ylim = [-0.1, 0]
fig_plot_2, ax_plot_2 = plt.subplots()
for k in range(0, len(folder)):
    ax_plot_2.plot(data, collective_variable[k], 'o-', fillstyle='none', color=plotting_colors[k], label=labels[k])
ax_plot_2.set_xlim([xlim[0], xlim[1]])
ax_plot_2.set_ylim([ylim[0], ylim[1]])
ax_plot_2.set_xlabel('Bias / V')
ax_plot_2.set_ylabel('Mean displacement / Å')
ax_plot_2.legend(frameon=False)
fig_plot_2.tight_layout()
for k in range(0, len(folder)):
    fig_plot_2.savefig('{}/collective_variable.png'.format(folder[k]), dpi=param.save_dpi)

if __name__ == "__main__":
    print('Finished.')
    plt.show()
