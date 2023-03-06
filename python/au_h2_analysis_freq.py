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
          '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/AuH2/transport/siesta-smeagol/phonon']
folder_save = folder
# data = [0.0, 0.1, 0.3, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.45, 1.5, 1.55, 1.6, 1.7, 1.8, 1.9]
data = [0.0, 0.1, 0.3, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3]
plotting_colors = ['r', 'g', 'b', 'm', 'grey', 'orange', 'y']
trans = [0.4, 1]

phonon_modes = np.zeros((len(folder), len(data)*6))
for i in range(len(folder)):
    phonon_modes[i] = np.genfromtxt('{}/data_freq.out'.format(folder[i]))

# Plot all
rows, cols = 2, 2
xlim = [0, 2]
ylim = [0.9, 1.20]
fig_plot_all, ax_plot_all = plt.subplots(rows, cols, sharex='col', sharey='row', figsize=(10, 8))

for i in range(len(folder)):
    ax_plot_all[0, 0].plot(data, phonon_modes[i, ::6]/phonon_modes[i, 0], 'bo-', alpha=trans[i], fillstyle='none', label='{} Tx'.format(labels[i]))
    ax_plot_all[0, 0].plot(data, phonon_modes[i, 1::6]/phonon_modes[i, 1], 'ro-', alpha=trans[i], fillstyle='none', label='{} Ty'.format(labels[i]))
    ax_plot_all[0, 1].plot(data, phonon_modes[i, 5::6]/phonon_modes[i, 5], 'bo-', alpha=trans[i], fillstyle='none', label='{}'.format(labels[i]))
    ax_plot_all[1, 0].plot(data, phonon_modes[i, 2::6]/phonon_modes[i, 2], 'bo-', alpha=trans[i], fillstyle='none', label='{} Rx'.format(labels[i]))
    ax_plot_all[1, 0].plot(data, phonon_modes[i, 3::6]/phonon_modes[i, 3], 'ro-', alpha=trans[i], fillstyle='none', label='{} Ry'.format(labels[i]))
    ax_plot_all[1, 1].plot(data, phonon_modes[i, 4::6]/phonon_modes[i, 4], 'bo-', alpha=trans[i], fillstyle='none', label='{}'.format(labels[i]))

ax_plot_all[0, 0].set_xlim([xlim[0], xlim[1]])
ax_plot_all[0, 0].set_ylim([ylim[0], ylim[1]])
ax_plot_all[0, 0].legend(frameon=False)
ax_plot_all[1, 0].set_xlim([xlim[0], xlim[1]])
ax_plot_all[1, 0].set_ylim([ylim[0], ylim[1]])
ax_plot_all[1, 0].legend(frameon=False)
ax_plot_all[0, 1].set_xlim([xlim[0], xlim[1]])
ax_plot_all[0, 1].set_ylim([ylim[0], ylim[1]])
ax_plot_all[0, 1].legend(frameon=False)
ax_plot_all[1, 1].set_xlim([xlim[0], xlim[1]])
ax_plot_all[1, 1].set_ylim([ylim[0], ylim[1]])
ax_plot_all[1, 1].legend(frameon=False)

ax_plot_all[1, 0].set_xlabel('Bias / V')
ax_plot_all[1, 1].set_xlabel('Bias / V')
ax_plot_all[1, 0].set_ylabel(r'$\mathrm{\omega}$(V)/$\mathrm{\omega}$(0)')
ax_plot_all[0, 0].set_ylabel(r'$\mathrm{\omega}$(V)/$\mathrm{\omega}$(0)')

fig_plot_all.tight_layout()
# fig_plot_all.subplots_adjust(hspace=0)
for i in range(len(folder)):
    fig_plot_all.savefig('{}/phonon_modes.png'.format(folder_save[i]), dpi=param.save_dpi)

if __name__ == "__main__":
    print('Finished.')
    plt.show()
