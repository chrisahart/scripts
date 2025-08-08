import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import MDAnalysis as mda
from MDAnalysis.analysis import distances
from general import parameters as param
from ase.io import read, write
# from general import load_coordinates
# from general import print_xyz
from collections import OrderedDict
import csv

"""
    Plot energy and forces for bulk hematite
"""


def read_output(folder, filename):
    """
        Return CP2K MD .ener file as re-structured Numpy array.
    """

    files = ['{}/{}'.format(folder, filename)]
    cols = ['Folder', 'Energy_first', 'Energy_last', 'Gap_alpha', 'Gap_beta', 
            'HOMO_alpha_last', 'HOMO_beta_last', 'HOMO_alpha_first', 'HOMO_beta_first', 
            'LUMO_alpha_last', 'LUMO_beta_last', 'LUMO_alpha_first', 'LUMO_beta_first']
    file_energy = pd.read_csv(files[0], delim_whitespace=True, names=cols, skiprows=1)

    return file_energy


files = ['output.out'] * 5
num_atoms = 324
box_size = [13.8, 13.8, 17.76, 90, 90, 90]

# Hubbard U
folder_use_hubbard = ['pbe-u-o']
colors_hubbard = ['k']
labels_hubbard = ['']
hubbard_values_polaron = [[0, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0]]
skip_folders_hubbard_start = 1
skip_folders_hubbard_end = -8
data_hubbard_neutral = []
data_hubbard_vertical_from_neutral = []
data_hubbard_electron_offset = []
data_hubbard_vertical_from_electron = []
trapping_energy_hubbard = []
folder_files_hubbard = []
for i in range(len(folder_use_hubbard)):
    folder_files_hubbard.append(['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/resources/wiki/wiki_cp2k/tio2/anatase/cell-441/mp-390/neutral/{}'.format(folder_use_hubbard[i]),
                                 '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/resources/wiki/wiki_cp2k/tio2/anatase/cell-441/mp-390/hole/from-neutral/{}'.format(folder_use_hubbard[i]),
                                 '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/resources/wiki/wiki_cp2k/tio2/anatase/cell-441/mp-390/hole/from-md-hole/{}'.format(folder_use_hubbard[i]),
                                 '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/resources/wiki/wiki_cp2k/tio2/anatase/cell-441/mp-390/neutral-from-hole/{}'.format(folder_use_hubbard[i])])
    data_hubbard_neutral.append(read_output(folder_files_hubbard[i][0], files[0]))
    data_hubbard_vertical_from_neutral.append(read_output(folder_files_hubbard[i][1], files[0]))
    data_hubbard_electron_offset.append(read_output(folder_files_hubbard[i][2], files[0]))
    data_hubbard_vertical_from_electron.append(read_output(folder_files_hubbard[i][3], files[0]))
    trapping_energy_hubbard.append(data_hubbard_vertical_from_neutral[i]['Energy_first'][:np.shape(hubbard_values_polaron[i])[0]] - data_hubbard_electron_offset[i]['Energy_last'][:np.shape(hubbard_values_polaron[i])[0]])
folder_save = folder_files_hubbard[0]

# Plot band gap Hubbard U
draw_legend = True
set_ylim = [2, 3.4]
draw_gap = 3.2
fig_band_gap1, ax_band_gap1 = plt.subplots()
ax_band_gap1.axhline(y=draw_gap, color='k', alpha=0.5)
for i in range(len(folder_use_hubbard)):
    ax_band_gap1.plot(hubbard_values_polaron[i], data_hubbard_neutral[i]['Gap_alpha'], 'x-', color=colors_hubbard[i], label='{}'.format(labels_hubbard[i]))
ax_band_gap1.set_ylabel('Band gap / eV')
ax_band_gap1.set_xlabel('Hubbard U / eV')
ax_band_gap1.set_ylim(set_ylim[0], set_ylim[1])
if draw_legend: ax_band_gap1.legend(frameon=False)
fig_band_gap1.tight_layout()
for i in range(np.shape(folder_save)[0]):
    fig_band_gap1.savefig('{}/band_gap_hubbard.png'.format(folder_save[i]), dpi=300)

# Plot trapping energy Hubbard U
fig_trapping1, ax_trapping1 = plt.subplots()
for i in range(len(folder_use_hubbard)):
    ax_trapping1.plot(hubbard_values_polaron[i][skip_folders_hubbard_start:skip_folders_hubbard_end],
                      trapping_energy_hubbard[i][skip_folders_hubbard_start:skip_folders_hubbard_end]*param.hartree_to_ev*1e3, 'x-', color=colors_hubbard[i], label='{}'.format(labels_hubbard[i]))
ax_trapping1.set_ylabel('Trapping energy / meV')
ax_trapping1.set_xlabel('Hubbard U / eV')
if draw_legend: ax_trapping1.legend(frameon=False)
fig_trapping1.tight_layout()
for i in range(np.shape(folder_save)[0]):
    fig_trapping1.savefig('{}/trapping_hubbard.png'.format(folder_save[i]), dpi=300)


if __name__ == "__main__":
    print('Finished.')
    plt.show()

