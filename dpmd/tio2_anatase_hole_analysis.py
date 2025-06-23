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
# folder_use_hubbard = ['pbe-u-o']
folder_use_hubbard = ['pbe-u-o', 'pbe-u-o-ti-6']
colors_hubbard = ['r', 'b']
labels_hubbard = ['Ti U = 0', 'Ti U = 6']
hubbard_values_polaron = [[0, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0],
                          [0, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0]]
# folder_use_hubbard = ['pbe-u-o', 'pbe-u-o-ti-3', 'pbe-u-o-ti-6']
# colors_hubbard = ['r', 'g', 'b']
# labels_hubbard = ['Ti U = 0', 'Ti U = 3', 'Ti U = 6']
# hubbard_values_polaron = [[0, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0],
# [0, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0],
#                           [0, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0]]
skip_folders_hubbard_start = 1
skip_folders_hubbard_end = -8
data_hubbard_neutral = []
data_hubbard_vertical_from_neutral = []
data_hubbard_electron_offset = []
data_hubbard_vertical_from_electron = []
trapping_energy_hubbard = []
non_linearity_alpha_hubbard = []
folder_files_hubbard = []
ionisation_potential_hubbard_e = []
non_linearity_alpha_hubbard_e = []
non_linearity_beta_hubbard_e = []
ionisation_potential_hubbard_h = []
non_linearity_alpha_hubbard_h = []
non_linearity_beta_hubbard_h = []
for i in range(len(folder_use_hubbard)):
    folder_files_hubbard.append(['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/anatase/archer/anatase/cell-441/geo-opt/mp-390/neutral/{}'.format(folder_use_hubbard[i]),
                                 '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/anatase/archer/anatase/cell-441/geo-opt/mp-390/hole/from-neutral/{}'.format(folder_use_hubbard[i]),
                                 '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/anatase/archer/anatase/cell-441/geo-opt/mp-390/hole/from-md-hole/{}'.format(folder_use_hubbard[i]),
                                 '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/anatase/archer/anatase/cell-441/geo-opt/mp-390/neutral-from-hole/{}'.format(folder_use_hubbard[i])])
    data_hubbard_neutral.append(read_output(folder_files_hubbard[i][0], files[0]))
    data_hubbard_vertical_from_neutral.append(read_output(folder_files_hubbard[i][1], files[0]))
    data_hubbard_electron_offset.append(read_output(folder_files_hubbard[i][2], files[0]))
    data_hubbard_vertical_from_electron.append(read_output(folder_files_hubbard[i][3], files[0]))
    trapping_energy_hubbard.append(data_hubbard_vertical_from_neutral[i]['Energy_first'][:np.shape(hubbard_values_polaron[i])[0]] - data_hubbard_electron_offset[i]['Energy_last'][:np.shape(hubbard_values_polaron[i])[0]])
    ionisation_potential_hubbard_e.append(-(data_hubbard_electron_offset[i]['Energy_last'][:np.shape(hubbard_values_polaron[i])[0]] - data_hubbard_vertical_from_electron[i]['Energy_first'][:np.shape(hubbard_values_polaron[i])[0]]))
    non_linearity_alpha_hubbard_e.append(ionisation_potential_hubbard_e[i] - data_hubbard_electron_offset[i]['HOMO_alpha_last'][:np.shape(hubbard_values_polaron[i])[0]])
    non_linearity_beta_hubbard_e.append(ionisation_potential_hubbard_e[i] - data_hubbard_electron_offset[i]['HOMO_beta_last'][:np.shape(hubbard_values_polaron[i])[0]])
    ionisation_potential_hubbard_h.append(data_hubbard_neutral[i]['Energy_last'][:np.shape(hubbard_values_polaron[i])[0]] - (data_hubbard_vertical_from_neutral[i]['Energy_first'][:np.shape(hubbard_values_polaron[i])[0]]))
    # non_linearity_alpha_hubbard_h.append(ionisation_potential_hubbard_h[i] - data_hubbard_vertical_from_electron[i]['HOMO_alpha_first'][:np.shape(hubbard_values_polaron[i])[0]])
    # non_linearity_beta_hubbard_h.append(ionisation_potential_hubbard_h[i] - data_hubbard_vertical_from_electron[i]['HOMO_beta_first'][:np.shape(hubbard_values_polaron[i])[0]])
    non_linearity_alpha_hubbard_h.append(ionisation_potential_hubbard_h[i] - data_hubbard_neutral[i]['HOMO_alpha_last'][:np.shape(hubbard_values_polaron[i])[0]])
    non_linearity_beta_hubbard_h.append(ionisation_potential_hubbard_h[i] - data_hubbard_neutral[i]['HOMO_beta_last'][:np.shape(hubbard_values_polaron[i])[0]])
# print(ionisation_potential_hse_h[0])
# print(data_hse_vertical_from_electron[0]['HOMO_alpha_first'][:np.shape(hse_values_polaron[0])[0]])
# print(non_linearity_alpha_hse[0])

# HSE
# folder_use_hse = ['hfx-schwarz-e-1e-4-f-1e-6-cfit11-cpfit3', 'hfx-schwarz-e-1e-4-f-1e-6-fit11-cpfit3',
#                   'hfx-schwarz-e-1e-6-f-1e-6-fit9-cpfit3', 'hfx-schwarz-e-1e-6-f-1e-6-fit9-pfit3']
# colors_hse = ['r', 'g', 'b', 'm']
# labels_hse = ['1', '2', '3', '4']
folder_use_hse = ['hfx-schwarz-e-1e-6-f-1e-6-fit9-pfit3']
colors_hse = ['r']
labels_hse = ['HSE']
hse_values_polaron = [[12, 15, 18, 21, 22, 23, 24, 25, 30]]
skip_folders_hse_start = 0
skip_folders_hse_end = -1
data_hse_neutral = []
data_hse_vertical_from_neutral = []
data_hse_electron_offset = []
data_hse_vertical_from_electron = []
trapping_energy_hse = []
non_linearity_alpha_hse = []
folder_files_hse = []
ionisation_potential_hse_e = []
non_linearity_alpha_hse_e = []
non_linearity_beta_hse_e = []
ionisation_potential_hse_h = []
non_linearity_alpha_hse_h = []
non_linearity_beta_hse_h = []
for i in range(len(folder_use_hse)):
    folder_files_hse.append(['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/anatase/archer/anatase/cell-441/geo-opt/mp-390/neutral/{}'.format(folder_use_hse[i]),
                                 '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/anatase/archer/anatase/cell-441/geo-opt/mp-390/hole/from-neutral/{}'.format(folder_use_hse[i]),
                                 '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/anatase/archer/anatase/cell-441/geo-opt/mp-390/hole/from-md-hole/{}'.format(folder_use_hse[i]),
                                 '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/anatase/archer/anatase/cell-441/geo-opt/mp-390/neutral-from-hole/{}'.format(folder_use_hse[i])])
    data_hse_neutral.append(read_output(folder_files_hse[i][0], files[0]))
    data_hse_vertical_from_neutral.append(read_output(folder_files_hse[i][1], files[0]))
    data_hse_electron_offset.append(read_output(folder_files_hse[i][2], files[0]))
    data_hse_vertical_from_electron.append(read_output(folder_files_hse[i][3], files[0]))
    trapping_energy_hse.append(data_hse_vertical_from_neutral[i]['Energy_first'][:np.shape(hse_values_polaron[i])[0]] - data_hse_electron_offset[i]['Energy_last'][:np.shape(hse_values_polaron[i])[0]])
    ionisation_potential_hse_e.append(-(data_hse_electron_offset[i]['Energy_last'][:np.shape(hse_values_polaron[i])[0]] - data_hse_vertical_from_electron[i]['Energy_first'][:np.shape(hse_values_polaron[i])[0]]))
    non_linearity_alpha_hse_e.append(ionisation_potential_hse_e[i] - data_hse_electron_offset[i]['HOMO_alpha_last'][:np.shape(hse_values_polaron[i])[0]])
    non_linearity_beta_hse_e.append(ionisation_potential_hse_e[i] - data_hse_electron_offset[i]['HOMO_beta_last'][:np.shape(hse_values_polaron[i])[0]])
    ionisation_potential_hse_h.append(data_hse_neutral[i]['Energy_last'][:np.shape(hse_values_polaron[i])[0]] - (data_hse_vertical_from_neutral[i]['Energy_first'][:np.shape(hse_values_polaron[i])[0]]))
    # non_linearity_alpha_hse_h.append(ionisation_potential_hse_h[i] - data_hse_vertical_from_electron[i]['HOMO_alpha_first'][:np.shape(hse_values_polaron[i])[0]])
    # non_linearity_beta_hse_h.append(ionisation_potential_hse_h[i] - data_hse_vertical_from_electron[i]['HOMO_beta_first'][:np.shape(hse_values_polaron[i])[0]])
    non_linearity_alpha_hse_h.append(ionisation_potential_hse_h[i] - data_hse_neutral[i]['HOMO_alpha_last'][:np.shape(hse_values_polaron[i])[0]])
    non_linearity_beta_hse_h.append(ionisation_potential_hse_h[i] - data_hse_neutral[i]['HOMO_beta_last'][:np.shape(hse_values_polaron[i])[0]])
# non_linearity_hse = non_linearity_alpha_hse_h
# print('non_linearity_alpha_hse_h', non_linearity_alpha_hse_h[0])
# print('ionisation_potential_hse_h', ionisation_potential_hse_h[0])
# print('data_hse_neutral', data_hse_neutral[0]['HOMO_beta_last'][:np.shape(hse_values_polaron[0])[0]])
# print('data_hse_vertical_from_electron', data_hse_vertical_from_electron[0]['HOMO_beta_first'][:np.shape(hse_values_polaron[0])[0]])
# print('data_hse_vertical_from_electron', data_hse_vertical_from_electron[0]['HOMO_alpha_first'][:np.shape(hse_values_polaron[0])[0]])

# pbe0
folder_use_pbe0 = ['pbe0-tz-tz2p-hfx-schwarz-e-1e-6-f-1e-6-fit9-pfit3-mckenna-input']
colors_pbe0 = ['g']
labels_pbe0 = ['PBE0']
pbe0_values_polaron = [[10, 10.5, 11, 11.5, 12, 13, 14, 15, 25]]
skip_folders_pbe0_start = 0
skip_folders_pbe0_end = -1
data_pbe0_neutral = []
data_pbe0_vertical_from_neutral = []
data_pbe0_electron_offset = []
data_pbe0_vertical_from_electron = []
trapping_energy_pbe0 = []
non_linearity_alpha_pbe0 = []
folder_files_pbe0 = []
ionisation_potential_pbe0_e = []
non_linearity_alpha_pbe0_e = []
non_linearity_beta_pbe0_e = []
ionisation_potential_pbe0_h = []
non_linearity_alpha_pbe0_h = []
non_linearity_beta_pbe0_h = []

electron_affinity_pbe0_h = []
non_linearity2_alpha_pbe0_h = []
non_linearity2_beta_pbe0_h = []

for i in range(len(folder_use_pbe0)):
    folder_files_pbe0.append(['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/anatase/archer/anatase/cell-441/geo-opt/mp-390/neutral/{}'.format(folder_use_pbe0[i]),
                                 '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/anatase/archer/anatase/cell-441/geo-opt/mp-390/hole/from-neutral/{}'.format(folder_use_pbe0[i]),
                                 '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/anatase/archer/anatase/cell-441/geo-opt/mp-390/hole/from-md-hole/{}'.format(folder_use_pbe0[i]),
                                 '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/anatase/archer/anatase/cell-441/geo-opt/mp-390/neutral-from-hole/{}'.format(folder_use_pbe0[i])])
    data_pbe0_neutral.append(read_output(folder_files_pbe0[i][0], files[0]))
    data_pbe0_vertical_from_neutral.append(read_output(folder_files_pbe0[i][1], files[0]))
    data_pbe0_electron_offset.append(read_output(folder_files_pbe0[i][2], files[0]))
    data_pbe0_vertical_from_electron.append(read_output(folder_files_pbe0[i][3], files[0]))
    trapping_energy_pbe0.append(data_pbe0_vertical_from_neutral[i]['Energy_first'][:np.shape(pbe0_values_polaron[i])[0]] - data_pbe0_electron_offset[i]['Energy_last'][:np.shape(pbe0_values_polaron[i])[0]])

    ionisation_potential_pbe0_e.append(-(data_pbe0_electron_offset[i]['Energy_last'][:np.shape(pbe0_values_polaron[i])[0]] - data_pbe0_vertical_from_electron[i]['Energy_first'][:np.shape(pbe0_values_polaron[i])[0]]))
    # non_linearity_alpha_pbe0_e.append(ionisation_potential_pbe0_e[i] - data_pbe0_electron_offset[i]['HOMO_alpha_last'][:np.shape(pbe0_values_polaron[i])[0]])
    # non_linearity_beta_pbe0_e.append(ionisation_potential_pbe0_e[i] - data_pbe0_electron_offset[i]['HOMO_beta_last'][:np.shape(pbe0_values_polaron[i])[0]])

    non_linearity_alpha_pbe0_e.append(ionisation_potential_pbe0_e[i] - data_pbe0_electron_offset[i]['LUMO_alpha_last'][:np.shape(pbe0_values_polaron[i])[0]])
    non_linearity_beta_pbe0_e.append(ionisation_potential_pbe0_e[i] - data_pbe0_electron_offset[i]['LUMO_beta_last'][:np.shape(pbe0_values_polaron[i])[0]])


    ionisation_potential_pbe0_h.append(data_pbe0_neutral[i]['Energy_last'][:np.shape(pbe0_values_polaron[i])[0]] - (data_pbe0_vertical_from_neutral[i]['Energy_first'][:np.shape(pbe0_values_polaron[i])[0]]))
    # non_linearity_alpha_pbe0_h.append(ionisation_potential_pbe0_h[i] - data_pbe0_vertical_from_electron[i]['HOMO_alpha_first'][:np.shape(pbe0_values_polaron[i])[0]])
    # non_linearity_beta_pbe0_h.append(ionisation_potential_pbe0_h[i] - data_pbe0_vertical_from_electron[i]['HOMO_beta_first'][:np.shape(pbe0_values_polaron[i])[0]])
    non_linearity_alpha_pbe0_h.append(ionisation_potential_pbe0_h[i] - data_pbe0_neutral[i]['HOMO_alpha_last'][:np.shape(pbe0_values_polaron[i])[0]])
    non_linearity_beta_pbe0_h.append(ionisation_potential_pbe0_h[i] - data_pbe0_neutral[i]['HOMO_beta_last'][:np.shape(pbe0_values_polaron[i])[0]])

# non_linearity_pbe0 = non_linearity_alpha_pbe0_h
# print('non_linearity_alpha_pbe0_h', non_linearity_alpha_pbe0_h[0])
# print('ionisation_potential_pbe0_h', ionisation_potential_pbe0_h[0])
# print('data_pbe0_neutral', data_pbe0_neutral[0]['HOMO_beta_last'][:np.shape(pbe0_values_polaron[0])[0]])
# print('data_pbe0_vertical_from_electron', data_pbe0_vertical_from_electron[0]['HOMO_beta_first'][:np.shape(pbe0_values_polaron[0])[0]])
# print('data_pbe0_vertical_from_electron', data_pbe0_vertical_from_electron[0]['HOMO_alpha_first'][:np.shape(pbe0_values_polaron[0])[0]])

# PBE0
# folder_use_pbe0 = ['pbe0-tz-tz2p-hfx-schwarz-e-1e-6-f-1e-6-fit9-pfit3',
#                    'pbe0-tz-tz2p-hfx-schwarz-e-1e-6-f-1e-6-fit9-pfit3-mckenna-input']
# folder_use_pbe0 = folder_use_pbe0[0]
# folder_files_pbe0 = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-336/geo-opt/convergence/neutral/neutral-opt/{}'.format(folder_use_pbe0),
#                      '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-336/geo-opt/convergence/electron/from-neutral-opt/{}'.format(folder_use_pbe0),
#                      '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-336/geo-opt/convergence/electron/from-u-4.0/{}'.format(folder_use_pbe0),
#                      '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-336/geo-opt/convergence/neutral/from-electron-opt/{}'.format(folder_use_pbe0)]
# pbe0_values_polaron = [10, 10.5, 11, 11.5, 12, 13, 14, 15, 25]
# skip_folders_pbe0 = 4
# plot_folders_pbe0 = [4, ]
# data_pbe0_neutral = read_output(folder_files_pbe0[0], files[0])
# data_pbe0_vertical_from_neutral = read_output(folder_files_pbe0[1], files[0])
# data_pbe0_electron_offset = read_output(folder_files_pbe0[2], files[0])
# data_pbe0_vertical_from_electron = read_output(folder_files_pbe0[3], files[0])
# trapping_energy_pbe0 = data_pbe0_vertical_from_neutral['Energy_first'][:np.shape(pbe0_values_polaron)[0]] - data_pbe0_electron_offset['Energy_last'][:np.shape(pbe0_values_polaron)[0]]
# ionisation_potential_pbe0 = (data_pbe0_electron_offset['Energy_last'][:np.shape(pbe0_values_polaron)[0]] - data_pbe0_vertical_from_electron['Energy_first'][:np.shape(pbe0_values_polaron)[0]])
# non_linearity_alpha_pbe0 = ionisation_potential_pbe0 - data_pbe0_electron_offset['HOMO_alpha_last'][:np.shape(pbe0_values_polaron)[0]]  # E(polaron) - E(neutral) - E_polaron(HOMO_alpha)
# non_linearity_beta_pbe0 = ionisation_potential_pbe0 - data_pbe0_electron_offset['HOMO_beta_last'][:np.shape(pbe0_values_polaron)[0]]  # E(polaron) - E(neutral) - E_polaron(HOMO_beta)

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

# Plot band gap HFX
set_ylim = [2.6, 3.8]
fig_band_gap2, ax_band_gap2 = plt.subplots()
ax_band_gap2.axhline(y=draw_gap, color='k', alpha=0.5)
for i in range(len(folder_use_hse)):
    ax_band_gap2.plot(hse_values_polaron[i], data_hse_neutral[i]['Gap_alpha'], 'x-', color=colors_hse[i], label='{}'.format(labels_hse[i]))
    ax_band_gap2.plot(pbe0_values_polaron[i], data_pbe0_neutral[i]['Gap_alpha'], 'x-', color=colors_pbe0[i], label='{}'.format(labels_pbe0[i]))
# ax_band_gap2.plot(pbe0_values_polaron, data_pbe0_neutral['Gap_alpha'], 'gx-', label='PBE0')
ax_band_gap2.set_ylabel('Band gap / eV')
ax_band_gap2.set_xlabel('% HFX')
ax_band_gap2.set_ylim(set_ylim[0], set_ylim[1])
fig_band_gap2.tight_layout()
if draw_legend: ax_band_gap2.legend(frameon=False)
for i in range(np.shape(folder_save)[0]):
    fig_band_gap2.savefig('{}/band_gap_hfx.png'.format(folder_save[i]), dpi=300)

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

# Plot trapping energy HFX
set_ylim = [100, 500]
fig_trapping2, ax_trapping2 = plt.subplots()
for i in range(len(folder_use_hse)):
    ax_trapping2.plot(hse_values_polaron[i][skip_folders_hse_start:skip_folders_hse_end],
                      trapping_energy_hse[i][skip_folders_hse_start:skip_folders_hse_end] * param.hartree_to_ev*1e3, 'x-', color=colors_hse[i], label='{}'.format(labels_hse[i]))
    ax_trapping2.plot(pbe0_values_polaron[i][skip_folders_pbe0_start:skip_folders_pbe0_end],
                      trapping_energy_pbe0[i][skip_folders_pbe0_start:skip_folders_pbe0_end] * param.hartree_to_ev*1e3, 'x-', color=colors_pbe0[i], label='{}'.format(labels_pbe0[i]))
# ax_trapping2.plot(pbe0_values_polaron[skip_folders_pbe0:], trapping_energy_pbe0[skip_folders_pbe0:] * param.hartree_to_ev*1e3, 'gx-', label='PBE0')
ax_trapping2.set_ylabel('Trapping energy / meV')
ax_trapping2.set_xlabel('% HFX')
ax_trapping2.set_ylim(set_ylim[0], set_ylim[1])
if draw_legend: ax_trapping2.legend(frameon=False)
fig_trapping2.tight_layout()
for i in range(np.shape(folder_save)[0]):
    fig_trapping2.savefig('{}/trapping_hfx.png'.format(folder_save[i]), dpi=300)

# Plot Hubbard U
draw_nonlinearity = 0.05
set_ylim = [-0.2, 0.5]
set_ylim = [0, 0.12]
fig_linearity1, ax_linearity1 = plt.subplots()
ax_linearity1.axhline(y=draw_nonlinearity, color='k', alpha=0.5)
for i in range(len(folder_use_hubbard)):
    # ax_linearity1.plot(hubbard_values_polaron[i][skip_folders_hubbard_start:skip_folders_hubbard_end],
    #                    non_linearity_alpha_hubbard[i][skip_folders_hubbard_start:skip_folders_hubbard_end] * param.hartree_to_ev, 'x-', color=colors_hubbard[i], label='{}'.format(labels_hubbard[i]))
    ax_linearity1.plot(hubbard_values_polaron[i][skip_folders_hubbard_start:skip_folders_hubbard_end],
                       non_linearity_beta_hubbard_h[i][skip_folders_hubbard_start:skip_folders_hubbard_end] * param.hartree_to_ev, 'x-', color=colors_hubbard[i], label='{}'.format(labels_hubbard[i]))
# ax_linearity1.plot(hubbard_values_polaron[skip_folders_hubbard:], non_linearity_beta_hubbard_h[skip_folders_hubbard:] * param.hartree_to_ev, 'kx-', label='PBE+U')
ax_linearity1.set_ylabel('Non-linearity / eV')
ax_linearity1.set_xlabel('Hubbard U / eV')
if draw_legend: ax_linearity1.legend(frameon=False)
ax_linearity1.set_ylim(set_ylim[0], set_ylim[1])
fig_linearity1.tight_layout()
for i in range(np.shape(folder_save)[0]):
    fig_linearity1.savefig('{}/linearity_hubbard.png'.format(folder_save[i]), dpi=300)

# Plot non-linearity HFX
set_ylim = [0, 0.12]
fig_linearity2, ax_linearity2 = plt.subplots()
ax_linearity2.axhline(y=draw_nonlinearity, color='k', alpha=0.5)
for i in range(len(folder_use_hse)):
    # ax_trapping2.plot(hse_values_polaron[i][skip_folders_hse_start:skip_folders_hse_end],
    #                   trapping_energy_hse[i][skip_folders_hse_start:skip_folders_hse_end] * param.hartree_to_ev*1e3, 'x-', color=colors_hse[i], label='{}'.format(labels_hse[i]))
    # ax_linearity2.plot(hse_values_polaron[i][skip_folders_hse_start:skip_folders_hse_end],
    #                    non_linearity_beta_hse_h[i][skip_folders_hse_start:skip_folders_hse_end] * param.hartree_to_ev, 'x-', color=colors_hse[i], label='{}'.format(labels_hse[i]))
    # ax_linearity2.plot(pbe0_values_polaron[i][skip_folders_pbe0_start:skip_folders_pbe0_end],
    #                    non_linearity_alpha_pbe0_h[i][skip_folders_pbe0_start:skip_folders_pbe0_end] * param.hartree_to_ev, 'x-', color=colors_pbe0[i], label='{}'.format(labels_pbe0[i]))
    # ax_linearity2.plot(pbe0_values_polaron[i][skip_folders_pbe0_start:skip_folders_pbe0_end],
    #                    non_linearity_beta_pbe0_h[i][skip_folders_pbe0_start:skip_folders_pbe0_end] * param.hartree_to_ev, 'x-', color=colors_pbe0[i], label='{}'.format(labels_pbe0[i]))
    ax_linearity2.plot(pbe0_values_polaron[i][skip_folders_pbe0_start:skip_folders_pbe0_end],
                       non_linearity_alpha_pbe0_e[i][skip_folders_pbe0_start:skip_folders_pbe0_end] * param.hartree_to_ev, 'r-',  label='{}'.format(labels_pbe0[i]))
    ax_linearity2.plot(pbe0_values_polaron[i][skip_folders_pbe0_start:skip_folders_pbe0_end],
                       non_linearity_beta_pbe0_e[i][skip_folders_pbe0_start:skip_folders_pbe0_end] * param.hartree_to_ev, 'g-', label='{}'.format(labels_pbe0[i]))
    ax_linearity2.plot(pbe0_values_polaron[i][skip_folders_pbe0_start:skip_folders_pbe0_end],
                       non_linearity_alpha_pbe0_h[i][skip_folders_pbe0_start:skip_folders_pbe0_end] * param.hartree_to_ev, 'b-', label='{}'.format(labels_pbe0[i]))
    ax_linearity2.plot(pbe0_values_polaron[i][skip_folders_pbe0_start:skip_folders_pbe0_end],
                       non_linearity_beta_pbe0_h[i][skip_folders_pbe0_start:skip_folders_pbe0_end] * param.hartree_to_ev, 'm-', label='{}'.format(labels_pbe0[i]))



# ax_linearity2.plot(pbe0_values_polaron[skip_folders_pbe0:], non_linearity_alpha_pbe0[skip_folders_pbe0:] * param.hartree_to_ev, 'gx-', label='PBE0')
ax_linearity2.set_ylabel('Non-linearity / eV')
ax_linearity2.set_ylim(set_ylim[0], set_ylim[1])
ax_linearity2.set_xlabel('% HFX')
if draw_legend: ax_linearity2.legend(frameon=False)
fig_linearity2.tight_layout()
for i in range(np.shape(folder_save)[0]):
    fig_linearity2.savefig('{}/linearity_hfx.png'.format(folder_save[i]), dpi=300)

if __name__ == "__main__":
    print('Finished.')
    plt.show()

