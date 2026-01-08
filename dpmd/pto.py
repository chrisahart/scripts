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

# Plotting
draw_nonlinearity = 0.1
set_ylim_linearity_hubbard = [-0.2, 0.2]
set_ylim_linearity_hfx = [-0.2, 0.2]
set_ylim_band_gap_hubbard = [1.7, 3.2]
# set_ylim_band_gap_hubbard = [2, 3.4]
draw_gap = 3.03  # rutile
# draw_gap = 3.2   # anatase
set_ylim_band_gap_hfx = [2.4, 3.4]
set_ylim_trapping_hfx = [0, 200]

# Plotting
draw_nonlinearity = 0.05
set_ylim_linearity_hubbard = [-0.3, 0.3]
set_ylim_linearity_hfx = [-0.07, 0.20]
set_ylim_linearity_hfx = [-0.001, 0.15]
set_ylim_band_gap_hubbard = [1.7, 3.2]
# set_ylim_band_gap_hubbard = [2, 3.4]
draw_gap = 3.03  # rutile
# draw_gap = 3.2   # anatase
set_ylim_band_gap_hfx = [2.4, 3.4]
set_ylim_band_gap_hfx = [2.9, 3.3]
set_ylim_trapping_hfx = [0, 210]

folder_use_hubbard = []
folder_files_hubbard = []
folder_files_pbe0 = []
folder_use_pbe0 = []
folder_files_hse = []
folder_use_hse = []

# cell 3x3x12
colors_hubbard = ['k']
labels_hubbard = ['']
skip_folders_hubbard_start = 0
skip_folders_hubbard_end = -6
folder_use_hubbard = []
folder_files_hubbard = []
colors_hse = ['k']
labels_hse = ['']
skip_folders_hse_start = 2
skip_folders_hse_end = -1
skip_folders_hse_end_gap = -1
folder_use_hse = ['hfx-schwarz-e-1e-6-f-1e-6-fit9-pfit3']
folder_files_hse = []
for i in range(len(folder_use_hse)):
    folder_files_hse.append(['/Volumes/Samsung/Data/Postdoc2/Data/Work/calculations/pto/mengyu/test/pristine-hse',
                        '/Volumes/Samsung/Data/Postdoc2/Data/Work/calculations/pto/mengyu/test/electron-from-neutral-todo',
                        '/Volumes/Samsung/Data/Postdoc2/Data/Work/calculations/pto/mengyu/test/electron',
                        '/Volumes/Samsung/Data/Postdoc2/Data/Work/calculations/pto/mengyu/test/vertical-from-electron'])
hse_values_polaron = [[25, 30, 33, 35]]

# HSE
data_hse_neutral = []
data_hse_vertical_from_neutral = []
data_hse_electron_offset = []
data_hse_vertical_from_electron = []
trapping_energy_hse = []
non_linearity_alpha_hse = []
ionisation_potential_hse_e = []
non_linearity_alpha_hse_e = []
non_linearity_beta_hse_e = []
ionisation_potential_hse_h = []
non_linearity_alpha_hse_h = []
non_linearity_beta_hse_h = []
non_linearity_homo_lumo_hse = []
for i in range(len(folder_use_hse)):
    data_hse_neutral.append(read_output(folder_files_hse[i][0], files[0]))
    data_hse_vertical_from_neutral.append(read_output(folder_files_hse[i][1], files[0]))
    data_hse_electron_offset.append(read_output(folder_files_hse[i][2], files[0]))
    data_hse_vertical_from_electron.append(read_output(folder_files_hse[i][3], files[0]))
    trapping_energy_hse.append(data_hse_vertical_from_neutral[i]['Energy_first'][:np.shape(hse_values_polaron[i])[0]] - data_hse_electron_offset[i]['Energy_last'][:np.shape(hse_values_polaron[i])[0]])
    # HOMO(N) = -IP = E(N) - E(N-1)
    ionisation_potential_hse_e.append((data_hse_electron_offset[i]['Energy_last'][:np.shape(hse_values_polaron[i])[0]] - data_hse_vertical_from_electron[i]['Energy_first'][:np.shape(hse_values_polaron[i])[0]]))
    non_linearity_alpha_hse_e.append(ionisation_potential_hse_e[i] - data_hse_electron_offset[i]['HOMO_alpha_last'][:np.shape(hse_values_polaron[i])[0]])
    non_linearity_beta_hse_e.append(ionisation_potential_hse_e[i] - data_hse_electron_offset[i]['HOMO_beta_last'][:np.shape(hse_values_polaron[i])[0]])
    # LUMO(N-1) = HOMO(N) on polaron geometry
    non_linearity_homo_lumo_hse.append(data_hse_vertical_from_electron[i]['LUMO_alpha_first'][:np.shape(hse_values_polaron[i])[0]] - data_hse_electron_offset[i]['HOMO_alpha_last'][:np.shape(hse_values_polaron[i])[0]])
    # non_linearity_homo_lumo_hse.append(data_hse_neutral[i]['LUMO_alpha_last'][:np.shape(hse_values_polaron[i])[0]] - data_hse_vertical_from_neutral[i]['HOMO_alpha_first'][:np.shape(hse_values_polaron[i])[0]])
    # LUMO(N-1) = -EA = E(N) - E(N-1)
    ionisation_potential_hse_h.append(data_hse_vertical_from_neutral[i]['Energy_first'][:np.shape(hse_values_polaron[i])[0]]-data_hse_neutral[i]['Energy_last'][:np.shape(hse_values_polaron[i])[0]])
    non_linearity_alpha_hse_h.append(ionisation_potential_hse_h[i] - data_hse_neutral[i]['LUMO_alpha_last'][:np.shape(hse_values_polaron[i])[0]])
    non_linearity_beta_hse_h.append(ionisation_potential_hse_h[i] - data_hse_neutral[i]['LUMO_beta_last'][:np.shape(hse_values_polaron[i])[0]])

# folder_save = folder_files_hse[0]
folder_save = []
for i in range(np.shape(folder_files_hse)[0]):
    folder_save.append(folder_files_hse[i][0])

# Plot band gap HFX
fig_band_gap2, ax_band_gap2 = plt.subplots()
for i in range(len(folder_use_hse)):
    ax_band_gap2.plot(hse_values_polaron[i], data_hse_neutral[i]['Gap_alpha'], 'x-', color=colors_hse[i], label='{}'.format(labels_hse[i]))
ax_band_gap2.set_ylabel('Band gap / eV')
ax_band_gap2.set_xlabel('% HFX')
# ax_band_gap2.set_ylim(set_ylim_band_gap_hfx[0], set_ylim_band_gap_hfx[1])
fig_band_gap2.tight_layout()
for i in range(np.shape(folder_save)[0]):
    fig_band_gap2.savefig('{}/band_gap_hfx.png'.format(folder_save[i]), dpi=300)

# Plot non-linearity HFX
fig_linearity2, ax_linearity2 = plt.subplots()
# ax_linearity2.axhline(y=draw_nonlinearity, color='k', alpha=0.5)
# ax_linearity2.axhline(y=-draw_nonlinearity, color='k', alpha=0.5)
# ax_linearity2.axhline(y=0.00, color='k', alpha=0.5, linestyle='--')
for i in range(len(folder_use_hse)):
    ax_linearity2.plot(hse_values_polaron[i],
                       non_linearity_homo_lumo_hse[i] * param.hartree_to_ev, 's-', fillstyle='none', color=colors_hse[i])
    ax_linearity2.plot(hse_values_polaron[i],
                       non_linearity_alpha_hse_e[i] * param.hartree_to_ev, 'o-', fillstyle='none', color=colors_hse[i], label='{}'.format(labels_hse[i]))
ax_linearity2.set_ylabel('Non-linearity / eV')
# ax_linearity2.set_ylim(set_ylim_linearity_hfx[0], set_ylim_linearity_hfx[1])
ax_linearity2.set_xlabel('% HFX')
fig_linearity2.tight_layout()
for i in range(np.shape(folder_save)[0]):
    fig_linearity2.savefig('{}/linearity_hfx.png'.format(folder_save[i]), dpi=300)

if __name__ == "__main__":
    print('Finished.')
    plt.show()

