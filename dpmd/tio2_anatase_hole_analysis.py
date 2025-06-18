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
folder_use_hubbard1 = 'pbe-u-ti'
folder_use_hubbard1 = 'pbe-u-o-ti-3'
# folder_use_hubbard1 = 'pbe-u-o'
folder_files_hubbard = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/anatase/archer/anatase/cell-441/geo-opt/mp-390/neutral/{}'.format(folder_use_hubbard1),
                        '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/anatase/archer/anatase/cell-441/geo-opt/mp-390/hole/from-neutral/{}'.format(folder_use_hubbard1),
                        '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/anatase/archer/anatase/cell-441/geo-opt/mp-390/hole/from-md-hole/{}'.format(folder_use_hubbard1),
                        '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/anatase/archer/anatase/cell-441/geo-opt/mp-390/neutral-from-hole/{}'.format(folder_use_hubbard1)]
hubbard_values_polaron = [0, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0]
skip_folders_hubbard = 1
data_hubbard_neutral = read_output(folder_files_hubbard[0], files[0])
data_hubbard_vertical_from_neutral = read_output(folder_files_hubbard[1], files[0])
data_hubbard_electron_offset = read_output(folder_files_hubbard[2], files[0])
data_hubbard_vertical_from_electron = read_output(folder_files_hubbard[3], files[0])
trapping_energy_hubbard = data_hubbard_vertical_from_neutral['Energy_first'][:np.shape(hubbard_values_polaron)[0]] - data_hubbard_electron_offset['Energy_last'][:np.shape(hubbard_values_polaron)[0]]
ionisation_potential_hubbard = -(data_hubbard_electron_offset['Energy_last'][:np.shape(hubbard_values_polaron)[0]] - data_hubbard_vertical_from_electron['Energy_first'][:np.shape(hubbard_values_polaron)[0]])  # Electron polaron IP
ionisation_potential_hubbard = -(data_hubbard_electron_offset['Energy_first'][:np.shape(hubbard_values_polaron)[0]] - data_hubbard_neutral['Energy_last'][:np.shape(hubbard_values_polaron)[0]])  # Hole polaron IP

# non_linearity_alpha_hubbard = -ionisation_potential_hubbard - data_hubbard_vertical_from_electron['HOMO_beta_first'][:np.shape(hubbard_values_polaron)[0]]  # E(polaron) - E(neutral) - E_polaron(HOMO_alpha)
# non_linearity_beta_hubbard = -ionisation_potential_hubbard - data_hubbard_vertical_from_electron['HOMO_beta_first'][:np.shape(hubbard_values_polaron)[0]]  # E(polaron) - E(neutral) - E_polaron(HOMO_beta)


# non_linearity_alpha_hubbard = ionisation_potential_hubbard - data_hubbard_electron_offset['HOMO_alpha_last'][:np.shape(hubbard_values_polaron)[0]]  # E(polaron) - E(neutral) - E_polaron(HOMO_alpha)
# non_linearity_beta_hubbard = ionisation_potential_hubbard - data_hubbard_electron_offset['HOMO_beta_last'][:np.shape(hubbard_values_polaron)[0]]  # E(polaron) - E(neutral) - E_polaron(HOMO_beta)
non_linearity_alpha_hubbard = ionisation_potential_hubbard - data_hubbard_vertical_from_electron['HOMO_alpha_first'][:np.shape(hubbard_values_polaron)[0]]  # E(polaron) - E(neutral) - E_polaron(HOMO_alpha)
non_linearity_beta_hubbard = ionisation_potential_hubbard - data_hubbard_vertical_from_electron['HOMO_beta_first'][:np.shape(hubbard_values_polaron)[0]]  # E(polaron) - E(neutral) - E_polaron(HOMO_beta)

print(ionisation_potential_hubbard)
# print(data_hubbard_electron_offset['HOMO_alpha_last'][:np.shape(hubbard_values_polaron)[0]])
print(data_hubbard_vertical_from_electron['HOMO_beta_first'][:np.shape(hubbard_values_polaron)[0]])

# HSE
# folder_use_hse = ['hfx-schwarz-e-1e-4-f-1e-6-cfit11-cpfit3', 'hfx-schwarz-e-1e-4-f-1e-6-fit11-cpfit3',
#                   'hfx-schwarz-e-1e-6-f-1e-6-fit9-cpfit3', 'hfx-schwarz-e-1e-6-f-1e-6-fit9-pfit3']
# folder_use_hse = folder_use_hse[3]
# folder_files_hse = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-336/geo-opt/convergence/neutral/neutral-opt/{}'.format(folder_use_hse),
#                     '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-336/geo-opt/convergence/electron/from-neutral-opt/{}'.format(folder_use_hse),
#                     '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-336/geo-opt/convergence/electron/from-u-4.0/{}'.format(folder_use_hse),
#                     '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-336/geo-opt/convergence/neutral/from-electron-opt/{}'.format(folder_use_hse)]
# hse_values_polaron = [12, 15, 18, 21, 22, 23, 24, 25, 30]
# skip_folders_hse_start = 4
# skip_folders_hse_end = 1
# data_hse_neutral = read_output(folder_files_hse[0], files[0])
# data_hse_vertical_from_neutral = read_output(folder_files_hse[1], files[0])
# data_hse_electron_offset = read_output(folder_files_hse[2], files[0])
# data_hse_vertical_from_electron = read_output(folder_files_hse[3], files[0])
# trapping_energy_hse = data_hse_vertical_from_neutral['Energy_first'][:np.shape(hse_values_polaron)[0]] - data_hse_electron_offset['Energy_last'][:np.shape(hse_values_polaron)[0]]
# ionisation_potential_hse = (data_hse_electron_offset['Energy_last'][:np.shape(hse_values_polaron)[0]] - data_hse_vertical_from_electron['Energy_first'][:np.shape(hse_values_polaron)[0]])
# non_linearity_alpha_hse = ionisation_potential_hse - data_hse_electron_offset['HOMO_alpha_last'][:np.shape(hse_values_polaron)[0]]  # E(polaron) - E(neutral) - E_polaron(HOMO_alpha)
# non_linearity_beta_hse = ionisation_potential_hse - data_hse_electron_offset['HOMO_beta_last'][:np.shape(hse_values_polaron)[0]]  # E(polaron) - E(neutral) - E_polaron(HOMO_beta)

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

folder_save = folder_files_hubbard

# Plot band gap Hubbard U
set_ylim = [2, 3.4]
draw_gap = 3.2
fig_band_gap1, ax_band_gap1 = plt.subplots()
ax_band_gap1.axhline(y=draw_gap, color='k', alpha=0.5)
ax_band_gap1.plot(hubbard_values_polaron, data_hubbard_neutral['Gap_alpha'], 'kx-', label='PBE+U')
ax_band_gap1.set_ylabel('Band gap / eV')
ax_band_gap1.set_xlabel('Hubbard U / eV')
ax_band_gap1.set_ylim(set_ylim[0], set_ylim[1])
fig_band_gap1.tight_layout()
for i in range(np.shape(folder_save)[0]):
    fig_band_gap1.savefig('{}/band_gap_hubbard.png'.format(folder_save[i]), dpi=300)

# Plot band gap HFX
# set_ylim = [2.4, 3.4]
# draw_legend = True
# fig_band_gap2, ax_band_gap2 = plt.subplots()
# ax_band_gap2.axhline(y=draw_gap, color='k', alpha=0.5)
# ax_band_gap2.plot(hse_values_polaron, data_hse_neutral['Gap_alpha'], 'rx-', label='HSE')
# ax_band_gap2.plot(pbe0_values_polaron, data_pbe0_neutral['Gap_alpha'], 'gx-', label='PBE0')
# ax_band_gap2.set_ylabel('Band gap / eV')
# ax_band_gap2.set_xlabel('% HFX')
# ax_band_gap2.set_ylim(set_ylim[0], set_ylim[1])
# fig_band_gap2.tight_layout()
# if draw_legend: ax_band_gap2.legend(frameon=False)
# for i in range(np.shape(folder_save)[0]):
#     fig_band_gap2.savefig('{}/band_gap_hfx.png'.format(folder_save[i]), dpi=300)

# Plot trapping energy Hubbard U
fig_trapping1, ax_trapping1 = plt.subplots()
ax_trapping1.plot(hubbard_values_polaron[skip_folders_hubbard:], trapping_energy_hubbard[skip_folders_hubbard:]*param.hartree_to_ev*1e3, 'kx-')
ax_trapping1.set_ylabel('Trapping energy / meV')
ax_trapping1.set_xlabel('Hubbard U / eV')
fig_trapping1.tight_layout()
for i in range(np.shape(folder_save)[0]):
    fig_trapping1.savefig('{}/trapping_hubbard.png'.format(folder_save[i]), dpi=300)

# Plot trapping energy HFX
# set_ylim = [0, 160]
# fig_trapping2, ax_trapping2 = plt.subplots()
# ax_trapping2.plot(hse_values_polaron[skip_folders_hse_start:-skip_folders_hse_end], trapping_energy_hse[skip_folders_hse_start:-skip_folders_hse_end] * param.hartree_to_ev*1e3, 'rx-', label='HSE')
# ax_trapping2.plot(pbe0_values_polaron[skip_folders_pbe0:], trapping_energy_pbe0[skip_folders_pbe0:] * param.hartree_to_ev*1e3, 'gx-', label='PBE0')
# ax_trapping2.set_ylabel('Trapping energy / meV')
# ax_trapping2.set_xlabel('% HFX')
# ax_trapping2.set_ylim(set_ylim[0], set_ylim[1])
# if draw_legend: ax_trapping2.legend(frameon=False)
# fig_trapping2.tight_layout()
# for i in range(np.shape(folder_save)[0]):
#     fig_trapping2.savefig('{}/trapping_hfx.png'.format(folder_save[i]), dpi=300)

# Plot Hubbard U
set_ylim = [-0.1, 1]
fig_linearity1, ax_linearity1 = plt.subplots()
ax_linearity1.plot(hubbard_values_polaron[skip_folders_hubbard:], non_linearity_alpha_hubbard[skip_folders_hubbard:] * param.hartree_to_ev, 'kx-', label='PBE+U')
# ax_linearity1.plot(hubbard_values_polaron[skip_folders_hubbard:], non_linearity_beta_hubbard[skip_folders_hubbard:] * param.hartree_to_ev, 'kx-', label='PBE+U')
ax_linearity1.set_ylabel('Non-linearity / eV')
ax_linearity1.set_xlabel('Hubbard U / eV')
ax_linearity1.set_ylim(set_ylim[0], set_ylim[1])
fig_linearity1.tight_layout()
for i in range(np.shape(folder_save)[0]):
    fig_linearity1.savefig('{}/linearity_hubbard.png'.format(folder_save[i]), dpi=300)

# Plot non-linearity HFX
# set_ylim = [0, 0.16]
# fig_linearity2, ax_linearity2 = plt.subplots()
# ax_linearity2.plot(hse_values_polaron[skip_folders_hse_start:-skip_folders_hse_end], non_linearity_alpha_hse[skip_folders_hse_start:-skip_folders_hse_end] * param.hartree_to_ev, 'rx-', label='HSE')
# ax_linearity2.plot(pbe0_values_polaron[skip_folders_pbe0:], non_linearity_alpha_pbe0[skip_folders_pbe0:] * param.hartree_to_ev, 'gx-', label='PBE0')
# ax_linearity2.set_ylabel('Non-linearity / eV')
# ax_linearity2.set_ylim(set_ylim[0], set_ylim[1])
# ax_linearity2.set_xlabel('% HFX')
# if draw_legend: ax_linearity2.legend(frameon=False)
# fig_linearity2.tight_layout()
# for i in range(np.shape(folder_save)[0]):
#     fig_linearity2.savefig('{}/linearity_hfx.png'.format(folder_save[i]), dpi=300)

if __name__ == "__main__":
    print('Finished.')
    plt.show()

