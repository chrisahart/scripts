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


# Hubbard U
# folder_files = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-336/geo-opt/md-struct/electron-opt',
#                 '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-336/geo-opt/md-struct/neutral-from-electron-opt',
#                 '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-336/geo-opt/md-struct/electron-from-neutral-from-electron-opt']
# folder_save = folder_files[1]
# files = ['output.out'] * 5
# num_atoms = 324
# u_values_polaron = [0, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0]
# u_values_neutral = [0, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 8.0, 9.0, 10.0]
# box_size = [13.8, 13.8, 17.76, 90, 90, 90]
# skip_folders = 1
# data_electron_md = read_output(folder_files[0], files[0])
# data_neutral = read_output(folder_files[1], files[0])
# data_electron = read_output(folder_files[2], files[0])

# HSE
folder_files = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-336/geo-opt/convergence/neutral/neutral-opt/hfx-schwarz-e-1e-4-f-1e-6-fit11-cpfit3',
                '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-336/geo-opt/convergence/neutral/neutral-opt/hfx-schwarz-e-1e-6-f-1e-6-fit9-cpfit3',
                '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-336/geo-opt/convergence/neutral/neutral-opt/hfx-schwarz-e-1e-6-f-1e-6-fit9-pfit3',
                '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-336/geo-opt/convergence/electron/from-neutral-opt/hfx-schwarz-e-1e-4-f-1e-6-fit11-cpfit3',
                '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-336/geo-opt/convergence/electron/from-neutral-opt/hfx-schwarz-e-1e-6-f-1e-6-fit9-cpfit3',
                '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-336/geo-opt/convergence/electron/from-neutral-opt/hfx-schwarz-e-1e-6-f-1e-6-fit9-pfit3'
                '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-336/geo-opt/convergence/electron/from-u-4.0/hfx-schwarz-e-1e-4-f-1e-6-fit11-cpfit3',
                '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-336/geo-opt/convergence/electron/from-u-4.0/hfx-schwarz-e-1e-6-f-1e-6-fit9-cpfit3',
                '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-336/geo-opt/convergence/electron/from-u-4.0/hfx-schwarz-e-1e-6-f-1e-6-fit9-pfit3'
                ]
folder_save = folder_files[1]
files = ['output.out'] * 5
hfx_values = [12, 15, 18, 22, 25]
skip_folders = 0

num_atoms = 324
box_size = [13.8, 13.8, 17.76, 90, 90, 90]

test = 0
data_neutral = read_output(folder_files[test+3*0], files[0])
data_electron_vertical = read_output(folder_files[test+3*1], files[0])
data_electron_offset = read_output(folder_files[test+3*2], files[0])

band_gap_alpha = (data_neutral['LUMO_alpha_last'] - data_neutral['HOMO_alpha_last']) * param.hartree_to_ev
band_gap_beta = (data_neutral['LUMO_beta_last'] - data_neutral['HOMO_beta_last']) * param.hartree_to_ev

electron_affinity = data_electron_vertical['Energy_first'][:np.shape(hfx_values)[0]] - data_neutral['Energy_last'][:np.shape(hfx_values)[0]]
# print(electron_affinity)
# print(data_neutral['HOMO_alpha_last'][:np.shape(u_values_polaron)[0]])

# non-linearity as E(N+1) - E(N) - HOMO(N+1), meaning E(polaron) - E(neutral) - E_polaron(HOMO)
non_linearity = electron_affinity - data_neutral['HOMO_alpha_last'][:np.shape(hfx_values)[0]]
# print(non_linearity)

trapping_energy = data_electron_vertical['Energy_first'][:np.shape(hfx_values)[0]] - data_electron_offset['Energy_last'][:np.shape(hfx_values)[0]]
# trapping_energy = data_electron['Energy_first'][:np.shape(u_values_polaron)[0]] - data_electron['Energy_last'][:np.shape(u_values_polaron)[0]]
print(data_electron_vertical['Energy_first'][:np.shape(hfx_values)[0]])
print(data_electron_offset['Energy_last'][:np.shape(hfx_values)[0]])
print(trapping_energy)

# emission = data_electron_md['Energy_last'][:np.shape(u_values_polaron)[0]] - data_neutral['Energy_first'][:np.shape(u_values_polaron)[0]]
# print(emission*param.hartree_to_ev)

# Plot band gap
fig_band_gap1, ax_band_gap1 = plt.subplots()
ax_band_gap1.plot(hfx_values, data_neutral['Gap_alpha'], 'kx-')
ax_band_gap1.plot(hfx_values, data_neutral['Gap_beta'], 'kx-')
ax_band_gap1.plot(hfx_values, band_gap_alpha, 'kx-')
ax_band_gap1.plot(hfx_values, band_gap_beta, 'kx-')
ax_band_gap1.set_ylabel('Band gap / eV')
# ax_band_gap1.set_xlabel('Hubbard U / eV')
ax_band_gap1.set_xlabel('% HFX')
fig_band_gap1.tight_layout()
fig_band_gap1.savefig('{}/band_gap.png'.format(folder_save), dpi=300)

# Plot electron affinity
# fig_electron_affinity, ax_electron_affinity = plt.subplots()
# ax_electron_affinity.plot(u_values_polaron, electron_affinity*param.hartree_to_ev, 'kx-')
# ax_electron_affinity.set_ylabel('Electron affinity / eV')
# ax_electron_affinity.set_xlabel('Hubbard U / eV')
# fig_electron_affinity.tight_layout()
# fig_electron_affinity.savefig('{}/electron_affinity.png'.format(folder_save), dpi=300)

# Plot trapping energy
fig_trapping, ax_trapping = plt.subplots()
ax_trapping.plot(hfx_values[skip_folders:], trapping_energy[skip_folders:]*param.hartree_to_ev, 'kx-')
ax_trapping.set_ylabel('Trapping energy / eV')
# ax_trapping.set_xlabel('Hubbard U / eV')
ax_trapping.set_xlabel('% HFX')
fig_trapping.tight_layout()
fig_trapping.savefig('{}/trapping.png'.format(folder_save), dpi=300)

# Plot non-linearity
# fig_nonlinearity, ax_nonlinearity = plt.subplots()
# ax_nonlinearity.plot(u_values_polaron[1:], non_linearity[1:]*param.hartree_to_ev, 'kx-')
# ax_nonlinearity.set_ylabel('Non-linearity / eV')
# ax_nonlinearity.set_xlabel('Hubbard U / eV')
# fig_nonlinearity.tight_layout()
# fig_nonlinearity.savefig('{}/nonlinearity.png'.format(folder_save), dpi=300)

if __name__ == "__main__":
    print('Finished.')
    plt.show()

