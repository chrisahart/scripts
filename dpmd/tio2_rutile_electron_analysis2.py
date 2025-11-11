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

# # cell 3x3x6 OLD cell vectors a = b = 4.60, c = 2.96
# cell vectors from materials project https://next-gen.materialsproject.org/materials/mp-2657/
# num_atoms = 324
# box_size = [13.8, 13.8, 17.76, 90, 90, 90]
# colors_hubbard = ['k']
# labels_hubbard = ['']
# skip_folders_hubbard_start = 4  # starting at 3 eV
# # skip_folders_hubbard_start = 7  # starting at 3 eV
# # skip_folders_hubbard_end = NaN
# folder_use_hubbard = ['pbe-u-ti']
# # folder_use_hubbard = []
# folder_files_hubbard = []
# for i in range(len(folder_use_hubbard)):
#     folder_files_hubbard.append(['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-336/geo-opt/md-struct/neutral-from-electron-opt',
#                         '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-336/geo-opt/md-struct/electron-from-neutral-from-electron-opt',
#                         '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-336/geo-opt/md-struct/electron-opt',
#                         '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-336/geo-opt/md-struct/neutral-from-electron-opt'])
# hubbard_values_polaron = [[0, 2.0, 2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 8.0, 9.0]]

# # cell 2x2x3 HSE only (PBE0 gives wierd results, likely due to rcut = 6 > cell/2)
# num_atoms = 72
# box_size = [9.18, 9.18, 8.88, 90, 90, 90]
# folder_use_hubbard = []
# colors_hse = ['k']
# labels_hse = ['']
# skip_folders_hse_start = 3
# skip_folders_hse_end = -1
# skip_folders_hse_end_gap = -1
# folder_use_hse = ['hfx-schwarz-e-1e-6-f-1e-6-fit9-pfit3']
# folder_files_hse = []
# for i in range(len(folder_use_hse)):
#     folder_files_hse.append(['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-223/geo-opt-cell-opt/neutral/neutral-opt/{}'.format(folder_use_hse[i]),
#                         '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-223/geo-opt-cell-opt/electron/from-neutral/{}'.format(folder_use_hse[i]),
#                         '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-223/geo-opt-cell-opt/electron/from-md-1/{}'.format(folder_use_hse[i]),
#                         '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-223/geo-opt-cell-opt/neutral/from-electron//{}'.format(folder_use_hse[i])])
# hse_values_polaron = [[12, 15, 18, 20, 21, 22, 23, 24, 25, 30]]
# folder_use_pbe0 = []

# # cell 3x3x4 HSE only (PBE0 gives wierd results, likely due to rcut = 6 > cell/2)
# num_atoms = 216
# box_size = [13.77, 13.77, 11.84, 90, 90, 90]
# colors_hubbard = ['k']
# labels_hubbard = ['']
# skip_folders_hubbard_start = 0
# skip_folders_hubbard_end = -6
# folder_use_hubbard = []
# colors_hse = ['k']
# labels_hse = ['']
# skip_folders_hse_start = 3
# skip_folders_hse_end = -1
# skip_folders_hse_end_gap = -1
# folder_use_hse = ['hfx-schwarz-e-1e-6-f-1e-6-fit9-pfit3']
# folder_files_hse = []
# for i in range(len(folder_use_hse)):
#     folder_files_hse.append(['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-334/geo-opt-cell-opt/neutral/neutral-opt/{}'.format(folder_use_hse[i]),
#                         '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-334/geo-opt-cell-opt/electron/from-neutral/{}'.format(folder_use_hse[i]),
#                         '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-334/geo-opt-cell-opt/electron/from-md-1/{}'.format(folder_use_hse[i]),
#                         '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-334/geo-opt-cell-opt/neutral/from-electron//{}'.format(folder_use_hse[i])])
# hse_values_polaron = [[12, 15, 18, 20, 21, 22, 23, 24, 25, 30]]
# folder_use_pbe0 = []


# cell 3x3x5
# num_atoms = 270
# box_size = [13.77, 13.77, 14.80, 90, 90, 90]
# colors_hubbard = ['k']
# labels_hubbard = ['']
# skip_folders_hubbard_start = 0
# skip_folders_hubbard_end = -6
# folder_use_hubbard = ['pbe-u-ti']
# folder_use_hubbard = []
# folder_files_hubbard = []
# for i in range(len(folder_use_hubbard)):
#     folder_files_hubbard.append(['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-335/geo-opt-cell-opt/neutral/neutral-opt/{}'.format(folder_use_hubbard[i]),
#                         '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-335/geo-opt-cell-opt/electron/from-neutral/{}'.format(folder_use_hubbard[i]),
#                         '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-335/geo-opt-cell-opt/electron/from-md-1/{}'.format(folder_use_hubbard[i]),
#                         '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-335/geo-opt-cell-opt/neutral/from-electron//{}'.format(folder_use_hubbard[i])])
# hubbard_values_polaron = [[2.5, 2.6, 2.7, 2.8, 2.9, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 7.0, 8.0, 9.0]]
# colors_hse = ['k']
# labels_hse = ['']
# skip_folders_hse_start = 3
# skip_folders_hse_end = -1
# skip_folders_hse_end_gap = -1
# folder_use_hse = ['hfx-schwarz-e-1e-6-f-1e-6-fit9-pfit3']
# folder_files_hse = []
# for i in range(len(folder_use_hse)):
#     folder_files_hse.append(['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-335/geo-opt-cell-opt/neutral/neutral-opt/{}'.format(folder_use_hse[i]),
#                         '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-335/geo-opt-cell-opt/electron/from-neutral/{}'.format(folder_use_hse[i]),
#                         '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-335/geo-opt-cell-opt/electron/from-md-1/{}'.format(folder_use_hse[i]),
#                         '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-335/geo-opt-cell-opt/neutral/from-electron//{}'.format(folder_use_hse[i])])
# hse_values_polaron = [[12, 15, 18, 20, 21, 22, 23, 24, 25, 30]]
# folder_use_pbe0 = ['pbe0-tz-tz2p-hfx-schwarz-e-1e-6-f-1e-6-fit9-pfit3-mckenna-input']
# folder_use_pbe0 = []
# folder_files_pbe0 = []
# for i in range(len(folder_use_pbe0)):
#     folder_files_pbe0.append(['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-335/geo-opt-cell-opt/neutral/neutral-opt/{}'.format(folder_use_pbe0[i]),
#                         '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-335/geo-opt-cell-opt/electron/from-neutral/{}'.format(folder_use_pbe0[i]),
#                         '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-335/geo-opt-cell-opt/electron/from-md-1/{}'.format(folder_use_pbe0[i]),
#                         '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-335/geo-opt-cell-opt/neutral/from-electron//{}'.format(folder_use_pbe0[i])])
# pbe0_values_polaron = [[10, 10.5, 11, 11.5, 12, 13, 14, 15, 25]]
# colors_pbe0 = ['g']
# labels_pbe0 = ['PBE0']
# skip_folders_pbe0_start = 0
# skip_folders_pbe0_end = -1
# skip_folders_pbe0_end_gap = -1

# # cell 3x3x6 HSE only
num_atoms = 324
box_size = [13.77, 13.77, 17.76, 90, 90, 90]
colors_hse = ['k']
labels_hse = ['']
skip_folders_hse_start = 3
skip_folders_hse_end = -1
skip_folders_hse_start_gap = 3
skip_folders_hse_end_gap = -1
folder_use_hse = ['hfx-schwarz-e-1e-6-f-1e-6-fit9-pfit3']
folder_files_hse = []
for i in range(len(folder_use_hse)):
    folder_files_hse.append(['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-336/geo-opt-cell-opt/neutral/neutral-opt/{}'.format(folder_use_hse[i]),
                        '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-336/geo-opt-cell-opt/electron/from-neutral/hfx-schwarz-e-1e-6-f-1e-6-fit9-pfit3-rs-scf-1e-7',
                        '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-336/geo-opt-cell-opt/electron/from-pbe0-14/hfx-schwarz-e-1e-6-f-1e-6-fit9-pfit3-tight',
                        '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-336/geo-opt-cell-opt/neutral/from-electron-pbe0-14/hfx-schwarz-e-1e-6-f-1e-6-fit9-pfit3-tight'])
hse_values_polaron = [[12, 15, 18, 20, 21, 22, 23, 24, 25, 30]]

# # cell 4x4x6
# num_atoms = 576
# box_size = [18.36, 18.36, 17.76, 90, 90, 90]
# colors_hubbard = ['k']
# labels_hubbard = ['']
# skip_folders_hubbard_start = 0
# skip_folders_hubbard_end = -6
# folder_use_hubbard = ['pbe-u-ti']
# folder_use_hubbard = []
# folder_files_hubbard = []
# colors_hse = ['k']
# labels_hse = ['']
# skip_folders_hse_start = 2
# skip_folders_hse_end = -1
# skip_folders_hse_end_gap = -1
# folder_use_hse = ['hfx-schwarz-e-1e-6-f-1e-6-fit9-pfit3']
# folder_use_hse = ['hfx-schwarz-e-1e-6-f-1e-6-fit9-pfit3-rs-tz-derived-ti-tz2p-o']  # FAILED
# # folder_use_hse = ['hfx-schwarz-e-1e-6-f-1e-6-fit9-pfit3-rs-tz-ti']
# # folder_use_hse = ['hfx-schwarz-e-1e-6-f-1e-6-fit9-pfit3-rs-tz-ti-o']
# # folder_use_hse = ['hfx-schwarz-e-1e-6-f-1e-6-fit9-pfit3-rs-tz-ti-tz2p-o']  # FAILED
# folder_files_hse = []
# # for i in range(len(folder_use_hse)):
# #     folder_files_hse.append(['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-446/geo-opt-cell-opt/neutral/neutral-opt/{}'.format(folder_use_hse[i]),
# #                         '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-446/geo-opt-cell-opt/electron/from-neutral/{}'.format(folder_use_hse[i]),
# #                         '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-446/geo-opt-cell-opt/electron/from-offset/{}'.format(folder_use_hse[i]),
# #                         '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-446/geo-opt-cell-opt/neutral/from-electron-offset/{}'.format(folder_use_hse[i])])
# for i in range(len(folder_use_hse)):
#     folder_files_hse.append(['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-446/geo-opt-cell-opt/neutral/neutral-opt/{}'.format(folder_use_hse[i]),
#                         '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-446/geo-opt-cell-opt/electron/from-neutral/{}'.format(folder_use_hse[i]),
#                         '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-446/geo-opt-cell-opt/electron/from-pbe0-14/{}'.format(folder_use_hse[i]),
#                         '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-446/geo-opt-cell-opt/neutral/from-electron-pbe0-14/{}'.format(folder_use_hse[i])])
# hse_values_polaron = [[12, 15, 18, 20, 21, 22, 23, 24, 25, 30]]
# folder_use_pbe0 = ['pbe0-tz-tz2p-hfx-schwarz-e-1e-6-f-1e-6-fit9-pfit3-mckenna-input']
# folder_files_pbe0 = []
# # for i in range(len(folder_use_pbe0)):
# #     folder_files_pbe0.append(['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-446/geo-opt-cell-opt/neutral/neutral-opt/{}'.format(folder_use_pbe0[i]),
# #                         '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-446/geo-opt-cell-opt/electron/from-neutral/{}'.format(folder_use_pbe0[i]),
# #                         '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-446/geo-opt-cell-opt/electron/from-offset/{}'.format(folder_use_pbe0[i]),
# #                         '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-446/geo-opt-cell-opt/neutral/from-electron-offset/{}'.format(folder_use_pbe0[i])])
# for i in range(len(folder_use_pbe0)):
#     folder_files_pbe0.append(['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-446/geo-opt-cell-opt/neutral/neutral-opt/{}'.format(folder_use_pbe0[i]),
#                         '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-446/geo-opt-cell-opt/electron/from-neutral/{}'.format(folder_use_pbe0[i]),
#                         '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-446/geo-opt-cell-opt/electron/from-pbe0-14/{}'.format(folder_use_pbe0[i]),
#                         '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-446/geo-opt-cell-opt/neutral/from-electron-pbe0-14/{}'.format(folder_use_pbe0[i])])
# pbe0_values_polaron = [[10, 10.5, 11, 11.5, 12, 13, 14, 15, 25]]
# colors_pbe0 = ['g']
# labels_pbe0 = ['PBE0']
# skip_folders_pbe0_start = 3
# skip_folders_pbe0_end = -1
# skip_folders_pbe0_end_gap = -1

# cell 4x4x6 basis set
# num_atoms = 576
# box_size = [18.36, 18.36, 17.76, 90, 90, 90]
# folder_use_hubbard = []
# folder_files_hubbard = []
# colors_hse = ['r', 'g', 'b', 'm']
# labels_hse = ['DZ, O DZP', 'Ti TZP, O DZP', 'Ti TZP, O TZP', 'Ti TZP, O TZ2P']
# labels_hse = ['DZ', 'TZ']
# # labels_hse = ['DZ', 'DZ -> TZ', 'TZ']
# skip_folders_hse_start = 3
# skip_folders_hse_end = -1
# skip_folders_hse_end_gap = -1
# folder_use_hse = ['hfx-schwarz-e-1e-6-f-1e-6-fit9-pfit3', 'hfx-schwarz-e-1e-6-f-1e-6-fit9-pfit3-rs-tz-ti-tz2p-o']
# folder_files_hse = []
# # for i in range(len(folder_use_hse)):
# #     folder_files_hse.append(['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-446/geo-opt-cell-opt/neutral/neutral-opt/{}'.format(folder_use_hse[i]),
# #                         '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-446/geo-opt-cell-opt/electron/from-neutral/{}'.format(folder_use_hse[i]),
# #                         '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-446/geo-opt-cell-opt/electron/from-offset/{}'.format(folder_use_hse[i]),
# #                         '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-446/geo-opt-cell-opt/neutral/from-electron-offset/{}'.format(folder_use_hse[i])])
# for i in range(len(folder_use_hse)):
#     folder_files_hse.append(['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-446/geo-opt-cell-opt/neutral/neutral-opt/{}'.format(folder_use_hse[i]),
#                         '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-446/geo-opt-cell-opt/electron/from-neutral/{}-rs-scf-1e-7'.format(folder_use_hse[i]),
#                         '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-446/geo-opt-cell-opt/electron/from-pbe0-14/{}'.format(folder_use_hse[i]),
#                         '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-446/geo-opt-cell-opt/neutral/from-electron-pbe0-14/{}'.format(folder_use_hse[i])])
# # for i in range(len(folder_use_hse)):
# #     folder_files_hse.append(['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-446/geo-opt-cell-opt/neutral/neutral-opt/{}'.format(folder_use_hse[i]),
# #                         '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-446/geo-opt-cell-opt/electron/from-neutral/hfx-schwarz-e-1e-6-f-1e-6-fit9-pfit3-rs-scf-1e-7',
# #                         '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-446/geo-opt-cell-opt/electron/from-pbe0-14/{}'.format(folder_use_hse[i]),
# #                         '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-446/geo-opt-cell-opt/neutral/from-electron-pbe0-14/{}'.format(folder_use_hse[i])])
# hse_values_polaron = [[12, 15, 18, 20, 21, 22, 23, 24, 25, 30], [12, 15, 18, 20, 21, 22, 23, 24, 25, 30], [12, 15, 18, 20, 21, 22, 23, 24, 25, 30], [12, 15, 18, 20, 21, 22, 23, 24, 25, 30]]
# folder_use_pbe0 = []
# folder_files_pbe0 = []

# cell 3x3x8
# colors_hubbard = ['k']
# labels_hubbard = ['']
# skip_folders_hubbard_start = 0
# skip_folders_hubbard_end = -6
# folder_use_hubbard = []
# folder_files_hubbard = []
# colors_hse = ['k']
# labels_hse = ['']
# skip_folders_hse_start = 2
# skip_folders_hse_end = -1
# skip_folders_hse_end_gap = -1
# folder_use_hse = ['hfx-schwarz-e-1e-6-f-1e-6-fit9-pfit3']
# folder_files_hse = []
# for i in range(len(folder_use_hse)):
#     folder_files_hse.append(['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-338/geo-opt-cell-opt/neutral/neutral-opt/{}'.format(folder_use_hse[i]),
#                         '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-338/geo-opt-cell-opt/electron/from-neutral/{}'.format(folder_use_hse[i]),
#                         '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-338/geo-opt-cell-opt/electron/from-offset-107-hse-25/{}'.format(folder_use_hse[i]),
#                         '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-338/geo-opt-cell-opt/neutral/from-offset-107-hse-25/{}'.format(folder_use_hse[i])])
# hse_values_polaron = [[12, 15, 18, 20, 21, 22, 23, 24, 25, 30]]
# folder_use_pbe0 = []
# folder_files_pbe0 = []

# cell 3x3x10
# colors_hubbard = ['k']
# labels_hubbard = ['']
# skip_folders_hubbard_start = 0
# skip_folders_hubbard_end = -6
# folder_use_hubbard = []
# folder_files_hubbard = []
# colors_hse = ['k']
# labels_hse = ['']
# skip_folders_hse_start = 2
# skip_folders_hse_end = -1
# skip_folders_hse_end_gap = -1
# folder_use_hse = ['hfx-schwarz-e-1e-6-f-1e-6-fit9-pfit3']
# folder_files_hse = []
# for i in range(len(folder_use_hse)):
#     folder_files_hse.append(['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-3310/geo-opt-cell-opt/neutral/neutral-opt/{}'.format(folder_use_hse[i]),
#                         '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-3310/geo-opt-cell-opt/electron/from-neutral/{}'.format(folder_use_hse[i]),
#                         '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-3310/geo-opt-cell-opt/electron/from-offset-135-hse-25/{}'.format(folder_use_hse[i]),
#                         '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-3310/geo-opt-cell-opt/neutral/from-offset-135-hse-25/{}'.format(folder_use_hse[i])])
# hse_values_polaron = [[12, 15, 18, 20, 21, 22, 23, 24, 25, 30]]
# folder_use_pbe0 = []
# folder_files_pbe0 = []

# cell 3x3x12
# colors_hubbard = ['k']
# labels_hubbard = ['']
# skip_folders_hubbard_start = 0
# skip_folders_hubbard_end = -6
# folder_use_hubbard = []
# folder_files_hubbard = []
# colors_hse = ['k']
# labels_hse = ['']
# skip_folders_hse_start = 2
# skip_folders_hse_end = -1
# skip_folders_hse_end_gap = -1
# folder_use_hse = ['hfx-schwarz-e-1e-6-f-1e-6-fit9-pfit3']
# folder_files_hse = []
# for i in range(len(folder_use_hse)):
#     folder_files_hse.append(['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-3312/geo-opt-cell-opt/neutral/neutral-opt/{}'.format(folder_use_hse[i]),
#                         '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-3312/geo-opt-cell-opt/electron/from-neutral/{}'.format(folder_use_hse[i]),
#                         '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-3312/geo-opt-cell-opt/electron/from-offset-163-hse-25/{}'.format(folder_use_hse[i]),
#                         '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-3312/geo-opt-cell-opt/neutral/from-offset-163-hse-25/{}'.format(folder_use_hse[i])])
# hse_values_polaron = [[12, 15, 18, 20, 21, 22, 23, 24, 25, 30]]
# folder_use_pbe0 = []
# folder_files_pbe0 = []

# cell 3x3x6 4x4x6 HSE PBE0
# set_ylim_linearity_hfx = [-0.07, 0.20]
# set_ylim_band_gap_hfx = [2.8, 3.3]
# set_ylim_trapping_hfx = [0, 210]
#
# folder_use_hubbard = []
# colors_hse = ['r', 'g']
# labels_hse = ['HSE 3x3x6 (324 atoms)', 'HSE 4x4x6 (576 atoms)']
# colors_pbe0 = ['b', 'm']
# labels_pbe0 = ['PBE0 3x3x6 (324 atoms)', 'PBE0 4x4x6 (576 atoms)']
# skip_folders_hse_start = 3
# skip_folders_hse_end = -1
# skip_folders_hse_start_gap = 3
# skip_folders_hse_end_gap = -1
# skip_folders_pbe0_start_gap = 3
# skip_folders_pbe0_end_gap = -1
# skip_folders_pbe0_start = 3
# skip_folders_pbe0_end = -1
# skip_folders_pbe0_end_gap = -1
# folder_files_hse = []
# folder_use_pbe0 = []
# folder_use_hse = ['hfx-schwarz-e-1e-6-f-1e-6-fit9-pfit3']
# # folder_use_hse = ['hfx-schwarz-e-1e-6-f-1e-6-fit9-pfit3-tight']
# # for i in range(len(folder_use_hse)):
# #     folder_files_hse.append(['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-336/geo-opt-cell-opt/neutral/neutral-opt/{}'.format(folder_use_hse[i]),
# #                         '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-336/geo-opt-cell-opt/electron/from-neutral/hfx-schwarz-e-1e-6-f-1e-6-fit9-pfit3-rs-scf-1e-7',
# #                         '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-336/geo-opt-cell-opt/electron/from-pbe0-14/hfx-schwarz-e-1e-6-f-1e-6-fit9-pfit3',
# #                         '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-336/geo-opt-cell-opt/neutral/from-electron-pbe0-14/hfx-schwarz-e-1e-6-f-1e-6-fit9-pfit3'])
# for i in range(len(folder_use_hse)):
#     folder_files_hse.append(['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-336/geo-opt-cell-opt/neutral/neutral-opt/{}'.format(folder_use_hse[i]),
#                         '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-336/geo-opt-cell-opt/electron/from-neutral/hfx-schwarz-e-1e-6-f-1e-6-fit9-pfit3-rs-scf-1e-7',
#                         '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-336/geo-opt-cell-opt/electron/from-pbe0-14/hfx-schwarz-e-1e-6-f-1e-6-fit9-pfit3-tight',
#                         '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-336/geo-opt-cell-opt/neutral/from-electron-pbe0-14/hfx-schwarz-e-1e-6-f-1e-6-fit9-pfit3-tight'])
# folder_use_hse = ['hfx-schwarz-e-1e-6-f-1e-6-fit9-pfit3']
# # folder_use_hse = ['hfx-schwarz-e-1e-6-f-1e-6-fit9-pfit3-tight']
# # folder_use_hse = ['hfx-schwarz-e-1e-6-f-1e-6-fit9-pfit3-rs-tz-ti-o']
# for i in range(len(folder_use_hse)):
#     folder_files_hse.append(['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-446/geo-opt-cell-opt/neutral/neutral-opt/{}'.format(folder_use_hse[i]),
#                         '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-446/geo-opt-cell-opt/electron/from-neutral/hfx-schwarz-e-1e-6-f-1e-6-fit9-pfit3-rs-scf-1e-7',
#                         '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-446/geo-opt-cell-opt/electron/from-pbe0-14/{}'.format(folder_use_hse[i]),
#                         '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-446/geo-opt-cell-opt/neutral/from-electron-pbe0-14/{}'.format(folder_use_hse[i])])
# hse_values_polaron = [[12, 15, 18, 20, 21, 22, 23, 24, 25, 30], [12, 15, 18, 20, 21, 22, 23, 24, 25, 30]]
# folder_use_hse = ['hfx-schwarz-e-1e-6-f-1e-6-fit9-pfit3'] * 2
# folder_files_pbe0 = []
# folder_use_pbe0 = ['pbe0-tz-tz2p-hfx-schwarz-e-1e-6-f-1e-6-fit9-pfit3-mckenna-input']
# for i in range(len(folder_use_pbe0)):
#     folder_files_pbe0.append(['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-336/geo-opt-cell-opt/neutral/neutral-opt/{}'.format(folder_use_pbe0[i]),
#                         '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-336/geo-opt-cell-opt/electron/from-neutral/{}'.format(folder_use_pbe0[i]),
#                         '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-336/geo-opt-cell-opt/electron/from-pbe0-14/{}'.format(folder_use_pbe0[i]),
#                         '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-336/geo-opt-cell-opt/neutral/from-electron-pbe0-14/{}'.format(folder_use_pbe0[i])])
# for i in range(len(folder_use_pbe0)):
#     folder_files_pbe0.append(['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-446/geo-opt-cell-opt/neutral/neutral-opt/{}'.format(folder_use_pbe0[i]),
#                         '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-446/geo-opt-cell-opt/electron/from-neutral/{}'.format(folder_use_pbe0[i]),
#                         '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-446/geo-opt-cell-opt/electron/from-pbe0-14/{}'.format(folder_use_pbe0[i]),
#                         '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-446/geo-opt-cell-opt/neutral/from-electron-pbe0-14/{}'.format(folder_use_pbe0[i])])
# pbe0_values_polaron = [[10, 10.5, 11, 11.5, 12, 13, 14, 15, 25], [10, 10.5, 11, 11.5, 12, 13, 14, 15, 25]]
# folder_use_pbe0 = ['pbe0-tz-tz2p-hfx-schwarz-e-1e-6-f-1e-6-fit9-pfit3-mckenna-input'] * 2


# # # cell all
# folder_use_hubbard = []
# colors_hse = ['r', 'g', 'b', 'm', 'y', 'orange', 'k']
# labels_hse = ['334', '335', '336', '446']
# # labels_hse = ['3x3x4', '3x3x5', '3x3x6', '4x4x6', '3x3x8', '3x3x10', '3x3x12']
# # labels_hse = ['3x3x4 (216 atoms)', '3x3x5 (270 atoms)', '3x3x6 (324 atoms)', '4x4x6 (576 atoms)', '3x3x8 (432 atoms)', '3x3x10 (540 atoms)', '3x3x12 (648 atoms)']
# # labels_hse = ['3x3x4 (216 atoms)', '3x3x6 (324 atoms)', '4x4x6 (576 atoms)', '3x3x8 (432 atoms)', '3x3x10 (540 atoms)', '3x3x12 (648 atoms)']
# labels_hse = ['3x3x4 (216 atoms)', '3x3x6 (324 atoms)', '3x3x8 (432 atoms)', '4x4x6 (576 atoms)', '3x3x10 (540 atoms)', '3x3x12 (648 atoms)']
# skip_folders_hse_start = 3
# skip_folders_hse_end = -1
# skip_folders_hse_start_gap = 3
# skip_folders_hse_end_gap = -1
# folder_use_hse = ['hfx-schwarz-e-1e-6-f-1e-6-fit9-pfit3']
# folder_files_hse = []
# folder_use_pbe0 = []
# # for i in range(len(folder_use_hse)): # from-md-1 better for non-linearity
# #     folder_files_hse.append(['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-334/geo-opt-cell-opt/neutral/neutral-opt/{}'.format(folder_use_hse[i]),
# #                         '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-334/geo-opt-cell-opt/electron/from-neutral/{}'.format(folder_use_hse[i]),
# #                         '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-334/geo-opt-cell-opt/electron/from-md-1/{}'.format(folder_use_hse[i]),
# #                         '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-334/geo-opt-cell-opt/neutral/from-electron//{}'.format(folder_use_hse[i])])
# # folder_use_hse = ['hfx-schwarz-e-1e-6-f-1e-6-fit9-pfit3-tight']
# for i in range(len(folder_use_hse)): # from-hse-30 better for trapping energy
#     folder_files_hse.append(['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-334/geo-opt-cell-opt/neutral/neutral-opt/{}'.format(folder_use_hse[i]),
#                         '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-334/geo-opt-cell-opt/electron/from-neutral/hfx-schwarz-e-1e-6-f-1e-6-fit9-pfit3-rs-scf-1e-7',
#                         '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-334/geo-opt-cell-opt/electron/from-hse-30/{}'.format(folder_use_hse[i]),
#                         '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-334/geo-opt-cell-opt/neutral/from-electron-hse-30//{}'.format(folder_use_hse[i])])
# # for i in range(len(folder_use_hse)):
# #     folder_files_hse.append(['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-335/geo-opt-cell-opt/neutral/neutral-opt/{}'.format(folder_use_hse[i]),
# #                         '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-335/geo-opt-cell-opt/electron/from-neutral/{}'.format(folder_use_hse[i]),
# #                         '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-335/geo-opt-cell-opt/electron/from-md-1/{}'.format(folder_use_hse[i]),
# #                         '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-335/geo-opt-cell-opt/neutral/from-electron//{}'.format(folder_use_hse[i])])
# folder_use_hse = ['hfx-schwarz-e-1e-6-f-1e-6-fit9-pfit3']
# # folder_use_hse = ['hfx-schwarz-e-1e-6-f-1e-6-fit9-pfit3-tight']
# # for i in range(len(folder_use_hse)):
# #     folder_files_hse.append(['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-336/geo-opt-cell-opt/neutral/neutral-opt/{}'.format(folder_use_hse[i]),
# #                         '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-336/geo-opt-cell-opt/electron/from-neutral/{}'.format(folder_use_hse[i]),
# #                         '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-336/geo-opt-cell-opt/electron/from-pbe0-14/hfx-schwarz-e-1e-6-f-1e-6-fit9-pfit3',
# #                         '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-336/geo-opt-cell-opt/neutral/from-electron-pbe0-14/hfx-schwarz-e-1e-6-f-1e-6-fit9-pfit3'])
# # for i in range(len(folder_use_hse)):  # best trapping energy
# #     folder_files_hse.append(['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-336/geo-opt-cell-opt/neutral/neutral-opt/{}'.format(folder_use_hse[i]),
# #                         '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-336/geo-opt-cell-opt/electron/from-neutral/{}'.format(folder_use_hse[i]),
# #                         '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-336/geo-opt-cell-opt/electron/from-md-atom-78/{}'.format(folder_use_hse[i]),
# #                         '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-336/geo-opt-cell-opt/neutral/from-electron/{}'.format(folder_use_hse[i])])
# for i in range(len(folder_use_hse)):
#     folder_files_hse.append(['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-336/geo-opt-cell-opt/neutral/neutral-opt/{}'.format(folder_use_hse[i]),
#                         '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-336/geo-opt-cell-opt/electron/from-neutral/hfx-schwarz-e-1e-6-f-1e-6-fit9-pfit3-rs-scf-1e-7',
#                         '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-336/geo-opt-cell-opt/electron/from-pbe0-14/hfx-schwarz-e-1e-6-f-1e-6-fit9-pfit3-tight',
#                         '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-336/geo-opt-cell-opt/neutral/from-electron-pbe0-14/hfx-schwarz-e-1e-6-f-1e-6-fit9-pfit3-tight'])
# folder_use_hse = ['hfx-schwarz-e-1e-6-f-1e-6-fit9-pfit3']
# for i in range(len(folder_use_hse)):
#     folder_files_hse.append(['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-338/geo-opt-cell-opt/neutral/neutral-opt/hfx-schwarz-e-1e-6-f-1e-6-fit9-pfit3-from-hse-15',
#                         '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-338/geo-opt-cell-opt/electron/from-neutral/hfx-schwarz-e-1e-6-f-1e-6-fit9-pfit3-from-hse-15-rs-scf-1e-7',
#                         '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-338/geo-opt-cell-opt/electron/from-offset-107-hse-25/{}'.format(folder_use_hse[i]),
#                         '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-338/geo-opt-cell-opt/neutral/from-offset-107-hse-25/{}'.format(folder_use_hse[i])])
# folder_use_hse = ['hfx-schwarz-e-1e-6-f-1e-6-fit9-pfit3']
# # folder_use_hse = ['hfx-schwarz-e-1e-6-f-1e-6-fit9-pfit3-tight']
# # folder_use_hse = ['hfx-schwarz-e-1e-6-f-1e-6-fit9-pfit3-rs-tz-ti-o']
# for i in range(len(folder_use_hse)):
#     folder_files_hse.append(['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-446/geo-opt-cell-opt/neutral/neutral-opt/{}'.format(folder_use_hse[i]),
#                         '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-446/geo-opt-cell-opt/electron/from-neutral/hfx-schwarz-e-1e-6-f-1e-6-fit9-pfit3-rs-scf-1e-7',
#                         '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-446/geo-opt-cell-opt/electron/from-pbe0-14/{}'.format(folder_use_hse[i]),
#                         '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-446/geo-opt-cell-opt/neutral/from-electron-pbe0-14/{}'.format(folder_use_hse[i])])
# for i in range(len(folder_use_hse)):
#     folder_files_hse.append(['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-3310/geo-opt-cell-opt/neutral/neutral-opt/hfx-schwarz-e-1e-6-f-1e-6-fit9-pfit3-from-hse-22',
#                         '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-3310/geo-opt-cell-opt/electron/from-neutral/hfx-schwarz-e-1e-6-f-1e-6-fit9-pfit3-from-hse-22-rs-scf-1e-7',
#                         '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-3310/geo-opt-cell-opt/electron/from-offset-135-hse-25/{}'.format(folder_use_hse[i]),
#                         '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-3310/geo-opt-cell-opt/neutral/from-offset-135-hse-25/{}'.format(folder_use_hse[i])])
# for i in range(len(folder_use_hse)):
#     folder_files_hse.append(['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-3312/geo-opt-cell-opt/neutral/neutral-opt/hfx-schwarz-e-1e-6-f-1e-6-fit9-pfit3-from-hse-20',
#                         '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-3312/geo-opt-cell-opt/electron/from-neutral/hfx-schwarz-e-1e-6-f-1e-6-fit9-pfit3-from-hse-20-rs-scf-1e-7',
#                         '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-3312/geo-opt-cell-opt/electron/from-offset-163-hse-25/{}'.format(folder_use_hse[i]),
#                         '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-3312/geo-opt-cell-opt/neutral/from-offset-163-hse-25/{}'.format(folder_use_hse[i])])
# folder_use_hse = ['hfx-schwarz-e-1e-6-f-1e-6-fit9-pfit3'] * 4
# # folder_use_hse = ['hfx-schwarz-e-1e-6-f-1e-6-fit9-pfit3'] * 6
# # hse_values_polaron = [[12, 15, 18, 20, 21, 22, 23, 24, 25, 30], [12, 15, 18, 20, 21, 22, 23, 24, 25, 30], [12, 15, 18, 20, 21, 22, 23, 24, 25, 30],
# #                       [12, 15, 18, 20, 21, 22, 23, 24, 25, 30], [12, 15, 18, 20, 21, 22, 23, 24, 25, 30], [12, 15, 18, 20, 21, 22, 23, 24, 25, 30]]
# hse_values_polaron = [[12, 15, 18, 20, 21, 22, 23, 24, 25, 30], [12, 15, 18, 20, 21, 22, 23, 24, 25, 30], [12, 15, 18, 20, 21, 22, 23, 24, 25, 30],
#                       [12, 15, 18, 20, 21, 22, 23, 24, 25, 30], [12, 15, 18, 20, 21, 22, 23, 24, 25, 30], [12, 15, 18, 20, 21, 22, 23, 24, 25, 30], [12, 15, 18, 20, 21, 22, 23, 24, 25, 30]]

# Hubbard U
data_hubbard_neutral = []
data_hubbard_vertical_from_neutral = []
data_hubbard_electron_offset = []
data_hubbard_vertical_from_electron = []
trapping_energy_hubbard = []
non_linearity_alpha_hubbard = []
ionisation_potential_hubbard_e = []
non_linearity_alpha_hubbard_e = []
non_linearity_beta_hubbard_e = []
ionisation_potential_hubbard_h = []
non_linearity_alpha_hubbard_h = []
non_linearity_beta_hubbard_h = []
non_linearity_homo_lumo_hubbard = []
for i in range(len(folder_use_hubbard)):
    data_hubbard_neutral.append(read_output(folder_files_hubbard[i][0], files[0]))
    data_hubbard_vertical_from_neutral.append(read_output(folder_files_hubbard[i][1], files[0]))
    data_hubbard_electron_offset.append(read_output(folder_files_hubbard[i][2], files[0]))
    data_hubbard_vertical_from_electron.append(read_output(folder_files_hubbard[i][3], files[0]))
    trapping_energy_hubbard.append(data_hubbard_vertical_from_neutral[i]['Energy_first'][:np.shape(hubbard_values_polaron[i])[0]] - data_hubbard_electron_offset[i]['Energy_last'][:np.shape(hubbard_values_polaron[i])[0]])
    # HOMO(N) = -IP = E(N) - E(N-1)
    ionisation_potential_hubbard_e.append((data_hubbard_electron_offset[i]['Energy_last'][:np.shape(hubbard_values_polaron[i])[0]] - data_hubbard_vertical_from_electron[i]['Energy_first'][:np.shape(hubbard_values_polaron[i])[0]]))
    non_linearity_alpha_hubbard_e.append(ionisation_potential_hubbard_e[i] - data_hubbard_electron_offset[i]['HOMO_alpha_last'][:np.shape(hubbard_values_polaron[i])[0]])
    non_linearity_beta_hubbard_e.append(ionisation_potential_hubbard_e[i] - data_hubbard_electron_offset[i]['HOMO_beta_last'][:np.shape(hubbard_values_polaron[i])[0]])
    # LUMO(N-1) = HOMO(N) on polaron geometry
    non_linearity_homo_lumo_hubbard.append(data_hubbard_vertical_from_electron[i]['LUMO_alpha_first'][:np.shape(hubbard_values_polaron[i])[0]] - data_hubbard_electron_offset[i]['HOMO_alpha_last'][:np.shape(hubbard_values_polaron[i])[0]])
    # non_linearity_homo_lumo_hubbard.append(data_hubbard_neutral[i]['LUMO_alpha_last'][:np.shape(hubbard_values_polaron[i])[0]] - data_hubbard_vertical_from_neutral[i]['HOMO_alpha_first'][:np.shape(hubbard_values_polaron[i])[0]])
    # LUMO(N-1) = -EA = E(N) - E(N-1)
    ionisation_potential_hubbard_h.append(data_hubbard_vertical_from_neutral[i]['Energy_first'][:np.shape(hubbard_values_polaron[i])[0]]-data_hubbard_neutral[i]['Energy_last'][:np.shape(hubbard_values_polaron[i])[0]])
    non_linearity_alpha_hubbard_h.append(ionisation_potential_hubbard_h[i] - data_hubbard_neutral[i]['LUMO_alpha_last'][:np.shape(hubbard_values_polaron[i])[0]])
    non_linearity_beta_hubbard_h.append(ionisation_potential_hubbard_h[i] - data_hubbard_neutral[i]['LUMO_beta_last'][:np.shape(hubbard_values_polaron[i])[0]])

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

# pbe0
data_pbe0_neutral = []
data_pbe0_vertical_from_neutral = []
data_pbe0_electron_offset = []
data_pbe0_vertical_from_electron = []
trapping_energy_pbe0 = []
non_linearity_alpha_pbe0 = []
ionisation_potential_pbe0_e = []
non_linearity_alpha_pbe0_e = []
non_linearity_beta_pbe0_e = []
ionisation_potential_pbe0_h = []
non_linearity_alpha_pbe0_h = []
non_linearity_beta_pbe0_h = []
non_linearity_homo_lumo_pbe0 = []
for i in range(len(folder_use_pbe0)):
    data_pbe0_neutral.append(read_output(folder_files_pbe0[i][0], files[0]))
    data_pbe0_vertical_from_neutral.append(read_output(folder_files_pbe0[i][1], files[0]))
    data_pbe0_electron_offset.append(read_output(folder_files_pbe0[i][2], files[0]))
    data_pbe0_vertical_from_electron.append(read_output(folder_files_pbe0[i][3], files[0]))
    trapping_energy_pbe0.append(data_pbe0_vertical_from_neutral[i]['Energy_first'][:np.shape(pbe0_values_polaron[i])[0]] - data_pbe0_electron_offset[i]['Energy_last'][:np.shape(pbe0_values_polaron[i])[0]])
    # HOMO(N) = -IP = E(N) - E(N-1)
    ionisation_potential_pbe0_e.append((data_pbe0_electron_offset[i]['Energy_last'][:np.shape(pbe0_values_polaron[i])[0]] - data_pbe0_vertical_from_electron[i]['Energy_first'][:np.shape(pbe0_values_polaron[i])[0]]))
    non_linearity_alpha_pbe0_e.append(ionisation_potential_pbe0_e[i] - data_pbe0_electron_offset[i]['HOMO_alpha_last'][:np.shape(pbe0_values_polaron[i])[0]])
    non_linearity_beta_pbe0_e.append(ionisation_potential_pbe0_e[i] - data_pbe0_electron_offset[i]['HOMO_beta_last'][:np.shape(pbe0_values_polaron[i])[0]])
    # LUMO(N-1) = HOMO(N) on polaron geometry
    non_linearity_homo_lumo_pbe0.append(data_pbe0_vertical_from_electron[i]['LUMO_alpha_first'][:np.shape(pbe0_values_polaron[i])[0]] - data_pbe0_electron_offset[i]['HOMO_alpha_last'][:np.shape(pbe0_values_polaron[i])[0]])
    # non_linearity_homo_lumo_pbe0.append(data_pbe0_neutral[i]['LUMO_alpha_last'][:np.shape(pbe0_values_polaron[i])[0]] - data_pbe0_vertical_from_neutral[i]['HOMO_alpha_first'][:np.shape(pbe0_values_polaron[i])[0]])
    # LUMO(N-1) = -EA = E(N) - E(N-1)
    ionisation_potential_pbe0_h.append(data_pbe0_vertical_from_neutral[i]['Energy_first'][:np.shape(pbe0_values_polaron[i])[0]]-data_pbe0_neutral[i]['Energy_last'][:np.shape(pbe0_values_polaron[i])[0]])
    non_linearity_alpha_pbe0_h.append(ionisation_potential_pbe0_h[i] - data_pbe0_neutral[i]['LUMO_alpha_last'][:np.shape(pbe0_values_polaron[i])[0]])
    non_linearity_beta_pbe0_h.append(ionisation_potential_pbe0_h[i] - data_pbe0_neutral[i]['LUMO_beta_last'][:np.shape(pbe0_values_polaron[i])[0]])

# folder_save = folder_files_hse[0]
folder_save = []
for i in range(np.shape(folder_files_hse)[0]):
    folder_save.append(folder_files_hse[i][0])
for i in range(np.shape(folder_files_pbe0)[0]):
    folder_save.append(folder_files_pbe0[i][0])
for i in range(np.shape(folder_files_hubbard)[0]):
    folder_save.append(folder_files_hubbard[i][0])

# Plot band gap Hubbard U
draw_legend = True
fig_band_gap1, ax_band_gap1 = plt.subplots()
ax_band_gap1.axhline(y=draw_gap, color='k', alpha=0.5)
for i in range(len(folder_use_hubbard)):
    ax_band_gap1.plot(hubbard_values_polaron[i], data_hubbard_neutral[i]['Gap_alpha'], 'x-', color=colors_hubbard[i], label='{}'.format(labels_hubbard[i]))
ax_band_gap1.set_ylabel('Band gap / eV')
ax_band_gap1.set_xlabel('Hubbard U / eV')
ax_band_gap1.set_ylim(set_ylim_band_gap_hubbard[0], set_ylim_band_gap_hubbard[1])
if draw_legend: ax_band_gap1.legend(frameon=False)
fig_band_gap1.tight_layout()
for i in range(np.shape(folder_save)[0]):
    fig_band_gap1.savefig('{}/band_gap_hubbard.png'.format(folder_save[i]), dpi=300)

# Plot band gap HFX
fig_band_gap2, ax_band_gap2 = plt.subplots()
ax_band_gap2.axhline(y=draw_gap, color='k', alpha=0.5)
for i in range(len(folder_use_hse)):
    ax_band_gap2.plot(hse_values_polaron[i][skip_folders_hse_start_gap:skip_folders_hse_end_gap], data_hse_neutral[i]['Gap_alpha'][skip_folders_hse_start_gap:skip_folders_hse_end_gap], 'x-', color=colors_hse[i], label='{}'.format(labels_hse[i]))
for i in range(len(folder_use_pbe0)):
    ax_band_gap2.plot(pbe0_values_polaron[i][skip_folders_pbe0_start_gap:skip_folders_pbe0_end_gap], data_pbe0_neutral[i]['Gap_alpha'][skip_folders_pbe0_start_gap:skip_folders_pbe0_end_gap], 'x-', color=colors_pbe0[i], label='{}'.format(labels_pbe0[i]))
# ax_band_gap2.plot(pbe0_values_polaron, data_pbe0_neutral['Gap_alpha'], 'gx-', label='PBE0')
ax_band_gap2.set_ylabel('Band gap / eV')
ax_band_gap2.set_xlabel('% HFX')
ax_band_gap2.set_ylim(set_ylim_band_gap_hfx[0], set_ylim_band_gap_hfx[1])
fig_band_gap2.tight_layout()
if draw_legend: ax_band_gap2.legend(frameon=False)
for i in range(np.shape(folder_save)[0]):
    fig_band_gap2.savefig('{}/band_gap_hfx.png'.format(folder_save[i]), dpi=300)

# Plot trapping energy Hubbard U
fig_trapping1, ax_trapping1 = plt.subplots()
for i in range(len(folder_use_hubbard)):
    ax_trapping1.plot(hubbard_values_polaron[i][skip_folders_hubbard_start:],
                      trapping_energy_hubbard[i][skip_folders_hubbard_start:]*param.hartree_to_ev*1e3, 'x-', color=colors_hubbard[i], label='{}'.format(labels_hubbard[i]))
ax_trapping1.set_ylabel('Trapping energy / meV')
ax_trapping1.set_xlabel('Hubbard U / eV')
if draw_legend: ax_trapping1.legend(frameon=False)
fig_trapping1.tight_layout()
for i in range(np.shape(folder_save)[0]):
    fig_trapping1.savefig('{}/trapping_hubbard.png'.format(folder_save[i]), dpi=300)

# Plot trapping energy HFX
fig_trapping2, ax_trapping2 = plt.subplots()
for i in range( len(folder_use_hse)):
    ax_trapping2.plot(hse_values_polaron[i][skip_folders_hse_start:skip_folders_hse_end],
                      trapping_energy_hse[i][skip_folders_hse_start:skip_folders_hse_end] * param.hartree_to_ev*1e3, 'x-', color=colors_hse[i], label='{}'.format(labels_hse[i]))
# skip_folders_hse_start_temp = skip_folders_hse_start + 1
# for i in range(0, 1):
#     ax_trapping2.plot(hse_values_polaron[i][skip_folders_hse_start_temp:skip_folders_hse_end],
#                       trapping_energy_hse[i][skip_folders_hse_start_temp:skip_folders_hse_end] * param.hartree_to_ev*1e3, 'x-', color=colors_hse[i], label='{}'.format(labels_hse[i]))
# for i in range(1, len(folder_use_hse) - 2):
#     ax_trapping2.plot(hse_values_polaron[i][skip_folders_hse_start:skip_folders_hse_end],
#                       trapping_energy_hse[i][skip_folders_hse_start:skip_folders_hse_end] * param.hartree_to_ev*1e3, 'x-', color=colors_hse[i], label='{}'.format(labels_hse[i]))
# skip_folders_hse_start_temp = skip_folders_hse_start + 2
# skip_folders_hse_end_temp = skip_folders_hse_end - 2
# for i in range(2, 4):
#     ax_trapping2.plot(hse_values_polaron[i][skip_folders_hse_start_temp:skip_folders_hse_end_temp],
#                       trapping_energy_hse[i][skip_folders_hse_start_temp:skip_folders_hse_end_temp] * param.hartree_to_ev*1e3, 'x-', color=colors_hse[i], label='{}'.format(labels_hse[i]))
for i in range(len(folder_use_pbe0)):
    ax_trapping2.plot(pbe0_values_polaron[i][skip_folders_pbe0_start:skip_folders_pbe0_end],
                      trapping_energy_pbe0[i][skip_folders_pbe0_start:skip_folders_pbe0_end] * param.hartree_to_ev*1e3, 'x-', color=colors_pbe0[i], label='{}'.format(labels_pbe0[i]))
# ax_trapping2.plot(pbe0_values_polaron[skip_folders_pbe0:], trapping_energy_pbe0[skip_folders_pbe0:] * param.hartree_to_ev*1e3, 'gx-', label='PBE0')
ax_trapping2.set_ylabel('Trapping energy / meV')
ax_trapping2.set_xlabel('% HFX')
ax_trapping2.set_ylim(set_ylim_trapping_hfx[0], set_ylim_trapping_hfx[1])
if draw_legend: ax_trapping2.legend(frameon=False)
fig_trapping2.tight_layout()
for i in range(np.shape(folder_save)[0]):
    fig_trapping2.savefig('{}/trapping_hfx.png'.format(folder_save[i]), dpi=300)

# Plot Hubbard U
fig_linearity1, ax_linearity1 = plt.subplots()
ax_linearity1.axhline(y=-draw_nonlinearity, color='k', alpha=0.5)
ax_linearity1.axhline(y=draw_nonlinearity, color='k', alpha=0.5)
for i in range(len(folder_use_hubbard)):
    # ax_linearity1.plot(hubbard_values_polaron[i][skip_folders_hubbard_start:skip_folders_hubbard_end],
    #                    non_linearity_homo_lumo_hubbard[i][skip_folders_hubbard_start:skip_folders_hubbard_end] * param.hartree_to_ev, 's-', fillstyle='none', color=colors_hubbard[i])
    ax_linearity1.plot(hubbard_values_polaron[i][skip_folders_hubbard_start:],
                       non_linearity_beta_hubbard_e[i][skip_folders_hubbard_start:] * param.hartree_to_ev, 'o-', fillstyle='none', color=colors_hubbard[i])
    # ax_linearity1.plot(hubbard_values_polaron[i][skip_folders_hubbard_start:skip_folders_hubbard_end],
    #                    non_linearity_beta_hubbard_h[i][skip_folders_hubbard_start:skip_folders_hubbard_end] * param.hartree_to_ev, 'x-', fillstyle='none', color=colors_hubbard[i], label='{}'.format(labels_hubbard[i]))
ax_linearity1.set_ylabel('Non-linearity / eV')
ax_linearity1.set_xlabel('Hubbard U / eV')
if draw_legend: ax_linearity1.legend(frameon=False)
# ax_linearity1.set_ylim(set_ylim_linearity_hubbard[0], set_ylim_linearity_hubbard[1])
fig_linearity1.tight_layout()
for i in range(np.shape(folder_save)[0]):
    fig_linearity1.savefig('{}/linearity_hubbard.png'.format(folder_save[i]), dpi=300)

# Plot non-linearity HFX
fig_linearity2, ax_linearity2 = plt.subplots()
ax_linearity2.axhline(y=draw_nonlinearity, color='k', alpha=0.5)
ax_linearity2.axhline(y=-draw_nonlinearity, color='k', alpha=0.5)
# ax_linearity2.axhline(y=0.00, color='k', alpha=0.5, linestyle='--')
for i in range(len(folder_use_hse)):
    # ax_linearity2.plot(hse_values_polaron[i][skip_folders_hse_start:skip_folders_hse_end],
    #                    non_linearity_homo_lumo_hse[i][skip_folders_hse_start:skip_folders_hse_end] * param.hartree_to_ev, 's-', fillstyle='none', color=colors_hse[i])
    ax_linearity2.plot(hse_values_polaron[i][skip_folders_hse_start:skip_folders_hse_end],
                       non_linearity_alpha_hse_e[i][skip_folders_hse_start:skip_folders_hse_end] * param.hartree_to_ev, 'o-', fillstyle='none', color=colors_hse[i], label='{}'.format(labels_hse[i]))
    # ax_linearity2.plot(hse_values_polaron[i][skip_folders_hse_start:skip_folders_hse_end],
    #                    non_linearity_alpha_hse_h[i][skip_folders_hse_start:skip_folders_hse_end] * param.hartree_to_ev, 'x-', fillstyle='none', color=colors_hse[i], label='{}'.format(labels_hse[i]))
# skip_folders_hse_start_temp = skip_folders_hse_start + 1
# for i in range(0, 1):
#     # ax_linearity2.plot(hse_values_polaron[i][skip_folders_hse_start:skip_folders_hse_end],
#     #                    non_linearity_homo_lumo_hse[i][skip_folders_hse_start:skip_folders_hse_end] * param.hartree_to_ev, 's-', fillstyle='none', color=colors_hse[i])
#     ax_linearity2.plot(hse_values_polaron[i][skip_folders_hse_start_temp:skip_folders_hse_end],
#                        non_linearity_alpha_hse_e[i][skip_folders_hse_start_temp:skip_folders_hse_end] * param.hartree_to_ev, 'o-', fillstyle='none', color=colors_hse[i], label='{}'.format(labels_hse[i]))
#     # ax_linearity2.plot(hse_values_polaron[i][skip_folders_hse_start:skip_folders_hse_end],
#     #                    non_linearity_alpha_hse_h[i][skip_folders_hse_start:skip_folders_hse_end] * param.hartree_to_ev, 'x-', fillstyle='none', color=colors_hse[i], label='{}'.format(labels_hse[i]))
# for i in range(1, len(folder_use_hse)-2):
#     # ax_linearity2.plot(hse_values_polaron[i][skip_folders_hse_start:skip_folders_hse_end],
#     #                    non_linearity_homo_lumo_hse[i][skip_folders_hse_start:skip_folders_hse_end] * param.hartree_to_ev, 's-', fillstyle='none', color=colors_hse[i])
#     ax_linearity2.plot(hse_values_polaron[i][skip_folders_hse_start:skip_folders_hse_end],
#                        non_linearity_alpha_hse_e[i][skip_folders_hse_start:skip_folders_hse_end] * param.hartree_to_ev, 'o-', fillstyle='none', color=colors_hse[i], label='{}'.format(labels_hse[i]))
#     # ax_linearity2.plot(hse_values_polaron[i][skip_folders_hse_start:skip_folders_hse_end],
#     #                    non_linearity_alpha_hse_h[i][skip_folders_hse_start:skip_folders_hse_end] * param.hartree_to_ev, 'x-', fillstyle='none', color=colors_hse[i], label='{}'.format(labels_hse[i]))
# skip_folders_hse_start_temp = skip_folders_hse_start + 2
# skip_folders_hse_end_temp = skip_folders_hse_end - 2
# for i in range(2, 4):
#     # ax_linearity2.plot(hse_values_polaron[i][skip_folders_hse_start:skip_folders_hse_end],
#     #                    non_linearity_homo_lumo_hse[i][skip_folders_hse_start:skip_folders_hse_end] * param.hartree_to_ev, 's-', fillstyle='none', color=colors_hse[i])
#     ax_linearity2.plot(hse_values_polaron[i][skip_folders_hse_start_temp:skip_folders_hse_end_temp],
#                        non_linearity_alpha_hse_e[i][skip_folders_hse_start_temp:skip_folders_hse_end_temp] * param.hartree_to_ev, 'o-', fillstyle='none', color=colors_hse[i], label='{}'.format(labels_hse[i]))
#     # ax_linearity2.plot(hse_values_polaron[i][skip_folders_hse_start:skip_folders_hse_end],
#     #                    non_linearity_alpha_hse_h[i][skip_folders_hse_start:skip_folders_hse_end] * param.hartree_to_ev, 'x-', fillstyle='none', color=colors_hse[i], label='{}'.format(labels_hse[i]))
for i in range(len(folder_use_pbe0)):
    # ax_linearity2.plot(pbe0_values_polaron[i][skip_folders_pbe0_start:skip_folders_pbe0_end],
    #                    non_linearity_homo_lumo_pbe0[i][skip_folders_pbe0_start:skip_folders_pbe0_end] * param.hartree_to_ev, 's-', fillstyle='none', color=colors_pbe0[i])
    ax_linearity2.plot(pbe0_values_polaron[i][skip_folders_pbe0_start:skip_folders_pbe0_end],
                       non_linearity_alpha_pbe0_e[i][skip_folders_pbe0_start:skip_folders_pbe0_end] * param.hartree_to_ev, 'o-', fillstyle='none', color=colors_pbe0[i], label='{}'.format(labels_pbe0[i]))
    # ax_linearity2.plot(pbe0_values_polaron[i][skip_folders_pbe0_start:skip_folders_pbe0_end],
    #                    non_linearity_alpha_pbe0_h[i][skip_folders_pbe0_start:skip_folders_pbe0_end] * param.hartree_to_ev, 'x-', fillstyle='none', color=colors_pbe0[i], label='{}'.format(labels_pbe0[i]))
ax_linearity2.set_ylabel('Non-linearity / eV')
ax_linearity2.set_ylim(set_ylim_linearity_hfx[0], set_ylim_linearity_hfx[1])
ax_linearity2.set_xlabel('% HFX')
if draw_legend: ax_linearity2.legend(frameon=False)
fig_linearity2.tight_layout()
for i in range(np.shape(folder_save)[0]):
    fig_linearity2.savefig('{}/linearity_hfx.png'.format(folder_save[i]), dpi=300)


rows = 1
cols = 2
fig_plot_nonlinearity_gap, ax_plot_nonlinearity_gap = plt.subplots(rows, cols, figsize=(12, 6))
# fig_plot_nonlinearity_gap, ax_plot_nonlinearity_gap = plt.subplots(rows, cols, figsize=(12, 5))

ax_plot_nonlinearity_gap[0].axhline(y=draw_gap, color='k', alpha=0.5)
for i in range(len(folder_use_hse)):
    ax_plot_nonlinearity_gap[0].plot(hse_values_polaron[i][skip_folders_hse_start_gap:skip_folders_hse_end_gap], data_hse_neutral[i]['Gap_alpha'][skip_folders_hse_start_gap:skip_folders_hse_end_gap], 'x-', color=colors_hse[i], label='{}'.format(labels_hse[i]))
for i in range(len(folder_use_pbe0)):
    ax_plot_nonlinearity_gap[0].plot(pbe0_values_polaron[i][skip_folders_pbe0_start_gap:skip_folders_pbe0_end_gap], data_pbe0_neutral[i]['Gap_alpha'][skip_folders_pbe0_start_gap:skip_folders_pbe0_end_gap], 'x-', color=colors_pbe0[i], label='{}'.format(labels_pbe0[i]))
# ax_band_gap2.plot(pbe0_values_polaron, data_pbe0_neutral['Gap_alpha'], 'gx-', label='PBE0')
ax_plot_nonlinearity_gap[0].set_ylabel('Band gap / eV')
ax_plot_nonlinearity_gap[0].set_xlabel('% HFX')
ax_plot_nonlinearity_gap[0].set_ylim(set_ylim_band_gap_hfx[0], set_ylim_band_gap_hfx[1])
if draw_legend: ax_plot_nonlinearity_gap[0].legend(frameon=False)

ax_plot_nonlinearity_gap[1].axhline(y=draw_nonlinearity, color='k', alpha=0.5)
ax_plot_nonlinearity_gap[1].axhline(y=-draw_nonlinearity, color='k', alpha=0.5)
for i in range(len(folder_use_hse)):
    ax_plot_nonlinearity_gap[1].plot(hse_values_polaron[i][skip_folders_hse_start:skip_folders_hse_end],
                       non_linearity_alpha_hse_e[i][skip_folders_hse_start:skip_folders_hse_end] * param.hartree_to_ev, 'o-', fillstyle='none', color=colors_hse[i], label='{}'.format(labels_hse[i]))
for i in range(len(folder_use_pbe0)):
    ax_plot_nonlinearity_gap[1].plot(pbe0_values_polaron[i][skip_folders_pbe0_start:skip_folders_pbe0_end],
                       non_linearity_alpha_pbe0_e[i][skip_folders_pbe0_start:skip_folders_pbe0_end] * param.hartree_to_ev, 'o-', fillstyle='none', color=colors_pbe0[i], label='{}'.format(labels_pbe0[i]))
ax_plot_nonlinearity_gap[1].set_ylabel('Non-linearity / eV')
ax_plot_nonlinearity_gap[1].set_ylim(set_ylim_linearity_hfx[0], set_ylim_linearity_hfx[1])
ax_plot_nonlinearity_gap[1].set_xlabel('% HFX')
if draw_legend: ax_plot_nonlinearity_gap[1].legend(frameon=False)

fig_plot_nonlinearity_gap.tight_layout()
for i in range(np.shape(folder_save)[0]):
    fig_plot_nonlinearity_gap.savefig('{}/subplot_nonlinearity_gap.png'.format(folder_save[i]), dpi=300)

rows = 1
cols = 2
fig_plot_nonlinearity_trapping, ax_plot_nonlinearity_trapping = plt.subplots(rows, cols, figsize=(12, 6))
# fig_plot_nonlinearity_trapping, ax_plot_nonlinearity_trapping = plt.subplots(rows, cols, figsize=(12, 5))

for i in range( len(folder_use_hse)):
    ax_plot_nonlinearity_trapping[0].plot(hse_values_polaron[i][skip_folders_hse_start:skip_folders_hse_end],
                      trapping_energy_hse[i][skip_folders_hse_start:skip_folders_hse_end] * param.hartree_to_ev*1e3, 'x-', color=colors_hse[i], label='{}'.format(labels_hse[i]))
for i in range(len(folder_use_pbe0)):
    ax_plot_nonlinearity_trapping[0].plot(pbe0_values_polaron[i][skip_folders_pbe0_start:skip_folders_pbe0_end],
                      trapping_energy_pbe0[i][skip_folders_pbe0_start:skip_folders_pbe0_end] * param.hartree_to_ev*1e3, 'x-', color=colors_pbe0[i], label='{}'.format(labels_pbe0[i]))
ax_plot_nonlinearity_trapping[0].set_ylabel('Trapping energy / meV')
ax_plot_nonlinearity_trapping[0].set_xlabel('% HFX')
ax_plot_nonlinearity_trapping[0].set_ylim(set_ylim_trapping_hfx[0], set_ylim_trapping_hfx[1])
if draw_legend: ax_plot_nonlinearity_trapping[0].legend(frameon=False)

ax_plot_nonlinearity_trapping[1].axhline(y=draw_nonlinearity, color='k', alpha=0.5)
ax_plot_nonlinearity_trapping[1].axhline(y=-draw_nonlinearity, color='k', alpha=0.5)
for i in range(len(folder_use_hse)):
    ax_plot_nonlinearity_trapping[1].plot(hse_values_polaron[i][skip_folders_hse_start:skip_folders_hse_end],
                       non_linearity_alpha_hse_e[i][skip_folders_hse_start:skip_folders_hse_end] * param.hartree_to_ev, 'o-', fillstyle='none', color=colors_hse[i], label='{}'.format(labels_hse[i]))
for i in range(len(folder_use_pbe0)):
    ax_plot_nonlinearity_trapping[1].plot(pbe0_values_polaron[i][skip_folders_pbe0_start:skip_folders_pbe0_end],
                       non_linearity_alpha_pbe0_e[i][skip_folders_pbe0_start:skip_folders_pbe0_end] * param.hartree_to_ev, 'o-', fillstyle='none', color=colors_pbe0[i], label='{}'.format(labels_pbe0[i]))
ax_plot_nonlinearity_trapping[1].set_ylabel('Non-linearity / eV')
ax_plot_nonlinearity_trapping[1].set_ylim(set_ylim_linearity_hfx[0], set_ylim_linearity_hfx[1])
ax_plot_nonlinearity_trapping[1].set_xlabel('% HFX')
if draw_legend: ax_plot_nonlinearity_trapping[1].legend(frameon=False)

fig_plot_nonlinearity_trapping.tight_layout()
for i in range(np.shape(folder_save)[0]):
    fig_plot_nonlinearity_trapping.savefig('{}/subplot_nonlinearity_trapping.png'.format(folder_save[i]), dpi=300)


rows = 1
cols = 3
# fig_plot_nonlinearity_gap, ax_plot_nonlinearity_gap = plt.subplots(rows, cols, figsize=(12, 6))
fig_plot_nonlinearity_gap_trapping, ax_plot_nonlinearity_gap_trapping = plt.subplots(rows, cols, figsize=(18, 6))

ax_plot_nonlinearity_gap_trapping[0].axhline(y=draw_gap, color='k', alpha=0.5)
for i in range(len(folder_use_hse)):
    ax_plot_nonlinearity_gap_trapping[0].plot(hse_values_polaron[i][skip_folders_hse_start_gap:skip_folders_hse_end_gap], data_hse_neutral[i]['Gap_alpha'][skip_folders_hse_start_gap:skip_folders_hse_end_gap], 'x-', color=colors_hse[i], label='{}'.format(labels_hse[i]))
for i in range(len(folder_use_pbe0)):
    ax_plot_nonlinearity_gap_trapping[0].plot(pbe0_values_polaron[i][skip_folders_pbe0_start_gap:skip_folders_pbe0_end_gap], data_pbe0_neutral[i]['Gap_alpha'][skip_folders_pbe0_start_gap:skip_folders_pbe0_end_gap], 'x-', color=colors_pbe0[i], label='{}'.format(labels_pbe0[i]))
# ax_band_gap2.plot(pbe0_values_polaron, data_pbe0_neutral['Gap_alpha'], 'gx-', label='PBE0')
ax_plot_nonlinearity_gap_trapping[0].set_ylabel('Band gap / eV')
ax_plot_nonlinearity_gap_trapping[0].set_xlabel('% HFX')
ax_plot_nonlinearity_gap_trapping[0].set_ylim(set_ylim_band_gap_hfx[0], set_ylim_band_gap_hfx[1])
if draw_legend: ax_plot_nonlinearity_gap_trapping[0].legend(frameon=False)

ax_plot_nonlinearity_gap_trapping[1].axhline(y=draw_nonlinearity, color='k', alpha=0.5)
ax_plot_nonlinearity_gap_trapping[1].axhline(y=-draw_nonlinearity, color='k', alpha=0.5)
for i in range(len(folder_use_hse)):
    ax_plot_nonlinearity_gap_trapping[1].plot(hse_values_polaron[i][skip_folders_hse_start:skip_folders_hse_end],
                       non_linearity_alpha_hse_e[i][skip_folders_hse_start:skip_folders_hse_end] * param.hartree_to_ev, 'o-', fillstyle='none', color=colors_hse[i], label='{}'.format(labels_hse[i]))
for i in range(len(folder_use_pbe0)):
    ax_plot_nonlinearity_gap_trapping[1].plot(pbe0_values_polaron[i][skip_folders_pbe0_start:skip_folders_pbe0_end],
                       non_linearity_alpha_pbe0_e[i][skip_folders_pbe0_start:skip_folders_pbe0_end] * param.hartree_to_ev, 'o-', fillstyle='none', color=colors_pbe0[i], label='{}'.format(labels_pbe0[i]))
ax_plot_nonlinearity_gap_trapping[1].set_ylabel('Non-linearity / eV')
ax_plot_nonlinearity_gap_trapping[1].set_ylim(set_ylim_linearity_hfx[0], set_ylim_linearity_hfx[1])
ax_plot_nonlinearity_gap_trapping[1].set_xlabel('% HFX')
if draw_legend: ax_plot_nonlinearity_gap_trapping[1].legend(frameon=False)

for i in range( len(folder_use_hse)):
    ax_plot_nonlinearity_gap_trapping[2].plot(hse_values_polaron[i][skip_folders_hse_start:skip_folders_hse_end],
                      trapping_energy_hse[i][skip_folders_hse_start:skip_folders_hse_end] * param.hartree_to_ev*1e3, 'x-', color=colors_hse[i], label='{}'.format(labels_hse[i]))
for i in range(len(folder_use_pbe0)):
    ax_plot_nonlinearity_gap_trapping[2].plot(pbe0_values_polaron[i][skip_folders_pbe0_start:skip_folders_pbe0_end],
                      trapping_energy_pbe0[i][skip_folders_pbe0_start:skip_folders_pbe0_end] * param.hartree_to_ev*1e3, 'x-', color=colors_pbe0[i], label='{}'.format(labels_pbe0[i]))
ax_plot_nonlinearity_gap_trapping[2].set_ylabel('Trapping energy / meV')
ax_plot_nonlinearity_gap_trapping[2].set_xlabel('% HFX')
ax_plot_nonlinearity_gap_trapping[2].set_ylim(set_ylim_trapping_hfx[0], set_ylim_trapping_hfx[1])
if draw_legend: ax_plot_nonlinearity_gap_trapping[2].legend(frameon=False)

fig_plot_nonlinearity_gap_trapping.tight_layout()
for i in range(np.shape(folder_save)[0]):
    fig_plot_nonlinearity_gap_trapping.savefig('{}/subplot_nonlinearity_gap_trapping.png'.format(folder_save[i]), dpi=300)
    
if __name__ == "__main__":
    print('Finished.')
    plt.show()

