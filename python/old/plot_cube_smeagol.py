import numpy as np
from matplotlib import pyplot as plt
from ase.io.cube import read_cube_data
from general import parameters as param

""" Plotting of CP2K .cube files for Hartree potential and charge density """

params = {'axes.formatter.limits': [-4, 4],
          'axes.labelsize': 'x-large',
          'axes.titlesize': 'x-large',
          'legend.fontsize': 'large',
          'lines.markersize': '8',
          }
plt.rcParams.update(params)

# Defaults
plot_fermi = False
plot_dft = True
use_xlim = False
draw_mirror = False
mirror_scale = 1
draw_markers = False
plot_diff = True
plot_leads = False
fermi_dft = 0
labels = ['V=0.cube', 'V=1.cube', 'Bulk.cube']
diff_label = ['L', 'R', 'W(L, R)']
file_charge_em2 = ['1V-ELECTRON_DENSITY-1_0.cube'] * 5
file_hartree_em2 = ['1V-v_hartree-1_0.cube'] * 5
factor = [1] * 5
diff_color = ['r', 'g', 'k']

# au capacitor
# folder_em1 = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/testing/delete/cp2k-smeagol/kpoints-2-2-20_V-1_double-contour-rerun'] * 5
# folder_em2 = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/testing/delete/cp2k-smeagol/kpoints-2-2-20_V-0.001_double-contour-rerun',
#               '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/testing/delete/cp2k-smeagol/kpoints-2-2-20_V-0.01_double-contour-rerun',
#               '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/testing/delete/cp2k-smeagol/kpoints-2-2-20_V-0.1_double-contour-rerun',
#               '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/testing/delete/cp2k-smeagol/kpoints-2-2-20_V-1_double-contour-rerun',
#               '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/testing/delete/cp2k-smeagol/kpoints-2-2-20_V-4_double-contour-rerun']
# diff_label = ['0.001', '0.01', '0.1', '1', '4']
# factor = [float(s) for s in diff_label]
# file_charge_em2 = []
# file_hartree_em2 = []
# for i in range(len(diff_label)):
#     file_charge_em2.append('{}V-ELECTRON_DENSITY-1_0.cube'.format(diff_label[i]))
#     file_hartree_em2.append('{}V-v_hartree-1_0.cube'.format(diff_label[i]))
# print(file_charge_em2)
# diff_color = ['r', 'g', 'b', 'm', 'k']

# au capacitor
# folder_em1 = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/testing/delete/cp2k-smeagol/kpoints-2-2-20_V-1_double-contour-rerun'] * 5
# folder_em2 = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/testing/delete/cp2k-smeagol/kpoints-2-2-20_V-0.01_WeightRho-0-rerun',
#               '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/testing/delete/cp2k-smeagol/kpoints-2-2-20_V-0.01_WeightRho-1-rerun',
#               '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/testing/delete/cp2k-smeagol/kpoints-2-2-20_V-0.01_WeightRho-0.5-rerun',
#               '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/testing/delete/cp2k-smeagol/kpoints-2-2-20_V-0.01_double-contour-rerun']
# diff_label = ['L', 'R', '1/2 L + 1/2 R',  'W(L, R)']
# file_charge_em2 = ['0.01V-ELECTRON_DENSITY-1_0.cube'] * 5
# file_hartree_em2 = ['0.01V-v_hartree-1_0.cube'] * 5
# diff_color = ['r', 'g', 'b', 'm', 'k']

# cu chain
# folder_em1 = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/chain/centred/wire-q1-f3'] * 5
# folder_em2 = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/chain/centred/wire-q1-f3-WeightRho-0-cp2k2024',
#               '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/chain/centred/wire-q1-f3-WeightRho-1-cp2k2024',
#               '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/chain/centred/wire-q1-f3']
#
# folder_em1 = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/chain/centred/wire-q1-f3-two-cu-removed-middle'] * 5
# folder_em2 = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/chain/centred/wire-q1-f3-two-cu-removed-middle-WeightRho-0-cp2k2024',
#               '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/chain/centred/wire-q1-f3-two-cu-removed-middle-WeightRho-1-cp2k2024',
#               '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/chain/centred/wire-q1-f3-two-cu-removed-middle-WeightRho-0.5-cp2k2024',
#               '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/chain/centred/wire-q1-f3-two-cu-removed-middle']
# diff_label = ['L', 'R', '1/2 (L + R)',  'W(L, R)']
# diff_color = ['r', 'g', 'b', 'k']

# folder_em1 = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/chain/centred/wire-q1-f3-two-cu-removed-middle-V-0.01'] * 5
# folder_em2 = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/chain/centred/wire-q1-f3-two-cu-removed-middle-WeightRho-0-cp2k2024-V-0.01',
#               '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/chain/centred/wire-q1-f3-two-cu-removed-middle-WeightRho-0.5-cp2k2024-V-0.01',
#               '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/chain/centred/wire-q1-f3-two-cu-removed-middle-V-0.01']
# diff_label = ['L', 'R', '1/2 (L + R)',  'W(L, R)']
# diff_color = ['r', 'g', 'b', 'k']
# file_charge_em2 = ['0.01V-ELECTRON_DENSITY-1_0.cube'] * 5
# file_hartree_em2 = ['0.01V-v_hartree-1_0.cube'] * 5
#
# folder_em1 = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/chain/centred/wire-q1-f3-two-cu-removed-middle-V-0.01'] * 5
# folder_em2 = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/chain/centred/wire-q1-f3-two-cu-removed-middle-WeightRho-0-cp2k2024-V-0.01',
#               '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/chain/centred/wire-q1-f3-two-cu-removed-middle-WeightRho-1-cp2k2024-V-0.01',
#               '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/chain/centred/wire-q1-f3-two-cu-removed-middle-V-0.01']
# file_charge_em2 = ['0.01V-ELECTRON_DENSITY-1_0.cube'] * 5
# file_hartree_em2 = ['0.01V-v_hartree-1_0.cube'] * 5

# folder_em1 = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/chain/centred/wire-q1-f3-broken-left-right-square-6atoms'] * 5
# folder_em2 = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/chain/centred/wire-q1-f3-broken-left-right-sc-WeightRho-0-6atoms-cp2k2024',
#               '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/chain/centred/wire-q1-f3-broken-left-right-sc-WeightRho-1-6atoms-cp2k2024',
#               '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/chain/centred/wire-q1-f3-broken-left-right-square-6atoms']
#
# folder_em1 = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/chain/centred/wire-q1-f3-broken-left-square-6atoms'] * 5
# folder_em2 = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/chain/centred/wire-q1-f3-broken-left-sc-WeightRho-0-6atoms-cp2k2024',
#               '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/chain/centred/wire-q1-f3-broken-left-sc-WeightRho-1-6atoms-cp2k2024',
#               '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/chain/centred/wire-q1-f3-broken-left-sc-WeightRho-0.5-6atoms-cp2k2024',
#               '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/chain/centred/wire-q1-f3-broken-left-square-6atoms']
# diff_label = ['L', 'R', '1/2 (L + R)',  'W(L, R)']
# diff_color = ['r', 'g', 'b', 'k']

# folder_em1 = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/chain/centred/wire-q1-f3-broken-middle-square-6atoms'] * 5
# folder_em2 = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/chain/centred/wire-q1-f3-broken-middle-sc-WeightRho-0-6atoms-cp2k2024',
#               '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/chain/centred/wire-q1-f3-broken-middle-sc-WeightRho-1-6atoms-cp2k2024',
#               '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/chain/centred/wire-q1-f3-broken-middle-sc-WeightRho-0.5-6atoms-cp2k2024',
#               '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/chain/centred/wire-q1-f3-broken-middle-square-6atoms']
# diff_label = ['L', 'R', '1/2 (L + R)',  'W(L, R)']
# diff_color = ['r', 'g', 'b', 'k']

# folder_em1 = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/chain/centred/wire-q1-f3-broken-middle-square-6atoms'] * 5
# folder_em2 = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/chain/centred/wire-q1-f3-broken-middle-sc-WeightRho-0-6atoms-cp2k2024',
#               '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/chain/centred/wire-q1-f3-broken-middle-sc-WeightRho-1-6atoms-cp2k2024',
#               '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/chain/centred/wire-q1-f3-broken-middle-sc-WeightRho-0.5-6atoms-cp2k2024',
#               '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/chain/centred/wire-q1-f3-broken-middle-mulliken-6atoms',
#               '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/chain/centred/wire-q1-f3-broken-middle-square-6atoms']
# diff_label = ['L', 'R', '1/2 (L + R)', 'WM(L, R)', 'W(L, R)']
# diff_color = ['r', 'g', 'b', 'm', 'k']

# folder_em1 = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/chain/centred/wire-q1-f3-broken-middle-q11-square-6atoms'] * 5
# folder_em2 = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/chain/centred/wire-q1-f3-broken-middle-q11-sc-WeightRho-0-6atoms-cp2k2024',
#               '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/chain/centred/wire-q1-f3-broken-middle-q11-sc-WeightRho-1-6atoms-cp2k2024',
#               '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/chain/centred/wire-q1-f3-broken-middle-q11-square-6atoms']
#
# folder_em1 = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/chain/centred/capacitor-q1-f3-6atoms'] * 5
# folder_em2 = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/chain/centred/capacitor-q1-f3-6atoms-WeightRho-0-cp2k2024',
#               '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/chain/centred/capacitor-q1-f3-6atoms-WeightRho-1-cp2k2024',
#               '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/chain/centred/capacitor-q1-f3-6atoms-WeightRho-0-cp2k2024']

# supercell-1-1-3-bulk-6-cu-1.86
# xlim2 = [21.66000-0.5, 40.49800+0.5]
# folder_em1 = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/supercell-1-1-3-bulk-6-cu-1.86/junction/bias/kpoints-2-2-V-0'] * 5
# folder_em2 = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/supercell-1-1-3-bulk-6-cu-1.86/junction/bias/kpoints-2-2-V-0.01-WeightRho-0',
#               '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/supercell-1-1-3-bulk-6-cu-1.86/junction/bias/kpoints-2-2-V-0.01-WeightRho-1',
#               '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/supercell-1-1-3-bulk-6-cu-1.86/junction/bias/kpoints-2-2-V-0.01-dc']
# file_charge_em2 = ['0.01V-ELECTRON_DENSITY-1_0.cube'] * 5
# file_hartree_em2 = ['0.01V-v_hartree-1_0.cube'] * 5
# diff_label = ['L', 'R', 'W(L, R)']
# diff_color = ['r', 'g', 'k']

# folder_em1 = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/supercell-1-1-3-bulk-6-cu-1.86/junction/bias/kpoints-2-2-V-0'] * 5
# folder_em2 = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/supercell-1-1-3-bulk-6-cu-1.86/junction/bias/kpoints-2-2-V-0.1-WeightRho-0',
#               '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/supercell-1-1-3-bulk-6-cu-1.86/junction/bias/kpoints-2-2-V-0.1-WeightRho-1',
#               '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/supercell-1-1-3-bulk-6-cu-1.86/junction/bias/kpoints-2-2-V-0.1-dc']
# file_charge_em2 = ['0.1V-ELECTRON_DENSITY-1_0.cube'] * 5
# file_hartree_em2 = ['0.1V-v_hartree-1_0.cube'] * 5
# diff_label = ['L', 'R', 'W(L, R)']
# diff_color = ['r', 'g', 'k']

# folder_em1 = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/supercell-1-1-3-bulk-6-cu-1.86/junction/bias/kpoints-2-2-V-0'] * 5
# folder_em2 = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/supercell-1-1-3-bulk-6-cu-1.86/junction/bias/kpoints-2-2-V-1-WeightRho-0',
#               '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/supercell-1-1-3-bulk-6-cu-1.86/junction/bias/kpoints-2-2-V-1-WeightRho-1',
#               '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/supercell-1-1-3-bulk-6-cu-1.86/junction/bias/kpoints-2-2-V-1-dc']
# file_charge_em2 = ['1V-ELECTRON_DENSITY-1_0.cube'] * 5
# file_hartree_em2 = ['1V-v_hartree-1_0.cube'] * 5
# diff_label = ['L', 'R', 'W(L, R)']
# diff_color = ['r', 'g', 'k']

# supercell-1-1-4-bulk-6-cu-1.86
xlim2 = [21.66000-0.5, 45.54800+0.5]
folder_em1 = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/supercell-1-1-4-bulk-6-cu-1.86/junction/bias/kpoints-2-2-V-0'] * 5
folder_em2 = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/supercell-1-1-4-bulk-6-cu-1.86/junction/bias/kpoints-2-2-V-1-sc-WeightRho-0',
              '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/supercell-1-1-4-bulk-6-cu-1.86/junction/bias/kpoints-2-2-V-1-sc-WeightRho-1',
              '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/supercell-1-1-4-bulk-6-cu-1.86/junction/bias/kpoints-2-2-V-1-dc']
file_charge_em2 = ['1V-ELECTRON_DENSITY-1_0.cube'] * 5
file_hartree_em2 = ['1V-v_hartree-1_0.cube'] * 5
diff_label = ['L', 'R', 'W(L, R)']
diff_color = ['r', 'g', 'k']

# supercell-1-1-5-bulk-6-cu-1.86 V=1
xlim2 = [21.66000-0.5, 50.59700+0.5]
folder_em1 = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/supercell-1-1-5-bulk-6-cu-1.86/junction/bias/energy/kpoints-2-2-V-0'] * 5
folder_em2 = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/supercell-1-1-5-bulk-6-cu-1.86/junction/bias/energy/kpoints-2-2-V-0.1']
file_charge_em2 = ['0.1V-ELECTRON_DENSITY-1_0.cube'] * 5
file_hartree_em2 = ['0.1V-v_hartree-1_0.cube'] * 5
diff_label = ['L', 'R', 'W(L, R)']
diff_color = ['r', 'g', 'k']

# folder_em1 = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/supercell-1-1-5-bulk-6-cu-1.86/junction/bias/energy/kpoints-2-2-V-0'] * 5
# folder_em2 = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/supercell-1-1-5-bulk-6-cu-1.86/junction/bias/energy/kpoints-2-2-V-0-rs-V-1-rs-sc-WeightRho-0',
#               '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/supercell-1-1-5-bulk-6-cu-1.86/junction/bias/energy/kpoints-2-2-V-0-rs-V-1-rs-sc-WeightRho-1',
#               '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/supercell-1-1-5-bulk-6-cu-1.86/junction/bias/energy/kpoints-2-2-V-0-rs-V-1']
# file_charge_em2 = ['1V-ELECTRON_DENSITY-1_0.cube'] * 5
# file_hartree_em2 = ['1V-v_hartree-1_0.cube'] * 5
# diff_label = ['L', 'R', 'W(L, R)']
# diff_color = ['r', 'g', 'k']


# folder_em1 = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/supercell-1-1-5-bulk-6-cu-1.86/junction/bias/energy/kpoints-2-2-V-0'] * 5
# folder_em2 = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/supercell-1-1-5-bulk-6-cu-1.86/junction/bias/energy/kpoints-2-2-V-0-rs-V-1',
#               '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/supercell-1-1-5-bulk-6-cu-1.86/junction/bias/energy/kpoints-2-2-V-0-rs-V-1-rs-sc-WeightRho-0',
#               '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/supercell-1-1-5-bulk-6-cu-1.86/junction/bias/energy/kpoints-2-2-V-0-rs-V-1-rs-sc-WeightRho-1']
# file_charge_em2 = ['1V-ELECTRON_DENSITY-1_0.cube'] * 5
# file_hartree_em2 = ['1V-v_hartree-1_0.cube'] * 5

zoom = True

file_charge_em1 = ['0V-ELECTRON_DENSITY-1_0.cube'] * 5
file_hartree_em1 = ['0V-v_hartree-1_0.cube'] * 5
file_charge_dft = ['dft_wfn-ELECTRON_DENSITY-1_0.cube'] * 4
file_hartree_dft = ['dft_wfn-v_hartree-1_0.cube'] * 4
file_charge_bulk = ['bulk-ELECTRON_DENSITY-1_0.cube'] * 4
file_hartree_bulk = ['bulk-v_hartree-1_0.cube'] * 4
labels = ['Double contour', 'Left', 'Right']
plot_color = ['k', 'b', 'r', 'g']
read_labels = [False, False, True, True]
plot_labels = read_labels
print_label = 'average_z_all'
axis = 2
plot_diff_legend = [True, False]
draw_mirror = False
draw_mirror_diff = True
mirror_scale = -1
use_xlim = False
xlim_specify_left = [-0.1, 6]
xlim_specify_right = [78, 84.1]
xlim_specify = xlim_specify_left
draw_markers = False
z_bond_length = 2.084
z_num = 6
z_right = 21.422
mid_pos = z_bond_length*(z_num - 1)+(z_right-z_bond_length*(z_num - 1))/2
print(mid_pos)

# Read .cube using ASE
data_charge_em2 = []
data_hartree_em2 = []
data_charge_em1 = []
data_hartree_em1 = []
data_charge_bulk = []
data_charge_dft = []
data_hartree_dft = []
data_hartree_bulk = []
for i in range(len(folder_em2)):
    print('reading .cube files', folder_em2[i])
    if read_labels[3]: data_charge_em2.append(read_cube_data('{}/{}'.format(folder_em2[i], file_charge_em2[i])))
    if read_labels[3]: data_hartree_em2.append(read_cube_data('{}/{}'.format(folder_em2[i], file_hartree_em2[i])))
    if read_labels[2]: data_charge_em1.append(read_cube_data('{}/{}'.format(folder_em1[i], file_charge_em1[i])))
    if read_labels[2]: data_hartree_em1.append(read_cube_data('{}/{}'.format(folder_em1[i], file_hartree_em1[i])))
    if read_labels[1]: data_charge_dft.append(read_cube_data('{}/{}'.format(folder_em2[i], file_charge_dft[i])))
    if read_labels[1]: data_hartree_dft.append(read_cube_data('{}/{}'.format(folder_em2[i], file_hartree_dft[i])))
    if read_labels[0]: data_charge_bulk.append(read_cube_data('{}/{}'.format(folder_em2[i], file_charge_bulk[i])))
    if read_labels[0]: data_hartree_bulk.append(read_cube_data('{}/{}'.format(folder_em2[i], file_hartree_bulk[i])))
print('Finished reading .cube files')

# Calculate average along axis
average_charge_em2 = []
average_hartree_em2 = []
average_charge_em1 = []
average_hartree_em1 = []
average_charge_bulk = []
average_charge_dft = []
average_hartree_dft = []
average_hartree_bulk = []
for i in range(len(folder_em2)):
    if read_labels[3]: average_charge_em2.append(np.zeros(data_charge_em2[i][0].shape[axis]))
    if read_labels[3]: average_hartree_em2.append(np.zeros(data_hartree_em2[i][0].shape[axis]))
    if read_labels[2]: average_charge_em1.append(np.zeros(data_charge_em1[i][0].shape[axis]))
    if read_labels[2]: average_hartree_em1.append(np.zeros(data_hartree_em1[i][0].shape[axis]))
    if read_labels[1]: average_charge_dft.append(np.zeros(data_charge_dft[i][0].shape[axis]))
    if read_labels[1]: average_hartree_dft.append(np.zeros(data_hartree_dft[i][0].shape[axis]))
    if read_labels[0]: average_charge_bulk.append(np.zeros(data_charge_bulk[i][0].shape[axis]))
    if read_labels[0]: average_hartree_bulk.append(np.zeros(data_hartree_bulk[i][0].shape[axis]))

    if read_labels[0]:
        for j in range(data_charge_bulk[i][0].shape[axis]):
            average_charge_bulk[i][j] = np.mean(data_charge_bulk[i][0][:, :, j])
            average_hartree_bulk[i][j] = np.mean(data_hartree_bulk[i][0][:, :, j] * param.hartree_to_ev)
    if read_labels[1]:
        for j in range(average_charge_dft[i][0].shape[axis]):
            average_charge_dft[i][j] = np.mean(data_charge_dft[i][0][:, :, j])
            average_hartree_dft[i][j] = np.mean(data_hartree_dft[i][0][:, :, j] * param.hartree_to_ev)
    if read_labels[2]:
        for j in range(data_charge_em1[i][0].shape[axis]):
            average_charge_em1[i][j] = np.mean(data_charge_em1[i][0][:, :, j])
            average_hartree_em1[i][j] = np.mean(data_hartree_em1[i][0][:, :, j] * param.hartree_to_ev)
    if read_labels[3]:
        for j in range(data_charge_em2[i][0].shape[axis]):
            average_charge_em2[i][j] = np.mean(data_charge_em2[i][0][:, :, j])
            average_hartree_em2[i][j] = np.mean(data_hartree_em2[i][0][:, :, j] * param.hartree_to_ev)

# Setup axis
for i in range(len(folder_em2)):
    if read_labels[2]: print('size em1 folder', i,  data_charge_em1[i][1].get_cell()[axis][axis])
    if read_labels[3]: print('size em2 folder', i,  data_charge_em2[i][1].get_cell()[axis][axis])
    if read_labels[2]: print('size em1 folder', i, average_charge_em1[0].shape[0])
    if read_labels[3]: print('size em2 folder', i, average_charge_em2[0].shape[0])
if read_labels[3]: energy_grid_hartree_em2 = np.linspace(start=0, stop=data_charge_em2[0][1].get_cell()[axis][axis], num=average_charge_em2[0].shape[0])
if read_labels[2]: energy_grid_hartree_em1 = np.linspace(start=0, stop=data_charge_em1[0][1].get_cell()[axis][axis], num=average_charge_em1[0].shape[0])
if read_labels[1]: energy_grid_hartree_dft = np.linspace(start=0, stop=data_charge_dft[0][1].get_cell()[axis][axis], num=average_charge_dft[0].shape[0])
if read_labels[0]: energy_grid_hartree_bulk = np.linspace(start=0, stop=data_charge_bulk[0][1].get_cell()[axis][axis], num=average_charge_bulk[0].shape[0])
mid_index = np.argmin(abs(energy_grid_hartree_em1-mid_pos))
mid_pos_grid = energy_grid_hartree_em1[mid_index]

# Plot absolute
rows, cols = 2, 1
xlim = [0-1, data_charge_em1[0][1].get_cell()[axis][axis]+1]
# fig_cube_z, ax_cube_z = plt.subplots(rows, cols, sharex='col', sharey='row', figsize=(6, 8))
# if plot_fermi: ax_cube_z[0].axhline(y=fermi_dft, color='grey', linestyle='--', label='DFT Fermi energy', alpha=0.5)
# if draw_mirror: ax_cube_z[0].plot(energy_grid_hartree_em1[:mid_index]+mid_pos_grid, np.flip(average_hartree_em2[plot_folder][:mid_index]), 'm--')
# if draw_mirror: ax_cube_z[0].plot(energy_grid_hartree_em1[:mid_index]+mid_pos_grid, np.flip(average_hartree_em1[plot_folder][:mid_index]), 'm--')
# if plot_labels[3]: ax_cube_z[0].plot(energy_grid_hartree_em2, average_hartree_em2[0], '-', color=plot_color[2], label=labels[3])
# if plot_labels[3]: ax_cube_z[0].plot(energy_grid_hartree_em2, average_hartree_em2[1], '-', color=plot_color[3], label=labels[3])
# ax_cube_z[0].set_xlim([xlim[0], xlim[1]])
# ax_cube_z[0].legend(frameon=False)
# ax_cube_z[0].set_ylabel('Hartree potential / eV')
# if draw_mirror: ax_cube_z[1].plot(energy_grid_hartree_em1[:mid_index]+mid_pos_grid, np.flip(average_charge_em2[plot_folder][:mid_index]), 'm--')
# if draw_mirror: ax_cube_z[1].plot(energy_grid_hartree_em1[:mid_index]+mid_pos_grid, np.flip(average_charge_em1[plot_folder][:mid_index]), 'm--')
# if plot_labels[3]: ax_cube_z[1].plot(energy_grid_hartree_em2, average_charge_em2[0], '-', color=plot_color[2], label=labels[3])
# if plot_labels[3]: ax_cube_z[1].plot(energy_grid_hartree_em2, average_charge_em2[1], '-', color=plot_color[3], label=labels[3])
# if draw_markers: ax_cube_z[1].plot(markers, markers*0, 'o', color='orange', fillstyle='none')
# ax_cube_z[1].set_xlim([xlim[0], xlim[1]])
# ax_cube_z[1].legend(frameon=False)
# ax_cube_z[1].set_xlabel(r'Position / Å')
# ax_cube_z[1].set_ylabel('Charge density')
# fig_cube_z.tight_layout()
# fig_cube_z.savefig('{}/charge_hartree_cube_{}.png'.format(folder_em2[plot_folder], print_label), dpi=300)
# print('Finished plotting average')

# Plot Hartree and charge .cube difference
if plot_diff:
    fig_cube_both, ax_cube_both = plt.subplots(rows, cols, sharex='col', sharey='row', figsize=(6, 8))
    print(len(average_hartree_em2))
    print(average_hartree_em2)
    print(average_hartree_em1)
    for i in range(len(folder_em2)):
        ax_cube_both[0].plot(energy_grid_hartree_em1, (average_hartree_em2[i]-average_hartree_em1[i])/factor[i], '-', color=diff_color[i], label=diff_label[i])
        # if draw_mirror_diff: ax_cube_both[0].plot(energy_grid_hartree_em1[:mid_index] + mid_pos_grid, -1*np.flip(average_hartree_em2[i][:mid_index]-average_hartree_em1[i][:mid_index]), '--', color=diff_color[i])
    if plot_diff_legend[0]: ax_cube_both[0].legend(frameon=False)
    ax_cube_both[0].set_xlim([xlim[0], xlim[1]])
    # ax_cube_both[0].set_xlabel(r'Position / Å')
    ax_cube_both[0].set_ylabel('Hartree potential z / eV')
    for i in range(len(folder_em2)):
        ax_cube_both[1].plot(energy_grid_hartree_em1, (average_charge_em2[i]-average_charge_em1[i])/factor[i], '-', color=diff_color[i], label=diff_label[i])
        # if draw_mirror_diff: ax_cube_both[1].plot(energy_grid_hartree_em1[:mid_index] + mid_pos_grid, -1*np.flip(average_charge_em2[i][:mid_index]-average_charge_em1[i][:mid_index]), '--', color=diff_color[i])
    if plot_diff_legend[1]: ax_cube_both[1].legend(frameon=False)
    ax_cube_both[1].set_xlim([xlim[0], xlim[1]])
    ax_cube_both[1].set_xlabel(r'Position / Å')
    ax_cube_both[1].set_ylabel('Charge density z')
    fig_cube_both.tight_layout()
    for i in range(len(folder_em2)):
        fig_cube_both.savefig('{}/charge_hartree_cube_diff.png'.format(folder_em2[i]), dpi=300)
    if zoom:
        ax_cube_both[1].set_xlim([xlim2[0], xlim2[1]])
        fig_cube_both.tight_layout()
        for i in range(len(folder_em2)):
            fig_cube_both.savefig('{}/charge_hartree_cube_diff_zoom.png'.format(folder_em2[i]), dpi=300)
    print('Finished plotting difference Hartree and charge ')

# Plot Hartree
# fig_hartree, ax_hartree = plt.subplots()
# if plot_fermi: ax_hartree.axhline(y=fermi_dft, color='grey', linestyle='--', label='DFT Fermi energy', alpha=0.5)
# if plot_dft: ax_hartree.plot(energy_grid_hartree_dft, average_hartree_dft, 'r-', label=labels[0])
# ax_hartree.plot(energy_grid_hartree_em1, average_hartree_em1, 'g-', label=labels[1])
# if draw_mirror: ax_hartree.plot(energy_grid_hartree_em1[:mid_index]+mid_pos_grid, np.flip(average_hartree_em1[:mid_index]), 'm--')
# if plot_leads: ax_hartree.plot(energy_grid_hartree_bulk, average_hartree_bulk, 'k-', label=labels[2])
# ax_hartree.set_xlim([xlim[0], xlim[1]])
# ax_hartree.legend(frameon=False)
# ax_hartree.set_xlabel(r'Position / Å')
# ax_hartree.set_ylabel('Hartree potential z / eV')
# fig_hartree.tight_layout()
# fig_hartree.savefig('{}/hartree.png'.format(folder_save), dpi=300)

# Plot Hartree
# fig_hartree, ax_hartree = plt.subplots()
# if plot_diff:
#     for i in range(len(folder_em2)):
#         ax_hartree.plot(energy_grid_hartree_em1, average_hartree_em2[i] - average_hartree_em1[i], '-', color=diff_color[i], label=diff_label[i])
#         if draw_mirror_diff: ax_hartree.plot(energy_grid_hartree_em1[:mid_index] + mid_pos_grid, -1 * np.flip(average_hartree_em2[i][:mid_index] - average_hartree_em1[i][:mid_index]), '--', color=diff_color[i])
# ax_hartree.set_xlim([xlim[0], xlim[1]])
# ax_hartree.legend(frameon=False)
# ax_hartree.set_xlabel(r'Position / Å')
# ax_hartree.set_ylabel('Hartree potential z / eV')
# fig_hartree.tight_layout()
# for i in range(len(folder_em2)):
#     fig_hartree.savefig('{}/hartree_cube_diff_{}.png'.format(folder_em2[i], print_label), dpi=300)

# Plot charge
# fig_charge, ax_charge = plt.subplots()
# if plot_fermi: ax_charge.axhline(y=fermi_dft, color='grey', linestyle='--', label='DFT Fermi energy', alpha=0.5)
# if plot_dft: ax_charge.plot(energy_grid_charge_dft, average_charge_dft, 'r-', label=labels[0])
# ax_charge.plot(energy_grid_charge_em1, average_charge_em1, 'g-', label=labels[1])
# if draw_mirror: ax_charge.plot(energy_grid_charge_em1[:mid_index]+mid_pos_grid, np.flip(average_charge_em1[:mid_index]), 'm--')
# if plot_leads: ax_charge.plot(energy_grid_charge_bulk, average_charge_bulk, 'k-', label=labels[2])
# ax_charge.set_xlim([xlim[0], xlim[1]])
# ax_charge.legend(frameon=False)
# ax_charge.set_xlabel(r'Position / Å')
# ax_charge.set_ylabel('Charge density z / eV')
# fig_charge.tight_layout()
# fig_charge.savefig('{}/charge.png'.format(folder_save), dpi=300)

if __name__ == "__main__":
    print(folder_em2)
    print('Finished.')
    plt.show()
