import numpy as np
from matplotlib import pyplot as plt
from ase.io.cube import read_cube_data
from general import parameters as param

plot_fermi = False
plot_dft = True
use_xlim = False
draw_mirror = False
mirror_scale = 1
draw_markers = False
plot_dft_diff = True
fermi_dft = 0
labels = ['V=0.cube', 'V=1.cube', 'Bulk.cube']
print_label = 'V0_DFT'

# CP2K Li chain 1e density fixed HLB=F (Chris implementation)
# folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/iv/li/cp2k-smeagol/lda-q3_atoms-28/single-points/density-fix/chris/V-0_HLB-F_z-0-0'
# file_charge_em = '0V-ELECTRON_DENSITY-1_0.cube'
# file_hartree_em = '0V-v_hartree-1_0.cube'
# file_charge_dft = 'dft_wfn-ELECTRON_DENSITY-1_0.cube'
# file_hartree_dft = 'dft_wfn-v_hartree-1_0.cube'
# file_charge_bulk = 'bulk-ELECTRON_DENSITY-1_0.cube'
# file_hartree_bulk = 'bulk-v_hartree-1_0.cube'
# use_xlim = True
# xlim_specify_left = [-0.1, 6]
# xlim_specify_right = [78, 84.1]
# xlim_specify = xlim_specify_right

# # CP2K Li chain 1e density fixed HLB=T (Chris implementation)
# folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/iv/li/cp2k-smeagol/lda-q3_atoms-28/single-points/density-fix/chris/V-0_HLB-T-auto_z-0-0'
# folder_save = '{}/output'.format(folder)
# file_charge_em = '0V-ELECTRON_DENSITY-1_0.cube'
# file_hartree_em = '0V-v_hartree-1_0.cube'
# file_charge_dft = 'dft_wfn-ELECTRON_DENSITY-1_0.cube'
# file_hartree_dft = 'dft_wfn-v_hartree-1_0.cube'
# file_charge_bulk = 'bulk-ELECTRON_DENSITY-1_0.cube'
# file_hartree_bulk = 'bulk-v_hartree-1_0.cube'
# use_xlim = True
# xlim_specify_left = [-0.1, 6]
# xlim_specify_right = [78, 84.1]
# xlim_specify = xlim_specify_left

# CP2K Li chain 1e density fixed HLB=F (Sergey implementation)
# folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/iv/li/cp2k-smeagol/lda-q3_atoms-28/single-points/density-fix/sergey/V-0_HLB-F_z-0-0_sergey'
# file_charge_em = '0V-ELECTRON_DENSITY-1_0.cube'
# file_hartree_em = '0V-v_hartree-1_0.cube'
# file_charge_dft = 'dft_wfn-ELECTRON_DENSITY-1_0.cube'
# file_hartree_dft = 'dft_wfn-v_hartree-1_0.cube'
# file_charge_bulk = 'bulk-ELECTRON_DENSITY-1_0.cube'
# file_hartree_bulk = 'bulk-v_hartree-1_0.cube'
# use_xlim = True
# xlim_specify_left = [-0.1, 6]
# xlim_specify_right = [78, 84.1]
# xlim_specify = xlim_specify_right

# CP2K Li chain 1e density fixed HLB=T (Sergey implementation)
# folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/iv/li/cp2k-smeagol/lda-q3_atoms-28/single-points/density-fix/sergey/V-0_HLB-T-0.2022_z-0-0_sergey'
# file_charge_em = '0V-ELECTRON_DENSITY-1_0.cube'
# file_hartree_em = '0V-v_hartree-1_0.cube'
# file_charge_dft = 'dft_wfn-ELECTRON_DENSITY-1_0.cube'
# file_hartree_dft = 'dft_wfn-v_hartree-1_0.cube'
# file_charge_bulk = 'bulk-ELECTRON_DENSITY-1_0.cube'
# file_hartree_bulk = 'bulk-v_hartree-1_0.cube'
# use_xlim = True
# xlim_specify_left = [-0.1, 6]
# xlim_specify_right = [78, 84.1]
# xlim_specify = xlim_specify_right

# CP2K Au chain density fixed HLB=F (Chris implementation)
# folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/iv/au-chain/cp2k-smeagol/single-points/density-fix/chris/V-0_HLB-F_z-0-0'
# file_charge_em = '0V-ELECTRON_DENSITY-1_0.cube'
# file_hartree_em = '0V-v_hartree-1_0.cube'
# file_charge_dft = 'dft_wfn-ELECTRON_DENSITY-1_0.cube'
# file_hartree_dft = 'dft_wfn-v_hartree-1_0.cube'
# file_charge_bulk = 'bulk-ELECTRON_DENSITY-1_0.cube'
# file_hartree_bulk = 'bulk-v_hartree-1_0.cube'
# use_xlim = False
# xlim_specify_left = [-0.1, 6]
# xlim_specify_right = [78, 84.1]
# xlim_specify = xlim_specify_right

# CP2K Au chain density fixed HLB=T (Chris implementation)
# folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/iv/au-chain/cp2k-smeagol/single-points/density-fix/chris/V-0_HLB-F_z-0-0'
# file_charge_em = '0V-ELECTRON_DENSITY-1_0.cube'
# file_hartree_em = '0V-v_hartree-1_0.cube'
# file_charge_dft = 'dft_wfn-ELECTRON_DENSITY-1_0.cube'
# file_hartree_dft = 'dft_wfn-v_hartree-1_0.cube'
# file_charge_bulk = 'bulk-ELECTRON_DENSITY-1_0.cube'
# file_hartree_bulk = 'bulk-v_hartree-1_0.cube'
# use_xlim = False
# xlim_specify_left = [-0.1, 6]
# xlim_specify_right = [78, 84.1]
# xlim_specify = xlim_specify_right

# CP2K Au capacitor 1
# folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/au-capacitor/cp2k-smeagol/single-points/V-0_HLB-F_z-0-0_hirshfeld-gaussian_nodummy_hartree'
# file_charge_em = '0V-ELECTRON_DENSITY-1_0.cube'
# file_hartree_em = '0V-v_hartree-1_0.cube'
# file_charge_dft = 'dft_wfn-ELECTRON_DENSITY-1_0.cube'
# file_hartree_dft = 'dft_wfn-v_hartree-1_0.cube'
# file_charge_bulk = 'bulk-ELECTRON_DENSITY-1_0.cube'
# file_hartree_bulk = 'bulk-v_hartree-1_0.cube'
# fermi_dft = 0.20301667225796 * param.hartree_to_ev

# CP2K Au capacitor 2
# folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/au-capacitor/cp2k-smeagol/single-points/V-0_HLB-F_z-0-0_hirshfeld-gaussian_ordered_hartree'
# file_charge_em = '0V-ELECTRON_DENSITY-1_0.cube'
# file_hartree_em = '0V-v_hartree-1_0.cube'
# file_charge_dft = 'dft_wfn-ELECTRON_DENSITY-1_0.cube'
# file_hartree_dft = 'dft_wfn-v_hartree-1_0.cube'
# file_charge_bulk = 'bulk-ELECTRON_DENSITY-1_0.cube'
# file_hartree_bulk = 'bulk-v_hartree-1_0.cube'
# fermi_dft = 0.20359259837197 * param.hartree_to_ev

# CP2K Li chain 1e density fixed HLB=T (Chris implementation)
# folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/iv/li/cp2k-smeagol/lda-q3_atoms-28/single-points/density-fix/chris/V-0_HLB-T-auto_z-0-0'
# folder_save = '{}/output'.format(folder)
# file_charge_em = '0V-ELECTRON_DENSITY-1_0.cube'
# file_hartree_em = '0V-v_hartree-1_0.cube'
# file_charge_dft = 'dft_wfn-ELECTRON_DENSITY-1_0.cube'
# file_hartree_dft = 'dft_wfn-v_hartree-1_0.cube'
# file_charge_bulk = 'bulk-ELECTRON_DENSITY-1_0.cube'
# file_hartree_bulk = 'bulk-v_hartree-1_0.cube'
# use_xlim = True
# xlim_specify_left = [-0.1, 6]
# xlim_specify_right = [78, 84.1]
# xlim_specify = xlim_specify_left

# Au-BDT 001 CP2K
# folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/scarf/geometry-ordered/3x3-4/bulk_layers-4/kpoints_bulk-2-2-100_em-2-2-1_hlb-10.98783_scf-1'
# folder_save = '{}/output'.format(folder)
# file_charge_em = '0V-ELECTRON_DENSITY-1_0.cube'
# file_hartree_em = '0V-v_hartree-1_0.cube'
# file_charge_dft = 'dft_wfn-ELECTRON_DENSITY-1_0.cube'
# file_hartree_dft = 'dft_wfn-v_hartree-1_0.cube'
# file_charge_bulk = 'bulk-ELECTRON_DENSITY-1_0.cube'
# file_hartree_bulk = 'bulk-v_hartree-1_0.cube'
# use_xlim = False
# xlim_specify_left = [-0.1, 6]
# xlim_specify_right = [78, 84.1]
# xlim_specify = xlim_specify_left

# Li chain CP2K
# folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/li-chain/cp2k-smeagol/transmission/V-0_HLB-F_z-0-0_sergey'
# folder_save = '{}/output'.format(folder)
# file_charge_em = '0V-ELECTRON_DENSITY-1_0.cube'
# file_hartree_em = '0V-v_hartree-1_0.cube'
# file_charge_dft = 'dft_wfn-ELECTRON_DENSITY-1_0.cube'
# file_hartree_dft = 'dft_wfn-v_hartree-1_0.cube'
# file_charge_bulk = 'bulk-ELECTRON_DENSITY-1_0.cube'
# file_hartree_bulk = 'bulk-v_hartree-1_0.cube'
# use_xlim = False
# xlim_specify_left = [-0.1, 6]
# xlim_specify_right = [78, 84.1]
# xlim_specify = xlim_specify_left

# Au chain CP2K
# folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-chain/cp2k-smeagol/transmission/testing/chris/HLB-F_z-0-0'
# folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-chain/cp2k-smeagol/transmission/testing/chris/HLB-auto'
# folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-chain/cp2k-smeagol/transmission/testing/chris/HLB-auto_kpoints-1-1-20'
# folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-chain/cp2k-smeagol/transmission/testing/chris/HLB-auto_kpoints-4-4-20'
# folder_save = '{}/output'.format(folder)
# file_charge_em = '0V-ELECTRON_DENSITY-1_0.cube'
# file_hartree_em = '0V-v_hartree-1_0.cube'
# file_charge_em = 'V-ELECTRON_DENSITY-1_0.cube'
# file_hartree_em = 'V-v_hartree-1_0.cube'
# file_charge_dft = '0V-ELECTRON_DENSITY-1_0.cube'
# file_hartree_dft = '0V-v_hartree-1_0.cube'
# file_charge_dft = 'dft_wfn-ELECTRON_DENSITY-1_0.cube'
# file_hartree_dft = 'dft_wfn-v_hartree-1_0.cube'
# file_charge_bulk = 'bulk-ELECTRON_DENSITY-1_0.cube'
# file_hartree_bulk = 'bulk-v_hartree-1_0.cube'
# use_xlim = False
# plot_dft = True
# plot_dft_diff = True
# draw_mirror = False
# draw_markers = False

# Au capacitor testing
# folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/cp2k-smeagol/sz/transmission/exp/capacitor/sergey-equal/auto/kpoints-20_test/chris_kpoints-2-2-20_hlb-auto'
# folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/cp2k-smeagol/sz/transmission/exp/capacitor/sergey-2d-equal/layers-1/kpoints-4-4-20_hlb-auto'
# folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/cp2k-smeagol/sz/transmission/exp/capacitor/sergey-2d-equal/layers-1-2/kpoints-4-4-20_hlb-auto'
# folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/cp2k-smeagol/sz/transmission/exp/capacitor/sergey-2d-equal/layers-1-4/kpoints-4-4-20_hlb-auto'
# folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/cp2k-smeagol/sz/transmission/exp/capacitor/sergey-2d-equal/layers-1-2-3/kpoints-4-4-20_hlb-auto'
# folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/cp2k-smeagol/sz/transmission/exp/capacitor/sergey-2d-equal/layers-1-2-3-4/kpoints-1-1-20_hlb-auto'
# folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/cp2k-smeagol/sz/transmission/exp/capacitor/sergey-2d-equal/layers-1-2-3-4/kpoints-4-4-20_hlb-auto_vacuum'
# folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/cp2k-smeagol/sz/transmission/exp/capacitor/sergey-2d-equal/layers-1-2-3-4/NIMAGES/kpoints-4-4-20_hlb-auto_chris'
# folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/cp2k-smeagol/sz/transmission/exp/capacitor/sergey-equal/auto/kpoints-20_test/sergey_kpoints-2-2-20_hlb-auto-NIMAGES_IJ'
folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/cp2k-smeagol/sz/transmission/exp/capacitor/sergey-2d-equal/archer2/layers-1-2-3-4/kpoints-2-2-20_hlb-auto'
# folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/cp2k-smeagol/sz/transmission/exp/capacitor/sergey-2d-equal/archer2/layers-1-2-3-4/kpoints-4-4-20_hlb-auto-300K_CUTOFF-1500'
# folder_save = '{}/output'.format(folder)
# file_charge_em = '0V-ELECTRON_DENSITY-1_0.cube'
# file_hartree_em = '0V-v_hartree-1_0.cube'
# file_charge_em = 'V-ELECTRON_DENSITY-1_0.cube'
# file_hartree_em = 'V-v_hartree-1_0.cube'
# file_charge_dft = '0V-ELECTRON_DENSITY-1_0.cube'
# file_hartree_dft = '0V-v_hartree-1_0.cube'
# file_charge_dft = 'dft_wfn-ELECTRON_DENSITY-1_0.cube'
# file_hartree_dft = 'dft_wfn-v_hartree-1_0.cube'
# file_charge_bulk = 'bulk-ELECTRON_DENSITY-1_0.cube'
# file_hartree_bulk = 'bulk-v_hartree-1_0.cube'
# print_label = 'V0_DFT'
# labels = ['NEGF V=0', 'NEGF V=1', 'Leads']
# labels = ['DFT', 'NEGF', 'Leads']
# labels = ['DFT.cube', 'V-0.cube', 'Bulk.cube']
# mirror_scale = 1
# print_label = 'V1_DFT'
# labels = ['DFT.cube', 'V=1.cube', 'Bulk.cube']
# mirror_scale = 1
# print_label = 'V1_V0'
# labels = ['V=0.cube', 'V=1.cube', 'Bulk.cube']
# mirror_scale = -1
# use_xlim = False
# plot_dft = True
# plot_leads = False
# plot_dft_diff = True
# draw_mirror = False
# draw_markers = False

# Au capacitor
# folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/cp2k-smeagol/sz/transmission/exp/capacitor/sergey-equal/auto/kpoints-100/V-0_kpoints_bulk-4-4-20_em-4-4-1_hlb-auto'
# folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/cp2k-smeagol/sz/transmission/exp/capacitor/sergey-equal/auto/kpoints-20/V-1_kpoints_bulk-2-2-20_em-2-2-1_hlb-auto'
# folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/cp2k-smeagol/sz/transmission/exp/capacitor/sergey-equal/auto/kpoints-20_scf-1_atomic/V-1_kpoints_bulk-2-2-20_em-2-2-1_hlb-auto'
# folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/cp2k-smeagol/sz/transmission/exp/capacitor/sergey-equal/auto/kpoints-20_scf-1_restart-dft/V-0_kpoints_bulk-2-2-20_em-2-2-1_hlb-auto'
# # folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/cp2k-smeagol/sz/transmission/exp/capacitor/sergey-equal/auto/kpoints-20_scf-1_restart-0V/V-1_kpoints_bulk-4-4-20_em-4-4-1_hlb-auto'
# folder_save = '{}/output'.format(folder)
# file_charge_em = '0V-ELECTRON_DENSITY-1_0.cube'
# file_hartree_em = '0V-v_hartree-1_0.cube'
# file_charge_dft = 'dft_wfn-ELECTRON_DENSITY-1_0.cube'
# file_hartree_dft = 'dft_wfn-v_hartree-1_0.cube'
# file_charge_bulk = 'bulk-ELECTRON_DENSITY-1_0.cube'
# file_hartree_bulk = 'bulk-v_hartree-1_0.cube'
# use_xlim = False
# plot_dft = True
# plot_dft_diff = True
# draw_mirror = True
# draw_markers = False
# xlim_specify_left = [-0.1, 6]
# xlim_specify_right = [78, 84.1]
# xlim_specify = xlim_specify_left

# Au-BDT
# folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/cp2k-smeagol/sz/transmission/exp/sergey/kpoints_bulk-4-4-100_em-4-4-1_hlb-t-11.03197_scf-500'
# folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/cp2k-smeagol/sz/transmission/exp/dzvp/sergey/kpoints_bulk-4-4-100_em-4-4-1_hlb-t-auto'
# folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/cp2k-smeagol/sz/transmission/exp/capacitor/sergey-equal/auto/V-0_kpoints_bulk-2-2-100_em-2-2-1_hlb-auto'
# folder_save = '{}/output'.format(folder)
# file_charge_em = '0V-ELECTRON_DENSITY-1_0.cube'
# file_hartree_em = '0V-v_hartree-1_0.cube'
# file_charge_dft = 'dft_wfn-ELECTRON_DENSITY-1_0.cube'
# file_hartree_dft = 'dft_wfn-v_hartree-1_0.cube'
# file_charge_bulk = 'bulk-ELECTRON_DENSITY-1_0.cube'
# file_hartree_bulk = 'bulk-v_hartree-1_0.cube'
# use_xlim = False
# plot_dft = False
# xlim_specify_left = [-0.1, 6]
# xlim_specify_right = [78, 84.1]
# xlim_specify = xlim_specify_left

# Au-BDT Archer2
# folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/AuH2/transport/cp2k-smeagol/archer/transmission/kpoints-4-4-20_hlb-auto_sergey'
# folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/AuH2/transport/cp2k-smeagol/archer/transmission/kpoints-10-10-20_hlb-auto_sergey_nodes-4'
# folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/AuH2/transport/cp2k-smeagol/archer/transmission_capacitor/kpoints-5-5-20_hlb-auto_sergey_nodes-4'
# folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/AuH2/transport/cp2k-smeagol/archer/transmission_vacuum/kpoints-3-3-20_hlb-auto_sergey'
# folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/AuH2/transport/cp2k-smeagol/archer/transmission_capacitor/kpoints-5-5-20_hlb-auto_sergey_nodes-4'
# folder_save = '{}/output'.format(folder)
# # file_charge_em = '0V-ELECTRON_DENSITY-1_0.cube'
# # file_hartree_em = '0V-v_hartree-1_0.cube'
# file_charge_em = 'V-ELECTRON_DENSITY-1_0.cube'
# file_hartree_em = 'V-v_hartree-1_0.cube'
# file_charge_dft = '0V-ELECTRON_DENSITY-1_0.cube'
# file_hartree_dft = '0V-v_hartree-1_0.cube'
# # file_charge_dft = 'dft_wfn-ELECTRON_DENSITY-1_0.cube'
# # file_hartree_dft = 'dft_wfn-v_hartree-1_0.cube'
# file_charge_bulk = 'bulk-ELECTRON_DENSITY-1_0.cube'
# file_hartree_bulk = 'bulk-v_hartree-1_0.cube'
# # print_label = 'V0_DFT'
# # labels = ['DFT.cube', 'V-0.cube', 'Bulk.cube']
# # mirror_scale = 1
# # print_label = 'V1_DFT'
# # labels = ['DFT.cube', 'V=1.cube', 'Bulk.cube']
# # mirror_scale = 1
# print_label = 'V1_V0'
# labels = ['V=0.cube', 'V=1.cube', 'Bulk.cube']
# mirror_scale = -1
# use_xlim = False
# plot_dft = True
# plot_dft_diff = True
# draw_mirror = False
# draw_markers = False
# use_xlim = False
# plot_dft = False
# xlim_specify_left = [-0.1, 6]
# xlim_specify_right = [78, 84.1]
# xlim_specify = xlim_specify_left

# Read .cube using ASE
data_1, atoms_1 = read_cube_data('{}/{}'.format(folder, file_charge_em))
data_2, atoms_2 = read_cube_data('{}/{}'.format(folder, file_hartree_em))
if plot_dft or plot_dft_diff: data_3, atoms_3 = read_cube_data('{}/{}'.format(folder, file_charge_dft))
if plot_dft or plot_dft_diff: data_4, atoms_4 = read_cube_data('{}/{}'.format(folder, file_hartree_dft))
data_5, atoms_5 = read_cube_data('{}/{}'.format(folder, file_charge_bulk))
data_6, atoms_6 = read_cube_data('{}/{}'.format(folder, file_hartree_bulk))
print('Finished reading .cube files')

# Plot charge and hartree .cube x average
axis = 0
x_average_1 = np.zeros(data_1.shape[axis])
x_average_2 = np.zeros(data_2.shape[axis])
if plot_dft or plot_dft_diff: x_average_3 = np.zeros(data_3.shape[axis])
if plot_dft or plot_dft_diff: x_average_4 = np.zeros(data_4.shape[axis])
x_average_5 = np.zeros(data_5.shape[axis])
x_average_6 = np.zeros(data_6.shape[axis])
for i in range(data_1.shape[axis]):
    x_average_1[i] = np.mean(data_1[i, :, :])
    x_average_2[i] = np.mean(data_2[i, :, :] * param.hartree_to_ev)
    if plot_dft or plot_dft_diff: x_average_3[i] = np.mean(data_3[i, :, :])
    if plot_dft or plot_dft_diff: x_average_4[i] = np.mean(data_4[i, :, :] * param.hartree_to_ev)
for i in range(data_5.shape[axis]):
    x_average_5[i] = np.mean(data_5[i, :, :])
    x_average_6[i] = np.mean(data_6[i, :, :] * param.hartree_to_ev)
energx_grid_x_1 = np.linspace(start=0, stop=atoms_1.get_cell()[axis][axis], num=data_1.shape[axis])
energx_grid_x_2 = np.linspace(start=0, stop=atoms_1.get_cell()[axis][axis], num=data_2.shape[axis])
if plot_dft or plot_dft_diff: energx_grid_x_3 = np.linspace(start=0, stop=atoms_1.get_cell()[axis][axis], num=data_3.shape[axis])
if plot_dft or plot_dft_diff: energx_grid_x_4 = np.linspace(start=0, stop=atoms_1.get_cell()[axis][axis], num=data_4.shape[axis])
energx_grid_x_5 = np.linspace(start=0, stop=atoms_5.get_cell()[axis][axis], num=data_5.shape[axis])
energx_grid_x_6 = np.linspace(start=0, stop=atoms_5.get_cell()[axis][axis], num=data_6.shape[axis])
xlim = [0 - 1, atoms_1.get_cell()[axis][axis] + 1]

# rows, cols = 2, 1
# figcube_both, ax_cube_x = plt.subplots(rows, cols, sharex='col', sharey='row', figsize=(6, 8))
# if plot_dft: ax_cube_x[0].plot(energx_grid_x_4, x_average_4, 'r-', label='DFT .cube')
# ax_cube_x[0].plot(energx_grid_x_2, x_average_2, 'g-', label='EM .cube')
# ax_cube_x[0].plot(energx_grid_x_6, x_average_6, 'k-', label='Bulk .cube')
# ax_cube_x[0].set_xlim([xlim[0], xlim[1]])
# ax_cube_x[0].legend(frameon=False)
# ax_cube_x[0].set_xlabel(r'Position / Å')
# ax_cube_x[0].set_ylabel('Hartree potential x / eV')
# if plot_dft: ax_cube_x[1].plot(energx_grid_x_3, x_average_3, 'r-', label='DFT .cube')
# ax_cube_x[1].plot(energx_grid_x_1, x_average_1, 'g-', label='EM .cube')
# ax_cube_x[1].plot(energx_grid_x_5, x_average_5, 'k-', label='Bulk .cube')
# ax_cube_x[1].set_xlim([xlim[0], xlim[1]])
# ax_cube_x[1].legend(frameon=False)
# ax_cube_x[1].set_xlabel(r'Position / Å')
# ax_cube_x[1].set_ylabel('Charge density x')
# figcube_both.tight_layout()
# figcube_both.savefig('{}/charge_hartree_cube_x.png'.format(folder), dpi=300)
# print('Finished plotting x average')

# Plot charge and hartree .cube y average
axis = 1
y_average_1 = np.zeros(data_1.shape[axis])
y_average_2 = np.zeros(data_2.shape[axis])
if plot_dft or plot_dft_diff: y_average_3 = np.zeros(data_3.shape[axis])
if plot_dft or plot_dft_diff: y_average_4 = np.zeros(data_4.shape[axis])
y_average_5 = np.zeros(data_5.shape[axis])
y_average_6 = np.zeros(data_6.shape[axis])
for i in range(data_1.shape[axis]):
    y_average_1[i] = np.mean(data_1[:, i, :])
    y_average_2[i] = np.mean(data_2[:, i, :] * param.hartree_to_ev)
    if plot_dft or plot_dft_diff: y_average_3[i] = np.mean(data_3[:, i, :])
    if plot_dft or plot_dft_diff: y_average_4[i] = np.mean(data_4[:, i, :] * param.hartree_to_ev)
for i in range(data_5.shape[axis]):
    y_average_5[i] = np.mean(data_5[:, i, :])
    y_average_6[i] = np.mean(data_6[:, i, :] * param.hartree_to_ev)
energy_grid_y_1 = np.linspace(start=0, stop=atoms_1.get_cell()[axis][axis], num=data_1.shape[axis])
energy_grid_y_2 = np.linspace(start=0, stop=atoms_1.get_cell()[axis][axis], num=data_2.shape[axis])
if plot_dft or plot_dft_diff: energy_grid_y_3 = np.linspace(start=0, stop=atoms_1.get_cell()[axis][axis], num=data_3.shape[axis])
if plot_dft or plot_dft_diff: energy_grid_y_4 = np.linspace(start=0, stop=atoms_1.get_cell()[axis][axis], num=data_4.shape[axis])
energy_grid_y_5 = np.linspace(start=0, stop=atoms_5.get_cell()[axis][axis], num=data_5.shape[axis])
energy_grid_y_6 = np.linspace(start=0, stop=atoms_5.get_cell()[axis][axis], num=data_6.shape[axis])
xlim = [0 - 1, atoms_1.get_cell()[axis][axis] + 1]

# rows, cols = 2, 1
# figcube_both, ax_cube_y = plt.subplots(rows, cols, sharex='col', sharey='row', figsize=(6, 8))
# if plot_dft: ax_cube_y[0].plot(energy_grid_y_4, y_average_4, 'r-', label='DFT .cube')
# ax_cube_y[0].plot(energy_grid_y_2, y_average_2, 'g-', label='EM .cube')
# ax_cube_y[0].plot(energy_grid_y_6, y_average_6, 'k-', label='Bulk .cube')
# ax_cube_y[0].set_xlim([xlim[0], xlim[1]])
# ax_cube_y[0].legend(frameon=False)
# ax_cube_y[0].set_xlabel(r'Position / Å')
# ax_cube_y[0].set_ylabel('Hartree potential y / eV')
# if plot_dft: ax_cube_y[1].plot(energy_grid_y_3, y_average_3, 'r-', label='DFT .cube')
# ax_cube_y[1].plot(energy_grid_y_1, y_average_1, 'g-', label='EM .cube')
# ax_cube_y[1].plot(energy_grid_y_5, y_average_5, 'k-', label='Bulk .cube')
# ax_cube_y[1].set_xlim([xlim[0], xlim[1]])
# ax_cube_y[1].legend(frameon=False)
# ax_cube_y[1].set_xlabel(r'Position / Å')
# ax_cube_y[1].set_ylabel('Charge density y')
# figcube_both.tight_layout()
# figcube_both.savefig('{}/charge_hartree_cube_y.png'.format(folder), dpi=300)
# print('Finished plotting y average')

# Plot charge and hartree .cube z average
axis = 2
z_average_1 = np.zeros(data_1.shape[axis])
z_average_2 = np.zeros(data_2.shape[axis])
if plot_dft or plot_dft_diff: z_average_3 = np.zeros(data_3.shape[axis])
if plot_dft or plot_dft_diff: z_average_4 = np.zeros(data_4.shape[axis])
z_average_5 = np.zeros(data_5.shape[axis])
z_average_6 = np.zeros(data_6.shape[axis])
for i in range(data_1.shape[axis]):
    z_average_1[i] = np.mean(data_1[:, :, i])
    z_average_2[i] = np.mean(data_2[:, :, i] * param.hartree_to_ev)
    if plot_dft or plot_dft_diff: z_average_3[i] = np.mean(data_3[:, :, i])
    if plot_dft or plot_dft_diff: z_average_4[i] = np.mean(data_4[:, :, i] * param.hartree_to_ev) #+ (14.923291563435138-10.989786508975172)
for i in range(data_5.shape[axis]):
    z_average_5[i] = np.mean(data_5[:, :, i])
    z_average_6[i] = np.mean(data_6[:, :, i] * param.hartree_to_ev)

energy_grid_z_1 = np.linspace(start=0, stop=atoms_1.get_cell()[axis][axis], num=data_1.shape[axis])
energy_grid_z_2 = np.linspace(start=0, stop=atoms_1.get_cell()[axis][axis], num=data_2.shape[axis])
if plot_dft or plot_dft_diff: energy_grid_z_3 = np.linspace(start=0, stop=atoms_1.get_cell()[axis][axis], num=data_3.shape[axis])
if plot_dft or plot_dft_diff:energy_grid_z_4 = np.linspace(start=0, stop=atoms_1.get_cell()[axis][axis], num=data_4.shape[axis])
energy_grid_z_5 = np.linspace(start=0, stop=atoms_5.get_cell()[axis][axis], num=data_5.shape[axis])
energy_grid_z_6 = np.linspace(start=0, stop=atoms_5.get_cell()[axis][axis], num=data_6.shape[axis])
print('\nMin of Hartree potential z EM .cube', np.min(z_average_2))
if plot_dft or plot_dft_diff:print('Min of Hartree potential z DFT .cube', np.min(z_average_4))
print('Min of Hartree potential z bulk .cube', np.min(z_average_6))

print('\nValue(z=Min) of Hartree potential z EM .cube', z_average_2[0])
if plot_dft or plot_dft_diff: print('Value(z=Min) of Hartree potential z DFT .cube', z_average_4[0])
print('Value(z=Min) of Hartree potential z bulk .cube', z_average_6[0])

print('\nValue(z=Max) of Hartree potential z EM .cube', z_average_2[-1])
if plot_dft or plot_dft_diff: print('Value(z=Max) of Hartree potential z DFT .cube', z_average_4[-1])
print('Value(z=Max) of Hartree potential z bulk .cube', z_average_6[-1])
xlim = [0-1, atoms_1.get_cell()[axis][axis]+1]

mid_pos = 2.08400*5+(21.42200-2.08400*5)/2
mid_index = np.argmin(abs(energy_grid_z_2-mid_pos))
mid_pos_grid = energy_grid_z_2[mid_index]
print(mid_pos)
print(mid_pos_grid)

# Plot all
rows, cols = 2, 1
figcube_both, ax_cube_z = plt.subplots(rows, cols, sharex='col', sharey='row', figsize=(6, 8))
if plot_fermi: ax_cube_z[0].axhline(y=fermi_dft, color='grey', linestyle='--', label='DFT Fermi energy', alpha=0.5)
if plot_dft: ax_cube_z[0].plot(energy_grid_z_4, z_average_4, 'r-', label=labels[0])
ax_cube_z[0].plot(energy_grid_z_2, z_average_2, 'g-', label=labels[1])
if draw_mirror: ax_cube_z[0].plot(energy_grid_z_2[:mid_index]+mid_pos_grid, np.flip(z_average_2[:mid_index]), 'm--')
if plot_leads: ax_cube_z[0].plot(energy_grid_z_6, z_average_6, 'k-', label=labels[2])
ax_cube_z[0].set_xlim([xlim[0], xlim[1]])
ax_cube_z[0].legend(frameon=False)
ax_cube_z[0].set_ylabel('Hartree potential z / eV')
if plot_dft: ax_cube_z[1].plot(energy_grid_z_3, z_average_3, 'r-', label=labels[0])
ax_cube_z[1].plot(energy_grid_z_1, z_average_1, 'g-', label=labels[1])
if draw_mirror: ax_cube_z[1].plot(energy_grid_z_1[:mid_index]+mid_pos_grid, np.flip(z_average_1[:mid_index]), 'm--')
if plot_leads: ax_cube_z[1].plot(energy_grid_z_5, z_average_5, 'k-', label=labels[2])
ax_cube_z[1].set_xlim([xlim[0], xlim[1]])
ax_cube_z[1].legend(frameon=False)
ax_cube_z[1].set_xlabel(r'Position / Å')
ax_cube_z[1].set_ylabel('Charge density z')
figcube_both.tight_layout()
figcube_both.savefig('{}/charge_hartree_cube_z_{}.png'.format(folder, print_label), dpi=300)
print('Finished plotting z average')

# Plot with markers
if draw_markers:
    rows, cols = 2, 1
    figcube_both_markers, ax_cube_z_markers = plt.subplots(rows, cols, sharex='col', sharey='row', figsize=(6, 8))
    if plot_fermi: ax_cube_z_markers[0].axhline(y=fermi_dft, color='grey', linestyle='--', label='DFT Fermi energy', alpha=0.5)
    if plot_dft: ax_cube_z_markers[0].plot(energy_grid_z_4, z_average_4, 'r.-', label=labels[0])
    ax_cube_z_markers[0].plot(energy_grid_z_2, z_average_2, 'g.-', label=labels[1])
    if draw_mirror: ax_cube_z_markers[0].plot(energy_grid_z_2[:mid_index]+mid_pos_grid, np.flip(z_average_2[:mid_index]), 'm.--')
    ax_cube_z_markers[0].plot(energy_grid_z_6, z_average_6, 'k.-', label=labels[2])
    ax_cube_z_markers[0].set_xlim([xlim[0], xlim[1]])
    ax_cube_z_markers[0].legend(frameon=False)
    ax_cube_z_markers[0].set_ylabel('Hartree potential z / eV')
    if plot_dft: ax_cube_z_markers[1].plot(energy_grid_z_3, z_average_3, 'r.-', label=labels[0])
    ax_cube_z_markers[1].plot(energy_grid_z_1, z_average_1, 'g.-', label=labels[1])
    if draw_mirror: ax_cube_z_markers[1].plot(energy_grid_z_1[:mid_index]+mid_pos_grid, np.flip(z_average_1[:mid_index]), 'm.--')
    ax_cube_z_markers[1].plot(energy_grid_z_5, z_average_5, 'k.-', label=labels[2])
    ax_cube_z_markers[1].set_xlim([xlim[0], xlim[1]])
    ax_cube_z_markers[1].legend(frameon=False)
    ax_cube_z_markers[1].set_xlabel(r'Position / Å')
    ax_cube_z_markers[1].set_ylabel('Charge density z')
    figcube_both_markers.tight_layout()
    figcube_both_markers.savefig('{}/charge_hartree_cube_z_markers_{}.png'.format(folder, print_label), dpi=300)

# Plot Hartree .cube difference
# if plot_dft:
#     figcube_both, ax_cube_z = plt.subplots()
#     ax_cube_z.plot(energy_grid_z_2, z_average_2-z_average_4, 'k-', label=labels[0])
#     if draw_mirror: ax_cube_z.plot(energy_grid_z_2[:mid_index] + mid_pos_grid, np.flip(z_average_2[:mid_index])-np.flip(z_average_4[:mid_index]), 'm--')
#     ax_cube_z.set_xlim([xlim[0], xlim[1]])
#     # if use_xlim:ax_cube_z.set_xlim(xlim_specify)
#     ax_cube_z.set_xlabel(r'Position / Å')
#     ax_cube_z.set_ylabel('Hartree potential z / eV')
#     figcube_both.tight_layout()
#     figcube_both.savefig('{}/hartree_cube_z_diff.png'.format(folder), dpi=300)
#     print('Finished plotting z difference Hartree')

# Plot Hartree and charge .cube difference
if plot_dft_diff:
    figcube_both, ax_cube_z = plt.subplots(rows, cols, sharex='col', sharey='row', figsize=(6, 8))
    if plot_fermi: ax_cube_z[0].axhline(y=fermi_dft, color='grey', linestyle='--', label='DFT Fermi energy', alpha=0.5)
    ax_cube_z[0].plot(energy_grid_z_2, z_average_2-z_average_4, 'k-', label=labels[0])
    if draw_mirror: ax_cube_z[0].plot(energy_grid_z_2[:mid_index] + mid_pos_grid, (np.flip(z_average_2[:mid_index])-np.flip(z_average_4[:mid_index]))*mirror_scale, 'm--')
    ax_cube_z[0].set_xlim([xlim[0], xlim[1]])
    # ax_cube_z[0].legend(frameon=False)
    ax_cube_z[0].set_xlabel(r'Position / Å')
    ax_cube_z[0].set_ylabel('Hartree potential z / eV')
    ax_cube_z[1].plot(energy_grid_z_3, z_average_1-z_average_3, 'k-', label=labels[0])
    if draw_mirror: ax_cube_z[1].plot(energy_grid_z_2[:mid_index] + mid_pos_grid, (np.flip(z_average_1[:mid_index])-np.flip(z_average_3[:mid_index]))*mirror_scale, 'm--')
    ax_cube_z[1].set_xlim([xlim[0], xlim[1]])
    # ax_cube_z[1].legend(frameon=False)
    ax_cube_z[1].set_xlabel(r'Position / Å')
    ax_cube_z[1].set_ylabel('Charge density z')
    figcube_both.tight_layout()
    figcube_both.savefig('{}/charge_hartree_cube_z_diff_{}.png'.format(folder, print_label), dpi=300)
    print('Finished plotting z difference Hartree and charge ')

# Plot Hartree
fig_hartree, ax_hartree = plt.subplots()
if plot_fermi: ax_hartree.axhline(y=fermi_dft, color='grey', linestyle='--', label='DFT Fermi energy', alpha=0.5)
if plot_dft: ax_hartree.plot(energy_grid_z_4, z_average_4, 'r-', label=labels[0])
ax_hartree.plot(energy_grid_z_2, z_average_2, 'g-', label=labels[1])
if draw_mirror: ax_hartree.plot(energy_grid_z_2[:mid_index]+mid_pos_grid, np.flip(z_average_2[:mid_index]), 'm--')
if plot_leads: ax_hartree.plot(energy_grid_z_6, z_average_6, 'k-', label=labels[2])
ax_hartree.set_xlim([xlim[0], xlim[1]])
ax_hartree.legend(frameon=False)
ax_hartree.set_xlabel(r'Position / Å')
ax_hartree.set_ylabel('Hartree potential z / eV')
fig_hartree.tight_layout()
fig_hartree.savefig('{}/hartree.png'.format(folder), dpi=300)

# Plot charge
fig_charge, ax_charge = plt.subplots()
if plot_fermi: ax_charge.axhline(y=fermi_dft, color='grey', linestyle='--', label='DFT Fermi energy', alpha=0.5)
if plot_dft: ax_charge.plot(energy_grid_z_3, z_average_3, 'r-', label=labels[0])
ax_charge.plot(energy_grid_z_1, z_average_1, 'g-', label=labels[1])
if draw_mirror: ax_charge.plot(energy_grid_z_1[:mid_index]+mid_pos_grid, np.flip(z_average_1[:mid_index]), 'm--')
if plot_leads: ax_charge.plot(energy_grid_z_5, z_average_5, 'k-', label=labels[2])
ax_charge.set_xlim([xlim[0], xlim[1]])
ax_charge.legend(frameon=False)
ax_charge.set_xlabel(r'Position / Å')
ax_charge.set_ylabel('Charge density z / eV')
fig_charge.tight_layout()
fig_charge.savefig('{}/charge.png'.format(folder), dpi=300)

if __name__ == "__main__":
    print(folder)
    print('Finished.')
    plt.show()
