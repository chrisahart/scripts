import numpy as np
from matplotlib import pyplot as plt
from ase.io.cube import read_cube_data
from general import parameters as param

plot_fermi = False
fermi_dft = 0

# Au capacitor
# folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/cp2k-smeagol/sz/transmission/exp/capacitor/sergey-2d-equal/NEnergReal-640/layers-1-2-3-4/kpoints-2-2-20_hlb-auto'
# folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/cp2k-smeagol/sz/transmission/exp/capacitor/sergey-2d-equal/NEnergReal-640/layers-1-2-3-4/kpoints-4-4-20_hlb-auto'
# folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/cp2k-smeagol/sz/transmission/exp/capacitor/sergey-2d-equal/NEnergReal-80/layers-1-2-3-4/kpoints-4-4-20_hlb-auto_cores-120'
# folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/cp2k-smeagol/sz/transmission/exp/capacitor/sergey-2d-equal/NEnergReal-80-memory/layers-1-2-3-4/kpoints-4-4-20_hlb-auto_cores-120'
# folder_2 = folder_1
# file_charge_1 = '0V-ELECTRON_DENSITY-1_0.cube'
# file_hartree_1 = '0V-v_hartree-1_0.cube'
# file_charge_2 = 'V-ELECTRON_DENSITY-1_0.cube'
# file_hartree_2 = 'V-v_hartree-1_0.cube'
# plot_markers = False
# draw_mirror = False
# markers = np.array([2.08400*0, 2.08400*1, 2.08400*2, 2.08400*3, 2.08400*4, 2.08400*5, 21.42200+2.08400*0,
#                     21.42200+2.08400*1, 21.42200+2.08400*2, 21.42200+2.08400*3, 21.42200+2.08400*4, 21.42200+2.08400*5])
# use_xlim = False
# xlim_specify_left = [-0.1, 6]
# xlim_specify_right = [78, 84.1]
# xlim_specify = xlim_specify_left

# # Au Au-4TPA-C60-Au
# folder_2 = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/mengxuan/archer/archer/Au-4TPA-C60-Au/iv_parralel/kpoints-1-1-20_hlb-auto_cores-128_OrderN-F_finer/iv_curve/V_0.25'
# folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/mengxuan/archer/archer/Au-4TPA-C60-Au/iv_parralel/kpoints-1-1-20_hlb-auto_cores-128_OrderN-F_finer/iv_curve/V_0.00'
# file_charge_2 = 'V-0.25-ELECTRON_DENSITY-1_0.cube'
# file_hartree_2 = 'V-0.25-v_hartree-1_0.cube'
# # file_charge_2 = 'V--1.00-ELECTRON_DENSITY-1_0.cube'
# # file_hartree_2 = 'V--1.00-v_hartree-1_0.cube'
# file_charge_1 = 'V-0.00-ELECTRON_DENSITY-1_0.cube'
# file_hartree_1 = 'V-0.00-v_hartree-1_0.cube'
# plot_markers = False
# draw_mirror = False
# markers = np.array(
#     [2.08400 * 0, 2.08400 * 1, 2.08400 * 2, 2.08400 * 3, 2.08400 * 4, 2.08400 * 5, 21.42200 + 2.08400 * 0,
#      21.42200 + 2.08400 * 1, 21.42200 + 2.08400 * 2, 21.42200 + 2.08400 * 3, 21.42200 + 2.08400 * 4,
#      21.42200 + 2.08400 * 5])
# use_xlim = False
# xlim_specify_left = [-0.1, 6]
# xlim_specify_right = [78, 84.1]
# xlim_specify = xlim_specify_left

# Au ben_ant_chris
# folder_2 = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/mengxuan/archer/archer/ben_ant_chris/iv_parralel/kpoints-1-1-20_hlb-auto_cores-128_OrderN-F_finer/iv_curve/V_-1.00'
# folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/mengxuan/archer/archer/ben_ant_chris/iv_parralel/kpoints-1-1-20_hlb-auto_cores-128_OrderN-F_finer/iv_curve/V_0.00'
# file_charge_2 = 'V--1.00-ELECTRON_DENSITY-1_0.cube'
# file_hartree_2 = 'V--1.00-v_hartree-1_0.cube'
# file_charge_1 = 'V-0.00-ELECTRON_DENSITY-1_0.cube'
# file_hartree_1 = 'V-0.00-v_hartree-1_0.cube'
# plot_markers = False
# draw_mirror = False
# markers = np.array([2.08400*0, 2.08400*1, 2.08400*2, 2.08400*3, 2.08400*4, 2.08400*5, 21.42200+2.08400*0,
#                     21.42200+2.08400*1, 21.42200+2.08400*2, 21.42200+2.08400*3, 21.42200+2.08400*4, 21.42200+2.08400*5])
# use_xlim = False
# xlim_specify_left = [-0.1, 6]
# xlim_specify_right = [78, 84.1]
# xlim_specify = xlim_specify_left

# Au capacitor
folder_2 = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-capacitor/archer/layers-1-2-3-4/iv_parralel/kpoints-4-4-20_hlb-auto_vacuum/iv_curve/V_-1.0'
folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-capacitor/archer/layers-1-2-3-4/iv_parralel/kpoints-4-4-20_hlb-auto_vacuum/iv_curve/V_0.0'
# # folder_1 = folder_2
# # folder_2 = folder_1
# file_charge_2 = 'V-ELECTRON_DENSITY-1_0.cube'
# file_hartree_2 = 'V-v_hartree-1_0.cube'
# # file_charge_1 = 'dft_wfn-ELECTRON_DENSITY-1_0.cube'
# # file_hartree_1 = 'dft_wfn-v_hartree-1_0.cube'
# file_charge_1 = 'V-ELECTRON_DENSITY-1_0.cube'
# file_hartree_1 = 'V-v_hartree-1_0.cube'
# plot_markers = False
# draw_mirror = True
# mirror_factor = -1
# markers = np.array(
#     [2.08400 * 0, 2.08400 * 1, 2.08400 * 2, 2.08400 * 3, 2.08400 * 4, 2.08400 * 5, 21.42200 + 2.08400 * 0,
#      21.42200 + 2.08400 * 1, 21.42200 + 2.08400 * 2, 21.42200 + 2.08400 * 3, 21.42200 + 2.08400 * 4,
#      21.42200 + 2.08400 * 5])
# use_xlim = False
# xlim_specify_left = [-0.1, 6]
# xlim_specify_right = [78, 84.1]
# xlim_specify = xlim_specify_left

# folder_2 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/delete/old/au-q1-au-capacitor-kpoints-1-1-20-bulk-4layers-lattice-cu-10layers-NEnergReal-64-omp-2'
# folder_2 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/delete/2-2-10-layers/capacitor/testing/au-s-NIMAGES_IJ-7-EPS_PGF_ORB-1E-5-omp-1-kpoints-3-3_pbe-dft-rs-dc'
# folder_2 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/delete/2-2-10-layers/capacitor/testing/au-s-NIMAGES_IJ-7-EPS_PGF_ORB-1E-5-omp-1-kpoints-3-3_pbe-dft-rs-dc-symmetric'
# folder_2 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/delete/2-2-10-layers/capacitor/testing/au-s-NIMAGES_IJ-7-EPS_PGF_ORB-1E-5-omp-1-kpoints-4-4_pbe-dft-rs-dc'
# folder_2 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/delete/2-2-10-layers/capacitor/testing/au-s-NIMAGES_IJ-7-EPS_PGF_ORB-1E-5-omp-1-kpoints-3-3_pbe-dft-rs-dc-symmetric-au-1-dz-V-1'
# folder_2 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/supercell-1-1-5/capacitor/testing/NIMAGES_IJ-7-kpoints-3-3-V-0-screen-symmetric'
# folder_2 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/supercell-1-1-5/capacitor/testing/NIMAGES_IJ-7-kpoints-3-3-V-0-screen'
# folder_2 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/supercell-1-1-5/capacitor/testing/NIMAGES_IJ-7-kpoints-3-3-V-0'
# folder_2 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/supercell-1-1-5-bulk-6/junction/testing/NIMAGES_IJ-7-kpoints-2-2-SZV-PBE-CUSTOM-q1-f3-mpi-32-omp-4-V-0'
# folder_2 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/supercell-1-1-5-bulk-6-cu-1.86/capacitor/testing/NIMAGES_IJ-7-kpoints-2-2-V-0-SZV-PBE-CUSTOM-q1-f3-screen-symmetric'
# folder_2 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/supercell-1-1-5-bulk-6-cu-1.86/capacitor/testing/NIMAGES_IJ-7-kpoints-2-2-V-0-SZV-PBE-CUSTOM-q1-f3-screen'
# folder_2 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/supercell-1-1-5-bulk-6-cu-1.86/capacitor/testing/NIMAGES_IJ-7-kpoints-2-2-V-0-SZV-PBE-CUSTOM-q1-f3-V-0'
# folder_2 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/supercell-1-1-5-bulk-6-cu-1.86/capacitor/testing/NIMAGES_IJ-7-kpoints-2-2-V-0-SZV-PBE-CUSTOM-q1-f3'
# folder_2 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/supercell-1-1-5-bulk-6-cu-1.86/capacitor/testing/NIMAGES_IJ-7-kpoints-2-2-V-0-SZV-PBE-CUSTOM-q1-f3-all'
# folder_2 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/supercell-1-1-5-bulk-6-cu-1.86/capacitor/testing/NIMAGES_IJ-7-kpoints-2-2-V-0-SZV-PBE-CUSTOM-q1-f3-au'
# folder_2 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/supercell-1-1-5-bulk-6-cu-1.86/capacitor/testing/NIMAGES_IJ-7-kpoints-2-2-V-0-SZV-PBE-CUSTOM-q1-f3-au-symmetric'
# folder_2 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/test/kpoints-2-2-20_V-1_double-contour-rerun'
# folder_2 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/delete/2-2-10-layers/capacitor/testing/au-s-NIMAGES_IJ-7-EPS_PGF_ORB-1E-5-omp-1-kpoints-3-3_pbe-dft-rs-dc-symmetric-au-1-dz-V-1'
# folder_2 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/au/2-2-10-layers/capacitor/testing/au-s-NIMAGES_IJ-7-EPS_PGF_ORB-1E-5-omp-1-kpoints-3-3_pbe-dft-rs-dc-symmetric-au-1-dz-cu-tight'
folder_2 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/au/2-2-10-layers/junction/testing/au-s-NIMAGES_IJ-7-EPS_PGF_ORB-1E-5-omp-1-kpoints-2-2_pbe-dft-rs-V-1-1hr'
# folder_2 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/au/2-2-10-layers/capacitor/testing/au-s-NIMAGES_IJ-7-EPS_PGF_ORB-1E-5-omp-1-kpoints-3-3_pbe-dft-rs-dc-symmetric-au-1-dz-cu-1-1-5'
# folder_2 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/supercell-1-1-5-bulk-6-cu-1.86/capacitor/testing/NIMAGES_IJ-7-kpoints-3-3-V-0-SZV-PBE-CUSTOM-q1-f3'
# folder_2 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/supercell-1-1-5/capacitor/testing/NIMAGES_IJ-7-kpoints-3-3-V-0-test-symmetric'
# folder_2 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/supercell-1-1-5/capacitor/testing/NIMAGES_IJ-7-kpoints-3-3-V-0-test-symmetric-2-2-10'
# folder_2 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/supercell-1-1-5/capacitor/testing/NIMAGES_IJ-7-kpoints-3-3-V-0-test-symmetric-293K-3'
# folder_2 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/supercell-1-1-5-bulk-6/junction/testing/NIMAGES_IJ-7-kpoints-2-2-SZV-PBE-CUSTOM-q1-f3-mpi-32-omp-4-V-0'
# folder_2 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/supercell-1-1-5-bulk-6-cu-1.86/capacitor/testing/NIMAGES_IJ-7-kpoints-2-2-V-0-SZV-PBE-CUSTOM-q1-f3-V-1-293K'
# folder_2 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/supercell-1-1-5-bulk-6-cu-1.86/capacitor/bias/kpoints-3-3-symmetric-293K-600mgrid-6layers-V-0'
# folder_2 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/au/au-s-NIMAGES_IJ-7-EPS_PGF_ORB-1E-5-omp-1-kpoints-2-2_pbe-dft-rs-V-1-1hr-dc'
# folder_2 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/supercell-1-1-5-bulk-6-cu-1.86/capacitor/bias/kpoints-2-2-V-1-sc-WeightRho-0.5'
# folder_2 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/supercell-1-1-5-bulk-6-cu-1.86/capacitor/bias/kpoints-2-2-V-1-sc-WeightRho-0.5'

# folder_2 = "/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/supercell-1-1-5-bulk-6-cu-1.86/junction/bias/energy/kpoints-2-2-V-0-au-3pd-2-rs"
# folder_2 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/supercell-1-1-5-bulk-6-cu-1.86/junction/bias/energy/kpoints-2-2-V-0-au'
# folder_2 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/supercell-1-1-5-bulk-6-cu-1.86/junction/bias/energy/kpoints-2-2-V-0-rs-V-1-rs-sc'
# folder_2 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/supercell-1-1-5-bulk-6-cu-1.86/junction/bias/energy/kpoints-2-2-V-0-rs-V-1-rs-sc-WeightRho-1'

# folder_2 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/supercell-1-1-5-bulk-6-cu-1.86/capacitor/testing/NIMAGES_IJ-7-kpoints-2-2-V-0-SZV-PBE-CUSTOM-q1-f3-screen'
# folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/supercell-1-1-5-bulk-6-cu-1.86/capacitor/dft/NIMAGES_IJ-7-kpoints-2-2-SZV-PBE-CUSTOM-q1-f3-screen'
# folder_2 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/au-s-NIMAGES_IJ-7-EPS_PGF_ORB-1E-5-omp-1-kpoints-2-2_pbe-dft-rs-V-1-1hr-dc-test'

# folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/supercell-1-1-3-bulk-6-cu-1.86/junction/bias/kpoints-2-2-V-0'
# folder_2 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/supercell-1-1-3-bulk-6-cu-1.86/junction/bias/kpoints-2-2-V-0.1-dc'
# folder_2 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/supercell-1-1-3-bulk-6-cu-1.86/junction/bias/kpoints-2-2-V-1-dc'

# folder_2 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/supercell-1-1-3-bulk-6-cu-1.86/capacitor/bias/kpoints-2-2-dc-all-q1'
# folder_1 = folder_2

# folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/supercell-1-1-3-bulk-6-cu-1.86/junction/bias/kpoints-2-2-V-0-all-q1'
# folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/supercell-1-1-3-bulk-6-cu-1.86/junction/bias/kpoints-2-2-V-0'
# folder_2 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/supercell-1-1-3-bulk-6-cu-1.86/junction/bias/kpoints-2-2-V-0.1-dc-all-q1'
# folder_2 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/supercell-1-1-3-bulk-6-cu-1.86/junction/bias/kpoints-2-2-V-1-dc-all-q1'
# folder_2 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/supercell-1-1-3-bulk-6-cu-1.86/junction/bias/kpoints-2-2-V-1-dc-all-q1'
# folder_2 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/supercell-1-1-3-bulk-6-cu-1.86/junction/bias/kpoints-2-2-V-0.01-dc'
# folder_2 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/supercell-1-1-3-bulk-6-cu-1.86/junction/bias/kpoints-2-2-V-0.01-dc-all-q1'

folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/supercell-1-1-5-bulk-6-cu-1.86/junction/bias/energy/kpoints-2-2-V-0'
folder_2 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/supercell-1-1-5-bulk-6-cu-1.86/junction/bias/energy/kpoints-2-2-V-0.1'
folder_2 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/supercell-1-1-5-bulk-6-cu-1.86/junction/bias/energy/kpoints-2-2-V-0.1-sc-WeightRho-1'
# folder_2 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/supercell-1-1-5-bulk-6-cu-1.86/junction/bias/energy/kpoints-2-2-V-0.01'
# folder_2 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/supercell-1-1-5-bulk-6-cu-1.86/junction/bias/energy/kpoints-2-2-V-0.01-sc-WeightRho-1'
# folder_2 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/supercell-1-1-5-bulk-6-cu-1.86/junction/bias/energy/kpoints-2-2-V-0-rs-V-1-rs-sc-WeightRho-1'


# cu hafnia interface
folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/supercell-1-1-3-bulk-6-cu-1.86/junction/bias/kpoints-2-2-V-0-all-q1'
folder_2 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/supercell-1-1-3-bulk-6-cu-1.86/junction/bias/kpoints-2-2-V-0.1-dc-all-q1'
folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/supercell-1-1-3-bulk-6-cu-1.86/junction/bias/kpoints-2-2-V-0'
folder_2 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/supercell-1-1-3-bulk-6-cu-1.86/junction/bias/kpoints-2-2-V-1-dc'
folder_2 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/supercell-1-1-3-bulk-6-cu-1.86/junction/bias/kpoints-2-2-V-1-WeightRho-1'
folder_2 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/supercell-1-1-3-bulk-6-cu-1.86/junction/bias/kpoints-2-2-V-0.1-dc'
folder_2 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/supercell-1-1-3-bulk-6-cu-1.86/junction/bias/kpoints-2-2-V-0.1-WeightRho-1'
folder_2 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/supercell-1-1-3-bulk-6-cu-1.86/junction/bias/kpoints-2-2-V-0.01-dc'
folder_2 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/supercell-1-1-3-bulk-6-cu-1.86/junction/bias/kpoints-2-2-V-0.01-WeightRho-1'

folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/supercell-1-1-5-bulk-6-cu-1.86/junction/bias/energy/kpoints-2-2-V-0'
folder_2 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/supercell-1-1-5-bulk-6-cu-1.86/junction/bias/energy/kpoints-2-2-V-0-rs-V-1'
# folder_2 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/supercell-1-1-5-bulk-6-cu-1.86/junction/bias/energy/kpoints-2-2-V-0-rs-V-1-rs'
# folder_2 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/supercell-1-1-5-bulk-6-cu-1.86/junction/bias/energy/kpoints-2-2-V-0-rs-V-1-rs-square'
# folder_2 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/supercell-1-1-5-bulk-6-cu-1.86/junction/bias/energy/kpoints-2-2-V-0-rs-V-1-rs-mulliken'

# cu wire
# folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/chain/centred/capacitor-q1-f3-6atoms'
# folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/chain/centred/wire-q1-f3-broken-left-square-6atoms'
folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/chain/centred/capacitor-q1-f3-6atoms'
folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/chain/centred/wire-q1-f3-broken-left-square-6atoms'
folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/chain/centred/wire-q1-f3-broken-left-sc-WeightRho-0-6atoms'
folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/chain/centred/wire-q1-f3-broken-left-sc-WeightRho-1-6atoms'
folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/chain/centred/wire-q1-f3-broken-middle-square-6atoms'
folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/chain/centred/wire-q1-f3-broken-left-right-square-6atoms'
folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/chain/centred/wire-q1-f3-broken-left-right-sc-WeightRho-0-6atoms-cp2k2024'
# folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/chain/centred/wire-q1-f3-broken-middle-sc-WeightRho-0-6atoms-cp2k2024'
# folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/chain/centred/wire-q1-f3-broken-middle-sc-WeightRho-1-6atoms-cp2k2024'
# folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/chain/centred/wire-q1-f3-broken-middle-q11-square-6atoms'
# folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/chain/centred/wire-q1-f3-broken-middle-q11-sc-WeightRho-0-6atoms-cp2k2024'
folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/chain/centred/wire-q1-f3-one-cu-removed-middle'
folder_2 = folder_1

#
# folder_1 = folder_2
# file_charge_2 = '4V-ELECTRON_DENSITY-1_0.cube'
# file_hartree_2 = '4V-v_hartree-1_0.cube'
file_charge_2 = '1V-ELECTRON_DENSITY-1_0.cube'
file_hartree_2 = '1V-v_hartree-1_0.cube'
# file_charge_2 = '0.1V-ELECTRON_DENSITY-1_0.cube'
# file_hartree_2 = '0.1V-v_hartree-1_0.cube'
# file_charge_2 = '0.01V-ELECTRON_DENSITY-1_0.cube'
# file_hartree_2 = '0.01V-v_hartree-1_0.cube'
# file_charge_2 = '-0.1V-ELECTRON_DENSITY-1_0.cube'
# file_hartree_2 = '-0.1V-v_hartree-1_0.cube'
# file_charge_2 = '0V-ELECTRON_DENSITY-1_0.cube'
# file_hartree_2 = '0V-v_hartree-1_0.cube'
file_charge_1 = '0V-ELECTRON_DENSITY-1_0.cube'
file_hartree_1 = '0V-v_hartree-1_0.cube'
# file_charge_1 = 'dft_wfn-ELECTRON_DENSITY-1_0.cube'
# file_hartree_1 = 'dft_wfn-v_hartree-1_0.cube'
plot_markers = False
draw_mirror = False
# draw_mirror = False
mirror_factor = -1
marker_z = 1.80500
marker_left = np.array(
    [marker_z * 0, marker_z * 1, marker_z * 2, marker_z * 3, marker_z * 4, marker_z * 5, marker_z * 6, marker_z * 7,
     marker_z * 8, marker_z * 9])
marker_right = np.array(
    [marker_z * 0, marker_z * 1, marker_z * 2, marker_z * 3, marker_z * 4, marker_z * 5, marker_z * 6, marker_z * 7,
     marker_z * 8, marker_z * 9]) + 47.23200
markers = np.concatenate([marker_left, marker_right])
use_xlim = False
xlim_specify_left = [-0.1, 6]
xlim_specify_right = [78, 84.1]
xlim_specify = xlim_specify_left

# Read .cube using ASE
print('Start reading .cube files')
data_1, atoms_1 = read_cube_data('{}/{}'.format(folder_1, file_charge_1))
data_2, atoms_2 = read_cube_data('{}/{}'.format(folder_1, file_hartree_1))
data_3, atoms_3 = read_cube_data('{}/{}'.format(folder_2, file_charge_2))
data_4, atoms_4 = read_cube_data('{}/{}'.format(folder_2, file_hartree_2))
print('Finished reading .cube files')

# Plot charge and hartree .cube z average
axis = 2
z_average_1 = np.zeros(data_1.shape[axis])
z_average_2 = np.zeros(data_2.shape[axis])
z_average_3 = np.zeros(data_3.shape[axis])
z_average_4 = np.zeros(data_4.shape[axis])
for i in range(data_1.shape[axis]):
    z_average_1[i] = np.mean(data_1[:, :, i])
    z_average_2[i] = np.mean(data_2[:, :, i] * param.hartree_to_ev)
    z_average_3[i] = np.mean(data_3[:, :, i])
    z_average_4[i] = np.mean(data_4[:, :, i] * param.hartree_to_ev)

energy_grid_z_1 = np.linspace(start=0, stop=atoms_1.get_cell()[axis][axis], num=data_1.shape[axis])
energy_grid_z_2 = np.linspace(start=0, stop=atoms_1.get_cell()[axis][axis], num=data_2.shape[axis])
energy_grid_z_3 = np.linspace(start=0, stop=atoms_1.get_cell()[axis][axis], num=data_3.shape[axis])
energy_grid_z_4 = np.linspace(start=0, stop=atoms_1.get_cell()[axis][axis], num=data_4.shape[axis])
xlim = [0-1, atoms_1.get_cell()[axis][axis]+1]
# mid_pos = 18.04900  +(37.04600 -18.04900)/2
mid_pos = 16.24500  +(47.23200 -16.24500)/2
mid_pos = 14.43900 + (57.81600-14.43900)/2

# symmetric
# mid_pos = 19.85400 + (50.59700-19.85400)/2

# non-symmetric
# mid_pos = 21.65900 + (50.59700-21.65900)/2

# supercell-1-1-3-bulk-6-cu-1.86 non-symmetric
mid_pos = 21.66000 + (40.49800-21.66000)/2

# mid_pos = 65.282/2 # 2.08400*5+(21.42200-2.08400*5)/2
mid_index = np.argmin(abs(energy_grid_z_2-mid_pos))
mid_pos_grid = energy_grid_z_2[mid_index]
print('mid_pos', mid_pos)
# print(mid_pos_grid)

rows, cols = 2, 1
figcube_both, ax_cube_z = plt.subplots(rows, cols, sharex='col', sharey='row', figsize=(6, 8))
ax_cube_z[0].plot(energy_grid_z_4, z_average_4-z_average_2, 'k-', label='1V - 0 V')
if plot_markers: ax_cube_z[0].plot(markers, markers*0, 'o', color='orange', fillstyle='none')
if draw_mirror: ax_cube_z[0].plot(energy_grid_z_2[:mid_index] + mid_pos_grid, mirror_factor*(np.flip(z_average_4[:mid_index]) - np.flip(z_average_2[:mid_index])), 'm--', label='EM .cube mirror')
ax_cube_z[0].set_xlim([xlim[0], xlim[1]])
# ax_cube_z[0].legend(frameon=False)
# ax_cube_z[0].set_xlabel(r'Position / Å')
ax_cube_z[0].set_ylabel('Hartree potential z / eV')
ax_cube_z[1].plot(energy_grid_z_3, z_average_3-z_average_1, 'k-', label='1V - 0 V')
if plot_markers: ax_cube_z[1].plot(markers, markers*0, 'o', color='orange', fillstyle='none')
if draw_mirror: ax_cube_z[1].plot(energy_grid_z_2[:mid_index] + mid_pos_grid, mirror_factor*(np.flip(z_average_3[:mid_index]) - np.flip(z_average_1[:mid_index])), 'm--', label='EM .cube mirror')
ax_cube_z[1].set_xlim([xlim[0], xlim[1]])
# ax_cube_z[1].legend(frameon=False)
ax_cube_z[1].set_xlabel(r'Position / Å')
ax_cube_z[1].set_ylabel('Charge density z')
figcube_both.tight_layout()
figcube_both.savefig('{}/charge_hartree_cube_z_diff.png'.format(folder_1), dpi=300)
figcube_both.savefig('{}/charge_hartree_cube_z_diff.png'.format(folder_2), dpi=300)
print('Finished plotting z average')

if __name__ == "__main__":
    print('folder_1', folder_1)
    print('folder_1', folder_2)
    print('Finished.')
    plt.show()
