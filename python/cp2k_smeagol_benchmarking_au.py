import numpy as np
from matplotlib import pyplot as plt
from ase.io.cube import read_cube_data
from general import parameters as param
import pandas as pd
from scipy.optimize import curve_fit


def f_1(x, A, B):
    return A * x ** 1 + B


def f_2(x, A, B):
    return A * x ** 2 + B


def f_3(x, A, B):
    return A * x ** 3 + B


siesta_scale = 1
cp2k_scale = 1
scale_cores_siesta = False
scale_cores_cp2k = False
plot_dft = True
labels_cp2k_1 = ['CP2K', 'CP2K-SMEAGOL V=0', 'CP2K-SMEAGOL V=1']
labels_cp2k_2 = ['CP2K 2', 'CP2K-SMEAGOL 2 V=0', 'CP2K-SMEAGOL 2 V=1']
labels_siesta = ['SIESTA', 'SIESTA-SMEAGOL V=0', 'SIESTA-SMEAGOL V=1']
sub_folders = ['layers-1', 'layers-1-2', 'layers-1-2-3', 'layers-1-2-3-4']
num_kpoints = 4
kpoints = num_kpoints * 0
size_cp2k = [24, 24*2, 24*3, 24*4]
size_siesta = [24, 24*2, 24*3, 24*4]

# Compare CP2K and SIESTA Au capacitor
# folder_cp2k = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/cp2k-smeagol/sz/transmission/exp/capacitor/sergey-2d-equal/benchmarking/NEnergReal-640'
# folder_siesta = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/siesta-smeagol/capacitor/benchmarking/NEnergReal-640'
# cols_cp2k = ['folder', 'File', 'SCF', 'A', 'B', 'C', 'D', 'Time']
# cols_siesta = ['folder', 'File', 'A', 'SCF', 'B', 'Time', 'C']
# cores = 1
# scale_cores_siesta = True
# scale_cores_cp2k = False

# Compare CP2K and SIESTA Au capacitor
# folder_cp2k = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/cp2k-smeagol/sz/transmission/exp/capacitor/sergey-2d-equal/benchmarking/NEnergReal-80'
# folder_siesta = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/siesta-smeagol/capacitor/benchmarking/NEnergReal-80_OrderN'
# cols_cp2k = ['folder', 'File', 'SCF', 'A', 'B', 'C', 'D', 'Time']
# cols_siesta = ['folder', 'File', 'A', 'SCF', 'B', 'Time', 'C']
# cores = 1
# scale_cores_siesta = True
# scale_cores_cp2k = False
# siesta_offset = -1
# axis_limits = True
# ylim_1 = [0, 40]
# ylim_2 = [0, 135]
# ylim_3 = [0, 70]

# Compare CP2K EPS_PGF_ORB
# folder_cp2k = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/cp2k-smeagol/sz/transmission/exp/capacitor/sergey-2d-equal/benchmarking/NEnergReal-80_EPS_PGF_ORB-0'
# folder_siesta = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/cp2k-smeagol/sz/transmission/exp/capacitor/sergey-2d-equal/benchmarking/NEnergReal-80_EPS_PGF_ORB-1e-6'
# folder_siesta = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/cp2k-smeagol/sz/transmission/exp/capacitor/sergey-2d-equal/benchmarking/NEnergReal-80_EPS_PGF_ORB-1e-12'
# cols_cp2k = ['folder', 'File', 'SCF', 'A', 'B', 'C', 'D', 'Time']
# cols_siesta = ['folder', 'File', 'SCF', 'A', 'B', 'C', 'D', 'Time']
# labels_cp2k_1 = ['CP2K EPS_PGF_ORB=0', 'CP2K-SMEAGOL V=0 EPS_PGF_ORB=0', 'CP2K-SMEAGOL V=1 EPS_PGF_ORB=0']
# labels_siesta = ['CP2K EPS_PGF_ORB=1e-6', 'CP2K-SMEAGOL V=0 EPS_PGF_ORB=1e-6', 'CP2K-SMEAGOL V=1 EPS_PGF_ORB=1e-6']
# labels_siesta = ['CP2K EPS_PGF_ORB=1e-12', 'CP2K-SMEAGOL V=0 EPS_PGF_ORB=1e-12', 'CP2K-SMEAGOL V=1 EPS_PGF_ORB=1e-12']
# cores = 1
# scale_cores_siesta = False
# scale_cores_cp2k = False
# plot_dft = False
# axis_limits = True
# ylim_1 = [0, 40]
# ylim_2 = [0, 135]
# ylim_3 = [0, 70]
# siesta_offset = 0

# Compare CP2K smeagol-interface-cp2k-20220907 and CP2K smeagol-interface-cp2k-20220907-memory
# folder_cp2k = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/cp2k-smeagol/sz/transmission/exp/capacitor/sergey-2d-equal/benchmarking/NEnergReal-80'
# folder_siesta = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/cp2k-smeagol/sz/transmission/exp/capacitor/sergey-2d-equal/benchmarking/NEnergReal-80-memory'
# cols_cp2k = ['folder', 'File', 'SCF', 'A', 'B', 'C', 'D', 'Time']
# cols_siesta = ['folder', 'File', 'SCF', 'A', 'B', 'C', 'D', 'Time']
# labels_cp2k_1 = ['CP2K 20220907', 'CP2K-SMEAGOL V=0 20220907', 'CP2K-SMEAGOL V=1 20220907']
# labels_siesta = ['CP2K 20220907-memory', 'CP2K-SMEAGOL V=0 20220907-memory', 'CP2K-SMEAGOL V=1 20220907-memory']
# cores = 1
# scale_cores_siesta = False
# scale_cores_cp2k = False
# plot_dft = False
# axis_limits = True
# ylim_1 = [0, 40]
# ylim_2 = [0, 135]
# ylim_3 = [0, 70]
# siesta_offset = 0

# Compare CP2K smeagol-interface-cp2k-20220907-memory OrderN T and SIESTA EM.OrderN T
# folder_cp2k = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/cp2k-smeagol/sz/transmission/exp/capacitor/sergey-2d-equal/benchmarking/NEnergReal-80-memory'
# folder_siesta = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/siesta-smeagol/capacitor/benchmarking/NEnergReal-80_OrderN'
# cols_cp2k = ['folder', 'File', 'SCF', 'A', 'B', 'C', 'D', 'Time']
# cols_siesta = ['folder', 'File', 'A', 'SCF', 'B', 'Time', 'C']
# labels_cp2k_1 = ['CP2K EM.OrderN T', 'CP2K-SMEAGOL V=0 EM.OrderN T', 'CP2K-SMEAGOL V=1 EM.OrderN T']
# labels_siesta = ['SIESTA EM.OrderN T', 'SIESTA-SMEAGOL V=0 EM.OrderN T', 'SIESTA-SMEAGOL V=1 EM.OrderN T']
# cores = 1
# scale_cores_siesta = True
# scale_cores_cp2k = False
# plot_dft = True
# axis_limits = True
# ylim_1 = [0, 40]
# ylim_2 = [0, 135]
# ylim_3 = [0, 70]
# siesta_offset = -1

# # Compare CP2K smeagol-interface-cp2k-20220907-memory OrderN F and SIESTA EM.OrderN F
# folder_cp2k = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/cp2k-smeagol/sz/transmission/exp/capacitor/sergey-2d-equal/benchmarking/NEnergReal-80-memory_OrderN-F'
# folder_siesta = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/siesta-smeagol/capacitor/benchmarking/NEnergReal-80'
# cols_cp2k = ['folder', 'File', 'SCF', 'A', 'B', 'C', 'D', 'Time']
# cols_siesta = ['folder', 'File', 'A', 'SCF', 'B', 'Time', 'C']
# labels_cp2k_1 = ['CP2K EM.OrderN F', 'CP2K-SMEAGOL V=0 EM.OrderN F', 'CP2K-SMEAGOL V=1 EM.OrderN F']
# labels_siesta = ['SIESTA EM.OrderN F', 'SIESTA-SMEAGOL V=0 EM.OrderN F', 'SIESTA-SMEAGOL V=1 EM.OrderN F']
# cores = 1
# scale_cores_siesta = True
# scale_cores_cp2k = False
# plot_dft = True
# axis_limits = True
# ylim_1 = [0, 40]
# ylim_2 = [0, 135]
# ylim_3 = [0, 70]
# siesta_offset = -1

# Compare CP2K NEnergReal-80-memory and NEnergReal-80-memory_OutputInfo-F
# folder_cp2k = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/cp2k-smeagol/sz/transmission/exp/capacitor/sergey-2d-equal/benchmarking/NEnergReal-80-memory'
# folder_siesta = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/cp2k-smeagol/sz/transmission/exp/capacitor/sergey-2d-equal/benchmarking/NEnergReal-80-memory_OutputInfo-F'
# cols_cp2k = ['folder', 'File', 'SCF', 'A', 'B', 'C', 'D', 'Time']
# cols_siesta = ['folder', 'File', 'SCF', 'A', 'B', 'C', 'D', 'Time']
# labels_cp2k_1 = ['CP2K OutputInfo T', 'CP2K-SMEAGOL V=0 OutputInfo T', 'CP2K-SMEAGOL V=1 OutputInfo T']
# labels_siesta = ['CP2K OutputInfo F', 'CP2K-SMEAGOL V=0 OutputInfo F', 'CP2K-SMEAGOL V=1 OutputInfo F']
# cores = 1
# scale_cores_siesta = False
# scale_cores_cp2k = False
# plot_dft = False
# axis_limits = True
# ylim_1 = [0, 40]
# ylim_2 = [0, 135]
# ylim_3 = [0, 70]
# siesta_offset = 0

# Compare CP2K NEnergReal-80-memory and NEnergReal-80-memory_OutputInfo-F
# folder_cp2k = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/cp2k-smeagol/sz/transmission/exp/capacitor/sergey-2d-equal/benchmarking/NEnergReal-80-memory'
# folder_siesta = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/cp2k-smeagol/sz/transmission/exp/capacitor/sergey-2d-equal/benchmarking/NEnergReal-80-memory_ParallelOverKNum-auto'
# cols_cp2k = ['folder', 'File', 'SCF', 'A', 'B', 'C', 'D', 'Time']
# cols_siesta = ['folder', 'File', 'SCF', 'A', 'B', 'C', 'D', 'Time']
# labels_cp2k_1 = ['CP2K ParallelOverKNum 1', 'CP2K-SMEAGOL V=0 ParallelOverKNum 1', 'CP2K-SMEAGOL V=1 ParallelOverKNum 1']
# labels_siesta = ['CP2K ParallelOverKNum -1', 'CP2K-SMEAGOL V=0 ParallelOverKNum -1', 'CP2K-SMEAGOL V=1 ParallelOverKNum -1']
# cores = 1
# scale_cores_siesta = False
# scale_cores_cp2k = False
# plot_dft = False
# axis_limits = True
# ylim_1 = [0, 40]
# ylim_2 = [0, 135]
# ylim_3 = [0, 70]
# siesta_offset = 0

# Compare CP2K NEnergReal-80-memory_NTransmPoints-0 and  NEnergReal-80-memory_OutputInfo-F_NTransmPoints-0_printing-off
# folder_cp2k = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/cp2k-smeagol/sz/transmission/exp/capacitor/sergey-2d-equal/benchmarking/NEnergReal-80-memory'
# folder_siesta = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/cp2k-smeagol/sz/transmission/exp/capacitor/sergey-2d-equal/benchmarking/NEnergReal-80-memory_OutputInfo-F_NTransmPoints-0_printing-off'
# cols_cp2k = ['folder', 'File', 'SCF', 'A', 'B', 'C', 'D', 'Time']
# cols_siesta = ['folder', 'SCF', 'Time']
# labels_cp2k_1 = ['CP2K debug T', 'CP2K-SMEAGOL V=0 debug T', 'CP2K-SMEAGOL V=1 debug T']
# labels_siesta = ['CP2K debug F', 'CP2K-SMEAGOL V=0 debug F', 'CP2K-SMEAGOL V=1 debug F']
# cores = 1
# scale_cores_siesta = False
# scale_cores_cp2k = False
# plot_dft = False
# axis_limits = True
# ylim_1 = [0, 40]
# ylim_2 = [0, 135]
# ylim_3 = [0, 70]
# siesta_offset = 0

# Compare CP2K and CP2K with shorter Au basis set
# folder_cp2k = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/cp2k-smeagol/sz/transmission/exp/capacitor/sergey-2d-equal/benchmarking/NEnergReal-80'
# folder_siesta = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/cp2k-smeagol/sz/transmission/exp/capacitor/sergey-2d-equal/benchmarking/NEnergReal-80_Au-SR'
# cols_cp2k = ['folder', 'File', 'SCF', 'A', 'B', 'C', 'D', 'Time']
# cols_siesta = ['folder', 'File', 'SCF', 'A', 'B', 'C', 'D', 'Time']
# labels_cp2k_1 = ['CP2K ζ=2.95', 'CP2K-SMEAGOL V=0 ζ=2.95', 'CP2K-SMEAGOL V=1 ζ=2.95']
# labels_siesta = ['CP2K ζ=1.95', 'CP2K-SMEAGOL V=0 ζ=1.95', 'CP2K-SMEAGOL V=1 ζ=1.95']
# cores = 1
# scale_cores_siesta = False
# scale_cores_cp2k = False

# Compare SIESTA EM.OrderN F and EM.OrderN T
# folder_cp2k = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/siesta-smeagol/capacitor/benchmarking/NEnergReal-80'
# folder_siesta = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/siesta-smeagol/capacitor/benchmarking/NEnergReal-80_OrderN'
# cols_cp2k = ['folder', 'File', 'A', 'SCF', 'B', 'Time', 'C']
# cols_siesta = ['folder', 'File', 'A', 'SCF', 'B', 'Time', 'C']
# labels_cp2k_1 = ['SIESTA EM.OrderN F', 'SIESTA-SMEAGOL V=0 EM.OrderN F', 'SIESTA-SMEAGOL V=1 EM.OrderN F']
# labels_siesta = ['SIESTA EM.OrderN T', 'SIESTA-SMEAGOL V=0 EM.OrderN T', 'SIESTA-SMEAGOL V=1 EM.OrderN T']
# cores = 1
# scale_cores_siesta = True
# scale_cores_cp2k = True

# Compare SIESTA NSlices 1 and NSlices 5
# folder_cp2k = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/siesta-smeagol/capacitor/benchmarking/NEnergReal-80_OrderN'
# folder_siesta = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/siesta-smeagol/capacitor/benchmarking/NEnergReal-80_OrderN_NSlices-5'
# cols_cp2k = ['folder', 'File', 'A', 'SCF', 'B', 'Time', 'C']
# cols_siesta = ['folder', 'File', 'A', 'SCF', 'B', 'Time', 'C']
# labels_cp2k_1 = ['SIESTA NSlices 1', 'SIESTA-SMEAGOL V=0 NSlices 1', 'SIESTA-SMEAGOL V=1 NSlices 1']
# labels_siesta = ['SIESTA NSlices 5', 'SIESTA-SMEAGOL V=0 NSlices 5', 'SIESTA-SMEAGOL V=1 NSlices 5']
# cores = 1
# scale_cores_siesta = True
# scale_cores_cp2k = True

# Compare SIESTA ParallelOverKNum=1 and ParallelOverKNum=-1
# folder_cp2k = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/siesta-smeagol/capacitor/benchmarking/NEnergReal-40'
# folder_siesta = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/siesta-smeagol/capacitor/benchmarking/NEnergReal-40_ParallelOverKNum-2'
# cols_cp2k = ['folder', 'File', 'A', 'SCF', 'B', 'Time', 'C']
# cols_siesta = ['folder', 'File', 'A', 'SCF', 'B', 'Time', 'C']
# labels_cp2k_1 = ['SIESTA ParallelOverKNum=1', 'SIESTA-SMEAGOL V=0 ParallelOverKNum=1', 'SIESTA-SMEAGOL V=1 ParallelOverKNum=1']
# # labels_siesta = ['SIESTA ParallelOverKNum=-1', 'SIESTA-SMEAGOL V=0 ParallelOverKNum=-1', 'SIESTA-SMEAGOL V=1 ParallelOverKNum=-1']
# labels_siesta = ['SIESTA ParallelOverKNum=2', 'SIESTA-SMEAGOL V=0 ParallelOverKNum=2', 'SIESTA-SMEAGOL V=1 ParallelOverKNum=2']
# cores = 1
# scale_cores_siesta = True
# scale_cores_cp2k = True

# Au-BDT-Au SIESTA EM.OrderN F and SIESTA EM.OrderN T
# folder_cp2k = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/benchmarking/siesta/benchmarking/NEnergReal-80_OrderN-F'
# folder_siesta = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/benchmarking/siesta/benchmarking/NEnergReal-80_OrderN-T'
# sub_folders = ['atoms-116', 'atoms-246', 'atoms-428', 'atoms-948']
# cols_cp2k = ['folder', 'File', 'A', 'SCF', 'B', 'Time', 'C']
# cols_siesta = ['folder', 'File', 'A', 'SCF', 'B', 'Time', 'C']
# labels_cp2k_1 = ['SIESTA EM.OrderN F', 'SIESTA-SMEAGOL V=0 EM.OrderN F', 'SIESTA-SMEAGOL V=1 EM.OrderN F']
# labels_siesta = ['SIESTA EM.OrderN T', 'SIESTA-SMEAGOL V=0 EM.OrderN T', 'SIESTA-SMEAGOL V=1 EM.OrderN T']
# size_cp2k = [116, 246, 428, 948]
# size_siesta = [116, 246, 428, 948]
# cores = 1
# scale_cores_siesta = True
# scale_cores_cp2k = True
# plot_dft = False
# axis_limits = False
# ylim_1 = [0, 40]
# ylim_2 = [0, 135]
# ylim_3 = [0, 70]
# siesta_offset = 0
# num_kpoints = 3
# kpoints = num_kpoints * 0

#  Au-BDT-Au CP2K NEnergReal-80-memory_OrderF-F and NEnergReal-80-memory_OrderF-T
# folder_cp2k = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/benchmarking/cp2k/archer/benchmarking/NEnergReal-80-memory_OrderF-F'
# folder_siesta = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/benchmarking/cp2k/archer/benchmarking/NEnergReal-80-memory_OrderN-T'
# folder_cp2k = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/benchmarking/cp2k/archer/benchmarking/NEnergReal-80-memory_OrderN-F_tasks-per-node-64'
# folder_siesta = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/benchmarking/cp2k/archer/benchmarking/NEnergReal-80-memory_OrderN-T_tasks-per-node-64'
# folder_cp2k = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/benchmarking/cp2k/archer/benchmarking/NEnergReal-80-memory_OrderF-F'
# folder_siesta = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/benchmarking/cp2k/archer/benchmarking/NEnergReal-80-memory_OrderN-T-highmem'
# sub_folders = ['atoms-116', 'atoms-246', 'atoms-428', 'atoms-948']
# size_cp2k = [116, 246, 428, 948]
# size_siesta = [116, 246, 428, 948]
# cols_cp2k = ['folder', 'File', 'SCF', 'A', 'B', 'C', 'D', 'Time']
# cols_siesta = ['folder', 'File', 'SCF', 'A', 'B', 'C', 'D', 'Time']
# labels_cp2k_1 = ['CP2K', 'CP2K+SMEAGOL OrderN F', 'CP2K+SMEAGOL OrderN F']
# labels_siesta = ['CP2K', 'CP2K+SMEAGOL OrderN T', 'CP2K+SMEAGOL OrderN T']
# cores = 0
# scale_cores_siesta = False
# scale_cores_cp2k = False
# plot_dft = False
# axis_limits = False
# ylim_1 = [0, 40]
# ylim_2 = [0, 135]
# ylim_3 = [0, 70]
# siesta_offset = 0
# num_kpoints = 3
# kpoints = num_kpoints * 0

#  Au-BDT-Au CP2K NEnergReal-80-memory and SIESTA NEnergReal-80
# folder_cp2k = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/benchmarking/cp2k/archer/benchmarking/NEnergReal-80-memory_OrderF-F'
# folder_cp2k = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/benchmarking/cp2k/archer/benchmarking/NEnergReal-80-memory_OrderN-F_tasks-per-node-64'
# # folder_cp2k = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/benchmarking/cp2k/archer/benchmarking/NEnergReal-80-memory_OrderN-T'
# # folder_cp2k = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/benchmarking/cp2k/archer/benchmarking/NEnergReal-80-memory_OrderN-T_tasks-per-node-64'
# folder_siesta = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/benchmarking/siesta/benchmarking/scarf/NEnergReal-80_OrderN-F'
# # folder_siesta = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/benchmarking/siesta/benchmarking/scarf/NEnergReal-80_OrderN-F_bulk-EM.ParallelOverKNum'
# # folder_siesta = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/benchmarking/siesta/benchmarking/scarf/NEnergReal-80_OrderN-T_bulk-EM.ParallelOverKNum'
# sub_folders = ['atoms-116', 'atoms-246', 'atoms-428', 'atoms-948']
# cols_cp2k = ['folder', 'File', 'SCF', 'A', 'B', 'C', 'D', 'Time']
# cols_siesta = ['folder', 'File', 'A', 'SCF', 'B', 'Time', 'C']
# labels_cp2k_1 = ['CP2K', 'CP2K+SMEAGOL V=0 OrderN F', 'CP2K+SMEAGOL V=1 OrderN F']
# labels_siesta = ['SIESTA', 'SIESTA+SMEAGOL V=0 OrderN F', 'SIESTA+SMEAGOL V=1 OrderN F']
# # labels_cp2k_1 = ['CP2K', 'CP2K+SMEAGOL V=0 OrderN T', 'CP2K+SMEAGOL V=1 OrderN T']
# # labels_siesta = ['SIESTA', 'SIESTA+SMEAGOL V=0 OrderN T', 'SIESTA+SMEAGOL V=1 OrderN T']
# # size_cp2k = [116, 246, 428, 948]
# # size_siesta = [116, 246, 428, 948]
# size_cp2k = [1184, 2614, 4616, 10336]
# size_siesta = [1184, 2614, 4616, 10336]
# cores_list_cp2k = np.array([128, 512, 1024]) / 2
# cores_list_siesta = [40, 80, 120]
# cores_cp2k = 0
# cores_siesta = 1
# siesta_offset = -1
# num_kpoints = 3
# kpoints = num_kpoints * 0
# scale_cores_siesta = True
# scale_cores_cp2k = False
# plot_dft = False
# plot_fit = True
# axis_limits = False
# ylim_1 = [0, 40]
# ylim_2 = [0, 135]
# ylim_3 = [0, 70]

#  Au-BDT-Au CP2K omp-threads
# folders = ['/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/benchmarking/cp2k/archer/benchmarking/omp-threads/NEnergReal-128-memory_OrderF-F']
folders = ['/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/benchmarking/cp2k/archer/benchmarking/omp-threads/NEnergReal-80-memory_OrderF-F']
# folders = ['/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/benchmarking/cp2k/archer/benchmarking/omp-threads/NEnergReal-64-memory_OrderF-F-2']
sub_folders = ['/atoms-116/kpoints-1-1-20_hlb-auto', '/atoms-428/kpoints-1-1-20_hlb-auto', '/atoms-948/kpoints-1-1-20_hlb-auto']
files = ['scf-log_2.out', 'scf-log_3.out', 'scf-log_4.out']
plot_files = [True, True, True]
cols = [['folder', 'SCF', 'Time']]
labels_legend = [['CP2K', 'CP2K+SMEAGOL V=0 OrderN F', 'CP2K+SMEAGOL V=1 OrderN F']]
labels_axis_x = ['Number of electrons', 'Number of threads']
electrons = [1184, 4616, 10336]
cores = np.array([1, 2, 4, 8, 32])

bool_plot_atoms = True
plot_atoms = 0
bool_plot_cores = True
plot_cores = 0
bool_scale_time_by_cores = [False]

# cores_list_cp2k = np.array([1, 2, 4, 8, 32])
# cores_cp2k = 0
# num_kpoints = 0
# kpoints = num_kpoints * 0
# scale_cores_siesta = False
# scale_cores_cp2k = False
# plot_dft = False
# plot_fit = False
# axis_limits = False
# plot_siesta = False


data_file_1 = []
data_file_2 = []
data_file_3 = []
for i in range(len(folders)):
    data_file_1.append([pd.read_csv('{}/{}/{}'.format(folders[i], sub_folders[0], files[0]), names=cols[i], delim_whitespace=True).fillna(np.NaN),
                        pd.read_csv('{}/{}/{}'.format(folders[i], sub_folders[1], files[0]), names=cols[i], delim_whitespace=True).fillna(np.NaN),
                        pd.read_csv('{}/{}/{}'.format(folders[i], sub_folders[2], files[0]), names=cols[i], delim_whitespace=True).fillna(np.NaN)])
    data_file_2.append([pd.read_csv('{}/{}/{}'.format(folders[i], sub_folders[0], files[1]), names=cols[i], delim_whitespace=True).fillna(np.NaN),
                        pd.read_csv('{}/{}/{}'.format(folders[i], sub_folders[1], files[1]), names=cols[i], delim_whitespace=True).fillna(np.NaN),
                        pd.read_csv('{}/{}/{}'.format(folders[i], sub_folders[2], files[1]), names=cols[i], delim_whitespace=True).fillna(np.NaN)])
    data_file_3.append([pd.read_csv('{}/{}/{}'.format(folders[i], sub_folders[0], files[2]), names=cols[i], delim_whitespace=True).fillna(np.NaN),
                        pd.read_csv('{}/{}/{}'.format(folders[i], sub_folders[1], files[2]), names=cols[i], delim_whitespace=True).fillna(np.NaN),
                        pd.read_csv('{}/{}/{}'.format(folders[i], sub_folders[2], files[2]), names=cols[i], delim_whitespace=True).fillna(np.NaN)])

# Plot time vs atoms for kpoints 'kpoints'
# row_cp2k = cores_cp2k + kpoints
# if plot_siesta: row_siesta = cores_siesta + kpoints
# print('row_cp2k', row_cp2k)
# if plot_siesta: print('row_siesta', row_siesta)

# print(data_file_1[:][0]['Time'])
# print(data_file_2[0])
# print(data_file_3[0])

time_electrons_file_1 = []
time_electrons_file_2 = []
time_electrons_file_3 = []
for i in range(len(folders)):
    row = 0
    scale_factor = 1
    if bool_scale_time_by_cores: scale_factor = cores[plot_cores]

    time_electrons_file_1.append([data_file_1[i][0]['Time'][row] / data_file_1[i][0]['SCF'][row],
                                  data_file_1[i][1]['Time'][row] / data_file_1[i][1]['SCF'][row],
                                  data_file_1[i][2]['Time'][row] / data_file_1[i][2]['SCF'][row]])
    time_electrons_file_2.append([data_file_2[i][0]['Time'][row] / data_file_2[i][0]['SCF'][row],
                                  data_file_2[i][1]['Time'][row] / data_file_2[i][1]['SCF'][row],
                                  data_file_2[i][2]['Time'][row] / data_file_2[i][2]['SCF'][row]])
    time_electrons_file_3.append([data_file_3[i][0]['Time'][row] / data_file_3[i][0]['SCF'][row],
                                  data_file_3[i][1]['Time'][row] / data_file_3[i][1]['SCF'][row],
                                  data_file_3[i][2]['Time'][row] / data_file_3[i][2]['SCF'][row]])

time_cores_file_1 = []
time_cores_file_2 = []
time_cores_file_3 = []
for i in range(len(folders)):
    row = 0
    scale_factor = 1
    if bool_scale_time_by_cores: scale_factor = cores[plot_cores]

    time_cores_file_1.append([data_file_1[i][row]['Time'][0] / data_file_1[i][0]['SCF'][0],
                              data_file_1[i][row]['Time'][1] / data_file_1[i][1]['SCF'][1],
                              data_file_1[i][row]['Time'][2] / data_file_1[i][2]['SCF'][2],
                              data_file_1[i][row]['Time'][3] / data_file_1[i][2]['SCF'][3],
                              data_file_1[i][row]['Time'][4] / data_file_1[i][2]['SCF'][4]])
                            
    time_cores_file_2.append([data_file_2[i][row]['Time'][0] / data_file_2[i][0]['SCF'][0],
                              data_file_2[i][row]['Time'][1] / data_file_2[i][1]['SCF'][1],
                              data_file_2[i][row]['Time'][2] / data_file_2[i][2]['SCF'][2],
                              data_file_2[i][row]['Time'][3] / data_file_2[i][2]['SCF'][3],
                              data_file_2[i][row]['Time'][4] / data_file_2[i][2]['SCF'][4]])

    time_cores_file_3.append([data_file_3[i][row]['Time'][0] / data_file_3[i][0]['SCF'][0],
                              data_file_3[i][row]['Time'][1] / data_file_3[i][1]['SCF'][1],
                              data_file_3[i][row]['Time'][2] / data_file_3[i][2]['SCF'][2],
                              data_file_3[i][row]['Time'][3] / data_file_3[i][2]['SCF'][3],
                              data_file_3[i][row]['Time'][4] / data_file_3[i][2]['SCF'][4]])


# Fitting DFT
# if plot_fit: atoms_array = np.linspace(size_siesta[0], size_siesta[-1]*1, num=int(1e3))
# if plot_fit: cp2k_1_fit_2, _ = curve_fit(f_2, size_cp2k[:-1], cp2k_1[:-1] / cp2k_scale)
# if plot_fit and plot_siesta: siesta_1_fit_3, _ = curve_fit(f_3, size_siesta[:-1], siesta_1[:-1] / siesta_scale)

# Fitting SMEAGOL
# if plot_fit: atoms_array = np.linspace(size_siesta[0], size_siesta[-2]*1, num=int(1e3))
# if plot_fit: cp2k_2_fit_3, _ = curve_fit(f_3, size_cp2k[:-1], cp2k_2[:-1] / cp2k_scale)
# if plot_fit: cp2k_3_fit_3, _ = curve_fit(f_3, size_cp2k[:-1], cp2k_3[:-1] / cp2k_scale)

# Plot time taken vs electrons
print(time_electrons_file_3)
fig_plot_all_1, ax_plot_all_1 = plt.subplots()
for i in range(len(folders)):
    if plot_files[0]: ax_plot_all_1.plot(electrons, time_electrons_file_1[i], 'kx-', label=labels_legend[i][0])
    if plot_files[1]: ax_plot_all_1.plot(electrons, time_electrons_file_2[i], 'gx-', label=labels_legend[i][1])
    if plot_files[2]: ax_plot_all_1.plot(electrons, time_electrons_file_3[i], 'bx-', label=labels_legend[i][2])
# if plot_fit: ax_plot_all_1.plot(atoms_array, f_2(atoms_array, cp2k_1_fit_2[0], cp2k_1_fit_2[1]), '--', color='grey', alpha=0.5)
ax_plot_all_1.set_xlabel(labels_axis_x[0])
ax_plot_all_1.set_ylabel('Time taken per scf step / s')
ax_plot_all_1.get_xaxis().get_major_formatter().set_scientific(False)
ax_plot_all_1.legend(frameon=False)
fig_plot_all_1.tight_layout()
for i in range(len(folders)):
    fig_plot_all_1.savefig('{}/time_electrons.png'.format(folders[i]), dpi=param.save_dpi)

# Plot time taken vs MPI cores or OpenMP threads
print(time_cores_file_3)
fig_plot_all_2, ax_plot_all_2 = plt.subplots()
for i in range(len(folders)):
    if plot_files[0]: ax_plot_all_2.plot(cores, time_cores_file_1[i], 'kx-', label=labels_legend[i][0])
    if plot_files[1]: ax_plot_all_2.plot(cores, time_cores_file_2[i], 'gx-', label=labels_legend[i][1])
    if plot_files[2]: ax_plot_all_2.plot(cores, time_cores_file_3[i], 'bx-', label=labels_legend[i][2])
# if plot_fit: ax_plot_all_2.plot(atoms_array, f_2(atoms_array, cp2k_1_fit_2[0], cp2k_1_fit_2[1]), '--', color='grey', alpha=0.5)
ax_plot_all_2.set_xlabel(labels_axis_x[1])
ax_plot_all_2.set_ylabel('Time taken per scf step / s')
ax_plot_all_2.get_xaxis().get_major_formatter().set_scientific(False)
ax_plot_all_2.legend(frameon=False)
fig_plot_all_2.tight_layout()
for i in range(len(folders)):
    fig_plot_all_2.savefig('{}/time_cores.png'.format(folders[i]), dpi=param.save_dpi)
    
# Plot time vs kpoints for largest system size
# row_cp2k = cores_cp2k + kpoints
# if plot_siesta: row_siesta = cores_siesta + kpoints
# size_cp2k = ['1x1x1', '2x2x1', '4x4x1']
# size_siesta = ['1x1x1', '2x2x1', '4x4x1']
#
# cp2k_1 = np.array([cp2k_4_2['Time'][row_cp2k+num_kpoints*0]/cp2k_4_2['SCF'][row_cp2k+num_kpoints*0], cp2k_4_2['Time'][row_cp2k+num_kpoints*1]/cp2k_4_2['SCF'][row_cp2k+num_kpoints*1], cp2k_4_2['Time'][row_cp2k+num_kpoints*2]/cp2k_4_2['SCF'][row_cp2k+num_kpoints*2]])
# cp2k_2 = np.array([cp2k_4_3['Time'][row_cp2k+num_kpoints*0]/cp2k_4_3['SCF'][row_cp2k+num_kpoints*0], cp2k_4_3['Time'][row_cp2k+num_kpoints*1]/cp2k_4_3['SCF'][row_cp2k+num_kpoints*1], cp2k_4_3['Time'][row_cp2k+num_kpoints*2]/cp2k_4_3['SCF'][row_cp2k+num_kpoints*2]])
# cp2k_3 = np.array([cp2k_4_4['Time'][row_cp2k+num_kpoints*0]/cp2k_4_4['SCF'][row_cp2k+num_kpoints*0], cp2k_4_4['Time'][row_cp2k+num_kpoints*1]/cp2k_4_4['SCF'][row_cp2k+num_kpoints*1], cp2k_4_4['Time'][row_cp2k+num_kpoints*2]/cp2k_4_4['SCF'][row_cp2k+num_kpoints*2]])
#
# siesta_1 = np.array([siesta_4_2['Time'][row_siesta+num_kpoints*0]/(siesta_4_2['SCF'][row_siesta+num_kpoints*0]+siesta_offset),
#                      siesta_4_2['Time'][row_siesta+num_kpoints*1]/(siesta_4_2['SCF'][row_siesta+num_kpoints*1]+siesta_offset),
#                      siesta_4_2['Time'][row_siesta+num_kpoints*2]/(siesta_4_2['SCF'][row_siesta+num_kpoints*2]+siesta_offset)])
# if plot_siesta: siesta_2 = np.array([siesta_4_3['Time'][row_siesta+num_kpoints*0]/(siesta_4_3['SCF'][row_siesta+num_kpoints*0]+siesta_offset),
#                      siesta_4_3['Time'][row_siesta+num_kpoints*1]/(siesta_4_3['SCF'][row_siesta+num_kpoints*1]+siesta_offset),
#                      siesta_4_3['Time'][row_siesta+num_kpoints*2]/(siesta_4_3['SCF'][row_siesta+num_kpoints*2]+siesta_offset)])
# if plot_siesta: siesta_3 = np.array([siesta_4_4['Time'][row_siesta+num_kpoints*0]/(siesta_4_4['SCF'][row_siesta+num_kpoints*0]+siesta_offset),
#                      siesta_4_4['Time'][row_siesta+num_kpoints*1]/(siesta_4_4['SCF'][row_siesta+num_kpoints*1]+siesta_offset),
#                      siesta_4_4['Time'][row_siesta+num_kpoints*2]/(siesta_4_4['SCF'][row_siesta+num_kpoints*2]+siesta_offset)])
#
# fig_plot_dft_2, ax_plot_dft_2 = plt.subplots()
# ax_plot_dft_2.plot(size_cp2k, cp2k_1/ cp2k_scale, 'kx-', label=labels_cp2k_1[0])
# ax_plot_dft_2.plot(size_siesta, siesta_1/siesta_scale, 'kx--', label=labels_siesta[0])
# ax_plot_dft_2.set_xlabel('Kpoints')
# ax_plot_dft_2.set_ylabel('Time taken per scf step / s')
# ax_plot_dft_2.legend(frameon=False)
# fig_plot_dft_2.tight_layout()
# fig_plot_dft_2.savefig('{}/time_kpoints_dft.png'.format(folder_cp2k), dpi=param.save_dpi)
# fig_plot_dft_2.savefig('{}/time_kpoints_siesta.png'.format(folder_cp2k), dpi=param.save_dpi)
#
# fig_plot_all_2, ax_plot_all_2 = plt.subplots()
# if plot_dft: ax_plot_all_2.plot(size_cp2k, cp2k_1/ cp2k_scale, 'kx-', label=labels_cp2k_1[0])
# ax_plot_all_2.plot(size_cp2k, cp2k_2/ cp2k_scale, 'gx-', label=labels_cp2k_1[1])
# ax_plot_all_2.plot(size_cp2k, cp2k_3/ cp2k_scale, 'bx-', label=labels_cp2k_1[2])
# if plot_dft: ax_plot_all_2.plot(size_siesta, siesta_1/siesta_scale, 'kx--', label=labels_siesta[0])
# if plot_siesta: ax_plot_all_2.plot(size_siesta, siesta_2/siesta_scale, 'gx--', label=labels_siesta[1])
# if plot_siesta: ax_plot_all_2.plot(size_siesta, siesta_3/siesta_scale, 'bx--', label=labels_siesta[2])
# ax_plot_all_2.set_xlabel('Kpoints')
# ax_plot_all_2.set_ylabel('Time taken per scf step / s')
# ax_plot_all_2.legend(frameon=False)
# if axis_limits: ax_plot_all_2.set_ylim([ylim_2[0], ylim_2[1]])
# fig_plot_all_2.tight_layout()
# fig_plot_all_2.savefig('{}/time_kpoints_all.png'.format(folder_cp2k), dpi=param.save_dpi)
# if plot_siesta: fig_plot_all_2.savefig('{}/time_kpoints_all.png'.format(folder_siesta), dpi=param.save_dpi)
#
# # Plot time vs cores for largest system size and 2x2x1 kpoints
# row_cp2k = 0
# row_siesta = 0
#
# if scale_cores_cp2k:
#     cp2k_scale = cores_list_cp2k
# else:
#     # cp2k_scale = np.ones(cores_list_cp2k.shape[0])
#     cp2k_scale = np.ones(len(cores_list_cp2k))
# if scale_cores_siesta:
#     siesta_scale = cores_list_siesta
# else:
#     # siesta_scale = np.ones(cores_list_siesta.shape[0])
#     siesta_scale = np.ones(len(cores_list_siesta))
#
# cp2k_1 = np.array([cp2k_3_2['Time'][row_cp2k+0]/cp2k_3_2['SCF'][row_cp2k+0], cp2k_3_2['Time'][row_cp2k+1]/cp2k_3_2['SCF'][row_cp2k+1], cp2k_3_2['Time'][row_cp2k+2]/cp2k_3_2['SCF'][row_cp2k+2]])
# cp2k_2 = np.array([cp2k_3_3['Time'][row_cp2k+0]/cp2k_3_3['SCF'][row_cp2k+0], cp2k_3_3['Time'][row_cp2k+1]/cp2k_3_3['SCF'][row_cp2k+1], cp2k_3_3['Time'][row_cp2k+2]/cp2k_3_3['SCF'][row_cp2k+2]])
# cp2k_3 = np.array([cp2k_3_4['Time'][row_cp2k+0]/cp2k_3_4['SCF'][row_cp2k+0], cp2k_3_4['Time'][row_cp2k+1]/cp2k_3_4['SCF'][row_cp2k+1], cp2k_3_4['Time'][row_cp2k+2]/cp2k_3_4['SCF'][row_cp2k+2]])
#
# if plot_siesta: siesta_1 = np.array([siesta_3_2['Time'][row_siesta+0]/(siesta_3_2['SCF'][row_siesta+0]+ siesta_offset),
#                      siesta_3_2['Time'][row_siesta+1]/(siesta_3_2['SCF'][row_siesta+1]+ siesta_offset),
#                      siesta_3_2['Time'][row_siesta+2]/(siesta_3_2['SCF'][row_siesta+2]+ siesta_offset)])
# if plot_siesta: siesta_2 = np.array([siesta_3_3['Time'][row_siesta+0]/(siesta_3_3['SCF'][row_siesta+0]+ siesta_offset),
#                      siesta_3_3['Time'][row_siesta+1]/(siesta_3_3['SCF'][row_siesta+1]+ siesta_offset),
#                      siesta_3_3['Time'][row_siesta+2]/(siesta_3_3['SCF'][row_siesta+2]+ siesta_offset)])
# if plot_siesta: siesta_3 = np.array([siesta_3_4['Time'][row_siesta+0]/(siesta_3_4['SCF'][row_siesta+0]+ siesta_offset),
#                      siesta_3_4['Time'][row_siesta+1]/(siesta_3_4['SCF'][row_siesta+1]+ siesta_offset),
#                      siesta_3_4['Time'][row_siesta+2]/(siesta_3_4['SCF'][row_siesta+2]+ siesta_offset)])
#
# print('row_cp2k', row_cp2k)
# print('row_siesta', row_siesta)
#
# fig_plot_dft_3, ax_plot_dft_3 = plt.subplots()
# ax_plot_dft_3.plot(cores_list_cp2k, cp2k_1/ cp2k_scale, 'kx-', label=labels_cp2k_1[0])
# if plot_siesta: ax_plot_dft_3.plot(cores_list_siesta, siesta_1/siesta_scale, 'kx--', label=labels_siesta[0])
# ax_plot_dft_3.set_xlabel('Cores')
# ax_plot_dft_3.set_ylabel('Time taken per scf step / s')
# ax_plot_dft_3.legend(frameon=False)
# fig_plot_dft_3.tight_layout()
# fig_plot_dft_3.savefig('{}/time_cores_dft.png'.format(folder_cp2k), dpi=param.save_dpi)
# if plot_siesta: fig_plot_dft_3.savefig('{}/time_cores_dft.png'.format(folder_siesta), dpi=param.save_dpi)
#
# fig_plot_all_3, ax_plot_all_3 = plt.subplots()
# if plot_dft: ax_plot_all_3.plot(cores_list_cp2k, cp2k_1/ cp2k_scale, 'kx-', label=labels_cp2k_1[0])
# ax_plot_all_3.plot(cores_list_cp2k, cp2k_2/ cp2k_scale, 'gx-', label=labels_cp2k_1[1])
# ax_plot_all_3.plot(cores_list_cp2k, cp2k_3/ cp2k_scale, 'bx-', label=labels_cp2k_1[2])
# if plot_siesta and plot_dft: ax_plot_all_3.plot(cores_list_siesta, siesta_1/siesta_scale, 'kx--', label=labels_siesta[0])
# if plot_siesta: ax_plot_all_3.plot(cores_list_siesta, siesta_2/siesta_scale, 'gx--', label=labels_siesta[1])
# if plot_siesta: ax_plot_all_3.plot(cores_list_siesta, siesta_3/siesta_scale, 'bx--', label=labels_siesta[2])
# ax_plot_all_3.set_xlabel('Cores')
# ax_plot_all_3.set_ylabel('Time taken per scf step / s')
# ax_plot_all_3.legend(frameon=False)
# if axis_limits: ax_plot_all_3.set_ylim([ylim_3[0], ylim_3[1]])
# fig_plot_all_3.tight_layout()
# fig_plot_all_3.savefig('{}/time_cores_all.png'.format(folder_cp2k), dpi=param.save_dpi)
# if plot_siesta: fig_plot_all_3.savefig('{}/time_cores_all.png'.format(folder_siesta), dpi=param.save_dpi)

if __name__ == "__main__":
    print('Finished.')
    plt.show()
