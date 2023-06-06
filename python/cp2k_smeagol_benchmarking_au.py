import numpy as np
from matplotlib import pyplot as plt
from ase.io.cube import read_cube_data
from general import parameters as param
import pandas as pd

siesta_scale = 1
cp2k_scale = 1
scale_cores_siesta = False
scale_cores_cp2k = False
plot_dft = True
labels_cp2k_1 = ['CP2K', 'CP2K-SMEAGOL V=0', 'CP2K-SMEAGOL V=1']
labels_cp2k_2 = ['CP2K 2', 'CP2K-SMEAGOL 2 V=0', 'CP2K-SMEAGOL 2 V=1']
labels_siesta = ['SIESTA', 'SIESTA-SMEAGOL V=0', 'SIESTA-SMEAGOL V=1']

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

# Compare CP2K smeagol-interface-cp2k-20220907-memory OrderN F and SIESTA EM.OrderN F
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
folder_cp2k = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/cp2k-smeagol/sz/transmission/exp/capacitor/sergey-2d-equal/benchmarking/NEnergReal-80-memory'
folder_siesta = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/cp2k-smeagol/sz/transmission/exp/capacitor/sergey-2d-equal/benchmarking/NEnergReal-80-memory_ParallelOverKNum-auto'
cols_cp2k = ['folder', 'File', 'SCF', 'A', 'B', 'C', 'D', 'Time']
cols_siesta = ['folder', 'File', 'SCF', 'A', 'B', 'C', 'D', 'Time']
labels_cp2k_1 = ['CP2K ParallelOverKNum 1', 'CP2K-SMEAGOL V=0 ParallelOverKNum 1', 'CP2K-SMEAGOL V=1 ParallelOverKNum 1']
labels_siesta = ['CP2K ParallelOverKNum -1', 'CP2K-SMEAGOL V=0 ParallelOverKNum -1', 'CP2K-SMEAGOL V=1 ParallelOverKNum -1']
cores = 1
scale_cores_siesta = False
scale_cores_cp2k = False
plot_dft = False
axis_limits = True
ylim_1 = [0, 40]
ylim_2 = [0, 135]
ylim_3 = [0, 70]
siesta_offset = 0

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

if scale_cores_siesta: siesta_scale = 20 + cores * 20
if scale_cores_cp2k: cp2k_scale = 20 + cores * 20
replace_NaN = np.NaN
cp2k_1_2 = pd.read_csv('{}/layers-1/scf-log_2.out'.format(folder_cp2k), names=cols_cp2k, delim_whitespace=True).fillna(replace_NaN)
cp2k_1_3 = pd.read_csv('{}/layers-1/scf-log_3.out'.format(folder_cp2k), names=cols_cp2k, delim_whitespace=True).fillna(replace_NaN)
cp2k_1_4 = pd.read_csv('{}/layers-1/scf-log_4.out'.format(folder_cp2k), names=cols_cp2k, delim_whitespace=True).fillna(replace_NaN)
cp2k_2_2 = pd.read_csv('{}/layers-1-2/scf-log_2.out'.format(folder_cp2k), names=cols_cp2k, delim_whitespace=True).fillna(replace_NaN)
cp2k_2_3 = pd.read_csv('{}/layers-1-2/scf-log_3.out'.format(folder_cp2k), names=cols_cp2k, delim_whitespace=True).fillna(replace_NaN)
cp2k_2_4 = pd.read_csv('{}/layers-1-2/scf-log_4.out'.format(folder_cp2k), names=cols_cp2k, delim_whitespace=True).fillna(replace_NaN)
cp2k_3_2 = pd.read_csv('{}/layers-1-2-3/scf-log_2.out'.format(folder_cp2k), names=cols_cp2k, delim_whitespace=True).fillna(replace_NaN)
cp2k_3_3 = pd.read_csv('{}/layers-1-2-3/scf-log_3.out'.format(folder_cp2k), names=cols_cp2k, delim_whitespace=True).fillna(replace_NaN)
cp2k_3_4 = pd.read_csv('{}/layers-1-2-3/scf-log_4.out'.format(folder_cp2k), names=cols_cp2k, delim_whitespace=True).fillna(replace_NaN)
cp2k_4_2 = pd.read_csv('{}/layers-1-2-3-4/scf-log_2.out'.format(folder_cp2k), names=cols_cp2k, delim_whitespace=True).fillna(replace_NaN)
cp2k_4_3 = pd.read_csv('{}/layers-1-2-3-4/scf-log_3.out'.format(folder_cp2k), names=cols_cp2k, delim_whitespace=True).fillna(replace_NaN)
cp2k_4_4 = pd.read_csv('{}/layers-1-2-3-4/scf-log_4.out'.format(folder_cp2k), names=cols_cp2k, delim_whitespace=True).fillna(replace_NaN)
siesta_1_2 = pd.read_csv('{}/layers-1/scf-log_2.out'.format(folder_siesta), names=cols_siesta, delim_whitespace=True).fillna(replace_NaN)
siesta_1_3 = pd.read_csv('{}/layers-1/scf-log_3.out'.format(folder_siesta), names=cols_siesta, delim_whitespace=True).fillna(replace_NaN)
siesta_1_4 = pd.read_csv('{}/layers-1/scf-log_4.out'.format(folder_siesta), names=cols_siesta, delim_whitespace=True).fillna(replace_NaN)
siesta_2_2 = pd.read_csv('{}/layers-1-2/scf-log_2.out'.format(folder_siesta), names=cols_siesta, delim_whitespace=True).fillna(replace_NaN)
siesta_2_3 = pd.read_csv('{}/layers-1-2/scf-log_3.out'.format(folder_siesta), names=cols_siesta, delim_whitespace=True).fillna(replace_NaN)
siesta_2_4 = pd.read_csv('{}/layers-1-2/scf-log_4.out'.format(folder_siesta), names=cols_siesta, delim_whitespace=True).fillna(replace_NaN)
siesta_3_2 = pd.read_csv('{}/layers-1-2-3/scf-log_2.out'.format(folder_siesta), names=cols_siesta, delim_whitespace=True).fillna(replace_NaN)
siesta_3_3 = pd.read_csv('{}/layers-1-2-3/scf-log_3.out'.format(folder_siesta), names=cols_siesta, delim_whitespace=True).fillna(replace_NaN)
siesta_3_4 = pd.read_csv('{}/layers-1-2-3/scf-log_4.out'.format(folder_siesta), names=cols_siesta, delim_whitespace=True).fillna(replace_NaN)
siesta_4_2 = pd.read_csv('{}/layers-1-2-3-4/scf-log_2.out'.format(folder_siesta), names=cols_siesta, delim_whitespace=True).fillna(replace_NaN)
siesta_4_3 = pd.read_csv('{}/layers-1-2-3-4/scf-log_3.out'.format(folder_siesta), names=cols_siesta, delim_whitespace=True).fillna(replace_NaN)
siesta_4_4 = pd.read_csv('{}/layers-1-2-3-4/scf-log_4.out'.format(folder_siesta), names=cols_siesta, delim_whitespace=True).fillna(replace_NaN)

# Plot time vs atoms for kpoints 2x2
kpoints = 4 * 1
row = cores + kpoints
size_cp2k = [24, 24*2, 24*3, 24*4]
size_siesta = [24, 24*2, 24*3, 24*4]

cp2k_1 = np.array([cp2k_1_2['Time'][row]/cp2k_1_2['SCF'][row], cp2k_2_2['Time'][row]/cp2k_2_2['SCF'][row], cp2k_3_2['Time'][row]/cp2k_3_2['SCF'][row], cp2k_4_2['Time'][row]/cp2k_4_2['SCF'][row]])
cp2k_2 = np.array([cp2k_1_3['Time'][row]/cp2k_1_3['SCF'][row], cp2k_2_3['Time'][row]/cp2k_2_3['SCF'][row], cp2k_3_3['Time'][row]/cp2k_3_3['SCF'][row], cp2k_4_3['Time'][row]/cp2k_4_3['SCF'][row]])
cp2k_3 = np.array([cp2k_1_4['Time'][row]/cp2k_1_4['SCF'][row], cp2k_2_4['Time'][row]/cp2k_2_4['SCF'][row], cp2k_3_4['Time'][row]/cp2k_3_4['SCF'][row], cp2k_4_4['Time'][row]/cp2k_4_4['SCF'][row]])
siesta_1 = np.array([siesta_1_2['Time'][row]/(siesta_1_2['SCF'][row]+siesta_offset), siesta_2_2['Time'][row]/(siesta_2_2['SCF'][row]+siesta_offset),
                     siesta_3_2['Time'][row]/(siesta_3_2['SCF'][row]+siesta_offset), siesta_4_2['Time'][row]/(siesta_4_2['SCF'][row]+siesta_offset)])
siesta_2 = np.array([siesta_1_3['Time'][row]/(siesta_1_3['SCF'][row]+siesta_offset), siesta_2_3['Time'][row]/(siesta_2_3['SCF'][row]+siesta_offset),
                     siesta_3_3['Time'][row]/(siesta_3_3['SCF'][row]+siesta_offset), siesta_4_3['Time'][row]/(siesta_4_3['SCF'][row]+siesta_offset)])
siesta_3 = np.array([siesta_1_4['Time'][row]/(siesta_1_4['SCF'][row]+siesta_offset), siesta_2_4['Time'][row]/(siesta_2_4['SCF'][row]+siesta_offset), 
                     siesta_3_4['Time'][row]/(siesta_3_4['SCF'][row]+siesta_offset), siesta_4_4['Time'][row]/(siesta_4_4['SCF'][row]+siesta_offset)])

# Plot DFT only
fig_plot_dft_1, ax_plot_dft_1 = plt.subplots()
ax_plot_dft_1.plot(size_cp2k, cp2k_1 / cp2k_scale, 'kx-', label=labels_cp2k_1[0])
ax_plot_dft_1.plot(size_siesta, siesta_1 / siesta_scale, 'kx--', label=labels_siesta[0])
ax_plot_dft_1.set_xlabel('Number atoms')
ax_plot_dft_1.set_ylabel('Time taken per scf step / s')
ax_plot_dft_1.legend(frameon=False)
fig_plot_dft_1.tight_layout()
fig_plot_dft_1.savefig('{}/time_atoms_dft.png'.format(folder_cp2k), dpi=param.save_dpi)
fig_plot_dft_1.savefig('{}/time_atoms_dft.png'.format(folder_siesta), dpi=param.save_dpi)

# Plot all
fig_plot_all_1, ax_plot_all_1 = plt.subplots()
if plot_dft: ax_plot_all_1.plot(size_cp2k, cp2k_1/ cp2k_scale, 'kx-', label=labels_cp2k_1[0])
ax_plot_all_1.plot(size_cp2k, cp2k_2/ cp2k_scale, 'gx-', label=labels_cp2k_1[1])
ax_plot_all_1.plot(size_cp2k, cp2k_3/ cp2k_scale, 'bx-', label=labels_cp2k_1[2])
if plot_dft: ax_plot_all_1.plot(size_siesta, siesta_1/siesta_scale, 'kx--', label=labels_siesta[0])
ax_plot_all_1.plot(size_siesta, siesta_2/siesta_scale, 'gx--', label=labels_siesta[1])
ax_plot_all_1.plot(size_siesta, siesta_3/siesta_scale, 'bx--', label=labels_siesta[2])
ax_plot_all_1.set_xlabel('Number atoms')
ax_plot_all_1.set_ylabel('Time taken per scf step / s')
ax_plot_all_1.legend(frameon=False)
if axis_limits: ax_plot_all_1.set_ylim([ylim_1[0], ylim_1[1]])
fig_plot_all_1.tight_layout()
fig_plot_all_1.savefig('{}/time_atoms_all.png'.format(folder_cp2k), dpi=param.save_dpi)
fig_plot_all_1.savefig('{}/time_atoms_all.png'.format(folder_siesta), dpi=param.save_dpi)

# Plot time vs kpoints for largest system size
kpoints = 4 * 0
row = cores + kpoints
size_cp2k = ['1x1x1', '2x2x1', '4x4x1']
size_siesta = ['1x1x1', '2x2x1', '4x4x1']

cp2k_1 = np.array([cp2k_4_2['Time'][row+4*0]/cp2k_4_2['SCF'][row+4*0], cp2k_4_2['Time'][row+4*1]/cp2k_4_2['SCF'][row+4*1], cp2k_4_2['Time'][row+4*2]/cp2k_4_2['SCF'][row+4*2]])
cp2k_2 = np.array([cp2k_4_3['Time'][row+4*0]/cp2k_4_3['SCF'][row+4*0], cp2k_4_3['Time'][row+4*1]/cp2k_4_3['SCF'][row+4*1], cp2k_4_3['Time'][row+4*2]/cp2k_4_3['SCF'][row+4*2]])
cp2k_3 = np.array([cp2k_4_4['Time'][row+4*0]/cp2k_4_4['SCF'][row+4*0], cp2k_4_4['Time'][row+4*1]/cp2k_4_4['SCF'][row+4*1], cp2k_4_4['Time'][row+4*2]/cp2k_4_4['SCF'][row+4*2]])

siesta_1 = np.array([siesta_4_2['Time'][row+4*0]/(siesta_4_2['SCF'][row+4*0]+siesta_offset),
                     siesta_4_2['Time'][row+4*1]/(siesta_4_2['SCF'][row+4*1]+siesta_offset), 
                     siesta_4_2['Time'][row+4*2]/(siesta_4_2['SCF'][row+4*2]+siesta_offset)])
siesta_2 = np.array([siesta_4_3['Time'][row+4*0]/(siesta_4_3['SCF'][row+4*0]+siesta_offset), 
                     siesta_4_3['Time'][row+4*1]/(siesta_4_3['SCF'][row+4*1]+siesta_offset), 
                     siesta_4_3['Time'][row+4*2]/(siesta_4_3['SCF'][row+4*2]+siesta_offset)])
siesta_3 = np.array([siesta_4_4['Time'][row+4*0]/(siesta_4_4['SCF'][row+4*0]+siesta_offset), 
                     siesta_4_4['Time'][row+4*1]/(siesta_4_4['SCF'][row+4*1]+siesta_offset), 
                     siesta_4_4['Time'][row+4*2]/(siesta_4_4['SCF'][row+4*2]+siesta_offset)])

fig_plot_dft_2, ax_plot_dft_2 = plt.subplots()
ax_plot_dft_2.plot(size_cp2k, cp2k_1/ cp2k_scale, 'kx-', label=labels_cp2k_1[0])
ax_plot_dft_2.plot(size_siesta, siesta_1/siesta_scale, 'kx--', label=labels_siesta[0])
ax_plot_dft_2.set_xlabel('Kpoints')
ax_plot_dft_2.set_ylabel('Time taken per scf step / s')
ax_plot_dft_2.legend(frameon=False)
fig_plot_dft_2.tight_layout()
fig_plot_dft_2.savefig('{}/time_kpoints_dft.png'.format(folder_cp2k), dpi=param.save_dpi)
fig_plot_dft_2.savefig('{}/time_kpoints_siesta.png'.format(folder_cp2k), dpi=param.save_dpi)

fig_plot_all_2, ax_plot_all_2 = plt.subplots()
if plot_dft: ax_plot_all_2.plot(size_cp2k, cp2k_1/ cp2k_scale, 'kx-', label=labels_cp2k_1[0])
ax_plot_all_2.plot(size_cp2k, cp2k_2/ cp2k_scale, 'gx-', label=labels_cp2k_1[1])
ax_plot_all_2.plot(size_cp2k, cp2k_3/ cp2k_scale, 'bx-', label=labels_cp2k_1[2])
if plot_dft: ax_plot_all_2.plot(size_siesta, siesta_1/siesta_scale, 'kx--', label=labels_siesta[0])
ax_plot_all_2.plot(size_siesta, siesta_2/siesta_scale, 'gx--', label=labels_siesta[1])
ax_plot_all_2.plot(size_siesta, siesta_3/siesta_scale, 'bx--', label=labels_siesta[2])
ax_plot_all_2.set_xlabel('Kpoints')
ax_plot_all_2.set_ylabel('Time taken per scf step / s')
ax_plot_all_2.legend(frameon=False)
if axis_limits: ax_plot_all_2.set_ylim([ylim_2[0], ylim_2[1]])
fig_plot_all_2.tight_layout()
fig_plot_all_2.savefig('{}/time_kpoints_all.png'.format(folder_cp2k), dpi=param.save_dpi)
fig_plot_all_2.savefig('{}/time_kpoints_all.png'.format(folder_siesta), dpi=param.save_dpi)

# Plot time vs cores for largest system size and 2x2x1 kpoints
kpoints = 4 * 1
row = kpoints
size_cp2k = np.array([20, 40, 80, 120])
size_siesta = np.array([20, 40, 80, 120])

if scale_cores_cp2k: 
    cp2k_scale = size_cp2k
else:
    cp2k_scale = np.ones(size_cp2k.shape[0])
if scale_cores_siesta: 
    siesta_scale = size_siesta
else:
    siesta_scale = np.ones(size_siesta.shape[0])

cp2k_1 = np.array([cp2k_4_2['Time'][row+0]/cp2k_4_2['SCF'][row+0], cp2k_4_2['Time'][row+1]/cp2k_4_2['SCF'][row+1], cp2k_4_2['Time'][row+2]/cp2k_4_2['SCF'][row+2], cp2k_4_2['Time'][row+3]/cp2k_4_2['SCF'][row+3]])
cp2k_2 = np.array([cp2k_4_3['Time'][row+0]/cp2k_4_3['SCF'][row+0], cp2k_4_3['Time'][row+1]/cp2k_4_3['SCF'][row+1], cp2k_4_3['Time'][row+2]/cp2k_4_3['SCF'][row+2], cp2k_4_3['Time'][row+3]/cp2k_4_3['SCF'][row+3]])
cp2k_3 = np.array([cp2k_4_4['Time'][row+0]/cp2k_4_4['SCF'][row+0], cp2k_4_4['Time'][row+1]/cp2k_4_4['SCF'][row+1], cp2k_4_4['Time'][row+2]/cp2k_4_4['SCF'][row+2], cp2k_4_4['Time'][row+3]/cp2k_4_4['SCF'][row+3]])

siesta_1 = np.array([siesta_4_2['Time'][row+0]/(siesta_4_2['SCF'][row+0]+ siesta_offset), 
                     siesta_4_2['Time'][row+1]/(siesta_4_2['SCF'][row+1]+ siesta_offset), 
                     siesta_4_2['Time'][row+2]/(siesta_4_2['SCF'][row+2]+ siesta_offset),
                     siesta_4_2['Time'][row+3]/(siesta_4_2['SCF'][row+3]+ siesta_offset)])
siesta_2 = np.array([siesta_4_3['Time'][row+0]/(siesta_4_3['SCF'][row+0]+ siesta_offset), 
                     siesta_4_3['Time'][row+1]/(siesta_4_3['SCF'][row+1]+ siesta_offset), 
                     siesta_4_3['Time'][row+2]/(siesta_4_3['SCF'][row+2]+ siesta_offset),
                     siesta_4_3['Time'][row+3]/(siesta_4_3['SCF'][row+3]+ siesta_offset)])
siesta_3 = np.array([siesta_4_4['Time'][row+0]/(siesta_4_4['SCF'][row+0]+ siesta_offset), 
                     siesta_4_4['Time'][row+1]/(siesta_4_4['SCF'][row+1]+ siesta_offset), 
                     siesta_4_4['Time'][row+2]/(siesta_4_4['SCF'][row+2]+ siesta_offset),
                     siesta_4_4['Time'][row+3]/(siesta_4_4['SCF'][row+3]+ siesta_offset)])

fig_plot_dft_3, ax_plot_dft_3 = plt.subplots()
ax_plot_dft_3.plot(size_cp2k, cp2k_1/ cp2k_scale, 'kx-', label=labels_cp2k_1[0])
ax_plot_dft_3.plot(size_siesta, siesta_1/siesta_scale, 'kx--', label=labels_siesta[0])
ax_plot_dft_3.set_xlabel('Cores')
ax_plot_dft_3.set_ylabel('Time taken per scf step / s')
ax_plot_dft_3.legend(frameon=False)
fig_plot_dft_3.tight_layout()
fig_plot_dft_3.savefig('{}/time_cores_dft.png'.format(folder_cp2k), dpi=param.save_dpi)
fig_plot_dft_3.savefig('{}/time_cores_dft.png'.format(folder_siesta), dpi=param.save_dpi)

fig_plot_all_3, ax_plot_all_3 = plt.subplots()
if plot_dft: ax_plot_all_3.plot(size_cp2k, cp2k_1/ cp2k_scale, 'kx-', label=labels_cp2k_1[0])
ax_plot_all_3.plot(size_cp2k, cp2k_2/ cp2k_scale, 'gx-', label=labels_cp2k_1[1])
ax_plot_all_3.plot(size_cp2k, cp2k_3/ cp2k_scale, 'bx-', label=labels_cp2k_1[2])
if plot_dft: ax_plot_all_3.plot(size_siesta, siesta_1/siesta_scale, 'kx--', label=labels_siesta[0])
ax_plot_all_3.plot(size_siesta, siesta_2/siesta_scale, 'gx--', label=labels_siesta[1])
ax_plot_all_3.plot(size_siesta, siesta_3/siesta_scale, 'bx--', label=labels_siesta[2])
ax_plot_all_3.set_xlabel('Cores')
ax_plot_all_3.set_ylabel('Time taken per scf step / s')
ax_plot_all_3.legend(frameon=False)
if axis_limits: ax_plot_all_3.set_ylim([ylim_3[0], ylim_3[1]])
fig_plot_all_3.tight_layout()
fig_plot_all_3.savefig('{}/time_cores_all.png'.format(folder_cp2k), dpi=param.save_dpi)
fig_plot_all_3.savefig('{}/time_cores_all.png'.format(folder_siesta), dpi=param.save_dpi)

if __name__ == "__main__":
    print('Finished.')
    plt.show()
