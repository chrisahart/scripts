import numpy as np
from matplotlib import pyplot as plt
from ase.io.cube import read_cube_data
from general import parameters as param
import pandas as pd

folder_cp2k = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/cp2k-smeagol/sz/transmission/exp/capacitor/sergey-2d-equal/benchmarking'
# folder_cp2k = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/siesta-smeagol/capacitor/benchmarking/NEnergReal-640'
cols_cp2k = ['folder', 'File', 'SCF', 'A', 'B', 'C', 'D', 'Time']
# cols_cp2k = ['folder', 'File', 'A', 'SCF', 'B', 'Time', 'C']
folder_siesta = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/siesta-smeagol/capacitor/benchmarking/NEnergReal-80'
cols_siesta = ['folder', 'File', 'A', 'SCF', 'B', 'Time', 'C']
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
cores = 1
cp2k_scale = 40
siesta_scale = 40
row = cores + kpoints
size_cp2k = [24, 24*2, 24*3, 24*4]
size_siesta = [24, 24*2, 24*3, 24*4]

cp2k_1 = np.array([cp2k_1_2['Time'][row]/cp2k_1_2['SCF'][row], cp2k_2_2['Time'][row]/cp2k_2_2['SCF'][row], cp2k_3_2['Time'][row]/cp2k_3_2['SCF'][row], cp2k_4_2['Time'][row]/cp2k_4_2['SCF'][row]])
cp2k_2 = np.array([cp2k_1_3['Time'][row]/cp2k_1_3['SCF'][row], cp2k_2_3['Time'][row]/cp2k_2_3['SCF'][row], cp2k_3_3['Time'][row]/cp2k_3_3['SCF'][row], cp2k_4_3['Time'][row]/cp2k_4_3['SCF'][row]])
cp2k_3 = np.array([cp2k_1_4['Time'][row]/cp2k_1_4['SCF'][row], cp2k_2_4['Time'][row]/cp2k_2_4['SCF'][row], cp2k_3_4['Time'][row]/cp2k_3_4['SCF'][row], cp2k_4_4['Time'][row]/cp2k_4_4['SCF'][row]])
siesta_1 = np.array([siesta_1_2['Time'][row]/siesta_1_2['SCF'][row], siesta_2_2['Time'][row]/siesta_2_2['SCF'][row], siesta_3_2['Time'][row]/siesta_3_2['SCF'][row], siesta_4_2['Time'][row]/siesta_4_2['SCF'][row]])
siesta_2 = np.array([siesta_1_3['Time'][row]/siesta_1_3['SCF'][row], siesta_2_3['Time'][row]/siesta_2_3['SCF'][row], siesta_3_3['Time'][row]/siesta_3_3['SCF'][row], siesta_4_3['Time'][row]/siesta_4_3['SCF'][row]])
siesta_3 = np.array([siesta_1_4['Time'][row]/siesta_1_4['SCF'][row], siesta_2_4['Time'][row]/siesta_2_4['SCF'][row], siesta_3_4['Time'][row]/siesta_3_4['SCF'][row], siesta_4_4['Time'][row]/siesta_4_4['SCF'][row]])

# Plot DFT only
fig_plot_dft_1, ax_plot_dft_1 = plt.subplots()
ax_plot_dft_1.plot(size_cp2k, cp2k_1 / cp2k_scale, 'kx-', label='CP2K')
ax_plot_dft_1.plot(size_siesta, siesta_1 / siesta_scale, 'kx--', label='SIESTA')
ax_plot_dft_1.set_xlabel('Number atoms')
ax_plot_dft_1.set_ylabel('Time taken per scf step / s')
ax_plot_dft_1.legend(frameon=False)
fig_plot_dft_1.tight_layout()
fig_plot_dft_1.savefig('{}/time_atoms_dft.png'.format(folder_cp2k), dpi=param.save_dpi)
fig_plot_dft_1.savefig('{}/time_atoms_dft.png'.format(folder_siesta), dpi=param.save_dpi)

# Plot all
fig_plot_all_1, ax_plot_all_1 = plt.subplots()
ax_plot_all_1.plot(size_cp2k, cp2k_1/ cp2k_scale, 'kx-', label='CP2K')
ax_plot_all_1.plot(size_cp2k, cp2k_2/ cp2k_scale, 'gx-', label='CP2K-SMEAGOL V=0')
ax_plot_all_1.plot(size_cp2k, cp2k_3/ cp2k_scale, 'bx-', label='CP2K-SMEAGOL V=1')
ax_plot_all_1.plot(size_siesta, siesta_1/siesta_scale, 'kx--', label='SIESTA')
ax_plot_all_1.plot(size_siesta, siesta_2/siesta_scale, 'gx--', label='SIESTA-SMEAGOL V=0')
ax_plot_all_1.plot(size_siesta, siesta_3/siesta_scale, 'bx--', label='SIESTA-SMEAGOL V=1')
ax_plot_all_1.set_xlabel('Number atoms')
ax_plot_all_1.set_ylabel('Time taken per scf step / s')
ax_plot_all_1.legend(frameon=False)
fig_plot_all_1.tight_layout()
fig_plot_all_1.savefig('{}/time_atoms_all.png'.format(folder_cp2k), dpi=param.save_dpi)
fig_plot_all_1.savefig('{}/time_atoms_all.png'.format(folder_siesta), dpi=param.save_dpi)

# Plot time vs kpoints for largest system size
kpoints = 0 * 1
cores = 1
cp2k_scale = 40
siesta_scale = 40
row = cores + kpoints
size_cp2k = ['1x1x1', '2x2x1', '4x4x1']
size_siesta = ['1x1x1', '2x2x1', '4x4x1']

cp2k_1 = np.array([cp2k_4_2['Time'][row+4*0]/cp2k_4_2['SCF'][row+4*0], cp2k_4_2['Time'][row+4*1]/cp2k_4_2['SCF'][row+4*1], cp2k_4_2['Time'][row+4*2]/cp2k_4_2['SCF'][row+4*2]])
cp2k_2 = np.array([cp2k_4_3['Time'][row+4*0]/cp2k_4_3['SCF'][row+4*0], cp2k_4_3['Time'][row+4*1]/cp2k_4_3['SCF'][row+4*1], cp2k_4_3['Time'][row+4*2]/cp2k_4_3['SCF'][row+4*2]])
cp2k_3 = np.array([cp2k_4_4['Time'][row+4*0]/cp2k_4_4['SCF'][row+4*0], cp2k_4_4['Time'][row+4*1]/cp2k_4_4['SCF'][row+4*1], cp2k_4_4['Time'][row+4*2]/cp2k_4_4['SCF'][row+4*2]])
siesta_1 = np.array([siesta_4_2['Time'][row+4*0]/siesta_4_2['SCF'][row+4*0], siesta_4_2['Time'][row+4*1]/siesta_4_2['SCF'][row+4*1], siesta_4_2['Time'][row+4*2]/siesta_4_2['SCF'][row+4*2]])
siesta_2 = np.array([siesta_4_3['Time'][row+4*0]/siesta_4_3['SCF'][row+4*0], siesta_4_3['Time'][row+4*1]/siesta_4_3['SCF'][row+4*1], siesta_4_3['Time'][row+4*2]/siesta_4_3['SCF'][row+4*2]])
siesta_3 = np.array([siesta_4_4['Time'][row+4*0]/siesta_4_4['SCF'][row+4*0], siesta_4_4['Time'][row+4*1]/siesta_4_4['SCF'][row+4*1], siesta_4_4['Time'][row+4*2]/siesta_4_4['SCF'][row+4*2]])

fig_plot_dft_2, ax_plot_dft_2 = plt.subplots()
ax_plot_dft_2.plot(size_cp2k, cp2k_1/ cp2k_scale, 'kx-', label='CP2K')
ax_plot_dft_2.plot(size_siesta, siesta_1/siesta_scale, 'kx--', label='SIESTA')
ax_plot_dft_2.set_xlabel('Kpoints')
ax_plot_dft_2.set_ylabel('Time taken per scf step / s')
ax_plot_dft_2.legend(frameon=False)
fig_plot_dft_2.tight_layout()
fig_plot_dft_2.savefig('{}/time_kpoints_dft.png'.format(folder_cp2k), dpi=param.save_dpi)
fig_plot_dft_2.savefig('{}/time_kpoints_siesta.png'.format(folder_cp2k), dpi=param.save_dpi)

fig_plot_all_2, ax_plot_all_2 = plt.subplots()
ax_plot_all_2.plot(size_cp2k, cp2k_1/ cp2k_scale, 'kx-', label='CP2K')
ax_plot_all_2.plot(size_cp2k, cp2k_2/ cp2k_scale, 'gx-', label='CP2K-SMEAGOL V=0')
ax_plot_all_2.plot(size_cp2k, cp2k_3/ cp2k_scale, 'bx-', label='CP2K-SMEAGOL V=1')
ax_plot_all_2.plot(size_siesta, siesta_1/siesta_scale, 'kx--', label='SIESTA')
ax_plot_all_2.plot(size_siesta, siesta_2/siesta_scale, 'gx--', label='SIESTA-SMEAGOL V=0')
ax_plot_all_2.plot(size_siesta, siesta_3/siesta_scale, 'bx--', label='SIESTA-SMEAGOL V=1')
ax_plot_all_2.set_xlabel('Kpoints')
ax_plot_all_2.set_ylabel('Time taken per scf step / s')
ax_plot_all_2.legend(frameon=False)
fig_plot_all_2.tight_layout()
fig_plot_all_2.savefig('{}/time_kpoints_all.png'.format(folder_cp2k), dpi=param.save_dpi)
fig_plot_all_2.savefig('{}/time_kpoints_all.png'.format(folder_siesta), dpi=param.save_dpi)

# Plot time vs cores for largest system size and 2x2x1 kpoints
kpoints = 4 * 1
row = kpoints
size_cp2k = [20, 40, 80, 120]
size_siesta = [20, 40, 80, 120]

cp2k_1 = np.array([cp2k_4_2['Time'][row+0]/cp2k_4_2['SCF'][row+0], cp2k_4_2['Time'][row+1]/cp2k_4_2['SCF'][row+1], cp2k_4_2['Time'][row+2]/cp2k_4_2['SCF'][row+2], cp2k_4_2['Time'][row+3]/cp2k_4_2['SCF'][row+3]])
cp2k_2 = np.array([cp2k_4_3['Time'][row+0]/cp2k_4_3['SCF'][row+0], cp2k_4_3['Time'][row+1]/cp2k_4_3['SCF'][row+1], cp2k_4_3['Time'][row+2]/cp2k_4_3['SCF'][row+2], cp2k_4_3['Time'][row+3]/cp2k_4_3['SCF'][row+3]])
cp2k_3 = np.array([cp2k_4_4['Time'][row+0]/cp2k_4_4['SCF'][row+0], cp2k_4_4['Time'][row+1]/cp2k_4_4['SCF'][row+1], cp2k_4_4['Time'][row+2]/cp2k_4_4['SCF'][row+2], cp2k_4_4['Time'][row+3]/cp2k_4_4['SCF'][row+3]])
siesta_1 = np.array([siesta_4_2['Time'][row+0]/siesta_4_2['SCF'][row+0], siesta_4_2['Time'][row+1]/siesta_4_2['SCF'][row+1], siesta_4_2['Time'][row+2]/siesta_4_2['SCF'][row+2], siesta_4_2['Time'][row+3]/siesta_4_2['SCF'][row+3]])
siesta_2 = np.array([siesta_4_3['Time'][row+0]/siesta_4_3['SCF'][row+0], siesta_4_3['Time'][row+1]/siesta_4_3['SCF'][row+1], siesta_4_3['Time'][row+2]/siesta_4_3['SCF'][row+2], siesta_4_3['Time'][row+3]/siesta_4_3['SCF'][row+3]])
siesta_3 = np.array([siesta_4_4['Time'][row+0]/siesta_4_4['SCF'][row+0], siesta_4_4['Time'][row+1]/siesta_4_4['SCF'][row+1], siesta_4_4['Time'][row+2]/siesta_4_4['SCF'][row+2], siesta_4_4['Time'][row+3]/siesta_4_4['SCF'][row+3]])

fig_plot_dft_3, ax_plot_dft_3 = plt.subplots()
ax_plot_dft_3.plot(size_cp2k, cp2k_1/ cp2k_scale, 'kx-', label='CP2K')
ax_plot_dft_3.plot(size_siesta, siesta_1/size_siesta, 'kx--', label='SIESTA')
ax_plot_dft_3.set_xlabel('Cores')
ax_plot_dft_3.set_ylabel('Time taken per scf step / s')
ax_plot_dft_3.legend(frameon=False)
fig_plot_dft_3.tight_layout()
fig_plot_dft_3.savefig('{}/time_cores_dft.png'.format(folder_cp2k), dpi=param.save_dpi)
fig_plot_dft_3.savefig('{}/time_cores_dft.png'.format(folder_siesta), dpi=param.save_dpi)

fig_plot_all_3, ax_plot_all_3 = plt.subplots()
ax_plot_all_3.plot(size_cp2k, cp2k_1/ cp2k_scale, 'kx-', label='CP2K')
ax_plot_all_3.plot(size_cp2k, cp2k_2/ cp2k_scale, 'gx-', label='CP2K-SMEAGOL V=0')
ax_plot_all_3.plot(size_cp2k, cp2k_3/ cp2k_scale, 'bx-', label='CP2K-SMEAGOL V=1')
ax_plot_all_3.plot(size_siesta, siesta_1/size_siesta, 'kx--', label='SIESTA')
ax_plot_all_3.plot(size_siesta, siesta_2/size_siesta, 'gx--', label='SIESTA-SMEAGOL V=0')
ax_plot_all_3.plot(size_siesta, siesta_3/size_siesta, 'bx--', label='SIESTA-SMEAGOL V=1')
ax_plot_all_3.set_xlabel('Cores')
ax_plot_all_3.set_ylabel('Time taken per scf step / s')
ax_plot_all_3.legend(frameon=False)
fig_plot_all_3.tight_layout()
fig_plot_all_3.savefig('{}/time_cores_all.png'.format(folder_cp2k), dpi=param.save_dpi)
fig_plot_all_3.savefig('{}/time_cores_all.png'.format(folder_siesta), dpi=param.save_dpi)

if __name__ == "__main__":
    print('Finished.')
    plt.show()
