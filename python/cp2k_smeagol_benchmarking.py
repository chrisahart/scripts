import numpy as np
from matplotlib import pyplot as plt
from ase.io.cube import read_cube_data
from general import parameters as param
import pandas as pd

folder_cp2k = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/cp2k-smeagol/sz/transmission/exp/capacitor/sergey-2d-equal/benchmarking'
cols = ['folder_cp2k', 'File', 'SCF', 'A', 'B', 'C', 'D', 'E']

file_1_2 = '{}/layers-1/scf-log_2.out'.format(folder_cp2k)
file_1_3 = '{}/layers-1/scf-log_3.out'.format(folder_cp2k)
file_1_4 = '{}/layers-1/scf-log_4.out'.format(folder_cp2k)
test_1_2 = pd.read_csv(file_1_2, names=cols, delim_whitespace=True)
test_1_3 = pd.read_csv(file_1_3, names=cols, delim_whitespace=True)
test_1_4 = pd.read_csv(file_1_4, names=cols, delim_whitespace=True)

file_2_2 = '{}/layers-1-2/scf-log_2.out'.format(folder_cp2k)
file_2_3 = '{}/layers-1-2/scf-log_3.out'.format(folder_cp2k)
file_2_4 = '{}/layers-1-2/scf-log_4.out'.format(folder_cp2k)
test_2_2 = pd.read_csv(file_2_2, names=cols, delim_whitespace=True)
test_2_3 = pd.read_csv(file_2_3, names=cols, delim_whitespace=True)
test_2_4 = pd.read_csv(file_2_4, names=cols, delim_whitespace=True)

file_3_2 = '{}/layers-1-2-3/scf-log_2.out'.format(folder_cp2k)
file_3_3 = '{}/layers-1-2-3/scf-log_3.out'.format(folder_cp2k)
file_3_4 = '{}/layers-1-2-3/scf-log_4.out'.format(folder_cp2k)
test_3_2 = pd.read_csv(file_3_2, names=cols, delim_whitespace=True)
test_3_3 = pd.read_csv(file_3_3, names=cols, delim_whitespace=True)
test_3_4 = pd.read_csv(file_3_4, names=cols, delim_whitespace=True)

file_4_2 = '{}/layers-1-2-3-4/scf-log_2.out'.format(folder_cp2k)
file_4_3 = '{}/layers-1-2-3-4/scf-log_3.out'.format(folder_cp2k)
file_4_4 = '{}/layers-1-2-3-4/scf-log_4.out'.format(folder_cp2k)
test_4_2 = pd.read_csv(file_4_2, names=cols, delim_whitespace=True)
test_4_3 = pd.read_csv(file_4_3, names=cols, delim_whitespace=True)
test_4_4 = pd.read_csv(file_4_4, names=cols, delim_whitespace=True)

# Plot vs row
# row = [1, 2, 4]
# fig_plot_1, ax_plot_1 = plt.subplots()
# ax_plot_1.plot(row, test_1_4['E']/test_1_4['SCF'], 'rx-', label='CP2K-SMEAGOL V=1 1')
# ax_plot_1.plot(row, test_1_3['E']/test_1_3['SCF'], 'rx-.', label='CP2K-SMEAGOL V=0 1')
# ax_plot_1.plot(row, test_1_2['E']/test_1_2['SCF'], 'rx--', label='CP2K 1')
# ax_plot_1.plot(row, test_2_3['E']/test_2_3['SCF'], 'gx-', label='CP2K-SMEAGOL V=0 2')
# ax_plot_1.plot(row, test_2_2['E']/test_2_2['SCF'], 'gx--', label='CP2K 2')
# ax_plot_1.plot(row, test_4_3['E']/test_4_3['SCF'], 'mx-', label='CP2K-SMEAGOL V=0 4')
# ax_plot_1.plot(row, test_4_2['E']/test_4_2['SCF'], 'mx--', label='CP2K 4')
# ax_plot_1.set_xlabel('Kpoint mesh')
# ax_plot_1.set_ylabel('Time taken per scf step / s')
# ax_plot_1.legend(frameon=False)
# fig_plot_1.tight_layout()
# fig_plot_1.savefig('{}/scf_step_time.png'.format(folder_cp2k), dpi=param.save_dpi)

# Plot vs atoms for row 1x1
# row = 3
# size = [24, 24*2, 24*3, 24*4]
# data_1 = np.array([test_1_2['E'][row]/test_1_2['SCF'][row], test_2_2['E'][row]/test_2_2['SCF'][row], test_3_2['E'][row]/test_3_2['SCF'][row], test_4_2['E'][row]/test_4_2['SCF'][row]])
# data_2 = np.array([test_1_3['E'][row]/test_1_3['SCF'][row], test_2_3['E'][row]/test_2_3['SCF'][row], test_3_3['E'][row]/test_3_3['SCF'][row], test_4_3['E'][row]/test_4_3['SCF'][row]])
# data_3 = np.array([test_1_4['E'][row]/test_1_4['SCF'][row], test_2_4['E'][row]/test_2_4['SCF'][row], test_3_4['E'][row]/test_3_4['SCF'][row], test_4_4['E'][row]/test_4_4['SCF'][row]])
# fig_plot_2, ax_plot_2 = plt.subplots()
# ax_plot_2.plot(size, data_1, 'rx-', label='CP2K')
# ax_plot_2.plot(size, data_2, 'gx-', label='CP2K-SMEAGOL V=0')
# ax_plot_2.plot(size, data_3, 'bx-', label='CP2K-SMEAGOL V=1')
# ax_plot_2.set_xlabel('Number atoms')
# ax_plot_2.set_ylabel('Time taken per scf step / s')
# ax_plot_2.legend(frameon=False)
# fig_plot_2.tight_layout()
# fig_plot_2.savefig('{}/scf_step_time_row_1.png'.format(folder_cp2k), dpi=param.save_dpi)

# Plot vs atoms for row 2x2
row = 3 + 4
size = [24, 24*2, 24*3, 24*4]
data_1 = np.array([test_1_2['E'][row]/test_1_2['SCF'][row], test_2_2['E'][row]/test_2_2['SCF'][row], test_3_2['E'][row]/test_3_2['SCF'][row], test_4_2['E'][row]/test_4_2['SCF'][row]])
data_2 = np.array([test_1_3['E'][row]/test_1_3['SCF'][row], test_2_3['E'][row]/test_2_3['SCF'][row], test_3_3['E'][row]/test_3_3['SCF'][row], test_4_3['E'][row]/test_4_3['SCF'][row]])
data_3 = np.array([test_1_4['E'][row]/test_1_4['SCF'][row], test_2_4['E'][row]/test_2_4['SCF'][row], test_3_4['E'][row]/test_3_4['SCF'][row], test_4_4['E'][row]/test_4_4['SCF'][row]])
fig_plot_3, ax_plot_3 = plt.subplots()
ax_plot_3.plot(size, data_1, 'rx-', label='CP2K')
ax_plot_3.plot(size, data_2, 'gx-', label='CP2K-SMEAGOL V=0')
ax_plot_3.plot(size, data_3, 'bx-', label='CP2K-SMEAGOL V=1')
ax_plot_3.set_xlabel('Number atoms')
ax_plot_3.set_ylabel('Time taken per scf step / s')
ax_plot_3.legend(frameon=False)
fig_plot_3.tight_layout()
fig_plot_3.savefig('{}/scf_step_time_row_2.png'.format(folder_cp2k), dpi=param.save_dpi)

# Plot vs atoms for row 4x4
row = 3 + 4*2
size = [24, 24*2, 24*3, 24*4]
data_1 = np.array([test_1_2['E'][row]/test_1_2['SCF'][row], test_2_2['E'][row]/test_2_2['SCF'][row], test_3_2['E'][row]/test_3_2['SCF'][row], test_4_2['E'][row]/test_4_2['SCF'][row]])
data_2 = np.array([test_1_3['E'][row]/test_1_3['SCF'][row], test_2_3['E'][row]/test_2_3['SCF'][row], test_3_3['E'][row]/test_3_3['SCF'][row], test_4_3['E'][row]/test_4_3['SCF'][row]])
data_3 = np.array([test_1_4['E'][row]/test_1_4['SCF'][row], test_2_4['E'][row]/test_2_4['SCF'][row], test_3_4['E'][row]/test_3_4['SCF'][row], test_4_4['E'][row]/test_4_4['SCF'][row]])
fig_plot_4, ax_plot_4 = plt.subplots()
ax_plot_4.plot(size, data_1, 'rx-', label='CP2K')
ax_plot_4.plot(size, data_2, 'gx-', label='CP2K-SMEAGOL V=0')
ax_plot_4.plot(size, data_3, 'bx-', label='CP2K-SMEAGOL V=1')
ax_plot_4.set_xlabel('Number atoms')
ax_plot_4.set_ylabel('Time taken per scf step / s')
ax_plot_4.legend(frameon=False)
fig_plot_4.tight_layout()
fig_plot_4.savefig('{}/scf_step_time_row_4.png'.format(folder_cp2k), dpi=param.save_dpi)

if __name__ == "__main__":
    print('Finished.')
    plt.show()
