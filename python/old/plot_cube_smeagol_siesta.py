import numpy as np
from matplotlib import pyplot as plt
from ase.io.cube import read_cube_data
from general import parameters as param

""" Plotting of CP2K .cube files for Hartree potential and charge density """

params = {'axes.formatter.limits': [-4, 4],
          'axes.labelsize': 'large',
          'axes.titlesize': 'large',
          'legend.fontsize': 'medium',
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
print_label = 'V0_DFT'

folder_cp2k = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/supercell-1-1-5-bulk-6-cu-1.86/capacitor/bias/kpoints-2-2-symmetric-screen']
folder_cp2k = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/supercell-1-1-5-bulk-6-cu-1.86/capacitor/bias/kpoints-2-2']
folder_cp2k = ["/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/supercell-1-1-5-bulk-6-cu-1.86/junction/bias/energy/kpoints-2-2-V-0"]
# electrons = np.array([24, 0, 420, 420])
# electrons = np.array([24, 0, 64, 64])
# /np.sum(average_charge_bulk[plot_folder])*electrons[0]

file_charge_em2 = ['1V-ELECTRON_DENSITY-1_0.cube']
file_hartree_em2 = ['1V-v_hartree-1_0.cube']
diff_label = ['CP2K+SMEAGOL 1 V - 0 V']
file_charge_em1 = ['0V-ELECTRON_DENSITY-1_0.cube'] * 4
file_hartree_em1 = ['0V-v_hartree-1_0.cube'] * 4
file_charge_dft = ['dft_wfn-ELECTRON_DENSITY-1_0.cube'] * 4
file_hartree_dft = ['dft_wfn-v_hartree-1_0.cube'] * 4
file_charge_bulk = ['bulk-ELECTRON_DENSITY-1_0.cube'] * 4
file_hartree_bulk = ['bulk-v_hartree-1_0.cube'] * 4
labels = ['CP2K+SMEAGOL bulk', 'CP2K', 'CP2K+SMEAGOL V=0', 'CP2K+SMEAGOL V=1']
labels = ['Bulk', 'CP2K', 'V=0', 'V=1']
plot_folder = 0

plot_color = ['k', 'b', 'r', 'g']
read_labels = [True, False, True, True]

diff_color = ['g', 'b']
diff_color_siesta = ['m']
plot_labels = read_labels
print_label = 'average_z_all'
axis = 2
plot_diff_legend = [True, False]
draw_mirror = False
draw_mirror_diff = False
mirror_scale = -1
use_xlim = True
use_xlim = False
xlim_values = np.array([14.43900, 57.81600])

# Read .cube using ASE
data_charge_em2 = []
data_hartree_em2 = []
data_charge_em1 = []
data_hartree_em1 = []
data_charge_bulk = []
data_charge_dft = []
data_hartree_dft = []
data_hartree_bulk = []
for i in range(len(folder_cp2k)):
    if read_labels[3]: data_charge_em2.append(read_cube_data('{}/{}'.format(folder_cp2k[i], file_charge_em2[i])))
    if read_labels[3]: data_hartree_em2.append(read_cube_data('{}/{}'.format(folder_cp2k[i], file_hartree_em2[i])))
    if read_labels[2]: data_charge_em1.append(read_cube_data('{}/{}'.format(folder_cp2k[i], file_charge_em1[i])))
    if read_labels[2]: data_hartree_em1.append(read_cube_data('{}/{}'.format(folder_cp2k[i], file_hartree_em1[i])))
    if read_labels[1]: data_charge_dft.append(read_cube_data('{}/{}'.format(folder_cp2k[i], file_charge_dft[i])))
    if read_labels[1]: data_hartree_dft.append(read_cube_data('{}/{}'.format(folder_cp2k[i], file_hartree_dft[i])))
    if read_labels[0]: data_charge_bulk.append(read_cube_data('{}/{}'.format(folder_cp2k[i], file_charge_bulk[i])))
    if read_labels[0]: data_hartree_bulk.append(read_cube_data('{}/{}'.format(folder_cp2k[i], file_hartree_bulk[i])))
print('Finished reading .cube files')
# atoms, data = data_charge_em2[0]
# print(data_charge_em2[0][0])
# print(np.sum(data_charge_em2[0][0]))

# Calculate average along axis
average_charge_em2 = []
average_hartree_em2 = []
average_charge_em1 = []
average_hartree_em1 = []
average_charge_bulk = []
average_charge_dft = []
average_hartree_dft = []
average_hartree_bulk = []
for i in range(len(folder_cp2k)):
    if read_labels[3]: average_charge_em2.append(np.zeros(data_charge_em2[i][0].shape[axis]))
    if read_labels[3]: average_hartree_em2.append(np.zeros(data_hartree_em2[i][0].shape[axis]))
    if read_labels[2]: average_charge_em1.append(np.zeros(data_charge_em1[i][0].shape[axis]))
    if read_labels[2]: average_hartree_em1.append(np.zeros(data_hartree_em1[i][0].shape[axis]))
    if read_labels[1]: average_charge_dft.append(np.zeros(data_charge_dft[i][0].shape[axis]))
    if read_labels[1]: average_hartree_dft.append(np.zeros(data_hartree_dft[i][0].shape[axis]))
    if read_labels[0]: average_charge_bulk.append(np.zeros(data_charge_bulk[i][0].shape[axis]))
    if read_labels[0]: average_hartree_bulk.append(np.zeros(data_hartree_bulk[i][0].shape[axis]))
    
    for j in range(data_charge_em1[i][0].shape[axis]):
        if read_labels[3]: average_charge_em2[i][j] = np.mean(data_charge_em2[i][0][:, :, j])
        if read_labels[3]: average_hartree_em2[i][j] = np.mean(data_hartree_em2[i][0][:, :, j] * param.hartree_to_ev)
        if read_labels[2]: average_charge_em1[i][j] = np.mean(data_charge_em1[i][0][:, :, j])
        if read_labels[2]: average_hartree_em1[i][j] = np.mean(data_hartree_em1[i][0][:, :, j] * param.hartree_to_ev)
        if read_labels[1]: average_charge_dft[i][j] = np.mean(data_charge_dft[i][0][:, :, j])
        if read_labels[1]: average_hartree_dft[i][j] = np.mean(data_hartree_dft[i][0][:, :, j] * param.hartree_to_ev)
    for j in range(data_charge_bulk[i][0].shape[axis]):
        if read_labels[0]: average_charge_bulk[i][j] = np.mean(data_charge_bulk[i][0][:, :, j])
        if read_labels[0]: average_hartree_bulk[i][j] = np.mean(data_hartree_bulk[i][0][:, :, j] * param.hartree_to_ev)

# Setup axis
energy_grid_hartree_em1 = np.linspace(start=0, stop=data_charge_em1[0][1].get_cell()[axis][axis], num=average_charge_em1[0].shape[0])
energy_grid_hartree_bulk = np.linspace(start=0, stop=data_charge_bulk[0][1].get_cell()[axis][axis], num=average_charge_bulk[0].shape[0])

# Plot all
rows, cols = 2, 1
xlim = [0-1, data_charge_em1[0][1].get_cell()[axis][axis]+1]
if use_xlim:
    xlim = xlim_values
fig_cube_z, ax_cube_z = plt.subplots(rows, cols, sharex='col', sharey='row', figsize=(6, 8))
if plot_fermi: ax_cube_z[0].axhline(y=fermi_dft, color='grey', linestyle='--', label='DFT Fermi energy', alpha=0.5)
if plot_labels[1]: ax_cube_z[0].plot(energy_grid_hartree_em1, average_hartree_dft[plot_folder], '-', color=plot_color[1], label=labels[1])
if plot_labels[3]: ax_cube_z[0].plot(energy_grid_hartree_em1, average_hartree_em2[plot_folder], '-', color=plot_color[3], label=labels[3])
if plot_labels[2]: ax_cube_z[0].plot(energy_grid_hartree_em1, average_hartree_em1[plot_folder], '-', color=plot_color[2], label=labels[2])
if plot_labels[0]: ax_cube_z[0].plot(energy_grid_hartree_bulk, average_hartree_bulk[plot_folder], '-', color=plot_color[0], label=labels[0])
ax_cube_z[0].set_xlim([xlim[0], xlim[1]])
ax_cube_z[0].legend(frameon=False)
ax_cube_z[0].set_ylabel('Hartree potential / eV')
if plot_labels[1]: ax_cube_z[1].plot(energy_grid_hartree_em1, average_charge_dft[plot_folder], '-', color=plot_color[1], label=labels[1])
if plot_labels[3]: ax_cube_z[1].plot(energy_grid_hartree_em1, average_charge_em2[plot_folder], '-', color=plot_color[3], label=labels[3])
if plot_labels[2]: ax_cube_z[1].plot(energy_grid_hartree_em1, average_charge_em1[plot_folder], '-', color=plot_color[2], label=labels[2])
if plot_labels[0]: ax_cube_z[1].plot(energy_grid_hartree_bulk, average_charge_bulk[plot_folder], '-', color=plot_color[0], label=labels[0])
ax_cube_z[1].set_xlim([xlim[0], xlim[1]])
ax_cube_z[1].legend(frameon=False)
ax_cube_z[1].set_xlabel(r'Position / Å')
ax_cube_z[1].set_ylabel('Charge density')
fig_cube_z.tight_layout()
fig_cube_z.savefig('{}/charge_hartree_cube_{}.png'.format(folder_cp2k[plot_folder], print_label), dpi=300)
print('Finished plotting average')

# Plot Hartree and charge .cube difference
# if plot_diff:
#     fig_cube_both, ax_cube_both = plt.subplots(rows, cols, sharex='col', sharey='row', figsize=(6, 8))
#     # fig_cube_both, ax_cube_both = plt.subplots(rows, cols, sharex='col', sharey='row', figsize=(7, 8))
#     for i in range(len(folder_cp2k)):
#         ax_cube_both[0].plot(energy_grid_hartree_em1, average_hartree_em2[i]-average_hartree_em1[i], '-', color=diff_color[i], label=diff_label[i])
#     ax_cube_both[0].set_xlim([xlim[0], xlim[1]])
#     # ax_cube_both[0].set_xlabel(r'Position / Å')
#     ax_cube_both[0].set_ylabel('Hartree potential / eV')
#     for i in range(len(folder_cp2k)):
#         ax_cube_both[1].plot(energy_grid_hartree_em1, (average_charge_em2[i]-average_charge_em1[i])/np.abs(np.max(np.abs(average_charge_em2[i])-np.abs(average_charge_em1[i]))), '-', color=diff_color[i], label=diff_label[i])
#         # ax_cube_both[1].plot(energy_grid_hartree_em1, (average_charge_em2[i]-average_charge_em1[i])/np.max(average_charge_em2[1]-average_charge_em1[1]), '-', color=diff_color[i], label=diff_label[i])
#         # ax_cube_both[1].plot(energy_grid_hartree_em1, (average_charge_em2[i]-average_charge_em1[i]), '-', color=diff_color[i], label=diff_label[i])
#     if draw_markers: ax_cube_both[1].plot(markers, markers * 0, 'o', color='orange', fillstyle='none')
#     ax_cube_both[1].set_xlim([xlim[0], xlim[1]])
#     ax_cube_both[1].set_xlabel(r'Position / Å')
#     ax_cube_both[1].set_ylabel('Normalised charge density')
#     for i in range(len(folder_siesta)):
#         ax_cube_both[0].plot(data_hartree_3[i][:, 0], data_hartree_3[i][:, 1] - data_hartree_2[i][:, 1], '-',  label=diff_label_siesta[i], color=diff_color_siesta[i])
#         ax_cube_both[1].plot(data_charge_3[i][:, 0], (data_charge_3[i][:, 1] - data_charge_2[i][:, 1])/np.max(np.abs((data_charge_3[i][:, 1]) - (data_charge_2[i][:, 1]))), '-', label=diff_label_siesta[i], color=diff_color_siesta[i])
#         # ax_cube_both[1].plot(data_charge_3[i][:, 0], (data_charge_3[i][:, 1] - data_charge_2[i][:, 1])/np.max(data_charge_3[1][:, 1] - data_charge_2[1][:, 1]), '-', label=diff_label_siesta[i], color=diff_color_siesta[i])
#         # ax_cube_both[1].plot(data_charge_3[i][:, 0], (data_charge_3[i][:, 1] - data_charge_2[i][:, 1]), '-', label=diff_label_siesta[i], color=diff_color_siesta[i])
#     if plot_diff_legend[0]: ax_cube_both[0].legend(frameon=False)
#     if plot_diff_legend[1]: ax_cube_both[1].legend(frameon=False)
#     fig_cube_both.tight_layout()
#     for i in range(len(folder_cp2k)):
#         fig_cube_both.savefig('{}/charge_hartree_cube_diff_{}.png'.format(folder_cp2k[i], print_label), dpi=300)
#     print('Finished plotting difference Hartree and charge ')
#
# # Plot Hartree only
# if plot_diff:
#     fig_cube_hartree, ax_cube_hartree = plt.subplots()
#     for i in range(len(folder_cp2k)):
#         ax_cube_hartree.plot(energy_grid_hartree_em1, average_hartree_em2[i]-average_hartree_em1[i], '-', color=diff_color[i], label=diff_label[i])
#     for i in range(len(folder_siesta)):
#         ax_cube_hartree.plot(data_hartree_3[i][:, 0], data_hartree_3[i][:, 1] - data_hartree_2[i][:, 1], '-',  label=diff_label_siesta[i], color=diff_color_siesta[i])
#     ax_cube_hartree.set_xlim([xlim[0], xlim[1]])
#     ax_cube_hartree.set_xlabel(r'Position / Å')
#     ax_cube_hartree.set_ylabel('Hartree potential / eV')
#     if plot_diff_legend[0]: ax_cube_hartree.legend(frameon=False)
#     fig_cube_hartree.tight_layout()
#     for i in range(len(folder_cp2k)):
#         fig_cube_hartree.savefig('{}/hartree_cube_diff_{}.png'.format(folder_cp2k[i], print_label), dpi=300)
#     print('Finished plotting difference Hartree  ')
#
if __name__ == "__main__":
    print(folder_cp2k)
    print('Finished.')
    plt.show()
