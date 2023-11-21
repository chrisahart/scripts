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
print_label = 'V0_DFT'

# cp2k-smeagol-examples/examples/au-capacitor
# folder1 = [
#     '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/au-capacitor/cp2k-smeagol/kpoints-2-2-20_V-4_double-contour']
# folder1 = ['/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/au-capacitor/cp2k-smeagol/kpoints-2-2-20_V-1_WeightRho-0.5',
#            '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/au-capacitor/cp2k-smeagol/kpoints-2-2-20_V-4_WeightRho-0.5',
#            '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/au-capacitor/cp2k-smeagol/kpoints-2-2-20_V-4_double-contour']
# file_charge_em2 = ['4V-ELECTRON_DENSITY-1_0.cube', '4V-ELECTRON_DENSITY-1_0.cube', '4V-ELECTRON_DENSITY-1_0.cube']
# file_hartree_em2 = ['4V-v_hartree-1_0.cube', '4V-v_hartree-1_0.cube', '4V-v_hartree-1_0.cube']
# file_charge_em1 = ['0V-ELECTRON_DENSITY-1_0.cube'] * 4
# file_hartree_em1 = ['0V-v_hartree-1_0.cube'] * 4
# file_charge_dft = ['dft_wfn-ELECTRON_DENSITY-1_0.cube'] * 4
# file_hartree_dft = ['dft_wfn-v_hartree-1_0.cube'] * 4
# file_charge_bulk = ['bulk-ELECTRON_DENSITY-1_0.cube'] * 4
# file_hartree_bulk = ['bulk-v_hartree-1_0.cube'] * 4
# labels = ['CP2K+SMEAGOL bulk', 'CP2K', 'CP2K+SMEAGOL V=0', 'CP2K+SMEAGOL V=1']
# plot_folder = 0
# diff_label = ['V=1 - V=0 EM.WeightRho 0.5', 'V=4 - V=0 EM.WeightRho 0.5', 'V=4 - V=0 weighted double contour']
# plot_color = ['k', 'b', 'r', 'g']
# read_labels = [True, True, True, True]
# diff_color = ['r', 'g', 'b']
# plot_labels = read_labels
# print_label = 'average_z_all'
# axis = 2
# diff_color = ['r', 'g', 'b']
# plot_diff_legend = [True, False]
# draw_mirror = False
# draw_mirror_diff = True
# mirror_scale = -1
# use_xlim = False
# xlim_specify_left = [-0.1, 6]
# xlim_specify_right = [78, 84.1]
# xlim_specify = xlim_specify_left
# draw_markers = False
# z_bond_length = 2.084
# z_num = 6
# z_right = 21.422
# mid_pos = z_bond_length*(z_num - 1)+(z_right-z_bond_length*(z_num - 1))/2
# print(mid_pos)
# markers = np.array(
#     [z_bond_length * 0, z_bond_length * 1, z_bond_length * 2, z_bond_length * 3, z_bond_length * 4, z_bond_length * 5,
#      z_right + z_bond_length * 0, z_right + z_bond_length * 1, z_right + z_bond_length * 2,
#      z_right + z_bond_length * 3, z_right + z_bond_length * 4, z_right + z_bond_length * 5])

# cp2k-smeagol-examples/examples/au-h2
# folder1 = ['/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/au-h2/cp2k-smeagol/transmission/kpoints-3-3-20_omp-2_ParallelOverKNum-9_contour-single_EM.WeightRho-0',
#            '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/au-h2/cp2k-smeagol/transmission/kpoints-3-3-20_omp-2_ParallelOverKNum-9_contour-single_EM.WeightRho-1',
#            '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/au-h2/cp2k-smeagol/transmission/kpoints-3-3-20_omp-2_ParallelOverKNum-9_contour-double']
# folder1 = ['/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/au-h2/cp2k-smeagol/transmission/kpoints-3-3-20_omp-2_ParallelOverKNum-9_contour-single',
           # '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/au-h2/cp2k-smeagol/transmission/kpoints-3-3-20_omp-2_ParallelOverKNum-9_contour-double']
# file_charge_em2 = ['1.9V-ELECTRON_DENSITY-1_0.cube', '1.9V-ELECTRON_DENSITY-1_0.cube', '1.9V-ELECTRON_DENSITY-1_0.cube', '1.9V-ELECTRON_DENSITY-1_0.cube']
# file_hartree_em2 = ['1.9V-v_hartree-1_0.cube', '1.9V-v_hartree-1_0.cube', '1.9V-v_hartree-1_0.cube', '1.9V-v_hartree-1_0.cube']
# folder1 = ['/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/au-h2/cp2k-smeagol/transmission/kpoints-3-3-20_omp-2_ParallelOverKNum-9_contour-single_V-0.9',
#            '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/au-h2/cp2k-smeagol/transmission/kpoints-3-3-20_omp-2_ParallelOverKNum-9_contour-single_V-1.1',
#            '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/au-h2/cp2k-smeagol/transmission/kpoints-3-3-20_omp-2_ParallelOverKNum-9_contour-single']
# file_charge_em2 = ['0.9V-ELECTRON_DENSITY-1_0.cube', '1.1V-ELECTRON_DENSITY-1_0.cube', '1.9V-ELECTRON_DENSITY-1_0.cube', '1.9V-ELECTRON_DENSITY-1_0.cube']
# file_hartree_em2 = ['0.9V-v_hartree-1_0.cube', '1.1V-v_hartree-1_0.cube', '1.9V-v_hartree-1_0.cube', '1.9V-v_hartree-1_0.cube']
# folder1 = [
#     '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/au-h2/cp2k-smeagol/geo_opt/kpoints-3-3-20_omp-2_ParallelOverKNum-3_contour-single/1.3',
#     '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/au-h2/cp2k-smeagol/geo_opt/kpoints-3-3-20_omp-2_ParallelOverKNum-3_contour-single/1.5',
#     '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/au-h2/cp2k-smeagol/geo_opt/kpoints-3-3-20_omp-2_ParallelOverKNum-3_contour-single/1.7',
#     '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/au-h2/cp2k-smeagol/geo_opt/kpoints-3-3-20_omp-2_ParallelOverKNum-3_contour-single/1.9']
# file_charge_em2 = ['1.3V-ELECTRON_DENSITY-1_0.cube', '1.5V-ELECTRON_DENSITY-1_0.cube',  '1.7V-ELECTRON_DENSITY-1_0.cube', '1.9V-ELECTRON_DENSITY-1_0.cube']
# file_hartree_em2 = ['1.3V-v_hartree-1_0.cube', '1.5V-v_hartree-1_0.cube', '1.7V-v_hartree-1_0.cube', '1.9V-v_hartree-1_0.cube']
# folder1 = ['/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/au-h2/cp2k-smeagol/transmission_all/kpoints-3-3-20_omp-2_ParallelOverKNum-9_contour-single/kpoints-3-3-20_omp-2_ParallelOverKNum-9_contour-single/0.1']
# file_charge_em2 = ['0.1V-ELECTRON_DENSITY-1_0.cube', '1.9V-ELECTRON_DENSITY-1_0.cube', '1.9V-ELECTRON_DENSITY-1_0.cube', '1.9V-ELECTRON_DENSITY-1_0.cube']
# file_hartree_em2 = ['0.1V-v_hartree-1_0.cube', '1.9V-v_hartree-1_0.cube', '1.9V-v_hartree-1_0.cube', '1.9V-v_hartree-1_0.cube']
folder1 = [
    '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/au-h2/cp2k-smeagol/geo_opt/kpoints-3-3-20_omp-2_ParallelOverKNum-3_contour-single_dynamic-14/1.5'
    ]
file_charge_em2 = ['1.5V-ELECTRON_DENSITY-1_41.cube',  '1.7V-ELECTRON_DENSITY-1_0.cube', '1.9V-ELECTRON_DENSITY-1_0.cube']
file_hartree_em2 = ['1.5V-v_hartree-1_41.cube', '1.7V-v_hartree-1_0.cube', '1.9V-v_hartree-1_0.cube']
file_charge_em1 = ['0V-ELECTRON_DENSITY-1_0.cube'] * 10
file_hartree_em1 = ['0V-v_hartree-1_0.cube'] * 10
# file_charge_em1 = ['0.0V-ELECTRON_DENSITY-1_0.cube'] * 10
# file_hartree_em1 = ['0.0V-v_hartree-1_0.cube'] * 10
file_charge_dft = ['dft_wfn-ELECTRON_DENSITY-1_0.cube'] * 10
file_hartree_dft = ['dft_wfn-v_hartree-1_0.cube'] * 10
file_charge_bulk = ['bulk-ELECTRON_DENSITY-1_0.cube'] * 10
file_hartree_bulk = ['bulk-v_hartree-1_0.cube'] * 10
labels = ['CP2K+SMEAGOL bulk', 'CP2K', 'CP2K+SMEAGOL V=0', 'CP2K+SMEAGOL V=1']
plot_folder = 0
# diff_label = ['V=1.9 - V=0 EM.WeightRho 0.5', 'V=1.9 - V=0 weighted double contour']
# diff_label = ['CP2K+SMEAGOL EM.WeightRho 0.5', 'CP2K+SMEAGOL weighted double contour']
diff_label = ['V=1.3 - V=0', 'V=1.5 - V=0', 'V=1.7 - V=0', 'V=1.9 - V=0']
plot_color = ['k', 'b', 'r', 'g', 'm', 'grey', 'orange', 'y', 'brown', 'cyan', 'pink']
read_labels = [False, False, True, True]
plot_labels = read_labels
plot_diff = True
print_label = 'average_z_all'
axis = 2
diff_color = ['r', 'g', 'b', 'm', 'grey', 'orange', 'y', 'brown', 'cyan', 'pink']
plot_diff_legend = [True, False]
draw_mirror = False
draw_mirror_diff = True
mirror_scale = -1
use_xlim = False
xlim_specify_left = [-0.1, 6]
xlim_specify_right = [78, 84.1]
xlim_specify = xlim_specify_left
draw_markers = False
mid_pos = 19.3615

# Read .cube using ASE
data_charge_em2 = []
data_hartree_em2 = []
data_charge_em1 = []
data_hartree_em1 = []
data_charge_bulk = []
data_charge_dft = []
data_hartree_dft = []
data_hartree_bulk = []
for i in range(len(folder1)):
    print('reading .cube files', folder1[i])
    if read_labels[3]: data_charge_em2.append(read_cube_data('{}/{}'.format(folder1[i], file_charge_em2[i])))
    if read_labels[3]: data_hartree_em2.append(read_cube_data('{}/{}'.format(folder1[i], file_hartree_em2[i])))
    if read_labels[2]: data_charge_em1.append(read_cube_data('{}/{}'.format(folder1[i], file_charge_em1[i])))
    if read_labels[2]: data_hartree_em1.append(read_cube_data('{}/{}'.format(folder1[i], file_hartree_em1[i])))
    if read_labels[1]: data_charge_dft.append(read_cube_data('{}/{}'.format(folder1[i], file_charge_dft[i])))
    if read_labels[1]: data_hartree_dft.append(read_cube_data('{}/{}'.format(folder1[i], file_hartree_dft[i])))
    if read_labels[0]: data_charge_bulk.append(read_cube_data('{}/{}'.format(folder1[i], file_charge_bulk[i])))
    if read_labels[0]: data_hartree_bulk.append(read_cube_data('{}/{}'.format(folder1[i], file_hartree_bulk[i])))
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
for i in range(len(folder1)):
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
for i in range(len(folder1)):
    print('size em1 folder', i,  data_charge_em1[i][1].get_cell()[axis][axis])
    print('size em2 folder', i,  data_charge_em2[i][1].get_cell()[axis][axis])
    print('size em1 folder', i, average_charge_em1[0].shape[0])
    print('size em2 folder', i, average_charge_em2[0].shape[0])
if read_labels[3]: energy_grid_hartree_em2 = np.linspace(start=0, stop=data_charge_em2[0][1].get_cell()[axis][axis], num=average_charge_em2[0].shape[0])
if read_labels[2]: energy_grid_hartree_em1 = np.linspace(start=0, stop=data_charge_em1[0][1].get_cell()[axis][axis], num=average_charge_em1[0].shape[0])
if read_labels[1]: energy_grid_hartree_dft = np.linspace(start=0, stop=data_charge_dft[0][1].get_cell()[axis][axis], num=average_charge_dft[0].shape[0])
if read_labels[0]: energy_grid_hartree_bulk = np.linspace(start=0, stop=data_charge_bulk[0][1].get_cell()[axis][axis], num=average_charge_bulk[0].shape[0])
mid_index = np.argmin(abs(energy_grid_hartree_em1-mid_pos))
mid_pos_grid = energy_grid_hartree_em1[mid_index]

# Plot all
rows, cols = 2, 1
xlim = [0-1, data_charge_em1[0][1].get_cell()[axis][axis]+1]
fig_cube_z, ax_cube_z = plt.subplots(rows, cols, sharex='col', sharey='row', figsize=(6, 8))
if plot_fermi: ax_cube_z[0].axhline(y=fermi_dft, color='grey', linestyle='--', label='DFT Fermi energy', alpha=0.5)
if draw_mirror: ax_cube_z[0].plot(energy_grid_hartree_em1[:mid_index]+mid_pos_grid, np.flip(average_hartree_em2[plot_folder][:mid_index]), 'm--')
if draw_mirror: ax_cube_z[0].plot(energy_grid_hartree_em1[:mid_index]+mid_pos_grid, np.flip(average_hartree_em1[plot_folder][:mid_index]), 'm--')
if plot_labels[1]: ax_cube_z[0].plot(energy_grid_hartree_dft, average_hartree_dft[plot_folder], '-', color=plot_color[1], label=labels[1])
if plot_labels[3]: ax_cube_z[0].plot(energy_grid_hartree_em2, average_hartree_em2[plot_folder], '-', color=plot_color[3], label=labels[3])
if plot_labels[2]: ax_cube_z[0].plot(energy_grid_hartree_em1, average_hartree_em1[plot_folder], '-', color=plot_color[2], label=labels[2])
if plot_labels[0]: ax_cube_z[0].plot(energy_grid_hartree_bulk, average_hartree_bulk[plot_folder], '-', color=plot_color[0], label=labels[0])
if draw_markers: ax_cube_z[0].plot(markers, markers*0, 'o', color='orange', fillstyle='none')
ax_cube_z[0].set_xlim([xlim[0], xlim[1]])
ax_cube_z[0].legend(frameon=False)
ax_cube_z[0].set_ylabel('Hartree potential / eV')
if draw_mirror: ax_cube_z[1].plot(energy_grid_hartree_em1[:mid_index]+mid_pos_grid, np.flip(average_charge_em2[plot_folder][:mid_index]), 'm--')
if draw_mirror: ax_cube_z[1].plot(energy_grid_hartree_em1[:mid_index]+mid_pos_grid, np.flip(average_charge_em1[plot_folder][:mid_index]), 'm--')
if plot_labels[1]: ax_cube_z[1].plot(energy_grid_hartree_dft, average_charge_dft[plot_folder], '-', color=plot_color[1], label=labels[1])
if plot_labels[3]: ax_cube_z[1].plot(energy_grid_hartree_em2, average_charge_em2[plot_folder], '-', color=plot_color[3], label=labels[3])
if plot_labels[2]: ax_cube_z[1].plot(energy_grid_hartree_em1, average_charge_em1[plot_folder], '-', color=plot_color[2], label=labels[2])
if plot_labels[0]: ax_cube_z[1].plot(energy_grid_hartree_bulk, average_charge_bulk[plot_folder], '-', color=plot_color[0], label=labels[0])
if draw_markers: ax_cube_z[1].plot(markers, markers*0, 'o', color='orange', fillstyle='none')
ax_cube_z[1].set_xlim([xlim[0], xlim[1]])
ax_cube_z[1].legend(frameon=False)
ax_cube_z[1].set_xlabel(r'Position / Å')
ax_cube_z[1].set_ylabel('Charge density')
fig_cube_z.tight_layout()
fig_cube_z.savefig('{}/charge_hartree_cube_{}.png'.format(folder1[plot_folder], print_label), dpi=300)
print('Finished plotting average')

# Plot Hartree and charge .cube difference
if plot_diff:
    fig_cube_both, ax_cube_both = plt.subplots(rows, cols, sharex='col', sharey='row', figsize=(6, 8))
    for i in range(len(folder1)):
        ax_cube_both[0].plot(energy_grid_hartree_em1, average_hartree_em2[i]-average_hartree_em1[i], '-', color=diff_color[i], label=diff_label[i])
        if draw_mirror_diff: ax_cube_both[0].plot(energy_grid_hartree_em1[:mid_index] + mid_pos_grid, -1*np.flip(average_hartree_em2[i][:mid_index]-average_hartree_em1[i][:mid_index]), '--', color=diff_color[i])
    if draw_markers: ax_cube_both[0].plot(markers, markers * 0, 'o', color='orange', fillstyle='none')
    if plot_diff_legend[0]: ax_cube_both[0].legend(frameon=False)
    ax_cube_both[0].set_xlim([xlim[0], xlim[1]])
    ax_cube_both[0].set_xlabel(r'Position / Å')
    ax_cube_both[0].set_ylabel('Hartree potential z / eV')
    for i in range(len(folder1)):
        ax_cube_both[1].plot(energy_grid_hartree_em1, average_charge_em2[i]-average_charge_em1[i], '-', color=diff_color[i], label=diff_label[i])
        if draw_mirror_diff: ax_cube_both[1].plot(energy_grid_hartree_em1[:mid_index] + mid_pos_grid, -1*np.flip(average_charge_em2[i][:mid_index]-average_charge_em1[i][:mid_index]), '--', color=diff_color[i])
    if draw_markers: ax_cube_both[1].plot(markers, markers * 0, 'o', color='orange', fillstyle='none')
    if plot_diff_legend[1]: ax_cube_both[1].legend(frameon=False)
    ax_cube_both[1].set_xlim([xlim[0], xlim[1]])
    ax_cube_both[1].set_xlabel(r'Position / Å')
    ax_cube_both[1].set_ylabel('Charge density z')
    fig_cube_both.tight_layout()
    for i in range(len(folder1)):
        fig_cube_both.savefig('{}/charge_hartree_cube_diff_{}.png'.format(folder1[i], print_label), dpi=300)
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
fig_hartree, ax_hartree = plt.subplots()
if plot_diff:
    for i in range(len(folder1)):
        ax_hartree.plot(energy_grid_hartree_em1, average_hartree_em2[i] - average_hartree_em1[i], '-', color=diff_color[i], label=diff_label[i])
        if draw_mirror_diff: ax_hartree.plot(energy_grid_hartree_em1[:mid_index] + mid_pos_grid, -1 * np.flip(average_hartree_em2[i][:mid_index] - average_hartree_em1[i][:mid_index]), '--', color=diff_color[i])
ax_hartree.set_xlim([xlim[0], xlim[1]])
ax_hartree.legend(frameon=False)
ax_hartree.set_xlabel(r'Position / Å')
ax_hartree.set_ylabel('Hartree potential z / eV')
fig_hartree.tight_layout()
for i in range(len(folder1)):
    fig_hartree.savefig('{}/hartree_cube_diff_{}.png'.format(folder1[i], print_label), dpi=300)

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
    print(folder1)
    print('Finished.')
    plt.show()
