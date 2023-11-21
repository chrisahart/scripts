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

# CP2K+SMEAGOL
# folder_cp2k = ['/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/au-capacitor/cp2k-smeagol/kpoints-2-2-20_V-4_double-contour']
# file_charge_em2 = ['4V-ELECTRON_DENSITY-1_0.cube', '4V-ELECTRON_DENSITY-1_0.cube', '4V-ELECTRON_DENSITY-1_0.cube']
# file_hartree_em2 = ['4V-v_hartree-1_0.cube', '4V-v_hartree-1_0.cube', '4V-v_hartree-1_0.cube']
folder_cp2k = ['/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/au-capacitor/cp2k-smeagol/kpoints-4-4-20_V-1_double-contour']
file_charge_em2 = ['1V-ELECTRON_DENSITY-1_0.cube']
file_hartree_em2 = ['1V-v_hartree-1_0.cube']
file_charge_em1 = ['0V-ELECTRON_DENSITY-1_0.cube'] * 4
file_hartree_em1 = ['0V-v_hartree-1_0.cube'] * 4
file_charge_dft = ['dft_wfn-ELECTRON_DENSITY-1_0.cube'] * 4
file_hartree_dft = ['dft_wfn-v_hartree-1_0.cube'] * 4
file_charge_bulk = ['bulk-ELECTRON_DENSITY-1_0.cube'] * 4
file_hartree_bulk = ['bulk-v_hartree-1_0.cube'] * 4
labels = ['CP2K+SMEAGOL bulk', 'CP2K', 'CP2K+SMEAGOL V=0', 'CP2K+SMEAGOL V=4']
plot_folder = 0
diff_label = ['CP2K+SMEAGOL']
plot_color = ['k', 'b', 'r', 'g']
read_labels = [True, True, True, True]
diff_color = ['r', 'g', 'b']
diff_color_siesta = ['g', 'b']
plot_labels = read_labels
print_label = 'average_z_all'
axis = 2
diff_color = ['r', 'g', 'b']
plot_diff_legend = [True, True]
draw_mirror = False
draw_mirror_diff = False
mirror_scale = -1
use_xlim = False
xlim_specify_left = [-0.1, 6]
xlim_specify_right = [78, 84.1]
xlim_specify = xlim_specify_left

# SIESTA
# folder_siesta = ['/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/au-capacitor/siesta1-smeagol/kpoints-2-2-20_V-4_double-contour']
# folder_siesta = ['/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/au-capacitor/siesta1-smeagol/kpoints-2-2-20_V-4_double-contour',
#                  '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/au-capacitor/siesta3-smeagol/kpoints-2-2-20_V-4_cores-20']
# charge_3 = ['0.4V-RHO_AV.dat', '0.4V-RHO_AV.dat', '0.4V-RHO_AV.dat']
# hartree_3 = ['0.4V-VH_AV.dat', '0.4V-VH_AV.dat', '0.4V-VH_AV.dat']
folder_siesta = ['/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/au-capacitor/siesta1-smeagol/kpoints-4-4-20_V-1_WeightRho-0.5']
charge_3 = ['0.1V-RHO_AV.dat']
hartree_3 = ['0.1V-VH_AV.dat']
charge_1 = ['0.bulk-RHO_AV.dat'] * 3
hartree_1 = ['0.bulk-VH_AV.dat'] * 3
charge_2 = ['0.0V-RHO_AV.dat'] * 3
hartree_2 = ['0.0V-VH_AV.dat'] * 3
labels_siesta = ['V=0', 'V=1', 'Bulk']
ylim = [-3, 2]
plot_dft = True
plot_markers = False
plot_diff = True
plot_diff_multiple = True
plot_mirror = False
# diff_label_siesta = ['SIESTA1+SMEAGOL', 'SIESTA3+SMEAGOL']
diff_label_siesta = ['SIESTA+SMEAGOL']
# diff_label = ['V=1 - V=0 EM.WeightRho 0.5', 'V=4 - V=0 EM.WeightRho 0.5', 'V=4 - V=0 weighted double contour']
x = 78.4 - 11.2
z_bond_length = 2.084
z_num = 6
z_right = 21.422
mid_pos = z_bond_length*(z_num - 1)+(z_right-z_bond_length*(z_num - 1))/2
markers = np.array(
    [z_bond_length * 0, z_bond_length * 1, z_bond_length * 2, z_bond_length * 3, z_bond_length * 4, z_bond_length * 5,
     z_right + z_bond_length * 0, z_right + z_bond_length * 1, z_right + z_bond_length * 2,
     z_right + z_bond_length * 3, z_right + z_bond_length * 4, z_right + z_bond_length * 5])
diff_color = ['r', 'g', 'b']
data_charge_1 = []
data_hartree_1 = []
data_hartree_3 = []
data_charge_2 = []
data_hartree_2 = []
data_charge_3 = []
data_cp2k_hartree = []
data_cp2k_charge = []
for i in range(len(folder_siesta)):
    data_charge_1.append(np.genfromtxt('{}/{}'.format(folder_siesta[i], charge_1[i]), skip_header=0, skip_footer=0))
    data_hartree_1.append(np.genfromtxt('{}/{}'.format(folder_siesta[i], hartree_1[i]), skip_header=0, skip_footer=0))
    if plot_dft: data_hartree_2.append(np.genfromtxt('{}/{}'.format(folder_siesta[i], hartree_2[i]), skip_header=0, skip_footer=0))
    if plot_dft: data_charge_2.append(np.genfromtxt('{}/{}'.format(folder_siesta[i], charge_2[i]), skip_header=0, skip_footer=0))
    data_hartree_3.append(np.genfromtxt('{}/{}'.format(folder_siesta[i], hartree_3[i]), skip_header=0, skip_footer=0))
    data_charge_3.append(np.genfromtxt('{}/{}'.format(folder_siesta[i], charge_3[i]), skip_header=0, skip_footer=0))

# Setup axis
mid_index = np.argmin(abs(data_hartree_3[0][:, 0]-mid_pos))
mid_pos_grid = data_hartree_3[0][mid_index, 0]

# Markers for Au capacitor
draw_markers = False
z_bond_length = 2.084
z_num = 6
z_right = 21.422
mid_pos = z_bond_length*(z_num - 1)+(z_right-z_bond_length*(z_num - 1))/2
markers = np.array(
    [z_bond_length * 0, z_bond_length * 1, z_bond_length * 2, z_bond_length * 3, z_bond_length * 4, z_bond_length * 5,
     z_right + z_bond_length * 0, z_right + z_bond_length * 1, z_right + z_bond_length * 2,
     z_right + z_bond_length * 3, z_right + z_bond_length * 4, z_right + z_bond_length * 5])

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
mid_index = np.argmin(abs(energy_grid_hartree_em1-mid_pos))
mid_pos_grid = energy_grid_hartree_em1[mid_index]

# Plot all
rows, cols = 2, 1
xlim = [0-1, data_charge_em1[0][1].get_cell()[axis][axis]+1]
fig_cube_z, ax_cube_z = plt.subplots(rows, cols, sharex='col', sharey='row', figsize=(6, 8))
if plot_fermi: ax_cube_z[0].axhline(y=fermi_dft, color='grey', linestyle='--', label='DFT Fermi energy', alpha=0.5)
if draw_mirror: ax_cube_z[0].plot(energy_grid_hartree_em1[:mid_index]+mid_pos_grid, np.flip(average_hartree_em2[plot_folder][:mid_index]), 'm--')
if draw_mirror: ax_cube_z[0].plot(energy_grid_hartree_em1[:mid_index]+mid_pos_grid, np.flip(average_hartree_em1[plot_folder][:mid_index]), 'm--')
if plot_labels[1]: ax_cube_z[0].plot(energy_grid_hartree_em1, average_hartree_dft[plot_folder], '-', color=plot_color[1], label=labels[1])
if plot_labels[3]: ax_cube_z[0].plot(energy_grid_hartree_em1, average_hartree_em2[plot_folder], '-', color=plot_color[3], label=labels[3])
if plot_labels[2]: ax_cube_z[0].plot(energy_grid_hartree_em1, average_hartree_em1[plot_folder], '-', color=plot_color[2], label=labels[2])
if plot_labels[0]: ax_cube_z[0].plot(energy_grid_hartree_bulk, average_hartree_bulk[plot_folder], '-', color=plot_color[0], label=labels[0])
if draw_markers: ax_cube_z[0].plot(markers, markers*0, 'o', color='orange', fillstyle='none')
ax_cube_z[0].set_xlim([xlim[0], xlim[1]])
ax_cube_z[0].legend(frameon=False)
ax_cube_z[0].set_ylabel('Hartree potential / eV')
if draw_mirror: ax_cube_z[1].plot(energy_grid_hartree_em1[:mid_index]+mid_pos_grid, np.flip(average_charge_em2[plot_folder][:mid_index]), 'm--')
if draw_mirror: ax_cube_z[1].plot(energy_grid_hartree_em1[:mid_index]+mid_pos_grid, np.flip(average_charge_em1[plot_folder][:mid_index]), 'm--')
if plot_labels[1]: ax_cube_z[1].plot(energy_grid_hartree_em1, average_charge_dft[plot_folder], '-', color=plot_color[1], label=labels[1])
if plot_labels[3]: ax_cube_z[1].plot(energy_grid_hartree_em1, average_charge_em2[plot_folder], '-', color=plot_color[3], label=labels[3])
if plot_labels[2]: ax_cube_z[1].plot(energy_grid_hartree_em1, average_charge_em1[plot_folder], '-', color=plot_color[2], label=labels[2])
if plot_labels[0]: ax_cube_z[1].plot(energy_grid_hartree_bulk, average_charge_bulk[plot_folder], '-', color=plot_color[0], label=labels[0])
if draw_markers: ax_cube_z[1].plot(markers, markers*0, 'o', color='orange', fillstyle='none')
ax_cube_z[1].set_xlim([xlim[0], xlim[1]])
ax_cube_z[1].legend(frameon=False)
ax_cube_z[1].set_xlabel(r'Position / Å')
ax_cube_z[1].set_ylabel('Charge density')
fig_cube_z.tight_layout()
fig_cube_z.savefig('{}/charge_hartree_cube_{}.png'.format(folder_cp2k[plot_folder], print_label), dpi=300)
print('Finished plotting average')

# Plot Hartree and charge .cube difference
if plot_diff:
    # fig_cube_both, ax_cube_both = plt.subplots(rows, cols, sharex='col', sharey='row', figsize=(6, 8))
    fig_cube_both, ax_cube_both = plt.subplots(rows, cols, sharex='col', sharey='row', figsize=(7, 8))
    for i in range(len(folder_cp2k)):
        ax_cube_both[0].plot(energy_grid_hartree_em1, average_hartree_em2[i]-average_hartree_em1[i], '-', color=diff_color[i], label=diff_label[i])
    ax_cube_both[0].set_xlim([xlim[0], xlim[1]])
    # ax_cube_both[0].set_xlabel(r'Position / Å')
    ax_cube_both[0].set_ylabel('Hartree potential / eV')
    for i in range(len(folder_cp2k)):
        ax_cube_both[1].plot(energy_grid_hartree_em1, (average_charge_em2[i]-average_charge_em1[i])/np.max(average_charge_em2[i]-average_charge_em1[i]), '-', color=diff_color[i], label=diff_label[i])
    if draw_markers: ax_cube_both[1].plot(markers, markers * 0, 'o', color='orange', fillstyle='none')
    ax_cube_both[1].set_xlim([xlim[0], xlim[1]])
    ax_cube_both[1].set_xlabel(r'Position / Å')
    ax_cube_both[1].set_ylabel('Normalised charge density')
    for i in range(len(folder_siesta)):
        ax_cube_both[0].plot(data_hartree_3[i][:, 0], data_hartree_3[i][:, 1] - data_hartree_2[i][:, 1], '-',  label=diff_label_siesta[i], color=diff_color_siesta[i])
        ax_cube_both[1].plot(data_charge_3[i][:, 0], (data_charge_3[i][:, 1] - data_charge_2[i][:, 1])/np.max(data_charge_3[i][:, 1] - data_charge_2[i][:, 1]), '-', label=diff_label_siesta[i], color=diff_color_siesta[i])
    if plot_diff_legend[0]: ax_cube_both[0].legend(frameon=False)
    if plot_diff_legend[1]: ax_cube_both[1].legend(frameon=False)
    fig_cube_both.tight_layout()
    for i in range(len(folder_cp2k)):
        fig_cube_both.savefig('{}/charge_hartree_cube_diff_{}.png'.format(folder_cp2k[i], print_label), dpi=300)
    print('Finished plotting difference Hartree and charge ')

if __name__ == "__main__":
    print(folder_cp2k)
    print('Finished.')
    plt.show()
