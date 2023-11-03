import numpy as np
from matplotlib import pyplot as plt
from ase.io.cube import read_cube_data
from general import parameters as param

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


# CP2K+SMEAGOL testing
# folder1 = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-capacitor/layers-1-2-3-4/cp2k/weightrho/weightrho-0/kpoints-4-4-20_hlb-auto'
# folder1 = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-capacitor/layers-1-2-3-4/cp2k/weightrho/weightrho-0.5/kpoints-4-4-20_hlb-auto'
# folder1 = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-capacitor/layers-1-2-3-4/cp2k/weightrho/weightrho-1/kpoints-4-4-20_hlb-auto'
# folder1 = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-capacitor/layers-1-2-3-4/cp2k/weightrho/weightrho-default/kpoints-4-4-20_hlb-auto'
# folder1 = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-capacitor/layers-1-2-3-4/cp2k/weightrho/testing/kpoints-1-1-20_hlb-auto_NEnergReal-64_supercell-2-2-1'
# folder1 = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-capacitor/layers-1-2-3-4/cp2k/weightrho/weightrho/testing/li/kpoints-4-4-20_hlb-auto_cores-20_hartree_em10220907-memory_hash-143_OrderN-T'
# folder1 = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-capacitor/layers-1-2-3-4/cp2k/weightrho/weightrho/testing/li/kpoints-8-8-20_hlb-auto_cores-20_hartree_em10220907-memory_hash-143_OrderN-T_TransmissionOverk_em1.TRCChannels'
# folder1 = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-capacitor/layers-1-2-3-4/cp2k/cp2k-versions/kpoints-4-4-20_hlb-auto_cores-20_hartree_em10220907-memory_hash-143_OrderN-T_em1.WeightRho-0'
# folder1 = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-capacitor/layers-1-2-3-4/cp2k/weightrho/testing/kpoints-6-6-20_hlb-auto_NEnergReal-64'
# folder1 = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-capacitor/layers-1-2-3-4/cp2k/weightrho/testing/kpoints-1-1-20_hlb-auto_NEnergReal-64_supercell-2-2-1'
# folder1 = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-capacitor/layers-1-2-3-4/cp2k/weightrho/testing/kpoints-4-4-20_hlb-auto_NEnergReal-64_pt-sz'
# folder1 = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-capacitor/layers-1-2-3-4/cp2k/weightrho/testing/kpoints-8-8-20_hlb-auto_NEnergReal-64'
# folder1 = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-capacitor/layers-1-2-3-4/cp2k/weightrho/weightrho-1/kpoints-4-4-20_hlb-auto'
# folder1 = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-capacitor/layers-1-2-3-4/cp2k/weightrho/testing/kpoints-4-4-20_hlb-auto_NEnergReal-64_li_chris'
folder1 = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-capacitor/layers-1-2-3-4/cp2k/cp2k-versions/kpoints-4-4-20_hlb-auto_NEnergReal-64_double-contour'
# folder1 = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-capacitor/layers-1-2-3-4/cp2k/weightrho/testing/kpoints-4-4-20_hlb-auto_NEnergReal-64_li_chris_5V'
# folder1 = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-capacitor/layers-1-2-3-4/cp2k/weightrho/weightrho-1/kpoints-2-2-20_hlb-auto'
# folder1 = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-capacitor/layers-1-2-3-4/cp2k/weightrho/testing/kpoints-2-2-20_hlb-auto_NEnergReal-64'
folder2 = folder1
# folder2 = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-capacitor/layers-1-2-3-4/cp2k/weightrho/weightrho-1/scf-steps/kpoints-4-4-20_hlb-auto_dft-rs-none_V-1/max_scf/scf_50'
# folder2 = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-capacitor/layers-1-2-3-4/cp2k/weightrho/weightrho-1/scf-steps/kpoints-4-4-20_hlb-auto_dft-rs-none_V-1/max_scf/scf_1'
# folder2 = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-capacitor/layers-1-2-3-4/cp2k/weightrho/weightrho-1/scf-steps/kpoints-4-4-20_hlb-auto_V0-rs-all_V-1/max_scf/scf_1'
# folder2 = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-capacitor/layers-1-2-3-4/cp2k/weightrho/weightrho-1/scf-steps/kpoints-2-2-20_hlb-auto_V0-rs-all_V-1/max_scf/scf_40'
# folder2 = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-capacitor/layers-1-2-3-4/cp2k/weightrho/weightrho-1/scf-steps/kpoints-2-2-20_hlb-auto_dft-rs-none_V-1/max_scf/scf_1'
print('folder1', folder1)
print('folder2', folder2)
folder_save = '{}'.format(folder1)
# file_charge_em2 = '5V-ELECTRON_DENSITY-1_0.cube'
# file_hartree_em2 = '5V-v_hartree-1_0.cube'
file_charge_em2 = '1V-ELECTRON_DENSITY-1_0.cube'
file_hartree_em2 = '1V-v_hartree-1_0.cube'
file_charge_em1 = '0V-ELECTRON_DENSITY-1_0.cube'
file_hartree_em1 = '0V-v_hartree-1_0.cube'
file_charge_dft = 'dft_wfn-ELECTRON_DENSITY-1_0.cube'
file_hartree_dft = 'dft_wfn-v_hartree-1_0.cube'
file_charge_bulk = 'bulk-ELECTRON_DENSITY-1_0.cube'
file_hartree_bulk = 'bulk-v_hartree-1_0.cube'
labels = ['Bulk.cube', 'DFT.cube', 'V=0.cube', 'V=1.cube']
plot_color = ['k', 'b', 'r', 'g']
read_labels = [True, False, True, True]
plot_labels = read_labels
print_label = 'average_z_all'
axis = 2
plot_diff_em1_em2 = True
plot_diff_dft_em2 = False
plot_diff_dft_em1 = False
diff_color = ['r', 'g', 'k']
plot_diff_legend = False
draw_mirror_em2 = False
draw_mirror_em1 = False
mirror_scale = -1
use_xlim = False
xlim_specify_left = [-0.1, 6]
xlim_specify_right = [78, 84.1]
xlim_specify = xlim_specify_left
draw_markers = True
z_bond_length = 2.084
z_num = 6
z_right = 21.422
mid_pos = z_bond_length*(z_num - 1)+(z_right-z_bond_length*(z_num - 1))/2
markers = np.array(
    [z_bond_length * 0, z_bond_length * 1, z_bond_length * 2, z_bond_length * 3, z_bond_length * 4, z_bond_length * 5,
     z_right + z_bond_length * 0, z_right + z_bond_length * 1, z_right + z_bond_length * 2,
     z_right + z_bond_length * 3, z_right + z_bond_length * 4, z_right + z_bond_length * 5])

# Read .cube using ASE
if read_labels[3]: data_charge_em2, atoms_charge_em2 = read_cube_data('{}/{}'.format(folder2, file_charge_em2))
if read_labels[3]: data_hartree_em2, atoms_hartree_em2 = read_cube_data('{}/{}'.format(folder2, file_hartree_em2))
if read_labels[2]: data_charge_em1, atoms_charge_em1 = read_cube_data('{}/{}'.format(folder1, file_charge_em1))
if read_labels[2]: data_hartree_em1, atoms_hartree_em1 = read_cube_data('{}/{}'.format(folder1, file_hartree_em1))
if read_labels[1]: data_charge_dft, atoms_charge_dft = read_cube_data('{}/{}'.format(folder1, file_charge_dft))
if read_labels[1]: data_hartree_dft, atoms_hartree_dft = read_cube_data('{}/{}'.format(folder1, file_hartree_dft))
if read_labels[0]: data_charge_bulk, atoms_charge_bulk = read_cube_data('{}/{}'.format(folder1, file_charge_bulk))
if read_labels[0]: data_hartree_bulk, atoms_hartree_bulk = read_cube_data('{}/{}'.format(folder1, file_hartree_bulk))
print('Finished reading .cube files')

# Calculate average along axis
# if read_labels[3]: print(data_charge_em2.shape[axis])
# if read_labels[2]: print(data_charge_em1.shape[axis])
# if read_labels[1]: print(data_charge_dft.shape[axis])
# if read_labels[0]: print(data_charge_bulk.shape[axis])
# if read_labels[3]: print(data_hartree_em2.shape[axis])
# if read_labels[2]: print(data_hartree_em1.shape[axis])
# if read_labels[1]: print(data_hartree_dft.shape[axis])
# if read_labels[0]: print(data_hartree_bulk.shape[axis])
if read_labels[3]: average_charge_em2 = np.zeros(data_charge_em2.shape[axis])
if read_labels[3]: average_hartree_em2 = np.zeros(data_hartree_em2.shape[axis])
if read_labels[2]: average_charge_em1 = np.zeros(data_charge_em1.shape[axis])
if read_labels[2]: average_hartree_em1 = np.zeros(data_hartree_em1.shape[axis])
if read_labels[1]: average_charge_dft = np.zeros(data_charge_dft.shape[axis])
if read_labels[1]: average_hartree_dft = np.zeros(data_hartree_dft.shape[axis])
if read_labels[0]: average_charge_bulk = np.zeros(data_charge_bulk.shape[axis])
if read_labels[0]: average_hartree_bulk = np.zeros(data_hartree_bulk.shape[axis])
for i in range(data_charge_em1.shape[axis]):
    if read_labels[3]: average_charge_em2[i] = np.mean(data_charge_em2[:, :, i])
    if read_labels[3]: average_hartree_em2[i] = np.mean(data_hartree_em2[:, :, i] * param.hartree_to_ev)
    if read_labels[2]: average_charge_em1[i] = np.mean(data_charge_em1[:, :, i])
    if read_labels[2]: average_hartree_em1[i] = np.mean(data_hartree_em1[:, :, i] * param.hartree_to_ev)
    if read_labels[1]: average_charge_dft[i] = np.mean(data_charge_dft[:, :, i])
    if read_labels[1]: average_hartree_dft[i] = np.mean(data_hartree_dft[:, :, i] * param.hartree_to_ev)
for i in range(data_charge_bulk.shape[axis]):
    if read_labels[0]: average_charge_bulk[i] = np.mean(data_charge_bulk[:, :, i])
    if read_labels[0]: average_hartree_bulk[i] = np.mean(data_hartree_bulk[:, :, i] * param.hartree_to_ev)
if read_labels[3]: energgrid_charge_em2 = np.linspace(start=0, stop=atoms_charge_em2.get_cell()[axis][axis], num=data_charge_em2.shape[axis])
if read_labels[3]: energgrid_hartree_em2 = np.linspace(start=0, stop=atoms_charge_em2.get_cell()[axis][axis], num=data_hartree_em2.shape[axis])
if read_labels[2]: energgrid_charge_em1 = np.linspace(start=0, stop=atoms_charge_em1.get_cell()[axis][axis], num=data_charge_em1.shape[axis])
if read_labels[2]: energgrid_hartree_em1 = np.linspace(start=0, stop=atoms_charge_em1.get_cell()[axis][axis], num=data_hartree_em1.shape[axis])
if read_labels[1]: energgrid_charge_dft = np.linspace(start=0, stop=atoms_charge_dft.get_cell()[axis][axis], num=data_charge_dft.shape[axis])
if read_labels[1]: energgrid_hartree_dft = np.linspace(start=0, stop=atoms_charge_dft.get_cell()[axis][axis], num=data_hartree_dft.shape[axis])
if read_labels[0]: energgrid_charge_bulk = np.linspace(start=0, stop=atoms_charge_bulk.get_cell()[axis][axis], num=data_charge_bulk.shape[axis])
if read_labels[0]: energgrid_hartree_bulk = np.linspace(start=0, stop=atoms_charge_bulk.get_cell()[axis][axis], num=data_hartree_bulk.shape[axis])

# Setup axis
if read_labels[3]: energy_grid_charge_em2 = np.linspace(start=0, stop=atoms_charge_em2.get_cell()[axis][axis], num=data_charge_em2.shape[axis])
if read_labels[3]: energy_grid_hartree_em2 = np.linspace(start=0, stop=atoms_charge_em2.get_cell()[axis][axis], num=data_hartree_em2.shape[axis])
if read_labels[2]: energy_grid_charge_em1 = np.linspace(start=0, stop=atoms_charge_em1.get_cell()[axis][axis], num=data_charge_em1.shape[axis])
if read_labels[2]: energy_grid_hartree_em1 = np.linspace(start=0, stop=atoms_charge_em1.get_cell()[axis][axis], num=data_hartree_em1.shape[axis])
if read_labels[1]: energy_grid_charge_dft = np.linspace(start=0, stop=atoms_charge_dft.get_cell()[axis][axis], num=data_charge_dft.shape[axis])
if read_labels[1]: energy_grid_hartree_dft = np.linspace(start=0, stop=atoms_charge_dft.get_cell()[axis][axis], num=data_hartree_dft.shape[axis])
if read_labels[0]: energy_grid_charge_bulk = np.linspace(start=0, stop=atoms_charge_bulk.get_cell()[axis][axis], num=data_charge_bulk.shape[axis])
if read_labels[0]: energy_grid_hartree_bulk = np.linspace(start=0, stop=atoms_charge_bulk.get_cell()[axis][axis], num=data_hartree_bulk.shape[axis])

# Setup axis
mid_index = np.argmin(abs(energy_grid_hartree_em1-mid_pos))
mid_pos_grid = energy_grid_hartree_em1[mid_index]

# Print value of potential
if read_labels[3]: print('\nMin of Hartree potential z EM2 .cube', np.min(average_hartree_em2))
if read_labels[2]: print('Min of Hartree potential z EM1 .cube', np.min(average_hartree_em1))
if read_labels[1]: print('Min of Hartree potential z DFT .cube', np.min(average_hartree_dft))
if read_labels[0]: print('Min of Hartree potential z bulk .cube', np.min(average_hartree_bulk))
if read_labels[3]: print('\nValue(z=Min) of Hartree potential z EM2 .cube', average_hartree_em2[0])
if read_labels[2]: print('Value(z=Min) of Hartree potential z EM1 .cube', average_hartree_em1[0])
if read_labels[1]: print('Value(z=Min) of Hartree potential z DFT .cube', average_hartree_dft[0])
if read_labels[0]: print('Value(z=Min) of Hartree potential z bulk .cube', average_hartree_bulk[0])
if read_labels[3]: print('Value(z=Max) of Hartree potential z EM2 .cube', average_hartree_em2[-1])
if read_labels[2]: print('\nValue(z=Max) of Hartree potential z EM1 .cube', average_hartree_em1[-1])
if read_labels[1]: print('Value(z=Max) of Hartree potential z DFT .cube', average_hartree_dft[-1])
if read_labels[0]: print('Value(z=Max) of Hartree potential z bulk .cube', average_hartree_bulk[-1])

# Plot all
rows, cols = 2, 1
xlim = [0-1, atoms_charge_em1.get_cell()[axis][axis]+1]
fig_cube_z, ax_cube_z = plt.subplots(rows, cols, sharex='col', sharey='row', figsize=(6, 8))
if plot_fermi: ax_cube_z[0].axhline(y=fermi_dft, color='grey', linestyle='--', label='DFT Fermi energy', alpha=0.5)
if draw_mirror_em2: ax_cube_z[0].plot(energy_grid_hartree_em2[:mid_index]+mid_pos_grid, np.flip(average_hartree_em2[:mid_index]), 'm--')
if draw_mirror_em1: ax_cube_z[0].plot(energy_grid_hartree_em1[:mid_index]+mid_pos_grid, np.flip(average_hartree_em1[:mid_index]), 'm--')
if plot_labels[1]: ax_cube_z[0].plot(energy_grid_hartree_dft, average_hartree_dft, '-', color=plot_color[1], label=labels[1])
if plot_labels[3]: ax_cube_z[0].plot(energy_grid_hartree_em2, average_hartree_em2, '-', color=plot_color[3], label=labels[3])
if plot_labels[2]: ax_cube_z[0].plot(energy_grid_hartree_em1, average_hartree_em1, '-', color=plot_color[2], label=labels[2])
if plot_labels[0]: ax_cube_z[0].plot(energy_grid_hartree_bulk, average_hartree_bulk, '-', color=plot_color[0], label=labels[0])
if draw_markers: ax_cube_z[0].plot(markers, markers*0, 'o', color='orange', fillstyle='none')
ax_cube_z[0].set_xlim([xlim[0], xlim[1]])
ax_cube_z[0].legend(frameon=False)
ax_cube_z[0].set_ylabel('Hartree potential / eV')
if draw_mirror_em2: ax_cube_z[1].plot(energy_grid_charge_em2[:mid_index]+mid_pos_grid, np.flip(average_charge_em2[:mid_index]), 'm--')
if draw_mirror_em1: ax_cube_z[1].plot(energy_grid_charge_em1[:mid_index]+mid_pos_grid, np.flip(average_charge_em1[:mid_index]), 'm--')
if plot_labels[1]: ax_cube_z[1].plot(energy_grid_charge_dft, average_charge_dft, '-', color=plot_color[1], label=labels[1])
if plot_labels[3]: ax_cube_z[1].plot(energy_grid_charge_em2, average_charge_em2, '-', color=plot_color[3], label=labels[3])
if plot_labels[2]: ax_cube_z[1].plot(energy_grid_charge_em1, average_charge_em1, '-', color=plot_color[2], label=labels[2])
if plot_labels[0]: ax_cube_z[1].plot(energy_grid_charge_bulk, average_charge_bulk, '-', color=plot_color[0], label=labels[0])
if draw_markers: ax_cube_z[1].plot(markers, markers*0, 'o', color='orange', fillstyle='none')
ax_cube_z[1].set_xlim([xlim[0], xlim[1]])
ax_cube_z[1].legend(frameon=False)
ax_cube_z[1].set_xlabel(r'Position / Å')
ax_cube_z[1].set_ylabel('Charge density')
fig_cube_z.tight_layout()
fig_cube_z.savefig('{}/charge_hartree_cube_{}.png'.format(folder_save, print_label), dpi=300)
print('Finished plotting average')

# Plot Hartree and charge .cube difference
if plot_diff:
    fig_cube_both, ax_cube_both = plt.subplots(rows, cols, sharex='col', sharey='row', figsize=(6, 8))
    if plot_diff_em1_em2: ax_cube_both[0].plot(energy_grid_hartree_em1, average_hartree_em2-average_hartree_em1, '-', color=diff_color[2], label='V=1 - V=0')
    if plot_diff_dft_em2: ax_cube_both[0].plot(energy_grid_hartree_em1, average_hartree_em2-average_hartree_dft, '-', color=diff_color[1], label='V=1 - DFT')
    if plot_diff_dft_em1: ax_cube_both[0].plot(energy_grid_hartree_em1, average_hartree_em1-average_hartree_dft, '-', color=diff_color[0], label='V=0 - DFT')
    if draw_markers: ax_cube_both[0].plot(markers, markers * 0, 'o', color='orange', fillstyle='none')
    if plot_diff_legend: ax_cube_both[0].legend(frameon=False)
    ax_cube_both[0].set_xlim([xlim[0], xlim[1]])
    ax_cube_both[0].set_xlabel(r'Position / Å')
    ax_cube_both[0].set_ylabel('Hartree potential z / eV')
    if plot_diff_em1_em2: ax_cube_both[1].plot(energy_grid_charge_em1, average_charge_em2-average_charge_em1, '-', color=diff_color[2], label='V=1 - V=0')
    if plot_diff_dft_em2: ax_cube_both[1].plot(energy_grid_charge_em1, average_charge_em2-average_charge_dft, '-', color=diff_color[1], label='V=1 - DFT')
    if plot_diff_dft_em1: ax_cube_both[1].plot(energy_grid_charge_em1, average_charge_em1-average_charge_dft, '-', color=diff_color[0], label='V=0 - DFT')
    if draw_markers: ax_cube_both[1].plot(markers, markers * 0, 'o', color='orange', fillstyle='none')
    if plot_diff_legend: ax_cube_both[1].legend(frameon=False)
    ax_cube_both[1].set_xlim([xlim[0], xlim[1]])
    ax_cube_both[1].set_xlabel(r'Position / Å')
    ax_cube_both[1].set_ylabel('Charge density z')
    fig_cube_both.tight_layout()
    fig_cube_both.savefig('{}/charge_hartree_cube_diff_{}.png'.format(folder_save, print_label), dpi=300)
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
