import numpy as np
import matplotlib.pyplot as plt
from general import parameters as param
from ase.io.cube import read_cube_data

""" Plotting of SMEAGOL output _TRC.agr by filename"""

# Au capacitor
# CP2K
# folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/cp2k-smeagol/sz/transmission/exp/capacitor/sergey-equal/kpoints_bulk-4-4-100_em-4-4-1_hlb-t-10.99872_scf-500'
# folder_2 = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/cp2k-smeagol/sz/transmission/exp/capacitor/sergey-equal/V-1_kpoints_bulk-4-4-100_em-4-4-1_hlb-t-10.99872_scf-500'
# folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/cp2k-smeagol/sz/transmission/exp/capacitor/sergey-equal/kpoints_bulk-31_em-1-1-1_hlb-t-11.03197_scf-500'
# folder_2 = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/cp2k-smeagol/sz/transmission/exp/capacitor/sergey-equal/V-1_kpoints_bulk-31_em-1-1-1_hlb-t-11.03197_scf-500'
# folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/cp2k-smeagol/sz/transmission/exp/capacitor/sergey-equal/auto/V-0_kpoints_bulk-2-2-100_em-2-2-1_hlb-auto'
# folder_2 = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/cp2k-smeagol/sz/transmission/exp/capacitor/sergey-equal/auto/V-1_kpoints_bulk-2-2-100_em-2-2-1_hlb-auto'
# folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/cp2k-smeagol/sz/transmission/exp/capacitor/sergey-equal/auto/V-0_kpoints_bulk-4-4-100_em-4-4-1_hlb-auto'
# folder_2 = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/cp2k-smeagol/sz/transmission/exp/capacitor/sergey-equal/auto/V-1_kpoints_bulk-4-4-100_em-4-4-1_hlb-auto'
# folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/cp2k-smeagol/sz/transmission/exp/capacitor/sergey-equal/auto/kpoints-20/V-0_kpoints_bulk-4-4-20_em-4-4-1_hlb-auto'
# folder_2 = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/cp2k-smeagol/sz/transmission/exp/capacitor/sergey-equal/auto/kpoints-20/V-1_kpoints_bulk-4-4-20_em-4-4-1_hlb-auto'
# file_charge_1 = '0V-ELECTRON_DENSITY-1_0.cube'
# file_hartree_1 = '0V-v_hartree-1_0.cube'
# file_charge_2 = '0V-ELECTRON_DENSITY-1_0.cube'
# file_hartree_2 = '0V-v_hartree-1_0.cube'
# plot_markers = True
# markers = np.array([2.08400*0, 2.08400*1, 2.08400*2, 2.08400*3, 2.08400*4, 2.08400*5, 21.42200+2.08400*0,
#                     21.42200+2.08400*1, 21.42200+2.08400*2, 21.42200+2.08400*3, 21.42200+2.08400*4, 21.42200+2.08400*5])
# use_xlim = False
# xlim_specify_left = [-0.1, 6]
# xlim_specify_right = [78, 84.1]
# xlim_specify = xlim_specify_left
#
# # SIESTA
# bulk_charge = '0.bulk-RHO_AV.dat'
# bulk_hartree = '0.bulk-VH_AV.dat'
# em_hartree = '0.transport-VH_AV.dat'
# em_charge = '0.transport-RHO_AV.dat'
# labels_em = ['EM']
# labels_bulk = ['Bulk']
# ylim = [-3, 2]
# plot_markers = True
# markers = np.array([2.08400*0, 2.08400*1, 2.08400*2, 2.08400*3, 2.08400*4, 2.08400*5, 21.42200+2.08400*0,
#                     21.42200+2.08400*1, 21.42200+2.08400*2, 21.42200+2.08400*3, 21.42200+2.08400*4, 21.42200+2.08400*5])
# x = 33.926 - 8.336
# folder = ['/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/siesta-smeagol/capacitor/bulk-4-4-100-em-4-4-1_hlb-15.2496_0-0',
#           '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/siesta-smeagol/capacitor/V-1_bulk-4-4-100-em-4-4-1_hlb-15.2496_0-0']
# folder_save = []
# for i in range(len(folder)):
#     folder_save.append('{}/output'.format(folder[i]))
# print('Saving to', folder_save)

# Au chain
# CP2K
folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-chain/cp2k-smeagol/transmission/testing/chris/HLB-auto_kpoints-4-4-20'
folder_2 = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-chain/cp2k-smeagol/transmission/testing/chris/HLB-auto_kpoints-4-4-20'
file_charge_1 = '0V-ELECTRON_DENSITY-1_0.cube'
file_hartree_1 = '0V-v_hartree-1_0.cube'
file_charge_2 = 'V-ELECTRON_DENSITY-1_0.cube'
file_hartree_2 = 'V-v_hartree-1_0.cube'
plot_markers = True
markers = np.array([2.08400*0, 2.08400*1, 2.08400*2, 2.08400*3, 2.08400*4, 2.08400*5, 21.42200+2.08400*0,
                    21.42200+2.08400*1, 21.42200+2.08400*2, 21.42200+2.08400*3, 21.42200+2.08400*4, 21.42200+2.08400*5])
use_xlim = False
xlim_specify_left = [-0.1, 6]
xlim_specify_right = [78, 84.1]
xlim_specify = xlim_specify_left

# SIESTA
bulk_charge = '0.bulk-RHO_AV.dat'
bulk_hartree = '0.bulk-VH_AV.dat'
em_hartree = '0.transport-VH_AV.dat'
em_charge = '0.transport-RHO_AV.dat'
labels_em = ['EM']
labels_bulk = ['Bulk']
ylim = [-3, 2]
plot_markers = False
markers = np.array([2.08400*0, 2.08400*1, 2.08400*2, 2.08400*3, 2.08400*4, 2.08400*5, 21.42200+2.08400*0,
                    21.42200+2.08400*1, 21.42200+2.08400*2, 21.42200+2.08400*3, 21.42200+2.08400*4, 21.42200+2.08400*5])
x = 33.926 - 8.336
folder = ['/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-chain/siesta/transmission/testing/single-points/V-0_HLB-auto_kpoints-1-1-20',
          '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-chain/siesta/transmission/testing/single-points/V-1_HLB-auto_kpoints-1-1-20']
folder_save = []
for i in range(len(folder)):
    folder_save.append('{}/output'.format(folder[i]))
print('Saving to', folder_save)


data_bulk_charge = []
data_bulk_hartree = []
data_em_hartree = []
data_em_charge = []
data_cp2k_hartree = []
data_cp2k_charge = []
for i in range(len(folder)):
    data_bulk_charge.append(np.genfromtxt('{}/{}'.format(folder[i], bulk_charge), skip_header=0, skip_footer=0))
    data_bulk_hartree.append(np.genfromtxt('{}/{}'.format(folder[i], bulk_hartree), skip_header=0, skip_footer=0))
    data_em_hartree.append(np.genfromtxt('{}/{}'.format(folder[i], em_hartree), skip_header=0, skip_footer=0))
    data_em_charge.append(np.genfromtxt('{}/{}'.format(folder[i], em_charge), skip_header=0, skip_footer=0))


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

# Plot Hartree potential and charge density difference
rows, cols = 2, 1
fig_plot_both_diff, ax_plot_both_diff = plt.subplots(rows, cols, sharex='col', sharey='row', figsize=(6, 8))
ax_plot_both_diff[0].plot(data_em_hartree[0][:, 0], data_em_hartree[1][:, 1]-data_em_hartree[0][:, 1], 'r-', label='SIESTA-SMEAGOL')
ax_plot_both_diff[0].plot(energy_grid_z_4, z_average_4-z_average_2, 'g-', label='CP2K-SMEAGOL')
if plot_markers: ax_plot_both_diff[0].plot(markers, markers*0, 'o', color='orange', fillstyle='none')
ax_plot_both_diff[0].legend(frameon=False)
ax_plot_both_diff[0].set_ylabel('Hartree potential / eV')
ax_plot_both_diff[1].plot(data_em_charge[0][:, 0], (data_em_charge[1][:, 1]-data_em_charge[0][:, 1])/np.max(abs((data_em_charge[1][:, 1]-data_em_charge[0][:, 1]),)), 'r-', label='SIESTA-SMEAGOL')
ax_plot_both_diff[1].plot(energy_grid_z_3, (z_average_3-z_average_1)/np.max(abs((z_average_3-z_average_1))), 'g-', label='CP2K-SMEAGOL')
if plot_markers: ax_plot_both_diff[1].plot(markers, markers*0, 'o', color='orange', fillstyle='none')
# ax_plot_both_diff[1].legend(frameon=False)
ax_plot_both_diff[1].set_xlabel(r'Position / Å')
ax_plot_both_diff[1].set_ylabel('Normalised charge density')
fig_plot_both_diff.tight_layout()
for i in range(len(folder)):
    fig_plot_both_diff.savefig('{}/siesta_cp2k_compare.png'.format(folder_save[i]), dpi=param.save_dpi)
fig_plot_both_diff.savefig('{}/siesta_cp2k_compare.png'.format(folder_1), dpi=300)
fig_plot_both_diff.savefig('{}/siesta_cp2k_compare.png'.format(folder_2), dpi=300)

if __name__ == "__main__":
    print('Finished.')
    plt.show()
