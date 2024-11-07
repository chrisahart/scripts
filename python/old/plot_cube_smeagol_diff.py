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
# folder_2 = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-capacitor/archer/layers-1-2-3-4/iv_parralel/kpoints-4-4-20_hlb-auto_vacuum/iv_curve/V_-1.0'
# folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-capacitor/archer/layers-1-2-3-4/iv_parralel/kpoints-4-4-20_hlb-auto_vacuum/iv_curve/V_0.0'
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

folder_2 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/delete/au-q1-capacitor-kpoints-1-1-20-bulk-4layers-lattice-cu'
# folder_2 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/au-q1-au-capacitor-kpoints-1-1-20-bulk-4layers-lattice-cu-10layers-NEnergReal-64-omp-2'
# folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/au-q1-capacitor-kpoints-1-1-20-bulk-4layers-lattice-cu'
folder_1 = folder_2
# folder_2 = folder_1
file_charge_2 = '1V-ELECTRON_DENSITY-1_0.cube'
file_hartree_2 = '1V-v_hartree-1_0.cube'
# file_charge_1 = 'dft_wfn-ELECTRON_DENSITY-1_0.cube'
# file_hartree_1 = 'dft_wfn-v_hartree-1_0.cube'
file_charge_1 = '0V-ELECTRON_DENSITY-1_0.cube'
file_hartree_1 = '0V-v_hartree-1_0.cube'
plot_markers = False
draw_mirror = False
mirror_factor = -1
markers = np.array(
    [2.08400 * 0, 2.08400 * 1, 2.08400 * 2, 2.08400 * 3, 2.08400 * 4, 2.08400 * 5, 21.42200 + 2.08400 * 0,
     21.42200 + 2.08400 * 1, 21.42200 + 2.08400 * 2, 21.42200 + 2.08400 * 3, 21.42200 + 2.08400 * 4,
     21.42200 + 2.08400 * 5])
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
mid_pos = 2.08400*5+(21.42200-2.08400*5)/2
mid_index = np.argmin(abs(energy_grid_z_2-mid_pos))
mid_pos_grid = energy_grid_z_2[mid_index]
print(mid_pos)
print(mid_pos_grid)

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
