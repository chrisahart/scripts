import numpy as np
from matplotlib import pyplot as plt
from ase.io.cube import read_cube_data
from general import parameters as param

# CP2K Li chain 3e
# folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/iv/li/cp2k/transmission/dev-chris/single-points/V-0_HLB-F_z-0-0'
# file_charge_em = '0V-ELECTRON_DENSITY-1_0.cube'
# file_hartree_em = '0V-v_hartree-1_0.cube'
# file_charge_dft = 'dft_wfn-ELECTRON_DENSITY-1_0.cube'
# file_hartree_dft = 'dft_wfn-v_hartree-1_0.cube'
# file_charge_bulk = 'bulk-ELECTRON_DENSITY-1_0.cube'
# file_hartree_bulk = 'bulk-v_hartree-1_0.cube'

# CP2K Li chain 1e
# folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/iv/li/cp2k/transmission/dev-chris/single-points/V-0_HLB-F_z-0-0-pbe-DZVP-MOLOPT-SR-GTH-q1'
# file_charge_em = '0V-ELECTRON_DENSITY-1_0.cube'
# file_hartree_em = '0V-v_hartree-1_0.cube'
# file_charge_dft = 'dft_wfn-ELECTRON_DENSITY-1_0.cube'
# file_hartree_dft = 'dft_wfn-v_hartree-1_0.cube'
# file_charge_bulk = 'bulk-ELECTRON_DENSITY-1_0.cube'
# file_hartree_bulk = 'bulk-v_hartree-1_0.cube'

# CP2K Au capacitor 1
folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/au-capacitor/cp2k-smeagol/single-points/V-0_HLB-F_z-0-0_hirshfeld-gaussian_nodummy_hartree'
file_charge_em = '0V-ELECTRON_DENSITY-1_0.cube'
file_hartree_em = '0V-v_hartree-1_0.cube'
file_charge_dft = 'dft_wfn-ELECTRON_DENSITY-1_0.cube'
file_hartree_dft = 'dft_wfn-v_hartree-1_0.cube'
file_charge_bulk = 'bulk-ELECTRON_DENSITY-1_0.cube'
file_hartree_bulk = 'bulk-v_hartree-1_0.cube'
fermi_dft = 0.20301667225796 * param.hartree_to_ev

# CP2K Au capacitor 2
# folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/au-capacitor/cp2k-smeagol/single-points/V-0_HLB-F_z-0-0_hirshfeld-gaussian_ordered_hartree'
# file_charge_em = '0V-ELECTRON_DENSITY-1_0.cube'
# file_hartree_em = '0V-v_hartree-1_0.cube'
# file_charge_dft = 'dft_wfn-ELECTRON_DENSITY-1_0.cube'
# file_hartree_dft = 'dft_wfn-v_hartree-1_0.cube'
# file_charge_bulk = 'bulk-ELECTRON_DENSITY-1_0.cube'
# file_hartree_bulk = 'bulk-v_hartree-1_0.cube'
# fermi_dft = 0.20359259837197 * param.hartree_to_ev

# Read .cube using ASE
data_1, atoms_1 = read_cube_data('{}/{}'.format(folder, file_charge_em))
data_2, atoms_2 = read_cube_data('{}/{}'.format(folder, file_hartree_em))
data_3, atoms_3 = read_cube_data('{}/{}'.format(folder, file_charge_dft))
data_4, atoms_4 = read_cube_data('{}/{}'.format(folder, file_hartree_dft))
data_5, atoms_5 = read_cube_data('{}/{}'.format(folder, file_charge_bulk))
data_6, atoms_6 = read_cube_data('{}/{}'.format(folder, file_hartree_bulk))

# Plot charge and hartree .cube x average
axis = 0
x_average_1 = np.zeros(data_1.shape[axis])
x_average_2 = np.zeros(data_2.shape[axis])
x_average_3 = np.zeros(data_3.shape[axis])
x_average_4 = np.zeros(data_4.shape[axis])
x_average_5 = np.zeros(data_5.shape[axis])
x_average_6 = np.zeros(data_6.shape[axis])
for i in range(data_1.shape[axis]):
    x_average_1[i] = np.mean(data_1[:, i, :])
    x_average_2[i] = np.mean(data_2[:, i, :] * param.hartree_to_ev)
    x_average_3[i] = np.mean(data_3[:, i, :])
    x_average_4[i] = np.mean(data_4[:, i, :] * param.hartree_to_ev)
for i in range(data_5.shape[axis]):
    x_average_5[i] = np.mean(data_5[:, i, :])
    x_average_6[i] = np.mean(data_6[:, i, :] * param.hartree_to_ev)
energx_grid_x_1 = np.linspace(start=0, stop=atoms_1.get_cell()[axis][axis], num=data_1.shape[axis])
energx_grid_x_2 = np.linspace(start=0, stop=atoms_5.get_cell()[axis][axis], num=data_5.shape[axis])
xlim = [0 - 1, atoms_1.get_cell()[axis][axis] + 1]
print('\nMax of Hartree potential x EM .cube', np.max(x_average_2))
print('Max of Hartree potential x DFT .cube', np.max(x_average_4))
print('Max of Hartree potential x bulk .cube', np.max(x_average_6))

rows, cols = 2, 1
figcube_both, ax_cube_x = plt.subplots(rows, cols, sharex='col', sharey='row', figsize=(6, 8))
ax_cube_x[0].plot(energx_grid_x_1, x_average_4, 'r-', label='DFT .cube')
ax_cube_x[0].plot(energx_grid_x_1, x_average_2, 'g-', label='EM .cube')
ax_cube_x[0].plot(energx_grid_x_2, x_average_6, 'k-', label='Bulk .cube')
ax_cube_x[0].set_xlim([xlim[0], xlim[1]])
ax_cube_x[0].legend(frameon=False)
ax_cube_x[0].set_xlabel(r'Position / Å')
ax_cube_x[0].set_ylabel('Hartree potential x / eV')
ax_cube_x[1].plot(energx_grid_x_1, x_average_3, 'r-', label='DFT .cube')
ax_cube_x[1].plot(energx_grid_x_1, x_average_1, 'g-', label='EM .cube')
ax_cube_x[1].plot(energx_grid_x_2, x_average_5, 'k-', label='Bulk .cube')
ax_cube_x[1].set_xlim([xlim[0], xlim[1]])
ax_cube_x[1].legend(frameon=False)
ax_cube_x[1].set_xlabel(r'Position / Å')
ax_cube_x[1].set_ylabel('Charge density x')
figcube_both.tight_layout()
figcube_both.savefig('{}/charge_hartree_cube_x.png'.format(folder), dpi=300)

# Plot charge and hartree .cube y average
axis = 1
y_average_1 = np.zeros(data_1.shape[axis])
y_average_2 = np.zeros(data_2.shape[axis])
y_average_3 = np.zeros(data_3.shape[axis])
y_average_4 = np.zeros(data_4.shape[axis])
y_average_5 = np.zeros(data_5.shape[axis])
y_average_6 = np.zeros(data_6.shape[axis])
for i in range(data_1.shape[axis]):
    y_average_1[i] = np.mean(data_1[:, i, :])
    y_average_2[i] = np.mean(data_2[:, i, :] * param.hartree_to_ev)
    y_average_3[i] = np.mean(data_3[:, i, :])
    y_average_4[i] = np.mean(data_4[:, i, :] * param.hartree_to_ev)
for i in range(data_5.shape[axis]):
    y_average_5[i] = np.mean(data_5[:, i, :])
    y_average_6[i] = np.mean(data_6[:, i, :] * param.hartree_to_ev)
energy_grid_y_1 = np.linspace(start=0, stop=atoms_1.get_cell()[axis][axis], num=data_1.shape[axis])
energy_grid_y_2 = np.linspace(start=0, stop=atoms_5.get_cell()[axis][axis], num=data_5.shape[axis])

xlim = [0 - 1, atoms_1.get_cell()[axis][axis] + 1]
print('\nMax of Hartree potential y EM .cube', np.max(y_average_2))
print('Max of Hartree potential y DFT .cube', np.max(y_average_4))
print('Max of Hartree potential y bulk .cube', np.max(y_average_6))
rows, cols = 2, 1
figcube_both, ax_cube_y = plt.subplots(rows, cols, sharex='col', sharey='row', figsize=(6, 8))
ax_cube_y[0].plot(energy_grid_y_1, y_average_4, 'r-', label='DFT .cube')
ax_cube_y[0].plot(energy_grid_y_1, y_average_2, 'g-', label='EM .cube')
ax_cube_y[0].plot(energy_grid_y_2, y_average_6, 'k-', label='Bulk .cube')
ax_cube_y[0].set_xlim([xlim[0], xlim[1]])
ax_cube_y[0].legend(frameon=False)
ax_cube_y[0].set_xlabel(r'Position / Å')
ax_cube_y[0].set_ylabel('Hartree potential y / eV')
ax_cube_y[1].plot(energy_grid_y_1, y_average_3, 'r-', label='DFT .cube')
ax_cube_y[1].plot(energy_grid_y_1, y_average_1, 'g-', label='EM .cube')
ax_cube_y[1].plot(energy_grid_y_2, y_average_5, 'k-', label='Bulk .cube')
ax_cube_y[1].set_xlim([xlim[0], xlim[1]])
ax_cube_y[1].legend(frameon=False)
ax_cube_y[1].set_xlabel(r'Position / Å')
ax_cube_y[1].set_ylabel('Charge density y')
figcube_both.tight_layout()
figcube_both.savefig('{}/charge_hartree_cube_y.png'.format(folder), dpi=300)

# Plot charge and hartree .cube z average
axis = 2
z_average_1 = np.zeros(data_1.shape[axis])
z_average_2 = np.zeros(data_2.shape[axis])
z_average_3 = np.zeros(data_3.shape[axis])
z_average_4 = np.zeros(data_4.shape[axis])
z_average_5 = np.zeros(data_5.shape[axis])
z_average_6 = np.zeros(data_6.shape[axis])
for i in range(data_1.shape[axis]):
    z_average_1[i] = np.mean(data_1[:, :, i])
    z_average_2[i] = np.mean(data_2[:, :, i] * param.hartree_to_ev)
    z_average_3[i] = np.mean(data_3[:, :, i])
    z_average_4[i] = np.mean(data_4[:, :, i] * param.hartree_to_ev)
for i in range(data_5.shape[axis]):
    z_average_5[i] = np.mean(data_5[:, :, i])
    z_average_6[i] = np.mean(data_6[:, :, i] * param.hartree_to_ev)

energy_grid_z_1 = np.linspace(start=0, stop=atoms_1.get_cell()[axis][axis], num=data_1.shape[axis])
energy_grid_z_2 = np.linspace(start=0, stop=atoms_5.get_cell()[axis][axis], num=data_5.shape[axis])
print('\nMax of Hartree potential z EM .cube', np.max(z_average_2))
print('Max of Hartree potential z DFT .cube', np.max(z_average_4))
print('Max of Hartree potential z bulk .cube', np.max(z_average_6))
xlim = [0-1, atoms_1.get_cell()[axis][axis]+1]
rows, cols = 2, 1
figcube_both, ax_cube_z = plt.subplots(rows, cols, sharex='col', sharey='row', figsize=(6, 8))
ax_cube_z[0].axhline(y=fermi_dft, color='grey', linestyle='--', label='DFT Fermi energy', alpha=0.5)
ax_cube_z[0].plot(energy_grid_z_1, z_average_4, 'r-', label='DFT .cube')
ax_cube_z[0].plot(energy_grid_z_1, z_average_2, 'g-', label='EM .cube')
ax_cube_z[0].plot(energy_grid_z_2, z_average_6, 'k-', label='Bulk .cube')
ax_cube_z[0].set_xlim([xlim[0], xlim[1]])
ax_cube_z[0].legend(frameon=False)
ax_cube_z[0].set_xlabel(r'Position / Å')
ax_cube_z[0].set_ylabel('Hartree potential z / eV')
ax_cube_z[1].plot(energy_grid_z_1, z_average_3, 'r-', label='DFT .cube')
ax_cube_z[1].plot(energy_grid_z_1, z_average_1, 'g-', label='EM .cube')
ax_cube_z[1].plot(energy_grid_z_2, z_average_5, 'k-', label='Bulk .cube')
ax_cube_z[1].set_xlim([xlim[0], xlim[1]])
ax_cube_z[1].legend(frameon=False)
ax_cube_z[1].set_xlabel(r'Position / Å')
ax_cube_z[1].set_ylabel('Charge density z')
figcube_both.tight_layout()
figcube_both.savefig('{}/charge_hartree_cube_z.png'.format(folder), dpi=300)

# Plot DFT .cube Hartree
fig_dft_hartree, ax_dft_hartree = plt.subplots()
ax_dft_hartree.axhline(y=fermi_dft, color='grey', linestyle='--', label='DFT Fermi energy', alpha=0.5)
ax_dft_hartree.plot(energy_grid_z_1, z_average_4, 'k-')
ax_dft_hartree.set_xlim([xlim[0], xlim[1]])
ax_dft_hartree.set_xlabel(r'Position / Å')
ax_dft_hartree.set_ylabel('Hartree potential z / eV')
fig_dft_hartree.tight_layout()
fig_dft_hartree.savefig('{}/dft_hartree.png'.format(folder), dpi=300)

wf_dft = np.max(z_average_4) - fermi_dft
print('Work function DFT', wf_dft)

if __name__ == "__main__":
    print('Finished.')
    plt.show()
