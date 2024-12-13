import numpy as np
import matplotlib.pyplot as plt
from general import parameters as param

""" Plotting of SMEAGOL output _TRC.agr by filename"""

# Au capacitor bias dependence (bound states?)
charge_1 = ['0.bulk-RHO_AV.dat'] * 3
hartree_1 = ['0.bulk-VH_AV.dat'] * 3
charge_2 = ['0.0V-RHO_AV.dat'] * 3
hartree_2 = ['0.0V-VH_AV.dat'] * 3
charge_3 = ['0.1V-RHO_AV.dat', '0.4V-RHO_AV.dat', '0.4V-RHO_AV.dat']
hartree_3 = ['0.1V-VH_AV.dat', '0.4V-VH_AV.dat', '0.4V-VH_AV.dat']
labels = ['V=0', 'V=1', 'Bulk']
ylim = [-3, 2]
plot_dft = True
plot_markers = True
plot_diff = False
plot_diff_multiple = True
plot_mirror = True
diff_label = ['V=1 - V=0', 'V=4 - V=0', 'V=4 - V=0 weighted double contour']
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
folder = ['/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/au-capacitor/siesta3-smeagol/kpoints-2-2-20_V-1',
          '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/au-capacitor/siesta3-smeagol/kpoints-4-4-20_V-4']
# folder = ['/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/au-capacitor/siesta1-smeagol/kpoints-2-2-20_V-1_WeightRho-0.5',
#           '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/au-capacitor/siesta1-smeagol/kpoints-2-2-20_V-4_WeightRho-0.5',
#           '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/au-capacitor/siesta1-smeagol/kpoints-2-2-20_V-4_double-contour']

folder1 = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/supercell-1-1-5-bulk-6-cu-1.86/junction/bias/energy/kpoints-2-2-V-0',
           '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/supercell-1-1-5-bulk-6-cu-1.86/junction/bias/energy/kpoints-2-2-V-0-rs-V-1']

folder_save = []
for i in range(len(folder)):
    folder_save.append('{}/output'.format(folder[i]))
print('Saving to', folder_save)

data_charge_1 = []
data_hartree_1 = []
data_hartree_3 = []
data_charge_2 = []
data_hartree_2 = []
data_charge_3 = []
data_cp2k_hartree = []
data_cp2k_charge = []
for i in range(len(folder)):
    data_charge_1.append(np.genfromtxt('{}/{}'.format(folder[i], charge_1[i]), skip_header=0, skip_footer=0))
    data_hartree_1.append(np.genfromtxt('{}/{}'.format(folder[i], hartree_1[i]), skip_header=0, skip_footer=0))
    if plot_dft: data_hartree_2.append(np.genfromtxt('{}/{}'.format(folder[i], hartree_2[i]), skip_header=0, skip_footer=0))
    if plot_dft: data_charge_2.append(np.genfromtxt('{}/{}'.format(folder[i], charge_2[i]), skip_header=0, skip_footer=0))
    data_hartree_3.append(np.genfromtxt('{}/{}'.format(folder[i], hartree_3[i]), skip_header=0, skip_footer=0))
    data_charge_3.append(np.genfromtxt('{}/{}'.format(folder[i], charge_3[i]), skip_header=0, skip_footer=0))

# Setup axis
mid_index = np.argmin(abs(data_hartree_3[0][:, 0]-mid_pos))
mid_pos_grid = data_hartree_3[0][mid_index, 0]

# Hartree potential
# fig_plot_1, ax_plot_1 = plt.subplots()
# ax_plot_1.plot(data_hartree_3[0][:, 0], data_hartree_3[0][:, 1], 'r-', label='HLB 0')
# ax_plot_1.plot(data_hartree_3[1][:, 0], data_hartree_3[1][:, 1], 'g-', label='HLB 1 eV')
# ax_plot_1.legend(frameon=False)
# ax_plot_1.set_xlabel(r'Position / Å')
# ax_plot_1.set_ylabel('Hartree potential / eV')
# fig_plot_1.tight_layout()
# for i in range(len(folder)):
#     fig_plot_1.savefig('{}/hartree_potential.png'.format(folder_save[i]), dpi=param.save_dpi)

# Plot charge density
# fig_plot_2, ax_plot_2 = plt.subplots()
# ax_plot_2.plot(data_charge_3[0][:, 0], data_charge_3[0][:, 1], 'r-', label='HLB 0')
# ax_plot_2.plot(data_charge_3[1][:, 0], data_charge_3[1][:, 1], 'g-', label='HLB 1 eV')
# print(np.sum(np.abs(data_charge_3[0][:, 1])))
# ax_plot_2.legend(frameon=False)
# ax_plot_2.set_xlabel(r'Position / Å')
# ax_plot_2.set_ylabel('Charge density')
# fig_plot_2.tight_layout()
# for i in range(len(folder)):
#     fig_plot_2.savefig('{}/charge_density.png'.format(folder_save[i]), dpi=param.save_dpi)

# Plot charge density difference
# fig_plot_3, ax_plot_3 = plt.subplots()
# ax_plot_3.plot(data_charge_3[0][:, 0], data_charge_3[0][:, 1]-data_charge_3[1][:, 1], 'k-')
# ax_plot_3.set_xlabel(r'Position / Å')
# ax_plot_3.set_ylabel('Charge density difference (0 - 1)')
# fig_plot_3.tight_layout()
# for i in range(len(folder)):
#     fig_plot_3.savefig('{}/charge_density_difference.png'.format(folder_save[i]), dpi=param.save_dpi)

# Plot Hartree potential and charge density
# rows, cols = 2, 1
# fig_plot_both, ax_plot_both = plt.subplots(rows, cols, sharex='col', sharey='row', figsize=(6, 8))
# ax_plot_both[0].plot(data_hartree_3[0][:, 0], data_hartree_3[0][:, 1], 'g-', label='EM HLB 0=F')
# ax_plot_both[0].plot(data_hartree_3[1][:, 0], data_hartree_3[1][:, 1], 'r-', label='EM HLB 1 eV')
# ax_plot_both[0].plot(data_hartree_1[0][:, 0], data_hartree_1[0][:, 1], 'k-', label='Bulk')
# ax_plot_both[0].legend(frameon=False)
# ax_plot_both[0].set_ylabel('Hartree potential / eV')
# ax_plot_both[1].plot(data_charge_3[0][:, 0], data_charge_3[0][:, 1], 'g-', label='EM HLB 0=F')
# ax_plot_both[1].plot(data_charge_3[1][:, 0], data_charge_3[1][:, 1], 'r-', label='EM HLB 1 eV')
# ax_plot_both[1].plot(data_charge_1[0][:, 0], data_charge_1[0][:, 1], 'k-', label='Bulk')
# ax_plot_both[1].legend(frameon=False)
# ax_plot_both[1].set_xlabel(r'Position / Å')
# ax_plot_both[1].set_ylabel('Charge density')
# fig_plot_both.tight_layout()
# for i in range(len(folder)):
#     fig_plot_both.savefig('{}/hartree_potential_charge.png'.format(folder_save[i]), dpi=param.save_dpi)

# Plot Hartree potential and charge density
rows, cols = 2, 1
fig_plot_both, ax_plot_both = plt.subplots(rows, cols, sharex='col', sharey='row', figsize=(6, 8))
if plot_dft: ax_plot_both[0].plot(data_hartree_2[0][:, 0], data_hartree_2[0][:, 1], 'r-', label=labels[0])
ax_plot_both[0].plot(data_hartree_3[0][:, 0], data_hartree_3[0][:, 1], 'g-', label=labels[1])
ax_plot_both[0].plot(data_hartree_1[0][:, 0], data_hartree_1[0][:, 1], 'k-', label=labels[2])
# ax_plot_both[0].plot(x+data_hartree_1[0][:, 0], data_hartree_1[0][:, 1], 'k-')
if plot_markers: ax_plot_both[0].plot(markers, markers*0, 'o', color='orange', fillstyle='none')
ax_plot_both[0].legend(frameon=False)
ax_plot_both[0].set_ylabel('Hartree potential / eV')
if plot_dft: ax_plot_both[1].plot(data_charge_2[0][:, 0], data_charge_2[0][:, 1], 'r-', label=labels[0])
ax_plot_both[1].plot(data_charge_3[0][:, 0], data_charge_3[0][:, 1], 'g-', label=labels[1])
ax_plot_both[1].plot(data_charge_1[0][:, 0], data_charge_1[0][:, 1], 'k-', label=labels[2])
# ax_plot_both[1].plot(x+data_charge_1[0][:, 0], data_charge_1[0][:, 1], 'k-')
if plot_markers: ax_plot_both[1].plot(markers, markers*0, 'o', color='orange', fillstyle='none')
ax_plot_both[1].legend(frameon=False)
ax_plot_both[1].set_xlabel(r'Position / Å')
ax_plot_both[1].set_ylabel('Charge density')
fig_plot_both.tight_layout()
for i in range(len(folder)):
    fig_plot_both.savefig('{}/hartree_potential_charge.png'.format(folder_save[i]), dpi=param.save_dpi)

# Plot Hartree potential and charge density difference
if plot_diff:
    rows, cols = 2, 1
    fig_plot_both_diff, ax_plot_both_diff = plt.subplots(rows, cols, sharex='col', sharey='row', figsize=(6, 8))
    # ax_plot_both_diff[0].plot(data_hartree_3[0][:, 0], data_hartree_3[1][:, 1]-data_hartree_3[0][:, 1], 'k-', label='EM')
    ax_plot_both_diff[0].plot(data_hartree_3[0][:, 0], data_hartree_3[0][:, 1]-data_hartree_2[0][:, 1], 'k-', label='EM')
    if plot_markers: ax_plot_both_diff[0].plot(markers, markers*0, 'o', color='orange', fillstyle='none')
    # ax_plot_both_diff[0].legend(frameon=False)
    ax_plot_both_diff[0].set_ylabel('Hartree potential / eV')
    # ax_plot_both_diff[1].plot(data_charge_3[0][:, 0], data_charge_3[1][:, 1]-data_charge_3[0][:, 1], 'k-', label='EM')
    ax_plot_both_diff[1].plot(data_charge_3[0][:, 0], data_charge_3[0][:, 1]-data_charge_2[0][:, 1], 'k-', label='EM')
    if plot_markers: ax_plot_both_diff[1].plot(markers, markers*0, 'o', color='orange', fillstyle='none')
    ax_plot_both_diff[1].set_xlabel(r'Position / Å')
    ax_plot_both_diff[1].set_ylabel('Charge density')
    fig_plot_both_diff.tight_layout()
    for i in range(len(folder)):
        fig_plot_both_diff.savefig('{}/hartree_potential_charge_diff.png'.format(folder_save[i]), dpi=param.save_dpi)

# Plot Hartree potential and charge density difference
rows, cols = 2, 1
fig_plot_both_diff_multiple, ax_plot_both_diff_multiple = plt.subplots(rows, cols, sharex='col', sharey='row', figsize=(6, 8))
if plot_diff_multiple:
    for i in range(len(folder)):
        ax_plot_both_diff_multiple[0].plot(data_hartree_3[i][:, 0], data_hartree_3[i][:, 1]-data_hartree_2[i][:, 1], '-', label=diff_label[i], color=diff_color[i])
        ax_plot_both_diff_multiple[1].plot(data_charge_3[i][:, 0], data_charge_3[i][:, 1]-data_charge_2[i][:, 1], '-', label=diff_label[i], color=diff_color[i])
        if plot_mirror: ax_plot_both_diff_multiple[0].plot(data_hartree_3[i][:mid_index,0] + mid_pos_grid,-np.flip(data_hartree_3[i][:mid_index, 1]-data_hartree_2[i][:mid_index, 1]), '--', color=diff_color[i])
        if plot_mirror: ax_plot_both_diff_multiple[1].plot(data_charge_3[i][:mid_index,0] + mid_pos_grid,-np.flip(data_charge_3[i][:mid_index, 1]-data_charge_2[i][:mid_index, 1]), '--', color=diff_color[i])
if plot_markers: ax_plot_both_diff_multiple[0].plot(markers, markers * 0, 'o', color='orange', fillstyle='none')
if plot_markers: ax_plot_both_diff_multiple[1].plot(markers, markers*0, 'o', color='orange', fillstyle='none')
ax_plot_both_diff_multiple[0].legend(frameon=False)
ax_plot_both_diff_multiple[1].set_xlabel(r'Position / Å')
ax_plot_both_diff_multiple[1].set_ylabel('Charge density')
ax_plot_both_diff_multiple[0].set_ylabel('Hartree potential / eV')
fig_plot_both_diff_multiple.tight_layout()
for i in range(len(folder)):
    fig_plot_both_diff_multiple.savefig('{}/hartree_potential_charge_diff_multiple.png'.format(folder_save[i]), dpi=param.save_dpi)
        
print('Value(z=Min) of Hartree potential z EM .DAT', data_hartree_3[0][0, 1])
print('Value(z=Min) of Hartree potential z bulk .DAT', data_hartree_1[0][0, 1])

print('Value(z=Max) of Hartree potential z EM .DAT', data_hartree_3[0][-1, 1])
print('Value(z=Max) of Hartree potential z bulk .DAT', data_hartree_1[0][-1, 1])

print('Average of z=min, z=max of Hartree potential z EM .DAT', np.average([data_hartree_3[0][0, 1], data_hartree_3[0][-1, 1]]))
print('Average of z=min, z=max of Hartree potential z bulk .DAT', np.average([data_hartree_1[0][0, 1], data_hartree_1[0][-1, 1]]))

if __name__ == "__main__":
    print('Finished.')
    plt.show()
