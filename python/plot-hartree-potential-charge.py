import numpy as np
import matplotlib.pyplot as plt
from general import parameters as param

""" Plotting of SMEAGOL output _TRC.agr by filename"""

# SIESTA Lithium 1D wire
# bulk_charge = '0.Li.leads-RHO_AV.dat'
# bulk_hartree = '0.Li.leads-VH_AV.dat'
# em_hartree = '0.Liwire-VH_AV.dat'
# em_charge = '0.Liwire-RHO_AV.dat'
# labels_em = ['EM']
# labels_bulk = ['Bulk']
# ylim = [-3, 2]
# folder = [
#     '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/iv/li/siesta/single-points/v-0_bottom-0',
#     '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/iv/li/siesta/single-points/v-0_bottom-1']

# CP2K Lithium 1D wire 3e
bulk_charge = 'bulk-RHO_z.dat'
bulk_hartree = 'bulk-VH_AV.dat'
em_hartree = '0V-VH_z.dat'
em_charge = '0V-RHO_z.dat'
labels_em = ['EM']
labels_bulk = ['Bulk']
ylim = [-3, 2]
folder = [
    '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/iv/li/cp2k/transmission/dev-chris/single-points/V-0_HLB-F_z-0-0/',
    '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/iv/li/cp2k/transmission/dev-chris/single-points/V-0_HLB-1_z-0-0']

# CP2K Lithium 1D wire 1e
# bulk_charge = 'bulk-RHO_z.dat'
# bulk_hartree = 'bulk-VH_AV.dat'
# em_hartree = '0V-VH_z.dat'
# em_charge = '0V-RHO_z.dat'
# labels_em = ['EM']
# labels_bulk = ['Bulk']
# ylim = [-3, 2]
# folder = [
#     '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/iv/li/cp2k/transmission/dev-chris/single-points/V-0_HLB-F_z-0-0-pbe-DZVP-MOLOPT-SR-GTH-q1',
#     '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/iv/li/cp2k/transmission/dev-chris/single-points/V-0_HLB-1_z-0-0-pbe-DZVP-MOLOPT-SR-GTH-q1']

# CP2K Au 1D wire
# bulk_charge = 'bulk-RHO_z.dat'
# bulk_hartree = 'bulk-VH_AV.dat'
# em_hartree = '0V-VH_z.dat'
# em_charge = '0V-RHO_z.dat'
# labels_em = ['EM']
# labels_bulk = ['Bulk']
# ylim = [-3, 2]
# folder = [
#     '/Volumes/Storage/Data/delete/Work/Postdoc/Work/calculations/transport/iv/au/cp2k/dev/single-points/V-0_HLB-0_z-0-0',
#     '/Volumes/Storage/Data/delete/Work/Postdoc/Work/calculations/transport/iv/au/cp2k/dev/single-points/V-0_HLB-1_z-0-0']

# SIESTA Au(s,d) 1D wire
# bulk_charge = '0.Au_wire_leads-RHO_AV.dat'
# bulk_hartree = '0.Au_wire_leads-VH_AV.dat'
# em_hartree = '0.Au-VH_AV.dat'
# em_charge = '0.Au-RHO_AV.dat'
# labels_em = ['EM']
# labels_bulk = ['Bulk']
# ylim = [-3, 2]
# folder = [
#     '/Volumes/Storage/Data/delete/Work/Postdoc/Work/calculations/transport/iv/au/siesta/V-0_HLB-0_z-0-0',
#     '/Volumes/Storage/Data/delete/Work/Postdoc/Work/calculations/transport/iv/au/siesta/V-0_HLB-1_z-0-0']

# SIESTA Au(s,p,d) 1D wire
# bulk_charge = '0.Au_wire_leads-RHO_AV.dat'
# bulk_hartree = '0.Au_wire_leads-VH_AV.dat'
# em_hartree = '0.Au-VH_AV.dat'
# em_charge = '0.Au-RHO_AV.dat'
# labels_em = ['EM']
# labels_bulk = ['Bulk']
# ylim = [-3, 2]
# folder = [
#     '/Volumes/Storage/Data/delete/Work/Postdoc/Work/calculations/transport/iv/au/siesta/V-0_HLB-0_z-0-0-sdp',
#     '/Volumes/Storage/Data/delete/Work/Postdoc/Work/calculations/transport/iv/au/siesta/V-0_HLB-1_z-0-0-sdp']

# CP2K Lithium 1D wire 1e
# bulk_charge = 'bulk-RHO_z.dat'
# bulk_hartree = 'bulk-VH_AV.dat'
# em_hartree = '0V-VH_z_ORIG-0001.dat'
# em_charge = 'RHO_AO-0001.txt'
# labels_em = ['EM']
# labels_bulk = ['Bulk']
# ylim = [-3, 2]
# folder = [
#     '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/au-capacitor-100/struct_5-4/cp2k-smeagol/single-points/old/V-0_HLB-F_z-0-0_hirshfeld-gaussian_ordered']

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

# Hartree potential
# fig_plot_1, ax_plot_1 = plt.subplots()
# ax_plot_1.plot(data_em_hartree[0][:, 0], data_em_hartree[0][:, 1], 'r-', label='HLB 0')
# ax_plot_1.plot(data_em_hartree[1][:, 0], data_em_hartree[1][:, 1], 'g-', label='HLB 1 eV')
# ax_plot_1.legend(frameon=False)
# ax_plot_1.set_xlabel(r'Position / Å')
# ax_plot_1.set_ylabel('Hartree potential / eV')
# fig_plot_1.tight_layout()
# for i in range(len(folder)):
#     fig_plot_1.savefig('{}/hartree_potential.png'.format(folder[i]), dpi=param.save_dpi)

# Plot charge density
# fig_plot_2, ax_plot_2 = plt.subplots()
# ax_plot_2.plot(data_em_charge[0][:, 0], data_em_charge[0][:, 1], 'r-', label='HLB 0')
# ax_plot_2.plot(data_em_charge[1][:, 0], data_em_charge[1][:, 1], 'g-', label='HLB 1 eV')
# print(np.sum(np.abs(data_em_charge[0][:, 1])))
# ax_plot_2.legend(frameon=False)
# ax_plot_2.set_xlabel(r'Position / Å')
# ax_plot_2.set_ylabel('Charge density')
# fig_plot_2.tight_layout()
# for i in range(len(folder)):
#     fig_plot_2.savefig('{}/charge_density.png'.format(folder[i]), dpi=param.save_dpi)

# Plot charge density difference
# fig_plot_3, ax_plot_3 = plt.subplots()
# ax_plot_3.plot(data_em_charge[0][:, 0], data_em_charge[0][:, 1]-data_em_charge[1][:, 1], 'k-')
# ax_plot_3.set_xlabel(r'Position / Å')
# ax_plot_3.set_ylabel('Charge density difference (0 - 1)')
# fig_plot_3.tight_layout()
# for i in range(len(folder)):
#     fig_plot_3.savefig('{}/charge_density_difference.png'.format(folder[i]), dpi=param.save_dpi)

# Plot Hartree potential and charge density
# rows, cols = 2, 1
# fig_plot_both, ax_plot_both = plt.subplots(rows, cols, sharex='col', sharey='row', figsize=(6, 8))
# ax_plot_both[0].plot(data_em_hartree[0][:, 0], data_em_hartree[0][:, 1], 'g-', label='EM HLB 0=F')
# ax_plot_both[0].plot(data_em_hartree[1][:, 0], data_em_hartree[1][:, 1], 'r-', label='EM HLB 1 eV')
# ax_plot_both[0].plot(data_bulk_hartree[0][:, 0], data_bulk_hartree[0][:, 1], 'k-', label='Bulk')
# ax_plot_both[0].legend(frameon=False)
# ax_plot_both[0].set_ylabel('Hartree potential / eV')
# ax_plot_both[1].plot(data_em_charge[0][:, 0], data_em_charge[0][:, 1], 'g-', label='EM HLB 0=F')
# ax_plot_both[1].plot(data_em_charge[1][:, 0], data_em_charge[1][:, 1], 'r-', label='EM HLB 1 eV')
# ax_plot_both[1].plot(data_bulk_charge[0][:, 0], data_bulk_charge[0][:, 1], 'k-', label='Bulk')
# ax_plot_both[1].legend(frameon=False)
# ax_plot_both[1].set_xlabel(r'Position / Å')
# ax_plot_both[1].set_ylabel('Charge density')
# fig_plot_both.tight_layout()
# for i in range(len(folder)):
#     fig_plot_both.savefig('{}/hartree_potential_charge.png'.format(folder[i]), dpi=param.save_dpi)

# Plot Hartree potential and charge density
rows, cols = 2, 1
fig_plot_both, ax_plot_both = plt.subplots(rows, cols, sharex='col', sharey='row', figsize=(6, 8))
ax_plot_both[0].plot(data_em_hartree[0][:, 0], data_em_hartree[0][:, 1], 'g-', label='EM')
ax_plot_both[0].plot(data_bulk_hartree[0][:, 0], data_bulk_hartree[0][:, 1], 'k-', label='Bulk')
ax_plot_both[0].legend(frameon=False)
ax_plot_both[0].set_ylabel('Hartree potential / eV')
ax_plot_both[1].plot(data_em_charge[0][:, 0], data_em_charge[0][:, 1], 'g-', label='EM')
ax_plot_both[1].plot(data_bulk_charge[0][:, 0], data_bulk_charge[0][:, 1], 'k-', label='Bulk')
ax_plot_both[1].legend(frameon=False)
ax_plot_both[1].set_xlabel(r'Position / Å')
ax_plot_both[1].set_ylabel('Charge density')
fig_plot_both.tight_layout()
for i in range(len(folder)):
    fig_plot_both.savefig('{}/hartree_potential_charge.png'.format(folder[i]), dpi=param.save_dpi)

if __name__ == "__main__":
    print('Finished.')
    plt.show()
