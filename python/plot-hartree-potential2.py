import numpy as np
import matplotlib.pyplot as plt
from general import parameters as param
""" Plotting of SMEAGOL output _TRC.agr by filename"""

# Ivan Au 4s wire
bulk = 'bulk-VH_AV.dat'
em = '0V-VH_z.dat'
labels_em = ['EM V=0 Δ=0', 'EM V=0 Δ=-1.8', 'EM V=2 Δ=0', 'EM V=2 Δ=-1.8', 'EM V=0 Δ=-2.8', 'EM V=0 Δ=-0.8']
labels_bulk = ['Bulk']
ylim = [-3, 2]
folder = [
    '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/iv/li/cp2k-smeagol/lda-q3_atoms-28/single-points/V-0_HLB-F_z-0-0']

# Ivan Au 4s wire vs translated
# bulk = '0.Au_wire_leads-VH_AV.dat'
# em = '0.Au-VH_AV.dat'
# labels_em = ['EM V=0 Au=0 Δ=0', 'EM V=0 Au=0 Δ=-1.8, Z=0', 'EM V=0 Au=1 Δ=0', 'EM V=0 Au=1 Δ=-0.28, Z=0', 'EM V=0 Au=1 Δ=-1.8, Z=0', 'EM V=0 Au=1 Δ=-1.8 Z=1']
# labels_bulk = ['Bulk Au=0', 'Bulk Au=0', 'Bulk Au=1', 'Bulk Au=1', 'Bulk Au=1', 'Bulk Au=1']
# # ylim = [-4, 0]
# ylim = [-3, 2]
# folder = [
#     '/Volumes/Storage/Data/Work/Postdoc/Work/calculations/testing/smeagol/scarf/Smeagol_Tutorial/Smeagol_Tutorial_Files/DFT_NEGF_Transport/Day_1/Inputs/example1/HartreeLeadsBottom/single-points/v-0_bottom-0',
#     '/Volumes/Storage/Data/Work/Postdoc/Work/calculations/testing/smeagol/scarf/Smeagol_Tutorial/Smeagol_Tutorial_Files/DFT_NEGF_Transport/Day_1/Inputs/example1/HartreeLeadsBottom/single-points/v-0_bottom-leads',
#     '/Volumes/Storage/Data/Work/Postdoc/Work/calculations/testing/smeagol/scarf/Smeagol_Tutorial/Smeagol_Tutorial_Files/DFT_NEGF_Transport/Day_1/Inputs/example1/HartreeLeadsBottom/translated/single-points/v-0_bottom-0',
#     '/Volumes/Storage/Data/Work/Postdoc/Work/calculations/testing/smeagol/scarf/Smeagol_Tutorial/Smeagol_Tutorial_Files/DFT_NEGF_Transport/Day_1/Inputs/example1/HartreeLeadsBottom/translated/single-points/v-0_bottom-0.28-shiftlr',
#     '/Volumes/Storage/Data/Work/Postdoc/Work/calculations/testing/smeagol/scarf/Smeagol_Tutorial/Smeagol_Tutorial_Files/DFT_NEGF_Transport/Day_1/Inputs/example1/HartreeLeadsBottom/translated/single-points/v-0_bottom-leads-minimum-shiftlr',
#     '/Volumes/Storage/Data/Work/Postdoc/Work/calculations/testing/smeagol/scarf/Smeagol_Tutorial/Smeagol_Tutorial_Files/DFT_NEGF_Transport/Day_1/Inputs/example1/HartreeLeadsBottom/translated/single-points/v-0_bottom-leads-minimum-z-1']

data_bulk = []
data_em = []
for i in range(len(folder)):
    data_bulk.append(np.genfromtxt('{}/{}'.format(folder[i], bulk), skip_header=0, skip_footer=0))
    data_em.append(np.genfromtxt('{}/{}'.format(folder[i], em), skip_header=0, skip_footer=0))

# Hartree potential
# fig_plot_1, ax_plot_1 = plt.subplots()
# ax_plot_1.plot(data_em[0][:, 0], data_em[0][:, 1], '-', color='grey', label=labels_em[0], linewidth=6, alpha=0.5)
# ax_plot_1.plot(data_bulk[0][:, 0], data_bulk[0][:, 1], 'k-', label=labels_bulk[0])
# ax_plot_1.legend(frameon=False)
# ax_plot_1.set_ylim([ylim[0], ylim[1]])
# ax_plot_1.set_xlabel(r'Position / Å')
# ax_plot_1.set_ylabel('Hartree potential / eV')
# fig_plot_1.tight_layout()
# for i in range(len(folder)):
#     fig_plot_1.savefig('{}/hartree_potential_basic.png'.format(folder[i]), dpi=param.save_dpi)

# Hartree potential
fig_plot_2, ax_plot_2 = plt.subplots()
ax_plot_2.plot(data_em[0][:, 0], data_em[0][:, 1], '-', color='grey', label=labels_em[0], linewidth=6, alpha=0.5)
# ax_plot_2.plot(data_em[1][:, 0], data_em[1][:, 1], 'g-', label=labels_em[1])
# ax_plot_2.plot(data_em[4][:, 0], data_em[4][:, 1], 'b-', label=labels_em[4])
# ax_plot_2.plot(data_em[5][:, 0], data_em[5][:, 1], 'm-', label=labels_em[5])
ax_plot_2.plot(data_bulk[0][:, 0], data_bulk[0][:, 1], 'k-', label=labels_bulk[0])
ax_plot_2.legend(frameon=False)
ax_plot_2.set_ylim([ylim[0], ylim[1]])
ax_plot_2.set_xlabel(r'Position / Å')
ax_plot_2.set_ylabel('Hartree potential / eV')
fig_plot_2.tight_layout()
for i in range(len(folder)):
    fig_plot_2.savefig('{}/hartree_potential.png'.format(folder[i]), dpi=param.save_dpi)
    
# Hartree potential under bias
# fig_plot_3, ax_plot_3 = plt.subplots()
# ax_plot_3.plot(data_em[2][:, 0], data_em[2][:, 1], 'r-', label=labels_em[2])
# ax_plot_3.plot(data_em[3][:, 0], data_em[3][:, 1], 'g-', label=labels_em[3])
# ax_plot_3.plot(data_bulk[0][:, 0], data_bulk[0][:, 1], 'k-', label=labels_bulk[0])
# ax_plot_3.legend(frameon=False)
# ax_plot_3.set_xlabel(r'Position / Å')
# ax_plot_3.set_ylabel('Hartree potential / eV')
# fig_plot_3.tight_layout()
# for i in range(len(folder)):
#     fig_plot_3.savefig('{}/hartree_potential_bias.png'.format(folder[i]), dpi=param.save_dpi)

# Translated
# fig_plot_1, ax_plot_1 = plt.subplots()
# ax_plot_1.plot(data_em[0][:, 0], data_em[0][:, 1], '-', color='grey', label=labels_em[0], linewidth=6, alpha=0.5)
# ax_plot_1.plot(data_em[2][:, 0], data_em[2][:, 1], '-', color='brown', label=labels_em[2], linewidth=6, alpha=0.5)
# ax_plot_1.plot(data_em[1][:, 0], data_em[1][:, 1], 'r-', label=labels_em[1])
# ax_plot_1.plot(data_em[3][:, 0], data_em[3][:, 1], 'g-', label=labels_em[3])
# ax_plot_1.plot(data_em[5][:, 0], data_em[5][:, 1], '-', color='orange', label=labels_em[5])
# ax_plot_1.plot(data_bulk[0][:, 0], data_bulk[0][:, 1], 'b-', label=labels_bulk[0])
# ax_plot_1.plot(data_bulk[2][:, 0], data_bulk[2][:, 1], 'm-', label=labels_bulk[2])
# ax_plot_1.legend(frameon=False)
# ax_plot_1.set_ylim([ylim[0], ylim[1]])
# ax_plot_1.set_xlabel(r'Position / Å')
# ax_plot_1.set_ylabel('Hartree potential / eV')
# fig_plot_1.tight_layout()
# for i in range(len(folder)):
#     fig_plot_1.savefig('{}/hartree_potential_translate.png'.format(folder[i]), dpi=param.save_dpi)

if __name__ == "__main__":
    print('Finished.')
    plt.show()
