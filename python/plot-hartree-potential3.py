import numpy as np
import matplotlib.pyplot as plt
from scripts.main import parameters as param

""" Plotting of SMEAGOL output _TRC.agr by filename"""

# CP2K Lithium 1D wire
bulk = 'bulk-VH_AV.dat'
em_1 = '0V-VH_z-0001.dat'
em_2 = '0V-VH_z-0002.dat'
em_3 = '0V-VH_z-0003.dat'
em_4 = '0V-VH_z-0004.dat'
em_7 = '0V-VH_z-0007.dat'
em_30 = '0V-VH_z-0030.dat'
em_orig_1 = '0V-VH_z_ORIG-0001.dat'
em_orig_2 = '0V-VH_z_ORIG-0002.dat'
em_orig_3 = '0V-VH_z_ORIG-0003.dat'
em_orig_4 = '0V-VH_z_ORIG-0004.dat'
em_orig_7 = '0V-VH_z_ORIG-0007.dat'
em_orig_30 = '0V-VH_z_ORIG-0030.dat'

labels_em = ['EM']
labels_bulk = ['Bulk']
ylim = [-3, 2]
folder = [
    '/Volumes/Storage/Data/delete/Work/Postdoc/Work/calculations/transport/iv/cp2k/transmission/dev/V-0_HLB-1_z-0-0']

data_bulk = []
data_em_1 = []
data_em_2 = []
data_em_3 = []
data_em_4 = []
data_em_7 = []
data_em_30 = []
data_em_orig_1 = []
data_em_orig_2 = []
data_em_orig_3 = []
data_em_orig_4 = []
data_em_orig_7 = []
data_em_orig_30 = []
for i in range(len(folder)):
    data_bulk.append(np.genfromtxt('{}/{}'.format(folder[i], bulk), skip_header=0, skip_footer=0))
    data_em_1.append(np.genfromtxt('{}/{}'.format(folder[i], em_1), skip_header=0, skip_footer=0))
    data_em_2.append(np.genfromtxt('{}/{}'.format(folder[i], em_2), skip_header=0, skip_footer=0))
    data_em_3.append(np.genfromtxt('{}/{}'.format(folder[i], em_3), skip_header=0, skip_footer=0))
    data_em_4.append(np.genfromtxt('{}/{}'.format(folder[i], em_4), skip_header=0, skip_footer=0))
    data_em_7.append(np.genfromtxt('{}/{}'.format(folder[i], em_7), skip_header=0, skip_footer=0))
    data_em_30.append(np.genfromtxt('{}/{}'.format(folder[i], em_30), skip_header=0, skip_footer=0))
    data_em_orig_1.append(np.genfromtxt('{}/{}'.format(folder[i], em_orig_1), skip_header=0, skip_footer=0))
    data_em_orig_2.append(np.genfromtxt('{}/{}'.format(folder[i], em_orig_2), skip_header=0, skip_footer=0))
    data_em_orig_3.append(np.genfromtxt('{}/{}'.format(folder[i], em_orig_3), skip_header=0, skip_footer=0))
    data_em_orig_4.append(np.genfromtxt('{}/{}'.format(folder[i], em_orig_4), skip_header=0, skip_footer=0))
    data_em_orig_7.append(np.genfromtxt('{}/{}'.format(folder[i], em_orig_7), skip_header=0, skip_footer=0))
    data_em_orig_30.append(np.genfromtxt('{}/{}'.format(folder[i], em_orig_30), skip_header=0, skip_footer=0))

# Hartree potential
fig_plot_1, ax_plot_1 = plt.subplots()
ax_plot_1.plot(data_em_2[0][:, 0], data_em_2[0][:, 1], 'b-', label='EM NEGF 1')
# ax_plot_1.plot(data_em_3[0][:, 0], data_em_3[0][:, 1], 'y-', label='EM NEGF 2')
ax_plot_1.plot(data_em_4[0][:, 0], data_em_4[0][:, 1], 'c-', label='EM NEGF 3')
ax_plot_1.plot(data_em_7[0][:, 0], data_em_7[0][:, 1], 'm-', label='EM NEGF 7')
ax_plot_1.plot(data_em_30[0][:, 0], data_em_30[0][:, 1], 'r-', label='EM NEGF 30')
ax_plot_1.plot(data_em_orig_1[0][:, 0], data_em_orig_1[0][:, 1], 'g-', label='EM DFT')
ax_plot_1.plot(data_bulk[0][:, 0], data_bulk[0][:, 1], 'k-', label='Bulk')
ax_plot_1.legend(frameon=False)
# ax_plot_1.set_ylim([ylim[0], ylim[1]])
ax_plot_1.set_xlabel(r'Position / Å')
ax_plot_1.set_ylabel('Hartree potential / eV')
fig_plot_1.tight_layout()
for i in range(len(folder)):
    fig_plot_1.savefig('{}/hartree_potential_basic.png'.format(folder[i]), dpi=param.save_dpi)

# Hartree potential
# fig_plot_2, ax_plot_2 = plt.subplots()
# ax_plot_2.plot(data_em[0][:, 0], data_em[0][:, 1], '-', color='grey', label=labels_em[0], linewidth=6, alpha=0.5)
# ax_plot_2.plot(data_em[1][:, 0], data_em[1][:, 1], 'g-', label=labels_em[1])
# ax_plot_2.plot(data_em[4][:, 0], data_em[4][:, 1], 'b-', label=labels_em[4])
# ax_plot_2.plot(data_em[5][:, 0], data_em[5][:, 1], 'm-', label=labels_em[5])
# ax_plot_2.plot(data_bulk[0][:, 0], data_bulk[0][:, 1], 'k-', label=labels_bulk[0])
# ax_plot_2.legend(frameon=False)
# ax_plot_2.set_ylim([ylim[0], ylim[1]])
# ax_plot_2.set_xlabel(r'Position / Å')
# ax_plot_2.set_ylabel('Hartree potential / eV')
# fig_plot_2.tight_layout()
# for i in range(len(folder)):
#     fig_plot_2.savefig('{}/hartree_potential.png'.format(folder[i]), dpi=param.save_dpi)
    
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
