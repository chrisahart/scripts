import numpy as np
import matplotlib.pyplot as plt
from general import parameters as param

""" Plotting of SMEAGOL output _TRC.agr by filename"""

# plotting_colors = ['r', 'g', 'b', 'm', 'grey']
plotting_colors = ['b', 'r', 'g', 'm', 'grey']
n = 1

# xlim = [-1, 1]  # Li and Au chain
# ylim = [0, 8e-5]

# labels = ['CP2K Li', 'SIESTA Li']
# cp2k_folder = '/Volumes/Storage/Data/Work/Postdoc/Work/calculations/transport/iv/cp2k/12-atoms-a9b-tidy-iv'
# siesta_folder = '/Volumes/Storage/Data/Work/Postdoc/Work/calculations/transport/iv/siesta/0V-tidy-iv'
# cp2k = np.genfromtxt('{}/IV.log'.format(cp2k_folder), skip_header=1, skip_footer=0)
# siesta = np.genfromtxt('{}/Liwire.CUR'.format(siesta_folder), skip_header=0, skip_footer=0)

# Old
# labels = ['CP2K Li', 'SIESTA Li', 'CP2K Li align_Vhartree=F']
# cp2k_folder1 = '/Volumes/Storage/Data/delete/Work/Postdoc/Work/calculations/transport/iv/cp2k/12-atoms-a9b-tidy-iv'
# cp2k1 = np.genfromtxt('{}/IV.log'.format(cp2k_folder1), skip_header=1, skip_footer=0)
# siesta_folder1 = '/Volumes/Storage/Data/delete/Work/Postdoc/Work/calculations/transport/iv/siesta/iv-bottom-0/'
# siesta1 = np.genfromtxt('{}/Liwire.CUR'.format(siesta_folder1), skip_header=0, skip_footer=0)

# Plot Li chain LDA SIESTA:q1 CP2K:q3 27 atoms
# labels = ['CP2K-NEGF', 'CP2K-SMEAGOL', 'SIESTA-SMEAGOL', 'CP2K-SMEAGOL-old']
# cp2k_folder1 = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/iv/li/cp2k-smeagol/lda-q3_atoms-28/iv/HLB-F_z-0-0'
# cp2k1 = np.genfromtxt('{}/IV.log'.format(cp2k_folder1), skip_header=1, skip_footer=0)
# cp2k_folder2 = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/iv/li/cp2k-smeagol/lda-q3_atoms-28/iv/HLB-F_z-0-0_sergey'
# cp2k2 = np.genfromtxt('{}/IV.log'.format(cp2k_folder2), skip_header=1, skip_footer=0)
# cp2k_negf_folder1 = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/iv/li/cp2k-negf/lda-q3_atoms-28/iv/V-0_fermi-calc_xy-12'
# cp2k_negf1 = np.genfromtxt('{}/IV.log'.format(cp2k_negf_folder1), skip_header=1, skip_footer=0)
# siesta_folder1 = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/iv/li/siesta/lda-q1_atoms-27/iv/bottom-0/'
# siesta1 = np.genfromtxt('{}/Liwire.CUR'.format(siesta_folder1), skip_header=0, skip_footer=0)
# xlim = [-1, 1]
# ylim = [-0.75e-4, 0.75e-4]

# # Plot Au chain LDA 27 atoms
# labels = ['CP2K-NEGF', 'CP2K-SMEAGOL', 'SIESTA-SMEAGOL', 'CP2K-SMEAGOL-old']
# cp2k_folder1 = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/iv/au-chain/cp2k-smeagol/iv/HLB-F_z-0-0'
# cp2k1 = np.genfromtxt('{}/IV.log'.format(cp2k_folder1), skip_header=1, skip_footer=0)
# cp2k_folder2 = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/iv/au-chain/cp2k-smeagol/iv/HLB-F_z-0-0'
# cp2k2 = np.genfromtxt('{}/IV.log'.format(cp2k_folder1), skip_header=1, skip_footer=0)
# cp2k_negf_folder1 = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/iv/au-chain/cp2k-negf/iv/multiple-force-eval_fermi-calc_kpoints-1-1-31_xy-12'
# cp2k_negf1 = np.genfromtxt('{}/IV.log'.format(cp2k_negf_folder1), skip_header=1, skip_footer=0)
# siesta_folder1 = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/iv/au-chain/siesta/iv/HLB-0_z-0-0_positive'
# siesta1 = np.genfromtxt('{}/Au.CUR'.format(siesta_folder1), skip_header=0, skip_footer=0)
# xlim = [-1, 1]
# ylim = [-0.75e-4, 0.75e-4]

# IV curve cp2k
# fig_plot_1, ax_plot_1 = plt.subplots()
# ax_plot_1.plot(cp2k_negf1[:, 0], cp2k_negf1[:, 1], '.-', color=plotting_colors[0], label=labels[0])
# ax_plot_1.plot(-cp2k_negf1[:, 0], -cp2k_negf1[:, 1], '.-', color=plotting_colors[0])
# ax_plot_1.plot(cp2k1[:, 0], cp2k1[:, 1], 'k-', label=labels[1])
# ax_plot_1.plot(-cp2k1[:, 0], -cp2k1[:, 1], 'k-')
# ax_plot_1.plot(siesta1[:, 0], siesta1[:, 1], '.-', color=plotting_colors[2], label=labels[2])
# ax_plot_1.plot(-siesta1[:, 0], -siesta1[:, 1], '.-', color=plotting_colors[2])
# ax_plot_1.set_xlim([xlim[0], xlim[1]])
# ax_plot_1.set_ylim([ylim[0], ylim[1]])
# ax_plot_1.legend(frameon=False)
# ax_plot_1.set_xlabel('Bias voltage / eV')
# ax_plot_1.set_ylabel('Current / A')
# fig_plot_1.tight_layout()
# fig_plot_1.savefig('{}/IV_cp2k.png'.format(cp2k_folder1), dpi=param.save_dpi)
# fig_plot_1.savefig('{}/IV.png'.format(cp2k_negf_folder1), dpi=param.save_dpi)
# fig_plot_1.savefig('{}/IV.png'.format(siesta_folder1), dpi=param.save_dpi)

# IV curve Au-BDT
labels = ['CP2K-SMEAGOL', 'SIESTA-SMEAGOL']
cp2k_folder1 = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/cp2k-smeagol/sz/iv/kpoints_bulk-2-2-100_em-2-2-1'
cp2k1 = np.genfromtxt('{}/IV.log'.format(cp2k_folder1), skip_header=1, skip_footer=0)
siesta_folder1 = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/siesta-smeagol/iv/exp/bulk-4-4-100-em-4-4-1_hlb-15.2496_0-25.59'
siesta1 = np.genfromtxt('{}/transport.CUR'.format(siesta_folder1), skip_header=0, skip_footer=0)
print(cp2k1)
print(siesta1)
xlim = [-5, 5]
ylim = [-60, 60]

# IV curve all
fig_plot_2, ax_plot_2 = plt.subplots()
# ax_plot_2.plot(cp2k_negf1[:, 0], cp2k_negf1[:, 1], '.-', color=plotting_colors[0], label=labels[0])
# ax_plot_2.plot(-cp2k_negf1[:, 0], -cp2k_negf1[:, 1], '.-', color=plotting_colors[0])
# ax_plot_2.plot(cp2k2[:, 0], cp2k2[:, 1], '.-', color=plotting_colors[3], label=labels[3])
# ax_plot_2.plot(-cp2k2[:, 0], -cp2k2[:, 1], '.-', color='k')
ax_plot_2.plot(cp2k1[:, 0], cp2k1[:, 1]* 1e6, '.-', color=plotting_colors[1], label=labels[0])
# ax_plot_2.plot(-cp2k1[:, 0], -cp2k1[:, 1], '.-', color=plotting_colors[1])
ax_plot_2.plot(siesta1[:, 0], siesta1[:, 1]* 1e6, '.-', color=plotting_colors[2], label=labels[1])
# ax_plot_2.plot(-siesta1[:, 0], -siesta1[:, 1], '.-', color=plotting_colors[2])
ax_plot_2.set_xlim([xlim[0], xlim[1]])
ax_plot_2.set_ylim([ylim[0], ylim[1]])
ax_plot_2.legend(frameon=False)
ax_plot_2.set_xlabel('Bias voltage / V')
# ax_plot_2.set_ylabel('Current / A')
ax_plot_2.set_ylabel('Current / μA')
fig_plot_2.tight_layout()
# fig_plot_2.savefig('{}/IV.png'.format(cp2k_folder2), dpi=param.save_dpi)
fig_plot_2.savefig('{}/IV.png'.format(cp2k_folder1), dpi=param.save_dpi)
# fig_plot_2.savefig('{}/IV.png'.format(cp2k_negf_folder1), dpi=param.save_dpi)
fig_plot_2.savefig('{}/IV.png'.format(siesta_folder1), dpi=param.save_dpi)


if __name__ == "__main__":
    print('Finished.')
    plt.show()
