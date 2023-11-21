import numpy as np
import matplotlib.pyplot as plt
from general import parameters as param

""" Plotting of SMEAGOL output _TRC.agr by filename"""

# plotting_colors = ['k', 'r', 'g', 'b', 'm', 'grey', 'orange', 'y']
plotting_colors = ['r', 'g', 'b', 'm', 'grey', 'orange', 'y']
# plotting_colors_2 = ['orange', 'm', 'grey', 'b', 'g', 'r']
plotting_colors_2 = ['b', 'orange', 'm', 'grey',  'g', 'r']
n = 1

# Plot Li chain LDA SIESTA:q1 CP2K:q3 27 atoms
# xlim = [-5, 10]
# ylim = [-0.1, 1.1]
# ylim_log = [0.008, 1.1]
# ylim_dos = [0, 50]
# labels = ['CP2K-SMEAGOL', 'SIESTA-SMEAGOL']
# fermi = np.zeros(len(labels))
# folder = ['/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/li-chain/cp2k-smeagol/transmission/V-0_HLB-F_z-0-0/output',
#           '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/li-chain/siesta/transmission/v-0_bottom-0/output']
# labels_cp2k_negf = ['CP2K-NEGF']
# fermi_cp2k_negf = [-0.09983271872231]
# dos_norm_cp2k_negf = 28  # Number of atoms
# folder_cp2k_negf = ['/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/li-chain/cp2k-negf/transmission/V-0_fermi-calc_xy-12/output']

# Plot Au chain LDA q11 27 atoms
# xlim = [-3, 3]
# ylim = [0.0, 6.1]
# ylim_log = [0.008, 6.1]
# ylim_dos = [0, 500]
# labels = ['CP2K-SMEAGOL', 'SIESTA-SMEAGOL']
# fermi = np.zeros(len(labels))
# folder = ['/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-chain/cp2k-smeagol/transmission/V-0_HLB-F_z-0-0_atoms-28_noprint/output',
#           '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-chain/siesta/transmission/V-0_HLB-0_z-0-0/output']
# labels_cp2k_negf = ['CP2K-NEGF']
# fermi_cp2k_negf = [-0.22124430588176]
# dos_norm_cp2k_negf = 28
# folder_cp2k_negf = ['/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-chain/cp2k-negf/transmission/V-0_multiple-force-eval_fermi-calc_kpoints-1-1-31_xy-12/output']

# # Plot Au-BDT for experimental CP2K HLB=F and SIESTA HLB=T
# xlim = [-5, 5]
# ylim = [0.0, 1.0]
# ylim_log = [0.008, 1.2]
# ylim_dos = [0, 610]
# labels = ['CP2K', 'SIESTA']
# fermi = np.zeros(len(labels))
# folder = ['/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/cp2k-smeagol/sz/transmission/exp/kpoints_bulk-2-2-100_em-2-2-1/output',
#           '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/siesta-smeagol/transmission/exp/bulk-4-4-100-em-4-4-1_hlb-15.2496/output']
# folder_cp2k_negf = []
# fermi_cp2k_negf = 0

# # Plot Au-BDT for SIESTA experimental HLB=F
# xlim = [-5, 5]
# ylim = [0.0, 1.0]
# ylim_log = [0.008, 1.2]
# ylim_dos = [0, 610]
# labels = ['SIESTA 2x2', 'SIESTA 3x3', 'SIESTA 4x4']
# fermi = np.zeros(len(labels))
# folder = ['/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/siesta-smeagol/transmission/exp/bulk-2-2-100-em-2-2-1_4/output',
#           '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/siesta-smeagol/transmission/exp/bulk-3-3-100-em-3-3-1/output',
#           '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/siesta-smeagol/transmission/exp/bulk-4-4-100-em-4-4-1/output']
# folder_cp2k_negf = []
# fermi_cp2k_negf = 0

# Plot Au-BDT for SIESTA experimental, optimised HLB=T
# xlim = [-3, 3]
# ylim = [0.0, 1.0]
# ylim_log = [0.009, 1.3]
# ylim_dos = [0, 610]
# labels = ['SIESTA experimental', 'SIESTA optimised']
# fermi = np.zeros(len(labels))
# folder = ['/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/siesta-smeagol/transmission/exp/bulk-4-4-100-em-4-4-1_hlb-15.2496/output',
#           '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/siesta-smeagol/transmission/opt-cg/bulk-4-4-100-em-4-4-1_hlb-15.412-0-0/output']
# folder_cp2k_negf = []
# fermi_cp2k_negf = 0

# Plot Au-BDT for experimental CP2K HLB=F and SIESTA HLB=T
xlim = [-4, 4]
ylim = [0.0, 1.0]
ylim_log = [0.008, 1.2]
ylim_dos = [0, 300]
labels = ['CP2K', 'SIESTA']
fermi = np.zeros(len(labels))
# fermi = [0, 0.0, 0.7]
folder = ['/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/cp2k-smeagol/sz/transmission/exp/dzvp/sergey/kpoints_bulk-4-4-100_em-4-4-1_hlb-t-11.03197_scf-500/output',
          '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/siesta-smeagol/transmission/exp/bulk-4-4-100-em-4-4-1_hlb-15.2496/output']
folder_cp2k_negf = []
fermi_cp2k_negf = 0
plot_lengend = True

# Au capacitor
# xlim = [-4, 4]
# ylim = [0.0, 1.0]
# ylim_log = [0.008, 1.2]
# ylim_dos = [0, 300]
# labels = ['CP2K', 'SIESTA']
# fermi = np.zeros(len(labels))
# # fermi = [0, 0.0, 0.7]
# folder = ['/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/cp2k-smeagol/sz/transmission/exp/capacitor/sergey-equal/kpoints_bulk-4-4-100_em-4-4-1_hlb-t-10.99872_scf-500/output',
#           '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/siesta-smeagol/capacitor/bulk-4-4-100-em-4-4-1_hlb-15.2496_0-0/output']
# folder_cp2k_negf = []
# fermi_cp2k_negf = 0

# Au capacitor
# xlim = [-4, 4]
# ylim = [0.0, 1.0]
# ylim_log = [0.008, 1.2]
# ylim_dos = [0, 300]
# labels = ['Delta 0', 'Delta 1e-4']
# fermi = np.zeros(len(labels))
# # fermi = [0, 0.0, 0.7]
# folder = ['/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/siesta-smeagol/capacitor/testing/au-c3/kpoints-4-4-20_hlb-auto_cores-64_restricted_delta-0/output',
#           '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/siesta-smeagol/capacitor/testing/au-c3/kpoints-4-4-20_hlb-auto_cores-64_restricted_delta-1e-4/output']
# folder_cp2k_negf = []
# fermi_cp2k_negf = 0
# plot_lengend = True

# cp2k-smeagol-examples/examples/li-chain
# xlim = [-4, 6]
# ylim = [-0.1, 1.1]
# ylim_log = [0.008, 1.1]
# ylim_dos = [0, 50]
# labels = ['CP2K+SMEAGOL', 'SIESTA1+SMEAGOL']
# fermi = np.zeros(len(labels))
# # fermi = [0, 0.0, 0.7]
# folder = ['/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/li-chain/cp2k-smeagol/transmission/kpoints-1-1-20_20220907-memory_hash-143/output',
#           '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/li-chain/siesta1-smeagol/transmission/kpoints-1-1-20/output']
# folder_cp2k_negf = []
# fermi_cp2k_negf = 0
# plot_lengend = True

file_1 = []
file_2 = []
file_3 = []
file_4 = []
file_5 = []
for i in range(len(folder)):
    file_1.append(np.genfromtxt('{}/G0-S0.out'.format(folder[i]), skip_header=2, skip_footer=3) * n)
    file_2.append(np.genfromtxt('{}/G1-S0.out'.format(folder[i]), skip_header=2, skip_footer=3) * n)
    file_3.append(np.genfromtxt('{}/G1-S2.out'.format(folder[i]), skip_header=2, skip_footer=3) * n)
    file_4.append(np.genfromtxt('{}/G2-S0.out'.format(folder[i]), skip_header=2, skip_footer=3) * n)
    file_5.append(np.genfromtxt('{}/G3-S0.out'.format(folder[i]), skip_header=2, skip_footer=1) * n)
print('number of SMEAGOL files', len(folder))

file_1_negf = []
file_4_negf = []
for i in range(len(folder_cp2k_negf)):
    file_1_negf.append(np.genfromtxt('{}/G0-S0.out'.format(folder_cp2k_negf[i]), skip_header=2, skip_footer=3) * n)
    file_4_negf.append(np.genfromtxt('{}/G2-S0.out'.format(folder_cp2k_negf[i]), skip_header=2, skip_footer=3) * n)
print('number of CP2K-NEGF files', len(folder_cp2k_negf))

# Transmission
fig_plot_1, ax_plot_1 = plt.subplots()
for i in range(len(folder_cp2k_negf)):
    ax_plot_1.plot((file_1_negf[i][:, 0] - fermi_cp2k_negf[i]) * param.hartree_to_ev, file_1_negf[i][:, 1], color=plotting_colors_2[i], label=labels_cp2k_negf[i])
for i in range(len(folder)):
    ax_plot_1.plot(file_1[i][:, 0] - fermi[i], file_1[i][:, 1], color=plotting_colors[i], label=labels[i])
ax_plot_1.set_xlim([xlim[0], xlim[1]])
ax_plot_1.set_ylim([ylim[0], ylim[1]])
ax_plot_1.legend(frameon=False)
ax_plot_1.set_xlabel(r'E-E$_{\mathrm{F}}$ (eV)')
ax_plot_1.set_ylabel('Transmission')
fig_plot_1.tight_layout()
for i in range(len(folder_cp2k_negf)):
    fig_plot_1.savefig('{}/compare-transmission.png'.format(folder_cp2k_negf[i]), dpi=param.save_dpi)
for i in range(len(folder)):
    fig_plot_1.savefig('{}/compare-transmission.png'.format(folder[i]), dpi=param.save_dpi)

# Transmission log
fig_plot_1_2, ax_plot_1_2 = plt.subplots()
for i in range(len(folder)):
    ax_plot_1_2.plot(file_1[i][:, 0]+fermi[i], file_1[i][:, 1], color=plotting_colors[i], label=labels[i])
ax_plot_1_2.set_xlim([xlim[0], xlim[1]])
ax_plot_1_2.set_ylim([ylim_log[0], ylim_log[1]])
ax_plot_1_2.set_yscale('log')
ax_plot_1_2.set_xlabel(r'E-E$_{\mathrm{F}}$ (eV)')
ax_plot_1_2.set_ylabel('Transmission')
ax_plot_1_2.legend(frameon=False)
fig_plot_1_2.tight_layout()
for i in range(len(folder)):
    fig_plot_1_2.savefig('{}/compare-transmission_log.png'.format(folder[i]), dpi=param.save_dpi)

# Plot transmission and transmission log
# rows, cols = 2, 1
# fig1_3, ax1_3 = plt.subplots(rows, cols, sharex='col', sharey='row', figsize=(6, 8))
# for i in range(len(folder)):
#     ax1_3[0].plot(file_1[i][:, 0]+fermi[i], file_1[i][:, 1], color=plotting_colors[i], label=labels[i])
# ax1_3[0].set_xlim([xlim[0], xlim[1]])
# ax1_3[0].set_ylim([ylim[0], ylim[1]])
# ax1_3[0].legend(frameon=False)
# ax1_3[0].set_ylabel('Transmission')
# for i in range(len(folder)):
#     ax1_3[1].plot(file_1[i][:, 0]+fermi[i], file_1[i][:, 1], color=plotting_colors[i], label=labels[i])
# ax1_3[1].set_xlim([xlim[0], xlim[1]])
# ax1_3[1].set_ylim([ylim_log[0], ylim_log[1]])
# ax1_3[1].legend(frameon=False)
# ax1_3[1].set_yscale('log')
# ax1_3[1].set_xlabel(r'E-E$_{\mathrm{F}}$ (eV)')
# ax1_3[1].set_ylabel('Log transmission')
# fig1_3.tight_layout()
# for i in range(len(folder)):
#     fig1_3.savefig('{}/compare_transmission_log_all.png'.format(folder[i]), dpi=param.save_dpi)

# Plot EM DOS
fig_plot_2, ax_plot_2 = plt.subplots()
for i in range(len(folder_cp2k_negf)):
    ax_plot_2.plot((file_4_negf[i][:, 0] - fermi_cp2k_negf[i]) * param.hartree_to_ev, file_4_negf[i][:, 1] / dos_norm_cp2k_negf, color=plotting_colors_2[i], label=labels_cp2k_negf[i])
for i in range(len(folder)):
    ax_plot_2.plot(file_4[i][:, 0] + fermi[i], file_4[i][:, 1], color=plotting_colors[i], label=labels[i])
ax_plot_2.set_xlim([xlim[0], xlim[1]])
ax_plot_2.set_ylim([ylim_dos[0], ylim_dos[1]])
ax_plot_2.legend(frameon=False)
ax_plot_2.set_xlabel(r'E-E$_{\mathrm{F}}$ (eV)')
# ax_plot_2.set_ylabel('Density of states (a.u.)')
ax_plot_2.set_ylabel(r'Density of states (atom$^{-1}$ eV$^{-1}$)')
fig_plot_2.tight_layout()
# fig_trans_dos.subplots_adjust(hspace=0)
for i in range(len(folder_cp2k_negf)):
    fig_plot_2.savefig('{}/compare-em_dos.png'.format(folder_cp2k_negf[i]), dpi=param.save_dpi)
for i in range(len(folder)):
    fig_plot_2.savefig('{}/compare-em_dos.png'.format(folder[i]), dpi=param.save_dpi)
    
# Plot transmission and EM DOS
rows, cols = 2, 1
fig_trans_dos, ax_trans_dos = plt.subplots(rows, cols, sharex='col', sharey='row', figsize=(5, 8))
for i in range(len(folder_cp2k_negf)):
    ax_trans_dos[0].plot((file_1_negf[i][:, 0] - fermi_cp2k_negf[i]) * param.hartree_to_ev, file_1_negf[i][:, 1], color=plotting_colors_2[i], label=labels_cp2k_negf[i])
for i in range(len(folder)):
    ax_trans_dos[0].plot(file_1[i][:, 0]+fermi[i], file_1[i][:, 1], color=plotting_colors[i], label=labels[i])
ax_trans_dos[0].set_xlim([xlim[0], xlim[1]])
# ax_trans_dos[0].set_ylim([ylim[0], ylim[1]])
ax_trans_dos[0].set_ylabel('Transmission')
if plot_lengend: ax_trans_dos[0].legend(frameon=False)
for i in range(len(folder_cp2k_negf)):
    ax_trans_dos[1].plot((file_4_negf[i][:, 0] - fermi_cp2k_negf[i]) * param.hartree_to_ev, file_4_negf[i][:, 1] / dos_norm_cp2k_negf, color=plotting_colors_2[i], label=labels_cp2k_negf[i])
for i in range(len(folder)):
    ax_trans_dos[1].plot(file_4[i][:, 0] + fermi[i], file_4[i][:, 1], color=plotting_colors[i], label=labels[i])
ax_trans_dos[1].set_xlim([xlim[0], xlim[1]])
# ax_trans_dos[1].set_ylim([ylim_dos[0], ylim_dos[1]])
if plot_lengend: ax_trans_dos[1].legend(frameon=False)
ax_trans_dos[1].set_xlabel(r'E-E$_{\mathrm{F}}$ (eV)')
# ax_trans_dos[1].set_ylabel('Density of states (a.u.)')
ax_trans_dos[1].set_ylabel(r'Density of states (atom$^{-1}$ eV$^{-1}$)')
fig_trans_dos.tight_layout()
# fig_trans_dos.subplots_adjust(hspace=0)
for i in range(len(folder_cp2k_negf)):
    fig_trans_dos.savefig('{}/compare-transmission.png'.format(folder_cp2k_negf[i]), dpi=param.save_dpi)
for i in range(len(folder)):
    fig_trans_dos.savefig('{}/compare-transmission_emdos.png'.format(folder[i]), dpi=param.save_dpi)

# Plot all
rows, cols = 2, 2
fig_plot_all, ax_plot_all = plt.subplots(rows, cols, figsize=(10, 8))
# fig_plot_all, ax_plot_all = plt.subplots(rows, cols, sharex='col', sharey='row', figsize=(10, 8))
for i in range(len(folder)):
    ax_plot_all[0, 0].plot(file_1[i][:, 0]+fermi[i], file_1[i][:, 1], color=plotting_colors[i], label=labels[i])
ax_plot_all[0, 0].set_xlim([xlim[0], xlim[1]])
ax_plot_all[0, 0].set_ylim([ylim[0], ylim[1]])
ax_plot_all[0, 0].legend(frameon=False)
# ax_plot_all[0, 0].set_xlabel(r'E-E$_{\mathrm{F}}$ (eV)')
ax_plot_all[0, 0].set_ylabel('Transmission')
for i in range(len(folder)):
    ax_plot_all[0, 1].plot(file_2[i][:, 0]+fermi[i], file_2[i][:, 1], color=plotting_colors[i], label=labels[i])
ax_plot_all[0, 1].set_xlim([xlim[0], xlim[1]])
ax_plot_all[0, 1].legend(frameon=False)
# ax_plot_all[0, 1].set_xlabel(r'E-E$_{\mathrm{F}}$ (eV)')
ax_plot_all[0, 1].set_ylabel('Number of channels')
for i in range(len(folder)):
    ax_plot_all[1, 0].plot(file_4[i][:, 0]+fermi[i], file_4[i][:, 1], color=plotting_colors[i], label=labels[i])
ax_plot_all[1, 0].set_xlim([xlim[0], xlim[1]])
ax_plot_all[1, 0].legend(frameon=False)
ax_plot_all[1, 0].set_xlabel(r'E-E$_{\mathrm{F}}$ (eV)')
ax_plot_all[1, 0].set_ylabel('EM DOS')
# for i in range(len(folder)):
#     ax_plot_all[1, 1].plot(file_5[i][:, 0]+fermi[i], file_5[i][:, 1], color=plotting_colors[i], label=labels[i])
# ax_plot_all[1, 1].set_xlim([xlim[0], xlim[1]])
# ax_plot_all[1, 1].legend(frameon=False)
# ax_plot_all[1, 1].set_xlabel(r'E-E$_{\mathrm{F}}$ (eV)')
# ax_plot_all[1, 1].set_ylabel('Leads DOS')
fig_plot_all.tight_layout()
# fig_plot_all.subplots_adjust(hspace=0)
for i in range(len(folder)):
    fig_plot_all.savefig('{}/compare-plot_all.png'.format(folder[i]), dpi=param.save_dpi)

if __name__ == "__main__":
    print('Finished.')
    plt.show()
