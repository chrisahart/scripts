import numpy as np
import matplotlib.pyplot as plt
from general import parameters as param

""" Plotting of SMEAGOL output _TRC.agr by filename"""

plotting_colors = ['r', 'g', 'b', 'm', 'grey', 'orange', 'y']
n = 1

# Plot Li comparison for HartreeBottomLeads 0V
# xlim = [-5, 10]  # Li and Au chain
# ylim = [-0.1, 1.1]
# ylim_log = [0.008, 1.1]
# labels = ['SIESTA Δ=0', 'SIESTA Δ=-0.38', 'CP2K align_Vhartree=F ', 'CP2K Δ=0', 'CP2K  Δ=-0.20']
# fermi_offset = np.zeros(len(labels))
# folder = ['/Volumes/Storage/Data/Work/Postdoc/Work/calculations/transport/iv/siesta/single-points/v-0_bottom-0/output',
#           '/Volumes/Storage/Data/Work/Postdoc/Work/calculations/transport/iv/siesta/single-points/v-0_bottom-0.38727754/output',
#           '/Volumes/Storage/Data/Work/Postdoc/Work/calculations/transport/iv/cp2k/transmission/0V-tidy-a9b-testing-1/output',
#           '/Volumes/Storage/Data/Work/Postdoc/Work/calculations/transport/iv/cp2k/transmission/0V-tidy-a9b_bottom-0_z-0-0/output',
#           '/Volumes/Storage/Data/Work/Postdoc/Work/calculations/transport/iv/cp2k/transmission/0V-tidy-a9b_bottom-0.2_z-0-0-ev/output']

# Plot Li comparison for HartreeBottomLeads 1V
xlim = [-5, 10]  # Li and Au chain
ylim = [-0.1, 1.1]
ylim_log = [0.008, 1.1]
labels = ['SIESTA Δ=0', 'SIESTA Δ=-0.38', 'CP2K align_Vhartree=F ', 'CP2K Δ=0', 'CP2K  Δ=-0.20']
fermi_offset = np.zeros(len(labels))
folder = ['/Volumes/Storage/Data/Work/Postdoc/Work/calculations/transport/iv/siesta/single-points/v-1_bottom-0/output',
          '/Volumes/Storage/Data/Work/Postdoc/Work/calculations/transport/iv/siesta/single-points/v-1_bottom-0.38727754/output',
          '/Volumes/Storage/Data/Work/Postdoc/Work/calculations/transport/iv/cp2k/transmission/1V-tidy-a9b-testing-1/output',
          '/Volumes/Storage/Data/Work/Postdoc/Work/calculations/transport/iv/cp2k/transmission/1V-tidy-a9b_bottom-0_z-0-0/output',
          '/Volumes/Storage/Data/Work/Postdoc/Work/calculations/transport/iv/cp2k/transmission/1V-tidy-a9b_bottom-0.2_z-0-0-ev/output']

# Plot Li comparison for CP2K and SIESTA
# xlim = [-5, 10]  # Li and Au chain
# ylim = [-0.1, 1.1]
# ylim_log = [0.008, 1.1]
# labels = ['CP2K', 'SIESTA']
# fermi_offset = np.zeros(len(labels))
# folder = ['/Volumes/Elements/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/1D-chains/li-chain/cp2k/0V/output',
#           '/Volumes/Elements/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/1D-chains/li-chain/siesta/0V/output']
# folder = ['/Volumes/Elements/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/testing/smeagol/cx1/Smeagol_Tutorial_Files/DFT_NEGF_Transport/Day_1/Inputs/testing/li_chain/li-sz-explicit-2s-dos-12atoms-2/output',
#           '/Volumes/Elements/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/testing/cp2k-smeagol/examples/SMEAGOL_TEST/Au_chain/transport-12-wfnrs-li-bulk-4-tidy-3/output']

# Plot Li and Au comparison for CP2K and SIESTA (12 atoms, SIESTA Li 2s valence vs Au 6s valence)
# xlim = [-5, 10]  # Li and Au chain
# ylim = [-0.1, 1.1]
# ylim_log = [0.008, 1.1]
# labels = ['CP2K Li', 'SIESTA Li', 'SIESTA Au (6s)']
# fermi_offset = np.zeros(len(labels))
# folder = ['/Volumes/Elements/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/1D-chains/li-chain/cp2k/0V/output',
#           '/Volumes/Elements/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/1D-chains/li-chain/siesta/0V/output']
# folder = ['/Volumes/Elements/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/testing/smeagol/cx1/Smeagol_Tutorial_Files/DFT_NEGF_Transport/Day_1/Inputs/testing/li_chain/li-sz-explicit-2s-dos-12atoms-2/output',
#           '/Volumes/Elements/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/testing/cp2k-smeagol/examples/SMEAGOL_TEST/Au_chain/transport-12-wfnrs-li-bulk-4-tidy-3/output',
#           '/Volumes/Elements/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/testing/smeagol/cx1/Smeagol_Tutorial_Files/DFT_NEGF_Transport/Day_1/Inputs/testing/au_chain-12-4atoms/run/Au-lda-ivan-sz-4s-def-lda/output']

# Plot Au-Au comparison for CP2K and SIESTA (12 atoms, SIESTA short and lone range pseudopotential)
# xlim = [-5, 10]  # Li and Au chain
# ylim = [-0.1, 6.1]
# ylim_log = [0.008, 6.1]
# labels = ['CP2K', 'SIESTA (2)', 'SIESTA (8)']
# fermi_offset = np.zeros(len(labels))
# folder = ['/Volumes/Elements/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/testing/cp2k-smeagol/examples/SMEAGOL_TEST/Au_chain/transport-12-AuAu-2p8/output',
#           '/Volumes/Elements/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/testing/smeagol/cx1/Smeagol_Tutorial_Files/DFT_NEGF_Transport/Day_1/Inputs/testing/au_chain-12-4atoms/run/Au-lda-leem-sz-4s-3d-cutoff-lda/output',
#           '/Volumes/Elements/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/testing/smeagol/cx1/Smeagol_Tutorial_Files/DFT_NEGF_Transport/Day_1/Inputs/testing/au_chain-12-4atoms/run/Au-lda-clop-sz-4s-3d-cutoff-lda/output']

# Plot Au-Au comparison for CP2K and SIESTA (12 atoms, SIESTA basis set cutoff)
# xlim = [-5, 10]  # Li and Au chain
# ylim = [-0.1, 6.1]
# ylim_log = [0.008, 6.1]
# labels = ['CP2K', 'SIESTA (8) default Rc', 'SIESTA (8) longer Rc']
# fermi_offset = np.zeros(len(labels))
# folder = ['/Volumes/Elements/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/testing/cp2k-smeagol/examples/SMEAGOL_TEST/Au_chain/transport-12-AuAu-2p8/output',
#           '/Volumes/Elements/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/testing/smeagol/cx1/Smeagol_Tutorial_Files/DFT_NEGF_Transport/Day_1/Inputs/testing/au_chain-12-4atoms/run/Au-lda-clop-sz-4s-3d-def-lda/output',
#           '/Volumes/Elements/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/testing/smeagol/cx1/Smeagol_Tutorial_Files/DFT_NEGF_Transport/Day_1/Inputs/testing/au_chain-12-4atoms/run/Au-lda-clop-sz-4s-3d-cutoff-lda/output']

# Plot Au-Au comparison for CP2K and SIESTA (12 atoms, CP2K q11 and q19 valence)
# xlim = [-5, 10]  # Li and Au chain
# ylim = [-0.1, 6.1]
# ylim_log = [0.008, 6.1]
# labels = ['CP2K LDA q11', 'CP2K PBE q11', 'CP2K PBE q19']
# fermi_offset = np.zeros(len(labels))
# folder = ['/Volumes/Elements/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/testing/cp2k-smeagol/examples/SMEAGOL_TEST/Au_chain/transport-12-AuAu-2p8-lda/output',
#           '/Volumes/Elements/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/testing/cp2k-smeagol/examples/SMEAGOL_TEST/Au_chain/transport-12-AuAu-2p8-pbe/output',
#           '/Volumes/Elements/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/testing/cp2k-smeagol/examples/SMEAGOL_TEST/Au_chain/transport-12-AuAu-2p8-pbe-q19/output']

# Plot Au-Au comparison for CP2K and SIESTA (12 atoms, SIESTA sd and sdp valence)
# xlim = [-5, 10]
# ylim = [-0.1, 6.1]
# ylim_log = [0.008, 6.1]
# labels = ['CP2K', 'SIESTA 6s, 5d', 'SIESTA 6s, 5d, 6p']
# fermi_offset = np.zeros(len(labels))
# folder = ['/Volumes/Elements/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/testing/cp2k-smeagol/examples/SMEAGOL_TEST/Au_chain/transport-12-AuAu-2p8-lda/output',
#           '/Volumes/Elements/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/1D-chains/au-chain/scarf/siesta/au_chain-12-4atoms/run/Au-lda-clop-sz-6s-5d-cutoff-lda/output',
#           '/Volumes/Elements/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/1D-chains/au-chain/scarf/siesta/au_chain-12-4atoms/run/Au-lda-clop-sz-6s-5d-6p-cutoff-lda/output']

# Plot Au-Au comparison for CP2K and SIESTA (12 atoms)
# xlim = [-3, 3]
# ylim = [0.0, 6.1]
# ylim_log = [0.008, 6.1]
# labels = ['CP2K', 'SIESTA']
# fermi_offset = np.zeros(len(labels))
# folder = ['/Volumes/Elements/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/testing/cp2k-smeagol/examples/SMEAGOL_TEST/Au_chain/transport-12-AuAu-2p8/output',
#           '/Volumes/Elements/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/testing/smeagol/cx1/Smeagol_Tutorial_Files/DFT_NEGF_Transport/Day_1/Inputs/testing/au_chain-12-4atoms/run/Au-lda-clop-sz-4s-3d-cutoff-lda/output']

# Plot Au-Au comparison for CP2K and SIESTA (21 atoms)
# xlim = [-3, 3]  # Au chain
# ylim = [0.0, 6.1]
# ylim_log = [0.008, 6.1]
# labels = ['CP2K', 'SIESTA']
# fermi_offset = np.zeros(len(labels))
# folder = ['/Volumes/Elements/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/testing/cp2k-smeagol/examples/SMEAGOL_TEST/Au_chain/transport-21-AuAu-2p8/output',
#           '/Volumes/Elements/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/testing/smeagol/cx1/Smeagol_Tutorial_Files/DFT_NEGF_Transport/Day_1/Inputs/testing/au_chain-21-7atoms/run/Au-lda-clop-sz-4s-3d-cutoff-lda/output']

# Plot Au-BDT old (4x4 Au, 7-6 layers)
# xlim = [-3, 3]
# ylim = [0.0, 1.0]
# ylim_log = [0.008, 1.2]
# labels = ['1x1x31, 1x1x1']
# fermi_offset = np.zeros(len(labels))
# folder = ['/Volumes/Storage/Data/Work/Postdoc/Work/calculations/testing/cp2k-smeagol/bdt/siesta/hollow-site/7-6/transport-siesta/bulk-5-4/dzp-lda-pp8/output']

# Plot Au-BDT with kpoint convergence (3x3 Au, 5-4 layers)
# xlim = [-3, 3]
# ylim = [0.0, 1.0]
# ylim_log = [0.008, 1.2]
# labels = ['1x1x31, 1x1x1', '2x2x100, 2x2x1', '4x4x100, 4x4x1']
# fermi_offset = np.zeros(len(labels))
# folder = ['/Volumes/Elements/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/new/au-bdt/young/au-bdt/3x3/kpoints/bulk-1-1-31-em-1-1-1/output',
#           '/Volumes/Elements/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/new/au-bdt/young/au-bdt/3x3/kpoints/bulk-2-2-100-em-2-2-1/output',
#           '/Volumes/Elements/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/new/au-bdt/young/au-bdt/3x3/kpoints/bulk-4-4-100-em-4-4-1/output']

# Plot Au-BDT with kpoint convergence (3x3 Au, 5-4 layers)
# xlim = [-3, 3]
# ylim = [0.0, 1.0]
# ylim_log = [0.008, 1.2]
# labels = ['1x1x31, 1x1x1', '1x1x100, 1x1x1', '2x2x100, 2x2x1', '4x4x100, 4x4x1']
# fermi_offset = np.zeros(len(labels))
# folder = ['/Volumes/Elements/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/new/au-bdt/young/au-bdt/3x3/kpoints/bulk-1-1-31-em-1-1-1/output',
#           '/Volumes/Elements/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/new/au-bdt/young/au-bdt/3x3/kpoints/bulk-2-2-100-em-2-2-1/output',
#           '/Volumes/Elements/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/new/au-bdt/young/au-bdt/3x3/kpoints/bulk-4-4-100-em-4-4-1/output']

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
print('number of files', len(folder))

# Transmission
fig_plot_1, ax_plot_1 = plt.subplots()
for i in range(len(folder)):
    ax_plot_1.plot(file_1[i][:, 0]+fermi_offset[i], file_1[i][:, 1], color=plotting_colors[i], label=labels[i])
ax_plot_1.set_xlim([xlim[0], xlim[1]])
ax_plot_1.set_ylim([ylim[0], ylim[1]])
ax_plot_1.legend(frameon=False)
ax_plot_1.set_xlabel(r'E-E$_{\mathrm{F}}$ (eV)')
ax_plot_1.set_ylabel('Transmission')
fig_plot_1.tight_layout()
for i in range(len(folder)):
    fig_plot_1.savefig('{}/compare-transmission.png'.format(folder[i]), dpi=param.save_dpi)

# Transmission log
fig_plot_1_2, ax_plot_1_2 = plt.subplots()
for i in range(len(folder)):
    ax_plot_1_2.plot(file_1[i][:, 0]+fermi_offset[i], file_1[i][:, 1], color=plotting_colors[i], label=labels[i])
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
rows, cols = 2, 1
fig1_3, ax1_3 = plt.subplots(rows, cols, sharex='col', sharey='row', figsize=(6, 8))
for i in range(len(folder)):
    ax1_3[0].plot(file_1[i][:, 0]+fermi_offset[i], file_1[i][:, 1], color=plotting_colors[i], label=labels[i])
ax1_3[0].set_xlim([xlim[0], xlim[1]])
ax1_3[0].set_ylim([ylim[0], ylim[1]])
ax1_3[0].legend(frameon=False)
ax1_3[0].set_ylabel('Transmission')
for i in range(len(folder)):
    ax1_3[1].plot(file_1[i][:, 0]+fermi_offset[i], file_1[i][:, 1], color=plotting_colors[i], label=labels[i])
ax1_3[1].set_xlim([xlim[0], xlim[1]])
ax1_3[1].set_ylim([ylim_log[0], ylim_log[1]])
ax1_3[1].legend(frameon=False)
ax1_3[1].set_yscale('log')
ax1_3[1].set_xlabel(r'E-E$_{\mathrm{F}}$ (eV)')
ax1_3[1].set_ylabel('Log transmission')
fig1_3.tight_layout()
for i in range(len(folder)):
    fig1_3.savefig('{}/compare_transmission_log_all.png'.format(folder[i]), dpi=param.save_dpi)

# Plot transmission and EM DOS
rows, cols = 2, 1
fig_trans_dos, ax_trans_dos = plt.subplots(rows, cols, sharex='col', sharey='row', figsize=(5, 8))
for i in range(len(folder)):
    ax_trans_dos[0].plot(file_1[i][:, 0]+fermi_offset[i], file_1[i][:, 1], color=plotting_colors[i], label=labels[i])
ax_trans_dos[0].set_xlim([xlim[0], xlim[1]])
ax_trans_dos[0].set_ylim([ylim[0], ylim[1]])
ax_trans_dos[0].set_ylabel('Transmission')
ax_trans_dos[0].legend(frameon=False)
for i in range(len(folder)):
    ax_trans_dos[1].plot(file_4[i][:, 0]+fermi_offset[i], file_4[i][:, 1], color=plotting_colors[i], label=labels[i])
ax_trans_dos[1].set_xlim([xlim[0], xlim[1]])
ax_trans_dos[1].legend(frameon=False)
ax_trans_dos[1].set_xlabel(r'E-E$_{\mathrm{F}}$ (eV)')
ax_trans_dos[1].set_ylabel('Density of states (a.u.)')
fig_trans_dos.tight_layout()
# fig_trans_dos.subplots_adjust(hspace=0)
for i in range(len(folder)):
    fig_trans_dos.savefig('{}/compare-transmission_emdos.png'.format(folder[i]), dpi=param.save_dpi)

# Plot all
rows, cols = 2, 2
fig_plot_all, ax_plot_all = plt.subplots(rows, cols, figsize=(10, 8))
# fig_plot_all, ax_plot_all = plt.subplots(rows, cols, sharex='col', sharey='row', figsize=(10, 8))
for i in range(len(folder)):
    ax_plot_all[0, 0].plot(file_1[i][:, 0]+fermi_offset[i], file_1[i][:, 1], color=plotting_colors[i], label=labels[i])
ax_plot_all[0, 0].set_xlim([xlim[0], xlim[1]])
ax_plot_all[0, 0].set_ylim([ylim[0], ylim[1]])
ax_plot_all[0, 0].legend(frameon=False)
# ax_plot_all[0, 0].set_xlabel(r'E-E$_{\mathrm{F}}$ (eV)')
ax_plot_all[0, 0].set_ylabel('Transmission')
for i in range(len(folder)):
    ax_plot_all[0, 1].plot(file_2[i][:, 0]+fermi_offset[i], file_2[i][:, 1], color=plotting_colors[i], label=labels[i])
ax_plot_all[0, 1].set_xlim([xlim[0], xlim[1]])
ax_plot_all[0, 1].legend(frameon=False)
# ax_plot_all[0, 1].set_xlabel(r'E-E$_{\mathrm{F}}$ (eV)')
ax_plot_all[0, 1].set_ylabel('Number of channels')
for i in range(len(folder)):
    ax_plot_all[1, 0].plot(file_4[i][:, 0]+fermi_offset[i], file_4[i][:, 1], color=plotting_colors[i], label=labels[i])
ax_plot_all[1, 0].set_xlim([xlim[0], xlim[1]])
ax_plot_all[1, 0].legend(frameon=False)
ax_plot_all[1, 0].set_xlabel(r'E-E$_{\mathrm{F}}$ (eV)')
ax_plot_all[1, 0].set_ylabel('EM DOS')
# for i in range(len(folder)):
#     ax_plot_all[1, 1].plot(file_5[i][:, 0]+fermi_offset[i], file_5[i][:, 1], color=plotting_colors[i], label=labels[i])
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
