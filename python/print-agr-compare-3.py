import numpy as np
import matplotlib.pyplot as plt
from general import parameters as param

""" Plotting of SMEAGOL output _TRC.agr by filename"""


# Plotting parameters
# xlim = [-5, 10]  # Li and Au chain
# labels = ['CP2K Li', 'SIESTA Li', 'SIESTA Au (4s)']
# folder_1 = '/Volumes/Storage/Data/Work/Postdoc/Work/testing/smeagol/cx1/Smeagol_Tutorial_Files/DFT_NEGF_Transport/Day_1/Inputs/testing/li_chain/li-sz-explicit-2s-dos-12atoms-2/output'
# folder_2 = '/Volumes/Storage/Data/Work/Postdoc/Work/testing/cp2k-smeagol/examples/SMEAGOL_TEST/Au_chain/transport-12-wfnrs-li-bulk-4-tidy-3/output'
# folder_3 = '/Volumes/Storage/Data/Work/Postdoc/Work/testing/smeagol/cx1/Smeagol_Tutorial_Files/DFT_NEGF_Transport/Day_1/Inputs/testing/au_chain/run/ivan-sz-4s/output'
# folder_3 = '/Volumes/Storage/Data/Work/Postdoc/Work/testing/smeagol/cx1/Smeagol_Tutorial_Files/DFT_NEGF_Transport/Day_1/Inputs/testing/li_chain/li-sz-explicit-2s-dos-12atoms-2-au-4s/output'

# xlim = [-3, 3]  # Au-BDT
# folder_3 = '/Volumes/Storage/Data/Work/Postdoc/Work/testing/smeagol/cx1/Smeagol_Tutorial_Files/DFT_NEGF_Transport/Day_1/Inputs/testing/au_chain/run/leem-sz-4s/output'
# folder_2 = '/Volumes/Storage/Data/Work/Postdoc/Work/testing/smeagol/cx1/Smeagol_Tutorial_Files/DFT_NEGF_Transport/Day_1/Inputs/testing/au_chain/run/leem-sz-default/output'
# folder_1 = '/Volumes/Storage/Data/Work/Postdoc/Work/testing/smeagol/cx1/Smeagol_Tutorial_Files/DFT_NEGF_Transport/Day_1/Inputs/testing/au_chain/run/ivan-sz-4s/output'
# folder_1 = '/Volumes/Storage/Data/Work/Postdoc/Work/testing/cp2k-smeagol/examples/SMEAGOL_TEST/Au_chain/transport-12-wfnrs-bulk-4-tidy-2/output'
# folder_2 = '/Volumes/Storage/Data/Work/Postdoc/Work/testing/smeagol/cx1/Smeagol_Tutorial_Files/DFT_NEGF_Transport/Day_1/Inputs/testing/li_chain/li-sz-explicit-2s-dos-12atoms-2/output'
# folder_2 = '/Volumes/Storage/Data/Work/Postdoc/Work/testing/cp2k-smeagol/examples/SMEAGOL_TEST/Au_chain/transport-12-wfnrs-li-bulk-4-tidy-3/output'
# folder_2 = '/Volumes/Storage/Data/Work/Postdoc/Work/testing/smeagol/cx1/Smeagol_Tutorial_Files/DFT_NEGF_Transport/Day_1/Inputs/testing/au_chain/run/clot-sz-4s/output'
# folder_2 = '/Volumes/Storage/Data/Work/Postdoc/Work/testing/smeagol/cx1/Smeagol_Tutorial_Files/DFT_NEGF_Transport/Day_1/Inputs/testing/au_chain/run/clot-sz-default/output'

# xlim = [-5, 10]  # Li and Au chain
# labels = ['vps', 'psf', 'leem']
# fermi_offset_1 = 0
# fermi_offset_2 = 0
# fermi_offset_3 = 0
# folder_1 = '/Volumes/Storage/Data/Work/Postdoc/Work/testing/smeagol/cx1/Smeagol_Tutorial_Files/DFT_NEGF_Transport/Day_1/Inputs/testing/au_chain-12-4atoms/run/clov-sz-4s-3d-cutoff/output'
# folder_2 = '/Volumes/Storage/Data/Work/Postdoc/Work/testing/smeagol/cx1/Smeagol_Tutorial_Files/DFT_NEGF_Transport/Day_1/Inputs/testing/au_chain-12-4atoms/run/clop-sz-4s-3d-cutoff/output'
# folder_3 = '/Volumes/Storage/Data/Work/Postdoc/Work/testing/smeagol/cx1/Smeagol_Tutorial_Files/DFT_NEGF_Transport/Day_1/Inputs/testing/au_chain-12-4atoms/run/leem-sz-4s-3d-cutoff/output'

# xlim = [-5, 10]  # Li and Au chain
# labels = ['CP2K', 'SIESTA (2)', 'SIESTA (8)']
# fermi_offset_1 = 0
# fermi_offset_2 = 0
# fermi_offset_3 = 0
# folder_1 = '/Volumes/Storage/Data/Work/Postdoc/Work/testing/cp2k-smeagol/examples/SMEAGOL_TEST/Au_chain/transport-12-AuAu-2p8/output'
# folder_2 = '/Volumes/Storage/Data/Work/Postdoc/Work/testing/smeagol/cx1/Smeagol_Tutorial_Files/DFT_NEGF_Transport/Day_1/Inputs/testing/au_chain-12-4atoms/run/Au-lda-leem-sz-4s-3d-cutoff-lda/output'
# folder_3 = '/Volumes/Storage/Data/Work/Postdoc/Work/testing/smeagol/cx1/Smeagol_Tutorial_Files/DFT_NEGF_Transport/Day_1/Inputs/testing/au_chain-12-4atoms/run/Au-lda-clop-sz-4s-3d-cutoff-lda/output'

# xlim = [-5, 10]  # Li and Au chain
# labels = ['CP2K', 'SIESTA (8) default Rc', 'SIESTA (8) longer Rc']
# fermi_offset_1 = 0
# fermi_offset_2 = 0
# fermi_offset_3 = 0
# folder_1 = '/Volumes/Elements/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/testing/cp2k-smeagol/examples/SMEAGOL_TEST/Au_chain/transport-12-AuAu-2p8/output'
# folder_2 = '/Volumes/Elements/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/testing/smeagol/cx1/Smeagol_Tutorial_Files/DFT_NEGF_Transport/Day_1/Inputs/testing/au_chain-12-4atoms/run/Au-lda-clop-sz-4s-3d-def-lda/output'
# folder_3 = '/Volumes/Elements/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/testing/smeagol/cx1/Smeagol_Tutorial_Files/DFT_NEGF_Transport/Day_1/Inputs/testing/au_chain-12-4atoms/run/Au-lda-clop-sz-4s-3d-cutoff-lda/output'

# xlim = [-5, 10]  # Li and Au chain
# ylim = [0.0, 1.0]
# ylim_log = [0.008, 1.2]
# labels = ['CP2K LDA q11', 'CP2K PBE q11', 'CP2K PBE q19']
# fermi_offset_1 = 0
# fermi_offset_2 = 0
# fermi_offset_3 = 0
# folder_1 = '/Volumes/Elements/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/testing/cp2k-smeagol/examples/SMEAGOL_TEST/Au_chain/transport-12-AuAu-2p8-lda/output'
# folder_2 = '/Volumes/Elements/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/testing/cp2k-smeagol/examples/SMEAGOL_TEST/Au_chain/transport-12-AuAu-2p8-pbe/output'
# folder_3 = '/Volumes/Elements/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/testing/cp2k-smeagol/examples/SMEAGOL_TEST/Au_chain/transport-12-AuAu-2p8-pbe-q19/output'

# xlim = [-3, 3]  # Au chain
# xlim = [-5, 10]  # Li chain
# ylim = [0.0, 1.0]
# ylim_log = [0.008, 1.2]
# labels = ['CP2K', 'SIESTA 6s, 5d', 'SIESTA 6s, 5d, 6p']
# fermi_offset_1 = 0
# fermi_offset_2 = 0
# fermi_offset_3 = 0
# folder_1 = '/Volumes/Elements/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/testing/cp2k-smeagol/examples/SMEAGOL_TEST/Au_chain/transport-12-AuAu-2p8-lda/output'
# folder_2 = '/Volumes/Elements/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/1D-chains/au-chain/scarf/siesta/au_chain-12-4atoms/run/Au-lda-clop-sz-6s-5d-cutoff-lda/output'
# folder_3 = '/Volumes/Elements/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/1D-chains/au-chain/scarf/siesta/au_chain-12-4atoms/run/Au-lda-clop-sz-6s-5d-6p-cutoff-lda/output'

# xlim = [-3, 3]  # Au chain
# labels = ['1x1x31, 1x1x1', '2x2x100, 2x2x1', '4x4x100, 4x4x1']
# fermi_offset_1 = 0
# fermi_offset_2 = 0
# fermi_offset_3 = 0
# ylim = [0.0, 1.0]
# ylim_log = [0.008, 1.2]
# folder_1 = '/Volumes/Elements/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/new/au-bdt/young/au-bdt/3x3/kpoints/bulk-1-1-31-em-1-1-1/output'
# folder_2 = '/Volumes/Elements/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/new/au-bdt/young/au-bdt/3x3/kpoints/bulk-2-2-100-em-2-2-1/output'
# folder_3 = '/Volumes/Elements/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/new/au-bdt/young/au-bdt/3x3/kpoints/bulk-4-4-100-em-4-4-1/output'

n = 1
file_1_1 = np.genfromtxt('{}/G0-S0.out'.format(folder_1), skip_header=2, skip_footer=3) * n
file_1_2 = np.genfromtxt('{}/G1-S0.out'.format(folder_1), skip_header=2, skip_footer=3) * n
file_1_3 = np.genfromtxt('{}/G1-S2.out'.format(folder_1), skip_header=2, skip_footer=3) * n
file_1_4 = np.genfromtxt('{}/G2-S0.out'.format(folder_1), skip_header=2, skip_footer=3) * n
file_1_5 = np.genfromtxt('{}/G3-S0.out'.format(folder_1), skip_header=2, skip_footer=1) * n
file_2_1 = np.genfromtxt('{}/G0-S0.out'.format(folder_2), skip_header=2, skip_footer=3) * n
file_2_2 = np.genfromtxt('{}/G1-S0.out'.format(folder_2), skip_header=2, skip_footer=3) * n
file_2_3 = np.genfromtxt('{}/G1-S2.out'.format(folder_2), skip_header=2, skip_footer=3) * n
file_2_4 = np.genfromtxt('{}/G2-S0.out'.format(folder_2), skip_header=2, skip_footer=3) * n
file_2_5 = np.genfromtxt('{}/G3-S0.out'.format(folder_2), skip_header=2, skip_footer=1) * n
file_3_1 = np.genfromtxt('{}/G0-S0.out'.format(folder_3), skip_header=2, skip_footer=3) * n
file_3_2 = np.genfromtxt('{}/G1-S0.out'.format(folder_3), skip_header=2, skip_footer=3) * n
file_3_3 = np.genfromtxt('{}/G1-S2.out'.format(folder_3), skip_header=2, skip_footer=3) * n
file_3_4 = np.genfromtxt('{}/G2-S0.out'.format(folder_3), skip_header=2, skip_footer=3) * n
file_3_5 = np.genfromtxt('{}/G3-S0.out'.format(folder_3), skip_header=2, skip_footer=1) * n

# Transmission
fig_plot_1, ax_plot_1 = plt.subplots()
ax_plot_1.plot(file_1_1[:, 0]+fermi_offset_1, file_1_1[:, 1], 'r-', label=labels[0])
ax_plot_1.plot(file_2_1[:, 0]+fermi_offset_2, file_2_1[:, 1], 'g-', label=labels[1])
ax_plot_1.plot(file_3_1[:, 0]+fermi_offset_3, file_3_1[:, 1], 'b-', label=labels[2])
ax_plot_1.set_xlim([xlim[0], xlim[1]])
ax_plot_1.legend(frameon=False)
ax_plot_1.set_xlabel(r'E-E$_{\mathrm{F}}$ (eV)')
ax_plot_1.set_ylabel('Transmission')
fig_plot_1.tight_layout()
fig_plot_1.savefig('{}/compare-transmission.png'.format(folder_1), dpi=param.save_dpi)
fig_plot_1.savefig('{}/compare-transmission.png'.format(folder_2), dpi=param.save_dpi)
fig_plot_1.savefig('{}/compare-transmission.png'.format(folder_3), dpi=param.save_dpi)

# Transmission log
fig_plot_1_2, ax_plot_1_2 = plt.subplots()
ax_plot_1_2.plot(file_1_1[:, 0]+fermi_offset_1, file_1_1[:, 1], 'r-', label=labels[0])
ax_plot_1_2.plot(file_2_1[:, 0]+fermi_offset_2, file_2_1[:, 1], 'g-', label=labels[1])
ax_plot_1_2.plot(file_3_1[:, 0]+fermi_offset_3, file_3_1[:, 1], 'b-', label=labels[2])
ax_plot_1_2.set_xlim([xlim[0], xlim[1]])
ax_plot_1_2.set_ylim([ylim_log[0], ylim_log[1]])
ax_plot_1_2.set_yscale('log')
ax_plot_1_2.set_xlabel(r'E-E$_{\mathrm{F}}$ (eV)')
ax_plot_1_2.set_ylabel('Transmission')
ax_plot_1_2.legend(frameon=False)
fig_plot_1_2.tight_layout()
fig_plot_1_2.savefig('{}/compare_transmission_log.png'.format(folder_1), dpi=param.save_dpi)
fig_plot_1_2.savefig('{}/compare_transmission_log.png'.format(folder_2), dpi=param.save_dpi)
fig_plot_1_2.savefig('{}/compare_transmission_log.png'.format(folder_3), dpi=param.save_dpi)

# Plot transmission and transmission log
rows, cols = 2, 1
fig1_3, ax1_3 = plt.subplots(rows, cols, sharex='col', sharey='row', figsize=(6, 8))
ax1_3[0].plot(file_1_1[:, 0]+fermi_offset_1, file_1_1[:, 1], 'r-', label=labels[0])
ax1_3[0].plot(file_2_1[:, 0]+fermi_offset_2, file_2_1[:, 1], 'g-', label=labels[1])
ax1_3[0].plot(file_3_1[:, 0]+fermi_offset_2, file_3_1[:, 1], 'b-', label=labels[2])
ax1_3[0].set_xlim([xlim[0], xlim[1]])
ax1_3[0].set_ylim([ylim[0], ylim[1]])
ax1_3[0].legend(frameon=False)
ax1_3[0].set_ylabel('Transmission')
ax1_3[1].plot(file_1_1[:, 0]+fermi_offset_1, file_1_1[:, 1], 'r-', label=labels[0])
ax1_3[1].plot(file_2_1[:, 0]+fermi_offset_2, file_2_1[:, 1], 'g-', label=labels[1])
ax1_3[1].plot(file_3_1[:, 0]+fermi_offset_2, file_3_1[:, 1], 'b-', label=labels[2])
ax1_3[1].set_xlim([xlim[0], xlim[1]])
ax1_3[1].set_ylim([ylim_log[0], ylim_log[1]])
ax1_3[1].legend(frameon=False)
ax1_3[1].set_yscale('log')
ax1_3[1].set_xlabel(r'E-E$_{\mathrm{F}}$ (eV)')
ax1_3[1].set_ylabel('Transmission')
fig1_3.tight_layout()
fig1_3.savefig('{}/compare_transmission_log_all.png'.format(folder_1), dpi=param.save_dpi)
fig1_3.savefig('{}/compare_transmission_log_all.png'.format(folder_2), dpi=param.save_dpi)
fig1_3.savefig('{}/compare_transmission_log_all.png'.format(folder_3), dpi=param.save_dpi)

# EM DOS
fig_plot_4, ax_plot_4 = plt.subplots()
ax_plot_4.plot(file_1_4[:, 0]+fermi_offset_1, file_1_4[:, 1], 'r-', label=labels[0])
ax_plot_4.plot(file_2_4[:, 0]+fermi_offset_2, file_2_4[:, 1], 'g-', label=labels[1])
ax_plot_4.plot(file_3_4[:, 0]+fermi_offset_3, file_3_4[:, 1], 'b-', label=labels[2])
ax_plot_4.set_xlim([xlim[0], xlim[1]])
ax_plot_4.legend(frameon=False)
ax_plot_4.set_xlabel(r'E-E$_{\mathrm{F}}$ (eV)')
ax_plot_4.set_ylabel('EM DOS')
fig_plot_4.tight_layout()
fig_plot_4.savefig('{}/compare-dos_em.png'.format(folder_1), dpi=param.save_dpi)
fig_plot_4.savefig('{}/compare-dos_em.png'.format(folder_2), dpi=param.save_dpi)
fig_plot_4.savefig('{}/compare-dos_em.png'.format(folder_3), dpi=param.save_dpi)

# Leads DOS
fig_plot_5, ax_plot_5 = plt.subplots()
ax_plot_5.plot(file_1_5[:, 0]+fermi_offset_1, file_1_5[:, 1], 'r-', label=labels[0])
ax_plot_5.plot(file_2_5[:, 0]+fermi_offset_2, file_2_5[:, 1], 'g-', label=labels[2])
ax_plot_5.plot(file_3_5[:, 0]+fermi_offset_3, file_3_5[:, 1], 'b-', label=labels[1])
ax_plot_5.set_xlim([xlim[0], xlim[1]])
ax_plot_5.legend(frameon=False)
ax_plot_5.set_xlabel(r'E-E$_{\mathrm{F}}$ (eV)')
ax_plot_5.set_ylabel('Leads DOS')
fig_plot_5.tight_layout()
fig_plot_5.savefig('{}/compare-dos_leads.png'.format(folder_1), dpi=param.save_dpi)
fig_plot_5.savefig('{}/compare-dos_leads.png'.format(folder_2), dpi=param.save_dpi)
fig_plot_5.savefig('{}/compare-dos_leads.png'.format(folder_3), dpi=param.save_dpi)

# Plot transmission and EM DOS
rows, cols = 2, 1
fig_plot_all, ax_plot_all = plt.subplots(rows, cols, sharex='col', sharey='row', figsize=(5, 8))
ax_plot_all[0].plot(file_1_1[:, 0]+fermi_offset_1, file_1_1[:, 1], 'r-', label=labels[0])
ax_plot_all[0].plot(file_2_1[:, 0]+fermi_offset_2, file_2_1[:, 1], 'g-', label=labels[1])
ax_plot_all[0].plot(file_3_1[:, 0]+fermi_offset_3, file_3_1[:, 1], 'b-', label=labels[2])
ax_plot_all[0].set_xlim([xlim[0], xlim[1]])
ax_plot_all[0].set_ylabel('Transmission')
ax_plot_all[0].legend(frameon=False)
ax_plot_all[1].plot(file_1_4[:, 0]+fermi_offset_1, file_1_4[:, 1], 'r-', label=labels[0])
ax_plot_all[1].plot(file_2_4[:, 0]+fermi_offset_2, file_2_4[:, 1], 'g-', label=labels[1])
ax_plot_all[1].plot(file_3_4[:, 0]+fermi_offset_3, file_3_4[:, 1], 'b-', label=labels[2])
ax_plot_all[1].set_xlim([xlim[0], xlim[1]])
ax_plot_all[1].legend(frameon=False)
ax_plot_all[1].set_xlabel(r'E-E$_{\mathrm{F}}$ (eV)')
ax_plot_all[1].set_ylabel('Density of states (a.u.)')
fig_plot_all.tight_layout()
# fig_plot_all.subplots_adjust(hspace=0)
fig_plot_all.savefig('{}/compare-transmission_emdos.png'.format(folder_1), dpi=param.save_dpi)
fig_plot_all.savefig('{}/compare-transmission_emdos.png'.format(folder_2), dpi=param.save_dpi)
fig_plot_all.savefig('{}/compare-transmission_emdos.png'.format(folder_3), dpi=param.save_dpi)

# Plot all
rows, cols = 2, 2
fig_plot_all, ax_plot_all = plt.subplots(rows, cols, figsize=(10, 8))
# fig_plot_all, ax_plot_all = plt.subplots(rows, cols, sharex='col', sharey='row', figsize=(10, 8))
ax_plot_all[0, 0].plot(file_1_1[:, 0]+fermi_offset_1, file_1_1[:, 1], 'r-', label=labels[0])
ax_plot_all[0, 0].plot(file_2_1[:, 0]+fermi_offset_2, file_2_1[:, 1], 'g-', label=labels[1])
ax_plot_all[0, 0].plot(file_3_1[:, 0]+fermi_offset_3, file_3_1[:, 1], 'b-', label=labels[2])
ax_plot_all[0, 0].set_xlim([xlim[0], xlim[1]])
ax_plot_all[0, 0].set_ylim([ylim[0], ylim[1]])
ax_plot_all[0, 0].legend(frameon=False)
# ax_plot_all[0, 0].set_xlabel(r'E-E$_{\mathrm{F}}$ (eV)')
ax_plot_all[0, 0].set_ylabel('Transmission')
ax_plot_all[0, 1].plot(file_1_2[:, 0]+fermi_offset_1, file_1_2[:, 1], 'r-', label=labels[0])
ax_plot_all[0, 1].plot(file_2_2[:, 0]+fermi_offset_2, file_2_2[:, 1], 'g-', label=labels[1])
ax_plot_all[0, 1].plot(file_3_2[:, 0]+fermi_offset_3, file_3_2[:, 1], 'b-', label=labels[2])
ax_plot_all[0, 1].set_xlim([xlim[0], xlim[1]])
ax_plot_all[0, 1].legend(frameon=False)
# ax_plot_all[0, 1].set_xlabel(r'E-E$_{\mathrm{F}}$ (eV)')
ax_plot_all[0, 1].set_ylabel('Number of channels')
ax_plot_all[1, 0].plot(file_1_4[:, 0]+fermi_offset_1, file_1_4[:, 1], 'r-', label=labels[0])
ax_plot_all[1, 0].plot(file_2_4[:, 0]+fermi_offset_2, file_2_4[:, 1], 'g-', label=labels[1])
ax_plot_all[1, 0].plot(file_3_4[:, 0]+fermi_offset_3, file_3_4[:, 1], 'b-', label=labels[2])
ax_plot_all[1, 0].set_xlim([xlim[0], xlim[1]])
ax_plot_all[1, 0].legend(frameon=False)
ax_plot_all[1, 0].set_xlabel(r'E-E$_{\mathrm{F}}$ (eV)')
ax_plot_all[1, 0].set_ylabel('EM DOS')
ax_plot_all[1, 1].plot(file_1_5[:, 0]+fermi_offset_1, file_1_5[:, 1], 'r-', label=labels[0])
ax_plot_all[1, 1].plot(file_2_5[:, 0]+fermi_offset_2, file_2_5[:, 1], 'g-', label=labels[1])
ax_plot_all[1, 1].plot(file_3_5[:, 0]+fermi_offset_3, file_3_5[:, 1], 'b-', label=labels[2])
ax_plot_all[1, 1].set_xlim([xlim[0], xlim[1]])
ax_plot_all[1, 1].legend(frameon=False)
ax_plot_all[1, 1].set_xlabel(r'E-E$_{\mathrm{F}}$ (eV)')
ax_plot_all[1, 1].set_ylabel('Leads DOS')
fig_plot_all.tight_layout()
# fig_plot_all.subplots_adjust(hspace=0)
fig_plot_all.savefig('{}/compare-plot_all.png'.format(folder_1), dpi=param.save_dpi)
fig_plot_all.savefig('{}/compare-plot_all.png'.format(folder_2), dpi=param.save_dpi)
fig_plot_all.savefig('{}/compare-plot_all.png'.format(folder_3), dpi=param.save_dpi)

if __name__ == "__main__":
    print('Finished.')
    plt.show()
