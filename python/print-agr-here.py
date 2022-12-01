import numpy as np
import matplotlib.pyplot as plt
import os

""" Plotting of SMEAGOL output _TRC.agr from command line"""

# 1D Li chain
xlim = [-5, 10]
ylim = [-0.1, 1.1]
ylim_log = [0.008, 1.0]
fermi_offset_1 = 0

# Au-BDT
# xlim = [-3, 3]
# ylim = [0.0, 1.2]
# ylim_log = [0.008, 1.2]
# fermi_offset_1 = 1.347003

# Au-BDT
# xlim = [-3, 3]
# ylim = [0.0, 1.0]
# ylim_log = [0.008, 1.0]
# fermi_offset_1 = 0

save_dpi = 500  # DPI to save figure
params = {'axes.formatter.limits': [-4, 4],
          'axes.labelsize': 'large',
          'axes.titlesize': 'large',
          'lines.markersize': '8',
          }
plt.rcParams.update(params)

# Get files
folder_1 = os.getcwd()
file_1_1 = np.genfromtxt('{}/G0-S0.out'.format(folder_1), skip_header=2, skip_footer=3)
file_1_2 = np.genfromtxt('{}/G1-S0.out'.format(folder_1), skip_header=2, skip_footer=3)
file_1_3 = np.genfromtxt('{}/G1-S2.out'.format(folder_1), skip_header=2, skip_footer=3)
file_1_4 = np.genfromtxt('{}/G2-S0.out'.format(folder_1), skip_header=2, skip_footer=3)
file_1_5 = np.genfromtxt('{}/G3-S0.out'.format(folder_1), skip_header=2, skip_footer=1)

# Transmission
fig_plot_1, ax_plot_1 = plt.subplots()
ax_plot_1.plot(file_1_1[:, 0]+fermi_offset_1, file_1_1[:, 1], 'k-')
ax_plot_1.set_xlim([xlim[0], xlim[1]])
ax_plot_1.set_ylim([ylim[0], ylim[1]])
ax_plot_1.set_xlabel(r'E-E$_{\mathrm{F}}$ (eV)')
ax_plot_1.set_ylabel('Transmission')
fig_plot_1.tight_layout()
fig_plot_1.savefig('{}/transmission.png'.format(folder_1), dpi=save_dpi)

# Transmission log
fig_plot_1_2, ax_plot_1_2 = plt.subplots()
ax_plot_1_2.plot(file_1_1[:, 0]+fermi_offset_1, file_1_1[:, 1], 'k-')
ax_plot_1_2.set_xlim([xlim[0], xlim[1]])
ax_plot_1_2.set_ylim([ylim_log[0], ylim_log[1]])
ax_plot_1_2.set_yscale('log')
ax_plot_1_2.set_xlabel(r'E-E$_{\mathrm{F}}$ (eV)')
ax_plot_1_2.set_ylabel('Transmission')
fig_plot_1_2.tight_layout()
fig_plot_1_2.savefig('{}/transmission_log.png'.format(folder_1), dpi=save_dpi)

# Plot transmission and transmission log
rows, cols = 2, 1
fig1_3, ax1_3 = plt.subplots(rows, cols, sharex='col', sharey='row', figsize=(6, 8))
ax1_3[0].plot(file_1_1[:, 0]+fermi_offset_1, file_1_1[:, 1], 'k-')
ax1_3[0].set_xlim([xlim[0], xlim[1]])
ax1_3[0].set_ylim([ylim[0], ylim[1]])
ax1_3[0].set_ylabel('Transmission')
ax1_3[1].plot(file_1_1[:, 0]+fermi_offset_1, file_1_1[:, 1], 'k-')
ax1_3[1].set_xlim([xlim[0], xlim[1]])
ax1_3[1].set_ylim([ylim_log[0], ylim_log[1]])
ax1_3[1].set_yscale('log')
ax1_3[1].set_xlabel(r'E-E$_{\mathrm{F}}$ (eV)')
ax1_3[1].set_ylabel('Transmission')
fig1_3.tight_layout()
fig1_3.savefig('{}/transmission_log_all.png'.format(folder_1), dpi=save_dpi)

# Number of channels 1
fig_plot_2, ax_plot_2 = plt.subplots()
ax_plot_2.plot(file_1_2[:, 0]+fermi_offset_1, file_1_2[:, 1], 'k-')
ax_plot_2.set_xlim([xlim[0], xlim[1]])
# ax_plot_2.set_ylim([ylim[0], ylim[1]])
ax_plot_2.set_xlabel(r'E-E$_{\mathrm{F}}$ (eV)')
ax_plot_2.set_ylabel('Number of channels')
fig_plot_2.tight_layout()
fig_plot_2.savefig('{}/channels_1.png'.format(folder_1), dpi=save_dpi)

# Number of channels 2
# fig_plot_3, ax_plot_3 = plt.subplots()
# ax_plot_3.plot(file_1_3[:, 0], file_1_3[:, 1], 'k-')
# ax_plot_3.set_xlim([xlim[0], xlim[1]])
# ax_plot_3.set_xlabel(r'E-E$_{\mathrm{F}}$ (eV)')
# ax_plot_3.set_ylabel('Number of channels')
# fig_plot_3.tight_layout()
# fig_plot_3.savefig('{}/channels_2.png'.format(folder_1), dpi=save_dpi)

# EM DOS
fig_plot_4, ax_plot_4 = plt.subplots()
ax_plot_4.plot(file_1_4[:, 0]+fermi_offset_1, file_1_4[:, 1], 'k-')
ax_plot_4.set_xlim([xlim[0], xlim[1]])
ax_plot_4.set_xlabel(r'E-E$_{\mathrm{F}}$ (eV)')
ax_plot_4.set_ylabel('EM DOS')
fig_plot_4.tight_layout()
fig_plot_4.savefig('{}/dos_em.png'.format(folder_1), dpi=save_dpi)

# Leads DOS
fig_plot_5, ax_plot_5 = plt.subplots()
ax_plot_5.plot(file_1_5[:, 0]+fermi_offset_1, file_1_5[:, 1], 'k-')
ax_plot_5.set_xlim([xlim[0], xlim[1]])
ax_plot_5.set_xlabel(r'E-E$_{\mathrm{F}}$ (eV)')
ax_plot_5.set_ylabel('Leads DOS')
fig_plot_5.tight_layout()
fig_plot_5.savefig('{}/dos_leads.png'.format(folder_1), dpi=save_dpi)

# Plot transmission and EM DOS
rows, cols = 2, 1
fig_all_trans_dos, ax_all_trans_dos = plt.subplots(rows, cols, sharex='col', sharey='row', figsize=(5, 8))
ax_all_trans_dos[0].plot(file_1_1[:, 0]+fermi_offset_1, file_1_1[:, 1], 'k-')
ax_all_trans_dos[0].set_xlim([xlim[0], xlim[1]])
ax_all_trans_dos[0].set_ylabel('Transmission')
ax_all_trans_dos[1].plot(file_1_4[:, 0], file_1_4[:, 1], 'k-')
ax_all_trans_dos[1].set_xlim([xlim[0], xlim[1]])
ax_all_trans_dos[1].set_xlabel(r'E-E$_{\mathrm{F}}$ (eV)')
ax_all_trans_dos[1].set_ylabel('Density of states (a.u.)')
fig_all_trans_dos.tight_layout()
fig_all_trans_dos.savefig('{}/transmission_emdos.png'.format(folder_1), dpi=save_dpi)

# Plot all
rows, cols = 2, 2
fig_plot_all, ax_plot_all = plt.subplots(rows, cols, sharex='col', sharey='row', figsize=(10, 8))
ax_plot_all[0, 0].plot(file_1_1[:, 0], file_1_1[:, 1], 'k-')
ax_plot_all[0, 0].set_xlim([xlim[0], xlim[1]])
ax_plot_all[0, 0].set_ylim([ylim[0], ylim[1]])
# ax_plot_all[0, 0].set_xlabel(r'E-E$_{\mathrm{F}}$ (eV)')
ax_plot_all[0, 0].set_ylabel('Transmission')
ax_plot_all[0, 1].plot(file_1_2[:, 0], file_1_2[:, 1], 'k-')
ax_plot_all[0, 1].set_xlim([xlim[0], xlim[1]])
# ax_plot_all[0, 1].set_ylim([ylim[0], ylim[1]])
# ax_plot_all[0, 1].set_xlabel(r'E-E$_{\mathrm{F}}$ (eV)')
ax_plot_all[0, 1].set_ylabel('Number of channels')
ax_plot_all[1, 0].plot(file_1_4[:, 0], file_1_4[:, 1], 'k-')
ax_plot_all[1, 0].set_xlim([xlim[0], xlim[1]])
ax_plot_all[1, 0].set_xlabel(r'E-E$_{\mathrm{F}}$ (eV)')
ax_plot_all[1, 0].set_ylabel('EM DOS')
ax_plot_all[1, 1].plot(file_1_5[:, 0], file_1_5[:, 1], 'k-')
ax_plot_all[1, 1].set_xlim([xlim[0], xlim[1]])
ax_plot_all[1, 1].set_xlabel(r'E-E$_{\mathrm{F}}$ (eV)')
ax_plot_all[1, 1].set_ylabel('Leads DOS')
fig_plot_all.tight_layout()
# fig_plot_all.subplots_adjust(hspace=0)
fig_plot_all.savefig('{}/plot_all.png'.format(folder_1), dpi=save_dpi)
