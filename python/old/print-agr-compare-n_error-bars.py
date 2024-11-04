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
# plot_lengend = True
# labels = ['CP2K-SMEAGOL', 'SIESTA-SMEAGOL']
# fermi = np.zeros(len(labels))
# folder = ['/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-chain/cp2k-smeagol/transmission/V-0_HLB-F_z-0-0_atoms-28_noprint/output',
#           '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-chain/siesta/transmission/V-0_HLB-0_z-0-0/output']
# labels_cp2k_negf = ['CP2K-NEGF']
# fermi_cp2k_negf = [-0.22124430588176]
# dos_norm_cp2k_negf = 28
# folder_cp2k_negf = []
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
# xlim = [-4, 4]
# ylim = [0.0, 1.0]
# ylim_log = [0.008, 1.2]
# ylim_dos = [0, 300]
# labels = ['CP2K', 'SIESTA']
# fermi = np.zeros(len(labels))
# # fermi = [0, 0.0, 0.7]
# folder = ['/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/cp2k-smeagol/sz/transmission/exp/dzvp/sergey/kpoints_bulk-4-4-100_em-4-4-1_hlb-t-11.03197_scf-500/output',
#           '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/siesta-smeagol/transmission/exp/bulk-4-4-100-em-4-4-1_hlb-15.2496/output']
# folder_cp2k_negf = []
# fermi_cp2k_negf = 0
# plot_lengend = True

# Plot Au-BDT for experimental CP2K and SIESTA
# xlim = [-4, 4]
# ylim = [0.0, 1.02]
# ylim_log = [0.008, 1.2]
# ylim_dos = [0, 300]
# labels = ['CP2K+SMEAGOL', 'SIESTA+SMEAGOL']
# # labels = ['2x2x20', '4x4x20', '4x4x100']
# fermi = np.zeros(len(labels))
# # fermi = [0, 0.0, 0.7]
# folder = [
#     '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/au-bdt/cp2k-smeagol/exp/transmission/kpoints-4-4-20/output',
#     '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/au-bdt/siesta1-smeagol/transmission/exp_kpoints-4-4-20/output',
#     ]
# # folder = [
# #     '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/au-bdt/cp2k-smeagol/exp/transmission/kpoints-1-1-20/output',
# #     '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/au-bdt/cp2k-smeagol/exp/transmission/kpoints-4-4-20/output',
# #     '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/au-bdt/cp2k-smeagol/exp/transmission/kpoints-4-4-100/output']
# # folder = ['/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/au-bdt/siesta1-smeagol/transmission/exp_kpoints-2-2-20/output',
# #           '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/au-bdt/siesta1-smeagol/transmission/exp_kpoints-4-4-20/output',
# #           '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/au-bdt/siesta1-smeagol/transmission/exp_kpoints-4-4-100/output']
# folder_cp2k_negf = []
# fermi_cp2k_negf = 0
# plot_lengend = True

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

# Plot Au wire solvated
# xlim = [-2, 3]
xlim = [-5, 4]
# ylim = [0.0, 2.2]
ylim = [-0.1, 3]
ylim_log = [0.008, 1.2]
ylim_dos = [0, 300]
# labels = ['0 fs', '50 fs']
# labels = ['0 fs no water', '50 fs no water']
# labels = ['0 fs dummy water', '50 fs dummy water']
# labels = ['0 fs dummy water', '50 fs dummy water', '100 fs dummy water', '150 fs dummy water', '200 fs dummy water']
# labels = ['Solvated Au wire', 'Solvated Au wire (water removed)', 'Ideal Au wire']
# labels = ['Au wire (solvated)', 'Au wire (vacuum)', 'Water', 'Au wire (ideal)']
labels = ['Au wire (solvated)', 'Au wire (vacuum)', 'Au wire (ideal)']
# labels = ['Ideal Au wire DZVP 11e', 'Ideal Au wire SZ 11e', 'Ideal Au wire SZ 1e']
fermi = np.zeros(len(labels))
# fermi = [0, 0.0, 0.7]
# folder = [
#     '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-wire-ismael/jad/calculations/single-points/snapshots/restart-auto-long_au-frozen/md_dft_long_LINEAR_P/jobs/step-2001/V-0/output',
#     '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-wire-ismael/jad/calculations/single-points/snapshots/restart-auto-long_au-frozen/md_dft_long_LINEAR_P/jobs/step-2101/V-0/output',
# ]
# folder = [
#     '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-wire-ismael/jad/calculations/single-points/snapshots/restart-auto-long_au-frozen/md_dft_long_LINEAR_P_nowater/jobs/step-2001/V-0/output',
#     '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-wire-ismael/jad/calculations/single-points/snapshots/restart-auto-long_au-frozen/md_dft_long_LINEAR_P_nowater/jobs/step-2101/V-0/output',
# ]
# folder = [
#     '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-wire-ismael/jad/calculations/single-points/snapshots/restart-auto-long_au-frozen/md_dft_long_LINEAR_P_dummy/jobs/step-2001/V-0/output',
#     '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-wire-ismael/jad/calculations/single-points/snapshots/restart-auto-long_au-frozen/md_dft_long_LINEAR_P_dummy/jobs/step-2101/V-0/output',
# ]
# folder = [
#     '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-wire-ismael/jad/calculations/single-points/snapshots/restart-auto-long_au-frozen/md_dft_long_LINEAR_P_dummy/jobs/step-2001/V-0/output',
#     '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-wire-ismael/jad/calculations/single-points/snapshots/restart-auto-long_au-frozen/md_dft_long_LINEAR_P_dummy/jobs/step-2101/V-0/output',
#     '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-wire-ismael/jad/calculations/single-points/snapshots/restart-auto-long_au-frozen/md_dft_long_LINEAR_P_dummy/jobs/step-2201/V-0/output',
#     '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-wire-ismael/jad/calculations/single-points/snapshots/restart-auto-long_au-frozen/md_dft_long_LINEAR_P_dummy/jobs/step-2301/V-0/output',
#     '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-wire-ismael/jad/calculations/single-points/snapshots/restart-auto-long_au-frozen/md_dft_long_LINEAR_P_dummy/jobs/step-2401/V-0/output',
# ]
# folder = [
#     '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-wire-ismael/jad/calculations/single-points/snapshots/restart-auto-long_au-frozen/md_dft_long_LINEAR_P/jobs/step-2001/V-0/output',
#     '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-wire-ismael/jad/calculations/single-points/snapshots/restart-auto-long_au-frozen/md_dft_long_LINEAR_P_dummy/jobs/step-2001/V-0/output',
#     '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-chain/cp2k-smeagol/transmission/delete/V-0_HLB-F_z-0-0_atoms-28_noprint/output',
# ]
# folder = [
#     '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-chain/cp2k-smeagol/transmission/au-dzvp-q11_143b/output',
#     '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-chain/cp2k-smeagol/transmission/au-sz-q11_143b/output',
#     '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-chain/cp2k-smeagol/transmission/au-sz-q1_143b/output',
# ]
folder = [
    '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-wire-ismael/jad/calculations/single-points/snapshots/restart-auto-long_au-frozen/md_dft_long_LINEAR_P/jobs/step-2001/V-0/output',
    '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-wire-ismael/jad/calculations/single-points/snapshots/restart-auto-long_au-frozen/md_dft_long_LINEAR_P/jobs/step-2101/V-0/output',
    '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-wire-ismael/jad/calculations/single-points/snapshots/restart-auto-long_au-frozen/md_dft_long_LINEAR_P/jobs/step-2201/V-0/output',
    '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-wire-ismael/jad/calculations/single-points/snapshots/restart-auto-long_au-frozen/md_dft_long_LINEAR_P/jobs/step-2301/V-0/output',
    '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-wire-ismael/jad/calculations/single-points/snapshots/restart-auto-long_au-frozen/md_dft_long_LINEAR_P/jobs/step-2401/V-0/output',
    '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-wire-ismael/jad/calculations/single-points/snapshots/restart-auto-long_au-frozen/md_dft_long_LINEAR_P_dummy/jobs/step-2001/V-0/output',
    '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-wire-ismael/jad/calculations/single-points/snapshots/restart-auto-long_au-frozen/md_dft_long_LINEAR_P_dummy/jobs/step-2101/V-0/output',
    '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-wire-ismael/jad/calculations/single-points/snapshots/restart-auto-long_au-frozen/md_dft_long_LINEAR_P_dummy/jobs/step-2201/V-0/output',
    '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-wire-ismael/jad/calculations/single-points/snapshots/restart-auto-long_au-frozen/md_dft_long_LINEAR_P_dummy/jobs/step-2301/V-0/output',
    '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-wire-ismael/jad/calculations/single-points/snapshots/restart-auto-long_au-frozen/md_dft_long_LINEAR_P_dummy/jobs/step-2401/V-0/output',
    '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-wire-ismael/jad/calculations/single-points/snapshots/restart-auto-long_au-frozen/md_dft_long_LINEAR_P_nowire/jobs/step-2001/V-0/output',
    '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-chain/cp2k-smeagol/transmission/delete/V-0_HLB-F_z-0-0_atoms-28_noprint/output',
]
# plotting_colors = ['r', 'g', 'b', 'k']
plotting_colors = ['r', 'g', 'k']
folder_cp2k_negf = []
fermi_cp2k_negf = 0
plot_lengend = True

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

data_1 = np.array([file_1[0][:, 1], file_1[1][:, 1], file_1[2][:, 1], file_1[3][:, 1], file_1[4][:, 1]])
data_2 = np.array([file_1[5][:, 1], file_1[6][:, 1], file_1[7][:, 1], file_1[8][:, 1], file_1[9][:, 1]])

# Transmission and error bards
fig_plot_1_error, ax_plot_1_error = plt.subplots()

# Mean and standard deviation 1
data_1_average = np.average(data_1, axis=0)
data_1_std = np.std(data_1, axis=0)
ax_plot_1_error.fill_between(file_1[0][:, 0], data_1_average - data_1_std, data_1_average + data_1_std, alpha=0.2, color=plotting_colors[0])
ax_plot_1_error.plot(file_1[0][:, 0], data_1_average, color=plotting_colors[0], label=labels[0])

# Mean and standard deviation 2
data_2_average = np.average(data_2, axis=0)
data_2_std = np.std(data_2, axis=0)
ax_plot_1_error.fill_between(file_1[0][:, 0], data_2_average - data_2_std, data_2_average + data_2_std, alpha=0.2, color=plotting_colors[1])
ax_plot_1_error.plot(file_1[0][:, 0], data_2_average, color=plotting_colors[1], label=labels[1])

# No wire
# ax_plot_1_error.plot(file_1[-2][:, 0], file_1[-2][:, 1], color=plotting_colors[-2], label=labels[-2])

# Model wire
ax_plot_1_error.plot(file_1[-1][:, 0], file_1[-1][:, 1], color=plotting_colors[-1], label=labels[-1])

ax_plot_1_error.set_xlim([xlim[0], xlim[1]])
ax_plot_1_error.set_ylim([ylim[0], ylim[1]])
# ax_plot_1_error.legend(frameon=False)
ax_plot_1_error.set_xlabel(r'E-E$_{\mathrm{F}}$ (eV)')
ax_plot_1_error.set_ylabel('Transmission')
fig_plot_1_error.tight_layout()
for i in range(len(folder_cp2k_negf)):
    fig_plot_1_error.savefig('{}/compare-transmission_no_label.png'.format(folder_cp2k_negf[i]), dpi=param.save_dpi)
for i in range(len(folder)):
    fig_plot_1_error.savefig('{}/compare-transmission_no_label.png'.format(folder[i]), dpi=param.save_dpi)

if __name__ == "__main__":
    print('Finished.')
    plt.show()
