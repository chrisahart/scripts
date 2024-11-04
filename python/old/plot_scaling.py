import numpy as np
import matplotlib.pyplot as plt
from general import parameters as param

""" Scaling """

folder_1 = '/Volumes/Storage/Data/Work/Postdoc/Work/calculations/testing/cp2k-smeagol/bdt/siesta/hollow-site/7-6/transport-siesta/bulk-4-4/new/scaling-pp8/output'

# MPI
cx1_cores = np.array([16, 32, 64, 96, 128])
cx1_scf_step0_1_mpi = np.genfromtxt('{}/scf-step0_1_mpi.out'.format(folder_1), skip_header=0, skip_footer=0)
cx1_scf_step0_2_mpi = np.genfromtxt('{}/scf-step0_2_mpi.out'.format(folder_1), skip_header=0, skip_footer=0)
cx1_scf_step0_3_mpi = np.genfromtxt('{}/scf-step0_3_mpi.out'.format(folder_1), skip_header=0, skip_footer=0)
cx1_scf_1_mpi = np.genfromtxt('{}/scf_1_mpi.out'.format(folder_1), skip_header=0, skip_footer=0)
cx1_scf_2_mpi = np.genfromtxt('{}/scf_2_mpi.out'.format(folder_1), skip_header=0, skip_footer=0)
cx1_scf_3_mpi = np.genfromtxt('{}/scf_3_mpi.out'.format(folder_1), skip_header=0, skip_footer=0)
cx1_siesta_1_mpi = np.genfromtxt('{}/siesta_1_mpi.out'.format(folder_1), skip_header=0, skip_footer=0)
cx1_siesta_2_mpi = np.genfromtxt('{}/siesta_2_mpi.out'.format(folder_1), skip_header=0, skip_footer=0)
cx1_siesta_3_mpi = np.genfromtxt('{}/siesta_3_mpi.out'.format(folder_1), skip_header=0, skip_footer=0)
cx1_transm_3_mpi = np.genfromtxt('{}/transm_3_mpi.out'.format(folder_1), skip_header=0, skip_footer=0)

# OMP
cx1_threads_omp = np.array([1, 2, 4, 8])
cx1_mpi_omp = np.array([128, 64, 32, 16])
cx1_scf_step0_1_omp = np.genfromtxt('{}/scf-step0_1_omp.out'.format(folder_1), skip_header=0, skip_footer=0)
cx1_scf_step0_2_omp = np.genfromtxt('{}/scf-step0_2_omp.out'.format(folder_1), skip_header=0, skip_footer=0)
cx1_scf_step0_3_omp = np.genfromtxt('{}/scf-step0_3_omp.out'.format(folder_1), skip_header=0, skip_footer=0)
cx1_scf_1_omp = np.genfromtxt('{}/scf_1_omp.out'.format(folder_1), skip_header=0, skip_footer=0)
cx1_scf_2_omp = np.genfromtxt('{}/scf_2_omp.out'.format(folder_1), skip_header=0, skip_footer=0)
cx1_scf_3_omp = np.genfromtxt('{}/scf_3_omp.out'.format(folder_1), skip_header=0, skip_footer=0)
cx1_siesta_1_omp = np.genfromtxt('{}/siesta_1_omp.out'.format(folder_1), skip_header=0, skip_footer=0)
cx1_siesta_2_omp = np.genfromtxt('{}/siesta_2_omp.out'.format(folder_1), skip_header=0, skip_footer=0)
cx1_siesta_3_omp = np.genfromtxt('{}/siesta_3_omp.out'.format(folder_1), skip_header=0, skip_footer=0)
cx1_transm_3_omp = np.genfromtxt('{}/transm_3_omp.out'.format(folder_1), skip_header=0, skip_footer=0)

folder_2 = '/Volumes/Storage/Data/Work/Postdoc/Work/calculations/transport/au-bdt/scarf/new/scaling-pp8/output'

# MPI
scarf_cores = np.array([16, 24, 48, 72, 96, 120])
scarf_scf_step0_1_mpi = np.genfromtxt('{}/scf-step0_1_mpi.out'.format(folder_2), skip_header=0, skip_footer=0)
scarf_scf_step0_2_mpi = np.genfromtxt('{}/scf-step0_2_mpi.out'.format(folder_2), skip_header=0, skip_footer=0)
scarf_scf_step0_3_mpi = np.genfromtxt('{}/scf-step0_3_mpi.out'.format(folder_2), skip_header=0, skip_footer=0)
scarf_scf_1_mpi = np.genfromtxt('{}/scf_1_mpi.out'.format(folder_2), skip_header=0, skip_footer=0)
scarf_scf_2_mpi = np.genfromtxt('{}/scf_2_mpi.out'.format(folder_2), skip_header=0, skip_footer=0)
scarf_scf_3_mpi = np.genfromtxt('{}/scf_3_mpi.out'.format(folder_2), skip_header=0, skip_footer=0)
scarf_siesta_1_mpi = np.genfromtxt('{}/siesta_1_mpi.out'.format(folder_2), skip_header=0, skip_footer=0)
scarf_siesta_2_mpi = np.genfromtxt('{}/siesta_2_mpi.out'.format(folder_2), skip_header=0, skip_footer=0)
scarf_siesta_3_mpi = np.genfromtxt('{}/siesta_3_mpi.out'.format(folder_2), skip_header=0, skip_footer=0)
scarf_transm_3_mpi = np.genfromtxt('{}/transm_3_mpi.out'.format(folder_2), skip_header=0, skip_footer=0)


# CX1 MPI SCF-step0 all 
fig_plot_1, ax_plot_1 = plt.subplots()
# ax_plot_1.plot(cx1_cores, cx1_scf_step0_1_mpi/cx1_cores/60, 'r.-', label='CX1 BulkLeads')
# ax_plot_1.plot(cx1_cores, cx1_scf_step0_2_mpi/cx1_cores/60, 'g.-', label='CX1 EM DFT')
# ax_plot_1.plot(scarf_cores, scarf_scf_step0_3_mpi/scarf_cores/60, 'b.-', label='CX1 EMTransport')
ax_plot_1.plot(scarf_cores, scarf_scf_step0_1_mpi/scarf_cores/60, 'r.-', label='SCARF BulkLeads')
ax_plot_1.plot(scarf_cores, scarf_scf_step0_2_mpi/scarf_cores/60, 'g.-', label='SCARF EM DFT')
ax_plot_1.plot(scarf_cores, scarf_scf_step0_3_mpi/scarf_cores/60, 'b.-', label='SCARF EMTransport')
ax_plot_1.legend(frameon=False)
ax_plot_1.set_xlabel('Number of cores (threads=1)')
ax_plot_1.set_ylabel('Time for first SCF step / minutes')
fig_plot_1.tight_layout()
# fig_plot_1.savefig('{}/cx1_step0_all_mpi.png'.format(folder_1), dpi=param.save_dpi)
fig_plot_1.savefig('{}/scarf_step0_all_mpi.png'.format(folder_2), dpi=param.save_dpi)

# # CX1 MPI SCF-step0 all
# fig_plot_2, ax_plot_2 = plt.subplots()
# ax_plot_2.plot(cx1_threads_omp, cx1_scf_step0_1_omp/cx1_mpi_omp/60, 'r.-', label='BulkLeads')
# ax_plot_2.plot(cx1_threads_omp, cx1_scf_step0_2_omp/cx1_mpi_omp/60, 'g.-', label='EM DFT')
# ax_plot_2.plot(cx1_threads_omp, cx1_scf_step0_3_omp/cx1_mpi_omp/60, 'b.-', label='EMTransport')
# ax_plot_2.legend(frameon=False)
# ax_plot_2.set_xlabel('Number of OMP threads (cores=128)')
# ax_plot_2.set_ylabel('Time for first SCF step / minutes')
# fig_plot_2.tight_layout()
# fig_plot_2.savefig('{}/cx1_step0_all_omp.png'.format(folder_1), dpi=param.save_dpi)

# CX1 MPI total time
fig_plot_3, ax_plot_3 = plt.subplots()
# ax_plot_3.plot(cx1_cores, cx1_siesta_1_mpi/cx1_cores/60, 'r.-', label='CX1 BulkLeads')
# ax_plot_3.plot(cx1_cores, cx1_siesta_2_mpi/cx1_cores/60, 'g.-', label='CX1 EM DFT')
# ax_plot_3.plot(cx1_cores, cx1_siesta_3_mpi/cx1_cores/60, 'b.-', label='CX1 EMTransport')
ax_plot_3.plot(scarf_cores, scarf_siesta_1_mpi/scarf_cores/60, 'r.-', label='BulkLeads')
ax_plot_3.plot(scarf_cores, scarf_siesta_2_mpi/scarf_cores/60, 'g.-', label='EM DFT')
ax_plot_3.plot(scarf_cores, scarf_siesta_3_mpi/scarf_cores/60, 'b.-', label='EMTransport')
ax_plot_3.legend(frameon=False)
ax_plot_3.set_xlabel('Number of cores (threads=1)')
ax_plot_3.set_ylabel('Total time / minutes')
fig_plot_3.tight_layout()
# fig_plot_3.savefig('{}/cx1_total_mpi.png'.format(folder_1), dpi=param.save_dpi)
fig_plot_3.savefig('{}/scarf_total_mpi.png'.format(folder_2), dpi=param.save_dpi)

# CX1 OMP total time
# fig_plot_4, ax_plot_4 = plt.subplots()
# ax_plot_4.plot(cx1_threads_omp, cx1_siesta_1_omp/cx1_mpi_omp/60, 'r.-', label='BulkLeads')
# ax_plot_4.plot(cx1_threads_omp, cx1_siesta_2_omp/cx1_mpi_omp/60, 'g.-', label='EM DFT')
# ax_plot_4.plot(cx1_threads_omp, cx1_siesta_3_omp/cx1_mpi_omp/60, 'b.-', label='EMTransport')
# ax_plot_4.legend(frameon=False)
# ax_plot_4.set_xlabel('Number of OMP threads (cores=128)')
# ax_plot_4.set_ylabel('Total time / minutes')
# fig_plot_4.tight_layout()
# fig_plot_4.savefig('{}/cx1_total_omp.png'.format(folder_1), dpi=param.save_dpi)

# CX1 MPI TRANSM
# fig_plot_5, ax_plot_5 = plt.subplots()
# ax_plot_5.plot(cx1_cores, cx1_siesta_3_mpi/cx1_cores/60, 'r.-', label='EMTransport Total')
# ax_plot_5.plot(cx1_cores, cx1_transm_3_mpi/cx1_cores/60, 'g.-', label='EMTransport TRANSM')
# ax_plot_5.plot(cx1_cores, cx1_siesta_3_mpi/cx1_cores/60-cx1_transm_3_mpi/cx1_cores/60, 'b.-', label='EMTransport Total-TRANSM')
# ax_plot_5.legend(frameon=False)
# ax_plot_5.set_xlabel('Number of cores (threads=1)')
# ax_plot_5.set_ylabel('Time / minutes')
# fig_plot_5.tight_layout()
# fig_plot_5.savefig('{}/cx1_step0_all.png'.format(folder_1), dpi=param.save_dpi)

# CX1 OMP TRANSM
# fig_plot_6, ax_plot_6 = plt.subplots()
# ax_plot_6.plot(cx1_threads_omp, cx1_siesta_3_omp/cx1_mpi_omp/60, 'r.-', label='EMTransport Total')
# ax_plot_6.plot(cx1_threads_omp, cx1_transm_3_omp/cx1_mpi_omp/60, 'g.-', label='EMTransport TRANSM')
# ax_plot_6.plot(cx1_threads_omp, cx1_siesta_3_omp/cx1_mpi_omp/60-cx1_transm_3_omp/cx1_mpi_omp/60, 'b.-', label='EMTransport Total-TRANSM')
# ax_plot_6.legend(frameon=False)
# ax_plot_6.set_xlabel('Number of OMP threads (cores=128)')
# ax_plot_6.set_ylabel('Time / minutes')
# fig_plot_6.tight_layout()
# fig_plot_6.savefig('{}/cx1_step0_all.png'.format(folder_1), dpi=param.save_dpi)

# Plot all
# rows, cols = 2, 2
# fig_plot_all, ax_plot_all = plt.subplots(rows, cols,  figsize=(10, 8))
#
# ax_plot_all[0, 0].plot(cx1_cores_1, cx1_scf_3/cx1_cores_1, 'r.-', label='CX1 siesta')
# ax_plot_all[0, 0].plot(cx1_cores_1, cx1_scf_3/cx1_cores_1, 'r.-', label='CX1 IterSCF')
# ax_plot_all[0, 0].plot(cx1_cores_1, cx1_scf_3/cx1_cores_1, 'g.-', label='SCARF siesta')
# ax_plot_all[0, 0].plot(cx1_cores_1, cx1_scf_3/cx1_cores_1, 'g.-', label='SCARF IterSCF')
#
# ax_plot_all[0, 1].plot(cx1_cores_1, cx1_scf_3/cx1_cores_1, 'r.-', label='CX1 siesta')
# ax_plot_all[0, 1].plot(cx1_cores_1, cx1_scf_3/cx1_cores_1, 'r.-', label='CX1 IterSCF')
# ax_plot_all[0, 1].plot(cx1_cores_1, cx1_scf_3/cx1_cores_1, 'g.-', label='SCARF siesta')
# ax_plot_all[0, 1].plot(cx1_cores_1, cx1_scf_3/cx1_cores_1, 'g.-', label='SCARF IterSCF')
#
# ax_plot_all[1, 0].plot(cx1_cores_1, cx1_scf_3/cx1_cores_1, 'r.-', label='SCARF IterSCF')
# ax_plot_all[1, 0].plot(cx1_cores_1, cx1_scf_3/cx1_cores_1, 'g.-', label='SCARF IterSCF')
#
# ax_plot_all[1, 1].plot(cx1_cores_2, cx1_siesta/cx1_cores_2, 'k.-', label='CX1 siesta')
# ax_plot_all[1, 1].plot(cx1_cores_2, cx1_TRANSM/cx1_cores_2, 'b.-', label='CX1 TRANSM')
# ax_plot_all[1, 1].plot(cx1_cores_2, cx1_IterSCF/cx1_cores_2, 'r.-', label='CX1 IterSCF')
# ax_plot_all[1, 1].plot(cx1_cores_2, cx1_emtrans/cx1_cores_2, 'g.-', label='CX1 emtrans')

# ax_plot_all[0, 1].plot(file_1_2[:, 0]+fermi_offset_1, file_1_2[:, 1], 'r-', label=labels[0])
# ax_plot_all[0, 1].plot(file_2_2[:, 0]+fermi_offset_2, file_2_2[:, 1], 'g-', label=labels[1])
# ax_plot_all[0, 1].set_xlim([xlim[0], xlim[1]])
# ax_plot_all[0, 1].legend(frameon=False)
# # ax_plot_all[0, 1].set_xlabel(r'E-E$_{\mathrm{F}}$ (eV)')
# ax_plot_all[0, 1].set_ylabel('Number of channels')
# ax_plot_all[1, 0].plot(file_1_4[:, 0]+fermi_offset_1, file_1_4[:, 1], 'r-', label=labels[0])
# ax_plot_all[1, 0].plot(file_2_4[:, 0]+fermi_offset_2, file_2_4[:, 1], 'g-', label=labels[1])
# ax_plot_all[1, 0].set_xlim([xlim[0], xlim[1]])
# ax_plot_all[1, 0].legend(frameon=False)
# ax_plot_all[1, 0].set_xlabel(r'E-E$_{\mathrm{F}}$ (eV)')
# ax_plot_all[1, 0].set_ylabel('EM DOS')
# ax_plot_all[1, 1].plot(file_1_5[:, 0]+fermi_offset_1, file_1_5[:, 1], 'r-', label=labels[0])
# ax_plot_all[1, 1].plot(file_2_5[:, 0]+fermi_offset_2, file_2_5[:, 1], 'g-', label=labels[1])
# ax_plot_all[1, 1].set_xlim([xlim[0], xlim[1]])
# ax_plot_all[1, 1].legend(frameon=False)
# ax_plot_all[1, 1].set_xlabel(r'E-E$_{\mathrm{F}}$ (eV)')
# ax_plot_all[1, 1].set_ylabel('Leads DOS')
# fig_plot_all.tight_layout()
# # fig_plot_all.subplots_adjust(hspace=0)
# fig_plot_all.savefig('{}/compare-plot_all.png'.format(folder_1), dpi=param.save_dpi)
# fig_plot_all.savefig('{}/compare-plot_all.png'.format(folder_2), dpi=param.save_dpi)

if __name__ == "__main__":
    print('Finished.')
    plt.show()
