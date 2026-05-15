import numpy as np
from matplotlib import pyplot as plt
from general import parameters as param

"""
    Plot energy and forces for bulk hematite
"""

plotting_colors = ['r', 'g', 'b', 'm', 'grey', 'orange', 'brown', 'hotpink']

# Rutile 336
# folder_save = '/Volumes/Elements/Data/Postdoc2/Data/Work/calculations/yungu/images/336'
#
# intel_cores = np.array([96, 192, 288, 384])
# intel_time_1 = np.array([150, 77.1, 59.8, 46.3])
# intel_time_2 = np.array([54.6, 17.7, 13.1, 10.7])
#
# yungyu_cores = np.array([192, 192*2, 192*3, 192*4])
# yungyu_time_1 = np.array([92.4, 55.5, 66.4, 93.3])
# yungyu_time_2 = np.array([17.4, 11.6, 15.6, 12.3])
#
# yungyu_cores_omp = np.array([192, 192*2, 192*3, 192*4])
# yungyu_time_1_omp = np.array([179.6, 95.8, 76.4, 73.6])
# yungyu_time_2_omp = np.array([59.8, 16.5, 14.2, 12.4])
#
# archer_cores = np.array([128*2, 128*3, 128*4, 128*5, 128*6])
# archer_time_1 = np.array([130, 98.1, 88.1, 89, 84.9])
# archer_time_2 = np.array([28, 19.9, 15.8, 13.4, 11.8])
#
# archer_cores_omp = np.array([128*2, 128*3, 128*4, 128*5, 128*6])
# archer_time_1_omp = np.array([128.3, 87.6, 70.2, 60.6, 53.7])
# archer_time_2_omp = np.array([27.1, 19.2, 14.9, 12.7, 10.8])
#
# fugaku_omp_12_cores = np.array([48*16, 48*32, 48*64, 48*96])
# # fugaku_omp_12_time_1 = np.array([148.0, 99.4, 60.0, 49.9])
# # fugaku_omp_12_time_2 = np.array([73.1, 24.9, 18.9, 17.2])
# fugaku_omp_12_time_1 = np.array([154.6, 99.3, 60.3, 49.6])
# fugaku_omp_12_time_2 = np.array([33.0, 24.6, 19.0, 18.2])

# Rutile 446
# folder_save = '/Volumes/Elements/Data/Postdoc2/Data/Work/calculations/yungu/images/446'
#
# yungyu_cores = np.array([192*1, 192*2, 192*3, 192*4])
# yungyu_time_1 = np.array([121.9, 122.7, 83.4, 99.3])
# yungyu_time_2 = np.array([50.2, 33.5, 20.7, 26.6])
#
# yungyu_cores_omp_2 = np.array([192*2, 192*3, 192*4, 192*5])
# yungyu_time_1_omp_2 = np.array([165.7, 115.8, 84, 77.6])
# yungyu_time_2_omp_2 = np.array([57.4, 31.6, 23.3, 27.4])
#
# yungyu_cores_omp_12 = np.array([192*2, 192*3, 192*4, 192*5])
# yungyu_time_1_omp_12 = np.array([165.7, 115.8, 84, 77.6])
# yungyu_time_2_omp_12 = np.array([57.4, 31.6, 23.3, 27.4])
#
# archer_cores = np.array([128*5, 128*6])
# archer_time_1 = np.array([101.8, 107.7])
# archer_time_2 = np.array([23.9, 24.0])
#
# archer_cores_omp_2 = np.array([3*128, 4*128, 5*128, 6*128])
# archer_time_1_omp_2 = np.array([115.8, 89, 76.4, 76.3])
# archer_time_2_omp_2 = np.array([43.6, 27.7, 22.5, 23.1])
#
# archer_cores_omp_8 = np.array([2*128, 3*128, 4*128, 5*128, 6*128])
# archer_time_1_omp_8 = np.array([180.2, 130.4, 100.1, 84.8, 80.9])
# archer_time_2_omp_8 = np.array([124.5, 82.5, 60.1, 46.9, 41.1])
#
# # benchmarking
# fugaku_omp_12_cores = np.array([48*8, 48*16, 48*32, 48*48, 48*64, 48*80, 48*96])
# fugaku_omp_12_time_1 = np.array([428.5, 251.1, 138.9, 138.9, 110.5, 97.2, 90.4])
# fugaku_omp_12_cores = np.array([48*16, 48*32, 48*64, 48*96])
# fugaku_omp_12_time_1 = np.array([251.1, 138.9, 110.5, 90.4])
# fugaku_omp_12_time_2 = np.array([163.3, 74.7, 39.5, 32.9])
#
# # benchmarking-memory-8000
# fugaku_omp_12_cores = np.array([48*32, 48*64, 48*96])
# fugaku_omp_12_time_1 = np.array([142.8, 112.6, 102.7])
# fugaku_omp_12_time_2 = np.array([36.6, 36.4, 35.5])
#
# # benchmarking
# # fugaku_omp_8_cores = np.array([48*16, 48*32, 48*64, 48*96])
# # fugaku_omp_8_time_1 = np.array([253.1, 142.2, 116.1, 90.6])
#
# # benchmarking-memory-8000
# fugaku_omp_8_cores = np.array([48*64, 48*96])
# fugaku_omp_8_time_1 = np.array([118.8, 95.7])
# fugaku_omp_8_time_2 = np.array([41.5, 38.6])

# Rutile 446 PBE
# folder_save = '/Volumes/Elements/Data/Postdoc2/Data/Work/calculations/yungu/images/446-pbe'
# plotting_colors = ['r', 'g', 'b', 'm', 'grey', 'orange', 'brown', 'hotpink']
#
# archer_cores = np.array([128*1, 128*2, 128*3, 128*4, 128*5])
# archer_time_1_omp_1 = np.array([15.1, 8.8, 12.0, 17.2, 21.2])
# archer_time_1_omp_2 = np.array([16.3, 9.4, 7.2, 5.7, 7.8])
# archer_time_1_omp_8 = np.array([25.9, 14.9, 11.4, 8.9, 9.6])
#
# fugaku_omp_12_cores = np.array([48*4, 48*8, 48*16, 48*32])
# fugaku_omp_12_time_1 = np.array([33.9, 19.6, 12.6, 9.2])
# fugaku_omp_4_cores = np.array([48*4, 48*8, 48*16, 48*32])
# fugaku_omp_4_time_1 = np.array([25.4, 15.9, 12.0, 10.9])
# fugaku_omp_1_cores = np.array([48*4, 48*8, 48*16, 48*32])
# fugaku_omp_1_time_1 = np.array([25.9, 20.8, 21.6, 26.9])
#
# fugaku_omp_12_cores = np.array([48*4, 48*8, 48*16])
# fugaku_omp_12_time_1 = np.array([33.9, 19.6, 12.6])
# fugaku_omp_4_cores = np.array([48*4, 48*8, 48*16])
# fugaku_omp_4_time_1 = np.array([25.4, 15.9, 12.0])
# fugaku_omp_1_cores = np.array([48*4, 48*8, 48*16])
# fugaku_omp_1_time_1 = np.array([25.9, 20.8, 21.6])

# rutile(110)/water HSE
folder_save = '/Volumes/Elements/Data/Postdoc2/Data/Work/calculations/tio2-h2o/archer/ahart/110/benchmarking/pbe'
plotting_colors = ['r', 'g', 'b', 'm', 'grey', 'orange', 'brown', 'hotpink']

archer_cores = np.array([128*3, 128*4, 128*5])
archer_time_1_omp_1 = np.array([116.9, 101.8, 98.5])
archer_time_1_omp_2 = np.array([112.5, 86.9, 73.9])
archer_time_1_omp_8 = np.array([123.1, 94.7, 80.4])
archer_time_2_omp_1 = np.array([30.9, 24.2, 20.9])
archer_time_2_omp_2 = np.array([29.5, 23.2, 22.3])
archer_time_2_omp_8 = np.array([65.5, 45.0, 32.5])

# MAX_MEMORY 2400
fugaku_omp_12_cores = np.array([48*16, 48*32, 48*48, 48*64])
fugaku_omp_8_cores = np.array([48*32, 48*48, 48*64])
fugaku_omp_6_cores = np.array([48*32, 48*48, 48*64])
fugaku_omp_4_cores = np.array([48*64])

fugaku_omp_12_time_1 = np.array([366.5, 197.9, 167.9, 134.8])
fugaku_omp_8_time_1 = np.array([202.3, 172.1, 142.4])
fugaku_omp_6_time_1 = np.array([151.2, 133.7, 114.4])
fugaku_omp_4_time_1 = np.array([104.0])

fugaku_omp_12_time_2 = np.array([162.5, 44.0, 42.7, 37.3])
fugaku_omp_8_time_2 = np.array([47.5, 46.1, 41.2])
fugaku_omp_6_time_2 = np.array([43.0, 44.2, 42.3])
fugaku_omp_4_time_2 = np.array([37.3])

# MAX_MEMORY 8000 freq-2200
# fugaku_omp_12_cores = np.array([48*32, 48*48, 48*64])
fugaku_omp_8_cores = np.array([48*32, 48*48, 48*64])
fugaku_omp_6_cores = np.array([48*32, 48*48, 48*64])
fugaku_omp_4_cores = np.array([48*48, 48*64])

# fugaku_omp_12_time_1 = np.array([200.0, 172.1, 136.1])
fugaku_omp_8_time_1 = np.array([200.2, 172.3, 142.6])
fugaku_omp_6_time_1 = np.array([148.8, 130.4, 109.8])
fugaku_omp_4_time_1 = np.array([118.4, 103.0])

# fugaku_omp_12_time_2 = np.array([45.3, 45.7, 37.8])
fugaku_omp_8_time_2 = np.array([46.3, 48.1, 41.9])
fugaku_omp_6_time_2 = np.array([41.4, 42.2, 40.2])
fugaku_omp_4_time_2 = np.array([42.0, 43.5])

# first scf step
# fig_benchmark_1, ax_benchmark_1 = plt.subplots()
fig_benchmark_1, ax_benchmark_1 = plt.subplots(figsize=(8, 6))
# ax_benchmark_1.plot(archer_cores, archer_time_1, 'kx-', label='Archer AMD 7742 OMP = 1')
# ax_benchmark_1.plot(archer_cores_omp, archer_time_1_omp, 'ko-', label='Archer AMD 7742 OMP = 2', fillstyle='none')
# ax_benchmark_1.plot(intel_cores, intel_time_1, 'x-', label='Yungu Intel')
# ax_benchmark_1.plot(yungyu_cores, yungyu_time_1, 'rx-', label='Yungu AMD 9654 OMP = 1')
# ax_benchmark_1.plot(yungyu_cores_omp, yungyu_time_1_omp, 'go-', label='Yungu AMD 9654 OMP = 2', fillstyle='none')
# ax_benchmark_1.plot(fugaku_cores, fugaku_time_1, 'bo-', label='Fugaku AMD OMP = 12', fillstyle='none')

# ax_benchmark_1.plot(archer_cores_omp, archer_time_1_omp, 'ro-', label='ARCHER2 AMD 7742', fillstyle='none')
# ax_benchmark_1.plot(yungyu_cores_omp_2, yungyu_time_1_omp_2, 'go-', label='Yungyu AMD 9654', fillstyle='none')
# ax_benchmark_1.plot(fugaku_omp_12_cores, fugaku_omp_12_time_1, 'bo-', label='Fugaku A64FX', fillstyle='none')

# ax_benchmark_1.plot(archer_cores_omp_2, archer_time_1_omp_2, 'o-', color=plotting_colors[0], label='ARCHER2 AMD 7742 OMP = 2', fillstyle='none')
# ax_benchmark_1.plot(archer_cores_omp_8, archer_time_1_omp_8, 'o-', color=plotting_colors[1], label='ARCHER2 AMD 7742 OMP = 8', fillstyle='none')
# ax_benchmark_1.plot(yungyu_cores_omp_2, yungyu_time_1_omp_2, 'o-', color=plotting_colors[2], label='Yungyu AMD 9654 OMP = 2', fillstyle='none')
# ax_benchmark_1.plot(fugaku_omp_12_cores, fugaku_omp_12_time_1, 'o-', color=plotting_colors[3], label='Fugaku A64FX OMP = 12', fillstyle='none')
# ax_benchmark_1.plot(fugaku_omp_8_cores, fugaku_omp_8_time_1, 'o-', color=plotting_colors[4], label='Fugaku A64FX OMP = 8', fillstyle='none')

# ax_benchmark_1.plot(archer_cores, archer_time_1_omp_1, 'o-', color=plotting_colors[0], label='ARCHER2 AMD 7742 OMP = 1', fillstyle='none')
ax_benchmark_1.plot(archer_cores, archer_time_1_omp_2, 's-', color=plotting_colors[0], label='ARCHER2 AMD 7742 OMP = 2', fillstyle='none')
# ax_benchmark_1.plot(archer_cores, archer_time_1_omp_8, 'o-', color=plotting_colors[2], label='ARCHER2 AMD 7742 OMP = 8', fillstyle='none')
# ax_benchmark_1.plot(fugaku_omp_1_cores, fugaku_omp_1_time_1, 'o-', color=plotting_colors[3], label='Fugaku A64FX OMP = 1', fillstyle='none')
ax_benchmark_1.plot(fugaku_omp_12_cores, fugaku_omp_12_time_1, 'o-', color=plotting_colors[4], label='Fugaku A64FX OMP = 12', fillstyle='none')
# ax_benchmark_1.plot(fugaku_omp_8_cores, fugaku_omp_8_time_1, 'o-', color=plotting_colors[3], label='Fugaku A64FX OMP = 8', fillstyle='none')
ax_benchmark_1.plot(fugaku_omp_6_cores, fugaku_omp_6_time_1, 'o-', color=plotting_colors[2], label='Fugaku A64FX OMP = 6', fillstyle='none')
ax_benchmark_1.plot(fugaku_omp_4_cores, fugaku_omp_4_time_1, 'o-', color=plotting_colors[1], label='Fugaku A64FX OMP = 4', fillstyle='none')

ax_benchmark_1.set_ylabel('Time for first SCF step / s')
ax_benchmark_1.set_xlabel('Number of CPU')
ax_benchmark_1.legend(frameon=False)
ax_benchmark_1.set_ylim([0, 380])
fig_benchmark_1.tight_layout()
fig_benchmark_1.savefig('{}/benchmark_1.png'.format(folder_save), dpi=300)

# second scf step
# fig_benchmark_2, ax_benchmark_2 = plt.subplots()
fig_benchmark_2, ax_benchmark_2 = plt.subplots(figsize=(8, 6))
# ax_benchmark_2.plot(archer_cores, archer_time_2, 'kx-', label='Archer AMD 7742 OMP = 1')
# ax_benchmark_2.plot(archer_cores_omp, archer_time_2_omp, 'ko-', label='Archer AMD 7742 OMP = 2', fillstyle='none')
# ax_benchmark_2.plot(intel_cores, intel_time_2, 'x-', label='Yungu Intel')
# ax_benchmark_2.plot(yungyu_cores, yungyu_time_2, 'rx-', label='Yungu AMD 9654 OMP = 1')
# ax_benchmark_2.plot(yungyu_cores_omp, yungyu_time_2_omp, 'go-', label='Yungu AMD 9654 OMP = 2', fillstyle='none')
# ax_benchmark_2.plot(fugaku_cores, fugaku_time_2, 'bo-', label='Fugaku AMD OMP = 12', fillstyle='none')

# ax_benchmark_2.plot(archer_cores_omp, archer_time_2_omp, 'ro-', label='ARCHER2 AMD 7742', fillstyle='none')
# ax_benchmark_2.plot(yungyu_cores_omp_2, yungyu_time_2_omp_2, 'go-', label='Yungyu AMD 9654', fillstyle='none')
# ax_benchmark_2.plot(fugaku_omp_12_cores, fugaku_omp_12_time_2, 'bo-', label='Fugaku A64FX', fillstyle='none')

ax_benchmark_2.plot(archer_cores, archer_time_2_omp_2, 's-', color=plotting_colors[0], label='ARCHER2 AMD 7742 OMP = 2', fillstyle='none')
# ax_benchmark_2.plot(archer_cores_omp_8, archer_time_2_omp_8, 'o-', color=plotting_colors[1], label='ARCHER2 AMD 7742 OMP = 8', fillstyle='none')
# ax_benchmark_2.plot(yungyu_cores_omp_2, yungyu_time_2_omp_2, 'o-', color=plotting_colors[2], label='Yungyu AMD 9654 OMP = 2', fillstyle='none')
ax_benchmark_2.plot(fugaku_omp_12_cores, fugaku_omp_12_time_2, 'o-', color=plotting_colors[4], label='Fugaku A64FX OMP = 12', fillstyle='none')
# ax_benchmark_2.plot(fugaku_omp_8_cores, fugaku_omp_8_time_2, 'o-', color=plotting_colors[4], label='Fugaku A64FX OMP = 8', fillstyle='none')
ax_benchmark_2.plot(fugaku_omp_6_cores, fugaku_omp_6_time_2, 'o-', color=plotting_colors[2], label='Fugaku A64FX OMP = 6', fillstyle='none')
ax_benchmark_2.plot(fugaku_omp_4_cores, fugaku_omp_4_time_2, 'o-', color=plotting_colors[1], label='Fugaku A64FX OMP = 4', fillstyle='none')
#
ax_benchmark_2.set_ylabel('Time for second SCF step / s')
ax_benchmark_2.set_xlabel('Number of CPU')
ax_benchmark_2.legend(frameon=False)
ax_benchmark_2.set_ylim([0, 170])
fig_benchmark_2.tight_layout()
fig_benchmark_2.savefig('{}/benchmark_2.png'.format(folder_save), dpi=300)

# rows, cols = 1, 2
# fig_plot_all, ax_plot_all = plt.subplots(rows, cols, sharex='col', sharey='row', figsize=(10, 6))
# ax_plot_all[0].plot(intel_cores, intel_time_1, 'x-', label='Yungu Intel')
# ax_plot_all[0].plot(yungyu_cores, yungyu_time_1, 'x-', label='Yungu AMD')
# ax_plot_all[0].set_ylabel('Time for first HSE SCF step / s')
# ax_plot_all[0].set_xlabel('Number of CPU')
# ax_plot_all[0].legend(frameon=False)
# ax_plot_all[1].plot(intel_cores, intel_time_2, 'x-', label='Yungu Intel')
# ax_plot_all[1].plot(yungyu_cores, yungyu_time_2, 'x-', label='Yungu AMD')
# ax_plot_all[1].set_ylabel('Time for second HSE SCF step / s')
# ax_plot_all[1].set_xlabel('Number of CPU')
# ax_plot_all[1].legend(frameon=False)
# fig_plot_all.tight_layout()
# fig_plot_all.savefig('{}/benchmark.png'.format(folder_save), dpi=300)

if __name__ == "__main__":
    print('Finished.')
    plt.show()
