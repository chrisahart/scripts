import numpy as np
from matplotlib import pyplot as plt
from general import parameters as param

"""
    Plot energy and forces for bulk hematite
"""

# Rutile 336
# folder_save = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/yungu/images/336'

# intel_cores = np.array([96, 192, 288, 384])
# intel_time_1 = np.array([150, 77.1, 59.8, 46.3])
# intel_time_2 = np.array([54.6, 17.7, 13.1, 10.7])

# amd_cores = np.array([192, 192*2, 192*3, 192*4])
# amd_time_1 = np.array([92.4, 55.5, 66.4, 93.3])
# amd_time_2 = np.array([17.4, 11.6, 15.6, 12.3])
#
# amd_cores_omp = np.array([192, 192*2, 192*3, 192*4])
# amd_time_1_omp = np.array([179.6, 95.8, 76.4, 73.6])
# amd_time_2_omp = np.array([59.8, 16.5, 14.2, 12.4])
#
# archer_cores = np.array([128*2, 128*3, 128*4, 128*5, 128*6])
# archer_time_1 = np.array([130, 98.1, 88.1, 89, 84.9])
# archer_time_2 = np.array([28, 19.9, 15.8, 13.4, 11.8])
#
# archer_cores_omp = np.array([128*2, 128*3, 128*4, 128*5, 128*6])
# archer_time_1_omp = np.array([128.3, 87.6, 70.2, 60.6, 53.7])
# archer_time_2_omp = np.array([27.1, 19.2, 14.9, 12.7, 10.8])

# Rutile 446
folder_save = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/yungu/images/446'

amd_cores = np.array([192*1, 192*2, 192*3, 192*4, 192*5])
amd_time_1 = np.array([121.9, 122.7, 83.4, 99.3, 81.3])
amd_time_2 = np.array([50.2, 33.5, 20.7, 26.6, 28.7])

amd_cores = np.array([192*1, 192*2, 192*3, 192*4])
amd_time_1 = np.array([121.9, 122.7, 83.4, 99.3])
amd_time_2 = np.array([50.2, 33.5, 20.7, 26.6])

amd_cores_omp = np.array([192*2, 192*3, 192*4, 192*5])
amd_time_1_omp = np.array([165.7, 115.8, 84, 77.6])
amd_time_2_omp = np.array([57.4, 31.6, 23.3, 27.4])

archer_cores = np.array([128*5, 128*6])
archer_time_1 = np.array([101.8, 107.7])
archer_time_2 = np.array([23.9, 24.0])

archer_cores_omp = np.array([3*128, 4*128, 5*128, 6*128])
archer_time_1_omp = np.array([115.8, 89, 76.4, 76.3])
archer_time_2_omp = np.array([43.6, 27.7, 22.5, 23.1])

# first scf step
fig_benchmark_1, ax_benchmark_1 = plt.subplots()
# ax_benchmark_1.plot(archer_cores, archer_time_1, 'kx-', label='Archer AMD 7742 OMP = 1')
ax_benchmark_1.plot(archer_cores_omp, archer_time_1_omp, 'ko-', label='Archer AMD 7742 OMP = 2', fillstyle='none')
# ax_benchmark_1.plot(intel_cores, intel_time_1, 'x-', label='Yungu Intel')
ax_benchmark_1.plot(amd_cores, amd_time_1, 'rx-', label='Yungu AMD 9654 OMP = 1')
ax_benchmark_1.plot(amd_cores_omp, amd_time_1_omp, 'go-', label='Yungu AMD 9654 OMP = 2', fillstyle='none')

ax_benchmark_1.set_ylabel('Time for first HSE SCF step / s')
ax_benchmark_1.set_xlabel('Number of CPU')
ax_benchmark_1.legend(frameon=False)
fig_benchmark_1.tight_layout()
fig_benchmark_1.savefig('{}/benchmark_1.png'.format(folder_save), dpi=300)

# second scf step
fig_benchmark_2, ax_benchmark_2 = plt.subplots()
ax_benchmark_2.plot(archer_cores, archer_time_2, 'kx-', label='Archer AMD 7742 OMP = 1')
ax_benchmark_2.plot(archer_cores_omp, archer_time_2_omp, 'ko-', label='Archer AMD 7742 OMP = 2', fillstyle='none')
# ax_benchmark_2.plot(intel_cores, intel_time_2, 'x-', label='Yungu Intel')
ax_benchmark_2.plot(amd_cores, amd_time_2, 'rx-', label='Yungu AMD 9654 OMP = 1')
ax_benchmark_2.plot(amd_cores_omp, amd_time_2_omp, 'go-', label='Yungu AMD 9654 OMP = 2', fillstyle='none')

ax_benchmark_2.set_ylabel('Time for second HSE SCF step / s')
ax_benchmark_2.set_xlabel('Number of CPU')
ax_benchmark_2.legend(frameon=False)
fig_benchmark_2.tight_layout()
fig_benchmark_2.savefig('{}/benchmark_2.png'.format(folder_save), dpi=300)

# rows, cols = 1, 2
# fig_plot_all, ax_plot_all = plt.subplots(rows, cols, sharex='col', sharey='row', figsize=(10, 6))
# ax_plot_all[0].plot(intel_cores, intel_time_1, 'x-', label='Yungu Intel')
# ax_plot_all[0].plot(amd_cores, amd_time_1, 'x-', label='Yungu AMD')
# ax_plot_all[0].set_ylabel('Time for first HSE SCF step / s')
# ax_plot_all[0].set_xlabel('Number of CPU')
# ax_plot_all[0].legend(frameon=False)
# ax_plot_all[1].plot(intel_cores, intel_time_2, 'x-', label='Yungu Intel')
# ax_plot_all[1].plot(amd_cores, amd_time_2, 'x-', label='Yungu AMD')
# ax_plot_all[1].set_ylabel('Time for second HSE SCF step / s')
# ax_plot_all[1].set_xlabel('Number of CPU')
# ax_plot_all[1].legend(frameon=False)
# fig_plot_all.tight_layout()
# fig_plot_all.savefig('{}/benchmark.png'.format(folder_save), dpi=300)

if __name__ == "__main__":
    print('Finished.')
    plt.show()
