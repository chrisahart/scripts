import numpy as np
from matplotlib import pyplot as plt
from ase.io.cube import read_cube_data
from general import parameters as param
import pandas as pd
from scipy.optimize import curve_fit


def f_1(x, A, B):
    return A * x ** 1 + B


def f_2(x, A, B):
    return A * x ** 2 + B


def f_3(x, A, B):
    return A * x ** 3 + B


# bulk-sz-q11_em-sz-q11
# folders = ['/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/performance/benchmarking/cp2k-smeagol/au-bdt/NEnergReal-64_kpoints-1-1-20/bulk-sz-q11_em-sz-q11/atoms-116',
#            '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/performance/benchmarking/cp2k-smeagol/au-bdt/NEnergReal-64_kpoints-1-1-20/bulk-sz-q11_em-sz-q11/atoms-246',
#            '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/performance/benchmarking/cp2k-smeagol/au-bdt/NEnergReal-64_kpoints-1-1-20/bulk-sz-q11_em-sz-q11/atoms-428',
#            '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/performance/benchmarking/cp2k-smeagol/au-bdt/NEnergReal-64_kpoints-1-1-20/bulk-sz-q11_em-sz-q11/atoms-662',
#            ]
# filename = '/log.out'
# cols = ['Filename', 'Time']
# data_atoms = [116, 246, 428, 662]
# data_basis = [1060, 2230, 3868, 5974]
# threads = np.array([1, 2, 4, 8, 16, 32, 64])

# bulk-sz-q1_em-sz-q11
folders = ['/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/performance/benchmarking/cp2k-smeagol/au-bdt/NEnergReal-64_kpoints-1-1-20/bulk-sz-q1_em-sz-q11/atoms-662',
           '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/performance/benchmarking/cp2k-smeagol/au-bdt/NEnergReal-64_kpoints-1-1-20/bulk-sz-q1_em-sz-q11/atoms-948',
           '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/performance/benchmarking/cp2k-smeagol/au-bdt/NEnergReal-64_kpoints-1-1-20/bulk-sz-q1_em-sz-q11/atoms-1286',
           '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/performance/benchmarking/cp2k-smeagol/au-bdt/NEnergReal-64_kpoints-1-1-20/bulk-sz-q1_em-sz-q11/atoms-1676',
           '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/performance/benchmarking/cp2k-smeagol/au-bdt/NEnergReal-64_kpoints-1-1-20/bulk-sz-q1_em-sz-q11/atoms-2118',
           '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/performance/benchmarking/cp2k-smeagol/au-bdt/NEnergReal-64_kpoints-1-1-20/bulk-sz-q1_em-sz-q11/atoms-2612',
           '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/performance/benchmarking/cp2k-smeagol/au-bdt/NEnergReal-64_kpoints-1-1-20/bulk-sz-q1_em-sz-q11/atoms-3756',
           ]
filename = '/log.out'
cols = ['Filename', 'Time']
data_atoms = [662, 948, 1286, 1678, 2118, 2612, 3756]
data_basis = [2774, 3940, 5318, 6908, 8710, 10724, 15388]
threads = np.array([1, 2, 4, 8, 16, 32, 64])
# threads = np.array([0.5, 1, 2, 4, 8, 16, 32]) * 128

diff_color = ['r', 'g', 'b', 'm', 'orange', 'y', 'brown', 'cyan', 'pink']

data = []
for i in range(len(folders)):
    print('Reading', folders[i])
    data.append(pd.read_csv('{}{}'.format(folders[i], filename), names=cols, delim_whitespace=True))

data_dft = []
data_threads_1 = []
data_threads_2 = []
data_threads_4 = []
data_threads_8 = []
data_threads_16 = []
data_threads_32 = []
data_threads_64 = []
for i in range(len(folders)):
    data_dft.append(data[i]['Time'][1])
    data_threads_1.append(data[i]['Time'][7 + 6 * 0])
    data_threads_2.append(data[i]['Time'][7 + 6 * 1])
    data_threads_4.append(data[i]['Time'][7 + 6 * 2])
    data_threads_8.append(data[i]['Time'][7 + 6 * 3])
    data_threads_16.append(data[i]['Time'][7 + 6 * 4])
    data_threads_32.append(data[i]['Time'][7 + 6 * 5])
    data_threads_64.append(data[i]['Time'][7 + 6 * 6])

# Plot threads vs time 1  for V=0
# fig_plot_1, ax_plot_1 = plt.subplots()
# for i in range(len(folders)):
#     ax_plot_1.plot(threads, data[i]['Time'][3::6], 'x-', label=data_basis[i], color=diff_color[i])
# ax_plot_1.set_xlabel('OMP threads')
# ax_plot_1.set_ylabel('Time for first SCF step / s')
# ax_plot_1.legend(frameon=False)
# fig_plot_1.tight_layout()
# for i in range(len(folders)):
#     fig_plot_1.savefig('{}/threads_time_V-0.png'.format(folders[i]), dpi=param.save_dpi)

# Plot threads vs time 1 for V=1
# fig_plot_2, ax_plot_2 = plt.subplots()
# for i in range(len(folders)):
#     ax_plot_2.plot(threads, data[i]['Time'][6::6], 'x-', label=data_basis[i], color=diff_color[i])
# ax_plot_2.set_xlabel('OMP threads')
# ax_plot_2.set_ylabel('Time for first SCF step / s')
# ax_plot_2.legend(frameon=False)
# fig_plot_2.tight_layout()
# for i in range(len(folders)):
#     fig_plot_2.savefig('{}/threads_time_V-1.png'.format(folders[i]), dpi=param.save_dpi)

# Plot threads vs time 1 for V=0 and V=1
fig_plot_3, ax_plot_3 = plt.subplots()
ax_plot_3.axhline(y=data_dft[0], color='k', linestyle='-', label='DFT')
for i in range(len(folders)):
    print(threads)
    print(data[i]['Time'][3::6])
    ax_plot_3.plot(threads, data[i]['Time'][3::6], '+--', color=diff_color[i])
    ax_plot_3.plot(threads, data[i]['Time'][6::6], 'x-', label=data_basis[i], color=diff_color[i])
ax_plot_3.set_xlabel('OMP threads')
ax_plot_3.set_ylabel('Time for first SCF step / s')
ax_plot_3.legend(frameon=False)
ax_plot_3.ticklabel_format(useOffset=False, style='plain')
fig_plot_3.tight_layout()
for i in range(len(folders)):
    fig_plot_3.savefig('{}/threads_time-1_V-0_V-1.png'.format(folders[i]), dpi=param.save_dpi)

# Plot threads vs time 2 for V=0 and V=1
fig_plot_4, ax_plot_4 = plt.subplots()
ax_plot_4.axhline(y=data_dft[1], color='k', linestyle='-', label='DFT')
for i in range(len(folders)):
    ax_plot_4.plot(threads, data[i]['Time'][4::6], '+--', color=diff_color[i])
    ax_plot_4.plot(threads, data[i]['Time'][7::6], 'x-', label=data_basis[i], color=diff_color[i])
ax_plot_4.set_xlabel('OMP threads')
ax_plot_4.set_ylabel('Time for second SCF step / s')
ax_plot_4.legend(frameon=False)
ax_plot_4.ticklabel_format(useOffset=False, style='plain')
fig_plot_4.tight_layout()
for i in range(len(folders)):
    fig_plot_4.savefig('{}/threads_time-2_V-0_V-1.png'.format(folders[i]), dpi=param.save_dpi)

# Plot threads vs time 1 for V=0 and V=1
# fig_plot_5, ax_plot_5 = plt.subplots()
# ax_plot_5.plot(data_basis[:-1], data_dft[:-1], 'x-', color=diff_color[0], label='DFT')
# ax_plot_5.plot(data_basis, data_threads_1, 'x-', color=diff_color[1], label='OMP threads 1')
# ax_plot_5.plot(data_basis, data_threads_2, 'x-', color=diff_color[2], label='OMP threads 2')
# ax_plot_5.plot(data_basis, data_threads_4, 'x-', color=diff_color[3], label='OMP threads 4')
# ax_plot_5.plot(data_basis, data_threads_8, 'x-', color=diff_color[4], label='OMP threads 8')
# ax_plot_5.plot(data_basis, data_threads_16, 'x-', color=diff_color[5], label='OMP threads 16')
# ax_plot_5.plot(data_basis, data_threads_32, 'x-', color=diff_color[6], label='OMP threads 32')
# ax_plot_5.plot(data_basis, data_threads_64, 'x-', color=diff_color[7], label='OMP threads 64')
# ax_plot_5.set_xlabel('Number basis functions')
# ax_plot_5.set_ylabel('Time for first SCF step / s')
# ax_plot_5.legend(frameon=False)
# fig_plot_5.tight_layout()
# for i in range(len(folders)):
#     fig_plot_5.savefig('{}/basis_time-1.png'.format(folders[i]), dpi=param.save_dpi)

# Plot threads vs time 2 for V=0 and V=1
fig_plot_6, ax_plot_6 = plt.subplots()
ax_plot_6.plot(data_basis, data_dft, 'k-', label='DFT')
# ax_plot_6.plot(data_basis[0:2], data_dft[0:2], 'k+')
# ax_plot_6.plot(data_basis[2:], data_dft[2:], 'kx')
ax_plot_6.plot(data_basis[0], data_dft[0], 'k+')
ax_plot_6.plot(data_basis[1:-1], data_dft[1:-1], 'kx')
ax_plot_6.plot(data_basis[-1], data_dft[-1], 'k', marker='|')
ax_plot_6.plot(data_basis, data_threads_1, 'x-', color=diff_color[1], label='OMP threads 1')
ax_plot_6.plot(data_basis, data_threads_2, 'x-', color=diff_color[2], label='OMP threads 2')
ax_plot_6.plot(data_basis, data_threads_4, 'x-', color=diff_color[3], label='OMP threads 4')
ax_plot_6.plot(data_basis, data_threads_8, 'x-', color=diff_color[4], label='OMP threads 8')
ax_plot_6.plot(data_basis, data_threads_16, 'x-', color=diff_color[5], label='OMP threads 16')
ax_plot_6.plot(data_basis, data_threads_32, 'x-', color=diff_color[6], label='OMP threads 32')
ax_plot_6.plot(data_basis, data_threads_64, 'x-', color=diff_color[7], label='OMP threads 64')
ax_plot_6.set_xlabel('Number basis functions')
ax_plot_6.set_ylabel('Time for second SCF step / s')
ax_plot_6.ticklabel_format(useOffset=False, style='plain')
ax_plot_6.legend(frameon=False)
fig_plot_6.tight_layout()
for i in range(len(folders)):
    fig_plot_6.savefig('{}/basis_time-2.png'.format(folders[i]), dpi=param.save_dpi)

# Speedup
fig_plot_7, ax_plot_7 = plt.subplots()
ax_plot_7.plot(threads[2:], threads[2:]/threads[2], 'k-', label='Ideal', alpha=0.5)
# i = 3
# ax_plot_7.plot(threads[2:], np.max(data[i]['Time'][7 + 6 * 2::6]) / data[i]['Time'][7 + 6 * 2::6], 'rx-')
# temp = np.max(data[i]['Time'][7 + 6 * 2::6]) / data[i]['Time'][7 + 6 * 2::6]
# speedup = temp / np.array([1, 2, 4, 8, 16])
# print(temp)
# print(speedup)
folders_plot = range(len(folders))
print(folders_plot)
folders_plot = folders_plot[1:]
for i in folders_plot:
    # ax_plot_7.plot(threads, np.max(data[i]['Time'][4::6])/data[i]['Time'][4::6], '+--', color=diff_color[i])
    # ax_plot_7.plot(threads, np.max(data[i]['Time'][7::6])/data[i]['Time'][7::6], 'x-', label=data_basis[i], color=diff_color[i])
    # ax_plot_7.plot(threads[2:], np.max(data[i]['Time'][4+6*2::6])/data[i]['Time'][4+6*2::6], '+--', color=diff_color[i])
    ax_plot_7.plot(threads[2:], np.max(data[i]['Time'][7+6*2::6])/data[i]['Time'][7+6*2::6], 'x-', label=data_basis[i], color=diff_color[i])
# ax_plot_7.set_xlabel('Archer2 cores')
ax_plot_7.set_xlabel('OMP threads')
ax_plot_7.set_ylabel('Speedup')
# ax_plot_7.set_ylabel('Speedup for second SCF step')
ax_plot_7.set_ylim([0.9, 3.6])
ax_plot_7.legend(frameon=False)
fig_plot_7.tight_layout()
for i in range(len(folders)):
    fig_plot_7.savefig('{}/threads_scaling_time-2_V-0_V-1.png'.format(folders[i]), dpi=param.save_dpi)
    
if __name__ == "__main__":
    print('Finished.')
    plt.show()
