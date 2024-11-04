from __future__ import division, print_function
import numpy as np
import matplotlib.pyplot as plt

"""
    Plot convergence
"""

folder_1 = '/Volumes/Storage/Data/Work/Postdoc/Work/testing/cp2k-smeagol/bdt/cp2k/bulk/convergence'
file_1_1 = np.genfromtxt('{}/Cutoff/output_cutoff.dat'.format(folder_1))
file_1_2 = np.genfromtxt('{}/Kpoints/output_kpoints.dat'.format(folder_1))

# Plotting parameters
save_dpi = 500  # DPI to save figure
params = {'axes.formatter.limits': [-4, 4],
          'axes.labelsize': 'large',
          'axes.titlesize': 'large',
          'lines.markersize': '8',
          }
plt.rcParams.update(params)

# Cutoff against energy
start = 3
fig_plot_1, ax_plot_1 = plt.subplots()
ax_plot_1.plot(file_1_1[start:, 0], 1*(file_1_1[start:, 2]-file_1_1[start:, 2][-1]), 'kx-')
ax_plot_1.set_xlabel('Cutoff')
ax_plot_1.set_ylabel('Energy / eV')
fig_plot_1.tight_layout()
fig_plot_1.savefig('{}/plot_cutoff_energy.png'.format(folder_1), dpi=save_dpi)

# Cutoff against time
fig_plot_2, ax_plot_2 = plt.subplots()
ax_plot_2.plot(file_1_1[start:, 0], file_1_1[start:, 3], 'kx-')
ax_plot_2.set_xlabel('Cutoff')
ax_plot_2.set_ylabel('Time / s')
fig_plot_2.tight_layout()
fig_plot_2.savefig('{}/plot_cutoff_time.png'.format(folder_1), dpi=save_dpi)

# Cutoff against energy
start = 1
fig_plot_3, ax_plot_3 = plt.subplots()
ax_plot_3.plot(file_1_2[start:, 0], 1*(file_1_2[start:, 2]-file_1_2[start:, 2][-1]), 'kx-')
ax_plot_3.set_xlabel('Kpoints')
ax_plot_3.set_ylabel('Energy / eV')
fig_plot_3.tight_layout()
fig_plot_3.savefig('{}/plot_kpoints_energy.png'.format(folder_1), dpi=save_dpi)

# Cutoff against time
fig_plot_4, ax_plot_4 = plt.subplots()
ax_plot_4.plot(file_1_2[start:, 0], file_1_2[start:, 3], 'kx-')
ax_plot_4.set_xlabel('Kpoints')
ax_plot_4.set_ylabel('Time / s')
fig_plot_4.tight_layout()
fig_plot_4.savefig('{}/plot_kpoints_time.png'.format(folder_1), dpi=save_dpi)

if __name__ == "__main__":
    print('Finished.')
    plt.show()
