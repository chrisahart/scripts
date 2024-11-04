from __future__ import division, print_function
import pandas as pd
import numpy as np
from general import parameters as param
import matplotlib.pyplot as plt

"""
    Sort .xyz sequentially along desired axes (useful for transport calculations)
"""

folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/group/denan_li/archer/l75/pto/convergence/kpoints/bs-ti4-pb2-o2'
file_input_1 = '{}/data.out'.format(folder_1)
file_plot_energy_1 = '{}/plot_energy.png'.format(folder_1)
file_plot_scf_1 = '{}/plot_scf.png'.format(folder_1)
file_plot_time_1 = '{}/plot_time.png'.format(folder_1)

folder_2 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/group/denan_li/archer/l75/pto/convergence/kpoints/bs-ti4-pb2-o2'
file_input_2 = '{}/data.out'.format(folder_2)
file_plot_energy_2 = '{}/plot_energy.png'.format(folder_2)
file_plot_scf_2 = '{}/plot_scf.png'.format(folder_2)
file_plot_time_2 = '{}/plot_time.png'.format(folder_2)

x_label = 'k-points'
save_plot = True

cols = ['ii', 'energy', 'scf', 'time']
file_coord_1 = pd.read_csv(file_input_1, names=cols, delim_whitespace=True)
file_coord_2 = pd.read_csv(file_input_2, names=cols, delim_whitespace=True)

# Delete rows SCF=0
file_coord_1 = file_coord_1.drop(fileß_coord_1[file_coord_1.scf == 0].index)
file_coord_1 = file_coord_1.reset_index(drop=True)
# print(file_coord)

# Plot ii vs energy
fig_plot_energy, ax_plot_energy = plt.subplots()
ax_plot_energy.plot(file_coord_1['ii'], (file_coord_1['energy']-file_coord_1['energy'].values[-1])*param.hartree_to_ev*1e3, 'kx-')
ax_plot_energy.set_xlabel(x_label)
ax_plot_energy.set_ylabel('Energy / meV')
fig_plot_energy.tight_layout()
if save_plot: fig_plot_energy.savefig(file_plot_energy_1, dpi=param.save_dpi)
if save_plot: fig_plot_energy.savefig(file_plot_energy_2, dpi=param.save_dpi)

# # Plot ii vs scf steps
fig_plot_scf, ax_plot_scf = plt.subplots()
ax_plot_scf.plot(file_coord_1['ii'], file_coord_1['scf'], 'kx-')
ax_plot_scf.set_xlabel(x_label)
ax_plot_scf.set_ylabel('SCF steps')
fig_plot_scf.tight_layout()
if save_plot: fig_plot_scf.savefig(file_plot_scf_1, dpi=param.save_dpi)
if save_plot: fig_plot_scf.savefig(file_plot_scf_2, dpi=param.save_dpi)

# Plot ii vs time
fig_plot_time, ax_plot_time = plt.subplots()
ax_plot_time.plot(file_coord_1['ii'], file_coord_1['time'], 'kx-')
ax_plot_time.set_xlabel(x_label)
ax_plot_time.set_ylabel('Time / s')
fig_plot_time.tight_layout()
if save_plot: fig_plot_time.savefig(file_plot_time_1, dpi=param.save_dpi)
if save_plot: fig_plot_time.savefig(file_plot_time_2, dpi=param.save_dpi)

if __name__ == "__main__":
    print('Finished.')
    plt.show()
