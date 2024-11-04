from __future__ import division, print_function
import pandas as pd
import numpy as np
from general import parameters as param
import matplotlib.pyplot as plt

"""
    Sort .xyz sequentially along desired axes (useful for transport calculations)
"""

folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/bulk/monoclinic/convergence/monoclinic-cell-ref-large-rs-equal-cell/basis-sr'
file_input_1 = '{}/data.out'.format(folder_1)
file_plot_energy_1 = '{}/plot_energy.png'.format(folder_1)
file_plot_scf_1 = '{}/plot_scf.png'.format(folder_1)
file_plot_time_1 = '{}/plot_time.png'.format(folder_1)
x_label_1 = ['DZP SR', 'TZP SR', 'TZ2P SR']

folder_2 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/bulk/monoclinic/convergence/monoclinic-cell-ref-large-rs-equal-cell/basis'
file_input_2 = '{}/data.out'.format(folder_2)
file_plot_energy_2 = '{}/plot_energy.png'.format(folder_2)
file_plot_scf_2 = '{}/plot_scf.png'.format(folder_2)
file_plot_time_2 = '{}/plot_time.png'.format(folder_2)
x_label_2 = ['DZP', 'TZP', 'TZ2P']

save_plot = True
cols = ['ii', 'energy', 'scf', 'time']
file_coord_1 = pd.read_csv(file_input_1, names=cols, delim_whitespace=True)
file_coord_2 = pd.read_csv(file_input_2, names=cols, delim_whitespace=True)

# Delete rows SCF=0
file_coord_1 = file_coord_1.drop(file_coord_1[file_coord_1.scf == 0].index)
file_coord_1 = file_coord_1.reset_index(drop=True)
# print(file_coord)

# Plot ii vs energy
fig_plot_energy, ax_plot_energy = plt.subplots()
ax_plot_energy.plot(x_label_1, (file_coord_1['energy']-file_coord_1['energy'].values[-1])*param.hartree_to_ev*1e3, 'rx-')
ax_plot_energy.plot(x_label_2, (file_coord_2['energy']-file_coord_2['energy'].values[-1])*param.hartree_to_ev*1e3, 'bx-')
ax_plot_energy.set_ylabel('Energy / meV')
fig_plot_energy.tight_layout()
if save_plot: fig_plot_energy.savefig(file_plot_energy_1, dpi=param.save_dpi)
if save_plot: fig_plot_energy.savefig(file_plot_energy_2, dpi=param.save_dpi)

# # Plot ii vs scf steps
fig_plot_scf, ax_plot_scf = plt.subplots()
ax_plot_scf.plot(x_label_1, file_coord_1['scf'], 'rx-')
ax_plot_scf.plot(x_label_2, file_coord_2['scf'], 'bx-')
ax_plot_scf.set_ylabel('SCF steps')
fig_plot_scf.tight_layout()
if save_plot: fig_plot_scf.savefig(file_plot_scf_1, dpi=param.save_dpi)
if save_plot: fig_plot_scf.savefig(file_plot_scf_2, dpi=param.save_dpi)

# Plot ii vs time
fig_plot_time, ax_plot_time = plt.subplots()
ax_plot_time.plot(x_label_1, file_coord_1['time'], 'rx-')
ax_plot_time.plot(x_label_2, file_coord_2['time'], 'bx-')
ax_plot_time.set_ylabel('Time / s')
fig_plot_time.tight_layout()
if save_plot: fig_plot_time.savefig(file_plot_time_1, dpi=param.save_dpi)
if save_plot: fig_plot_time.savefig(file_plot_time_2, dpi=param.save_dpi)

if __name__ == "__main__":
    print('Finished.')
    plt.show()
