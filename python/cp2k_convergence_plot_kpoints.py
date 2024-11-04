from __future__ import division, print_function
import pandas as pd
import numpy as np
from general import parameters as param
import matplotlib.pyplot as plt

"""
    Sort .xyz sequentially along desired axes (useful for transport calculations)
"""

folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/bulk/monoclinic/convergence/monoclinic-cell-ref-large-rs-equal-cell/kpoints/MOLOPT-SR-GTH_GTH_POTENTIALS'
file_input = '{}/data.out'.format(folder)
file_plot_energy = '{}/plot_energy.png'.format(folder)
file_plot_scf = '{}/plot_scf.png'.format(folder)
file_plot_time = '{}/plot_time.png'.format(folder)
x_label = 'k-points'
plot_xlim = [1.8, 8.2]
plot_ylim = [-1.1, 1.1]
save_plot = True

cols = ['ii', 'energy', 'scf', 'time']
file_coord = pd.read_csv(file_input, names=cols, delim_whitespace=True)

# Delete rows SCF=0
file_coord = file_coord.drop(file_coord[file_coord.scf == 0].index)
file_coord = file_coord.reset_index(drop=True)
# print(file_coord)

# Plot ii vs energy
# fig_plot_energy, ax_plot_energy = plt.subplots()
energy = (file_coord['energy']-file_coord['energy'].values[-1])*param.hartree_to_ev * 1e3 / 12
fig_plot_energy, ax_plot_energy = plt.subplots(figsize=(6, 6))
ax_plot_energy.plot(file_coord['ii'], energy , 'kx-')
ax_plot_energy.hlines(1, 0, 100, 'r', alpha=0.5)
ax_plot_energy.hlines(-1, 0, 100, 'r', alpha=0.5)
ax_plot_energy.set_xlabel(x_label)
ax_plot_energy.set_ylabel('Energy / meV per atom')
ax_plot_energy.set_xlim([plot_xlim[0], plot_xlim[1]])
ax_plot_energy.set_ylim([plot_ylim[0], plot_ylim[1]])
fig_plot_energy.tight_layout()
if save_plot: fig_plot_energy.savefig(file_plot_energy, dpi=param.save_dpi)

# # Plot ii vs scf steps
# fig_plot_scf, ax_plot_scf = plt.subplots()
# ax_plot_scf.plot(file_coord['ii'], file_coord['scf'], 'kx-')
# ax_plot_scf.set_xlabel(x_label)
# ax_plot_scf.set_ylabel('SCF steps')
# fig_plot_scf.tight_layout()
# if save_plot: fig_plot_scf.savefig(file_plot_scf, dpi=param.save_dpi)
#
# # Plot ii vs time
# fig_plot_time, ax_plot_time = plt.subplots()
# ax_plot_time.plot(file_coord['ii'], file_coord['time'], 'kx-')
# ax_plot_time.set_xlabel(x_label)
# ax_plot_time.set_ylabel('Time / s')
# fig_plot_time.tight_layout()
# if save_plot: fig_plot_time.savefig(file_plot_time, dpi=param.save_dpi)

if __name__ == "__main__":
    print('Finished.')
    plt.show()
