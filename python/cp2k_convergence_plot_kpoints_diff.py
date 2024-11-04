from __future__ import division, print_function
import pandas as pd
import numpy as np
from general import parameters as param
import matplotlib.pyplot as plt

"""
    Sort .xyz sequentially along desired axes (useful for transport calculations)
"""

params = {'axes.formatter.limits': [-4, 4],
          'axes.labelsize': 'x-large',
          'xtick.labelsize' : '12',
          'axes.titlesize': 'x-large',
          'legend.fontsize': 'medium',
          'lines.markersize': '10',
          }
plt.rcParams.update(params)

folders = [
    '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/bulk/monoclinic/convergence/monoclinic-cell-ref-large-rs-equal-cell/kpoints/MOLOPT-SR-GTH_GTH_POTENTIALS',
    '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/bulk/tetragonal/convergence/yudi-pbe-cell-ref-large-rs-equal-cell/kpoints/MOLOPT-SR-GTH_GTH_POTENTIALS',
    '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/bulk/po/convergence/yudi-pbe-cell-ref-large-rs-equal-cell/kpoints/MOLOPT-SR-GTH_GTH_POTENTIALS'
]

# Process each folder
data_frames = []
for folder in folders:
    file_input = f'{folder}/data.out'
    cols = ['ii', 'energy', 'scf', 'time']
    df = pd.read_csv(file_input, names=cols, delim_whitespace=True)
    data_frames.append(df)

energy_1 = np.abs((data_frames[0]['energy']-data_frames[1]['energy'])*param.hartree_to_ev) / 12 * 1e3
energy_2 = np.abs((data_frames[0]['energy']-data_frames[2]['energy'])*param.hartree_to_ev) / 12 * 1e3
print(data_frames[0]['ii'])
print(energy_1)
print(energy_2)

plot_xlim_1 = [2.8, 7.2]
plot_xlim_2 = plot_xlim_1
plot_ylim_1 = np.array([0.547, 0.553]) / 12 * 1e3
plot_ylim_2 = np.array([0.2736, 0.2746]) / 12 * 1e3

# Plot ii vs energy 1
fig_plot_energy_1, ax_plot_energy_1 = plt.subplots(figsize=(6, 6))
ax_plot_energy_1.plot(data_frames[0]['ii'], energy_1, 'kx-', label='A')
ax_plot_energy_1.hlines(energy_1.values[-1]+1, 0, 100, 'r', alpha=0.5)
ax_plot_energy_1.hlines(energy_1.values[-1]-1, 0, 100, 'r', alpha=0.5)
ax_plot_energy_1.set_xlim([plot_xlim_1[0], plot_xlim_1[1]])
ax_plot_energy_1.set_ylim([plot_ylim_1[0], plot_ylim_1[1]])
ax_plot_energy_1.set_ylabel('Energy / meV per atom')
fig_plot_energy_1.tight_layout()
fig_plot_energy_1.savefig('{}/energy_1_diff.png'.format(folders[2]), dpi=param.save_dpi)

# Plot ii vs energy
fig_plot_energy_2, ax_plot_energy_2 = plt.subplots(figsize=(6, 6))
ax_plot_energy_2.plot(data_frames[0]['ii'], energy_2, 'kx-', label='A')
ax_plot_energy_2.hlines(energy_2.values[-1]+1, 0, 100, 'r', alpha=0.5)
ax_plot_energy_2.hlines(energy_2.values[-1]-1, 0, 100, 'r', alpha=0.5)
ax_plot_energy_2.set_xlim([plot_xlim_2[0], plot_xlim_2[1]])
ax_plot_energy_2.set_ylim([plot_ylim_2[0], plot_ylim_2[1]])
ax_plot_energy_2.set_ylabel('Energy / meV per atom')
fig_plot_energy_2.tight_layout()
fig_plot_energy_2.savefig('{}/energy_2_diff.png'.format(folders[2]), dpi=param.save_dpi)

if __name__ == "__main__":
    print('Finished.')
    plt.show()
