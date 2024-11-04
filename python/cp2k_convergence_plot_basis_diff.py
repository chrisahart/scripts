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

params = {'axes.formatter.limits': [-4, 4],
          'axes.labelsize': 'x-large',
          'xtick.labelsize' : '12',
          'axes.titlesize': 'x-large',
          'legend.fontsize': 'medium',
          'lines.markersize': '10',
          }
plt.rcParams.update(params)


folders = [
    '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/bulk/monoclinic/convergence/monoclinic-cell-ref-large-rs-equal-cell/MOLOPT-SR-GTH_GTH_POTENTIALS',
    '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/bulk/tetragonal/convergence/yudi-pbe-cell-ref-large-rs-equal-cell/MOLOPT-SR-GTH_GTH_POTENTIALS',
    '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/bulk/po/convergence/yudi-pbe-cell-ref-large-rs-equal-cell/MOLOPT-SR-GTH_GTH_POTENTIALS'
]

folders = [
    '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/bulk/monoclinic/convergence/monoclinic-cell-ref-large-rs-equal-cell/MOLOPT-GTH_POTENTIAL_UZH',
    '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/bulk/tetragonal/convergence/yudi-pbe-cell-ref-large-rs-equal-cell/MOLOPT-GTH_POTENTIAL_UZH',
    '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/bulk/po/convergence/yudi-pbe-cell-ref-large-rs-equal-cell/MOLOPT-GTH_POTENTIAL_UZH'
]

folders = [
    '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/bulk/monoclinic/convergence/monoclinic-cell-ref-large-rs-equal-cell/MOLOPT_GAPW_POTENTIAL_UZH',
    '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/bulk/tetragonal/convergence/yudi-pbe-cell-ref-large-rs-equal-cell/MOLOPT_GAPW_POTENTIAL_UZH',
    '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/bulk/po/convergence/yudi-pbe-cell-ref-large-rs-equal-cell/MOLOPT_GAPW_POTENTIAL_UZH'
]

# Process each folder
data_frames = []
for folder in folders:
    file_input = f'{folder}/data.out'
    cols = ['ii', 'energy', 'scf', 'time']
    df = pd.read_csv(file_input, names=cols, delim_whitespace=True)
    data_frames.append(df)

x_label = ['DZP', 'TZP', 'TZ2P']
x_label = ['TV2P', 'QZ2P']

energy_1 = np.abs((data_frames[0]['energy']-data_frames[1]['energy'])*param.hartree_to_ev)
energy_2 = np.abs((data_frames[0]['energy']-data_frames[2]['energy'])*param.hartree_to_ev)
print(energy_1)
print(energy_2)

# Plot ii vs energy 1
fig_plot_energy_1, ax_plot_energy_1 = plt.subplots(figsize=(6, 6))
ax_plot_energy_1.plot(x_label, energy_1, 'kx-', label='A')
ax_plot_energy_1.set_ylabel('Energy / eV')
fig_plot_energy_1.tight_layout()
fig_plot_energy_1.savefig('{}/energy_1_diff.png'.format(folders[2]), dpi=param.save_dpi)

# Plot ii vs energy
fig_plot_energy_2, ax_plot_energy_2 = plt.subplots(figsize=(6, 6))
ax_plot_energy_2.plot(x_label, energy_2, 'kx-', label='A')
ax_plot_energy_2.set_ylabel('Energy / eV')
fig_plot_energy_2.tight_layout()
fig_plot_energy_2.savefig('{}/energy_2_diff.png'.format(folders[2]), dpi=param.save_dpi)

if __name__ == "__main__":
    print('Finished.')
    plt.show()
