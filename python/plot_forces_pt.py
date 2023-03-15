from __future__ import division, print_function
import pandas as pd
import numpy as np
import glob
from general import load_coordinates
from general import print_xyz
import matplotlib.pyplot as plt
from general import parameters as param

"""
    Plot forces
"""

folder = ['/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Demonstrating/Lab/VolcanoPlot-main/Pt-basis-2x2-input_archer2/lab/surface-DZV/kpoints-1-1-1',
          '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Demonstrating/Lab/VolcanoPlot-main/Pt-basis-2x2-input_archer2/lab/surface-DZV/kpoints-5-5-1',
          '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Demonstrating/Lab/VolcanoPlot-main/Pt-basis-2x2-input_archer2/lab/surface-DZVP-MOLOPT-SR/kpoints-1-1-1',
          '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Demonstrating/Lab/VolcanoPlot-main/Pt-basis-2x2-input_archer2/lab/surface-DZVP-MOLOPT-SR/kpoints-5-5-1',
          '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Demonstrating/Lab/VolcanoPlot-main/Pt-basis-2x2-input_archer2/lab/surface-DZVP-MOLOPT_UZH/kpoints-1-1-1',
          '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Demonstrating/Lab/VolcanoPlot-main/Pt-basis-2x2-input_archer2/lab/surface-DZVP-MOLOPT_UZH/kpoints-5-5-1',
          '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Demonstrating/Lab/VolcanoPlot-main/Pt-basis-2x2-input_archer2/lab/surface-TZVP-MOLOPT-SR_UCL/kpoints-1-1-1',
          '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Demonstrating/Lab/VolcanoPlot-main/Pt-basis-2x2-input_archer2/lab/surface-TZVP-MOLOPT-SR_UCL/kpoints-5-5-1',
          '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Demonstrating/Lab/VolcanoPlot-main/Pt-basis-2x2-input_archer2/lab/surface-TZV2P-MOLOPT-SR_UCL/kpoints-1-1-1',
          '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Demonstrating/Lab/VolcanoPlot-main/Pt-basis-2x2-input_archer2/lab/surface-TZV2P-MOLOPT-SR_UCL/kpoints-5-5-1',
          '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Demonstrating/Lab/VolcanoPlot-main/Pt-basis-2x2-input_archer2/lab/surface-TZV-conf/kpoints-1-1-1',
          '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Demonstrating/Lab/VolcanoPlot-main/Pt-basis-2x2-input_archer2/lab/surface-TZV-conf/kpoints-5-5-1'
          ]
labels = ['DZV-GTH-PADE 1-1-1', 'DZV-GTH-PADE 5-5-1',
          'DZVP-MOLOPT-SR 1-1-1', 'DZVP-MOLOPT-SR 5-5-1',
          'DZVP-MOLOPT 1-1-1', 'DZVP-MOLOPT 5-5-1',
          'TZVP-MOLOPT-SR 1-1-1', 'TZVP-MOLOPT-SR 5-5-1',
          'TZV2P-MOLOPT-SR 1-1-1', 'TZV2P-MOLOPT-SR 5-5-1',
          'TZV-GTH-LDA-q18-very-confined 1-1-1', 'TZV-GTH-LDA-q18-very-confined 5-5-1'
          ]

atoms = np.array([1, 4, 7, 10]) - 1
plotting_colors = ['r', 'r', 'g', 'g', 'b', 'b', 'm', 'm', 'orange', 'orange', 'y', 'y']
input_filename = 'Relax-frc-1.xyz'

# Read number of atoms and labels from .xyz file
cols = ['Species', 'X', 'Y', 'Z']
coord_z = []
for i in range(len(folder)):
    _, _, _, temp, _, _, _ = load_coordinates.load_values_coord(folder[i], input_filename, cols)
    coord_z.append(temp)

# Plot z force for step 0
fig_plot_1, ax_plot_1 = plt.subplots()
plotting = [0, 2, 4, 6, 10]
# plotting = [1, 3, 5, 7]
step = 0
for i in plotting:
    ax_plot_1.plot(atoms + 1, coord_z[i][step][atoms], 'x-', color=plotting_colors[i], label=labels[i])
ax_plot_1.legend(frameon=False)
ax_plot_1.set_xlabel('Atom index')
ax_plot_1.set_ylabel('Force (Hartree/bohr)')
fig_plot_1.tight_layout()
for i in plotting:
    fig_plot_1.savefig('{}/forces_step_{}.png'.format(folder[i], step), dpi=param.save_dpi)

# Plot z force for step n
fig_plot_2, ax_plot_2 = plt.subplots()
step = 3
for i in plotting:
    ax_plot_2.plot(atoms+1, coord_z[i][step][atoms], 'x-', color=plotting_colors[i], label=labels[i])
ax_plot_2.legend(frameon=False)
ax_plot_2.set_xlabel('Atom index')
ax_plot_2.set_ylabel('Force (Hartree/bohr)')
fig_plot_2.tight_layout()
for i in plotting:
    fig_plot_2.savefig('{}/forces_step_{}.png'.format(folder[i], step), dpi=param.save_dpi)

if __name__ == "__main__":
    print('Finished.')
    plt.show()
