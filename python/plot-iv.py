from __future__ import division, print_function
import pandas as pd
import numpy as np
from general import load_coordinates
from general import parameters as param
import matplotlib.pyplot as plt

""" Plotting of SMEAGOL IV filename"""

# Change matplotlib defaults to large title labels and standard form
params = {'axes.formatter.limits': [-4, 4],
          'axes.labelsize': 'x-large',
          'axes.titlesize': 'x-large',
          'legend.fontsize': 'large',
          'lines.markersize': '8',
          }
plt.rcParams.update(params)

plotting_colors = ['r', 'g', 'b']
n = 0

# cp2k-smeagol-examples/examples/li-chain
labels = ['CP2K+SMEAGOL', 'SIESTA1+SMEAGOL', 'SIESTA3+SMEAGOL']
cp2k_folder1 = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/li-chain/cp2k-smeagol/iv/kpoints-1-1-20'
cp2k1 = np.genfromtxt('{}/IV.log'.format(cp2k_folder1), skip_header=1, skip_footer=0)
siesta_folder1 = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/li-chain/siesta1-smeagol/iv/kpoints-1-1-20'
siesta_folder2 = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/li-chain/siesta3-smeagol/iv/kpoints-2-2-20'
siesta1 = np.genfromtxt('{}/IV.CUR'.format(siesta_folder1), skip_header=0, skip_footer=0)
siesta2 = np.genfromtxt('{}/IV.CUR'.format(siesta_folder2), skip_header=0, skip_footer=0)
plot_cp2k1 = True
plot_cp2k2 = False
plot_siesta1 = True
plot_siesta2 = True
xlim = [-0.5, 0.5]
ylim = np.array([-40, 40])
factor = 1e6
use_xlim = True
use_ylim = True

# IV curve all
fig_plot_2, ax_plot_2 = plt.subplots(figsize=(6, 4))
if plot_cp2k1: ax_plot_2.plot(cp2k1[:, 0], cp2k1[:, 1]* factor, '.-', color=plotting_colors[0], label=labels[0])
if plot_cp2k2: ax_plot_2.plot(cp2k2[:, 0], cp2k2[:, 1]* factor, '.-', color=plotting_colors[1], label=labels[1])
if plot_siesta1: ax_plot_2.plot(siesta1[:, 0], siesta1[:, 1]* factor, '.-', color=plotting_colors[1], label=labels[1])
if plot_siesta2: ax_plot_2.plot(siesta2[:, 0], siesta2[:, 1]* factor, '.-', color=plotting_colors[2], label=labels[2])
if use_xlim: ax_plot_2.set_xlim([xlim[0], xlim[1]])
if use_ylim: ax_plot_2.set_ylim([ylim[0], ylim[1]])
ax_plot_2.legend(frameon=False)
ax_plot_2.set_xlabel('Bias voltage / V')
ax_plot_2.set_ylabel('Current / μA')
fig_plot_2.tight_layout()
if plot_cp2k1: fig_plot_2.savefig('{}/IV.png'.format(cp2k_folder1), dpi=param.save_dpi)
if plot_cp2k2: fig_plot_2.savefig('{}/IV.png'.format(cp2k_folder2), dpi=param.save_dpi)
if plot_siesta1: fig_plot_2.savefig('{}/IV.png'.format(siesta_folder1), dpi=param.save_dpi)
if plot_siesta2: fig_plot_2.savefig('{}/IV.png'.format(siesta_folder2), dpi=param.save_dpi)

if __name__ == "__main__":
    print('Finished.')
    plt.show()
