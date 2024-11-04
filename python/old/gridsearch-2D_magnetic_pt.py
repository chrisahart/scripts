from __future__ import division, print_function
import pandas as pd
import numpy as np
from general import parameters as param
import matplotlib.pyplot as plt

"""
    Sort .xyz sequentially along desired axes (useful for transport calculations)
"""

folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Masters/2021-2022/Yike/Ni/ni/archer2/slab/optimise/magnetisation/scf_1e-4_temp-2000_mos-1000_alpha-0.08_beta-1.5_periodic-xyz'
file_input = '{}/data.out'.format(folder)
file_output_alpha_beta_contour_1 = '{}/alpha_beta_energy_contour_1.png'.format(folder)
file_output_alpha_beta_scatter_0 = '{}/alpha_beta_energy_scatter_0.png'.format(folder)
file_output_alpha_beta_scatter_1 = '{}/alpha_beta_energy_scatter_1.png'.format(folder)
file_output_alpha_beta_scatter_2 = '{}/alpha_beta_energy_scatter_2.png'.format(folder)
file_output_alpha_beta_scatter_3 = '{}/alpha_beta_energy_scatter_3.png'.format(folder)
save = True

cols = ['ii', 'jj', 'Energy', 'SCF', 'time', 'iasd']
file_coord = pd.read_csv(file_input, names=cols, delim_whitespace=True)

# Alpha beta
# plot_log = True
# x_label = 'Alpha'
# y_label = 'Beta'
# x_units = 1
# y_units = 1
# xlim = [1e-3-1e-4, 1e-1+1e-2]
# ylim = [1e-3-1e-4, 5e-1+1e-2]
# print(file_coord)

# Magnetisation
plot_log = False
x_label = 'Ni_b'
y_label = 'Ni_i'
x_units = 1
y_units = 1
file_coord['ii'] = file_coord['ii'] * x_units
file_coord['jj'] = file_coord['jj'] * y_units
xlim = np.array([np.min(file_coord['ii'])-0.05, np.max(file_coord['ii'])+0.05])
ylim = np.array([np.min(file_coord['jj'])-0.01, np.max(file_coord['jj'])+0.01])
print(file_coord)

# Plot ii, jj
fig_plot_scatter_0, ax_plot_scatter_0 = plt.subplots()
sc = ax_plot_scatter_0.scatter(file_coord['ii'], file_coord['jj'], cmap="copper")
ax_plot_scatter_0.set_xlim([xlim[0], xlim[1]])
ax_plot_scatter_0.set_ylim([ylim[0], ylim[1]])
if plot_log: ax_plot_scatter_0.set_xscale('log')
if plot_log: ax_plot_scatter_0.set_yscale('log')
ax_plot_scatter_0.set_xlabel(x_label)
ax_plot_scatter_0.set_ylabel(y_label)
if save: fig_plot_scatter_0.savefig(file_output_alpha_beta_scatter_0, dpi=param.save_dpi)

# Delete rows containing SCF=300 or SCF=0
# file_coord = file_coord.drop(file_coord[file_coord.SCF == 300].index)
file_coord = file_coord.drop(file_coord[file_coord.SCF == 0].index)
file_coord = file_coord.drop(file_coord[file_coord.SCF > 50].index)
file_coord = file_coord.drop(file_coord[file_coord.iasd < 23].index)
file_coord = file_coord.reset_index(drop=True)
# file_coord['Energy'] = file_coord['Energy'] - np.min(file_coord['Energy'])
print(file_coord)

# Plot ii, jj, energy as a colored scatter plot
# fig_plot_scatter_1, ax_plot_scatter_1 = plt.subplots()
# sc = ax_plot_scatter_1.scatter(file_coord['ii'], file_coord['jj'], c=file_coord['Energy'],  cmap="copper")
# cbar = fig_plot_scatter_1.colorbar(sc, ax=ax_plot_scatter_1)
# cbar.set_label('Energy difference / Ha')
# ax_plot_scatter_1.set_xlim([xlim[0], xlim[1]])
# ax_plot_scatter_1.set_ylim([ylim[0], ylim[1]])
# if plot_log: ax_plot_scatter_1.set_xscale('log')
# if plot_log: ax_plot_scatter_1.set_yscale('log')
# ax_plot_scatter_1.set_xlabel(x_label)
# ax_plot_scatter_1.set_ylabel(y_label)
# # fig_plot_scatter_1.tight_layout()
# if save: fig_plot_scatter_1.savefig(file_output_alpha_beta_scatter_1, dpi=param.save_dpi)

# Plot ii, jj, IASD as a colored scatter plot
fig_plot_scatter_2, ax_plot_scatter_2 = plt.subplots()
sc = ax_plot_scatter_2.scatter(file_coord['ii'], file_coord['jj'], c=file_coord['iasd'],  cmap="copper")
cbar = fig_plot_scatter_2.colorbar(sc, ax=ax_plot_scatter_2)
cbar.set_label('IASD')
ax_plot_scatter_2.set_xlim([xlim[0], xlim[1]])
ax_plot_scatter_2.set_ylim([ylim[0], ylim[1]])
if plot_log: ax_plot_scatter_2.set_xscale('log')
if plot_log: ax_plot_scatter_2.set_yscale('log')
ax_plot_scatter_2.set_xlabel(x_label)
ax_plot_scatter_2.set_ylabel(y_label)
fig_plot_scatter_2.tight_layout()
if save: fig_plot_scatter_2.savefig(file_output_alpha_beta_scatter_2, dpi=param.save_dpi)

# Plot ii, jj, IASD as a colored scatter plot
fig_plot_scatter_3, ax_plot_scatter_3 = plt.subplots()
sc = ax_plot_scatter_3.scatter(file_coord['ii'], file_coord['jj'], c=file_coord['SCF'],  cmap="copper")
cbar = fig_plot_scatter_3.colorbar(sc, ax=ax_plot_scatter_3)
cbar.set_label('SCF steps')
ax_plot_scatter_3.set_xlim([xlim[0], xlim[1]])
ax_plot_scatter_3.set_ylim([ylim[0], ylim[1]])
if plot_log: ax_plot_scatter_3.set_xscale('log')
if plot_log: ax_plot_scatter_3.set_yscale('log')
ax_plot_scatter_3.set_xlabel(x_label)
ax_plot_scatter_3.set_ylabel(y_label)
fig_plot_scatter_3.tight_layout()
if save: fig_plot_scatter_3.savefig(file_output_alpha_beta_scatter_3, dpi=param.save_dpi)

# Plot ii, jj, energy as a contour plot
# fig_plot_contour_1, ax_plot_contour_1 = plt.subplots()
# cntr2 = ax_plot_contour_1.tricontourf(file_coord['ii'], file_coord['jj'], file_coord['Energy'],  cmap="copper")
# fig_plot_contour_1.colorbar(cntr2, ax=ax_plot_contour_1)
# ax_plot_contour_1.set_xlim([xlim[0], xlim[1]])
# ax_plot_contour_1.set_ylim([ylim[0], ylim[1]])
# ax_plot_contour_1.set_xscale('log')
# ax_plot_contour_1.set_yscale('log')
# ax_plot_contour_1.set_xlabel(x_label)
# ax_plot_contour_1.set_ylabel(y_label)
# fig_plot_contour_1.tight_layout()
# fig_plot_contour_1.savefig(file_output_alpha_beta_contour_1, dpi=param.save_dpi)

if __name__ == "__main__":
    print('Finished.')
    plt.show()
