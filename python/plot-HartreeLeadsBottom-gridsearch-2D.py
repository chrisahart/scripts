import numpy as np
import matplotlib.pyplot as plt
from scripts.main import parameters as param

""" Plotting of SMEAGOL output _TRC.agr by filename"""

# Lithium chain SIESTA
folder_1 = '/Volumes/Storage/Data/Work/Postdoc/Work/calculations/transport/iv/siesta/bottom-grid-search-2D'
data_grid = np.genfromtxt('{}/grid.out'.format(folder_1), skip_header=0, skip_footer=0)
data_charge = np.genfromtxt('{}/charge.log'.format(folder_1), skip_header=0, skip_footer=0)
bias = np.round(np.arange(start=-2, stop=2+0.5, step=0.5), 3)
hartree = np.arange(start=-3, stop=0+0.01, step=0.01)
charge = 12
xlim = [-0.6, 0.01]
ylim = [-0.11, 0.03]

# Ivan Au 4s chain SIESTA coarse
# folder_1 = '/Volumes/Storage/Data/Work/Postdoc/Work/calculations/testing/smeagol/scarf/Smeagol_Tutorial/Smeagol_Tutorial_Files/DFT_NEGF_Transport/Day_1/Inputs/example1/HartreeLeadsBottom/bottom-grid-search-2D-coarse'
# data_grid = np.genfromtxt('{}/grid.out'.format(folder_1), skip_header=0, skip_footer=0)
# data_charge = np.genfromtxt('{}/charge.log'.format(folder_1), skip_header=0, skip_footer=0)
# bias = np.round(np.arange(start=-2, stop=2+0.5, step=0.5), 3)
# hartree = np.arange(start=-3, stop=0+0.01, step=0.01)
# charge = 8

# Ivan Au 4s chain SIESTA
# folder_1 = '/Volumes/Storage/Data/Work/Postdoc/Work/calculations/testing/smeagol/scarf/Smeagol_Tutorial/Smeagol_Tutorial_Files/DFT_NEGF_Transport/Day_1/Inputs/example1/HartreeLeadsBottom/bottom-grid-search-2D-fine'
# data_grid = np.genfromtxt('{}/grid.out'.format(folder_1), skip_header=0, skip_footer=0)
# data_charge = np.genfromtxt('{}/charge.log'.format(folder_1), skip_header=0, skip_footer=0)
# bias = np.round(np.arange(start=-0.5, stop=0.5+0.1, step=0.1), 3)
# hartree = np.arange(start=-2, stop=-1.5+0.01, step=0.01)
# charge = 8

# Find numerical solution
print('bias', bias)
print('hartree', hartree)
index = int(np.where(np.abs(bias) <= 1e-10)[0]) + 1
data_y = np.abs(data_charge[(index - 1) * hartree.shape[0]:index * hartree.shape[0]]-charge)[:, 2]
val = np.nanargmin(data_y)
print('Numerical solution:', hartree[val], 'with charge difference:', data_y[val])

# Charge
fig_plot_1, ax_plot_1 = plt.subplots()
# for i in range(1, bias.shape[0]+1):
# for i in range(3, bias.shape[0]-1):
for i in range(4, bias.shape[0]-2):
    data_x = data_grid[(i - 1) * hartree.shape[0]:i * hartree.shape[0]]
    data_y = data_charge[(i - 1) * hartree.shape[0]:i * hartree.shape[0]]
    ax_plot_1.plot(data_x[:, 1], data_y[:, 2]-charge, '-', label=str(bias[i-1]))
ax_plot_1.axhline(y=0, color='k', linestyle='-')
ax_plot_1.legend(frameon=False)
ax_plot_1.set_xlim([xlim[0], xlim[1]])
ax_plot_1.set_ylim([ylim[0], ylim[1]])
ax_plot_1.set_xlabel('HartreeLeadsBottom / eV')
ax_plot_1.set_ylabel('Charge difference / e')
fig_plot_1.tight_layout()
fig_plot_1.savefig('{}/gridsearch_2D_1.png'.format(folder_1), dpi=param.save_dpi)

# Charge
fig_plot_2, ax_plot_2 = plt.subplots()
for i in range(1, bias.shape[0]+1):
    data_x = data_grid[(i - 1) * hartree.shape[0]:i * hartree.shape[0]]
    data_y = data_charge[(i - 1) * hartree.shape[0]:i * hartree.shape[0]]
    ax_plot_2.plot(data_x[:, 1], data_y[:, 2]-charge, '.-', label=str(bias[i-1]), markersize=5)
ax_plot_2.axhline(y=0, color='k', linestyle='-')
# ax_plot_2.set_xlim([xlim[0], xlim[1]])
# ax_plot_2.set_ylim([ylim[0], ylim[1]])
ax_plot_2.legend(frameon=False)
ax_plot_2.set_xlabel('HartreeLeadsBottom / eV')
ax_plot_2.set_ylabel('Charge difference / e')
fig_plot_2.tight_layout()
fig_plot_2.savefig('{}/gridsearch_2D_2.png'.format(folder_1), dpi=param.save_dpi)

if __name__ == "__main__":
    print('Finished.')
    plt.show()
