import numpy as np
import matplotlib.pyplot as plt
from scripts.main import parameters as param

""" Plotting of SMEAGOL output _TRC.agr by filename"""

plotting_colors = ['r', 'g', 'b', 'm', 'grey', 'orange', 'cyan']

# Lithium chain SIESTA
# xlim = [-1, 1]
# charge = 12
# HartreeLeadsBottom = np.arange(start=-3, stop=0+0.01, step=0.01)
# labels = ['-2 eV', '-1 eV', '1 eV', '2 eV']
# folder = ['/Volumes/Storage/Data/Work/Postdoc/Work/calculations/transport/iv/siesta/bottom-grid-search-1D-V_-2',
#           '/Volumes/Storage/Data/Work/Postdoc/Work/calculations/transport/iv/siesta/bottom-grid-search-1D-V_-1',
#           '/Volumes/Storage/Data/Work/Postdoc/Work/calculations/transport/iv/siesta/bottom-grid-search-1D-V_1',
#           '/Volumes/Storage/Data/Work/Postdoc/Work/calculations/transport/iv/siesta/bottom-grid-search-1D-V_2']
    
# Lithium chain CP2K
xlim = [-1, 1]
charge = 12*3
HartreeLeadsBottom = np.arange(start=-3, stop=3+0.01, step=0.01)
xlim = [-3, -2.8]
ylim = [-0.04, 0.01]
labels = ['V = 0 eV']
folder = ['/Volumes/Storage/Data/delete/Work/Postdoc/Work/calculations/transport/iv/cp2k/transmission/bottom-grid-search-1D-V-0']
data_chr = []
for i in range(len(folder)):
    data_chr.append(np.genfromtxt('{}/charge.log'.format(folder[i]), skip_header=0, skip_footer=0))

# Find numerical solution
val = np.nanargmin(np.abs(data_chr[0]-36))
print(np.abs(data_chr[0])-36)
print('Numerical solution:', HartreeLeadsBottom[val])

# Charge CP2K
fig_plot_1, ax_plot_1 = plt.subplots()
for i in range(len(folder)):
    ax_plot_1.plot(HartreeLeadsBottom, data_chr[i]-charge, '-', color=plotting_colors[i], label=labels[i])
ax_plot_1.axhline(y=0, color='k', linestyle='-')
ax_plot_1.legend(frameon=False)
ax_plot_1.set_xlabel('HartreeLeadsBottom / eV')
ax_plot_1.set_ylabel('Charge difference / e')
fig_plot_1.tight_layout()
for i in range(len(folder)):
    fig_plot_1.savefig('{}/HartreeLeadsBottom_gridsearch_1D.png'.format(folder[i]), dpi=param.save_dpi)

# Charge CP2K
fig_plot_2, ax_plot_2 = plt.subplots()
for i in range(len(folder)):
    ax_plot_2.plot(HartreeLeadsBottom, data_chr[i]-charge, '.', color=plotting_colors[i], label=labels[i])
ax_plot_2.axhline(y=0, color='k', linestyle='-')
ax_plot_2.legend(frameon=False)
ax_plot_2.set_xlabel('HartreeLeadsBottom / eV')
ax_plot_2.set_ylabel('Charge difference / e')
# ax_plot_2.set_xlim([xlim[0], xlim[1]])
# ax_plot_2.set_ylim([ylim[0], ylim[1]])
fig_plot_2.tight_layout()
for i in range(len(folder)):
    fig_plot_2.savefig('{}/HartreeLeadsBottom_gridsearch_1D_2.png'.format(folder[i]), dpi=param.save_dpi)
    
# Charge SIESTA
# fig_plot_2, ax_plot_2 = plt.subplots()
# for i in range(len(folder)):
#     ax_plot_2.plot(HartreeLeadsBottom, data_chr[i][:, 2]-charge, '-', color=plotting_colors[i], label=labels[i])
# ax_plot_2.axhline(y=0, color='k', linestyle='-')
# ax_plot_2.legend(frameon=False)
# ax_plot_2.set_xlabel('HartreeLeadsBottom / eV')
# ax_plot_2.set_ylabel('Charge difference / e')
# fig_plot_2.tight_layout()
# for i in range(len(folder)):
#     fig_plot_2.savefig('{}/HartreeLeadsBottom_gridsearch_1D.png'.format(folder[i]), dpi=param.save_dpi)

if __name__ == "__main__":
    print('Finished.')
    plt.show()
