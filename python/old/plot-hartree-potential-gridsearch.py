import numpy as np
import matplotlib.pyplot as plt
from scripts.main import parameters as param

""" Plotting of SMEAGOL output _TRC.agr by filename"""

# CP2K Lithium 1D wire
folder_data = '/Volumes/Storage/Data/delete/Work/Postdoc/Work/calculations/transport/iv/cp2k/transmission/dev/bottom-grid-search-1D-V-0/output'
folder_ref = '/Volumes/Storage/Data/delete/Work/Postdoc/Work/calculations/transport/iv/cp2k/transmission/dev/V-0_HLB-0_z-0-0'
data_em = np.genfromtxt('{}/{}'.format(folder_ref, '0V-VH_z-0001.dat'), skip_header=0, skip_footer=0)
hlb = np.arange(start=-3, stop=3+0.01, step=0.01)
# hlb = np.arange(start=-1, stop=1+0.01, step=0.01)
# hlb = hlb[::20]

# Replace zero
hlb = np.round(hlb, 2)
hlb = np.delete(hlb, np.where(hlb == -0.))

str = []
hlb_str = []
data_hartree = []
for i in range(len(hlb)):
    temp = np.format_float_positional(hlb[i], unique=False, precision=2)
    hlb_str.append(temp)
    str.append('bottom_{}_VH_z.dat'.format(temp))
print(str)

# Hartree potential
fig_plot_1, ax_plot_1 = plt.subplots()
ax_plot_1.legend(frameon=False)
for i in range(len(hlb)):
    data_hartree = np.genfromtxt('{}/{}'.format(folder_data, str[i]), skip_header=0, skip_footer=0)
    ax_plot_1.plot(data_hartree[:, 0], data_hartree[:, 1], label=hlb_str[i])
ax_plot_1.plot(data_em[:, 0], data_em[:, 1], 'k-')
# ax_plot_1.legend(frameon=False)
ax_plot_1.set_xlabel(r'Position / Å')
ax_plot_1.set_ylabel('Hartree potential / eV')
fig_plot_1.tight_layout()
fig_plot_1.savefig('{}/hartree_potential_gridsearch.png'.format(folder_data), dpi=param.save_dpi)
# fig_plot_1.savefig('{}/hartree_potential_gridsearch_fine.png'.format(folder_data), dpi=param.save_dpi)

# Max min difference
fig_plot_2, ax_plot_2 = plt.subplots()
ax_plot_2.legend(frameon=False)
hartree_max = np.zeros(len(hlb))
hartree_min = np.zeros(len(hlb))
for i in range(len(hlb)):
    data_hartree = np.genfromtxt('{}/{}'.format(folder_data, str[i]), skip_header=0, skip_footer=0)
    hartree_max[i] = np.max(data_hartree[:, 1])
    hartree_min[i] = np.min(data_hartree[:, 1])
ax_plot_2.plot(hlb, hartree_max-hartree_min, 'r-', label='EM NEGF')
ax_plot_2.axhline(y=(np.max(data_em[:, 1])-np.min(data_em[:, 1])), color='k', linestyle='-', label='EM DFT')
ax_plot_2.set_xlabel(r'HartreeLeadsBottom / eV')
ax_plot_2.set_ylabel('Hartree potential max-min / eV')
ax_plot_2.legend(frameon=False)
fig_plot_2.tight_layout()
# fig_plot_2.savefig('{}/hartree_potential_difference_gridsearch.png'.format(folder_data), dpi=param.save_dpi)
fig_plot_2.savefig('{}/hartree_potential_difference_gridsearch_fine.png'.format(folder_data), dpi=param.save_dpi)

if __name__ == "__main__":
    print('Finished.')
    plt.show()
