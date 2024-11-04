import numpy as np
import matplotlib.pyplot as plt
from scripts.main import parameters as param

""" Plotting of SMEAGOL output _TRC.agr by filename"""

folder_1 = '/Volumes/Storage/Data/Work/Postdoc/Work/calculations/transport/iv/siesta/single-points/v-0_bottom-0'
data_hartree_1_bulk = np.genfromtxt('{}/0.Li.leads-VH_AV.dat'.format(folder_1), skip_header=0, skip_footer=0)
data_hartree_1_em = np.genfromtxt('{}/0.Liwire-VH_AV.dat'.format(folder_1), skip_header=0, skip_footer=0)
print('Min', np.min(data_hartree_1_bulk[:, 1]), 'Average', np.mean(data_hartree_1_bulk[:, 1]),)
print('Left', data_hartree_1_bulk[:, 1][0], 'Right', data_hartree_1_bulk[:, 1][-1])

folder_2 = '/Volumes/Storage/Data/Work/Postdoc/Work/calculations/transport/iv/siesta/single-points/v-0_bottom-0.38727754'
data_hartree_2_bulk = np.genfromtxt('{}/0.Li.leads-VH_AV.dat'.format(folder_2), skip_header=0, skip_footer=0)
data_hartree_2_em = np.genfromtxt('{}/0.Liwire-VH_AV.dat'.format(folder_2), skip_header=0, skip_footer=0)

folder_3 = '/Volumes/Storage/Data/Work/Postdoc/Work/calculations/transport/iv/siesta/single-points/v-0.5_bottom-0'
data_hartree_3_bulk = np.genfromtxt('{}/0.Li.leads-VH_AV.dat'.format(folder_3), skip_header=0, skip_footer=0)
data_hartree_3_em = np.genfromtxt('{}/0.Liwire-VH_AV.dat'.format(folder_3), skip_header=0, skip_footer=0)

folder_4 = '/Volumes/Storage/Data/Work/Postdoc/Work/calculations/transport/iv/siesta/single-points/v-0.5_bottom-0.38727754'
data_hartree_4_bulk = np.genfromtxt('{}/0.Li.leads-VH_AV.dat'.format(folder_4), skip_header=0, skip_footer=0)
data_hartree_4_em = np.genfromtxt('{}/0.Liwire-VH_AV.dat'.format(folder_4), skip_header=0, skip_footer=0)

# Voltage 0 Hartree 0 
fig_plot_1, ax_plot_1 = plt.subplots()
ax_plot_1.plot(data_hartree_1_em[:, 0], data_hartree_1_em[:, 1], 'g-', label='EM')
ax_plot_1.plot(data_hartree_1_bulk[:, 0], data_hartree_1_bulk[:, 1], 'k-', label='Bulk')
ax_plot_1.legend(frameon=False)
ax_plot_1.set_xlabel(r'Position / Å')
ax_plot_1.set_ylabel('Hartree potential / eV')
fig_plot_1.tight_layout()
fig_plot_1.savefig('{}/hartree_potential1.png'.format(folder_1), dpi=param.save_dpi)

# Voltage 0 Hartree -0.38727754
fig_plot_2, ax_plot_2 = plt.subplots()
ax_plot_2.plot(data_hartree_2_em[:, 0], data_hartree_2_em[:, 1], 'g-', label='EM')
ax_plot_2.plot(data_hartree_2_bulk[:, 0], data_hartree_2_bulk[:, 1], 'r-', label='Bulk')
ax_plot_2.legend(frameon=False)
ax_plot_2.set_xlabel(r'Position / Å')
ax_plot_2.set_ylabel('Hartree potential / eV')
fig_plot_2.tight_layout()
fig_plot_2.savefig('{}/hartree_potential2.png'.format(folder_1), dpi=param.save_dpi)

# Voltage 0.1 Hartree 0.0, -0.38727754
fig_plot_3, ax_plot_3 = plt.subplots()
ax_plot_3.plot(data_hartree_3_em[:, 0], data_hartree_3_em[:, 1], 'g-', label='EM HartreeBottomLeads=0.0')
ax_plot_3.plot(data_hartree_4_em[:, 0], data_hartree_4_em[:, 1], 'b-', label='EM HartreeBottomLeads=-0.38727754')
ax_plot_3.plot(data_hartree_4_bulk[:, 0], data_hartree_4_bulk[:, 1], 'y-', label='Bulk')
ax_plot_3.legend(frameon=False)
ax_plot_3.set_xlabel(r'Position / Å')
ax_plot_3.set_ylabel('Hartree potential / eV')
fig_plot_3.tight_layout()
fig_plot_3.savefig('{}/hartree_potential3.png'.format(folder_1), dpi=param.save_dpi)


if __name__ == "__main__":
    print('Finished.')
    plt.show()
