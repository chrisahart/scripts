import numpy as np
import matplotlib.pyplot as plt
from scripts.main import parameters as param

""" Plotting of SMEAGOL output _TRC.agr by filename"""

plotting_colors = ['r', 'g', 'b', 'm', 'grey']
n = 1

xlim = [-1, 1]  # Li and Au chain
# ylim = [0, 8e-5]

# labels = ['CP2K Li', 'SIESTA Li']
# cp2k_folder = '/Volumes/Storage/Data/Work/Postdoc/Work/calculations/transport/iv/cp2k/12-atoms-a9b-tidy-iv'
# siesta_folder = '/Volumes/Storage/Data/Work/Postdoc/Work/calculations/transport/iv/siesta/0V-tidy-iv'
# cp2k = np.genfromtxt('{}/IV.log'.format(cp2k_folder), skip_header=1, skip_footer=0)
# siesta = np.genfromtxt('{}/Liwire.CUR'.format(siesta_folder), skip_header=0, skip_footer=0)

labels = ['CP2K Li', 'SIESTA Li', 'CP2K Li align_Vhartree=F']
cp2k_folder1 = '/Volumes/Storage/Data/Work/Postdoc/Work/calculations/transport/iv/cp2k/12-atoms-a9b-tidy-iv'
cp2k1 = np.genfromtxt('{}/IV.log'.format(cp2k_folder1), skip_header=1, skip_footer=0)
cp2k_folder2 = '/Volumes/Storage/Data/Work/Postdoc/Work/calculations/transport/iv/cp2k/iv-bottom-0-VHartreeF'
cp2k2 = np.genfromtxt('{}/IV.log'.format(cp2k_folder2), skip_header=1, skip_footer=0)
siesta_folder = '/Volumes/Storage/Data/Work/Postdoc/Work/calculations/transport/iv/siesta/iv/iv-bottom-0.38727754-left-0-right-0'
siesta = np.genfromtxt('{}/Liwire.CUR'.format(siesta_folder), skip_header=0, skip_footer=0)

# IV curve CP2K vs SIESTA
fig_plot_1, ax_plot_1 = plt.subplots()
ax_plot_1.plot(cp2k1[:, 0], cp2k1[:, 1], '.-', color=plotting_colors[0], label=labels[0])
ax_plot_1.plot(cp2k2[:, 0], cp2k2[:, 1], '.-', color=plotting_colors[2], label=labels[2])
ax_plot_1.plot(siesta[:, 0], siesta[:, 1], '.-', color=plotting_colors[1], label=labels[1])
# ax_plot_1.set_xlim([xlim[0], xlim[1]])
# ax_plot_1.set_ylim([ylim[0], ylim[1]])
ax_plot_1.legend(frameon=False)
ax_plot_1.set_xlabel('Bias voltage / eV')
ax_plot_1.set_ylabel('Current / A')
fig_plot_1.tight_layout()
fig_plot_1.savefig('{}/IV.png'.format(cp2k_folder1), dpi=param.save_dpi)
# fig_plot_1.savefig('{}/IV.png'.format(siesta_folder), dpi=param.save_dpi)


if __name__ == "__main__":
    print('Finished.')
    plt.show()
