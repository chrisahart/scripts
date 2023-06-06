import numpy as np
from matplotlib import pyplot as plt
from ase.io.cube import read_cube_data
from general import parameters as param
import pandas as pd

folder_cp2k = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/stm/siesta'
size_siesta = np.array([120, 256])
siesta_scale = size_siesta
siesta_dft = np.array([1092862.351, 2286639.663]) / 5
siesta_V0 = np.array([455430.109, 660019.751]) / 5
siesta_V = np.array([709182.198, 895608.440]) / 5

# Plot DFT only
fig_plot_dft_1, ax_plot_dft_1 = plt.subplots()
ax_plot_dft_1.plot(size_siesta, siesta_dft / siesta_scale, 'kx-', label='SIESTA')
ax_plot_dft_1.plot(size_siesta, siesta_V0 / siesta_scale, 'gx-', label='SIESTA-SMEAGOL V=0')
ax_plot_dft_1.plot(size_siesta, siesta_V / siesta_scale, 'bx-', label='SIESTA-SMEAGOL V=0.2')
ax_plot_dft_1.set_xlabel('Number of cores')
ax_plot_dft_1.set_ylabel('Time taken per scf step / s')
ax_plot_dft_1.legend(frameon=False)
fig_plot_dft_1.tight_layout()
fig_plot_dft_1.savefig('{}/time_atoms_dft.png'.format(folder_cp2k), dpi=param.save_dpi)

# Plot all
# fig_plot_all_1, ax_plot_all_1 = plt.subplots()
# ax_plot_all_1.plot(size_cp2k, cp2k_1/ cp2k_scale, 'kx-', label=labels_cp2k_1[0])
# ax_plot_all_1.plot(size_cp2k, cp2k_2/ cp2k_scale, 'gx-', label=labels_cp2k_1[1])
# ax_plot_all_1.plot(size_cp2k, cp2k_3/ cp2k_scale, 'bx-', label=labels_cp2k_1[2])
# ax_plot_all_1.plot(size_siesta, siesta_1/siesta_scale, 'kx--', label=labels_siesta[0])
# ax_plot_all_1.plot(size_siesta, siesta_2/siesta_scale, 'gx--', label=labels_siesta[1])
# ax_plot_all_1.plot(size_siesta, siesta_3/siesta_scale, 'bx--', label=labels_siesta[2])
# ax_plot_all_1.set_xlabel('Number atoms')
# ax_plot_all_1.set_ylabel('Time taken per scf step / s')
# ax_plot_all_1.legend(frameon=False)
# fig_plot_all_1.tight_layout()
# fig_plot_all_1.savefig('{}/time_atoms_all.png'.format(folder_cp2k), dpi=param.save_dpi)
# fig_plot_all_1.savefig('{}/time_atoms_all.png'.format(folder_siesta), dpi=param.save_dpi)

if __name__ == "__main__":
    print('Finished.')
    plt.show()
