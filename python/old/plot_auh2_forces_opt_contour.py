from __future__ import division, print_function
import pandas as pd
import numpy as np
import glob
from general import load_coordinates
from general import load_forces
from general import parameters as param
from general import print_xyz
import matplotlib.pyplot as plt
import csv

"""
    Analysis script for au h2 as in clotilde phonon paper, to confirm SMEAGOL forces
"""

params = {'axes.formatter.limits': [-4, 4],
          'axes.labelsize': 'x-large',
          'axes.titlesize': 'x-large',
          'legend.fontsize': 'large',
          'lines.markersize': '8',
          }
plt.rcParams.update(params)


folders1 = ['/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/au-h2/cp2k-smeagol/transmission_all/kpoints-3-3-20_omp-2_ParallelOverKNum-9_contour-single_dynamic-14/0.0',
            '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/au-h2/cp2k-smeagol/transmission_all/kpoints-3-3-20_omp-2_ParallelOverKNum-9_contour-single_dynamic-14/1.5',
            '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/au-h2/cp2k-smeagol/transmission_all/kpoints-3-3-20_omp-2_ParallelOverKNum-9_contour-single_dynamic-14/1.9'
            ]
folders2 = ['/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/au-h2/cp2k-smeagol/transmission_all/kpoints-3-3-20_omp-2_ParallelOverKNum-9_contour-double_dynamic-14/0.0',
            '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/au-h2/cp2k-smeagol/transmission_all/kpoints-3-3-20_omp-2_ParallelOverKNum-9_contour-double_dynamic-14/1.5',
            '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/au-h2/cp2k-smeagol/transmission_all/kpoints-3-3-20_omp-2_ParallelOverKNum-9_contour-double_dynamic-14/1.9'
            ]
files1 = ['0.0V-frc-1.xyz', '1.5V-frc-1.xyz',  '1.9V-frc-1.xyz']
files2 = ['0.0V-frc-1.xyz', '1.5V-frc-1.xyz', '1.9V-frc-1.xyz']
labels = ['0.0V', '1.5V', '1.9V']
diff = 0.2
xlim = [58-diff, 65+diff]
step = 0
dimension = 2
colors = ['r', 'g', 'b', 'm', 'grey', 'orange', 'y', 'brown', 'cyan', 'pink']

# Read position 1
cols = ['Species', 'X', 'Y', 'Z']
coord1 = []
num_atoms1 = []
for i in range(len(folders1)):
    temp_coord, _, _, _, _, temp_num_atoms, _ = load_coordinates.load_values_coord(folders1[i], files1[i], cols)
    coord1.append(temp_coord)
    num_atoms1.append(temp_num_atoms)
index1 = np.arange(start=1, stop=num_atoms1[0]+1)

# Read position 2
cols = ['Species', 'X', 'Y', 'Z']
coord2 = []
num_atoms2 = []
for i in range(len(folders2)):
    temp_coord, _, _, _, _, temp_num_atoms, _ = load_coordinates.load_values_coord(folders2[i], files2[i], cols)
    coord2.append(temp_coord)
    num_atoms2.append(temp_num_atoms)
index2 = np.arange(start=1, stop=num_atoms2[0]+1)

# Plot force component for single timestep all atoms
test = ["Au"] * (58+1) + ["Au2", "Au1", "H1", "H1'", "Au1'", "Au2'"] + ["Au'"] * (113-65)
print(test)
fig_plot_1, ax_plot_1 = plt.subplots()
for i in range(len(folders1)):
    ax_plot_1.plot(np.arange(len(test)), (coord1[i][step, dimension, :]) * param.hartree_per_bohr_to_ev_per_angstrom, 'x--',
                   color=colors[i], fillstyle='none')
for i in range(len(folders2)):
    ax_plot_1.plot(np.arange(len(test)), (coord2[i][step, dimension, :]) * param.hartree_per_bohr_to_ev_per_angstrom, 'o-',
                   color=colors[i], fillstyle='none', label=labels[i])


plt.xticks(np.arange(len(test)), test, fontsize=13)
# ax_plot_1.set_xlabel('Atom index')
ax_plot_1.set_ylabel('Force z (eV / Å)')
ax_plot_1.set_xlim([xlim[0], xlim[1]])
ax_plot_1.legend(frameon=False)
fig_plot_1.tight_layout()
for i in range(len(folders1)):
    fig_plot_1.savefig('{}/force.png'.format(folders1[i]), dpi=param.save_dpi)

# Plot z force for single atom all timesteps
# atom = np.array([61]) - 1
# atoms = np.array([61, 62, 63]) - 1
# xlim = [1-diff, 18+diff]
# fig_plot_2, ax_plot_2 = plt.subplots()
# ax_plot_2.plot(dft_time, (coord_dft[:, 0, atom]) * param.hartree_per_bohr_to_ev_per_angstrom, 'ko-', fillstyle='none', label=labels[0])
# ax_plot_2.plot(V0_time, (coord_V0[:, 0, atom]) * param.hartree_per_bohr_to_ev_per_angstrom, 'ro-', fillstyle='none', label=labels[1])
# ax_plot_2.plot(dft_time, (coord_dft[:, 0, atom+1]) * param.hartree_per_bohr_to_ev_per_angstrom, 'ko--', fillstyle='none')
# ax_plot_2.plot(V0_time, (coord_V0[:, 0, atom+1]) * param.hartree_per_bohr_to_ev_per_angstrom, 'ro--', fillstyle='none')
# ax_plot_2.set_xlabel('GEO_OPT step')
# ax_plot_2.set_ylabel('Force x (eV / Å)')
# ax_plot_2.set_xlim([xlim[0], xlim[1]])
# ax_plot_2.legend(frameon=False)
# fig_plot_2.tight_layout()
# for i in range(len(folders1)):
#     fig_plot_2.savefig('{}/force_x_time.png'.format(folders1[i]), dpi=param.save_dpi)

# Plot all forces on atom
# rows, cols = 3, 1
# fig_plot_forces_all, ax_plot_forces_all = plt.subplots(rows, cols, sharex='col', sharey='row', figsize=(6, 8))
# ax_plot_forces_all[0].plot(dft_time, (coord_dft[:, 0, atom]) * param.hartree_per_bohr_to_ev_per_angstrom, 'ko-', fillstyle='none', label='CP2K')
# ax_plot_forces_all[0].plot(V0_time, (coord_V0[:, 0, atom]) * param.hartree_per_bohr_to_ev_per_angstrom, 'ro-', fillstyle='none', label='CP2K-SMEAGOL')
# ax_plot_forces_all[1].plot(dft_time, (coord_dft[:, 1, atom]) * param.hartree_per_bohr_to_ev_per_angstrom, 'ko-', fillstyle='none', label='CP2K')
# ax_plot_forces_all[1].plot(V0_time, (coord_V0[:, 1, atom]) * param.hartree_per_bohr_to_ev_per_angstrom, 'ro-', fillstyle='none', label='CP2K-SMEAGOL')
# ax_plot_forces_all[2].plot(dft_time, (coord_dft[:, 2, atom]) * param.hartree_per_bohr_to_ev_per_angstrom, 'ko-', fillstyle='none', label='CP2K')
# ax_plot_forces_all[2].plot(V0_time, (coord_V0[:, 2, atom]) * param.hartree_per_bohr_to_ev_per_angstrom, 'ro-', fillstyle='none', label='CP2K-SMEAGOL')
# ax_plot_forces_all[0].set_ylabel('Force x (eV / Å)')
# ax_plot_forces_all[1].set_ylabel('Force y (eV / Å)')
# ax_plot_forces_all[2].set_ylabel('Force z (eV / Å)')
# ax_plot_forces_all[2].set_xlabel('GEO_OPT step')
# ax_plot_forces_all[0].set_xlim([xlim[0], xlim[1]])
# ax_plot_forces_all[1].set_xlim([xlim[0], xlim[1]])
# ax_plot_forces_all[2].set_xlim([xlim[0], xlim[1]])
# ax_plot_forces_all[0].legend(frameon=False)
# fig_plot_forces_all.tight_layout()
# for i in range(len(folders1)):
#     fig_plot_forces_all.savefig('{}/force_all_time.png'.format(folders1[i]), dpi=param.save_dpi)

# Plot all x force on atoms
# rows, cols = 3, 1
# fig_plot_atoms_all, ax_plot_atoms_all = plt.subplots(rows, cols, sharex='col', sharey='row', figsize=(6, 8))
# ax_plot_atoms_all[0].plot(dft_time, (coord_dft[:, 0, atom-1]) * param.hartree_per_bohr_to_ev_per_angstrom, 'ko-', fillstyle='none', label='CP2K')
# ax_plot_atoms_all[0].plot(V0_time, (coord_V0[:, 0, atom-1]) * param.hartree_per_bohr_to_ev_per_angstrom, 'ro-', fillstyle='none', label='CP2K-SMEAGOL')
# ax_plot_atoms_all[1].plot(dft_time, (coord_dft[:, 0, atom]) * param.hartree_per_bohr_to_ev_per_angstrom, 'ko-', fillstyle='none', label='CP2K')
# ax_plot_atoms_all[1].plot(V0_time, (coord_V0[:, 0, atom]) * param.hartree_per_bohr_to_ev_per_angstrom, 'ro-', fillstyle='none', label='CP2K-SMEAGOL')
# ax_plot_atoms_all[2].plot(dft_time, (coord_dft[:, 0, atom-2]) * param.hartree_per_bohr_to_ev_per_angstrom, 'ko-', fillstyle='none', label='CP2K')
# ax_plot_atoms_all[2].plot(V0_time, (coord_V0[:, 0, atom-2]) * param.hartree_per_bohr_to_ev_per_angstrom, 'ro-', fillstyle='none', label='CP2K-SMEAGOL')
# ax_plot_atoms_all[0].set_ylabel('Force x atom 4 (eV / Å)')
# ax_plot_atoms_all[1].set_ylabel('Force x atom 5 (eV / Å)')
# ax_plot_atoms_all[2].set_ylabel('Force x atom 6 (eV / Å)')
# ax_plot_atoms_all[2].set_xlabel('GEO_OPT step')
# ax_plot_atoms_all[0].set_xlim([xlim[0], xlim[1]])
# ax_plot_atoms_all[1].set_xlim([xlim[0], xlim[1]])
# ax_plot_atoms_all[2].set_xlim([xlim[0], xlim[1]])
# ax_plot_atoms_all[0].legend(frameon=False)
# fig_plot_atoms_all.tight_layout()
# for i in range(len(folders1)):
#     fig_plot_atoms_all.savefig('{}/force_all_time.png'.format(folders1[i]), dpi=param.save_dpi)
    
if __name__ == "__main__":
    print('Finished.')
    plt.show()
