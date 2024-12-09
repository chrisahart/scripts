from __future__ import division, print_function
import pandas as pd
import numpy as np
from general import parameters as param
import matplotlib.pyplot as plt

"""
    Sort .xyz sequentially along desired axes (useful for transport calculations)
"""

folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/supercell-1-1-5-bulk-6/junction/dft/energy/NIMAGES_IJ-7-kpoints-2-2-SZV-PBE-CUSTOM-q1-f3/cu-o-test'
folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/supercell-1-1-5-bulk-6/junction/dft/energy/NIMAGES_IJ-7-kpoints-2-2-SZV-PBE-CUSTOM-q1-f3/cu-o-coarse'
folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/supercell-1-1-5-bulk-6/junction/dft/energy/NIMAGES_IJ-7-kpoints-2-2-SZV-PBE-CUSTOM-q1-f3/cu-o-fine'
file_input = '{}/data.out'.format(folder)
file_plot_energy = '{}/plot_energy.png'.format(folder)
file_plot_scf = '{}/plot_scf.png'.format(folder)
file_plot_time = '{}/plot_time.png'.format(folder)
x_label = 'k-points'
save_plot = True

cols = ['ii', 'energy', 'scf', 'time']
file_coord = pd.read_csv(file_input, names=cols, delim_whitespace=True)

# 1.790788 Vacuum optimised Cu-O z GEO_OPT only

# Delete rows SCF=0
file_coord = file_coord.drop(file_coord[file_coord.scf == 0].index)
file_coord = file_coord.reset_index(drop=True)
print(file_coord)
file_coord = file_coord.sort_values(by='ii')
file_coord = file_coord.reset_index(drop=True)
print(file_coord)

length_cu_o = []
for i in range(file_coord['ii'].shape[0]):
    print(file_coord['ii'][i])

    # Smaller
    atom_cu = np.array([2.55200, 2.80200, 21.65900])
    atom_o = np.array([2.55200, 2.59500, 23.65900 - file_coord['ii'][i]])
    bond_length_1 = np.sqrt((atom_cu[0] - atom_o[0]) ** 2 +
                          (atom_cu[1] - atom_o[1]) ** 2 +
                          (atom_cu[2] - atom_o[2]) ** 2)
    print(bond_length_1)

    # Larger
    atom_cu = np.array([0.00000, 2.80200, 21.65900])
    atom_o = np.array([0.02700, 3.05200,  23.65900 - file_coord['ii'][i]])
    bond_length_2 = np.sqrt((atom_cu[0] - atom_o[0]) ** 2 +
                          (atom_cu[1] - atom_o[1]) ** 2 +
                          (atom_cu[2] - atom_o[2]) ** 2)
    print(bond_length_2)

    bond_length = (bond_length_1 + bond_length_2) / 2
    length_cu_o.append(bond_length)

print(length_cu_o)
plot_xlim = [np.min(length_cu_o), np.max(length_cu_o)]
plot_ylim = [-0.005, 0.11]

# Plot ii vs energy
fig_plot_energy, ax_plot_energy = plt.subplots()
energy = (file_coord['energy'] - np.min(file_coord['energy'])) * param.hartree_to_ev
# energy = file_coord['energy']
fig_plot_energy, ax_plot_energy = plt.subplots(figsize=(6, 6))
ax_plot_energy.plot(length_cu_o, energy, 'kx-')
ax_plot_energy.vlines(1.790788, -1e3, 1e3, 'r', alpha=0.5)
# ax_plot_energy.hlines(1, 0, 100, 'r', alpha=0.5)
# ax_plot_energy.hlines(-1, 0, 100, 'r', alpha=0.5)
ax_plot_energy.set_xlabel('Cu-O / A')
ax_plot_energy.set_ylabel('Energy / eV')
ax_plot_energy.set_xlim([plot_xlim[0], plot_xlim[1]])
ax_plot_energy.set_ylim([plot_ylim[0], plot_ylim[1]])
fig_plot_energy.tight_layout()
if save_plot:
    fig_plot_energy.savefig(file_plot_energy, dpi=param.save_dpi)


opt = 0.156
# Smaller
atom_cu = np.array([2.55200, 2.80200, 21.65900])
atom_o = np.array([2.55200, 2.59500, 23.65900 - opt])
bond_length_1 = np.sqrt((atom_cu[0] - atom_o[0]) ** 2 +
                        (atom_cu[1] - atom_o[1]) ** 2 +
                        (atom_cu[2] - atom_o[2]) ** 2)
print(bond_length_1)

# Larger
atom_cu = np.array([0.00000, 2.80200, 21.65900])
atom_o = np.array([0.02700, 3.05200, 23.65900 - opt])
bond_length_2 = np.sqrt((atom_cu[0] - atom_o[0]) ** 2 +
                        (atom_cu[1] - atom_o[1]) ** 2 +
                        (atom_cu[2] - atom_o[2]) ** 2)
print(bond_length_2)

bond_length = (bond_length_1 + bond_length_2) / 2
print(bond_length)

if __name__ == "__main__":
    print('Finished.')
    plt.show()
