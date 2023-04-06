import numpy as np
import matplotlib.pyplot as plt
from general import parameters as param
from ase.io.cube import read_cube_data
import copy
import pandas as pd

""" Plot Mulliken and Hirshfeld charges"""


def read_mulliken(filename):
    """
    Read Mulliken analysis from CP2K output file (requires removal from CP2K output file to filename)
    """

    # Read number of atoms and labels from .xyz file
    # cols = ['Atom', 'Element', 'Kind', 'Ref Charge', 'Pop 1', 'Pop 2', 'Spin', 'Charge']
    cols = ['Atom', 'Element', 'Kind', 'Population', 'Charge']
    file_spec1 = pd.read_csv(filename, names=cols, delim_whitespace=True, skiprows=3)
    species = file_spec1['Element']

    Au_index = [i for i, e in enumerate(species) if e == 'Au']
    H_index = [i for i, e in enumerate(species) if e == 'H']
    O_index = [i for i, e in enumerate(species) if e == 'O']
    Fe_index = [i for i, e in enumerate(species) if e == 'Fe_a' or e == 'Fe_b' or e == 'Fe']

    # Force database to numeric, assigning any non-numeric as NaN
    file_spec1 = file_spec1.apply(pd.to_numeric, errors='coerce')

    Au_db1 = file_spec1.loc[Au_index]
    Fe_db1 = file_spec1.loc[Fe_index]
    O_db1 = file_spec1.loc[O_index]
    H_db1 = file_spec1.loc[H_index]

    return Au_db1, Fe_db1, O_db1, H_db1, file_spec1


def read_mulliken_siesta(filename):
    """
    Read Mulliken analysis from CP2K output file (requires removal from CP2K output file to filename)
    """

    # Read number of atoms and labels from .xyz file
    # cols = ['Atom', 'Element', 'Kind', 'Ref Charge', 'Pop 1', 'Pop 2', 'Spin', 'Charge']
    cols = ['Index', 'Total', 'A', 'B', 'C', 'D', 'E', 'F']
    file_spec1 = pd.read_csv(filename, names=cols, delim_whitespace=True, skiprows=3)

    return file_spec1

def read_hirsh(filename):
    """
    Read Hirshfeld analysis from CP2K output file (requires removal from CP2K output file to filename)
    """

    # Read number of atoms and labels from .xyz file
    # cols = ['Atom', 'Element', 'Kind', 'Ref Charge', 'Pop 1', 'Pop 2', 'Spin', 'Charge']
    cols = ['Atom', 'Element', 'Kind', 'Ref Charge', 'Pop', 'Charge']
    file_spec1 = pd.read_csv(filename, names=cols, delim_whitespace=True, skiprows=3)
    species = file_spec1['Element']

    Au_index = [i for i, e in enumerate(species) if e == 'Au']
    H_index = [i for i, e in enumerate(species) if e == 'H']
    O_index = [i for i, e in enumerate(species) if e == 'O']
    Fe_index = [i for i, e in enumerate(species) if e == 'Fe_a' or e == 'Fe_b' or e == 'Fe']

    # Force database to numeric, assigning any non-numeric as NaN
    file_spec1 = file_spec1.apply(pd.to_numeric, errors='coerce')

    Au_db1 = file_spec1.loc[Au_index]
    Fe_db1 = file_spec1.loc[Fe_index]
    O_db1 = file_spec1.loc[O_index]
    H_db1 = file_spec1.loc[H_index]

    return Au_db1, Fe_db1, O_db1, H_db1, file_spec1


# folder_cp2k = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/cp2k-smeagol/sz/transmission/exp/capacitor/sergey-equal/auto/kpoints-100/V-0_kpoints_bulk-1-1-100_em-1-1-1_hlb-auto'
# filename_hirsh_0 = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/cp2k-smeagol/sz/transmission/exp/capacitor/sergey-equal/auto/kpoints-100/V-0_kpoints_bulk-1-1-100_em-1-1-1_hlb-auto/hirsh_log_2_dft_wfn.out'
# filename_hirsh_1 = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/cp2k-smeagol/sz/transmission/exp/capacitor/sergey-equal/auto/kpoints-100/V-0_kpoints_bulk-1-1-100_em-1-1-1_hlb-auto/hirsh_log_3_0V.out'
# filename_hirsh_2 = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/cp2k-smeagol/sz/transmission/exp/capacitor/sergey-equal/auto/kpoints-100/V-1_kpoints_bulk-1-1-100_em-1-1-1_hlb-auto/hirsh_log_3_0V.out'
# filename_mulliken_0 = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/cp2k-smeagol/sz/transmission/exp/capacitor/sergey-equal/auto/kpoints-100/V-0_kpoints_bulk-1-1-100_em-1-1-1_hlb-auto/mulliken_log_2_dft_wfn.out'
# filename_mulliken_1 = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/cp2k-smeagol/sz/transmission/exp/capacitor/sergey-equal/auto/kpoints-100/V-0_kpoints_bulk-1-1-100_em-1-1-1_hlb-auto/mulliken_log_3_0V.out'
# filename_mulliken_2 = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/cp2k-smeagol/sz/transmission/exp/capacitor/sergey-equal/auto/kpoints-100/V-1_kpoints_bulk-1-1-100_em-1-1-1_hlb-auto/mulliken_log_3_0V.out'

folder_cp2k = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/cp2k-smeagol/sz/transmission/exp/capacitor/sergey-equal/auto/kpoints-100/V-0_kpoints_bulk-2-2-100_em-2-2-1_hlb-auto'
filename_hirsh_0 = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/cp2k-smeagol/sz/transmission/exp/capacitor/sergey-equal/auto/kpoints-100/V-0_kpoints_bulk-2-2-100_em-2-2-1_hlb-auto/hirsh_log_2_dft_wfn.out'
filename_hirsh_1 = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/cp2k-smeagol/sz/transmission/exp/capacitor/sergey-equal/auto/kpoints-100/V-0_kpoints_bulk-2-2-100_em-2-2-1_hlb-auto/hirsh_log_3_0V.out'
filename_hirsh_2 = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/cp2k-smeagol/sz/transmission/exp/capacitor/sergey-equal/auto/kpoints-100/V-1_kpoints_bulk-2-2-100_em-2-2-1_hlb-auto/hirsh_log_3_0V.out'
filename_mulliken_0 = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/cp2k-smeagol/sz/transmission/exp/capacitor/sergey-equal/auto/kpoints-100/V-0_kpoints_bulk-2-2-100_em-2-2-1_hlb-auto/mulliken_log_2_dft_wfn.out'
filename_mulliken_1 = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/cp2k-smeagol/sz/transmission/exp/capacitor/sergey-equal/auto/kpoints-100/V-0_kpoints_bulk-2-2-100_em-2-2-1_hlb-auto/mulliken_log_3_0V.out'
filename_mulliken_2 = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/cp2k-smeagol/sz/transmission/exp/capacitor/sergey-equal/auto/kpoints-100/V-1_kpoints_bulk-2-2-100_em-2-2-1_hlb-auto/mulliken_log_3_0V.out'

# folder_cp2k = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/cp2k-smeagol/sz/transmission/exp/capacitor/sergey-equal/auto/kpoints-100/V-0_kpoints_bulk-3-3-100_em-3-3-1_hlb-auto'
# filename_hirsh_0 = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/cp2k-smeagol/sz/transmission/exp/capacitor/sergey-equal/auto/kpoints-100/V-0_kpoints_bulk-3-3-100_em-3-3-1_hlb-auto/hirsh_log_2_dft_wfn.out'
# filename_hirsh_1 = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/cp2k-smeagol/sz/transmission/exp/capacitor/sergey-equal/auto/kpoints-100/V-0_kpoints_bulk-3-3-100_em-3-3-1_hlb-auto/hirsh_log_3_0V.out'
# filename_hirsh_2 = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/cp2k-smeagol/sz/transmission/exp/capacitor/sergey-equal/auto/kpoints-100/V-1_kpoints_bulk-3-3-100_em-3-3-1_hlb-auto/hirsh_log_3_0V.out'
# filename_mulliken_0 = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/cp2k-smeagol/sz/transmission/exp/capacitor/sergey-equal/auto/kpoints-100/V-0_kpoints_bulk-3-3-100_em-3-3-1_hlb-auto/mulliken_log_2_dft_wfn.out'
# filename_mulliken_1 = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/cp2k-smeagol/sz/transmission/exp/capacitor/sergey-equal/auto/kpoints-100/V-0_kpoints_bulk-3-3-100_em-3-3-1_hlb-auto/mulliken_log_3_0V.out'
# filename_mulliken_2 = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/cp2k-smeagol/sz/transmission/exp/capacitor/sergey-equal/auto/kpoints-100/V-1_kpoints_bulk-3-3-100_em-3-3-1_hlb-auto/mulliken_log_3_0V.out'

# folder_cp2k = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/cp2k-smeagol/sz/transmission/exp/capacitor/sergey-equal/auto/kpoints-100/V-0_kpoints_bulk-4-4-100_em-4-4-1_hlb-auto'
# filename_hirsh_0 = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/cp2k-smeagol/sz/transmission/exp/capacitor/sergey-equal/auto/kpoints-100/V-0_kpoints_bulk-4-4-100_em-4-4-1_hlb-auto/hirsh_log_2_dft_wfn.out'
# filename_hirsh_1 = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/cp2k-smeagol/sz/transmission/exp/capacitor/sergey-equal/auto/kpoints-100/V-0_kpoints_bulk-4-4-100_em-4-4-1_hlb-auto/hirsh_log_3_0V.out'
# filename_hirsh_2 = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/cp2k-smeagol/sz/transmission/exp/capacitor/sergey-equal/auto/kpoints-100/V-1_kpoints_bulk-4-4-100_em-4-4-1_hlb-auto/hirsh_log_3_0V.out'
# filename_mulliken_0 = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/cp2k-smeagol/sz/transmission/exp/capacitor/sergey-equal/auto/kpoints-100/V-0_kpoints_bulk-4-4-100_em-4-4-1_hlb-auto/mulliken_log_2_dft_wfn.out'
# filename_mulliken_1 = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/cp2k-smeagol/sz/transmission/exp/capacitor/sergey-equal/auto/kpoints-100/V-0_kpoints_bulk-4-4-100_em-4-4-1_hlb-auto/mulliken_log_3_0V.out'
# filename_mulliken_2 = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/cp2k-smeagol/sz/transmission/exp/capacitor/sergey-equal/auto/kpoints-100/V-1_kpoints_bulk-4-4-100_em-4-4-1_hlb-auto/mulliken_log_3_0V.out'

hirsh_0_Au, _, _, _, _ = read_hirsh(filename_hirsh_0)
hirsh_1_Au, _, _, _, _ = read_hirsh(filename_hirsh_1)
hirsh_2_Au, _, _, _, _ = read_hirsh(filename_hirsh_2)
mulliken_0_Au, _, _, _, _ = read_mulliken(filename_mulliken_0)
mulliken_1_Au, _, _, _, _ = read_mulliken(filename_mulliken_1)
mulliken_2_Au, _, _, _, _ = read_mulliken(filename_mulliken_2)

folder_siesta = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/siesta-smeagol/capacitor/bulk-2-2-100-em-2-2-1_hlb-15.2496_0-0_mulliken'
filename_mulliken_siesta_0 = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/siesta-smeagol/capacitor/bulk-2-2-100-em-2-2-1_hlb-15.2496_0-0_mulliken/mulliken_dft.out'
filename_mulliken_siesta_1 = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/siesta-smeagol/capacitor/bulk-2-2-100-em-2-2-1_hlb-15.2496_0-0_mulliken/mulliken_transport.out'
filename_mulliken_siesta_2 = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/siesta-smeagol/capacitor/V-1_bulk-2-2-100-em-2-2-1_hlb-15.2496_0-0_mulliken/mulliken_transport.out'
mulliken_siesta_0 = read_mulliken_siesta(filename_mulliken_siesta_0)
mulliken_siesta_1 = read_mulliken_siesta(filename_mulliken_siesta_1)
mulliken_siesta_2 = read_mulliken_siesta(filename_mulliken_siesta_2)
mulliken_siesta_0['Total'] = -1 * (mulliken_siesta_0['Total'] - 11)
mulliken_siesta_1['Total'] = -1 * (mulliken_siesta_1['Total'] - 11)
mulliken_siesta_2['Total'] = -1 * (mulliken_siesta_2['Total'] - 11)

markers = np.array([2.08400*0, 2.08400*1, 2.08400*2, 2.08400*3, 2.08400*4, 2.08400*5, 21.42200+2.08400*0,
                    21.42200+2.08400*1, 21.42200+2.08400*2, 21.42200+2.08400*3, 21.42200+2.08400*4, 21.42200+2.08400*5])
num_atoms = 96
atoms_per_layer = 8
num_layers = 96/atoms_per_layer
atoms_index = np.arange(start=1, stop=num_atoms+1)
layers = np.arange(start=1, stop=num_layers+1)
ylim = [0.085, -0.085]
ylim_siesta = [0.02, -0.02]

print(hirsh_0_Au['Charge'][::atoms_per_layer])
print(hirsh_1_Au['Charge'][::atoms_per_layer])
print(hirsh_2_Au['Charge'][::atoms_per_layer])

print(mulliken_0_Au['Charge'][::atoms_per_layer])
print(mulliken_1_Au['Charge'][::atoms_per_layer])
print(mulliken_2_Au['Charge'][::atoms_per_layer])

print(mulliken_siesta_0['Total'][::atoms_per_layer])
print(mulliken_siesta_1['Total'][::atoms_per_layer])
print(mulliken_siesta_2['Total'][::atoms_per_layer])

# Plot Hirsh
fig_plot_1, ax_plot_1 = plt.subplots()
ax_plot_1.plot(markers, hirsh_0_Au['Charge'][::atoms_per_layer], 'kx-', label='DFT')
ax_plot_1.plot(markers, hirsh_1_Au['Charge'][::atoms_per_layer], 'rx-', label='V=0')
ax_plot_1.plot(markers, hirsh_2_Au['Charge'][::atoms_per_layer], 'gx-', label='V=1')
ax_plot_1.legend(frameon=False)
ax_plot_1.set_xlabel('Z coordinate of layer / A')
ax_plot_1.set_ylabel('Average Hirshfeld charge')
fig_plot_1.tight_layout()
fig_plot_1.savefig('{}/hirsh.png'.format(folder_cp2k), dpi=param.save_dpi)

# Plot Mulliken CP2K
fig_plot_2, ax_plot_2 = plt.subplots()
ax_plot_2.plot(markers, mulliken_0_Au['Charge'][::atoms_per_layer], 'kx-', label='DFT')
ax_plot_2.plot(markers, mulliken_1_Au['Charge'][::atoms_per_layer], 'rx-', label='V=0')
ax_plot_2.plot(markers, mulliken_2_Au['Charge'][::atoms_per_layer], 'gx-', label='V=1')
ax_plot_2.set_ylim([ylim[0], ylim[1]])
ax_plot_2.legend(frameon=False)
ax_plot_2.set_xlabel('Z coordinate of layer / A')
ax_plot_2.set_ylabel('Average Mulliken charge')
fig_plot_2.tight_layout()
fig_plot_2.savefig('{}/mulliken.png'.format(folder_cp2k), dpi=param.save_dpi)

# Plot Mulliken SIESTA
fig_plot_siesta, ax_plot_siesta = plt.subplots()
ax_plot_siesta.plot(markers, mulliken_siesta_0['Total'][::atoms_per_layer], 'kx-', label='DFT')
ax_plot_siesta.plot(markers, mulliken_siesta_1['Total'][::atoms_per_layer], 'rx-', label='V=0')
ax_plot_siesta.plot(markers, mulliken_siesta_2['Total'][::atoms_per_layer], 'gx-', label='V=1')
ax_plot_siesta.set_ylim([ylim_siesta[0], ylim_siesta[1]])
ax_plot_siesta.legend(frameon=False)
ax_plot_siesta.set_xlabel('Z coordinate of layer / A')
ax_plot_siesta.set_ylabel('Average Mulliken charge')
fig_plot_siesta.tight_layout()
fig_plot_siesta.savefig('{}/mulliken_siesta.png'.format(folder_siesta), dpi=param.save_dpi)

# Plot Mulliken charge difference
fig_plot_3, ax_plot_3 = plt.subplots()
ax_plot_3.plot(markers, mulliken_1_Au['Charge'][::atoms_per_layer]-mulliken_2_Au['Charge'][::atoms_per_layer], 'kx-')
ax_plot_3.set_xlabel('Z coordinate of layer / A')
ax_plot_3.set_ylabel('Average Mulliken charge difference')
fig_plot_3.tight_layout()
fig_plot_3.savefig('{}/mulliken_difference.png'.format(folder_cp2k), dpi=param.save_dpi)

# Plot Hirsh charge difference
fig_plot_4, ax_plot_4 = plt.subplots()
ax_plot_4.plot(markers, hirsh_1_Au['Charge'][::atoms_per_layer]-hirsh_2_Au['Charge'][::atoms_per_layer], 'kx-')
ax_plot_4.set_xlabel('Z coordinate of layer / A')
ax_plot_4.set_ylabel('Average hirsh charge difference')
fig_plot_4.tight_layout()
fig_plot_4.savefig('{}/hirsh_difference.png'.format(folder_cp2k), dpi=param.save_dpi)

# Plot Mulliken charge SIESTA difference
fig_plotsiesta_diff, ax_plotsiesta_diff = plt.subplots()
ax_plotsiesta_diff.plot(markers, mulliken_siesta_1['Total'][::atoms_per_layer]-mulliken_siesta_2['Total'][::atoms_per_layer], 'kx-')
ax_plotsiesta_diff.set_xlabel('Z coordinate of layer / A')
ax_plotsiesta_diff.set_ylabel('Average Mulliken charge difference')
fig_plotsiesta_diff.tight_layout()
fig_plotsiesta_diff.savefig('{}/mulliken_difference_siesta.png'.format(folder_siesta), dpi=param.save_dpi)

# Plot Mulliken charge both difference
fig_plot_diff_all, ax_plot_diff_all = plt.subplots()
ax_plot_diff_all.plot(markers, mulliken_1_Au['Charge'][::atoms_per_layer]-mulliken_2_Au['Charge'][::atoms_per_layer], 'gx-', label='CP2K-SMEAGOL')
ax_plot_diff_all.plot(markers, mulliken_siesta_1['Total'][::atoms_per_layer]-mulliken_siesta_2['Total'][::atoms_per_layer], 'rx-', label='SIESTA-SMEAGOL')
ax_plot_diff_all.set_xlabel('Z coordinate of layer / A')
ax_plot_diff_all.set_ylabel('Average Mulliken charge difference')
ax_plot_diff_all.legend(frameon=False)
fig_plot_diff_all.tight_layout()
fig_plot_diff_all.savefig('{}/mulliken_difference_both.png'.format(folder_cp2k), dpi=param.save_dpi)
fig_plot_diff_all.savefig('{}/mulliken_difference_both.png'.format(folder_siesta), dpi=param.save_dpi)

if __name__ == "__main__":
    print('Finished.')
    plt.show()
