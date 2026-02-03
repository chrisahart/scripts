import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from general import parameters as param

"""
    Plot .pdos file from CP2K
"""


def gaussian(energy_grid, eigenvalue, width):
    """Return Gaussian centred on eigenvalue over energy_grid with width"""

    x = -((energy_grid - eigenvalue) / width) ** 2

    return np.exp(x) / (np.sqrt(np.pi) * width)


def smearing(eigenvalues, pdos, energy_grid, width):
    """ Calculate convoluted PDOS by summing Gaussian distribution centred on each eigenvalue"""

    cpdos = np.zeros(eigenvalues.shape[0])
    for i in range(eigenvalues.shape[0]):
        cpdos += pdos[i] * gaussian(energy_grid, eigenvalues[i], width)

    return cpdos


def sum_s(x, label):
    """Sum of s orbitals"""

    total = 0
    total += x[label[3]]

    return total


def sum_p(x, label):
    """Sum of p orbitals"""

    total = 0
    for i in range(4, 7):
        total += x[label[i]]

    return total


def sum_d(x, label):
    """Sum of d orbitals"""

    total = 0
    all = np.array([7, 8, 9, 10, 11], dtype=int) + 0
    t2g = np.array([0, 1, 3], dtype=int) + 7
    eg = np.array([2, 4], dtype=int) + 7
    for i in all:
        # print('sum_d', i)
        total += x[label[i]]

    return total


def sum_d_eg(x, label):
    """Sum of d orbitals"""

    total = 0
    all = np.array([7, 8, 9, 10, 11], dtype=int) + 0
    t2g = np.array([0, 1, 3], dtype=int) + 7
    eg = np.array([2, 4], dtype=int) + 7
    for i in eg:
        # print('sum_d eg', i)
        total += x[label[i]]

    return total


def sum_d_t2g(x, label):
    """Sum of d orbitals"""

    total = 0
    all = np.array([7, 8, 9, 10, 11], dtype=int) + 0
    t2g = np.array([0, 1, 3], dtype=int) + 7
    eg = np.array([2, 4], dtype=int) + 7
    for i in t2g:
        # print('sum_d t2g', i)
        total += x[label[i]]

    return total


def sum_f(x, label):
    """Sum of d orbitals"""

    total = 0
    for i in range(12, 18):
        total += x[label[i]]

    return total


def sum_all(x, label):
    """Sum of all orbitals"""

    total = 0
    for i in range(3, len(label)):
        total += x[label[i]]

    return total


def sum_n(x, label, start, end):
    """Sum of all orbitals"""

    total = 0
    for i in range(start, end):
        total += x[label[i]]

    return total

# Projected Density of states
width = 0.15
x_lim = [-4, 6]  # x axis limit eV
use_ylim = True
y_lim = [-1, 290]  # y axis limit
# folder = '/Volumes/Elements/Data/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-336/geo-opt-cell-opt/neutral/other/hse-22-amd-n-1-pdos-polaron'
# folder = '/Volumes/Elements/Data/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-336/geo-opt-cell-opt/electron/other/hse-22-amd-n-1-pdos-polaron'
folder = '/Volumes/Elements/Data/Postdoc2/Data/Work/calculations/tio2/anatase/archer/anatase/cell-442/geo-opt-cell-opt-hse-20/hole/other/hse-19'
# folder = '/Volumes/Elements/Data/Postdoc2/Data/Work/calculations/tio2/anatase/archer/anatase/cell-442/geo-opt-cell-opt-hse-20/neutral/neutral-opt/other/hse-19'

# filename_pdos_ti_1_polaron = '{}/tio2-ALPHA_k2-1.pdos'.format(folder)
# filename_pdos_ti_2_polaron = '{}/tio2-BETA_k-1.pdos'.format(folder)
# filename_pdos_ti_1 = '{}/tio2-ALPHA_k1-1.pdos'.format(folder)
# filename_pdos_ti_2 = '{}/tio2-BETA_k1-1.pdos'.format(folder)
# filename_pdos_o_1 = '{}/tio2-ALPHA_k3-1.pdos'.format(folder)
# filename_pdos_o_2 = '{}/tio2-BETA_k3-1.pdos'.format(folder)

filename_pdos_o_1_polaron = '{}/tio2-BETA_k3-1.pdos'.format(folder)
filename_pdos_o_2_polaron = '{}/tio2-BETA_k3-1.pdos'.format(folder)
filename_pdos_ti_1 = '{}/tio2-BETA_k1-1.pdos'.format(folder)
filename_pdos_ti_2 = '{}/tio2-BETA_k1-1.pdos'.format(folder)
filename_pdos_o_1 = '{}/tio2-BETA_k2-1.pdos'.format(folder)
filename_pdos_o_2 = '{}/tio2-BETA_k2-1.pdos'.format(folder)

# Read Fermi energy from file
input_file = open(filename_pdos_ti_1, 'r')
line_first = input_file.readline().strip().split()
fermi = float(line_first[15])
# fermi = 0.235319
fermi = 0.169388
print(fermi)

# Read energy and DOS from files
labels_ti = ['MO', 'Eigenvalue', 'Occupation', 's', 'py', 'pz', 'px', 'd-2', 'd-1', 'd0', 'd+1', 'd+2', 'f-3', 'f-2', 'f-1', 'f0', 'f+1', 'f+2', 'f+3']
labels_o = ['MO', 'Eigenvalue', 'Occupation', 's', 'py', 'pz', 'px', 'd-2', 'd-1', 'd0', 'd+1', 'd+2']

pdos_o_1_polaron = pd.read_csv(filename_pdos_o_1_polaron, names=labels_o, skiprows=[0, 1], delim_whitespace=True)
pdos_o_2_polaron = pd.read_csv(filename_pdos_o_2_polaron, names=labels_o, skiprows=[0, 1], delim_whitespace=True)
# pdos_ti_1_polaron = pd.read_csv(filename_pdos_ti_1_polaron, names=labels_ti, skiprows=[0, 1], delim_whitespace=True)
# pdos_ti_2_polaron = pd.read_csv(filename_pdos_ti_2_polaron, names=labels_ti, skiprows=[0, 1], delim_whitespace=True)
pdos_ti_1 = pd.read_csv(filename_pdos_ti_1, names=labels_ti, skiprows=[0, 1], delim_whitespace=True)
pdos_ti_2 = pd.read_csv(filename_pdos_ti_2, names=labels_ti, skiprows=[0, 1], delim_whitespace=True)
pdos_o_1 = pd.read_csv(filename_pdos_o_1, names=labels_o, skiprows=[0, 1], delim_whitespace=True)
pdos_o_2 = pd.read_csv(filename_pdos_o_2, names=labels_o, skiprows=[0, 1], delim_whitespace=True)

# pdos_ti_1_polaron_d = sum_d(pdos_ti_1_polaron, labels_ti)
pdos_o_1_polaron_p = sum_p(pdos_o_1_polaron, labels_o)
pdos_ti_1_d = sum_d(pdos_ti_1, labels_ti)
pdos_o_1_p = sum_p(pdos_o_1, labels_o)

# # Calculate convoluted PDOS
num_points = (pdos_ti_1['s']).shape[0]
eigenvalues = (pdos_ti_1['Eigenvalue'] - fermi) * param.hartree_to_ev
energy_grid = np.linspace(np.min(eigenvalues), np.max(eigenvalues), num=num_points)

# fermi_index = np.argmin(np.abs(energy_grid))
# pdos_ti_1_polaron = pdos_ti_1.copy()
# pdos_ti_1_d_polaron2 = sum_d(pdos_ti_1_polaron, labels_ti)
# pdos_ti_1_d_polaron2[fermi_index] = pdos_ti_1_d_polaron2[fermi_index] * 1e3

# # Plotting convoluted PDOS
fig_cpdos_1, ax_cpdos_1 = plt.subplots()
# ax_cpdos_1.plot(energy_grid, smearing(eigenvalues, pdos_o_1_polaron_p.values*10, energy_grid, width), 'k', label='Ti 3d polaron * 10')
ax_cpdos_1.plot(energy_grid, smearing(eigenvalues, pdos_o_1_polaron_p.values*10, energy_grid, width), 'k', label='O 2p polaron * 10')
ax_cpdos_1.plot(energy_grid, smearing(eigenvalues, pdos_ti_1_d.values, energy_grid, width), 'b', label='Ti 3d')
ax_cpdos_1.plot(energy_grid, smearing(eigenvalues, pdos_o_1_p.values, energy_grid, width), 'r', label='O 2p')
# ax_cpdos_1.plot(energy_grid, smearing(eigenvalues, pdos_ti_1_d_polaron2.values, energy_grid, width), 'k', label='Ti 3d polaron * 10')
ax_cpdos_1.set_xlim([x_lim[0], x_lim[1]])
ax_cpdos_1.set_ylim([y_lim[0], y_lim[1]])
ax_cpdos_1.set_xlabel(r'E - E$_\mathrm{f}$ (eV)')
ax_cpdos_1.set_ylabel('DOS (arb units)')
ax_cpdos_1.legend(frameon=True)
fig_cpdos_1.tight_layout()
fig_cpdos_1.savefig('{}/cpdos_width-{}.png'.format(folder, width), dpi=param.save_dpi)

if __name__ == "__main__":
    print('Finished.')
    plt.show()
