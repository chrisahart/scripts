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
        print('sum_d', i)
        total += x[label[i]]

    return total


def sum_d_eg(x, label):
    """Sum of d orbitals"""

    total = 0
    all = np.array([7, 8, 9, 10, 11], dtype=int) + 0
    t2g = np.array([0, 1, 3], dtype=int) + 7
    eg = np.array([2, 4], dtype=int) + 7
    for i in eg:
        print('sum_d eg', i)
        total += x[label[i]]

    return total


def sum_d_t2g(x, label):
    """Sum of d orbitals"""

    total = 0
    all = np.array([7, 8, 9, 10, 11], dtype=int) + 0
    t2g = np.array([0, 1, 3], dtype=int) + 7
    eg = np.array([2, 4], dtype=int) + 7
    for i in t2g:
        print('sum_d t2g', i)
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
width = 0.6
x_lim = [-3, 3]  # x axis limit eV
use_ylim = True
y_lim = [0, 0.3]  # y axis limit
folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/cp2k-smeagol/sz/pdos/dzvp/gamma'
filename_pdos_au = '{}/dft_wfn-k1-1.pdos'.format(folder)
filename_pdos_s = '{}/dft_wfn-k2-1.pdos'.format(folder)
filename_pdos_c = '{}/dft_wfn-k3-1.pdos'.format(folder)
filename_pdos_h = '{}/dft_wfn-k4-1.pdos'.format(folder)

# Read Fermi energy from file
input_file = open(filename_pdos_au, 'r')
line_first = input_file.readline().strip().split()
fermi = float(line_first[15])

# Read energy and DOS from files
labels_au = ['MO', 'Eigenvalue', 'Occupation', 's', 'py', 'pz', 'px', 'd-2', 'd-1', 'd0', 'd+1', 'd+2']
labels_s = ['MO', 'Eigenvalue', 'Occupation', 's', 'py', 'pz', 'px']
labels_c = ['MO', 'Eigenvalue', 'Occupation', 's', 'py', 'pz', 'px']
labels_h = ['MO', 'Eigenvalue', 'Occupation', 's']
pdos_au = pd.read_csv(filename_pdos_au, names=labels_au, skiprows=[0, 1], delim_whitespace=True)
pdos_s = pd.read_csv(filename_pdos_s, names=labels_s, skiprows=[0, 1], delim_whitespace=True)
pdos_c = pd.read_csv(filename_pdos_c, names=labels_c, skiprows=[0, 1], delim_whitespace=True)
pdos_h = pd.read_csv(filename_pdos_h, names=labels_h, skiprows=[0, 1], delim_whitespace=True)

# Calculate PDOS for Au
pdos_au_a_s = sum_s(pdos_au, labels_au)
pdos_au_a_p = sum_p(pdos_au, labels_au)
pdos_au_a_d = sum_d(pdos_au, labels_au)
pdos_au_a_all = sum_n(pdos_au, labels_au, 3, len(labels_au))

# Calculate PDOS for s
pdos_s_a_s = sum_s(pdos_s, labels_s)
pdos_s_a_p = sum_p(pdos_s, labels_s)
pdos_s_a_all = sum_n(pdos_s, labels_s, 3, len(labels_s))

# Calculate PDOS for s
pdos_c_a_s = sum_s(pdos_c, labels_c)
pdos_c_a_p = sum_p(pdos_c, labels_c)
pdos_c_a_all = sum_n(pdos_c, labels_c, 3, len(labels_c))

# Calculate PDOS for s
pdos_h_a_s = sum_s(pdos_h, labels_h)
pdos_h_a_all = sum_n(pdos_h, labels_h, 3, len(labels_h))

# Calculate PDOS for BDT
bdt = pdos_s_a_all.values + pdos_c_a_all.values + pdos_h_a_all.values

# Calculate convoluted PDOS
num_points = (pdos_au['s']).shape[0]
eigenvalues = (pdos_au['Eigenvalue'] - fermi) * param.hartree_to_ev
energy_grid = np.linspace(np.min(eigenvalues), np.max(eigenvalues), num=num_points)

# Plotting convoluted PDOS
fig_cpdos_1, ax_cpdos_1 = plt.subplots()
ax_cpdos_1.plot(energy_grid, smearing(eigenvalues, pdos_s_a_all.values, energy_grid, width), 'r', label='S')
ax_cpdos_1.plot(energy_grid, smearing(eigenvalues, pdos_c_a_all.values, energy_grid, width), 'g', label='C')
ax_cpdos_1.plot(energy_grid, smearing(eigenvalues, bdt, energy_grid, width), 'k', label='BDT DOS CP2K')
ax_cpdos_1.set_xlim([x_lim[0], x_lim[1]])
if use_ylim: ax_cpdos_1.set_ylim([y_lim[0], y_lim[1]])
ax_cpdos_1.set_xlabel(r'E - E$_\mathrm{f}$ (eV)')
ax_cpdos_1.set_ylabel('DOS (arb units)')
ax_cpdos_1.legend(frameon=True)
fig_cpdos_1.tight_layout()
fig_cpdos_1.savefig('{}/cpdos_width-{}.png'.format(folder, width), dpi=param.save_dpi, bbbox_inches='tight')

if __name__ == "__main__":
    print('Finished.')
    plt.show()
