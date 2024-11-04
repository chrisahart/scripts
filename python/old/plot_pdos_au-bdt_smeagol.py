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
y_lim = [0, 1]  # y axis limit
folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/siesta-smeagol/transmission/opt-cg/bulk-4-4-100-em-4-4-1_hlb-15.412-0-0_pdos'
filename_pdos_s = '{}/pdos_s.out'.format(folder)
filename_pdos_c = '{}/pdos_c.out'.format(folder)
labels = ['Eigenvalue', 'Occupation']
pdos_s = pd.read_csv(filename_pdos_s, names=labels, delim_whitespace=True)
pdos_c = pd.read_csv(filename_pdos_c, names=labels, delim_whitespace=True)

# Plotting convoluted PDOS including orbital labels
fig_cpdos_2, ax_cpdos_2 = plt.subplots()
ax_cpdos_2.plot(pdos_s['Eigenvalue'], pdos_s['Occupation'], 'y', label='S')
ax_cpdos_2.plot(pdos_c['Eigenvalue'], pdos_c['Occupation'],  color='brown', label='C')
ax_cpdos_2.set_xlim([x_lim[0], x_lim[1]])
if use_ylim: ax_cpdos_2.set_ylim([y_lim[0], y_lim[1]])
ax_cpdos_2.set_xlabel(r'E - E$_\mathrm{f}$ (eV)')
ax_cpdos_2.set_ylabel('DOS (arb units)')
ax_cpdos_2.legend(frameon=True)
fig_cpdos_2.tight_layout()
fig_cpdos_2.savefig('{}/pdos_smeagol.png'.format(folder, width), dpi=param.save_dpi, bbbox_inches='tight')

if __name__ == "__main__":
    print('Finished.')
    plt.show()
