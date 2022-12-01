import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scripts.general import parameters as param

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


# Projected Density of states
width = 0.2  # Gaussian width, 0.2 works for m400, 0.1 for m800, 0.05 for m1200
x_lim = [-5, 10]  # x axis limit eV
folder = '/Volumes/Storage/Data/Work/Postdoc/Work/testing/cp2k-smeagol/examples/SMEAGOL_TEST/Au_chain/dft-pdos-50atoms-q11'
filename_pdos_au = '{}/dft_wfn-k1-1.pdos'.format(folder)
labels_au = ['MO', 'Eigenvalue', 'Occupation', 's', 'py', 'pz', 'px', 'd-2', 'd-1', 'd0', 'd+1', 'd+2']

# width = 0.3  # Gaussian width, 0.2 works for m400, 0.1 for m800, 0.05 for m1200
# x_lim = [-8, 15]  # x axis limit eV
# folder = '/Volumes/Storage/Data/Work/Postdoc/Work/testing/cp2k-smeagol/bdt/cp2k/bulk/pdos'
# filename_pdos_au = '{}/dft_wfn-k1-1.pdos'.format(folder)

# Read Fermi energy from file
input_file = open(filename_pdos_au, 'r')
line_first = input_file.readline().strip().split()
fermi = float(line_first[15])

# Read energy and DOS from files
pdos_au = pd.read_csv(filename_pdos_au, names=labels_au, skiprows=[0, 1], delim_whitespace=True)

# Calculate PDOS for O
pdos_au_a_s = sum_s(pdos_au, labels_au)
pdos_au_a_p = sum_p(pdos_au, labels_au)
pdos_au_a_d = sum_d(pdos_au, labels_au)
pdos_au_a_d_eg = sum_d_eg(pdos_au, labels_au)
pdos_au_a_d_t2g = sum_d_t2g(pdos_au, labels_au)
pdos_au_a_all = sum_all(pdos_au, labels_au)

# Calculate convoluted PDOS
num_points = (pdos_au['s']).shape[0]
eigenvalues = (pdos_au['Eigenvalue'] - fermi) * param.hartree_to_ev
energy_grid = np.linspace(np.min(eigenvalues), np.max(eigenvalues), num=num_points)

# Plotting convoluted PDOS
fig_cpdos_1, ax_cpdos_1 = plt.subplots()
ax_cpdos_1.plot(energy_grid, smearing(eigenvalues, pdos_au_a_s.values, energy_grid, width), 'b', label='Au(s)')
ax_cpdos_1.plot(energy_grid, smearing(eigenvalues, pdos_au_a_p.values, energy_grid, width), 'r', label='Au(p)')
ax_cpdos_1.plot(energy_grid, smearing(eigenvalues, pdos_au_a_d.values, energy_grid, width), 'g', label='Au(d)')
ax_cpdos_1.set_xlim([x_lim[0], x_lim[1]])
ax_cpdos_1.set_xlabel(r'E - E$_\mathrm{f}$ (eV)')
ax_cpdos_1.set_ylabel('DOS (arb units)')
ax_cpdos_1.legend(frameon=True)
fig_cpdos_1.tight_layout()
fig_cpdos_1.savefig('{}/cpdos_width-{}.png'.format(folder, width), dpi=param.save_dpi, bbbox_inches='tight')

# Plotting convoluted PDOS including orbital labels
# t2g = [0, 1, 3]
# eg = [2, 4]
# fig_cpdos_2, ax_cpdos_2 = plt.subplots()
# ax_cpdos_2.plot(energy_grid, smearing(eigenvalues, pdos_au_a_s.values, energy_grid, width), 'b', label='Au(s)')
# ax_cpdos_2.plot(energy_grid, smearing(eigenvalues, pdos_au_a_p.values, energy_grid, width), 'r', label='Au(p)')
# ax_cpdos_2.plot(energy_grid, smearing(eigenvalues, pdos_au[labels_au[7]], energy_grid, width), 'r', label='d t2g')  # -2 does not seem to be dx2-y2
# ax_cpdos_2.plot(energy_grid, smearing(eigenvalues, pdos_au[labels_au[8]], energy_grid, width), 'g', label='d eg')  # -1 does not seem to be dyz?
# ax_cpdos_2.plot(energy_grid, smearing(eigenvalues, pdos_au[labels_au[9]], energy_grid, width), 'b', label=r'd$z^2$')  # 0 does not seem to be dz^2?
# ax_cpdos_2.plot(energy_grid, smearing(eigenvalues, pdos_au[labels_au[10]], energy_grid, width), 'y', label=r'd$xz$')  # 1 does not seem to be dxz?
# ax_cpdos_2.plot(energy_grid, smearing(eigenvalues, pdos_au[labels_au[11]], energy_grid, width), 'm', label=r'd$xy$')  # 1 does not seem to be dxy?
# ax_cpdos_2.plot(energy_grid, smearing(eigenvalues, pdos_au_a_d_eg.values, energy_grid, width), 'y--', label=r'Au(d) e$_{\mathrm{g}}$')
# ax_cpdos_2.plot(energy_grid, smearing(eigenvalues, pdos_au_a_d_t2g.values, energy_grid, width), 'c--', label=r'Au(d) t$_{\mathrm{2g}}$')
# ax_cpdos_2.plot(energy_grid, smearing(eigenvalues, pdos_au_a_d.values, energy_grid, width), 'g', label='Au(d)')
# ax_cpdos_2.set_xlim([-5, 10])
# ax_cpdos_2.set_xlabel(r'E - E$_\mathrm{f}$ (eV)')
# ax_cpdos_2.set_ylabel('DOS (arb units)')
# ax_cpdos_2.legend(frameon=True)
# fig_cpdos_2.tight_layout()
# fig_cpdos_2.savefig('{}/cpdos_width-{}.png'.format(folder, width), dpi=param.save_dpi, bbbox_inches='tight')

# Plot from smeared.dat https://wiki.wpi.edu/deskinsgroup/Density_of_States_-_CP2K
# file_1_1 = np.genfromtxt('{}/smeared.dat'.format(folder))
# fig_pdos_dat, ax_pdos_dat = plt.subplots()
# ax_pdos_dat.plot(file_1_1[:, 0], file_1_1[:, 1], 'b.-', label='Au(s)')
# ax_pdos_dat.set_xlim([-5, 10])

if __name__ == "__main__":
    print('Finished.')
    plt.show()
