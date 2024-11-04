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
width = 0.6
x_lim = [-5, 4]  # x axis limit eV
use_ylim = True
y_lim = [-1, 300]  # y axis limit
folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-wire-ismael/jad/calculations/single-points/snapshots/restart-auto-long_au-frozen-dos/md_dft_long_LINEAR_P/step-2001/dft-gamma'
# folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-wire-ismael/jad/calculations/single-points/snapshots/restart-auto-long_au-frozen-dos/md_dft_long_LINEAR_P_dummy/step-2001/dft-gamma'
filename_pdos_Au_lead = '{}/2_dft_wfn-k1-1.pdos'.format(folder)
filename_pdos_Au_screen = '{}/2_dft_wfn-k2-1.pdos'.format(folder)
filename_pdos_Au_al = '{}/2_dft_wfn-k3-1.pdos'.format(folder)
filename_pdos_Au_bl = '{}/2_dft_wfn-k4-1.pdos'.format(folder)
filename_pdos_Au_cl = '{}/2_dft_wfn-k5-1.pdos'.format(folder)
filename_pdos_h = '{}/2_dft_wfn-k6-1.pdos'.format(folder)
filename_pdos_Au_wire = '{}/2_dft_wfn-k7-1.pdos'.format(folder)
filename_pdos_o = '{}/2_dft_wfn-k8-1.pdos'.format(folder)
filename_pdos_Au_br = '{}/2_dft_wfn-k9-1.pdos'.format(folder)
filename_pdos_Au_cr = '{}/2_dft_wfn-k10-1.pdos'.format(folder)

# Read Fermi energy from file
input_file = open(filename_pdos_Au_wire, 'r')
line_first = input_file.readline().strip().split()
fermi = float(line_first[15])

# Read energy and DOS from files
labels_au_dz_q11 = ['MO', 'Eigenvalue', 'Occupation', 's', 'py', 'pz', 'px', 'd-2', 'd-1', 'd0', 'd+1', 'd+2', 'f-3', 'f-2', 'f-1', 'f0', 'f+1', 'f+2', 'f+3']
labels_au_sz_q11 = ['MO', 'Eigenvalue', 'Occupation', 's', 'py', 'pz', 'px', 'd-2', 'd-1', 'd0', 'd+1', 'd+2']
labels_au_sz_q1 = ['MO', 'Eigenvalue', 'Occupation', 's']
labels_h = ['MO', 'Eigenvalue', 'Occupation', 's', 'py', 'pz', 'px']
labels_o = ['MO', 'Eigenvalue', 'Occupation', 's', 'py', 'pz', 'px', 'd-2', 'd-1', 'd0', 'd+1', 'd+2']

pdos_Au_lead = pd.read_csv(filename_pdos_Au_lead, names=labels_au_sz_q1, skiprows=[0, 1], delim_whitespace=True)
pdos_Au_screen = pd.read_csv(filename_pdos_Au_screen, names=labels_au_sz_q1, skiprows=[0, 1], delim_whitespace=True)
pdos_Au_al = pd.read_csv(filename_pdos_Au_al, names=labels_au_dz_q11, skiprows=[0, 1], delim_whitespace=True)
pdos_Au_bl = pd.read_csv(filename_pdos_Au_bl, names=labels_au_sz_q11, skiprows=[0, 1], delim_whitespace=True)
pdos_Au_cl = pd.read_csv(filename_pdos_Au_cl, names=labels_au_dz_q11, skiprows=[0, 1], delim_whitespace=True)
pdos_Au_br = pd.read_csv(filename_pdos_Au_br, names=labels_au_dz_q11, skiprows=[0, 1], delim_whitespace=True)
pdos_Au_cr = pd.read_csv(filename_pdos_Au_cr, names=labels_au_sz_q11, skiprows=[0, 1], delim_whitespace=True)
pdos_Au_wire = pd.read_csv(filename_pdos_Au_wire, names=labels_au_dz_q11, skiprows=[0, 1], delim_whitespace=True)

pdos_h = pd.read_csv(filename_pdos_h, names=labels_h, skiprows=[0, 1], delim_whitespace=True)
pdos_o = pd.read_csv(filename_pdos_o, names=labels_o, skiprows=[0, 1], delim_whitespace=True)

# Calculate pdos_Au_lead
pdos_Au_lead_s = sum_s(pdos_Au_lead, labels_au_sz_q1)

# Calculate pdos_Au_screen
pdos_Au_screen_s = sum_s(pdos_Au_screen, labels_au_sz_q1)

# Calculate pdos_Au_al
# There was a mistake in input file, species Au_al was mis-labelled as Au_a meaning it was DZVP not SZ
pdos_Au_al_s = sum_s(pdos_Au_al, labels_au_dz_q11)
pdos_Au_al_p = sum_p(pdos_Au_al, labels_au_dz_q11)
pdos_Au_al_d = sum_d(pdos_Au_al, labels_au_dz_q11)
pdos_Au_al_f = sum_f(pdos_Au_al, labels_au_dz_q11)
pdos_Au_al_all = sum_n(pdos_Au_al, labels_au_dz_q11, 3, len(labels_au_dz_q11))

# Calculate pdos_Au_bl
pdos_Au_bl_s = sum_s(pdos_Au_bl, labels_au_sz_q11)
pdos_Au_bl_p = sum_p(pdos_Au_bl, labels_au_sz_q11)
pdos_Au_bl_d = sum_d(pdos_Au_bl, labels_au_sz_q11)
pdos_Au_bl_all = sum_n(pdos_Au_bl, labels_au_sz_q11, 3, len(labels_au_sz_q11))

# Calculate pdos_Au_cl
pdos_Au_cl_s = sum_s(pdos_Au_cl, labels_au_dz_q11)
pdos_Au_cl_p = sum_p(pdos_Au_cl, labels_au_dz_q11)
pdos_Au_cl_d = sum_d(pdos_Au_cl, labels_au_dz_q11)
pdos_Au_cl_f = sum_f(pdos_Au_cl, labels_au_dz_q11)
pdos_Au_cl_all = sum_n(pdos_Au_cl, labels_au_dz_q11, 3, len(labels_au_dz_q11))

# Calculate pdos_Au_wire
pdos_Au_wire_s = sum_s(pdos_Au_wire, labels_au_dz_q11)
pdos_Au_wire_p = sum_p(pdos_Au_wire, labels_au_dz_q11)
pdos_Au_wire_d = sum_d(pdos_Au_wire, labels_au_dz_q11)
pdos_Au_wire_f = sum_f(pdos_Au_wire, labels_au_dz_q11)
pdos_Au_wire_all = sum_n(pdos_Au_wire, labels_au_dz_q11, 3, len(labels_au_dz_q11))

# Calculate pdos_Au_br
pdos_Au_br_s = sum_s(pdos_Au_br, labels_au_dz_q11)
pdos_Au_br_p = sum_p(pdos_Au_br, labels_au_dz_q11)
pdos_Au_br_d = sum_d(pdos_Au_br, labels_au_dz_q11)
pdos_Au_br_f = sum_f(pdos_Au_br, labels_au_dz_q11)
pdos_Au_br_all = sum_n(pdos_Au_br, labels_au_dz_q11, 3, len(labels_au_dz_q11))

# Calculate pdos_Au_cr
pdos_Au_cr_s = sum_s(pdos_Au_cr, labels_au_sz_q11)
pdos_Au_cr_p = sum_p(pdos_Au_cr, labels_au_sz_q11)
pdos_Au_cr_d = sum_d(pdos_Au_cr, labels_au_sz_q11)
pdos_Au_cr_all = sum_n(pdos_Au_cr, labels_au_sz_q11, 3, len(labels_au_sz_q11))

# Calculate pdos_Au_all
pdos_Au_all = pdos_Au_lead_s.values + pdos_Au_screen_s.values + pdos_Au_al_all.values + pdos_Au_bl_all.values \
              + pdos_Au_cl_all.values + pdos_Au_wire_all.values + pdos_Au_br_all.values + pdos_Au_cr_all.values

# Calculate pdos_o
pdos_o_s = sum_s(pdos_o, labels_o)
pdos_o_p = sum_p(pdos_o, labels_o)
pdos_o_d = sum_d(pdos_o, labels_o)
pdos_o_all = sum_n(pdos_o, labels_o, 3, len(labels_o))

# Calculate pdos_h
pdos_h_s = sum_s(pdos_h, labels_h)
pdos_h_all = sum_n(pdos_h, labels_h, 3, len(labels_h))

# Calculate convoluted PDOS
num_points = (pdos_Au_wire['s']).shape[0]
eigenvalues = (pdos_Au_wire['Eigenvalue'] - fermi) * param.hartree_to_ev
energy_grid = np.linspace(np.min(eigenvalues), np.max(eigenvalues), num=num_points)

# Plotting convoluted PDOS
fig_cpdos_1, ax_cpdos_1 = plt.subplots()
ax_cpdos_1.plot(energy_grid, smearing(eigenvalues, pdos_Au_lead_s.values, energy_grid, width), 'r', label='Au leads')
ax_cpdos_1.plot(energy_grid, smearing(eigenvalues, pdos_Au_all, energy_grid, width), 'g', label='Au')
ax_cpdos_1.plot(energy_grid, smearing(eigenvalues, pdos_o_all.values, energy_grid, width), 'b', label='O')
ax_cpdos_1.plot(energy_grid, smearing(eigenvalues, pdos_h_s.values, energy_grid, width), 'm', label='H')

# ax_cpdos_1.plot(energy_grid, smearing(eigenvalues, pdos_Au_lead_s.values, energy_grid, width), 'r', label='Au_lead SZ 6s')
# ax_cpdos_1.plot(energy_grid, smearing(eigenvalues, pdos_Au_screen_s.values, energy_grid, width), 'g', label='Au_screen SZ 6s')
# ax_cpdos_1.plot(energy_grid, smearing(eigenvalues, pdos_Au_al_all.values, energy_grid, width), 'b', label='Au_al DZVP')
# ax_cpdos_1.plot(energy_grid, smearing(eigenvalues, pdos_Au_bl_all.values, energy_grid, width), 'y', label='Au_bl SZ')
# ax_cpdos_1.plot(energy_grid, smearing(eigenvalues, pdos_Au_cl_all.values, energy_grid, width), 'm', label='Au_cl DZVP')
# ax_cpdos_1.plot(energy_grid, smearing(eigenvalues, pdos_Au_wire_all.values, energy_grid, width), 'k', label='Au_wire DZVP')
# ax_cpdos_1.plot(energy_grid, smearing(eigenvalues, pdos_Au_br_all.values, energy_grid, width), 'b', label='Au_br DZVP')
# ax_cpdos_1.plot(energy_grid, smearing(eigenvalues, pdos_Au_cr_all.values, energy_grid, width), 'y', label='Au_cr SZ')

ax_cpdos_1.set_xlim([x_lim[0], x_lim[1]])
if use_ylim: ax_cpdos_1.set_ylim([y_lim[0], y_lim[1]])
ax_cpdos_1.set_xlabel(r'E - E$_\mathrm{f}$ (eV)')
ax_cpdos_1.set_ylabel('DOS (arb units)')
ax_cpdos_1.legend(frameon=True)
fig_cpdos_1.tight_layout()
fig_cpdos_1.savefig('{}/cpdos_width-{}.png'.format(folder, width), dpi=param.save_dpi)

if __name__ == "__main__":
    print('Finished.')
    plt.show()
