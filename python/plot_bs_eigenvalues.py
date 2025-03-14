from __future__ import division, print_function
import pandas as pd
import numpy as np
import glob
from matplotlib import pyplot as plt
from general import parameters as param
from matplotlib.colors import LogNorm

"""
    Plot eigenvalues of effective Hamiltonian to identify bound states
"""

# folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/supercell-1-1-2-bulk-6-cu-1.86/junction/bias/bound-states/delete/kpoints-2-2-V-0-rs-1-bs-1-2-short2_scf-2-6hrs'
# folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/supercell-1-1-2-bulk-6-cu-1.86/junction/bias/bound-states/delete/kpoints-2-2-V-0-rs-0.1-bs-1-2-short2_scf-2-6hrs'
# folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/supercell-1-1-2-bulk-6-cu-1.86/junction/iv-curve/rs-dft/square-bs-eig-mpi-128/iv_curve/V_{}'.format(0.4)
folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/supercell-1-1-2-bulk-6-cu-1.86/junction/iv-curve/rs-dft-alpha-0.1/square-bs-eig-mpi-128/iv_curve/V_{}'.format(1.0)

# folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/supercell-1-1-2-bulk-6-cu-1.86/capacitor/bias/kpoints-2-2-V-1-bs-1'
# folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/supercell-1-1-2-bulk-6-cu-1.86/capacitor/bias/kpoints-2-2-V-1-bs-1-rs-1-2-parralel4'
# folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/supercell-1-1-2-bulk-6-cu-1.86/capacitor/bias/kpoints-2-2-V-1-bs-1-rs-1-2-parralel4-nobs-conv-dftrs-cores128'
# folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/interface/cu/tetragonal/transport/cu/supercell-1-1-2-bulk-6-cu-1.86/capacitor/iv-curve/sym/mpi-128/square-eig-delta-1e-8/iv_curve/V_{}'.format(1.0)

fermi = -0.17105102842723646

ylim = [-1, 1]
# 1 V
xlim = [-0.62, 0.62]
bias = np.array([0.5, -0.5])
# 0.5 V
# ylim = [-0.3, 0.3]
# xlim = [-0.42, 0.42]
# bias = np.array([0.25, -0.25])
# 0.4 V
# ylim = [-0.3, 0.3]
# xlim = [-0.32, 0.32]
# bias = np.array([0.2, -0.2])
# 0.1 V
# xlim = [-0.12, 0.12]
# bias = np.array([0.05, -0.05])

# folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/testing/capacitor/cp2k-smeagol/testing/kpoints-2-2-20_V-0.1_double-contour-square-mpi-128-eig'
# xlim = [-0.062, 0.062]
# ylim = [-0.2, 0.2]
# bias = np.array([0.05, -0.05])
# folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/testing/capacitor/cp2k-smeagol/testing/kpoints-2-2-20_V-1_double-contour-square-mpi-128-eig'
# xlim = [-0.62, 0.62]
# bias = np.array([0.5, -0.5])
# folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/testing/capacitor/cp2k-smeagol/testing/kpoints-2-2-20_V-3_double-contour-square-mpi-128-eig'
# xlim = [-2, 2.]
# ylim = [-1.8, 1.8]
# bias = np.array([1.5, -1.5])
# fermi = 0.33887644568033  # H
# fermi = 0.67775289136065853  # Ry

file = '{}/eigenvalues.out'.format(folder)
folder_save = folder

# cols = ['Energy', 'Real', 'Im', 'Spin', 'kvalue']
cols = ['Energy', 'Real', 'Im', 'Spin', 'kvalue', 'Eig_index', 'Energy_index']
skip_rows = 500000
skip_rows = int(479136/3)
skip_rows = 0
test = pd.read_csv(file, names=cols, delim_whitespace=True, on_bad_lines='skip', skiprows=skip_rows)
test = test.apply(pd.to_numeric, errors='coerce')
test = test.reset_index(drop=True)
print(test.head())
k_point = 1
save_figure = True
filtered_df = test.loc[test['kvalue'] == k_point]
# filtered_df = filtered_df[(filtered_df['Im']) < 0.0]
# filtered_df = filtered_df[-(filtered_df['Im']) < 1e-5]
print(filtered_df.head())
rydberg_to_ev = 0.5 * param.hartree_to_ev

# Plot real
# fig_plot_real, ax_plot_real = plt.subplots()
# ax_plot_real.plot((filtered_df['Energy']-fermi)*rydberg_to_ev, (filtered_df['Real']-fermi)*rydberg_to_ev, 'k.')
# ax_plot_real.set_ylim([ylim[0], ylim[1]])
# ax_plot_real.set_xlim([xlim[0], xlim[1]])
# ax_plot_real.plot([-10, 10], [-10, 10], 'r--', label='y = x')
# fig_plot_real.tight_layout()
# if save_figure: fig_plot_real.savefig('{}/eigenvalues_real_kpoint_{}.png'.format(folder_save, k_point), dpi=300)

# Plot imaginary
# fig_plot_im, ax_plot_im = plt.subplots()
# ax_plot_im.plot((filtered_df['Energy']-fermi)*rydberg_to_ev, filtered_df['Im'], 'k.')
# fig_plot_im.tight_layout()
# if save_figure: fig_plot_im.savefig('{}/eigenvalues_im.png_kpoint_{}.png'.format(folder_save, k_point), dpi=300)

# Plot real with imaginary colorbar
fig_plot_colorbar, ax_plot_colorbar, = plt.subplots()
vmin = 1e-16
# vmin = 1e-7
vmax = 1e-3
log_im = np.log10((filtered_df['Im']))
sc = ax_plot_colorbar.scatter(
    (filtered_df['Energy'] - fermi) * rydberg_to_ev,
    (filtered_df['Real'] - fermi) * rydberg_to_ev,
    c=np.abs(filtered_df['Im']),
    cmap='viridis',
    norm=LogNorm(vmin=vmin, vmax=vmax,),
    s=20,
)
ax_plot_colorbar.set_ylim([ylim[0], ylim[1]])
ax_plot_colorbar.set_xlim([xlim[0], xlim[1]])
ax_plot_colorbar.plot([bias[0], bias[0]], [-10, 10], 'r--', label='left v')
ax_plot_colorbar.plot([bias[1], bias[1]], [-10, 10], 'r--', label='right v')
ax_plot_colorbar.plot([-10, 10], [-10, 10], 'r--', label='y = x')
# ax_plot_colorbar.set_xlabel('Energy (eV)')
ax_plot_colorbar.set_xlabel(r'E-E$_{\mathrm{F}}$ (eV)')
ax_plot_colorbar.set_ylabel('Real Part (eV)')
cbar = fig_plot_colorbar.colorbar(sc, ax=ax_plot_colorbar)
# cbar.set_ticks([1e-14, 1e-13, 1e-12, 1e-11, 1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2])
# cbar.set_ticklabels(['1e-14', '1e-13', '1e-12', '1e-11', '1e-10', '1e-9', '1e-8', '1e-7', '1e-6', '1e-5', '1e-4', '1e-3', '1e-2', '1e-1'])
cbar.set_label('-Imaginary Part (Ry)')
fig_plot_colorbar.tight_layout()
if save_figure: fig_plot_colorbar.savefig('{}/eigenvalues_colorbar_kpoint_{}.png'.format(folder_save, k_point), dpi=300)

if __name__ == "__main__":
    print('Finished.')
    plt.show()