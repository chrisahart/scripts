from __future__ import division, print_function
import pandas as pd
import numpy as np
import glob
from general import load_coordinates
from general import parameters as param
from general import print_xyz
import matplotlib.pyplot as plt
import math
import ase
from ase.build import surface
from ase.visualize import view
from ase import Atoms
from ase.build import bulk
from ase.io import read, write
from ase.build import make_supercell
from ase.build import fcc111
from pathlib import Path

"""
    Analysis for melamine
"""

params = {'axes.formatter.limits': [-4, 4],
          'axes.labelsize': 'x-large',
          'axes.titlesize': 'x-large',
          'legend.fontsize': 'large',
          'lines.markersize': '8',
          }
plt.rcParams.update(params)

# CP2K /Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/melamine/hx1/au_from-h2
v_bias_cp2k = np.array([-2, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2])
save_folder_cp2k = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/melamine/hx1/au_from-h2/iv-curve/single-point'

# without tip
# v_gs_without_tip_cp2k = np.loadtxt('/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/melamine/hx1/au_from-h2/iv-curve/single-point/gs-notip/energy.out')
v_c1_without_tip_cp2k = np.loadtxt('/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/melamine/hx1/au_from-h2/iv-curve/single-point/c1-notip/energy.out')
v_c2_without_tip_cp2k = np.loadtxt('/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/melamine/hx1/au_from-h2/iv-curve/single-point/c2-notip/energy.out')
v_ts2_without_tip_cp2k = np.loadtxt('/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/melamine/hx1/au_from-h2/iv-curve/single-point/ts2-guess-notip/energy.out')

# with tip 5A
v_gs_5A_cp2k = np.loadtxt('/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/melamine/hx1/au_from-h2/iv-curve/single-point/gs/energy.out')
v_c1_5A_cp2k = np.loadtxt('/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/melamine/hx1/au_from-h2/iv-curve/single-point/c1/energy.out')
v_c2_5A_cp2k = np.loadtxt('/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/melamine/hx1/au_from-h2/iv-curve/single-point/c2/energy.out')
v_ts2_5A_cp2k = np.loadtxt('/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/melamine/hx1/au_from-h2/iv-curve/single-point/ts2-guess/energy.out')

# with tip 3A
# v_gs_3A_cp2k = np.loadtxt('/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/melamine/hx1/au_from-h2/iv-curve/single-point/gs-3A/energy.out')
v_c1_3A_cp2k = np.loadtxt('/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/melamine/hx1/au_from-h2/iv-curve/single-point/c1-3A/energy.out')
v_c2_3A_cp2k = np.loadtxt('/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/melamine/hx1/au_from-h2/iv-curve/single-point/c2-3A/energy.out')
v_ts2_3A_cp2k = np.loadtxt('/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/melamine/hx1/au_from-h2/iv-curve/single-point/ts2-3A/energy.out')

v_gs_cp2k = v_gs_5A_cp2k
v_c1_cp2k = v_c1_5A_cp2k
v_c2_cp2k = v_c2_5A_cp2k
v_ts2_cp2k = v_ts2_5A_cp2k
print('cp2k C1', (v_gs_cp2k[2]-v_c1_cp2k[2])*param.hartree_to_ev)
print('cp2k C2', (v_gs_cp2k[2]-v_c2_cp2k[2])*param.hartree_to_ev)
print('cp2k TS2', (v_gs_cp2k[2]-v_ts2_cp2k[2])*param.hartree_to_ev)

# SIESTA /Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/melamine/hx1/au_from-h2/iv-curve/single-point
# v_bias_siesta = np.array([-2, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2])
v_bias_siesta = np.array([-1.0, -0.5, 0.0, 0.5, 1.0])
save_folder_siesta = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/melamine/young/siesta1-smeagol/au_from-h2/iv-curve'

v_gs_5A_siesta = np.loadtxt('/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/melamine/young/siesta1-smeagol/au_from-h2/iv-curve/gs_leads-sz-custom_surf-tip-dzp/energy.out')
v_c1_5A_siesta = np.loadtxt('/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/melamine/young/siesta1-smeagol/au_from-h2/iv-curve/c1_leads-sz-custom_surf-tip-dzp/energy.out')
v_c2_5A_siesta = np.loadtxt('/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/melamine/young/siesta1-smeagol/au_from-h2/iv-curve/c2_leads-sz-custom_surf-tip-dzp/energy.out')
v_ts1_5A_siesta = np.loadtxt('/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/melamine/young/siesta1-smeagol/au_from-h2/iv-curve/ts1_leads-sz-custom_surf-tip-dzp_from-cu/energy.out')
v_ts2_5A_siesta = np.loadtxt('/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/melamine/young/siesta1-smeagol/au_from-h2/iv-curve/ts2_leads-sz-custom_surf-tip-dzp_from-cu2/energy.out')

v_gs_3A_siesta = np.loadtxt('/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/melamine/young/siesta1-smeagol/au_from-h2/iv-curve-3A/gs_leads-sz-custom_surf-tip-dzp_from-cu/energy.out')
v_c1_3A_siesta = np.loadtxt('/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/melamine/young/siesta1-smeagol/au_from-h2/iv-curve-3A/c1_leads-sz-custom_surf-tip-dzp_from-cu/energy.out')
v_c2_3A_siesta = np.loadtxt('/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/melamine/young/siesta1-smeagol/au_from-h2/iv-curve-3A/c2_leads-sz-custom_surf-tip-dzp_from-cu/energy.out')
v_ts2_3A_siesta = np.loadtxt('/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/melamine/young/siesta1-smeagol/au_from-h2/iv-curve-3A/ts2_leads-sz-custom_surf-tip-dzp_from-cu/energy.out')

# v_gs_without_tip_siesta = np.loadtxt('/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/melamine/young/siesta1-smeagol/au_from-h2/iv-curve-no-tip/gs_leads-sz-custom_surf-tip-dzp_from-cu/energy.out')
# v_c1_without_tip_siesta = np.loadtxt('/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/melamine/young/siesta1-smeagol/au_from-h2/iv-curve-no-tip/c1_leads-sz-custom_surf-tip-dzp_from-cu/energy.out')
# v_c2_without_tip_siesta = np.loadtxt('/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/melamine/young/siesta1-smeagol/au_from-h2/iv-curve-no-tip/c2_leads-sz-custom_surf-tip-dzp_from-cu/energy.out')
# v_ts2_without_tip_siesta = np.loadtxt('/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/melamine/young/siesta1-smeagol/au_from-h2/iv-curve-no-tip/ts2_leads-sz-custom_surf-tip-dzp_from-cu/energy.out')

v_gs_siesta = v_gs_5A_siesta
v_c1_siesta = v_c1_5A_siesta
v_c2_siesta = v_c2_5A_siesta
v_ts2_siesta = v_ts2_5A_siesta
print('siesta C1', (v_gs_siesta[2]-v_c1_siesta[2])*1)
print('siesta C2', (v_gs_siesta[2]-v_c2_siesta[2])*1)
print('siesta TS2', (v_gs_siesta[2]-v_ts2_siesta[2])*1)

c1_ts2_without_tip_cp2k = []
c1_ts2_3A_cp2k = []
c1_ts2_5A_cp2k = []
c2_ts2_without_tip_cp2k = []
c2_ts2_3A_cp2k = []
c2_ts2_5A_cp2k = []
for i in range(v_bias_cp2k.shape[0] - 1):
    c1_ts2_without_tip_cp2k.append(v_c1_without_tip_cp2k[i] - v_ts2_without_tip_cp2k[i])
    c1_ts2_3A_cp2k.append(v_c1_3A_cp2k[i] - v_ts2_3A_cp2k[i])
    c1_ts2_5A_cp2k.append(v_c1_5A_cp2k[i] - v_ts2_5A_cp2k[i])
    c2_ts2_without_tip_cp2k.append(v_c2_without_tip_cp2k[i] - v_ts2_without_tip_cp2k[i])
    c2_ts2_3A_cp2k.append(v_c2_3A_cp2k[i] - v_ts2_3A_cp2k[i])
    c2_ts2_5A_cp2k.append(v_c2_5A_cp2k[i] - v_ts2_5A_cp2k[i])

c1_ts2_without_tip_siesta = []
c1_ts2_3A_siesta = []
c1_ts2_5A_siesta = []
c2_ts2_without_tip_siesta = []
c2_ts2_3A_siesta = []
c2_ts2_5A_siesta = []

for i in range(v_bias_siesta.shape[0]):
    # c1_ts2_without_tip_siesta.append(v_c1_without_tip_siesta[i] - v_ts2_without_tip_siesta[i])
    c1_ts2_3A_siesta.append(v_c1_3A_siesta[i] - v_ts2_3A_siesta[i])
    c1_ts2_5A_siesta.append(v_c1_5A_siesta[i] - v_ts2_5A_siesta[i])
    # c2_ts2_without_tip_siesta.append(v_c2_without_tip_siesta[i] - v_ts2_without_tip_siesta[i])
    c2_ts2_3A_siesta.append(v_c2_3A_siesta[i] - v_ts2_3A_siesta[i])
    c2_ts2_5A_siesta.append(v_c2_5A_siesta[i] - v_ts2_5A_siesta[i])
    
# fig_cp2k_c1, ax_cp2k_c1 = plt.subplots()
# ax_cp2k_c1.plot(v_bias[2:-2], (c1_ts2_5A_cp2k[2:-1] - c1_ts2_5A_cp2k[4]) * param.hartree_to_ev, 'ms-', label='5 A')
# ax_cp2k_c1.plot(v_bias[2:-2], (c1_ts2_3A_cp2k[2:-1] - c1_ts2_3A_cp2k[4]) * param.hartree_to_ev, 'r^-', label='3 A')
# ax_cp2k_c1.plot(v_bias[2:-2], (c1_ts2_without_tip_cp2k[2:-1] - c1_ts2_without_tip_cp2k[4]) * param.hartree_to_ev, 'bx-', label='Without tip')
# ax_cp2k_c1.legend(frameon=False)
# ax_cp2k_c1.set_xlabel('Bias / V')
# ax_cp2k_c1.set_ylabel('Barrier Change / eV')
# fig_cp2k_c1.tight_layout()

# fig_cp2k_c2, ax_cp2k_c2 = plt.subplots()
# ax_cp2k_c2.plot(v_bias[2:-2], (c2_ts2_5A_cp2k[2:-1] - c2_ts2_5A_cp2k[4]) * param.hartree_to_ev, 'ms-', label='5 A')
# ax_cp2k_c2.plot(v_bias[2:-2], (c2_ts2_3A_cp2k[2:-1] - c2_ts2_3A_cp2k[4]) * param.hartree_to_ev, 'r^-', label='3 A')
# ax_cp2k_c2.plot(v_bias[2:-2], (c2_ts2_without_tip_cp2k[2:-1] - c2_ts2_without_tip_cp2k[4]) * param.hartree_to_ev, 'bx-', label='Without tip')
# ax_cp2k_c2.legend(frameon=False)
# ax_cp2k_c2.set_xlabel('Bias / V')
# ax_cp2k_c2.set_ylabel('Barrier Change / eV')
# fig_cp2k_c2.tight_layout()

# rows, cols = 2, 1
# fig_cp2k_c1_c2, ax_cp2k_c1_c2 = plt.subplots(rows, cols, sharex='col', sharey='row', figsize=(6, 8))
# ax_cp2k_c1_c2[0].plot(v_bias_cp2k[2:-2], (c1_ts2_5A_cp2k[2:-1] - c1_ts2_5A_cp2k[4]) * param.hartree_to_ev, 'ms-', label='5 Å')
# ax_cp2k_c1_c2[0].plot(v_bias_cp2k[2:-2], (c1_ts2_3A_cp2k[2:-1] - c1_ts2_3A_cp2k[4]) * param.hartree_to_ev, 'r^-', label='3 Å')
# ax_cp2k_c1_c2[1].plot(v_bias_cp2k[2:-2], (c2_ts2_5A_cp2k[2:-1] - c2_ts2_5A_cp2k[4]) * param.hartree_to_ev, 'ms-', label='5 Å')
# ax_cp2k_c1_c2[1].plot(v_bias_cp2k[2:-2], (c2_ts2_3A_cp2k[2:-1] - c2_ts2_3A_cp2k[4]) * param.hartree_to_ev, 'r^-', label='3 Å')
# ax_cp2k_c1_c2[1].set_xlabel('Bias / V')
# ax_cp2k_c1_c2[0].set_ylabel('Barrier Change / eV')
# ax_cp2k_c1_c2[1].set_ylabel('Barrier Change / eV')
# ax_cp2k_c1_c2[0].set_ylim([-0.064, 0.007])
# ax_cp2k_c1_c2[1].set_ylim([-0.064, 0.007])
# ax_cp2k_c1_c2[0].legend(frameon=False)
# fig_cp2k_c1_c2.tight_layout()
# fig_cp2k_c1_c2.savefig('{}/cp2k_c1_c2.png'.format(save_folder_siesta), dpi=300)
# fig_cp2k_c1_c2.savefig('{}/cp2k_c1_c2.png'.format(save_folder_cp2k), dpi=300)

# rows, cols = 2, 1
# fig_siesta_c1_c2, ax_siesta_c1_c2 = plt.subplots(rows, cols, sharex='col', sharey='row', figsize=(6, 8))
# ax_siesta_c1_c2[0].plot(v_bias_siesta, (c1_ts2_5A_siesta - c1_ts2_5A_siesta[2]) * 1, 'ms-', label='5 Å')
# ax_siesta_c1_c2[0].plot(v_bias_siesta, (c1_ts2_3A_siesta - c1_ts2_3A_siesta[2]) * 1, 'r^-', label='3 Å')
# ax_siesta_c1_c2[1].plot(v_bias_siesta, (c2_ts2_5A_siesta - c2_ts2_5A_siesta[2]) * 1, 'ms-', label='5 Å')
# ax_siesta_c1_c2[1].plot(v_bias_siesta, (c2_ts2_3A_siesta - c2_ts2_3A_siesta[2]) * 1, 'r^-', label='3 Å')
# ax_siesta_c1_c2[1].set_xlabel('Bias / V')
# ax_siesta_c1_c2[0].set_ylabel('Barrier Change / eV')
# ax_siesta_c1_c2[1].set_ylabel('Barrier Change / eV')
# ax_siesta_c1_c2[0].set_ylim([-0.064, 0.007])
# ax_siesta_c1_c2[1].set_ylim([-0.064, 0.007])
# ax_siesta_c1_c2[0].legend(frameon=False)
# fig_siesta_c1_c2.tight_layout()
# fig_siesta_c1_c2.savefig('{}/siesta_c1_c2.png'.format(save_folder_siesta), dpi=300)
# fig_siesta_c1_c2.savefig('{}/siesta_c1_c2.png'.format(save_folder_cp2k), dpi=300)

rows, cols = 2, 2
fig_siesta_cp2k_c1_c2, ax_siesta_cp2k_c1_c2 = plt.subplots(rows, cols, sharex='col', sharey='row', figsize=(10, 8))
ax_siesta_cp2k_c1_c2[0, 0].plot(v_bias_cp2k[2:-2], (c1_ts2_5A_cp2k[2:-1] - c1_ts2_5A_cp2k[4]) * param.hartree_to_ev, 'ms-', label='CP2K+SMEAGOL 5 Å')
ax_siesta_cp2k_c1_c2[0, 0].plot(v_bias_cp2k[2:-2], (c1_ts2_3A_cp2k[2:-1] - c1_ts2_3A_cp2k[4]) * param.hartree_to_ev, 'r^-', label='CP2K+SMEAGOL 3 Å')
ax_siesta_cp2k_c1_c2[1, 0].plot(v_bias_cp2k[2:-2], (c2_ts2_5A_cp2k[2:-1] - c2_ts2_5A_cp2k[4]) * param.hartree_to_ev, 'ms-', label='CP2K+SMEAGOL 5 Å')
ax_siesta_cp2k_c1_c2[1, 0].plot(v_bias_cp2k[2:-2], (c2_ts2_3A_cp2k[2:-1] - c2_ts2_3A_cp2k[4]) * param.hartree_to_ev, 'r^-', label='CP2K+SMEAGOL 3 Å')
ax_siesta_cp2k_c1_c2[0, 1].plot(v_bias_siesta, (c1_ts2_5A_siesta - c1_ts2_5A_siesta[2]) * 1, 'ms-', label='SIESTA+SMEAGOL 5 Å')
ax_siesta_cp2k_c1_c2[0, 1].plot(v_bias_siesta, (c1_ts2_3A_siesta - c1_ts2_3A_siesta[2]) * 1, 'r^-', label='SIESTA+SMEAGOL 3 Å')
ax_siesta_cp2k_c1_c2[1, 1].plot(v_bias_siesta, (c2_ts2_5A_siesta - c2_ts2_5A_siesta[2]) * 1, 'ms-', label='SIESTA+SMEAGOL 5 Å')
ax_siesta_cp2k_c1_c2[1, 1].plot(v_bias_siesta, (c2_ts2_3A_siesta - c2_ts2_3A_siesta[2]) * 1, 'r^-', label='SIESTA+SMEAGOL 3 Å')
ax_siesta_cp2k_c1_c2[1, 0].set_xlabel('Bias / V')
ax_siesta_cp2k_c1_c2[1, 1].set_xlabel('Bias / V')
ax_siesta_cp2k_c1_c2[0, 0].set_ylabel('C1-TS2 Barrier Change / eV')
ax_siesta_cp2k_c1_c2[1, 0].set_ylabel('C2-TS2 Barrier Change / eV')
ax_siesta_cp2k_c1_c2[1, 0].set_ylim([-0.064, 0.007])
ax_siesta_cp2k_c1_c2[1, 1].set_ylim([-0.064, 0.007])
ax_siesta_cp2k_c1_c2[0, 0].legend(frameon=False)
ax_siesta_cp2k_c1_c2[0, 1].legend(frameon=False)
fig_siesta_cp2k_c1_c2.tight_layout()
fig_siesta_cp2k_c1_c2.savefig('{}/siesta_cp2k_c1_c2.png'.format(save_folder_siesta), dpi=300)
fig_siesta_cp2k_c1_c2.savefig('{}/siesta_cp2k_c1_c2.png'.format(save_folder_cp2k), dpi=300)

print('c1_ts2_5A_cp2k', (c1_ts2_5A_cp2k[2:-1] - c1_ts2_5A_cp2k[4]) * param.hartree_to_ev)
print('c1_ts2_3A_cp2k', (c1_ts2_5A_cp2k[2:-1] - c1_ts2_5A_cp2k[4]) * param.hartree_to_ev)
print('c2_ts2_5A_cp2k', (c2_ts2_5A_cp2k[2:-1] - c2_ts2_5A_cp2k[4]) * param.hartree_to_ev)
print('c2_ts2_3A_cp2k', (c2_ts2_3A_cp2k[2:-1] - c2_ts2_3A_cp2k[4]) * param.hartree_to_ev)

print('c1_ts2_5A_siesta', (c1_ts2_5A_siesta - c1_ts2_5A_siesta[2]) * 1)
print('c1_ts2_3A_siesta', (c1_ts2_3A_siesta - c1_ts2_3A_siesta[2]) * 1)
print('c2_ts2_5A_siesta', (c2_ts2_5A_siesta - c2_ts2_5A_siesta[2]) * 1)
print('c2_ts2_3A_siesta', (c2_ts2_3A_siesta - c2_ts2_3A_siesta[2]) * 1)

if __name__ == "__main__":
    print('Finished.')
    plt.show()
