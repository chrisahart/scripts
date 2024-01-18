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

# /work/e05/e05/cahart/postdoc/transport/cp2k-smeagol-examples/melamine/cp2k-smeagol/cu_size-6-6/iv
# gs_opt = np.array([-19370.900415185347811, -19370.945222938524239, -19370.910016315658140])
# c1_opt = np.array([-19370.879891856813629, -19370.923962915203447, -19370.888524923171644])
# c2_opt = np.array([-19370.869901778714848, -19370.913798194538685, -19370.878386992106243])
# ts1_linear = np.array([-19370.809433596892632, -19370.853919844092161, -19370.818485115851217])
# ts2_guess = np.array([-19370.838092091162252, -19370.880523932148208, -19370.845128781958920])
#
# c1_ts2 = (c1_opt - ts2_guess) * param.hartree_to_ev
# c1_ts2 = c1_ts2 - c1_ts2[1]
# print(c1_ts2)
#
# c2_ts2 = (c2_opt - ts2_guess) * param.hartree_to_ev
# c2_ts2 = c2_ts2 - c2_ts2[1]
# print(c2_ts2)

# /Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/melamine/hx1/au_from-h2
v_bias = np.array([-2, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2])

# without tip
# v_gs_without_tip = np.loadtxt('/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/melamine/hx1/au_from-h2/iv-curve/single-point/gs-notip/energy.out')
v_c1_without_tip = np.loadtxt('/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/melamine/hx1/au_from-h2/iv-curve/single-point/c1-notip/energy.out')
v_c2_without_tip = np.loadtxt('/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/melamine/hx1/au_from-h2/iv-curve/single-point/c2-notip/energy.out')
v_ts2_without_tip = np.loadtxt('/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/melamine/hx1/au_from-h2/iv-curve/single-point/ts2-guess-notip/energy.out')

# with tip 5A
v_gs_5A = np.loadtxt('/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/melamine/hx1/au_from-h2/iv-curve/single-point/gs/energy.out')
v_c1_5A = np.loadtxt('/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/melamine/hx1/au_from-h2/iv-curve/single-point/c1/energy.out')
v_c2_5A = np.loadtxt('/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/melamine/hx1/au_from-h2/iv-curve/single-point/c2/energy.out')
v_ts2_5A = np.loadtxt('/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/melamine/hx1/au_from-h2/iv-curve/single-point/ts2-guess/energy.out')

# with tip 3A
# v_gs_3A = np.loadtxt('/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/melamine/hx1/au_from-h2/iv-curve/single-point/gs-3A/energy.out')
v_c1_3A = np.loadtxt('/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/melamine/hx1/au_from-h2/iv-curve/single-point/c1-3A/energy.out')
v_c2_3A = np.loadtxt('/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/melamine/hx1/au_from-h2/iv-curve/single-point/c2-3A/energy.out')
v_ts2_3A = np.loadtxt('/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/melamine/hx1/au_from-h2/iv-curve/single-point/ts2-3A/energy.out')

v_gs = v_gs_5A
v_c1 = v_c1_5A
v_c2 = v_c2_5A
v_ts2 = v_ts2_5A

print('CP2K C1', (v_gs[2]-v_c1[2])*param.hartree_to_ev)
print('CP2K C2', (v_gs[2]-v_c2[2])*param.hartree_to_ev)
print('CP2K TS2', (v_gs[2]-v_ts2[2])*param.hartree_to_ev)

# v_gs = np.loadtxt('/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/melamine/hx1/au_from-h2/iv-curve/single-point/gs-notip/energy.out')
# v_c1 = np.loadtxt('/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/melamine/hx1/au_from-h2/iv-curve/single-point/c1-notip/energy.out')
# v_c2 = np.loadtxt('/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/melamine/hx1/au_from-h2/iv-curve/single-point/c2-notip/energy.out')
# v_ts2 = np.loadtxt('/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/melamine/hx1/au_from-h2/iv-curve/single-point/ts2-guess-notip/energy.out')

# v_c1 = np.loadtxt('/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/melamine/hx1/au_from-h2/iv-curve/single-point/c1-7A/energy.out')
# v_c2 = np.loadtxt('/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/melamine/hx1/au_from-h2/iv-curve/single-point/c2-7A/energy.out')
# v_ts2 = np.loadtxt('/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/melamine/hx1/au_from-h2/iv-curve/single-point/ts2-guess-7A/energy.out')

v_bias_siesta = np.array([-1.0, -0.5, 0.0, 0.5, 1.0])
v_gs_siesta = np.loadtxt('/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/melamine/young/siesta1-smeagol/au_from-h2/iv-curve/gs_leads-sz-custom_surf-tip-dzp/energy.out')
v_c1_siesta = np.loadtxt('/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/melamine/young/siesta1-smeagol/au_from-h2/iv-curve/c1_leads-sz-custom_surf-tip-dzp/energy.out')
v_c2_siesta = np.loadtxt('/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/melamine/young/siesta1-smeagol/au_from-h2/iv-curve/c2_leads-sz-custom_surf-tip-dzp/energy.out')
v_ts1_siesta = np.loadtxt('/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/melamine/young/siesta1-smeagol/au_from-h2/iv-curve/ts1_leads-sz-custom_surf-tip-dzp_from-cu/energy.out')
v_ts2_siesta = np.loadtxt('/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/melamine/young/siesta1-smeagol/au_from-h2/iv-curve/ts2_leads-sz-custom_surf-tip-dzp_from-cu2/energy.out')

# v_gs_siesta = np.loadtxt('/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/melamine/young/siesta1-smeagol/au_from-h2/iv-curve-3A/gs_leads-sz-custom_surf-tip-dzp_from-cu/energy.out')
# v_c1_siesta = np.loadtxt('/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/melamine/young/siesta1-smeagol/au_from-h2/iv-curve-3A/c1_leads-sz-custom_surf-tip-dzp_from-cu/energy.out')
# v_c2_siesta = np.loadtxt('/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/melamine/young/siesta1-smeagol/au_from-h2/iv-curve-3A/c2_leads-sz-custom_surf-tip-dzp_from-cu/energy.out')
# v_ts2_siesta = np.loadtxt('/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/melamine/young/siesta1-smeagol/au_from-h2/iv-curve-3A/ts2_leads-sz-custom_surf-tip-dzp_from-cu/energy.out')

v_gs_siesta = np.loadtxt('/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/melamine/young/siesta1-smeagol/au_from-h2/iv-curve-no-tip/gs_leads-sz-custom_surf-tip-dzp_from-cu/energy.out')
v_c1_siesta = np.loadtxt('/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/melamine/young/siesta1-smeagol/au_from-h2/iv-curve-no-tip/c1_leads-sz-custom_surf-tip-dzp_from-cu/energy.out')
v_c2_siesta = np.loadtxt('/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/melamine/young/siesta1-smeagol/au_from-h2/iv-curve-no-tip/c2_leads-sz-custom_surf-tip-dzp_from-cu/energy.out')
v_ts2_siesta = np.loadtxt('/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/melamine/young/siesta1-smeagol/au_from-h2/iv-curve-no-tip/ts2_leads-sz-custom_surf-tip-dzp_from-cu/energy.out')

# v_c1_siesta = np.loadtxt('/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/melamine/young/siesta1-smeagol/au_from-h2/iv-curve/c1_leads-sz-custom_surf-tip-dzp-sr/energy.out')
# v_c2_siesta = np.loadtxt('/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/melamine/young/siesta1-smeagol/au_from-h2/iv-curve/c2_leads-sz-custom_surf-tip-dzp-sr/energy.out')
# v_ts1_siesta = np.loadtxt('/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/melamine/young/siesta1-smeagol/au_from-h2/iv-curve/ts1_leads-sz-custom_surf-tip-dzp_from-cu/energy.out')
# v_ts2_siesta = np.loadtxt('/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/melamine/young/siesta1-smeagol/au_from-h2/iv-curve/ts2_leads-sz-custom_surf-tip-dzp-sr_from-cu/energy.out')

print('SIESTA C1', (v_gs_siesta[2]-v_c1_siesta[2])*1)
print('SIESTA C2', (v_gs_siesta[2]-v_c2_siesta[2])*1)
print('SIESTA TS1', (v_gs_siesta[2]-v_ts1_siesta[2])*1)
print('SIESTA TS2', (v_gs_siesta[2]-v_ts2_siesta[2])*1)

# v_c1_siesta = np.loadtxt('/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/melamine/young/siesta1-smeagol/au_from-h2/iv-curve/c1_leads-sz-custom_surf-tip-sz/energy.out')
# v_c2_siesta = np.loadtxt('/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/melamine/young/siesta1-smeagol/au_from-h2/iv-curve/c2_leads-sz-custom_surf-tip-sz/energy.out')
# v_ts2_siesta = np.loadtxt('/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/melamine/young/siesta1-smeagol/au_from-h2/iv-curve/ts2_leads-sz-custom_surf-tip-sz/energy.out')

# v_gs = np.array([-5729.588680033193668, -5729.643348350294218, -5729.614110315752441, -5729.621506063563174,
#                  -5729.623712865300149,
# #                  -5729.621315797976422, -5729.610264784145329, -5729.635332518842915, -5729.567002626033172])
# v_c1 = np.array([-5729.562016429242249, -5729.614947287638643, -5729.585309699680693, -5729.592769185197540,
#                  -5729.595002918055798,
#                  -5729.592639791269903, -5729.581606759504211, -5729.606425269940701, -5729.539703125408778])
# v_c2 = np.array([-5729.548446015055561, -5729.604502117189440, -5729.574869466943710, -5729.582400454502931,
#                  -5729.584537766557332,
#                  -5729.582247295762500, -5729.571170316562529, -5729.595991116643745, -5729.529300386735486])
# v_ts2 = np.array([-5729.505076646280941, -5729.569720795700050, -5729.540118150493981, -5729.547899988428981,
#                   -5729.549960014870521,
#                   -5729.547460663692618, -5729.536249907415367, -5729.561032570403768, -5729.494436481123557])

# without tip
# v_gs = np.array([-5563.504948362586219
# ,-5563.554226945822847
# ,-5563.527268125908449
# ,-5563.533097835539593
# ,-5563.533045286002562
# ,-5563.527588138555984
# ,-5563.498473253614065
# ,-5563.508801815541119
# ,-5563.491095289380610
# ])
# v_c1 = np.array([-5563.474383371479234
# ,-5563.523256205680809
# ,-5563.495832191755653
# ,-5563.502134403438504
# ,-5563.501797604520107
# ,-5563.496341717858741
# ,-5563.470577356697504
# ,-5563.481795163597781
# ,-5563.371056250286529
# ])
# v_c2 = np.array([-5563.463747781705933
# ,-5563.512759217190251
# ,-5563.485353718964689
# ,-5563.491676700714379
# ,-5563.491192605025390
# ,-5563.485938432986586
# ,-5563.460096313257964
# ,-5563.493936963448505
# ,-5563.469291015110684
# ])
# v_ts2 = np.array([-5563.434902609867095
# ,-5563.478362541952265
# ,-5563.450835632496819
# ,-5563.457426166355617
# ,-5563.456749982058682
# ,-5563.451267116628514
# ,-5563.425871900268248
# ,-5563.441081999412745
# ,-5563.342199693333896
# ])

c1_ts2_without_tip_cp2k = []
c1_ts2_3A_cp2k = []
c1_ts2_5A_cp2k = []

c2_ts2_without_tip_cp2k = []
c2_ts2_3A_cp2k = []
c2_ts2_5A_cp2k = []

for i in range(v_bias.shape[0] - 1):
    c1_ts2_without_tip_cp2k.append(v_c1_without_tip[i] - v_ts2_without_tip[i])
    c1_ts2_3A_cp2k.append(v_c1_3A[i] - v_ts2_3A[i])
    c1_ts2_5A_cp2k.append(v_c1_5A[i] - v_ts2_5A[i])
    c2_ts2_without_tip_cp2k.append(v_c2_without_tip[i] - v_ts2_without_tip[i])
    c2_ts2_3A_cp2k.append(v_c2_3A[i] - v_ts2_3A[i])
    c2_ts2_5A_cp2k.append(v_c2_5A[i] - v_ts2_5A[i])

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


# fig_cube_z, ax_cube_z = plt.subplots(rows, cols, sharex='col', sharey='row', figsize=(6, 8))
# if plot_fermi: ax_cube_z[0].axhline(y=fermi_dft, color='grey', linestyle='--', label='DFT Fermi energy', alpha=0.5)
# if draw_mirror: ax_cube_z[0].plot(energy_grid_hartree_em1[:mid_index]+mid_pos_grid, np.flip(average_hartree_em2[plot_folder][:mid_index]), 'm--')
# if draw_mirror: ax_cube_z[0].plot(energy_grid_hartree_em1[:mid_index]+mid_pos_grid, np.flip(average_hartree_em1[plot_folder][:mid_index]), 'm--')
# if plot_labels[1]: ax_cube_z[0].plot(energy_grid_hartree_dft, average_hartree_dft[plot_folder], '-', color=plot_color[1], label=labels[1])
# if plot_labels[3]: ax_cube_z[0].plot(energy_grid_hartree_em2, average_hartree_em2[plot_folder], '-', color=plot_color[3], label=labels[3])
# if plot_labels[2]: ax_cube_z[0].plot(energy_grid_hartree_em1, average_hartree_em1[plot_folder], '-', color=plot_color[2], label=labels[2])
# if plot_labels[0]: ax_cube_z[0].plot(energy_grid_hartree_bulk, average_hartree_bulk[plot_folder], '-', color=plot_color[0], label=labels[0])
# if draw_markers: ax_cube_z[0].plot(markers, markers*0, 'o', color='orange', fillstyle='none')
# ax_cube_z[0].set_xlim([xlim[0], xlim[1]])
# ax_cube_z[0].legend(frameon=False)
# ax_cube_z[0].set_ylabel('Hartree potential / eV')


rows, cols = 2, 1
fig_cp2k_c1_c2, ax_cp2k_c1_c2 = plt.subplots(rows, cols, sharex='col', sharey='row', figsize=(6, 8))
ax_cp2k_c1_c2[0].plot(v_bias[2:-2], (c1_ts2_5A_cp2k[2:-1] - c1_ts2_5A_cp2k[4]) * param.hartree_to_ev, 'ms-', label='5 A')
ax_cp2k_c1_c2[0].plot(v_bias[2:-2], (c1_ts2_3A_cp2k[2:-1] - c1_ts2_3A_cp2k[4]) * param.hartree_to_ev, 'r^-', label='3 A')
ax_cp2k_c1_c2[1].plot(v_bias[2:-2], (c2_ts2_5A_cp2k[2:-1] - c2_ts2_5A_cp2k[4]) * param.hartree_to_ev, 'ms-', label='5 A')
ax_cp2k_c1_c2[1].plot(v_bias[2:-2], (c2_ts2_3A_cp2k[2:-1] - c2_ts2_3A_cp2k[4]) * param.hartree_to_ev, 'r^-', label='3 A')
# ax_cp2k_c1_c2.plot(v_bias[2:-2], (c2_ts2_without_tip_cp2k[2:-1] - c2_ts2_without_tip_cp2k[4]) * param.hartree_to_ev, 'bx-', label='Without tip')
# ax_cp2k_c1_c2.legend(frameon=False)
ax_cp2k_c1_c2[1].set_xlabel('Bias / V')
ax_cp2k_c1_c2[0].set_ylabel('Barrier Change / eV')
ax_cp2k_c1_c2[1].set_ylabel('Barrier Change / eV')
ax_cp2k_c1_c2[0].set_ylim([-0.046, 0.007])
ax_cp2k_c1_c2[1].set_ylim([-0.046, 0.007])
fig_cp2k_c1_c2.tight_layout()

# fig_cp2k, ax_cp2k = plt.subplots()
# ax_cp2k.plot(v_bias[2:-2], (c1_ts2[2:-1] - c1_ts2[4])*param.hartree_to_ev, 'rx-', label='C1 - TS2')
# ax_cp2k.plot(v_bias[2:-2], (c2_ts2[2:-1] - c2_ts2[4])*param.hartree_to_ev, 'gx-', label='C2 - TS2')
# ax_cp2k.plot(v_bias_siesta , (c1_ts2_siesta  - c1_ts2_siesta [2])*1, 'rx-', label='C1 - TS2')
# ax_cp2k.plot(v_bias_siesta , (c2_ts2_siesta  - c2_ts2_siesta [2])*1, 'gx-', label='C2 - TS2')
# ax_cp2k.legend(frameon=False)
# ax_cp2k.set_xlabel('Bias / V')
# ax_cp2k.set_ylabel('Barrier Change / eV')
# fig_cp2k.tight_layout()

# fig_siesta, ax_siesta = plt.subplots()
# ax_siesta.plot(v_bias[2:-2], (c1_ts2[2:-1] - c1_ts2[4])*param.hartree_to_ev, 'rx-', label='C1 - TS2')
# ax_siesta.plot(v_bias[2:-2], (c2_ts2[2:-1] - c2_ts2[4])*param.hartree_to_ev, 'gx-', label='C2 - TS2')
# ax_siesta.plot(v_bias_siesta , (c1_ts2_siesta  - c1_ts2_siesta [2])*1, 'rx-', label='C1 - TS2')
# ax_siesta.plot(v_bias_siesta , (c2_ts2_siesta  - c2_ts2_siesta [2])*1, 'gx-', label='C2 - TS2')
# ax_siesta.legend(frameon=False)
# ax_siesta.set_xlabel('Bias / V')
# ax_siesta.set_ylabel('Barrier Change / eV')
# fig_siesta.tight_layout()

# for i in range(len(folder1)):
#     fig_charge.savefig('{}/charge_z.png'.format(folder1[i]), dpi=300)

if __name__ == "__main__":
    print('Finished.')
    plt.show()
