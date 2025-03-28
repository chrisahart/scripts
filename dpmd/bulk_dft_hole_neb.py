import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import MDAnalysis as mda
from MDAnalysis.analysis import distances
from general import parameters as param

"""
    Plot energy and forces for bulk hematite
"""


def read_energy(folder, filename):
    """
        Return CP2K MD .ener file as re-structured Numpy array.
    """

    files = ['{}/{}'.format(folder, filename)]
    cols = ['Step', 'Time', 'E_kin', 'Temp', 'E_pot', 'E_tot', 'Time_per_step']
    file_energy = pd.read_csv(files[0], delim_whitespace=True, names=cols, skiprows=[0])

    # Load energy data from Pandas database
    energy_kinetic = file_energy['E_kin'].values
    energy_potential = file_energy['E_pot'].values
    energy_total = file_energy['E_tot'].values
    temperature = file_energy['Temp'].values
    time = file_energy['Time'].values
    step = file_energy['Step'].values
    time_per_step = file_energy['Time_per_step'].values

    return energy_kinetic, energy_potential, energy_total, temperature, time, time_per_step, step


# def read_hirsh(folder, filename):
#     """
#     Read Hirshfeld
#     """
#
#     cols_hirsh = ['Atom', 'Element', 'Kind', 'Ref Charge', 'Pop 1', 'Pop 2', 'Spin', 'Charge']
#     data_hirsh = pd.read_csv('{}{}'.format(folder, filename), names=cols_hirsh, delim_whitespace=True)
#     species = data_hirsh['Element']
#
#     return data_hirsh, species


def read_hirsh(folder, filename):
    """
    Read Hirshfeld analysis from CP2K output file
    """

    # Read number of atoms and labels from .xyz file
    files = ['{}/{}'.format(folder, filename)]
    cols = ['Atom', 'Element', 'Kind', 'Ref Charge', 'Pop 1', 'Pop 2', 'Spin', 'Charge', 'A', 'B']
    file_spec1 = pd.read_csv(files[0], names=cols, delim_whitespace=True, skiprows=5)

    file_spec1 = file_spec1.drop(columns=['A'])
    file_spec1 = file_spec1.drop(columns=['B'])
    file_spec1 = file_spec1.drop(columns=['Element'])

    file_spec1 = file_spec1.apply(pd.to_numeric, errors='coerce')
    file_spec1 = file_spec1.dropna(axis='rows', thresh=2)
    file_spec1 = file_spec1.dropna(axis='columns', thresh=1)
    file_spec1 = file_spec1.reset_index(drop=True)

    return file_spec1


def func_metric(a, b, c):
    index = [0, 1, 2, 3, 4, 5]
    metric = np.average(a.flat[index])
    return metric


skip = 2
atoms = 120
value = 'Spin'

# MD
# run = '400K'
# folder_4 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/chris_phd/bulk/hole/hse/400k'
# folder_4 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/chris_phd/bulk/hole/hse/400k-b'
# folder_4 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/chris_phd/bulk/hole/hse/400k-f'
# folder_4 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/chris_phd/bulk/hole/hse/hops-0'
# folder_4 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/chris_phd/bulk/hole/hse/330k'

# NEB atom 0-1: highest coupling 0-a-1-a
# run = '0K'
# folder_4 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/archer/bulk/221_supercell/reftraj/neb/atom-0-1/0-a-1-a'
# folder_4 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/archer/bulk/221_supercell/reftraj/neb/atom-0-1/0-a-1-b'
# folder_4 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/archer/bulk/221_supercell/reftraj/neb/atom-0-1/0-a-1-c'
# folder_4 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/archer/bulk/221_supercell/reftraj/neb/atom-0-1/0-b-1-a'
# folder_4 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/archer/bulk/221_supercell/reftraj/neb/atom-0-1/0-b-1-b'
# folder_4 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/archer/bulk/221_supercell/reftraj/neb/atom-0-1/0-b-1-c'
# folder_4 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/archer/bulk/221_supercell/reftraj/neb/atom-0-1/0-c-1-a'
# folder_4 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/archer/bulk/221_supercell/reftraj/neb/atom-0-1/0-c-1-b'
# folder_4 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/archer/bulk/221_supercell/reftraj/neb/atom-0-1/0-c-1-c'
# atom 0-1
# atom 13: 50, 67, 61, 88, 95, 108
# atom 38: 106, 109, 63, 104, 108, 95
# o_plot = np.array([95, 108]) - 1
# o_plot = np.array([50, 67, 61, 88, 95, 108, 106, 109, 63, 104]) - 1

# NEB atom 0-3: highest coupling 0-b-3-c
# run = '0K'
# folder_4 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/archer/bulk/221_supercell/reftraj/neb/atom-0-3/0-a-3-a'
# folder_4 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/archer/bulk/221_supercell/reftraj/neb/atom-0-3/0-a-3-b'
# folder_4 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/archer/bulk/221_supercell/reftraj/neb/atom-0-3/0-a-3-c'
# folder_4 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/archer/bulk/221_supercell/reftraj/neb/atom-0-3/0-b-3-a'
# folder_4 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/archer/bulk/221_supercell/reftraj/neb/atom-0-3/0-b-3-b'
# folder_4 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/archer/bulk/221_supercell/reftraj/neb/atom-0-3/0-b-3-c'
# folder_4 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/archer/bulk/221_supercell/reftraj/neb/atom-0-3/0-c-3-a'
# folder_4 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/archer/bulk/221_supercell/reftraj/neb/atom-0-3/0-c-3-b'
# folder_4 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/archer/bulk/221_supercell/reftraj/neb/atom-0-3/0-c-3-c'
# atom 13: 50, 67, 61, 88, 95, 108
# atom 6: 54, 58, 88, 71, 67, 69
# o_plot = np.array([67, 88]) - 1

# REFTRAJ NEB all
run = '0K'
# folder_4 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/archer/bulk/221_supercell/reftraj/neb/atom-0-1/concat/all'
# o_plot = np.array([95, 108]) - 1
folder_4 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/archer/bulk/221_supercell/reftraj/neb/atom-0-3/database/geo-opt-all'
o_plot = np.array([67, 88]) - 1
# folder_4 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/archer/bulk/441_supercell/reftraj/neb/atom-0-1/database/geo-opt-441'
# o_plot = np.array([95, 108]) - 1

folder_save = folder_4
polaron_index_fe = {'0K': [13-1, 38-1],
                    '330K': [13-1, 27-1, 46-1, 5-1],
                    '400K': [13-1, 27-1, 38-1, 1-1, 28-1, 14-1, 29-1, 14-1, 42-1]}
labels = '0'
plot_color = 'r', 'b', 'g', 'c', 'm', 'orange', 'y', 'peru','yellowgreen', 'lightgreen'
index_fe_4 = polaron_index_fe[run]
energy_kinetic1_4, energy_potential1_4, energy_total1_4, temperature1_4, steps1_4, time_per_step1_4, steps1_4 = read_energy(folder_4, 'hematite-1.ener')
file_spec1_4 = read_hirsh(folder_4, '/hematite-charges-1-clean.hirshfeld')
print(file_spec1_4)
num_data1_4 = energy_kinetic1_4.shape[0]
topology_file = '{}/topology.xyz'.format(folder_4)
trajectory_file = '{}/hematite-pos-1.xyz'.format(folder_4)
timestep = 0.5
time_array = np.linspace(0, int(steps1_4.shape[0]*timestep), num=steps1_4.shape[0])
# steps1_4 = steps1_4 / timestep

transition_time = np.array([460, 117.5, 315, 656, 721, 1283])
zoom = False
transition_time_plot = 0
axis_lim_x_zoom = np.array([transition_time[transition_time_plot]-40, transition_time[transition_time_plot]+40])

# check sizes are consistent
print(int(file_spec1_4.shape[0]/atoms))
print(num_data1_4)
from ase.io import read, write
# ase_forces = read('{}/{}'.format(folder_4, 'hematite-frc-1.xyz'), format='xyz', index=':')
# print(all_frames)
# print('Number of atoms', len(all_frames))
# print('Number of frames', len(all_frames))
# extracted_frames = all_frames[1699:2900]
# write('{}/{}'.format(folder_4, 'hematite-frc-1-test.xyz'), extracted_frames, format='xyz')

# Plotting
draw_polaron = False
draw_legend = False
polaron_size = 4
polaron_alpha = 1
ylim_1 = [-4.0, -3.0]
ylim_2 = [2.16, 1.89]
ylim_2 = [2.10, 1.90]
ylim_3 = [-31.6, -30.4]
strength_limit = -31.0

# System
box_size = [10.071, 10.071, 13.747, 90, 90, 120]
h_all = np.NaN
water = np.NaN
fe_beta = np.array([1, 2, 3, 4, 5, 6, 13, 14, 15, 16, 17, 18, 25, 26, 27, 28, 29, 30, 37, 38, 41, 42, 45, 46]) - 1
fe_alpha = np.array([7, 8, 9, 10, 11, 12, 19, 20, 21, 22, 23, 24, 31, 32, 33, 34, 35, 36, 39, 40, 43, 44, 47, 48]) - 1
# fe_b = np.array([14, 16, 18, 42, 27, 45, 25, 29]) - 1
fe_b = np.array([27, 45, 18, 14, 25, 29, 42, 16]) - 1
# fe_d = np.array([6, 2, 13, 17, 38, 4, 15, 41]) - 1
fe_d = np.array([2, 6, 17, 13, 4, 38, 41, 15]) - 1
# fe_f = np.array([46, 28, 5, 1, 30, 26, 37, 3]) - 1
fe_f = np.array([28, 46, 1, 5, 26, 30, 3, 37]) - 1
o_neighbours = np.array([116, 76, 55, 89, 49, 102]) - 1
o_all = np.linspace(start=49, stop=120, num=120-49+1, dtype=int) - 1
num_species = np.array([len(o_all), len(fe_alpha), len(fe_beta)])
plotting_colors = ['r', 'g', 'b', 'm', 'grey', 'orange', 'brown', 'hotpink'] * 100

# Plot all iron spin 1
fig_spin2, ax_spin2 = plt.subplots()
x_end = steps1_4[-1]-steps1_4[0]
if draw_polaron:
    temp4 = np.zeros(num_data1_4)
    for j in range(len(index_fe_4)):
        for n in range(num_data1_4):
            temp4[n] = (file_spec1_4.loc[atoms * n + index_fe_4[j], 'Spin'])
        ax_spin2.plot(steps1_4-steps1_4[0], temp4, '--', color=plot_color[j], linewidth=polaron_size, alpha=polaron_alpha)
    ax_spin2.plot(steps1_4[0], temp4[0], '--', label='Polaron', color=plot_color[0], linewidth=polaron_size, alpha=polaron_alpha)
temp1 = np.zeros((8, num_data1_4))
temp2 = np.zeros((8, num_data1_4))
temp3 = np.zeros((8, num_data1_4))
temp4 = np.zeros((o_all.shape[0], num_data1_4))
temp5 = np.zeros((1, num_data1_4))
temp6 = np.zeros((o_neighbours.shape[0], num_data1_4))
for j in range(len(fe_b)):
    for n in range(num_data1_4):
        temp1[j, n] = (file_spec1_4.loc[atoms * n + fe_b[j], 'Spin'])
        temp2[j, n] = (file_spec1_4.loc[atoms * n + fe_d[j], 'Spin'])
        temp3[j, n] = (file_spec1_4.loc[atoms * n + fe_f[j], 'Spin'])
    ax_spin2.plot(steps1_4 - steps1_4[0], temp1[j, :], 'rx-')
    ax_spin2.plot(steps1_4 - steps1_4[0], temp2[j, :], 'gx-')
    ax_spin2.plot(steps1_4 - steps1_4[0], temp3[j, :], 'bx-')
print(temp1)
print(num_data1_4)
print(steps1_4)
print(steps1_4[0])
print(steps1_4 - steps1_4[0])
ax_spin2.plot(steps1_4[0] - steps1_4[0], temp1[0, 0], 'r-', label='Fe, B')
ax_spin2.plot(steps1_4[0] - steps1_4[0], temp2[0, 0], 'g-', label='Fe, D')
ax_spin2.plot(steps1_4[0] - steps1_4[0], temp3[0, 0], 'b-', label='Fe, F')
# ax_spin2.plot(steps1_4 - steps1_4[0], np.sum(temp1, axis=0)/8, 'r--')
# ax_spin2.plot(steps1_4 - steps1_4[0], np.sum(temp2, axis=0)/8, 'g--')
# ax_spin2.plot(steps1_4 - steps1_4[0], np.sum(temp3, axis=0)/8, 'b--')
if draw_legend: ax_spin2.legend(frameon=True)
# ax_spin2.plot([0, x_end], [-3.29, -3.29], '--', label='Bulk', color='grey')
# ax_spin2.set_xlabel('Timestep')
ax_spin2.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
ax_spin2.set_xlabel('Reaction coordinate')
# ax_spin2.set_xlabel('Time / fs')
ax_spin2.set_ylabel('Spin moment')
ax_spin2.set_ylim(ylim_1)
ax_spin2.set_xlim([0, x_end])
# ax_spin2.set_xlim([1500, 1600])
fig_spin2.tight_layout()
fig_spin2.savefig('{}/spin_all_{}.png'.format(folder_save, run), dpi=300)

# # Plot all iron spin layer
# fig_spin3, ax_spin3 = plt.subplots(figsize=(12, 6))
fig_spin3, ax_spin3 = plt.subplots()
x_end = steps1_4[-1]-steps1_4[0]
for j in range(len(fe_b)):
    # ax_spin3.plot(time_array - time_array[0], temp1[j, :], '-', color=plotting_colors[j], label='Fe {}'.format(j+1))
    ax_spin3.plot(time_array - time_array[0], temp2[j, :], 'x-', color=plotting_colors[j], label='Fe {}'.format(j+1))
    # ax_spin3.plot(time_array - time_array[0], temp3[j, :], '-', color=plotting_colors[j], label='Fe {}'.format(j+1))
if draw_legend: ax_spin3.legend(frameon=True)
ax_spin3.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
ax_spin3.set_xlabel('Reaction coordinate')
# ax_spin3.set_xlabel('Timestep')
# ax_spin3.set_xlabel('Time / fs')
ax_spin3.set_ylabel('Spin moment')
ax_spin3.set_ylim(ylim_1)
ax_spin3.set_xlim([0, x_end/2])
# ax_spin3.set_xlim([1500, 1600])
fig_spin3.tight_layout()
fig_spin3.savefig('{}/spin_layer_{}.png'.format(folder_save, run), dpi=300)
fig_spin3.tight_layout()
if zoom:
    ax_spin3.set_xlim(axis_lim_x_zoom[0], axis_lim_x_zoom[1])
    fig_spin3.tight_layout()
    fig_spin3.savefig('{}/spin_layer_zoom_{}.png'.format(folder_save, run), dpi=300)

# Plot all o spin
fig_spin_o, ax_spin_o = plt.subplots()
x_end = steps1_4[-1]-steps1_4[0]
temp1 = np.zeros((len(o_all), num_data1_4))
for j in range(len(o_plot)):
    for n in range(num_data1_4):
        temp1[j, n] = (file_spec1_4.loc[atoms * n + o_plot[j], 'Spin'])
    ax_spin_o.plot(steps1_4 - steps1_4[0], temp1[j, :], 'x-', color=plotting_colors[j], label=o_plot[j]+1)
# ax_spin_o.plot(steps1_4[0] - steps1_4[0], temp1[0, 0], '-', color=plotting_colors[j], label='O')
# ax_spin_o.legend(frameon=True)
if draw_legend: ax_spin_o.legend(frameon=True)
# ax_spin_o.set_xlabel('Timestep')
ax_spin_o.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
ax_spin_o.set_xlabel('Reaction coordinate')
# ax_spin_o.set_xlabel('Time / fs')
ax_spin_o.set_ylabel('Spin moment')
ax_spin_o.set_ylim([-0.06, 0.14])
# ax_spin_o.set_ylim([-0.05, 0.1])
# ax_spin_o.set_ylim(ylim_1)
ax_spin_o.set_xlim([0, x_end])
# ax_spin_o.set_xlim([1500, 1600])
fig_spin_o.tight_layout()
fig_spin_o.savefig('{}/spin_o_bridging_{}.png'.format(folder_save, run), dpi=300)

# Plot energy
fig_energy, ax_energy = plt.subplots()
# ax_energy.plot(steps1_4 - steps1_4[0], (energy_kinetic1_4 - energy_kinetic1_4[0]) * param.hartree_to_ev, 'r-', label='Kinetic')
ax_energy.plot(steps1_4 - steps1_4[0], (energy_potential1_4 - energy_potential1_4[0]) * param.hartree_to_ev, 'kx-', label='Potential')
# ax_energy.plot(steps1_4 - steps1_4[0], (energy_total1_4 - energy_total1_4[0]) * param.hartree_to_ev, 'b-', label='Total')
ax_energy.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
ax_energy.set_xlabel('Reaction coordinate')
# ax_energy.set_xlabel('Timestep')
ax_energy.set_ylabel('Energy / eV')
# ax_energy.legend(frameon=True)
fig_energy.tight_layout()
fig_energy.savefig('{}/energy{}.png'.format(folder_save, run), dpi=300)

# Plot force x
# dimension = 0
# fig_force_x, ax_force_x = plt.subplots()
# for j in range(len(fe_b)):
#     for i, atoms in enumerate(ase_forces):
#         temp1[j, i] = atoms[fe_b[j]].position[dimension]
#         temp2[j, i] = atoms[fe_d[j]].position[dimension]
#         temp3[j, i] = atoms[fe_f[j]].position[dimension]
#         temp4[j, i] = atoms[o_all[j]].position[dimension]
#     # ax_force_x.plot(steps1_4 - steps1_4[0], temp1[j, :], 'r-')
#     # ax_force_x.plot(steps1_4 - steps1_4[0], temp2[j, :], 'g-')
#     # ax_force_x.plot(steps1_4 - steps1_4[0], temp3[j, :], 'b-')
#     ax_force_x.plot(steps1_4 - steps1_4[0], temp4[j, :], 'm-')
# # ax_force_x.plot(steps1_4[0] - steps1_4[0], temp1[0, 0], 'r-', label='Fe, B')
# # ax_force_x.plot(steps1_4[0] - steps1_4[0], temp2[0, 0], 'g-', label='Fe, D')
# # ax_force_x.plot(steps1_4[0] - steps1_4[0], temp3[0, 0], 'b-', label='Fe, F')
# ax_force_x.plot(steps1_4[0] - steps1_4[0], temp4[0, 0], 'm-', label='O')
# ax_force_x.set_xlabel('Timestep')
# ax_force_x.set_ylabel('Force x / au')
# ax_force_x.legend(frameon=True)
# fig_force_x.tight_layout()
# fig_force_x.savefig('{}/force{}.png'.format(folder_save, run), dpi=300)

# Plot force y
# dimension = 1
# fig_force_y, ax_force_y = plt.subplots()
# for j in range(len(fe_b)):
#     for i, atoms in enumerate(ase_forces):
#         temp1[j, i] = atoms[fe_b[j]].position[dimension]
#         temp2[j, i] = atoms[fe_d[j]].position[dimension]
#         temp3[j, i] = atoms[fe_f[j]].position[dimension]
#         temp4[j, i] = atoms[o_all[j]].position[dimension]
#     # ax_force_y.plot(steps1_4 - steps1_4[0], temp1[j, :], 'r-')
#     # ax_force_y.plot(steps1_4 - steps1_4[0], temp2[j, :], 'g-')
#     # ax_force_y.plot(steps1_4 - steps1_4[0], temp3[j, :], 'b-')
#     ax_force_y.plot(steps1_4 - steps1_4[0], temp4[j, :], 'm-')
# # ax_force_y.plot(steps1_4[0] - steps1_4[0], temp1[0, 0], 'r-', label='Fe, B')
# # ax_force_y.plot(steps1_4[0] - steps1_4[0], temp2[0, 0], 'g-', label='Fe, D')
# # ax_force_y.plot(steps1_4[0] - steps1_4[0], temp3[0, 0], 'b-', label='Fe, F')
# ax_force_y.plot(steps1_4[0] - steps1_4[0], temp4[0, 0], 'm-', label='O')
# ax_force_x.set_xlabel('Timestep')
# ax_force_y.set_ylabel('Force y / au')
# ax_force_y.legend(frameon=True)
# fig_force_y.tight_layout()
# fig_force_y.savefig('{}/force{}.png'.format(folder_save, run), dpi=300)

# Plot force z
# dimension = 2
# fig_force_z, ax_force_z = plt.subplots()
# for j in range(len(fe_b)):
#     for i, atoms in enumerate(ase_forces):
#         temp1[j, i] = atoms[fe_b[j]].position[dimension]
#         temp2[j, i] = atoms[fe_d[j]].position[dimension]
#         temp3[j, i] = atoms[fe_f[j]].position[dimension]
#         temp4[j, i] = atoms[o_all[j]].position[dimension]
#     # ax_force_z.plot(steps1_4 - steps1_4[0], temp1[j, :], 'r-')
#     # ax_force_z.plot(steps1_4 - steps1_4[0], temp2[j, :], 'g-')
#     # ax_force_z.plot(steps1_4 - steps1_4[0], temp3[j, :], 'b-')
#     ax_force_z.plot(steps1_4 - steps1_4[0], temp4[j, :], 'm-')
# # ax_force_z.plot(steps1_4[0] - steps1_4[0], temp1[0, 0], 'r-', label='Fe, B')
# # ax_force_z.plot(steps1_4[0] - steps1_4[0], temp2[0, 0], 'g-', label='Fe, D')
# # ax_force_z.plot(steps1_4[0] - steps1_4[0], temp3[0, 0], 'b-', label='Fe, F')
# ax_force_z.plot(steps1_4[0] - steps1_4[0], temp4[0, 0], 'm-', label='O')
# ax_force_x.set_xlabel('Timestep')
# ax_force_z.set_ylabel('Force z / au')
# ax_force_z.legend(frameon=True)
# fig_force_z.tight_layout()
# fig_force_z.savefig('{}/force{}.png'.format(folder_save, run), dpi=300)

# Plot force
# fig_force_z, ax_force_z = plt.subplots()
# for j in range(len(fe_b)):
#     for i, atoms in enumerate(ase_forces):
#         temp1[j, i] = np.linalg.norm((atoms[fe_b[j]].position[0], atoms[fe_b[j]].position[1], atoms[fe_b[j]].position[2]))
#         temp2[j, i] = np.linalg.norm((atoms[fe_d[j]].position[0], atoms[fe_d[j]].position[1], atoms[fe_d[j]].position[2]))
#         temp3[j, i] = np.linalg.norm((atoms[fe_f[j]].position[0], atoms[fe_f[j]].position[1], atoms[fe_f[j]].position[2]))
#         # ax_force_z.plot(steps1_4 - steps1_4[0], temp1[j, :], 'r-')
#         # ax_force_z.plot(steps1_4 - steps1_4[0], temp2[j, :], 'g-')
#     # ax_force_z.plot(steps1_4 - steps1_4[0], temp3[j, :], 'b-')
# for j in range(len(o_all)):
#     for i, atoms in enumerate(ase_forces):
#         temp4[j, i] = np.linalg.norm((atoms[o_all[j]].position[0], atoms[o_all[j]].position[1], atoms[o_all[j]].position[2]))
#     # ax_force_z.plot(steps1_4 - steps1_4[0], temp4[j, :], 'm-')
# for j in range(len(o_neighbours)):
#     for i, atoms in enumerate(ase_forces):
#         temp6[j, i] = np.linalg.norm((atoms[o_neighbours[j]].position[0], atoms[o_neighbours[j]].position[1], atoms[o_neighbours[j]].position[2]))
# temp6_mean = np.zeros(num_data1_4)
# for i, atoms in enumerate(ase_forces):
#         temp6_mean[i] = np.mean(temp6[:, i])
# ax_force_z.plot(steps1_4 - steps1_4[0], temp6_mean, 'c-')
# for i, atoms in enumerate(ase_forces):
#     temp5[0, i] = np.linalg.norm((atoms[0].position[0], atoms[0].position[1], atoms[0].position[2]))
# ax_force_z.plot(steps1_4 - steps1_4[0], temp5[0, :], 'k-')
# # ax_force_z.plot(steps1_4[0] - steps1_4[0], temp1[0, 0], 'r-', label='Fe, B')
# # ax_force_z.plot(steps1_4[0] - steps1_4[0], temp2[0, 0], 'g-', label='Fe, D')
# ax_force_z.plot(steps1_4[0] - steps1_4[0], temp3[0, 0], 'b-', label='Fe, F')
# ax_force_z.plot(steps1_4[0] - steps1_4[0], temp4[0, 0], 'm-', label='O')
# ax_force_z.plot(steps1_4[0] - steps1_4[0], temp6_mean[0], 'c-', label='O polaron')
# ax_force_z.plot(steps1_4[0] - steps1_4[0], temp5[0, 0], 'k-', label='Fe polaron')
# ax_force_x.set_xlabel('Timestep')
# ax_force_z.set_ylabel('Force norm / au')
# ax_force_z.legend(frameon=True)
# fig_force_z.tight_layout()
# fig_force_z.savefig('{}/force{}.png'.format(folder_save, run), dpi=300)

# # Get indexes for mdanalysis
# species_numpy = species1_4[3:atoms+3].to_numpy()
# species_numpy_fe = np.where(species_numpy == 'Fe')
# fe_only_b = np.zeros(fe_b.shape[0])
# fe_only_d = np.zeros(fe_d.shape[0])
# fe_only_f = np.zeros(fe_f.shape[0])
# for i in range(fe_b.shape[0]):
#     fe_only_b[i] = np.count_nonzero(species_numpy[:fe_b[i]] == 'Fe_b')
#     fe_only_d[i] = np.count_nonzero(species_numpy[:fe_d[i]] == 'Fe_d')
#     fe_only_f[i] = np.count_nonzero(species_numpy[:fe_f[i]] == 'Fe_f')

# Setup md analysis environment
universe = mda.Universe(topology_file, trajectory_file)
atoms_fe = universe.select_atoms('name Fe')
atoms_o = universe.select_atoms('name O')
dist_arr = distances.distance_array(atoms_fe.positions, atoms_o.positions, box=box_size)
bond_lengths_time = np.zeros((len(universe.trajectory), len(atoms_fe), len(atoms_o)))
bond_lengths_mean_1 = np.zeros((len(universe.trajectory)))
bond_lengths_mean_2 = np.zeros((len(universe.trajectory)))
for ts in universe.trajectory:
    frame = universe.trajectory.frame
    bond_lengths_time[frame] = distances.distance_array(atoms_fe.positions, atoms_o.positions, box=box_size)
bond_lengths_time_sorted = np.zeros((len(universe.trajectory), len(atoms_fe), len(atoms_o)))
bond_lengths_time_sorted_mean = np.zeros((len(universe.trajectory), len(atoms_fe)))
bond_lengths_time_sorted_index = np.zeros((len(universe.trajectory), len(atoms_fe), 6))
for i in range(len(atoms_fe)):
    for j in range(len(universe.trajectory)):
        bond_lengths_time_sorted_index[j, i, :] = np.argsort(bond_lengths_time[j, i])[0:6]
        bond_lengths_time_sorted_mean[j, i] = np.mean(np.sort(bond_lengths_time[j, i])[0:6])
print(bond_lengths_time_sorted_mean[0, 13-1])
print(bond_lengths_time_sorted_index[0, 13-1, :])
# # Plot  metric (all)
time_plot = np.linspace(start=0, stop=len(universe.trajectory) * timestep, num=len(universe.trajectory))
metric = np.zeros((len(atoms_fe), len(universe.trajectory)))
fig_4, ax_4 = plt.subplots()
for i in range(len(fe_b)):
    # ax_4.plot(time_array - time_array[0], bond_lengths_time_sorted_mean[:, fe_b[i]], '-', color=plotting_colors[i], label='Fe {}'.format(i + 1))
    ax_4.plot(time_array - time_array[0], bond_lengths_time_sorted_mean[:, fe_d[i]], 'x-', color=plotting_colors[i], label='Fe {}'.format(i + 1))
    # ax_4.plot(time_array - time_array[0], bond_lengths_time_sorted_mean[:, fe_f[i]], '-', color=plotting_colors[i], label='Fe {}'.format(i + 1))
# ax_4.get_xaxis().set_visible(False)
ax_4.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
ax_4.set_xlabel('Reaction coordinate')
# ax_4.set_xlabel('Time / fs')
# ax_4.set_xlabel('Timestep')
ax_4.set_ylabel('Average Fe-O bond length / A')
# ax_4.set_xlim([0, len(universe.trajectory)])
ax_4.set_xlim([0, x_end/2])
ax_4.set_ylim(ylim_2)
if draw_legend: ax_4.legend(frameon=False)
fig_4.savefig('{}/bond_lengths_fe_o_{}.png'.format(folder_save, run), dpi=300)
fig_4.tight_layout()
if zoom:
    ax_4.set_xlim(axis_lim_x_zoom[0], axis_lim_x_zoom[1])
fig_4.tight_layout()
fig_4.savefig('{}/bond_lengths_fe_o_zoom_{}_{}.png'.format(folder_save, transition_time_plot, run), dpi=300)

# # Plot  metric (all)
fig_5, ax_5 = plt.subplots()
for i in range(len(fe_b)):
    ax_5.plot(time_array - time_array[0], bond_lengths_time_sorted_mean[:, fe_b[i]], 'rx-', label='Fe b')
    ax_5.plot(time_array - time_array[0], bond_lengths_time_sorted_mean[:, fe_d[i]], 'gx-', label='Fe d')
    ax_5.plot(time_array - time_array[0], bond_lengths_time_sorted_mean[:, fe_f[i]], 'bx-', label='Fe f')
ax_5.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
ax_5.set_xlabel('Reaction coordinate')
# ax_5.set_xlabel('Time / fs')
# ax_5.set_xlabel('Timestep')
ax_5.set_ylabel('Average Fe-O bond length / A')
# ax_5.set_xlim([0, len(universe.trajectory)])
ax_5.set_xlim([0, x_end/2])
ax_5.set_ylim(ylim_2)
if draw_legend: ax_5.legend(frameon=False)
fig_5.savefig('{}/bond_lengths_all_fe_o_{}.png'.format(folder_save, run), dpi=300)
fig_5.tight_layout()
if zoom:
    ax_5.set_xlim(axis_lim_x_zoom[0], axis_lim_x_zoom[1])
fig_5.tight_layout()
fig_5.savefig('{}/bond_lengths_all_fe_o_zoom_{}_{}.png'.format(folder_save, transition_time_plot, run), dpi=300)

# Plot  metric (color coded by layer)
# time_plot = np.linspace(start=0, stop=len(universe.trajectory)*timestep, num=len(universe.trajectory))
# metric = np.zeros((len(universe.trajectory)))
# fig_6, ax_6 = plt.subplots()
# temp1 = np.zeros((len(fe_only_b), len(universe.trajectory)))
# temp2 = np.zeros((len(fe_only_d), len(universe.trajectory)))
# temp3 = np.zeros((len(fe_only_f), len(universe.trajectory)))
# for i in range(len(fe_only_b)):
#     for j in range(len(universe.trajectory)):
#         sorted1 = np.sort(bond_lengths_time[j, int(fe_only_b[i])])[0:6]
#         sorted2 = np.sort(bond_lengths_time[j, int(fe_only_d[i])])[0:6]
#         sorted3 = np.sort(bond_lengths_time[j, int(fe_only_f[i])])[0:6]
#         temp1[i, j] = func_metric(sorted1, bond_lengths_mean_1[j], bond_lengths_mean_2[j])
#         temp2[i, j] = func_metric(sorted2, bond_lengths_mean_1[j], bond_lengths_mean_2[j])
#         temp3[i, j] = func_metric(sorted3, bond_lengths_mean_1[j], bond_lengths_mean_2[j])
#     ax_6.plot(time_plot, temp1[i, :], 'r')
#     ax_6.plot(time_plot, temp2[i, :], 'g')
#     ax_6.plot(time_plot, temp3[i, :], 'b')
# ax_6.plot(time_plot[0], temp1[0, 0], 'r-', label='Fe B')
# ax_6.plot(time_plot[0], temp2[0, 0], 'g-', label='Fe D')
# ax_6.plot(time_plot[0], temp3[0, 0], 'b-', label='Fe F')
# ax_6.set_xlabel('Time / fs')
# ax_6.set_ylabel('Average Fe-O bond length / A')
# ax_6.set_xlim([0, len(universe.trajectory)*timestep])
# ax_6.set_ylim(ylim_2)
# if draw_legend: ax_6.legend(frameon=False)
# fig_6.tight_layout()
# fig_6.savefig('{}/metric_color_layer_{}.png'.format(folder_save, run), dpi=300, bbbox_inches='tight')
#
# # Plot  metric (color coded by layer)
# time_plot = np.linspace(start=0, stop=len(universe.trajectory)*timestep, num=len(universe.trajectory))
# metric = np.zeros((len(universe.trajectory)))
# fig_6, ax_6 = plt.subplots()
# temp3 = np.zeros((len(fe_only_f), len(universe.trajectory)))
# for i in range(len(fe_only_f)):
#     for j in range(len(universe.trajectory)):
#         sorted3 = np.sort(bond_lengths_time[j, int(fe_only_f[i])])[0:6]
#         temp3[i, j] = func_metric(sorted3, bond_lengths_mean_1[j], bond_lengths_mean_2[j])
#     ax_6.plot(steps1_4 - steps1_4[0], temp3[i, :], '-', color=plot_color[i], label='Fe {}'.format(i + 1))
# ax_6.set_xlabel('Time / fs')
# ax_6.set_ylabel('Average Fe-O bond length / A')
# ax_6.set_xlim([0, len(universe.trajectory)*timestep])
# ax_6.set_ylim(ylim_2)
# # ax_6.set_xlim([700, 1500])
# ax_6.set_xlim([1000, 1500])
# if draw_legend: ax_6.legend(frameon=False)
# fig_6.tight_layout()
# fig_6.savefig('{}/metric_color_atom_{}.png'.format(folder_save, run), dpi=300)

# Plot all Fe-O (color coded by layer)
# time_plot = np.linspace(start=0, stop=len(universe.trajectory)*timestep, num=len(universe.trajectory))
# metric = np.zeros((len(universe.trajectory)))
# fig_6, ax_6 = plt.subplots()
# temp3 = np.zeros((len(fe_only_f), len(universe.trajectory), 6))
# for i in range(len(fe_only_f)):
#     for j in range(len(universe.trajectory)):
#         sorted3 = np.sort(bond_lengths_time[j, int(fe_only_f[i])])[0:6]
#         temp3[i, j, :] = sorted3
#     ax_6.plot(steps1_4 - steps1_4[0], temp3[i, :, 0], '-', color=plot_color[i], label='Fe {}'.format(i + 1))
#     ax_6.plot(steps1_4 - steps1_4[0], temp3[i, :, 1], '-', color=plot_color[i], label='Fe {}'.format(i + 1))
#     ax_6.plot(steps1_4 - steps1_4[0], temp3[i, :, 2], '-', color=plot_color[i], label='Fe {}'.format(i + 1))
#     ax_6.plot(steps1_4 - steps1_4[0], temp3[i, :, 3], '-', color=plot_color[i], label='Fe {}'.format(i + 1))
#     ax_6.plot(steps1_4 - steps1_4[0], temp3[i, :, 4], '-', color=plot_color[i], label='Fe {}'.format(i + 1))
#     ax_6.plot(steps1_4 - steps1_4[0], temp3[i, :, 5], '-', color=plot_color[i], label='Fe {}'.format(i + 1))
# ax_6.set_xlabel('Time / fs')
# ax_6.set_ylabel('Average Fe-O bond length / A')
# ax_6.set_xlim([0, len(universe.trajectory)*timestep])
# ax_6.set_ylim(ylim_2)
# ax_6.set_xlim([700, 1500])
# if draw_legend: ax_6.legend(frameon=False)
# fig_6.tight_layout()
# fig_6.savefig('{}/bonds_color_atom_{}.png'.format(folder_save, run), dpi=300)

if __name__ == "__main__":
    print('Finished.')
    plt.show()