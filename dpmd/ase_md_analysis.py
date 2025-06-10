import numpy as np
import matplotlib.pyplot as plt
from general import parameters as param
from ase.io import read

# folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/mgo/deepmd/md/cell-222/electron-u-6/temp-600-velocity-dft/single-task-se_e2_a'
# folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/mgo/deepmd/md/cell-222/electron-u-6/temp-600-velocity-dft/multi-task-se_e2_a'
# folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/mgo/deepmd/md/cell-222/electron-u-6/temp-600-velocity-dft/single-task-dpa3'
# folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/mgo/deepmd/md/cell-222/electron-u-6/temp-600-velocity-dft/multi-task-dpa3'

# folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/mgo/deepmd/md/cell-222/electron-u-6/temp-400-velocity-dft/single-task-dpa3'
# folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/mgo/deepmd/md/cell-222/electron-u-6/temp-500-velocity-dft/single-task-dpa3'
# folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/mgo/deepmd/md/cell-222/electron-u-6/temp-600-velocity-dft/single-task-dpa3'
folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/mgo/deepmd/md/cell-222/electron-u-6/temp-700-velocity-dft/single-task-dpa3'
folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/mgo/deepmd/md/cell-222/electron-u-6/temp-800-velocity-dft/single-task-dpa3'

spin = np.load("{}/spin_history.npy".format(folder))
traj = read("{}/md.traj".format(folder))
num_atoms = 64

num_timesteps = spin.shape[0]
time_array = np.linspace(0, int(num_timesteps / 2), num=num_timesteps)

draw_legend = True
plotting_colors = ['r', 'g', 'b', 'm', 'grey', 'orange', 'brown', 'hotpink'] * 100
# plotting_colors = ['b', 'g', 'r', 'm', 'grey', 'orange', 'brown', 'hotpink'] * 100
xlim_1 = [0, (5103+200)/2]
xlim_1 = [0, time_array[-1]]
xlim_2 = [0, time_array[-1]]
# ylim_1 = [0, 0.9]
ylim_1 = [-0.02, 0.9]


# Calculate distance between current timestep polaron atom and next timestep
# Then get all non-zero answers


# Plot average of 6 O-Mg bonds
metric = np.zeros((num_atoms_mg, num_timesteps2))
fig_bonds_1, ax_bonds_1 = plt.subplots()
for i in range(num_atoms_o):
    ax_bonds_1.plot(time_val_1 - time_val_1[0], bond_lengths_time_sorted_mean[:, i], '-')
for j in range(polaron_atoms.shape[0]):
    ax_bonds_1.plot(time_val_1 - time_val_1[0], bond_lengths_time_sorted_mean[:, polaron_atoms[j]-num_atoms_mg], '-', color=plotting_colors[j], label='{}'.format(polaron_atoms[j]+1))
ax_bonds_1.set_xlabel('Time / fs')
# ax_bonds_1.set_xlabel('Timestep')
ax_bonds_1.set_ylabel('Average of 6 O-Mg bond lengths / A')
if draw_legend: ax_bonds_1.legend(frameon=True)
# # ax_bonds_1.set_xlim([0, len(universe.trajectory)])
ax_bonds_1.set_xlim(xlim_1)
# ax_bonds_1.set_xlim([0, len(universe.trajectory) * timestep])
# ax_bonds_1.set_ylim([0.06, -0.10])
fig_bonds_1.savefig('{}/bond_lengths_average.png'.format(folder_save), dpi=300)
fig_bonds_1.tight_layout()





# Calculate polaron atom
polaron_atom_time = np.zeros(num_timesteps, dtype=int)
for j in range(num_timesteps):
    polaron_atom_time[j] = int(np.argmax(spin[j, :]))
polaron_atoms = np.unique(polaron_atom_time)
# print('polaron_atoms', polaron_atoms)

# Plot total spin
fig_spin2, ax_spin2 = plt.subplots(figsize=(10, 4))
# fig_spin2, ax_spin2 = plt.subplots()
ax_spin2.plot(time_array, np.sum(spin, axis=1), 'k-')
ax_spin2.set_xlim(0, time_array[-1])
ax_spin2.set_xlabel("Time / fs")
ax_spin2.set_ylabel("Spin")
ax_spin2.set_xlim(xlim_2)
fig_spin2.tight_layout()
fig_spin2.savefig("{}/dp_md_spin_sum.png".format(folder), dpi=600)
ax_spin2.set_xlim(xlim_1)
fig_spin2.savefig("{}/dp_md_spin_sum_lim.png".format(folder), dpi=600)
fig_spin2.tight_layout()

# Plot all spin
fig_spin1, ax_spin1 = plt.subplots(figsize=(10, 4))
# fig_spin1, ax_spin1 = plt.subplots()
for i in range(num_atoms):
    ax_spin1.plot(time_array, spin[:, i], '-', color=plotting_colors[i])
for j in range(polaron_atoms.shape[0]):
    print('(polaron_atoms[j]', polaron_atoms[j])
    ax_spin1.plot(time_array, spin[:, polaron_atoms[j]], '-', color=plotting_colors[j], label='{}'.format(polaron_atoms[j]+1))
# ax_spin1.plot(time_array, spin[:, 52], '-', color=plotting_colors[0], label='{}'.format(52))
# ax_spin1.plot(time_array, np.sum(spin, axis=1), 'k-', alpha=0.2)
if draw_legend: ax_spin1.legend(frameon=True)
ax_spin1.set_xlim(0, time_array[-1])
ax_spin1.set_xlabel("Time / fs")
ax_spin1.set_ylabel("Spin")
ax_spin1.set_ylim(ylim_1)
ax_spin1.set_xlim(xlim_2)
fig_spin1.tight_layout()
fig_spin1.savefig("{}/dp_md_spin.png".format(folder), dpi=600)
ax_spin1.set_xlim(xlim_1)
fig_spin1.savefig("{}/dp_md_spin_lim.png".format(folder), dpi=600)
fig_spin1.tight_layout()

if __name__ == "__main__":
    print('Finished.')
    plt.show()
