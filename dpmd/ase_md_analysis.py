import numpy as np
import matplotlib.pyplot as plt
from general import parameters as param
from ase.io import read

folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/mgo/deepmd/md/cell-222/electron-u-6/temp-600-velocity-dft/single-task-se_e2_a'
folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/mgo/deepmd/md/cell-222/electron-u-6/temp-600-velocity-dft/multi-task-se_e2_a'
folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/mgo/deepmd/md/cell-222/electron-u-6/temp-600-velocity-dft/single-task-dpa3'
# folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/mgo/deepmd/md/cell-222/electron-u-6/temp-600-velocity-dft/multi-task-dpa3'

spin = np.load("{}/spin_history.npy".format(folder))
traj = read("{}/md.traj".format(folder))
num_atoms = 64

num_timesteps = spin.shape[0]
time_array = np.linspace(0, int(num_timesteps / 2), num=num_timesteps)

draw_legend = True
plotting_colors = ['r', 'g', 'b', 'm', 'grey', 'orange', 'brown', 'hotpink'] * 100
xlim_1 = [0, (5103+200)/2]
xlim_2 = [0, time_array[-1]]

# Calculate polaron atom
polaron_atom_time = np.zeros(num_timesteps, dtype=int)
for j in range(num_timesteps):
    polaron_atom_time[j] = int(np.argmax(spin[j, :]))
polaron_atoms = np.unique(polaron_atom_time)
print('polaron_atoms', polaron_atoms)

# Plot all spin
fig_spin1, ax_spin1 = plt.subplots(figsize=(15, 6))
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
ax_spin1.set_xlim(xlim_2)
fig_spin1.tight_layout()
fig_spin1.savefig("{}/dp_md_spin.png".format(folder))
ax_spin1.set_xlim(xlim_1)
fig_spin1.savefig("{}/dp_md_spin_lim.png".format(folder))
fig_spin1.tight_layout()

# Plot total spin
fig_spin2, ax_spin2 = plt.subplots(figsize=(15, 6))
# fig_spin2, ax_spin2 = plt.subplots()
ax_spin2.plot(time_array, np.sum(spin, axis=1), 'k-')
ax_spin2.set_xlim(0, time_array[-1])
ax_spin2.set_xlabel("Time / fs")
ax_spin2.set_ylabel("Spin")
ax_spin2.set_xlim(xlim_2)
fig_spin2.tight_layout()
fig_spin2.savefig("{}/dp_md_spin_sum.png".format(folder))
ax_spin2.set_xlim(xlim_1)
fig_spin2.savefig("{}/dp_md_spin_sum_lim.png".format(folder))
fig_spin2.tight_layout()

# time = data_ase[:,0] * 1000 # ps to fs
# pe = data_ase[:,2]
# temp = data_ase[:,-1]
# fig, ax = plt.subplots(figsize=(4, 3))
# ax.plot(time, pe, color='blue')
# ax.set_xlabel('Time (fs)')
# ax.set_ylabel('Energy (eV)')
# for tl in ax.get_yticklabels():
#     tl.set_color('blue')
# ax2 = ax.twinx()
# ax2.plot(time, temp, color='red')
# ax2.set_ylabel('Temperature (K)')
# for tl in ax2.get_yticklabels():
#     tl.set_color('red')
# ax2.set_ylim(150, 400)
# plt.tight_layout()
# plt.savefig('thermo.png', dpi=600)

if __name__ == "__main__":
    print('Finished.')
    plt.show()
