import numpy as np
import matplotlib.pyplot as plt
from general import parameters as param

folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/chain/hse/frozen-fe-h2o/hops-1/ase_md'
folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/chain/hse/frozen-fe-h2o/hops-1-database-spin/dpmd/denan'
folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/chain/hematitese/frozen-fe-h2o/hops-1-database-spin/dpmd/denan-redo'
# folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/chain/hse/frozen-fe-h2o/hops-1-database-spin/dpmd/denan-script'
# folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/chain/hse/frozen-fe-h2o/hops-1-database-spin/dpmd/denan-script-xyz'

spin = np.load("{}/spin_history.npy".format(folder))
data_ase = np.loadtxt('{}/md.log'.format(folder), skiprows=1)
plotting_colors = ['r', 'g', 'b', 'm', 'grey', 'orange', 'brown', 'hotpink'] * 100

num_atoms_plot_spin = 7
num_timesteps = 10000
time_array = np.linspace(0, int(num_timesteps / 2), num=num_timesteps)
data = spin[:num_timesteps, :num_atoms_plot_spin]

# fig2, axes2 = plt.subplots(figsize=(15, 6))
fig2, axes2 = plt.subplots()
for i in range(num_atoms_plot_spin):
    axes2.plot(time_array, data[:, i], '-', color=plotting_colors[i], label='Fe {}'.format(i+1))
axes2.legend()
axes2.set_xlabel("Time / fs")
axes2.set_ylabel("Spin")
fig2.tight_layout()
fig2.savefig("{}/dp_md_spin.png".format(folder))

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
