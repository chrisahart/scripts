import numpy as np
import matplotlib.pyplot as plt
from general import parameters as param

folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/chain/hse/frozen-fe-h2o/hops-1/ase_md'

spin = np.load("{}/spin_history.npy".format(folder))
plotting_colors = ['r', 'g', 'b', 'm', 'grey', 'orange', 'cyan'] * 100

num_atoms_plot_spin = 7
num_timesteps = 20000
time_array = np.linspace(0, int(num_timesteps / 2), num=num_timesteps)
data = spin[:num_timesteps, :num_atoms_plot_spin]

fig2, axes2 = plt.subplots(figsize=(15, 6))
for i in range(num_atoms_plot_spin):
    axes2.plot(time_array, data[:, i], '-', color=plotting_colors[i], label='Fe {}'.format(i))
axes2.legend()
axes2.set_xlabel("Time / fs")
axes2.set_ylabel("Spin")
fig2.tight_layout()
fig2.savefig("{}/dp_md_spin.png".format(folder))

if __name__ == "__main__":
    print('Finished.')
    plt.show()
