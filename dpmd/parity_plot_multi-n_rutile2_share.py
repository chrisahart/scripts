import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import mean_absolute_error
from sklearn.metrics import root_mean_squared_error


def plot_ener(dft, dp, ax, color_plot, pos, text, num_atoms, title=None):
    ax.plot(dft.flatten(), dp.flatten(), '.', color=color_plot)
    ax.set_xlabel("DFT energy (eV)")
    ax.set_ylabel("DP energy (eV)")
    mae = mean_absolute_error(dft.flatten(), dp.flatten()) / num_atoms * 1000  # unit: meV/atom
    rmse = root_mean_squared_error(dft.flatten(), dp.flatten()) / num_atoms * 1000  # unit: meV/atom
    print('plot_ener')
    print('mean_absolute_error', mae)
    print('root_mean_squared_error', rmse)
    min_val, max_val = min(dft.min(), dp.min()), max(dft.max(), dp.max())
    ax.plot([min_val, max_val], [min_val, max_val], 'k--')
    ax.text(pos[0], pos[1], f"{text} RMSE: {rmse:.3f} meV/atom", transform=ax.transAxes, color=color_plot)
    if title is not None:
        ax.set_title(title)


def plot_force(dft, dp, ax, color_plot, pos, text, polaron_index=None, title=None):

    if polaron_index is None:
        ax.plot(dft.flatten(), dp.flatten(), '.', color=color_plot)
        ax.set_xlabel("DFT force (eV/Å)")
        ax.set_ylabel("DP force (eV/Å)")
        mae = mean_absolute_error(dft.flatten(), dp.flatten()) * 1000  # unit: meV/A
        rmse = root_mean_squared_error(dft.flatten(), dp.flatten()) * 1000  # unit: meV/A
        print('plot_force')
        print('mean_absolute_error', mae)
        print('root_mean_squared_error', rmse)
        ax.text(pos[0], pos[1], f"RMSE: {rmse:.3f} meV/Å", transform=ax.transAxes, color=color_plot)
        min_val, max_val = min(dft.min(), dp.min()), max(dft.max(), dp.max())
        ax.plot([min_val, max_val], [min_val, max_val], 'k--')
        if title is not None:
            ax.set_title(title)
    else:
        ax.plot([[-10, -10], [10, 10]], [[-10, -10], [10, 10]], 'k--')
        polaron_index = np.repeat(polaron_index, 3, axis=1)
        dft = dft.flatten()
        dp = dp.flatten()
        polaron_index = polaron_index.flatten()

        dft2 = []
        dp2 = []
        for i in range(len(polaron_index)):
            if polaron_index[i] == 1:
                dft2.append(dft[i])
                dp2.append(dp[i])
        ax.plot(dft2, dp2, '.', color=color_plot)
        ax.set_xlabel("DFT force (eV/Å)")
        ax.set_ylabel("DP force (eV/Å)")
        mae = mean_absolute_error(dft2, dp2) * 1000  # unit: meV/A
        rmse = root_mean_squared_error(dft2, dp2) * 1000  # unit: meV/A
        print('plot_force')
        print('mean_absolute_error', mae)
        print('root_mean_squared_error', rmse)
        # ax.text(pos[0], pos[1], f"{text} MAE: {mae:.2f} meV/$\AA$", transform=ax.transAxes)
        # ax.text(pos[0], pos[1], f"{text} RMSE: {rmse:.2f} meV/$\AA$", transform=ax.transAxes, color=color_plot)
        ax.text(pos[0], pos[1], f"RMSE: {rmse:.3f} meV/Å", transform=ax.transAxes, color=color_plot)
        min_val, max_val = min(dft.min(), dp.min()), max(dft.max(), dp.max())
        ax.plot([min_val, max_val], [min_val, max_val], 'k--')
        if title is not None:
            ax.set_title(title)


def plot_spin(dft, dp, ax, color_plot, pos, text, polaron_index=None, title=None):

    ax.plot([[-10, -10], [10, 10]], [[-10, -10], [10, 10]], 'k--')

    if polaron_index is None:
        ax.plot(dft.flatten(), dp.flatten(), '.', color=color_plot)
        ax.set_xlabel("DFT spin")
        ax.set_ylabel("DP spin")
        min_val, max_val = min(dft.min(), dp.min()), max(dft.max(), dp.max())
        mae = mean_absolute_error(dft.flatten(), dp.flatten())
        rmse = root_mean_squared_error(dft.flatten(), dp.flatten())
        print('plot_spin')
        print('mean_absolute_error', mae)
        print('root_mean_squared_error', rmse)
        ax.text(pos[0], pos[1], f"RMSE: {rmse:.3f}", transform=ax.transAxes, color=color_plot)
        ax.set_xlim([axis_lim_y[0], axis_lim_y[1]])
        ax.set_ylim([axis_lim_y[0], axis_lim_y[1]])

    else:
        ax.plot([[-10, -10], [10, 10]], [[-10, -10], [10, 10]], 'k--')
        dft = dft.flatten()
        dp = dp.flatten()
        polaron_index = polaron_index.flatten()
        dft2 = []
        dp2 = []
        for i in range(len(polaron_index)):
            if polaron_index[i] == 1:
                dft2.append(dft[i])
                dp2.append(dp[i])
        ax.plot(dft2, dp2, '.', color=color_plot)
        ax.set_xlabel("DFT spin")
        ax.set_ylabel("DP spin")
        mae = mean_absolute_error(dft2, dp2)
        rmse = root_mean_squared_error(dft2, dp2)
        print('plot_spin polaron_index is not None')
        print('mean_absolute_error', mae)
        print('root_mean_squared_error', rmse)
        ax.text(pos[0], pos[1], f"RMSE: {rmse:.3f}", transform=ax.transAxes, color=color_plot)
        ax.set_xlim([axis_lim_y[0], axis_lim_y[1]])
        ax.set_ylim([axis_lim_y[0], axis_lim_y[1]])

    if title is not None:
        ax.set_title(title)


def plot_spin_time1(dft, dp, ax, axis_lim_y, num_atoms, title=None):
    num_atoms = int(dp.shape[1])
    num_timesteps = int(dp.shape[0])
    time_array = np.linspace(0, int(num_timesteps), num=num_timesteps)
    num_atoms_plot_spin = num_atoms
    plotting_colors = ['r', 'g', 'b', 'm', 'grey', 'orange', 'brown', 'hotpink'] * 100

    for i in range(0, num_atoms_plot_spin):
        ax.plot(time_array, dft[:, i], '-', color=plotting_colors[i], label='Fe {}'.format(i))
        ax.plot(time_array, dp[:, i], '--', color=plotting_colors[i])

    ax.set_xlim(0, time_array.shape[0])
    ax.set_ylim([axis_lim_y[0], axis_lim_y[1]])
    ax.set_xlabel("Time / fs")
    ax.set_ylabel("Spin moment")


def plot_spin_time1_total(dft, dp, ax, axis_lim_y, num_atoms, title=None):
    num_timesteps = int(dp.shape[0])
    time_array = np.linspace(0, int(num_timesteps), num=num_timesteps)
    plotting_colors = ['r', 'g', 'b', 'm', 'grey', 'orange', 'brown', 'hotpink'] * 100

    dft_total_spin = np.sum(dft, axis=1)
    dp_total_spin = np.sum(dp, axis=1)

    mae = mean_absolute_error(dft_total_spin.flatten(), dp_total_spin.flatten())
    rmse = root_mean_squared_error(dft_total_spin.flatten(), dp_total_spin.flatten())
    print('plot_spin_time1_total')
    print('plot_spin')
    print('mean_absolute_error', mae)
    print('root_mean_squared_error', rmse)

    ax.plot(time_array, dft_total_spin, '-', color=plotting_colors[0], label='DFT')
    ax.plot(time_array, dp_total_spin, '-', color=plotting_colors[1], label='NNP')
    ax.set_xlim(0, time_array.shape[0])
    ax.set_ylim([axis_lim_y[0], axis_lim_y[1]])
    ax.set_xlabel("Time / fs")
    ax.legend()


# Bulk TiO2 leopold
model = ['ener-dpa3-start_pref-0.02-1000_limit_pref-1-1', 'pop-dpa3-pref_pop-1000-1']
folder = '/Volumes/Elements/Data/Postdoc2/Data/Work/calculations/tio2/rutile/deepmd/rutile/leopold-fixed3'
model_ener = ['{}/{}'.format(folder, model[0])] * 3
model_spin = ['{}/{}'.format(folder, model[1])] * 3
database = ['{}/database_train/'.format(folder),
            '{}/database_test/1'.format(folder),
            '{}/database_test/2'.format(folder)]
val = ['_0_', '_1_', '_2_']
axis_lim_y = np.array([-0.02, 0.95])
ylim_spin_time = [0, 1.6]
ylim_population_time = [239, 246]
pos_array_energy = np.array(([0.45, 0.25], [0.45, 0.15], [0.45, 0.05]))
pos_array_force = np.array(([0.6, 0.25], [0.6, 0.15], [0.6, 0.05]))
pos_array_spin = np.array(([0.72, 0.25], [0.72, 0.15], [0.72, 0.05]))
text_array = ['Train', 'Valid', 'Test']
num_atoms = 288
color_plot_array = ['r', 'g', 'b', 'm']

dft_e = []
dft_f = []
dft_s = []
dft_polaron = []
spin_1 = []
ener_1 = []
force_1 = []
for i in range(len(database)):
    dft_e.append(np.load("{}/set.000/energy.npy".format(database[i], val[i])))
    dft_f.append(np.load("{}/set.000/force.npy".format(database[i], val[i])))
    dft_polaron.append(np.load("{}/set.000/aparam.npy".format(database[i], val[i])))
    ener_1.append(np.load("{}/1{}ener.npy".format(model_ener[i], val[i])))
    force_1.append(np.load("{}/1{}force.npy".format(model_ener[i], val[i])))
    spin_1.append(np.load("{}/1{}spin.npy".format(model_spin[i], val[i])))
    dft_s.append(np.load("{}/set.000/atomic_population.npy".format(database[i], val[i])))

# Plot parity 1x3 subplot
print('Plot parity 1x3 subplot')
fig2, axes2 = plt.subplots(1, 3, figsize=(15, 5))
for i in range(len(database)):
    print(i)
    plot_ener(dft_e[i], ener_1[i], axes2[0], num_atoms=num_atoms, color_plot=color_plot_array[i],
              pos=pos_array_energy[i], text=text_array[i])
    plot_force(dft_f[i], force_1[i], axes2[1], color_plot=color_plot_array[i], pos=pos_array_force[i],
               text=text_array[i])
    plot_spin((dft_s[i][:, :, 0] - dft_s[i][:, :, 1]), (spin_1[i][:, :, 0] - spin_1[i][:, :, 1]), axes2[2],
                color_plot=color_plot_array[i], pos=pos_array_spin[i], text=text_array[i])
plt.tight_layout()
for i in range(len(database)):
    plt.savefig("{}/fit_1x3_folders_{}.png".format(model_spin[i], len(model_ener)), dpi=600)
    plt.savefig("{}/fit_1x3_folders_{}.png".format(model_ener[i], len(model_spin)), dpi=600)

# Plot energy parity
print('Plot energy parity')
fig_energy, axes_energy = plt.subplots(figsize=(5, 5))
for i in range(len(database)):
    print(i)
    plot_ener(dft_e[i], ener_1[i], axes_energy, num_atoms=num_atoms, color_plot=color_plot_array[i],
              pos=pos_array_energy[i], text=text_array[i])
plt.tight_layout()
for i in range(len(database)):
    plt.savefig("{}/energy_{}.png".format(model_spin[i], len(model_ener)), dpi=600)
    plt.savefig("{}/energy_{}.png".format(model_ener[i], len(model_spin)), dpi=600)

# Plot force parity
print('Plot force parity')
fig_force, axes_force = plt.subplots(figsize=(5, 5))
for i in range(len(database)):
    print(i)
    plot_force(dft_f[i], force_1[i], axes_force, color_plot=color_plot_array[i], pos=pos_array_force[i],
               text=text_array[i])
plt.tight_layout()
for i in range(len(database)):
    plt.savefig("{}/force_{}_{}.png".format(model_spin[i], len(model_ener), val[0]), dpi=600)
    plt.savefig("{}/force_{}_{}.png".format(model_ener[i], len(model_spin), val[0]), dpi=600)

# Plot force parity polaron
print('Plot force parity polaron')
fig_force2, axes_force2 = plt.subplots(figsize=(5, 5))
for i in range(len(database)):
    print(i)
    plot_force(dft_f[i], force_1[i], axes_force2, color_plot=color_plot_array[i], pos=pos_array_force[i],
               text=text_array[i], polaron_index=dft_polaron[i])
plt.tight_layout()
for i in range(len(database)):
    plt.savefig("{}/force_polaron_{}_{}.png".format(model_spin[i], len(model_ener), val[0]), dpi=600)
    plt.savefig("{}/force_polaron_{}_{}.png".format(model_ener[i], len(model_spin), val[0]), dpi=600)

# Plot spin polaron
print('Plot spin polaron')
fig3, axes3 = plt.subplots(figsize=(5, 5))
for i in range(len(database)):
    print(i)
    plot_spin((dft_s[i][:, :, 0] - dft_s[i][:, :, 1]), (spin_1[i][:, :, 0] - spin_1[i][:, :, 1]), axes3,
              color_plot=color_plot_array[i], pos=pos_array_spin[i], text=text_array[i], polaron_index=dft_polaron[i])
plt.tight_layout()
for i in range(len(database)):
    plt.savefig("{}/spin_moment_polaron_{}.png".format(model_spin[i], len(model_ener)), dpi=600)

# Plot spin
print('Plot spin')
fig4, axes4 = plt.subplots(figsize=(5, 5))
for i in range(len(database)):
    print(i)
    plot_spin((dft_s[i][:, :, 0] - dft_s[i][:, :, 1]), (spin_1[i][:, :, 0] - spin_1[i][:, :, 1]), axes4,
              color_plot=color_plot_array[i], pos=pos_array_spin[i], text=text_array[i])
plt.tight_layout()
for i in range(len(database)):
    plt.savefig("{}/spin_moment_all_{}.png".format(model_spin[i], len(model_ener)), dpi=600)

if __name__ == "__main__":
    print('Finished.')
    plt.show()
