import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import mean_absolute_error
from sklearn.metrics import root_mean_squared_error

# 1. Load data
# model = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/shared/denan_li/multi-task'
# model = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/shared/denan_li/testing/multi-task-copy'
# database = '{}/database_spin'.format(model)

# model = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/shared/denan_li/delete/denan_li/multi-task'
# database = '{}/database_spin'.format(model)

# model = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/shared/denan_li/delete/denan_li/multi-task-copy'
# database = '{}/database_spin'.format(model)

# model = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/chain/hf/frozen-fe/multi-task'
# database = '{}/database_spin'.format(model)

# HF with frozen Fe, H2O flies away
# model = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/chain/hf/frozen-fe/multi-task'
# database = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/chain/hf/frozen-fe/multi-task/database_spin'
# model = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/chain/hf/frozen-fe/multi-task2'
# database = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/chain/hf/frozen-fe/database'

# HF with frozen Fe and H2O
# model = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/chain/hf/frozen-fe-h2o/multi-task'
# database = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/chain/hf/frozen-fe-h2o/database'

# HSE06(35%) with frozen Fe and H2O, single hop
# model = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/chain/hse/frozen-fe-h2o/hops-1/multi-task'
# model = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/chain/hse/frozen-fe-h2o/hops-1-database-spin/multi-task'
# database = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/chain/hse/frozen-fe-h2o/database'
# axis_lim_y = np.array([3.7, 4.5])

# Bulk hematite, no hops with apram
# model = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/delete/hops-0/multi-task'
# model = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/delete/hops-0/multi-task-2-se_e3'
# model = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/delete/hops-0/multi-task-2-dpa1-se_atten_v2'
# model = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/delete/hops-0/multi-task-2-dpa2'
# database = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/delete/hops-0/database'

# Bulk hematite, no hops with species
# model = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/delete/hops-0-species/multi-task'  # force failed to train
# model = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/delete/hops-0-species/multi-task-2'  # force failed to train
# model = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/delete/hops-0-species/multi-task-2-se_e3'  # force failed to train
# # model = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/delete/hops-0-species/multi-task-2-dpa1'  # all failed to train
# model = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/delete/hops-0-species/multi-task-2-dpa1-se_atten_v2'  # force failed to train, but best of them
# model = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/delete/hops-0-species/multi-task-2-dpa2'  # force failed to train, but still perfect
# database = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/delete/hops-0-species/database'

# Bulk hematite, no hops with apram and database_ener_force, database_spin
# model = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/delete/hops-0-database-spin/multi-task-se_e2_a'
# model = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/delete/hops-0-database-spin/multi-task-se_e3'  # force failed to train, probably because axis_neuron=12 for all others while this uses default axis_neuron=4
# model = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/delete/hops-0-database-spin/multi-task-dpa1-se_atten_v2'
# model = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/delete/hops-0-database-spin/multi-task-dpa2'
# database = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/delete/hops-0-database-spin/database_spin'

# Bulk hematite, hops-1 with apram and database_ener_force, database_spin
# model = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/delete/hops-1-database-spin/multi-task-se_e2_a'
# model = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/delete/hops-1-database-spin/multi-task-dpa1-se_atten_v2'
# model = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/delete/hops-1-database-spin/multi-task-dpa2'
# database = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/delete/hops-1-database-spin/database_spin'

# Bulk hematite 400k f hops-1
# model_ener = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/400k-f/single-fit-ener-se_e2_a'
# model_spin = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/400k-f/single-fit-m-se_e2_a'
# model_ener = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/400k-f/single-fit-ener-dpa2'
# model_spin = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/400k-f/single-fit-m-dpa2'
# model_ener = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/400k-f/multi-task-se_e2_a'
# model_ener = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/400k-f/multi-task-se_e3'
# model_ener = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/400k-f/multi-task-dpa1-se_atten_v2'
# model_ener = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/400k-f/multi-task-dpa2'
# model_spin = model_ener
# database = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/400k-f/database_spin'

# Bulk hematite 400k b hops-5
# model_ener = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/400k-b/single-fit-ener-se_e2_a'
# model_spin = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/400k-b/single-fit-m-se_e2_a'
# model_ener = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/400k-b/multi-task-se_e2_a'
# model_ener = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/400k-b/multi-task-dpa1-se_atten_v2'
# model_ener = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/400k-b/multi-task-dpa2'
# model_spin = model_ener
# database = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/400k-b/database_spin'

# Bulk hematite 400k b hops-5 test set (400k-f)
# model_ener = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/400k-b/test-set/400k-f/multi-task-se_e2_a'
# model_ener = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/400k-b/test-set/400k-f/multi-task-dpa2'
# model_spin = model_ener
# database = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/400k-b/test-set/400k-f/database_spin'

# Bulk hematite 0K NEB
# model = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/geo-opt-all/multi-task-se_e2_a'
# model = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/geo-opt-all/multi-task-dpa2'
# database = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/geo-opt-all/database_spin'
# axis_lim_y = np.array([-4.05, -3.05])
# axis_lim_y = np.array([-4, -3])


# Bulk hematite 0K NEB + DFT-MD single-fit
model_ener = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/400k-b-geo-opt-all/single-fit-ener-se_e2_a'
model_spin = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/400k-b-geo-opt-all/single-fit-m-se_e2_a'
database = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/400k-b-geo-opt-all/database_spin'
axis_lim_y = np.array([-4.05, -3.05])

# Bulk hematite 0K NEB + DFT-MD multi-fit
# model_ener = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/400k-b-geo-opt-all/multi-task-se_e3'
# model_ener = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/400k-b-geo-opt-all/multi-task-dpa1-se_atten_v2'
# model_ener = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/400k-b-geo-opt-all/multi-task-dpa2'
# model_spin = model_ener
# database = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/400k-b-geo-opt-all/database_spin'
# axis_lim_y = np.array([-4.05, -3.05])

axis_lim_y = np.array([-4, -3])
transition_time = np.array([460, 117.5, 315, 656, 721, 1283])
zoom = False
transition_time_plot = 5
axis_lim_x_zoom = np.array([transition_time[transition_time_plot]-40, transition_time[transition_time_plot]+40])

dft_e = np.load("{}/set.000/energy.npy".format(database))
dft_f = np.load("{}/set.000/force.npy".format(database))
dft_s = np.load("{}/set.000/atom_ener.npy".format(database))
spin_0 = np.load("{}/0_spin.npy".format(model_spin))
ener_0 = np.load("{}/0_ener.npy".format(model_ener))
force_0 = np.load("{}/0_force.npy".format(model_ener))
spin_1 = np.load("{}/1_spin.npy".format(model_spin))
ener_1 = np.load("{}/1_ener.npy".format(model_ener))
force_1 = np.load("{}/1_force.npy".format(model_ener))


def plot_ener(dft, dp, ax, color_plot, pos, text, title=None):
    ax.plot(dft.flatten(), dp.flatten(), '.', color=color_plot)
    ax.set_xlabel("DFT energy (eV)")
    ax.set_ylabel("DP energy (eV)")
    mae = mean_absolute_error(dft.flatten(), dp.flatten()) / 120 * 1000  # unit: meV/atom
    rmse = root_mean_squared_error(dft.flatten(), dp.flatten()) / 120 * 1000  # unit: meV/atom
    print('plot_ener')
    print('mean_absolute_error', mae)
    print('root_mean_squared_error', rmse)
    min_val, max_val = min(dft.min(), dp.min()), max(dft.max(), dp.max())
    ax.plot([min_val, max_val], [min_val, max_val], 'k--')
    # ax.text(pos[0], pos[1], f"{text} MAE: {mae:.2f} meV/atom", transform=ax.transAxes)
    ax.text(pos[0], pos[1], f"{text} RMSE: {rmse:.2f} meV/atom", transform=ax.transAxes)
    if title is not None:
        ax.set_title(title)
    else:
        ax.set_title("Energy")


def plot_force(dft, dp, ax, color_plot, pos, text, title=None):
    ax.plot(dft.flatten(), dp.flatten(), '.', color=color_plot)
    ax.set_xlabel("DFT force (eV/$\AA$)")
    ax.set_ylabel("DP force (eV/$\AA$)")
    mae = mean_absolute_error(dft.flatten(), dp.flatten()) * 1000  # unit: meV/A
    rmse = root_mean_squared_error(dft.flatten(), dp.flatten()) * 1000  # unit: meV/A
    print('plot_force')
    print('mean_absolute_error', mae)
    print('root_mean_squared_error', rmse)
    # ax.text(pos[0], pos[1], f"{text} MAE: {mae:.2f} meV/$\AA$", transform=ax.transAxes)
    ax.text(pos[0], pos[1], f"{text} RMSE: {rmse:.2f} meV/$\AA$", transform=ax.transAxes)
    min_val, max_val = min(dft.min(), dp.min()), max(dft.max(), dp.max())
    ax.plot([min_val, max_val], [min_val, max_val], 'k--')
    if title is not None:
        ax.set_title(title)
    else:
        ax.set_title("Force")


def plot_spin(dft, dp, ax, color_plot, pos, text, title=None):
    ax.plot(dft.flatten(), dp.flatten(), '.', color=color_plot)
    ax.set_xlabel("DFT spin")
    ax.set_ylabel("DP spin")
    min_val, max_val = min(dft.min(), dp.min()), max(dft.max(), dp.max())
    mae = mean_absolute_error(dft.flatten(), dp.flatten())
    rmse = root_mean_squared_error(dft.flatten(), dp.flatten())
    print('plot_spin')
    print('mean_absolute_error', mae)
    print('root_mean_squared_error', rmse)
    # ax.text(pos[0], pos[1], f"{text} MAE: {mae:.2f}", transform=ax.transAxes)
    ax.text(pos[0], pos[1], f"{text} RMSE: {rmse:.2f}", transform=ax.transAxes)
    ax.plot([min_val, max_val], [min_val, max_val], 'k--')

    ax.set_xlim([axis_lim_y[0], axis_lim_y[1]])
    ax.set_ylim([axis_lim_y[0], axis_lim_y[1]])

    if title is not None:
        ax.set_title(title)
    else:
        ax.set_title("Spin")


def plot_spin_time1(dft, dp, ax, axis_lim_y, title=None):
    num_atoms = int(dp.shape[1])
    num_timesteps = int(dp.shape[0])
    time_array = np.linspace(0, int(num_timesteps / 2), num=num_timesteps)
    # num_atoms_plot_spin = 7
    num_atoms_plot_spin = 120
    plotting_colors = ['r', 'g', 'b', 'm', 'grey', 'orange', 'brown', 'hotpink'] * 100
    dft = np.reshape(dft, (num_timesteps, num_atoms))

    fe_b = np.array([27, 45, 18, 14, 25, 29, 42, 16]) - 1
    fe_d = np.array([2, 6, 17, 13, 4, 38, 41, 15]) - 1
    fe_f = np.array([28, 46, 1, 5, 26, 30, 3, 37]) - 1
    fe_plot = np.arange(0, 120)
    print(fe_plot)

    for i in range(fe_plot.shape[0]):
        # ax.plot(time_array, dft[:, fe_plot[i]], 'x-', color=plotting_colors[i], label='Fe {}'.format(i+1))
        # ax.plot(time_array, dp[:, fe_plot[i]], 'x--', color=plotting_colors[i])
        ax.plot(time_array, dft[:, fe_plot[i]], '-', color=plotting_colors[i], label='Fe {}'.format(i+1))
        ax.plot(time_array, dp[:, fe_plot[i]], '--', color=plotting_colors[i])

    # for i in range(num_atoms_plot_spin):
    #     ax.plot(time_array, dft[:, i], '-', color=plotting_colors[i], label='Fe {}'.format(i))
    #     ax.plot(time_array, dp[:, i], '--', color=plotting_colors[i])

    ax.set_xlim(0, time_array.shape[0]/2)
    ax.set_ylim([axis_lim_y[0], axis_lim_y[1]])
    ax.set_xlabel("Time / fs")
    ax.set_ylabel("Spin moment")
    # ax.legend()


# 2. Plot 2x3
pos_array = np.array(([0.5, 0.1], [0.5, 0.2]))
text_array = ['Test', 'Train']
color_plot_array = ['r', 'g']
i = 1
fig, axes = plt.subplots(2, 3, figsize=(15, 6))
plot_ener(dft_e, ener_0, axes[0, 0], color_plot=color_plot_array[i], pos=pos_array[i], text=text_array[i], title="Energy, No-aparam")
plot_ener(dft_e, ener_1, axes[1, 0], color_plot=color_plot_array[i], pos=pos_array[i], text=text_array[i], title="Energy, Yes-aparam")
plot_force(dft_f, force_0, axes[0, 1], color_plot=color_plot_array[i], pos=pos_array[i], text=text_array[i], title="Force, No-aparam")
plot_force(dft_f, force_1, axes[1, 1], color_plot=color_plot_array[i], pos=pos_array[i], text=text_array[i], title="Force, Yes-aparam")
plot_spin(dft_s, spin_0, axes[0, 2], color_plot=color_plot_array[i], pos=pos_array[i], text=text_array[i], title="Spin, No-aparam")
plot_spin(dft_s, spin_1, axes[1, 2], color_plot=color_plot_array[i], pos=pos_array[i], text=text_array[i], title="Spin, Yes-aparam")
plt.tight_layout()
plt.savefig("{}/fit_2x3_folders_{}.png".format(model_ener, len(model_ener)), dpi=600)

# 2. Plot 1x3
text_array = ['Test', 'Train']
color_plot_array = ['r', 'g']
i = 1
fig2, axes2 = plt.subplots(1, 3, figsize=(15, 3))
plot_ener(dft_e, ener_1, axes2[0], color_plot=color_plot_array[i], pos=pos_array[i], text=text_array[i], title="Energy, Yes-aparam")
plot_force(dft_f, force_1, axes2[1], color_plot=color_plot_array[i], pos=pos_array[i], text=text_array[i], title="Force, Yes-aparam")
plot_spin(dft_s, spin_1, axes2[2], color_plot=color_plot_array[i], pos=pos_array[i], text=text_array[i], title="Spin, Yes-aparam")
plt.tight_layout()
plt.savefig("{}/fit_1x3_folders_{}.png".format(model_ener, len(model_ener)), dpi=600)

# 2. Plot
# fig2, axes2 = plt.subplots(figsize=(15, 6))
fig2, axes2 = plt.subplots()
plot_spin_time1(dft_s, spin_1, axes2, axis_lim_y, title="Energy, No-aparam")
plt.tight_layout()
plt.savefig("{}/spin_time.png".format(model_ener), dpi=600)
if zoom:
    plt.xlim(axis_lim_x_zoom[0], axis_lim_x_zoom[1])
plt.savefig("{}/spin_time_zoom_{}.png".format(model_ener, transition_time_plot), dpi=600)

if __name__ == "__main__":
    print('Finished.')
    plt.show()
