import dpdata
import numpy as np
import matplotlib.pyplot as plt
from general import parameters as param

# training_systems = dpdata.LabeledSystem(
#     "/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/shared/denan_li/chris/database", fmt="deepmd/npy"
# )
# predict = training_systems.predict("/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/shared/denan_li/chris/single-fit-ener/graph.pb")
#
#
# plt.scatter(training_systems["energies"], predict["energies"])
#
# x_range = np.linspace(plt.xlim()[0], plt.xlim()[1])
#
# plt.plot(x_range, x_range, "r--", linewidth=0.25)
# plt.xlabel("Energy of DFT")
# plt.ylabel("Energy predicted by deep potential")
# plt.plot()

import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import mean_absolute_error

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
model = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/chain/hse/frozen-fe-h2o/hops-1-database-spin/multi-task'
database = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/chain/hse/frozen-fe-h2o/database'
axis_lim = np.array([3.7, 4.5])

# Bulk hematite, no hops
# model = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/hops-0/multi-task'
# model = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/hops-0/multi-task-2-se_e3'
# model = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/hops-0/multi-task-2-dpa1-se_atten_v2'
# model = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/hops-0/multi-task-2-dpa2'
# database = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/hops-0/database'

# Bulk hematite, no hops
# model = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/hops-0-species/multi-task'  # force failed to train
# model = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/hops-0-species/multi-task-2'  # force failed to train
# model = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/hops-0-species/multi-task-2-se_e3'  # force failed to train
# # model = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/hops-0-species/multi-task-2-dpa1'  # all failed to train
# model = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/hops-0-species/multi-task-2-dpa1-se_atten_v2'  # force failed to train, but best of them
# model = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/hops-0-species/multi-task-2-dpa2'  # force failed to train, but still perfect
# database = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/hops-0-species/database'
# axis_lim = np.array([-4.05, -3.05])

dft_e = np.load("{}/set.000/energy.npy".format(database))
dft_f = np.load("{}/set.000/force.npy".format(database))
dft_s = np.load("{}/set.000/atom_ener.npy".format(database))

spin_0 = np.load("{}/0_spin.npy".format(model))
ener_0 = np.load("{}/0_ener.npy".format(model))
force_0 = np.load("{}/0_force.npy".format(model))

spin_1 = np.load("{}/1_spin.npy".format(model))
ener_1 = np.load("{}/1_ener.npy".format(model))
force_1 = np.load("{}/1_force.npy".format(model))


def plot_ener(dft, dp, ax, title=None):
    ax.plot(dft.flatten(), dp.flatten(), 'o')
    ax.set_xlabel("DFT energy (eV)")
    ax.set_ylabel("DP energy (eV)")
    mae = mean_absolute_error(dft.flatten(), dp.flatten()) / 77 * 1000 # unit: meV/atom
    ax.text(0.5, 0.1, f"MAE: {mae:.2f} meV/atom", transform=ax.transAxes)
    min_val, max_val = min(dft.min(), dp.min()), max(dft.max(), dp.max())
    ax.plot([min_val, max_val], [min_val, max_val], 'k--')
    if title is not None:
        ax.set_title(title)
    else:
        ax.set_title("Energy")

def plot_force(dft, dp, ax, title=None):
    ax.plot(dft.flatten(), dp.flatten(), 'o')
    ax.set_xlabel("DFT force (eV/$\AA$)")
    ax.set_ylabel("DP force (eV/$\AA$)")
    mae = mean_absolute_error(dft.flatten(), dp.flatten())
    ax.text(0.5, 0.1, f"MAE: {mae:.2f} eV/$\AA$", transform=ax.transAxes)
    min_val, max_val = min(dft.min(), dp.min()), max(dft.max(), dp.max())
    ax.plot([min_val, max_val], [min_val, max_val], 'k--')
    if title is not None:
        ax.set_title(title)
    else:
        ax.set_title("Force")


def plot_spin_time1(dft, dp, ax, axis_lim, title=None):
    num_atoms = int(dp.shape[1])
    num_timesteps = int(dp.shape[0])
    time_array = np.linspace(0, int(num_timesteps / 2), num=num_timesteps)
    num_atoms_plot_spin = 7
    # num_atoms_plot_spin = 120
    plotting_colors = ['r', 'g', 'b', 'm', 'grey', 'orange', 'cyan'] * 100
    dft = np.reshape(dft, (num_timesteps, num_atoms))

    for i in range(num_atoms_plot_spin):
        ax.plot(time_array, dft[:, i], '-', color=plotting_colors[i], label='Fe {}'.format(i))
        ax.plot(time_array, dp[:, i], '--', color=plotting_colors[i])

    ax.set_ylim([axis_lim[0], axis_lim[1]])
    ax.set_xlabel("Time / fs")
    ax.set_ylabel("Spin")
    ax.legend()


def plot_spin(dft, dp, ax, axis_lim, title=None):
    ax.plot(dft.flatten(), dp.flatten(), 'o')
    ax.set_xlabel("DFT spin")
    ax.set_ylabel("DP spin")
    min_val, max_val = min(dft.min(), dp.min()), max(dft.max(), dp.max())
    mae = mean_absolute_error(dft.flatten(), dp.flatten())
    ax.text(0.5, 0.1, f"MAE: {mae:.2f}", transform=ax.transAxes)
    ax.plot([min_val, max_val], [min_val, max_val], 'k--')

    ax.set_xlim([axis_lim[0], axis_lim[1]])
    ax.set_ylim([axis_lim[0], axis_lim[1]])

    if title is not None:
        ax.set_title(title)
    else:
        ax.set_title("Spin")

# 2. Plot
fig, axes = plt.subplots(2, 3, figsize=(15, 6))
plot_ener(dft_e, ener_0, axes[0, 0], title="Energy, No-aparam")
plot_ener(dft_e, ener_1, axes[1, 0], title="Energy, Yes-aparam")
plot_force(dft_f, force_0, axes[0, 1], title="Force, No-aparam")
plot_force(dft_f, force_1, axes[1, 1], title="Force, Yes-aparam")
plot_spin(dft_s, spin_0, axes[0, 2], axis_lim, title="Spin, No-aparam", )
plot_spin(dft_s, spin_1, axes[1, 2], axis_lim, title="Spin, Yes-aparam")
plt.tight_layout()
plt.savefig("{}/fit.png".format(model), dpi=600)

# 2. Plot
fig2, axes2 = plt.subplots(figsize=(15, 6))
plot_spin_time1(dft_s, spin_1, axes2, axis_lim, title="Energy, No-aparam")
plt.tight_layout()
plt.savefig("{}/spin_time.png".format(model), dpi=600)

if __name__ == "__main__":
    print('Finished.')
    plt.show()
