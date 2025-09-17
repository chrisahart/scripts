import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import mean_absolute_error
from sklearn.metrics import root_mean_squared_error
from general import parameters as param


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
    # ax.text(pos[0], pos[1], f"{text} MAE: {mae:.2f} meV/atom", transform=ax.transAxes)
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
        # ax.text(pos[0], pos[1], f"{text} MAE: {mae:.2f} meV/$\AA$", transform=ax.transAxes)
        # ax.text(pos[0], pos[1], f"{text} RMSE: {rmse:.2f} meV/$\AA$", transform=ax.transAxes, color=color_plot)
        ax.text(pos[0], pos[1], f"RMSE: {rmse:.3f} meV/Å", transform=ax.transAxes, color=color_plot)
        min_val, max_val = min(dft.min(), dp.min()), max(dft.max(), dp.max())
        ax.plot([min_val, max_val], [min_val, max_val], 'k--')
        if title is not None:
            ax.set_title(title)
    else:
        ax.plot(dft.flatten(), dp.flatten(), '.', color=color_plot)


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
        # ax.text(pos[0], pos[1], f"{text} MAE: {mae:.2f}", transform=ax.transAxes)
        # ax.text(pos[0], pos[1], f"{text} RMSE: {rmse:.2f}", transform=ax.transAxes, color=color_plot)
        ax.text(pos[0], pos[1], f"RMSE: {rmse:.3f}", transform=ax.transAxes, color=color_plot)
        ax.set_xlim([axis_lim_y[0], axis_lim_y[1]])
        ax.set_ylim([axis_lim_y[0], axis_lim_y[1]])

    else:
        ax.plot([[-10, -10], [10, 10]], [[-10, -10], [10, 10]], 'k--')
        dft = dft.flatten()
        dp = dp.flatten()
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
    # num_atoms_plot_spin = 64
    num_atoms_plot_spin = int(num_atoms/3 * 2)
    plotting_colors = ['r', 'g', 'b', 'm', 'grey', 'orange', 'brown', 'hotpink'] * 100

    print('np.shape(dft)', np.shape(dft))
    # dft = np.reshape(dft, (num_timesteps, num_atoms))

    # fe_b = np.array([27, 45, 18, 14, 25, 29, 42, 16]) - 1
    # fe_d = np.array([2, 6, 17, 13, 4, 38, 41, 15]) - 1
    # fe_f = np.array([28, 46, 1, 5, 26, 30, 3, 37]) - 1
    # fe_plot = np.arange(0, 120)
    # print(fe_plot)
    # for i in range(fe_plot.shape[0]):
        # ax.plot(time_array, dft[:, fe_plot[i]], 'x-', color=plotting_colors[i], label='Fe {}'.format(i+1))
        # ax.plot(time_array, dp[:, fe_plot[i]], 'x--', color=plotting_colors[i])
        # ax.plot(time_array, dft[:, fe_plot[i]], '-', color=plotting_colors[i], label='Fe {}'.format(i+1))
        # ax.plot(time_array, dp[:, fe_plot[i]], '--', color=plotting_colors[i])

    # for i in range(num_atoms_plot_spin):
    for i in range(0, num_atoms_plot_spin):
        ax.plot(time_array, dft[:, i], '-', color=plotting_colors[i], label='Fe {}'.format(i))
        ax.plot(time_array, dp[:, i], '--', color=plotting_colors[i])

    ax.set_xlim(0, time_array.shape[0]/2)
    ax.set_ylim([axis_lim_y[0], axis_lim_y[1]])
    ax.set_xlabel("Time / fs")
    ax.set_ylabel("Spin moment")
    # ax.legend()


plot_aparam = False

# Bulk TiO2 leopold DELETE
# model = ['single-fit-ener-se_e2_a', 'single-fit-m-se_e2_a']
# model = ['single-fit-ener-dpa2', 'single-fit-m-dpa2']
# model = ['single-fit-ener-dpa2-finetune', 'single-fit-m-dpa2-finetune']
# model = ['single-fit-ener-dpa3-6-default', 'single-fit-m-dpa3-6-default']
# model = ['single-fit-ener-dpa3-16-official', 'single-fit-m-dpa3-16-official']
# model = ['single-fit-ener-dpa3-16-official-finetune', 'single-fit-m-dpa3-16-official-finetune']
# spin_is_population = True
# folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/deepmd/rutile/leopold'
# model_ener = ['{}/{}'.format(folder, model[0])] * 4
# model_spin = ['{}/{}'.format(folder, model[1])] * 4
# database = ['{}/database_train/'.format(folder),
#             '{}/database_test/1'.format(folder),
#             '{}/database_test/2'.format(folder)]
# val = ['_0_', '_1_', '_2_']
# axis_lim_y = np.array([0, 1])
# pos_array_energy = np.array(([0.45, 0.25], [0.45, 0.15], [0.45, 0.05]))
# pos_array_force = np.array(([0.55, 0.25], [0.55, 0.15], [0.55, 0.05]))
# pos_array_spin = np.array(([0.72, 0.25], [0.72, 0.15], [0.72, 0.05]))
# text_array = ['Train', 'Valid', 'Test']
# num_atoms = 288

# Bulk TiO2 leopold
# model = ['single-fit-ener-se_e2_a', 'single-fit-m-se_e2_a']
# # model = ['single-fit-ener-dpa2', 'single-fit-m-dpa2']
# # model = ['single-fit-ener-dpa2-official', 'single-fit-m-dpa2-official']
# # model = ['single-fit-ener-dpa2-finetune', 'single-fit-m-dpa2-finetune']
# # model = ['single-fit-ener-dpa3-6-default', 'single-fit-m-dpa3-6-default']
# # model = ['single-fit-ener-dpa3-16-official', 'single-fit-m-dpa3-16-official']
# # model = ['single-fit-ener-dpa3-16-official-finetune', 'single-fit-m-dpa3-16-official-finetune']
# spin_is_population = True
# folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/deepmd/rutile/leopold'
# database_folder = folder
# model_ener = ['{}/{}'.format(folder, model[0])] * 4
# model_spin = ['{}/{}'.format(folder, model[1])] * 4
# database = ['{}/database_train/'.format(database_folder),
#             '{}/database_test/1'.format(database_folder),
#             '{}/database_test/2'.format(database_folder)]
# val = ['_0_', '_1_', '_2_']
# axis_lim_y = np.array([0, 1])
# pos_array_energy = np.array(([0.45, 0.25], [0.45, 0.15], [0.45, 0.05]))
# pos_array_force = np.array(([0.6, 0.25], [0.6, 0.15], [0.6, 0.05]))
# pos_array_spin = np.array(([0.72, 0.25], [0.72, 0.15], [0.72, 0.05]))
# text_array = ['Train', 'Valid', 'Test']
# num_atoms = 288

# Bulk TiO2 PBE+U
database_folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/deepmd/rutile/336/md-cell-opt/database/pbe-u-4-ps'
# model = ['single-fit-ener-se_e2_a-v3.1.0a0-polaron', 'single-fit-pop-se_e2_a-v3.1.0a0-polaron']
# model = ['single-fit-ener-dpa3-nlayers-6-new-v3.1.0a0-polaron', 'single-fit-pop-dpa3-nlayers-6-new-v3.1.0a0-polaron']
# model = ['single-fit-ener-dpa3-nlayers-6-new2-v3.1.0a0-polaron', 'single-fit-pop-dpa3-nlayers-6-new2-v3.1.0a0-polaron']
# spin_is_population = True
# database = ['{}/database_population_train/'.format(database_folder),
#             '{}/database_population_test/1'.format(database_folder),
#             '{}/database_population_test/2'.format(database_folder)]
# model = ['single-fit-ener-se_e2_a-official-v3.1.0', 'single-fit-spin-se_e2_a-official-v3.1.0']
# model = ['single-fit-ener-se_e2_a-official-v3.1.0-sel-90', 'single-fit-spin-se_e2_a-official-v3.1.0-sel-90']
# model = ['single-fit-ener-dpa3-nlayers-6-new2-v3.1.0', 'single-fit-spin-dpa3-nlayers-6-new2-v3.1.0']
# model = ['single-fit-ener-dpa3-nlayers-6-new2-v3.1.0-4gpu', 'single-fit-spin-dpa3-nlayers-6-new2-v3.1.0-4gpu']
# model = ['single-fit-ener-dpa3-16-official-all-v3.1.0-4gpu', 'single-fit-spin-dpa3-16-official-all-v3.1.0-4gpu']
# model = ['single-fit-ener-dpa3-nlayers-6-new2-v3.1.0a0-polaron', 'single-fit-spin-dpa3-nlayers-6-new2-v3.1.0a0-polaron']
# model = ['single-fit-ener-dpa3-nlayers-6-new2-v3.1.0a0-polaron', 'single-fit-m-dpa3-freeze-test']
# model = ['single-fit-ener-dpa3-nlayers-6-new2-v3.1.0a0-polaron', 'single-fit-m-dpa3-freeze-test-2']
# model = ['single-fit-ener-dpa3-nlayers-6-new2-v3.1.0a0-polaron', 'single-fit-m-dpa3-freeze-test-3']
# model = ['single-fit-ener-dpa3-nlayers-6-new2-v3.1.0a0-polaron', 'single-fit-m-dpa3-freeze-test-4']
model = ['single-fit-ener-dpa3-nlayers-6-new2-v3.1.0a0-polaron', 'single-fit-m-dpa3']
spin_is_population = False
database = ['{}/database_spin_train/'.format(database_folder),
            '{}/database_spin_test/1'.format(database_folder),
            '{}/database_spin_test/2'.format(database_folder)]
folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/deepmd/rutile/336/md-cell-opt/pbe-u-4-ps-v2'
model_ener = ['{}/{}'.format(folder, model[0])] * 4
model_spin = ['{}/{}'.format(folder, model[1])] * 4
val = ['_0_', '_1_', '_2_']
axis_lim_y = np.array([0, 1])
pos_array_energy = np.array(([0.45, 0.25], [0.45, 0.15], [0.45, 0.05]))
pos_array_force = np.array(([0.6, 0.25], [0.6, 0.15], [0.6, 0.05]))
pos_array_spin = np.array(([0.72, 0.25], [0.72, 0.15], [0.72, 0.05]))
text_array = ['Train', 'Valid', 'Test']
num_atoms = 324

# Bulk TiO2 PBE+U no aparam
# database_folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/deepmd/rutile/336/md-cell-opt/database/pbe-u-4-ps'
# model = ['single-fit-ener-se_e2_a-official-v3.1.0-sel-90', 'single-fit-spin-se_e2_a-official-v3.1.0-sel-90']
# model = ['single-fit-ener-se_e2_a-official-v3.1.0a0-polaron-sel-90', 'single-fit-spin-se_e2_a-official-v3.1.0a0-polaron-sel-90']
# model = ['single-fit-ener-dpa3-nlayers-6-new2-v3.1.0', 'single-fit-spin-dpa3-nlayers-6-new2-v3.1.0']
# model = ['single-fit-ener-dpa3-nlayers-6-new2-v3.1.0a0-polaron', 'single-fit-spin-dpa3-nlayers-6-new2-v3.1.0a0-polaron']
# spin_is_population = False
# database = ['{}/database_spin_train/'.format(database_folder),
#             '{}/database_spin_test/1'.format(database_folder),
#             '{}/database_spin_test/2'.format(database_folder)]
# folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/deepmd/rutile/336/md-cell-opt/pbe-u-4-ps-v2-tip'
# model_ener = ['{}/{}'.format(folder, model[0])] * 4
# model_spin = ['{}/{}'.format(folder, model[1])] * 4
# val = ['_0_', '_1_', '_2_']
# axis_lim_y = np.array([0, 1])
# pos_array_energy = np.array(([0.45, 0.25], [0.45, 0.15], [0.45, 0.05]))
# pos_array_force = np.array(([0.6, 0.25], [0.6, 0.15], [0.6, 0.05]))
# pos_array_spin = np.array(([0.72, 0.25], [0.72, 0.15], [0.72, 0.05]))
# text_array = ['Train', 'Valid', 'Test']
# num_atoms = 324

# # Bulk TiO2 336 22% hse-22
# # model = ['single-fit-ener-se_e2_a', 'single-fit-m-se_e2_a']
# # model = ['single-fit-ener-dpa3', 'single-fit-m-dpa3']
# # model = ['multi-task-se_e2_a', 'multi-task-se_e2_a']
# model = ['multi-task-dpa3', 'multi-task-dpa3']
# spin_is_population = True
# folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/deepmd/rutile/hse-22'
# model_ener = ['{}/{}'.format(folder, model[0])] * 4
# model_spin = ['{}/{}'.format(folder, model[1])] * 4
# database = ['{}/database_population_train/'.format(folder),
#             '{}/database_population_test/1'.format(folder),
#             '{}/database_population_test/2'.format(folder)]
# val = ['_0_', '_1_', '_2_']
# axis_lim_y = np.array([0, 1])
# pos_array_energy = np.array(([0.45, 0.25], [0.45, 0.15], [0.45, 0.05]))
# pos_array_force = np.array(([0.6, 0.25], [0.6, 0.15], [0.6, 0.05]))
# pos_array_spin = np.array(([0.72, 0.25], [0.72, 0.15], [0.72, 0.05]))
# text_array = ['Train', 'Valid', 'Test']
# num_atoms = 324

# Bulk TiO2 336 25% hse-25-3-ps (no hops)
# model = ['single-fit-ener-se_e2_a', 'single-fit-m-se_e2_a']
# model = ['multi-task-se_e2_a', 'multi-task-se_e2_a']
# model = ['single-fit-ener-dpa3', 'single-fit-m-dpa3']
# model = ['multi-task-dpa3', 'multi-task-dpa3']
# spin_is_population = True
# folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/deepmd/rutile/hse-25-3-ps'
# model_ener = ['{}/{}'.format(folder, model[0])] * 4
# model_spin = ['{}/{}'.format(folder, model[1])] * 4
# database = ['{}/database_population_train/'.format(folder),
#             '{}/database_population_test/1'.format(folder),
#             '{}/database_population_test/2'.format(folder)]
# val = ['_0_', '_1_', '_2_']
# axis_lim_y = np.array([0, 1])
# pos_array_energy = np.array(([0.45, 0.25], [0.45, 0.15], [0.45, 0.05]))
# pos_array_force = np.array(([0.6, 0.25], [0.6, 0.15], [0.6, 0.05]))
# pos_array_spin = np.array(([0.72, 0.25], [0.72, 0.15], [0.72, 0.05]))
# text_array = ['Train', 'Valid', 'Test']
# num_atoms = 324

# Bulk TiO2 336 25% hse-25-4-ps (many hops)
# model = ['single-fit-ener-se_e2_a', 'single-fit-m-se_e2_a']
# model = ['multi-task-se_e2_a', 'multi-task-se_e2_a']
# model = ['single-fit-ener-dpa3', 'single-fit-m-dpa3']
# model = ['multi-task-dpa3', 'multi-task-dpa3']
# spin_is_population = True
# folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/deepmd/rutile/hse-25-4-ps'
# model_ener = ['{}/{}'.format(folder, model[0])] * 4
# model_spin = ['{}/{}'.format(folder, model[1])] * 4
# database = ['{}/database_population_train/'.format(folder),
#             '{}/database_population_test/1'.format(folder),
#             '{}/database_population_test/2'.format(folder)]
# val = ['_0_', '_1_', '_2_']
# axis_lim_y = np.array([0, 1])
# pos_array_energy = np.array(([0.45, 0.25], [0.45, 0.15], [0.45, 0.05]))
# pos_array_force = np.array(([0.6, 0.25], [0.6, 0.15], [0.6, 0.05]))
# pos_array_spin = np.array(([0.72, 0.25], [0.72, 0.15], [0.72, 0.05]))
# text_array = ['Train', 'Valid', 'Test']
# num_atoms = 324

# # Bulk TiO2 336 22% hse-22-v2
# model = ['single-fit-ener-se_e2_a', 'single-fit-m-se_e2_a']
# model = ['single-fit-ener-dpa3-nlayers-6-old', 'single-fit-m-dpa3-nlayers-6-old']
# # model = ['single-fit-ener-dpa3-nlayers-6-new', 'single-fit-m-dpa3-nlayers-6-new']
# # model = ['single-fit-ener-dpa3-nlayers-6-new2', 'single-fit-m-dpa3-nlayers-6-new2']
# # model = ['single-fit-ener-dpa3-nlayers-6-new2-dynamic', 'single-fit-m-dpa3-nlayers-6-new2-dynamic']
# # model = ['single-fit-ener-dpa3-16-official', 'single-fit-m-dpa3-16-official']
# # model = ['single-fit-ener-dpa3-16-official-all', 'single-fit-m-dpa3-16-official-all']
# spin_is_population = True
# folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/deepmd/rutile/hse-22-v2'
# model_ener = ['{}/{}'.format(folder, model[0])] * 4
# model_spin = ['{}/{}'.format(folder, model[1])] * 4
# database = ['{}/database_population_train/'.format(folder),
#             '{}/database_population_test/1'.format(folder),
#             '{}/database_population_test/2'.format(folder)]
# val = ['_0_', '_1_', '_2_']
# axis_lim_y = np.array([0, 1])
# pos_array_energy = np.array(([0.45, 0.25], [0.45, 0.15], [0.45, 0.05]))
# pos_array_force = np.array(([0.6, 0.25], [0.6, 0.15], [0.6, 0.05]))
# pos_array_spin = np.array(([0.72, 0.25], [0.72, 0.15], [0.72, 0.05]))
# text_array = ['Train', 'Valid', 'Test']
# num_atoms = 324

# # Bulk TiO2 336 22% multi-fit hse-22-25-4-ps-multi
# model = ['multi-task-ener-dpa3', 'multi-task-m-dpa3']
# spin_is_population = True
# folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/deepmd/rutile/hse-22-25-4-ps-multi'
# model_ener = ['{}/{}'.format(folder, model[0])] * 4
# model_spin = ['{}/{}'.format(folder, model[1])] * 4
# database_folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/deepmd/rutile/database/hse-22-backup-4-ps'
# database = ['{}/database_population_train/'.format(database_folder),
#             '{}/database_population_test/1'.format(database_folder),
#             '{}/database_population_test/2'.format(database_folder)]
# val = ['_0_', '_1_', '_2_']
# axis_lim_y = np.array([0, 1])
# pos_array_energy = np.array(([0.45, 0.25], [0.45, 0.15], [0.45, 0.05]))
# pos_array_force = np.array(([0.6, 0.25], [0.6, 0.15], [0.6, 0.05]))
# pos_array_spin = np.array(([0.72, 0.25], [0.72, 0.15], [0.72, 0.05]))
# text_array = ['Train', 'Valid', 'Test']
# num_atoms = 324

print('model', model)
color_plot_array = ['r', 'g', 'b', 'm']
# axis_lim_y = np.array([-4, -3])
transition_time = np.array([460, 117.5, 315, 656, 721, 1283])
zoom = False
transition_time_plot = 5
axis_lim_x_zoom = np.array([transition_time[transition_time_plot]-40, transition_time[transition_time_plot]+40])

dft_e = []
dft_f = []
dft_s = []
dft_polaron = []
spin_0 = []
ener_0 = []
force_0 = []
spin_1 = []
ener_1 = []
force_1 = []

# test = np.load("{}/database_population_train/set.000/atomic_spin.npy".format(folder))
# test = np.load("{}/database_population_test/1/set.000/atomic_spin.npy".format(folder))
# test = np.load("{}/database_train/set.000/atomic_spin.npy".format(folder))
# print('test', test[0, :, :])
# test = np.load("{}/database_train/set.000/force.npy".format(folder))
# print('test', test.shape)
# print('test', test[0, :])

for i in range(len(database)):
    dft_e.append(np.load("{}/set.000/energy.npy".format(database[i], val[i])))
    dft_f.append(np.load("{}/set.000/force.npy".format(database[i], val[i])))
    dft_polaron.append(np.load("{}/set.000/aparam.npy".format(database[i], val[i])))
    ener_0.append(np.load("{}/0{}ener.npy".format(model_ener[i], val[i])))
    force_0.append(np.load("{}/0{}force.npy".format(model_ener[i], val[i])))
    ener_1.append(np.load("{}/1{}ener.npy".format(model_ener[i], val[i])))
    force_1.append(np.load("{}/1{}force.npy".format(model_ener[i], val[i])))
    spin_0.append(np.load("{}/0{}spin.npy".format(model_spin[i], val[i])))
    spin_1.append(np.load("{}/1{}spin.npy".format(model_spin[i], val[i])))
    if not spin_is_population:
        dft_s.append(np.load("{}/set.000/atom_ener.npy".format(database[i], val[i])))
    if spin_is_population:
        dft_s.append(np.load("{}/set.000/atomic_spin.npy".format(database[i], val[i])))

# Plot parity 2x3 subplots
if plot_aparam:
    fig, axes = plt.subplots(2, 3, figsize=(15, 6))
    for i in range(len(database)):
        print(i)
        plot_ener(dft_e[i], ener_0[i], axes[0, 0], num_atoms=num_atoms, color_plot=color_plot_array[i], pos=pos_array_energy[i], text=text_array[i], title="Energy, No-aparam")
        plot_ener(dft_e[i], ener_1[i], axes[1, 0], num_atoms=num_atoms, color_plot=color_plot_array[i], pos=pos_array_energy[i], text=text_array[i], title="Energy, Yes-aparam")
        plot_force(dft_f[i], force_0[i], axes[0, 1], color_plot=color_plot_array[i], pos=pos_array_force[i], text=text_array[i], title="Force, No-aparam")
        plot_force(dft_f[i], force_1[i], axes[1, 1], color_plot=color_plot_array[i], pos=pos_array_force[i], text=text_array[i], title="Force, Yes-aparam")
        if not spin_is_population:
            plot_spin(dft_s[i], spin_0[i], axes[0, 2], color_plot=color_plot_array[i], pos=pos_array_spin[i], text=text_array[i], title="Spin, No-aparam")
            plot_spin(dft_s[i], spin_1[i], axes[1, 2], color_plot=color_plot_array[i], pos=pos_array_spin[i], text=text_array[i], title="Spin, Yes-aparam")
        if spin_is_population:
            plot_spin((dft_s[i][:, :, 0] - dft_s[i][:, :, 1]), (spin_0[i][:, :, 0] - spin_0[i][:, :, 1]), axes[0, 2], color_plot=color_plot_array[i], pos=pos_array_spin[i], text=text_array[i], title="Spin, No-aparam")
            plot_spin((dft_s[i][:, :, 0] - dft_s[i][:, :, 1]), (spin_1[i][:, :, 0] - spin_1[i][:, :, 1]), axes[1, 2], color_plot=color_plot_array[i], pos=pos_array_spin[i], text=text_array[i], title="Spin, Yes-aparam")
    plt.tight_layout()
    for i in range(len(database)):
        plt.savefig("{}/fit_2x3_folders_{}.png".format(model_spin[i], len(model_ener)), dpi=600)
        plt.savefig("{}/fit_2x3_folders_{}.png".format(model_ener[i], len(model_spin)), dpi=600)

# Plot parity 1x3 subplot
fig2, axes2 = plt.subplots(1, 3, figsize=(15, 3))
for i in range(len(database)):
    print(i)
    plot_ener(dft_e[i], ener_1[i], axes2[0], num_atoms=num_atoms, color_plot=color_plot_array[i], pos=pos_array_energy[i], text=text_array[i])
    plot_force(dft_f[i], force_1[i], axes2[1], color_plot=color_plot_array[i], pos=pos_array_force[i], text=text_array[i])
    if not spin_is_population:
        plot_spin(dft_s[i], spin_1[i], axes2[2], color_plot=color_plot_array[i], pos=pos_array_spin[i], text=text_array[i])
    if spin_is_population:
        plot_spin((dft_s[i][:, :, 0] - dft_s[i][:, :, 1]), (spin_1[i][:, :, 0] - spin_1[i][:, :, 1]), axes2[2], color_plot=color_plot_array[i], pos=pos_array_spin[i], text=text_array[i])
plt.tight_layout()
for i in range(len(database)):
    plt.savefig("{}/fit_1x3_folders_{}.png".format(model_spin[i], len(model_ener)), dpi=600)
    plt.savefig("{}/fit_1x3_folders_{}.png".format(model_ener[i], len(model_spin)), dpi=600)


# Plot energy parity
fig_energy, axes_energy = plt.subplots(figsize=(5, 5))
for i in range(len(database)):
    print(i)
    plot_ener(dft_e[i], ener_1[i], axes_energy, num_atoms=num_atoms, color_plot=color_plot_array[i], pos=pos_array_energy[i], text=text_array[i])
plt.tight_layout()
for i in range(len(database)):
    plt.savefig("{}/energy_{}.png".format(model_spin[i], len(model_ener)), dpi=600)
    plt.savefig("{}/energy_{}.png".format(model_ener[i], len(model_spin)), dpi=600)

# Plot force parity
fig_force, axes_force = plt.subplots(figsize=(5, 5))
for i in range(len(database)):
    print(i)
    plot_force(dft_f[i], force_1[i], axes_force, color_plot=color_plot_array[i], pos=pos_array_force[i], text=text_array[i])
plt.tight_layout()
for i in range(len(database)):
    plt.savefig("{}/force_{}.png".format(model_spin[i], len(model_ener)), dpi=600)
    plt.savefig("{}/force_{}.png".format(model_ener[i], len(model_spin)), dpi=600)

# Plot spin
if spin_is_population:
    fig3, axes3 = plt.subplots(figsize=(5, 5))
    for i in range(len(database)):
        print(i)
        plot_spin((dft_s[i][:, :, 0] - dft_s[i][:, :, 1]), (spin_1[i][:, :, 0] - spin_1[i][:, :, 1]), axes3, color_plot=color_plot_array[i], pos=pos_array_spin[i], text=text_array[i], polaron_index=dft_polaron[i])
    plt.tight_layout()
    for i in range(len(database)):
        plt.savefig("{}/spin_moment_polaron_{}.png".format(model_spin[i], len(model_ener)), dpi=600)

# Plot spin polaron
if spin_is_population:
    fig4, axes4 = plt.subplots(figsize=(5, 5))
    for i in range(len(database)):
        print(i)
        plot_spin((dft_s[i][:, :, 0] - dft_s[i][:, :, 1]), (spin_1[i][:, :, 0] - spin_1[i][:, :, 1]), axes4, color_plot=color_plot_array[i], pos=pos_array_spin[i], text=text_array[i])
    plt.tight_layout()
    for i in range(len(database)):
        plt.savefig("{}/spin_moment_all_{}.png".format(model_spin[i], len(model_ener)), dpi=600)

# Plot spin moment training
if spin_is_population:
    axis_lim_y = [-0.02, 1.0]
    plot_spin_time = True
    # plot_spin_time = False
    if plot_spin_time:
        fig_spin_train, axes_spin_train = plt.subplots()
        plot_spin_time1((dft_s[0][:, :, 0] - dft_s[0][:, :, 1]), (spin_1[0][:, :, 0] - spin_1[0][:, :, 1]),
                        axes_spin_train, axis_lim_y, num_atoms=num_atoms, title="Energy, No-aparam")
        plt.tight_layout()
        plt.savefig("{}/spin_train.png".format(model_spin[0]), dpi=600)
        if zoom:
            plt.xlim(axis_lim_x_zoom[0], axis_lim_x_zoom[1])
        plt.savefig("{}/spin_train_zoom{}.png".format(model_spin[0], transition_time_plot), dpi=600)

# Plot spin moment validation
if spin_is_population:
    if plot_spin_time:
        fig_spin_valid, axes_spin_valid = plt.subplots()
        plot_spin_time1((dft_s[1][:, :, 0] - dft_s[1][:, :, 1]), (spin_1[1][:, :, 0] - spin_1[1][:, :, 1]),
                        axes_spin_valid, axis_lim_y, num_atoms=num_atoms, title="Energy, No-aparam")
        plt.tight_layout()
        plt.savefig("{}/spin_valid.png".format(model_spin[0]), dpi=600)
        if zoom:
            plt.xlim(axis_lim_x_zoom[0], axis_lim_x_zoom[1])
        plt.savefig("{}/spin_valid_zoom{}.png".format(model_spin[0], transition_time_plot), dpi=600)

# Plot spin moment test
if spin_is_population:
    if plot_spin_time:
        fig_spin_test, axes_spin_test = plt.subplots()
        plot_spin_time1((dft_s[2][:, :, 0] - dft_s[2][:, :, 1]), (spin_1[2][:, :, 0] - spin_1[2][:, :, 1]),
                        axes_spin_test, axis_lim_y, num_atoms=num_atoms, title="Energy, No-aparam")
        plt.tight_layout()
        plt.savefig("{}/spin_test.png".format(model_spin[0]), dpi=600)
        if zoom:
            plt.xlim(axis_lim_x_zoom[0], axis_lim_x_zoom[1])
        plt.savefig("{}/spin_test_zoom{}.png".format(model_spin[0], transition_time_plot), dpi=600)

# Plot alpha population training
# axis_lim_y = [3.2, 3.52]
# fig_population_alpha_train, axes_population_alpha_train = plt.subplots()
# plot_spin_time1(dft_s[0][:, :, 0], spin_1[0][:, :, 0], axes_population_alpha_train, axis_lim_y, title="Energy, No-aparam")
# plt.tight_layout()
# plt.savefig("{}/population_alpha_train.png".format(model_spin[0]), dpi=600)
# if zoom:
#     plt.xlim(axis_lim_x_zoom[0], axis_lim_x_zoom[1])
# plt.savefig("{}/population_alpha_train_zoom_{}.png".format(model_spin[0], transition_time_plot), dpi=600)
#
# # Plot beta population training
# axis_lim_y = [3.25, 2.65]
# fig_population_beta_train, axes_population_beta_train = plt.subplots()
# plot_spin_time1(dft_s[0][:, :, 1], spin_1[0][:, :, 1], axes_population_beta_train, axis_lim_y, title="Energy, No-aparam")
# plt.tight_layout()
# plt.savefig("{}/population_beta_train.png".format(model_spin[0]), dpi=600)
# if zoom:
#     plt.xlim(axis_lim_x_zoom[0], axis_lim_x_zoom[1])
# plt.savefig("{}/population_beta_train_zoom_{}.png".format(model_spin[0], transition_time_plot), dpi=600)

if __name__ == "__main__":
    print('Finished.')
    plt.show()
