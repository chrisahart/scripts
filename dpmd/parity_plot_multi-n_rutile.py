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
        # ax.text(pos[0], pos[1], f"{text} MAE: {mae:.2f}", transform=ax.transAxes)
        # ax.text(pos[0], pos[1], f"{text} RMSE: {rmse:.2f}", transform=ax.transAxes, color=color_plot)
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
    time_array = np.linspace(0, int(num_timesteps / 2), num=num_timesteps)
    # num_atoms_plot_spin = 64
    # num_atoms_plot_spin = int(num_atoms/3 * 2)
    # num_atoms_plot_spin = 192
    num_atoms_plot_spin = num_atoms
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


def plot_spin_time1_total(dft, dp, ax, axis_lim_y, num_atoms, title=None):
    num_atoms = int(dp.shape[1])
    num_timesteps = int(dp.shape[0])
    time_array = np.linspace(0, int(num_timesteps / 2), num=num_timesteps)
    plotting_colors = ['r', 'g', 'b', 'm', 'grey', 'orange', 'brown', 'hotpink'] * 100

    dft_total_spin = np.sum(dft, axis=1)
    dp_total_spin = np.sum(dp, axis=1)

    mae = mean_absolute_error(dft_total_spin.flatten(), dp_total_spin.flatten())
    rmse = root_mean_squared_error(dft_total_spin.flatten(), dp_total_spin.flatten())
    print('plot_spin_time1_total')
    print('plot_spin')
    print('mean_absolute_error', mae)
    print('root_mean_squared_error', rmse)

    ax.plot(time_array, dft_total_spin, '-', color=plotting_colors[0])
    ax.plot(time_array, dp_total_spin, '-', color=plotting_colors[1])

    ax.set_xlim(0, time_array.shape[0]/2)
    ax.set_ylim([axis_lim_y[0], axis_lim_y[1]])
    ax.set_xlabel("Time / fs")
    ax.set_ylabel("Total spin moment")
    # ax.legend()


# do_concatenate = False

# Bulk TiO2 leopold leopold-fixed2
# model = ['single-fit-ener-dpa3-nlayers-6-official-v3.1.0-start_pref-0.02-1000_limit_pref-1-1',
#          'single-fit-pop-dpa3-nlayers-6-official-v3.1.0-dev-polaron-loss-mae-pref-1-pref_pop-1000-1']
# spin_is_population = True
# folder = '/Volumes/Elements/Data/Postdoc2/Data/Work/calculations/tio2/rutile/deepmd/rutile/leopold-fixed2'
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
# num_atoms = 288

# Bulk TiO2 leopold leopold-fixed3
# model = ['single-fit-ener-dpa3-nlayers-6-official-v3.1.0-start_pref-0.02-1000_limit_pref-1-1',
#          'single-fit-pop-dpa3-nlayers-6-official-v3.1.0-dev-polaron-loss-mae-pref-1-pref_pop-1000-1']
# spin_is_population = True
# folder = '/Volumes/Elements/Data/Postdoc2/Data/Work/calculations/tio2/rutile/deepmd/rutile/leopold-fixed3'
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
# num_atoms = 288

# Bulk TiO2 336 hse-22-10-ps-train-5-ps
# model = ['single-fit-ener-dpa3-nlayers-6-official-v3.1.0-start_pref-0.02-1000_limit_pref-1-1',
#          'single-fit-pop-dpa3-nlayers-6-official-v3.1.0-dev-polaron-loss-mae-pref-1-pref_pop-1000-1']  # energy: 0.205, 0.204, 0.205
# model = ['single-fit-ener-dpa3-nlayers-6-official-v3.1.0-start_pref-0.02-1000_limit_pref-1-1-data_stat_nbatch-2',
#          'single-fit-pop-dpa3-nlayers-6-official-v3.1.0-dev-polaron-loss-mae-pref-1-pref_pop-1000-1']  # energy: 0.161, 0.160, 0.161
# model = ['single-fit-ener-dpa3-nlayers-6-official-v3.1.0-start_pref-0.02-1000_limit_pref-1-1-data_stat_nbatch-3',
#          'single-fit-pop-dpa3-nlayers-6-official-v3.1.0-dev-polaron-loss-mae-pref-1-pref_pop-1000-1']  # energy: 0.060, 0.059, 0.062
# model = ['single-fit-ener-dpa3-nlayers-6-official-v3.1.0-start_pref-0.02-1000_limit_pref-1-1-data_stat_nbatch-4',
#          'single-fit-pop-dpa3-nlayers-6-official-v3.1.0-dev-polaron-loss-mae-pref-1-pref_pop-1000-1']  # energy: 0.191, 0.190, 0.189
# model = ['single-fit-ener-dpa3-nlayers-6-official-v3.1.0-start_pref-0.02-1000_limit_pref-1-1-data_stat_nbatch-5',
#          'single-fit-pop-dpa3-nlayers-6-official-v3.1.0-dev-polaron-loss-mae-pref-1-pref_pop-1000-1']  # energy: 0.1078, 0.077, 0.079
# model = ['single-fit-ener-dpa3-nlayers-6-official-v3.1.0-start_pref-0.02-1000_limit_pref-1-1-data_stat_nbatch-5000',
#          'single-fit-pop-dpa3-nlayers-6-official-v3.1.0-dev-polaron-loss-mae-pref-1-pref_pop-1000-1']  # energy: 0.279, 0.277, 0.280
# model = ['single-fit-ener-dpa3-nlayers-6-official-v3.1.0-start_pref-0.02-1000_limit_pref-1-1-twostep',
#          'single-fit-pop-dpa3-nlayers-6-official-v3.1.0-dev-polaron-loss-mae-pref-1-pref_pop-1000-1']  # energy: 0.209, 0.211, 0.217
# model = ['single-fit-ener-dpa3-nlayers-6-official-v3.1.0-start_pref-0.02-1000_limit_pref-1-1-twostep-10-1-10-1',
#          'single-fit-pop-dpa3-nlayers-6-official-v3.1.0-dev-polaron-loss-mae-pref-1-pref_pop-1000-1']  # energy: 0.057, 0.058, 0.074
# model = ['single-fit-ener-dpa3-nlayers-6-official-v3.1.0-start_pref-0.02-1000_limit_pref-1-1-twostep-lr-1e-5-1e-8',
#          'single-fit-pop-dpa3-nlayers-6-official-v3.1.0-dev-polaron-loss-mae-pref-1-pref_pop-1000-1']  # energy: 0.043, 0.044, 0.050. force: 40, 40, 42   ***
#
# model = ['single-fit-ener-dpa3-nlayers-6-official-v3.1.0-start_pref-0.2-100_limit_pref-20-60',
#          'single-fit-pop-dpa3-nlayers-6-official-v3.1.0-dev-polaron-loss-mae-pref-1-pref_pop-1000-1']  # energy: 0.115, 0.115, 0.121
# model = ['single-fit-ener-dpa3-nlayers-6-official-v3.1.0-start_pref-0.2-100_limit_pref-20-60-data_stat_nbatch-3',
#          'single-fit-pop-dpa3-nlayers-6-official-v3.1.0-dev-polaron-loss-mae-pref-1-pref_pop-1000-1']  # energy: 0.099, 0.101, 0.098
# model = ['single-fit-ener-dpa3-nlayers-6-official-v3.1.0-start_pref-0.2-100_limit_pref-20-60-twostep',
#          'single-fit-pop-dpa3-nlayers-6-official-v3.1.0-dev-polaron-loss-mae-pref-1-pref_pop-1000-1']  # energy: 0.075, 0.077, 0.086
# model = ['single-fit-ener-dpa3-nlayers-6-official-v3.1.0-start_pref-0.2-100_limit_pref-20-60-twostep-10-1-10-1',
#          'single-fit-pop-dpa3-nlayers-6-official-v3.1.0-dev-polaron-loss-mae-pref-1-pref_pop-1000-1']  # energy: 0.074, 0.075, 0.084
# model = ['single-fit-ener-dpa3-nlayers-6-official-v3.1.0-start_pref-0.2-100_limit_pref-20-60-twostep-lr-1e-5-1e-8',
#          'single-fit-pop-dpa3-nlayers-6-official-v3.1.0-dev-polaron-loss-mae-pref-1-pref_pop-1000-1']  # energy: 0.048, 0.049, 0.065. force 36, 37, 46

# model = ['single-fit-ener-dpa3-nlayers-6-official-v3.1.0-start_pref-1-1_limit_pref-1-1',
#          'single-fit-pop-dpa3-nlayers-6-official-v3.1.0-dev-polaron-loss-mae-pref-1-pref_pop-1000-1']  # energy: 0.230, 0.229, 0.234
# model = ['single-fit-ener-dpa3-nlayers-6-official-v3.1.0-start_pref-1-1_limit_pref-1-1-data_stat_nbatch-3',
#          'single-fit-pop-dpa3-nlayers-6-official-v3.1.0-dev-polaron-loss-mae-pref-1-pref_pop-1000-1']  # energy: 0.149, 0.147, 0.156
# model = ['single-fit-ener-dpa3-nlayers-6-official-v3.1.0-start_pref-1-1_limit_pref-1-1-twostep',
#          'single-fit-pop-dpa3-nlayers-6-official-v3.1.0-dev-polaron-loss-mae-pref-1-pref_pop-1000-1']  # energy: 0.121, 0.122, 0.139
# model = ['single-fit-ener-dpa3-nlayers-6-official-v3.1.0-start_pref-1-1_limit_pref-1-1-twostep-lr-1e-5-1e-8',
#          'single-fit-pop-dpa3-nlayers-6-official-v3.1.0-dev-polaron-loss-mae-pref-1-pref_pop-1000-1']  # energy: 0.119, 0.120, 0.123. force: 39, 39, 44

# model = ['single-fit-ener-dpa3-nlayers-6-official-v3.1.0-start_pref-0.02-1000_limit_pref-1-1-data_stat_nbatch-3',
#          'single-fit-pop-dpa3-nlayers-6-official-v3.1.0-dev-polaron-loss-mae-pref-1-pref_pop-1000-1']  # spin: 0.011, 0.012, 0.024
# model = ['single-fit-ener-dpa3-nlayers-6-official-v3.1.0-start_pref-0.02-1000_limit_pref-1-1-data_stat_nbatch-3',
         # 'single-fit-pop-dpa3-nlayers-6-official-v3.1.0-dev-polaron-loss-mae-pref-1-pref_pop-10000-10']  # spin: 0.010, 0.011, 0.024
# model = ['single-fit-ener-dpa3-nlayers-6-official-v3.1.0-start_pref-0.02-1000_limit_pref-1-1-data_stat_nbatch-3',
#          'single-fit-pop-dpa3-nlayers-6-official-v3.1.0-dev-polaron-loss-mae-pref-1-pref_pop-10000-10-twostep-lr-1e-5-1e-8']  # spin: 0.012, 0.014, 0.026 total spin slightly better

# spin_is_population = True
# folder = '/Volumes/Samsung/Data/Postdoc2/Data/Work/calculations/tio2/rutile/deepmd/rutile/336/md-cell-opt/hse-22-10-ps-train-5-ps'
# model_ener = ['{}/{}'.format(folder, model[0])] * 4
# model_spin = ['{}/{}'.format(folder, model[1])] * 4
# database = ['{}/database_population_train/'.format(folder),
#             '{}/database_population_test/1'.format(folder),
#             '{}/database_population_test/2'.format(folder)]
# val = ['_0_', '_1_', '_2_']
# database = ['{}/database_hse-22-10-ps-train-5-ps/database_population_train/'.format(folder),
#             '{}/database_hse-22-10-ps-train-5-ps/database_population_test/1'.format(folder),
#             '{}/database_hse-22-10-ps-train-5-ps/database_population_test/2'.format(folder)]
# val = ['_0_', '_1_', '_2_']
# database = ['{}/database_hse-22-10-ps-train-9-ps/database_population_train/'.format(folder),
#             '{}/database_hse-22-10-ps-train-9-ps/database_population_test/1'.format(folder),
#             '{}/database_hse-22-10-ps-train-9-ps/database_population_test/2'.format(folder)]
# val = ['_3_', '_4_', '_5_']
# axis_lim_y = np.array([-0.02, 0.805])
# pos_array_energy = np.array(([0.45, 0.25], [0.45, 0.15], [0.45, 0.05]))
# pos_array_force = np.array(([0.6, 0.25], [0.6, 0.15], [0.6, 0.05]))
# pos_array_spin = np.array(([0.72, 0.25], [0.72, 0.15], [0.72, 0.05]))
# text_array = ['Train', 'Valid', 'Test']
# num_atoms = 324

# Bulk TiO2 336 hse-22-10-ps-train-9-ps
# model = ['single-fit-ener-dpa3-nlayers-6-official-v3.1.0-start_pref-0.02-1000_limit_pref-1-1',
#          'single-fit-pop-dpa3-nlayers-6-official-v3.1.0-dev-polaron-loss-mae-pref-1-pref_pop-1000-1']  # decent, force 40, spin 0.014
# spin_is_population = True
# folder = '/Volumes/Samsung/Data/Postdoc2/Data/Work/calculations/tio2/rutile/deepmd/rutile/336/md-cell-opt/hse-22-10-ps-train-9-ps'
# model_ener = ['{}/{}'.format(folder, model[0])] * 4
# model_spin = ['{}/{}'.format(folder, model[1])] * 4
# database = ['{}/database_population_train/'.format(folder),
#             '{}/database_population_test/2'.format(folder)]
# val = ['_0_', '_1_']
# axis_lim_y = np.array([-0.02, 0.805])
# pos_array_energy = np.array(([0.45, 0.25], [0.45, 0.15], [0.45, 0.05]))
# pos_array_force = np.array(([0.6, 0.25], [0.6, 0.15], [0.6, 0.05]))
# pos_array_spin = np.array(([0.72, 0.25], [0.72, 0.15], [0.72, 0.05]))
# text_array = ['Train', 'Valid', 'Test']
# num_atoms = 324

# Bulk TiO2 336 hse-22-ts-md2
# model = ['single-fit-ener-dpa3-nlayers-6-official-v3.1.0-start_pref-0.02-1000_limit_pref-1-1-twostep-lr-1e-5-1e-8',
#          'single-fit-pop-dpa3-nlayers-6-official-v3.1.0-dev-polaron-loss-mae-pref-1-pref_pop-1000-1']  # Good, forces slightly overfit
# model = ['single-fit-ener-dpa3-nlayers-6-official-v3.1.0-start_pref-0.02-1000_limit_pref-1-1-rcut-4',
#          'single-fit-pop-dpa3-nlayers-6-official-v3.1.0-dev-polaron-loss-mae-pref-1-pref_pop-1000-1']  # md 5 ps: 0.125, 44
# model = ['single-fit-ener-dpa3-nlayers-6-official-v3.1.0-start_pref-1-1_limit_pref-1-1-twostep-lr-1e-5-1e-8',
#          'single-fit-pop-dpa3-nlayers-6-official-v3.1.0-dev-polaron-loss-mae-pref-1-pref_pop-1000-1']  # Good, forces slightly overfit

# model = ['single-fit-ener-dpa3-nlayers-6-official-v3.1.0-start_pref-0.02-1000_limit_pref-1-1-twostep-lr-1e-5-1e-8',
#          'single-fit-pop-dpa3-nlayers-6-official-v3.1.0-dev-polaron-loss-mae-pref-1-pref_pop-1000-1']  # md 5 ps: spin 0.003 polaron 0.027
# model = ['single-fit-ener-dpa3-nlayers-6-official-v3.1.0-start_pref-0.02-1000_limit_pref-1-1-twostep-lr-1e-5-1e-8',
#          'single-fit-pop-dpa3-nlayers-6-official-v3.1.0-dev-polaron-loss-mae-pref-1-pref_pop-1000-1-lr-1e-5-1e-8']  # md 5 ps: spin 0.003 polaron 0.027
# model = ['single-fit-ener-dpa3-nlayers-6-official-v3.1.0-start_pref-0.02-1000_limit_pref-1-1-twostep-lr-1e-5-1e-8',
#          'single-fit-pop-dpa3-nlayers-6-official-v3.1.0-dev-polaron-loss-mae-pref-1-pref_pop-1000-1-twostep-lr-1e-5-1e-8']  # md 5 ps: spin 0.003 polaron 0.027

# model = ['single-fit-ener-dpa3-nlayers-6-official-v3.1.0-start_pref-0.02-1000_limit_pref-1-1-twostep-lr-1e-5-1e-8',
#          'single-fit-pop-dpa3-nlayers-6-official-v3.1.0-dev-polaron-loss-mae-pref-1-pref_pop-10000-10']  # md 5 ps: spin 0.003 polaron 0.027
# model = ['single-fit-ener-dpa3-nlayers-6-official-v3.1.0-start_pref-0.02-1000_limit_pref-1-1-twostep-lr-1e-5-1e-8',
#          'single-fit-pop-dpa3-nlayers-6-official-v3.1.0-dev-polaron-loss-mae-pref-1-pref_pop-10000-10-lr-1e-5-1e-8']  # md 5 ps: spin 0.003 polaron 0.027
# model = ['single-fit-ener-dpa3-nlayers-6-official-v3.1.0-start_pref-0.02-1000_limit_pref-1-1-twostep-lr-1e-5-1e-8',
#          'single-fit-pop-dpa3-nlayers-6-official-v3.1.0-dev-polaron-loss-mae-pref-1-pref_pop-10000-10-twostep-lr-1e-5-1e-8']  # md 5 ps: spin 0.003 polaron 0.027

# model = ['single-fit-ener-dpa3-nlayers-6-official-v3.1.0-start_pref-0.02-1000_limit_pref-1-1-twostep-lr-1e-5-1e-8',
#          'single-fit-pop-dpa3-nlayers-6-official-v3.1.0-dev-polaron-loss-mae-pref-1-pref_pop-1000-1-twostep-lr-1e-5-1e-8']  # DP-MD 1
# model = ['single-fit-ener-dpa3-nlayers-6-official-v3.1.0-start_pref-0.02-1000_limit_pref-1-1',
#          'single-fit-pop-dpa3-nlayers-6-official-v3.1.0-dev-polaron-loss-mae-pref-1-pref_pop-1000-1-twostep-lr-1e-5-1e-8']  # DP-MD 2

# model = ['single-fit-ener-dpa3-nlayers-6-official-v3.1.0-start_pref-0.02-1000_limit_pref-1-1-rcut-4',
#          'single-fit-pop-dpa3-nlayers-6-official-v3.1.0-dev-polaron-loss-mae-pref-1-pref_pop-1000-1-twostep-lr-1e-5-1e-8']  # DP-MD 2

# model = ['single-fit-ener-dpa3-nlayers-6-official-v3.1.0-start_pref-0.02-1000_limit_pref-1-1-fixed',
#          'single-fit-pop-dpa3-nlayers-6-official-v3.1.0-dev-polaron-loss-mae-pref-1-pref_pop-1000-1-twostep-lr-1e-5-1e-8']  # good
# model = ['single-fit-ener-dpa3-nlayers-6-official-v3.1.0-start_pref-0.02-1000_limit_pref-1-1-fixed-decay_steps-50000',
#          'single-fit-pop-dpa3-nlayers-6-official-v3.1.0-dev-polaron-loss-mae-pref-1-pref_pop-1000-1-twostep-lr-1e-5-1e-8']  # bad
# model = ['single-fit-ener-dpa3-nlayers-6-official-v3.1.0-start_pref-0.02-1000_limit_pref-1-1-fixed-lr-1e-5-1e-8',
#          'single-fit-pop-dpa3-nlayers-6-official-v3.1.0-dev-polaron-loss-mae-pref-1-pref_pop-1000-1-twostep-lr-1e-5-1e-8']  # bad


# Bulk TiO2 336 hse-22-ts-md3
# folder = '/Volumes/Elements/Data/Postdoc2/Data/Work/calculations/tio2/rutile/deepmd/rutile/336/md-cell-opt/hse-22-ts-md3'
# model = ['ener-dpa3-start_pref-0.02-1000_limit_pref-1-1', 'pop-dpa3-pref_pop-1000-1']
# model = ['ener-dpa3-start_pref-0.02-1000_limit_pref-1-1-twostep-lr-1e-5-1e-8', 'pop-dpa3-pref_pop-1000-1-twostep-lr-1e-5-1e-8']
# model = ['ener-dpa3-start_pref-1-1_limit_pref-1-1', 'pop-dpa3-pref_pop-10000-10']
# model = ['ener-dpa3-start_pref-1-1_limit_pref-1-1-twostep-lr-1e-5-1e-8', 'pop-dpa3-pref_pop-10000-10-twostep-lr-1e-5-1e-8']
# model = ['ener-dpa3-start_pref-0.02-1000_limit_pref-1-1-twostep-lr-1e-5-1e-8', 'pop-dpa3-pref_pop-10000-10-twostep-lr-1e-5-1e-8']  # Best
# model_ener = ['{}/{}'.format(folder, model[0])] * 6
# model_spin = ['{}/{}'.format(folder, model[1])] * 6
# database = ['{}/database_ts/database_population_train/'.format(folder),
#             '{}/database_ts/database_population_test/1'.format(folder),
#             '{}/database_ts/database_population_test/2'.format(folder)]
# val = ['_0_', '_1_', '_2_']
# database += ['{}/database_md/database_population_train/'.format(folder),
#              '{}/database_md/database_population_test/1'.format(folder),
#              '{}/database_md/database_population_test/2'.format(folder)]
# val += ['_3_', '_4_', '_5_']
# spin_is_population = True
# do_concatenate = True
# axis_lim_y = np.array([-0.02, 0.805])
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

for i in range(len(database)):
    dft_e.append(np.load("{}/set.000/energy.npy".format(database[i], val[i])))
    dft_f.append(np.load("{}/set.000/force.npy".format(database[i], val[i])))
    dft_polaron.append(np.load("{}/set.000/aparam.npy".format(database[i], val[i])))
    # ener_0.append(np.load("{}/0{}ener.npy".format(model_ener[i], val[i])))
    # force_0.append(np.load("{}/0{}force.npy".format(model_ener[i], val[i])))
    ener_1.append(np.load("{}/1{}ener.npy".format(model_ener[i], val[i])))
    force_1.append(np.load("{}/1{}force.npy".format(model_ener[i], val[i])))
    # spin_0.append(np.load("{}/0{}spin.npy".format(model_spin[i], val[i])))
    spin_1.append(np.load("{}/1{}spin.npy".format(model_spin[i], val[i])))
    if not spin_is_population:
        dft_s.append(np.load("{}/set.000/atom_ener.npy".format(database[i], val[i])))
    if spin_is_population:
        dft_s.append(np.load("{}/set.000/atomic_population.npy".format(database[i], val[i])))

# Concatenate
if do_concatenate:
    merged_dft_e = []
    merged_dft_f = []
    merged_dft_polaron = []
    merged_ener_1 = []
    merged_force_1 = []
    merged_spin_1 = []
    merged_dft_s = []
    for i in range(3):
        merged_dft_e.append(np.concatenate((dft_e[i], dft_e[i+3])))
        merged_dft_f.append(np.concatenate((dft_f[i], dft_f[i+3])))
        merged_dft_polaron.append(np.concatenate((dft_polaron[i], dft_polaron[i+3])))
        merged_ener_1.append(np.concatenate((ener_1[i], ener_1[i+3])))
        merged_force_1.append(np.concatenate((force_1[i], force_1[i+3])))
        merged_spin_1.append(np.concatenate((spin_1[i], spin_1[i+3])))
        merged_dft_s.append(np.concatenate((dft_s[i], dft_s[i+3])))
    dft_e = merged_dft_e
    dft_f = merged_dft_f
    dft_polaron = merged_dft_polaron
    ener_1 = merged_ener_1
    force_1 = merged_force_1
    spin_1 = merged_spin_1
    dft_s = merged_dft_s
    database = ['{}/database_ts/database_population_train/'.format(folder),
                '{}/database_ts/database_population_test/1'.format(folder),
                '{}/database_ts/database_population_test/2'.format(folder)]
    val = ['_0_', '_1_', '_2_']

# Plot parity 2x3 subplots
# fig, axes = plt.subplots(2, 3, figsize=(15, 6))
# for i in range(len(database)):
#     print(i)
#     plot_ener(dft_e[i], ener_0[i], axes[0, 0], num_atoms=num_atoms, color_plot=color_plot_array[i], pos=pos_array_energy[i], text=text_array[i], title="Energy, No-aparam")
#     plot_ener(dft_e[i], ener_1[i], axes[1, 0], num_atoms=num_atoms, color_plot=color_plot_array[i], pos=pos_array_energy[i], text=text_array[i], title="Energy, Yes-aparam")
#     plot_force(dft_f[i], force_0[i], axes[0, 1], color_plot=color_plot_array[i], pos=pos_array_force[i], text=text_array[i], title="Force, No-aparam")
#     plot_force(dft_f[i], force_1[i], axes[1, 1], color_plot=color_plot_array[i], pos=pos_array_force[i], text=text_array[i], title="Force, Yes-aparam")
#     if not spin_is_population:
#         plot_spin(dft_s[i], spin_0[i], axes[0, 2], color_plot=color_plot_array[i], pos=pos_array_spin[i], text=text_array[i], title="Spin, No-aparam")
#         plot_spin(dft_s[i], spin_1[i], axes[1, 2], color_plot=color_plot_array[i], pos=pos_array_spin[i], text=text_array[i], title="Spin, Yes-aparam")
#     if spin_is_population:
#         plot_spin((dft_s[i][:, :, 0] - dft_s[i][:, :, 1]), (spin_0[i][:, :, 0] - spin_0[i][:, :, 1]), axes[0, 2], color_plot=color_plot_array[i], pos=pos_array_spin[i], text=text_array[i], title="Spin, No-aparam")
#         plot_spin((dft_s[i][:, :, 0] - dft_s[i][:, :, 1]), (spin_1[i][:, :, 0] - spin_1[i][:, :, 1]), axes[1, 2], color_plot=color_plot_array[i], pos=pos_array_spin[i], text=text_array[i], title="Spin, Yes-aparam")
# plt.tight_layout()
# for i in range(len(database)):
#     plt.savefig("{}/fit_2x3_folders_{}.png".format(model_spin[i], len(model_ener)), dpi=600)
#     plt.savefig("{}/fit_2x3_folders_{}.png".format(model_ener[i], len(model_spin)), dpi=600)

# Plot parity 1x3 subplot
print('Plot parity 1x3 subplot')
fig2, axes2 = plt.subplots(1, 3, figsize=(15, 5))
for i in range(len(database)):
    print(i)
    plot_ener(dft_e[i], ener_1[i], axes2[0], num_atoms=num_atoms, color_plot=color_plot_array[i],
              pos=pos_array_energy[i], text=text_array[i])
    plot_force(dft_f[i], force_1[i], axes2[1], color_plot=color_plot_array[i], pos=pos_array_force[i],
               text=text_array[i])
    if not spin_is_population:
        plot_spin(dft_s[i], spin_1[i], axes2[2], color_plot=color_plot_array[i], pos=pos_array_spin[i],
                  text=text_array[i])
    if spin_is_population:
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

# plot_energy_time = False
# plot_force_time = False
plot_energy_time = True
plot_force_time = True

if plot_energy_time:
    # Plot energy time training
    print('Plot energy time training')
    num_timesteps = int(dft_s[0][:, :, 0].shape[0])
    time_array = np.linspace(0, int(num_timesteps / 2), num=num_timesteps)
    fig_energy_time_train, axes_energy_time_train = plt.subplots()
    axes_energy_time_train.plot(time_array, (ener_1[0] - dft_e[0]) / num_atoms * 1000)
    axes_energy_time_train.set_xlim(0, time_array.shape[0] / 2)
    axes_energy_time_train.set_xlabel("Time / fs")
    axes_energy_time_train.set_ylabel("Energy DP - DFT / meV")
    plt.tight_layout()
    plt.savefig("{}/energy_time_train.png".format(model_spin[0]), dpi=600)

    # Plot energy time validation
    num_timesteps = int(dft_s[1][:, :, 0].shape[0])
    time_array = np.linspace(0, int(num_timesteps / 2), num=num_timesteps)
    fig_energy_time_valid, axes_energy_time_valid = plt.subplots()
    axes_energy_time_valid.plot(time_array, (ener_1[1] - dft_e[1]) / num_atoms * 1000)
    axes_energy_time_valid.set_xlim(0, time_array.shape[0] / 2)
    axes_energy_time_valid.set_xlabel("Time / fs")
    axes_energy_time_valid.set_ylabel("Energy DP - DFT / meV")
    plt.tight_layout()
    plt.savefig("{}/energy_time_valid.png".format(model_spin[0]), dpi=600)

    # Plot energy time testing
    num_timesteps = int(dft_s[2][:, :, 0].shape[0])
    time_array = np.linspace(0, int(num_timesteps / 2), num=num_timesteps)
    fig_energy_time_test, axes_energy_time_test = plt.subplots()
    axes_energy_time_test.plot(time_array, (ener_1[2] - dft_e[2]) / num_atoms * 1000)
    axes_energy_time_test.set_xlim(0, time_array.shape[0] / 2)
    axes_energy_time_test.set_xlabel("Time / fs")
    axes_energy_time_test.set_ylabel("Energy DP - DFT / meV")
    plt.tight_layout()
    plt.savefig("{}/energy_time_test.png".format(model_spin[0]), dpi=600)

if plot_force_time:
    # Plot force time training
    print('Plot force time training')
    dft_f_train = np.reshape(dft_f[0], (dft_f[0].shape[0], num_atoms, 3))
    force_error_train1 = force_1[0] - dft_f_train
    force_error_train2 = np.max(force_error_train1, axis=2)
    force_error_train3 = np.max(force_error_train2, axis=1)
    num_timesteps = int(dft_s[0][:, :, 0].shape[0])
    time_array = np.linspace(0, int(num_timesteps / 2), num=num_timesteps)
    fig_force_time_train, axes_force_time_train = plt.subplots()
    axes_force_time_train.plot(time_array, force_error_train3)
    axes_force_time_train.set_xlim(0, time_array.shape[0] / 2)
    axes_force_time_train.set_xlabel("Time / fs")
    axes_force_time_train.set_ylabel("Max force error / meV")
    plt.tight_layout()
    plt.savefig("{}/force_time_train.png".format(model_spin[0]), dpi=600)

    # Plot force time validing
    dft_f_valid = np.reshape(dft_f[1], (dft_f[1].shape[0], num_atoms, 3))
    force_error_valid1 = force_1[1] - dft_f_valid
    force_error_valid2 = np.max(force_error_valid1, axis=2)
    force_error_valid3 = np.max(force_error_valid2, axis=1)
    num_timesteps = int(dft_s[1][:, :, 0].shape[0])
    time_array = np.linspace(0, int(num_timesteps / 2), num=num_timesteps)
    fig_force_time_valid, axes_force_time_valid = plt.subplots()
    axes_force_time_valid.plot(time_array, force_error_valid3)
    axes_force_time_valid.set_xlim(0, time_array.shape[0] / 2)
    axes_force_time_valid.set_xlabel("Time / fs")
    axes_force_time_valid.set_ylabel("Max force error / meV")
    plt.tight_layout()
    plt.savefig("{}/force_time_valid.png".format(model_spin[0]), dpi=600)

    # Plot force time testing
    dft_f_test = np.reshape(dft_f[2], (dft_f[2].shape[0], num_atoms, 3))
    force_error_test1 = force_1[2] - dft_f_test
    force_error_test2 = np.max(force_error_test1, axis=2)
    force_error_test3 = np.max(force_error_test2, axis=1)
    num_timesteps = int(dft_s[2][:, :, 0].shape[0])
    time_array = np.linspace(0, int(num_timesteps / 2), num=num_timesteps)
    fig_force_time_test, axes_force_time_test = plt.subplots()
    axes_force_time_test.plot(time_array, force_error_test3)
    axes_force_time_test.set_xlim(0, time_array.shape[0] / 2)
    axes_force_time_test.set_xlabel("Time / fs")
    axes_force_time_test.set_ylabel("Max force error / meV")
    plt.tight_layout()
    plt.savefig("{}/force_time_test.png".format(model_spin[0]), dpi=600)

# Plot spin moment training
# axis_lim_y = [-0.02, 1.0]
plot_spin_time = True
# plot_spin_time = False
if plot_spin_time:
    print('fig_spin_train')
    fig_spin_train, axes_spin_train = plt.subplots()
    plot_spin_time1((dft_s[0][:, :, 0] - dft_s[0][:, :, 1]), (spin_1[0][:, :, 0] - spin_1[0][:, :, 1]),
                    axes_spin_train, axis_lim_y, num_atoms=num_atoms, title="Energy, No-aparam")
    plt.tight_layout()
    plt.savefig("{}/spin_train.png".format(model_spin[0]), dpi=600)
    if zoom:
        plt.xlim(axis_lim_x_zoom[0], axis_lim_x_zoom[1])
    plt.savefig("{}/spin_train_zoom{}.png".format(model_spin[0], transition_time_plot), dpi=600)

    print('fig_spin_train_total')
    fig_spin_train_total, axes_spin_train_total = plt.subplots()
    plot_spin_time1_total((dft_s[0][:, :, 0] - dft_s[0][:, :, 1]), (spin_1[0][:, :, 0] - spin_1[0][:, :, 1]),
                          axes_spin_train_total, [0, 1.2], num_atoms=num_atoms, title="Energy, No-aparam")
    plt.tight_layout()
    plt.savefig("{}/spin_train_total.png".format(model_spin[0]), dpi=600)

# Plot spin moment validation
if plot_spin_time:
    print('fig_spin_valid')
    fig_spin_valid, axes_spin_valid = plt.subplots()
    plot_spin_time1((dft_s[1][:, :, 0] - dft_s[1][:, :, 1]), (spin_1[1][:, :, 0] - spin_1[1][:, :, 1]),
                    axes_spin_valid, axis_lim_y, num_atoms=num_atoms, title="Energy, No-aparam")
    plt.tight_layout()
    plt.savefig("{}/spin_valid.png".format(model_spin[0]), dpi=600)
    if zoom:
        plt.xlim(axis_lim_x_zoom[0], axis_lim_x_zoom[1])
    plt.savefig("{}/spin_valid_zoom{}.png".format(model_spin[0], transition_time_plot), dpi=600)

    print('fig_spin_valid_total')
    fig_spin_valid_total, axes_spin_valid_total = plt.subplots()
    plot_spin_time1_total((dft_s[1][:, :, 0] - dft_s[1][:, :, 1]), (spin_1[1][:, :, 0] - spin_1[1][:, :, 1]),
                          axes_spin_valid_total, [0, 1.2], num_atoms=num_atoms, title="Energy, No-aparam")
    plt.tight_layout()
    plt.savefig("{}/spin_valid_total.png".format(model_spin[0]), dpi=600)

# Plot spin moment test
if plot_spin_time:
    print('fig_spin_test')
    fig_spin_test, axes_spin_test = plt.subplots()
    plot_spin_time1((dft_s[2][:, :, 0] - dft_s[2][:, :, 1]), (spin_1[2][:, :, 0] - spin_1[2][:, :, 1]),
                    axes_spin_test, axis_lim_y, num_atoms=num_atoms, title="Energy, No-aparam")
    plt.tight_layout()
    plt.savefig("{}/spin_test.png".format(model_spin[0]), dpi=600)
    if zoom:
        plt.xlim(axis_lim_x_zoom[0], axis_lim_x_zoom[1])
    plt.savefig("{}/spin_test_zoom{}.png".format(model_spin[0], transition_time_plot), dpi=600)

    print('fig_spin_train_total')
    fig_spin_train_total, axes_spin_train_total = plt.subplots()
    plot_spin_time1_total((dft_s[2][:, :, 0] - dft_s[2][:, :, 1]), (spin_1[2][:, :, 0] - spin_1[2][:, :, 1]),
                          axes_spin_train_total, [0, 1.2], num_atoms=num_atoms, title="Energy, No-aparam")
    plt.tight_layout()
    plt.savefig("{}/spin_train_total.png".format(model_spin[0]), dpi=600)

    # population
    print('fig_population_train_total')
    fig_population_train_total, axes_population_train_total = plt.subplots()
    plot_spin_time1_total((dft_s[0][:, :, 0] + dft_s[0][:, :, 1]), (spin_1[0][:, :, 0] + spin_1[0][:, :, 1]),
                          axes_population_train_total, [2593 - 0.5, 2593 + 0.5], num_atoms=num_atoms,
                          title="Energy, No-aparam")
    plt.tight_layout()
    plt.savefig("{}/population_train_total.png".format(model_spin[0]), dpi=600)

    print('fig_population_valid_total')
    fig_population_valid_total, axes_population_valid_total = plt.subplots()
    plot_spin_time1_total((dft_s[1][:, :, 0] + dft_s[1][:, :, 1]), (spin_1[1][:, :, 0] + spin_1[1][:, :, 1]),
                          axes_population_valid_total, [2593 - 0.5, 2593 + 0.5], num_atoms=num_atoms,
                          title="Energy, No-aparam")
    plt.tight_layout()
    plt.savefig("{}/population_valid_total.png".format(model_spin[0]), dpi=600)

    print('fig_population_test_total')
    fig_population_test_total, axes_population_test_total = plt.subplots()
    plot_spin_time1_total((dft_s[2][:, :, 0] + dft_s[2][:, :, 1]), (spin_1[2][:, :, 0] + spin_1[2][:, :, 1]),
                          axes_population_test_total, [2593 - 0.5, 2593 + 0.5], num_atoms=num_atoms,
                          title="Energy, No-aparam")
    plt.tight_layout()
    plt.savefig("{}/population_test_total.png".format(model_spin[0]), dpi=600)

if __name__ == "__main__":
    print('Finished.')
    plt.show()
