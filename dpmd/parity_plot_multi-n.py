import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import mean_absolute_error
from sklearn.metrics import root_mean_squared_error


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
    ax.text(pos[0], pos[1], f"{text} RMSE: {rmse:.2f} meV/atom", transform=ax.transAxes, color=color_plot)
    if title is not None:
        ax.set_title(title)


def plot_force(dft, dp, ax, color_plot, pos, text, title=None):
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
    ax.text(pos[0], pos[1], f"RMSE: {rmse:.2f} meV/Å", transform=ax.transAxes, color=color_plot)
    min_val, max_val = min(dft.min(), dp.min()), max(dft.max(), dp.max())
    ax.plot([min_val, max_val], [min_val, max_val], 'k--')
    if title is not None:
        ax.set_title(title)
    # else:
    #     ax.set_title("Force")


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
    # ax.text(pos[0], pos[1], f"{text} RMSE: {rmse:.2f}", transform=ax.transAxes, color=color_plot)
    ax.text(pos[0], pos[1], f"RMSE: {rmse:.2f}", transform=ax.transAxes, color=color_plot)
    ax.plot([min_val, max_val], [min_val, max_val], 'k--')

    ax.set_xlim([axis_lim_y[0], axis_lim_y[1]])
    ax.set_ylim([axis_lim_y[0], axis_lim_y[1]])

    if title is not None:
        ax.set_title(title)
    # else:
    #     ax.set_title("Spin")


def plot_spin_time1(dft, dp, ax, axis_lim_y, title=None):
    num_atoms = int(dp.shape[1])
    num_timesteps = int(dp.shape[0])
    time_array = np.linspace(0, int(num_timesteps / 2), num=num_timesteps)
    num_atoms_plot_spin = 64
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


# Bulk hematite MD 400k-f 5 hops training data with 400k-b 1 hop test data single-fit
# model_ener = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/400k-b/test-set/400k-f/single-fit-ener-se_e2_a',
#               '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/400k-b/single-fit-ener-se_e2_a']
# model_spin = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/400k-b/test-set/400k-f/single-fit-m-se_e2_a',
#               '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/400k-b/single-fit-m-se_e2_a']
# model_ener = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/400k-b/test-set/400k-f/single-fit-ener-dpa1-se_atten_v2',
#               '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/400k-b/single-fit-ener-dpa1-se_atten_v2']
# model_spin = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/400k-b/test-set/400k-f/single-fit-m-dpa1-se_atten_v2',
#               '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/400k-b/single-fit-m-dpa1-se_atten_v2']
# model_ener = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/400k-b/test-set/400k-f/single-fit-ener-dpa2',
#               '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/400k-b/single-fit-ener-dpa2']
# model_spin = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/400k-b/test-set/400k-f/single-fit-m-dpa2',
#               '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/400k-b/single-fit-m-dpa2']
# database = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/400k-b/test-set/400k-f/database_spin',
#             '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/400k-b/database_spin']
# axis_lim_y = np.array([-4.05, -3.05])
# pos_array = np.array(([0.45, 0.1], [0.45, 0.2]))

# Bulk hematite MD 400k-f 5 hops training data with 400k-b 1 hop test data multi-fit
# model_ener = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/400k-b/test-set/400k-f/multi-task-se_e2_a',
#               '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/400k-b/multi-task-se_e2_a']
# model_ener = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/400k-b/test-set/400k-f/multi-task-se_e3',
#               '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/400k-b/multi-task-se_e3']
# model_ener = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/400k-b/test-set/400k-f/multi-task-dpa1-se_atten_v2',
#               '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/400k-b/multi-task-dpa1-se_atten_v2']
# model_ener = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/400k-b/test-set/400k-f/multi-task-dpa2',
#               '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/400k-b/multi-task-dpa2']
# model_spin = model_ener
# database = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/400k-b/test-set/400k-f/database_spin',
#             '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/400k-b/database_spin']
# axis_lim_y = np.array([-4.05, -3.05])
# pos_array = np.array(([0.45, 0.1], [0.45, 0.2]))

# Bulk hematite 0K NEB multi-fit 221
# model_ener = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/geo-opt-all/test-set/atom-3-all/multi-task-se_e2_a',
#               '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/geo-opt-all/multi-task-se_e2_a']
# model_ener = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/geo-opt-all/test-set/atom-3-all/multi-task-se_e3',
#               '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/geo-opt-all/multi-task-se_e3']
# model_ener = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/geo-opt-all/test-set/atom-3-all/multi-task-dpa1-se_atten_v2',
#               '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/geo-opt-all/multi-task-dpa1-se_atten_v2']
# model_ener = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/geo-opt-all/test-set/atom-3-all/multi-task-dpa2',
#               '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/geo-opt-all/multi-task-dpa2']
# model_ener = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/geo-opt-all/test-set/atom-3-all/multi-task-dpa2-rcut-5',
#               '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/geo-opt-all/multi-task-dpa2-rcut-5']
# model_ener = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/geo-opt-all/multi-task-dpa2-rcut-5',
#               '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/geo-opt-all/test-set/atom-3-all/multi-task-dpa2-rcut-5']
# model_ener = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/geo-opt-all/test-set/atom-3-all/multi-task-dpa2-rcut-5-limit_pref_e-1',
#               '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/geo-opt-all/multi-task-dpa2-rcut-5-limit_pref_e-1']
# model_spin = model_ener
# database = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/geo-opt-all/test-set/atom-3-all/database_spin',
#             '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/geo-opt-all/database_spin']
# database = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/geo-opt-all/database_spin',
#             '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/geo-opt-all/test-set/atom-3-all/database_spin']
# axis_lim_y = np.array([-4.05, -3.05])
# pos_array = np.array(([0.45, 0.1], [0.45, 0.2]))

# Bulk hematite 0K NEB single fit 221
# model_ener = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/geo-opt-all/test-set/atom-3-all/single-fit-ener-se_e2_a-2',
#               '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/geo-opt-all/single-fit-ener-se_e2_a-2']
# model_spin = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/geo-opt-all/test-set/atom-3-all/single-fit-m-se_e2_a-2',
#               '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/geo-opt-all/single-fit-m-se_e2_a-2']
# model_ener = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/geo-opt-all/test-set/atom-3-all/single-fit-ener-se_e2_a-2-pt',
#               '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/geo-opt-all/single-fit-ener-se_e2_a-2-pt']
# model_spin = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/geo-opt-all/test-set/atom-3-all/single-fit-m-se_e2_a-2-pt',
#               '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/geo-opt-all/single-fit-m-se_e2_a-2-pt']
# model_ener = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/geo-opt-all/test-set/atom-3-all/single-fit-ener-se_e2_a-2-pt-limit_pref_e-1',
#               '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/geo-opt-all/single-fit-ener-se_e2_a-2-pt-limit_pref_e-1']
# model_spin = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/geo-opt-all/test-set/atom-3-all/single-fit-m-se_e2_a-2-pt',
#               '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/geo-opt-all/single-fit-m-se_e2_a-2-pt']
# model_ener = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/geo-opt-all/test-set/atom-3-all/single-fit-ener-se_e2_a-2-pt-rcut-20',
#               '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/geo-opt-all/single-fit-ener-se_e2_a-2-pt-rcut-20']
# model_spin = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/geo-opt-all/test-set/atom-3-all/single-fit-m-se_e2_a-2-pt-rcut-20',
#               '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/geo-opt-all/single-fit-m-se_e2_a-2-pt-rcut-20']
# model_ener = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/geo-opt-all/test-set/atom-3-all/single-fit-ener-dpa1-se_atten_v2',
#               '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/geo-opt-all/single-fit-ener-dpa1-se_atten_v2']
# model_spin = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/geo-opt-all/test-set/atom-3-all/single-fit-m-dpa1-se_atten_v2',
#               '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/geo-opt-all/single-fit-m-dpa1-se_atten_v2']
# model_ener = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/geo-opt-all/test-set/atom-3-all/single-fit-ener-dpa2',
#               '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/geo-opt-all/single-fit-ener-dpa2']
# model_spin = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/geo-opt-all/test-set/atom-3-all/single-fit-m-dpa2',
#               '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/geo-opt-all/single-fit-m-dpa2']
# model_ener = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/geo-opt-all/test-set/atom-3-all/multi-task-dpa2-limit_pref_e-1/',
#               '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/geo-opt-all/multi-task-dpa2-limit_pref_e-1/']
# model_spin = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/geo-opt-all/test-set/atom-3-all/single-fit-m-dpa2',
#               '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/geo-opt-all/single-fit-m-dpa2']
# database = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/geo-opt-all/test-set/atom-3-all/database_spin',
#             '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/geo-opt-all/database_spin']
# axis_lim_y = np.array([-4.05, -3.05])
# pos_array = np.array(([0.45, 0.1], [0.45, 0.2]))
# val = ['_', '_']

# Bulk hematite 0K NEB single fit 221 441
# model_ener = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/geo-opt-all/test-set/atom-3-all/single-fit-ener-se_e2_a-2',
#               '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/geo-opt-all/single-fit-ener-se_e2_a-2',
#               '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/geo-opt-all/test-set/441/single-fit-ener-se_e2_a-2']
# model_spin = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/geo-opt-all/test-set/atom-3-all/single-fit-m-se_e2_a-2',
#               '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/geo-opt-all/single-fit-m-se_e2_a-2',
#               '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/geo-opt-all/test-set/441/single-fit-m-se_e2_a-2']
# model_ener = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/geo-opt-all/test-set/atom-3-all/single-fit-ener-se_e2_a-2-pt',
#               '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/geo-opt-all/single-fit-ener-se_e2_a-2-pt',
#               '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/geo-opt-all/test-set/441/single-fit-ener-se_e2_a-2-pt']
# model_spin = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/geo-opt-all/test-set/atom-3-all/single-fit-m-se_e2_a-2-pt',
#               '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/geo-opt-all/single-fit-m-se_e2_a-2-pt',
#               '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/geo-opt-all/test-set/441/single-fit-m-se_e2_a-2-pt']
# model_ener = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/geo-opt-all/test-set/atom-3-all/single-fit-ener-dpa1-se_atten_v2',
#               '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/geo-opt-all/single-fit-ener-dpa1-se_atten_v2',
#               '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/geo-opt-all/test-set/441/single-fit-ener-dpa1-se_atten_v2']
# model_spin = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/geo-opt-all/test-set/atom-3-all/single-fit-m-dpa1-se_atten_v2',
#               '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/geo-opt-all/single-fit-m-dpa1-se_atten_v2',
#               '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/geo-opt-all/test-set/441/single-fit-m-dpa1-se_atten_v2']
# model_ener = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/geo-opt-all/test-set/atom-3-all/single-fit-ener-dpa2',
#               '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/geo-opt-all/single-fit-ener-dpa2',
#               '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/geo-opt-all/test-set/441/single-fit-ener-dpa2']
# model_spin = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/geo-opt-all/test-set/atom-3-all/single-fit-m-dpa2',
#               '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/geo-opt-all/single-fit-m-dpa2',
#               '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/geo-opt-all/test-set/441/single-fit-m-dpa2']
# database = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/geo-opt-all/test-set/atom-3-all/database_spin',
#             '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/geo-opt-all/database_spin',
#             '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/geo-opt-all/test-set/441/database_spin']
# axis_lim_y = np.array([-4.05, -3.05])
# pos_array = np.array(([0.38, 0.15], [0.38, 0.25], [0.38, 0.05]))

# Bulk hematite 0K NEB + DFT-MD
# model = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/400k-b-geo-opt-all/multi-task-dpa2'
# database = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/400k-b-geo-opt-all/database_spin'
# axis_lim_y = np.array([-4.05, -3.05])

# Bulk hematite 0K NEB 221 rcut 5 single fit
# model = ['single-fit-ener-se_e2_a-tf', 'single-fit-m-se_e2_a-tf']
# model = ['single-fit-ener-se_e2_a', 'single-fit-m-se_e2_a']
# model = ['single-fit-ener-se_e3', 'single-fit-m-se_e3']
# model = ['single-fit-ener-hybrid-se_e2_a-se_e3', 'single-fit-m-hybrid-se_e2_a-se_e3']
# model = ['single-fit-ener-dpa1-se_atten_v2-attn_layer-2', 'single-fit-m-dpa1-se_atten_v2-attn_layer-2']
# model = ['single-fit-ener-dpa1-se_atten_v2-attn_layer-3', 'single-fit-m-dpa1-se_atten_v2-attn_layer-3']
# model = ['single-fit-ener-dpa2-nlayers-6', 'single-fit-m-dpa2-nlayers-6']
# model_ener = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/geo-opt-all-rcut-5/{}'.format(model[0]),
#               '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/geo-opt-all-rcut-5/{}'.format(model[0])]
# model_spin = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/geo-opt-all-rcut-5/{}'.format(model[1]),
#               '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/geo-opt-all-rcut-5/{}'.format(model[1])]
# database = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/geo-opt-all-rcut-5/database_spin_test/1',
#             '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/geo-opt-all-rcut-5/database_spin_train']
# val = ['_1_', '_0_']
# text_array = ['0 K test', '0 K train', 'Test 441']
# axis_lim_y = np.array([-4.05, -3.05])
# pos_array_energy = np.array(([0.38, 0.05], [0.38, 0.15]))
# pos_array_force = np.array(([0.6, 0.05], [0.6, 0.15]))
# pos_array_spin = np.array(([0.78, 0.05], [0.78, 0.15]))

# Bulk hematite 0K NEB 221 rcut 5 multi fit
# model = 'multi-task-se_e2_a'
# model = 'multi-task-se_e3'
# model = 'multi-task-hybrid-se_e2_a-se_e3'
# model = 'multi-task-dpa1-se_atten_v2-attn_layer-1'
# model = 'multi-task-dpa1-se_atten_v2-attn_layer-2'
# model = 'multi-task-dpa1-se_atten_v2-attn_layer-3'
# model = 'multi-task-dpa2-nlayers-1'
# model = 'multi-task-dpa2-nlayers-2'
# model = 'multi-task-dpa2-nlayers-3'
# model = 'multi-task-dpa2-nlayers-3-update_g2_has_attn-f'
# model = 'multi-task-dpa2-nlayers-3-use_three_body-f'
# model = 'multi-task-dpa2-nlayers-6'
# model = 'multi-task-dpa2-nlayers-12'
# model_ener = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/geo-opt-all-rcut-5/{}'.format(model),
#               '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/geo-opt-all-rcut-5/{}'.format(model)]
# model_spin = model_ener
# database = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/geo-opt-all-rcut-5/database_spin_test/1',
#             '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/geo-opt-all-rcut-5/database_spin_train']
# val = ['_1_', '_0_']
# text_array = ['0 K test', '0 K train', 'Test 441']
# axis_lim_y = np.array([-4.05, -3.05])
# pos_array_energy = np.array(([0.38, 0.05], [0.38, 0.15]))
# pos_array_force = np.array(([0.6, 0.05], [0.6, 0.15]))
# pos_array_spin = np.array(([0.78, 0.05], [0.78, 0.15]))

# Bulk hematite MD 400k-f 5 hops training data with 400k-b 1 hop test data single-fit
# model = ['single-fit-ener-se_e2_a-tf', 'single-fit-m-se_e2_a-tf']
# model = ['single-fit-ener-se_e2_a', 'single-fit-m-se_e2_a']
# model = ['single-fit-ener-se_e3', 'single-fit-m-se_e3']
# model = ['single-fit-ener-hybrid-se_e2_a-se_e3', 'single-fit-m-hybrid-se_e2_a-se_e3']
# model = ['single-fit-ener-dpa1-se_atten_v2-attn_layer-2', 'single-fit-m-dpa1-se_atten_v2-attn_layer-2']
# model = ['single-fit-ener-dpa1-se_atten_v2-attn_layer-3', 'single-fit-m-dpa1-se_atten_v2-attn_layer-3']
# model = ['single-fit-ener-dpa2-nlayers-6', 'single-fit-m-dpa2-nlayers-6']
# model_ener = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/400k-b-rcut-5/{}'.format(model[0]),
#               '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/400k-b-rcut-5/{}'.format(model[0])]
# model_spin = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/400k-b-rcut-5/{}'.format(model[1]),
#               '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/400k-b-rcut-5/{}'.format(model[1])]
# database = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/400k-b-rcut-5/database_spin_test/1',
#             '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/400k-b-rcut-5/database_spin_train']
# val = ['_1_', '_0_']
# # axis_lim_y = np.array([-4.05, -3.05])
# # pos_array = np.array(([0.45, 0.1], [0.45, 0.2]))
# text_array = ['400 K test', '400 K train', 'Test 441']
# axis_lim_y = np.array([-4.05, -3.05])
# pos_array_energy = np.array(([0.38, 0.05], [0.38, 0.15]))
# pos_array_force = np.array(([0.6, 0.05], [0.6, 0.15]))
# pos_array_spin = np.array(([0.78, 0.05], [0.78, 0.15]))

# Bulk hematite MD 400k-f 5 hops training data with 400k-b 1 hop test data multi-fit
# model = 'multi-task-se_e2_a'
# model = 'multi-task-se_e3'
# model = 'multi-task-hybrid-se_e2_a-se_e3'
# model = 'multi-task-dpa1-se_atten_v2-attn_layer-0'
# model = 'multi-task-dpa1-se_atten_v2-attn_layer-1'
# model = 'multi-task-dpa1-se_atten_v2-attn_layer-2'
# model = 'multi-task-dpa1-se_atten_v2-attn_layer-3'
# model = 'multi-task-dpa2-nlayers-1'
# model = 'multi-task-dpa2-nlayers-2'
# model = 'multi-task-dpa2-nlayers-3'
# model = 'multi-task-dpa2-nlayers-3-update_g2_has_attn-f'
# model = 'multi-task-dpa2-nlayers-3-use_three_body-f'
# model = 'multi-task-dpa2-nlayers-6'
# model = 'multi-task-dpa2-nlayers-6-fe-o-only'
# model = 'multi-task-dpa2-nlayers-6-no-atomener'
# model = 'multi-task-dpa2-nlayers-6-rcut-6'
# model = 'multi-task-dpa2-nlayers-6-rcut-10'
# model = 'multi-task-dpa2-nlayers-6-rcut-20'
# model_ener = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/400k-b-rcut-5/{}'.format(model),
#               '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/400k-b-rcut-5/{}'.format(model)]
# model_spin = model_ener
# database = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/400k-b-rcut-5/database_spin_test/1',
#             '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/400k-b-rcut-5/database_spin_train']
# val = ['_1_', '_0_']
# text_array = ['400 K test', '400 K train', 'Test 441']
# axis_lim_y = np.array([-4.05, -3.05])
# pos_array_energy = np.array(([0.38, 0.05], [0.38, 0.15]))
# pos_array_force = np.array(([0.6, 0.05], [0.6, 0.15]))
# pos_array_spin = np.array(([0.78, 0.05], [0.78, 0.15]))

# Bulk MgO 222
# model = ['single-fit-ener-se_e2_a', 'single-fit-m-se_e2_a']
# # model = ['single-fit-ener-se_e3', 'single-fit-m-se_e3']
# # model = ['single-fit-ener-dpa1-se_atten_v2-attn_layer-2', 'single-fit-m-dpa1-se_atten_v2-attn_layer-2']
# # model = ['single-fit-ener-dpa1-se_atten_v2-attn_layer-3', 'single-fit-m-dpa1-se_atten_v2-attn_layer-3']
# # model = ['single-fit-ener-dpa2-nlayers-6', 'single-fit-m-dpa2-nlayers-6']
# # model = ['single-fit-ener-dpa3', 'single-fit-m-dpa3']
# # model = ['single-fit-ener-dpa3-prefactor-20-60', 'single-fit-m-dpa3-prefactor-20-60']
# # model = ['single-fit-ener-dpa3-prefactor-20-60-asel-48', 'single-fit-m-dpa3-prefactor-20-60-asel-48']
# # model = ['multi-task-se_e2_a', 'multi-task-se_e2_a']
# # model = ['multi-task-dpa1-se_atten_v2-attn_layer-2', 'multi-task-dpa1-se_atten_v2-attn_layer-2']
# # model = ['multi-task-dpa2-nlayers-6', 'multi-task-dpa2-nlayers-6']
# # model = ['multi-task-dpa3', 'multi-task-dpa3']
# # spin_is_population = False
# # folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/mgo/deepmd/cell-222/electron-u-6'
# # folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/mgo/deepmd/cell-222/electron-u-8-2'
# # folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/mgo/deepmd/cell-222/electron-u-8-rs'
# # database = ['{}/database_spin_train/'.format(folder),
# #             '{}/database_spin_test/1'.format(folder),
# #             '{}/database_spin_test/2'.format(folder)]
# spin_is_population = True
# folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/mgo/deepmd/cell-222/electron-u-6-population-fixed'
# model_ener = ['{}/{}'.format(folder, model[0])] * 4
# model_spin = ['{}/{}'.format(folder, model[1])] * 4
# database = ['{}/database_population_train/'.format(folder),
#             '{}/database_population_test/1'.format(folder),
#             '{}/database_population_test/2'.format(folder)]
# val = ['_0_', '_1_', '_2_']
# axis_lim_y = np.array([0, 1])
# pos_array_energy = np.array(([0.38, 0.25], [0.38, 0.15], [0.38, 0.05]))
# pos_array_force = np.array(([0.6, 0.25], [0.6, 0.15], [0.6, 0.05]))
# pos_array_spin = np.array(([0.78, 0.25], [0.78, 0.15], [0.78, 0.05]))
# text_array = ['400 K train', '400 K valid', '400 K test']
# num_atoms = 64

# Bulk TiO2 336
# model = ['single-fit-ener-se_e2_a', 'single-fit-m-se_e2_a']
model = ['single-fit-ener-dpa3', 'single-fit-m-dpa3']
spin_is_population = True
folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/deepmd/cell-336/hse-25-schwarz-e-1e-4-f-1e-6-cfit11-cpfit3'
model_ener = ['{}/{}'.format(folder, model[0])] * 4
model_spin = ['{}/{}'.format(folder, model[1])] * 4
database = ['{}/database_population_train/'.format(folder),
            '{}/database_population_test/1'.format(folder),
            '{}/database_population_test/2'.format(folder)]
val = ['_0_', '_1_', '_2_']
axis_lim_y = np.array([0, 1])
pos_array_energy = np.array(([0.38, 0.25], [0.38, 0.15], [0.38, 0.05]))
pos_array_force = np.array(([0.6, 0.25], [0.6, 0.15], [0.6, 0.05]))
pos_array_spin = np.array(([0.78, 0.25], [0.78, 0.15], [0.78, 0.05]))
text_array = ['400 K train', '400 K valid', '400 K test']
num_atoms = 324

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
spin_0 = []
ener_0 = []
force_0 = []
spin_1 = []
ener_1 = []
force_1 = []

test = np.load("{}/database_population_train/set.000/atomic_spin.npy".format(folder))
# test = np.load("{}/database_population_test/1/set.000/atomic_spin.npy".format(folder))
print('test', test.shape)
# print('test', test[0, :, :])

for i in range(len(database)):
    dft_e.append(np.load("{}/set.000/energy.npy".format(database[i], val[i])))
    dft_f.append(np.load("{}/set.000/force.npy".format(database[i], val[i])))
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

# 1. Plot parity
fig, axes = plt.subplots(2, 3, figsize=(15, 6))
for i in range(len(database)):
    plot_ener(dft_e[i], ener_0[i], axes[0, 0], color_plot=color_plot_array[i], pos=pos_array_energy[i], text=text_array[i], title="Energy, No-aparam")
    plot_ener(dft_e[i], ener_1[i], axes[1, 0], color_plot=color_plot_array[i], pos=pos_array_energy[i], text=text_array[i], title="Energy, Yes-aparam")
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

# 2. Plot parity
fig2, axes2 = plt.subplots(1, 3, figsize=(15, 3))
for i in range(len(database)):
    plot_ener(dft_e[i], ener_1[i], axes2[0], color_plot=color_plot_array[i], pos=pos_array_energy[i], text=text_array[i])
    plot_force(dft_f[i], force_1[i], axes2[1], color_plot=color_plot_array[i], pos=pos_array_force[i], text=text_array[i])
    if not spin_is_population:
        plot_spin(dft_s[i], spin_1[i], axes2[2], color_plot=color_plot_array[i], pos=pos_array_spin[i], text=text_array[i])
    if spin_is_population:
        plot_spin((dft_s[i][:, :, 0] - dft_s[i][:, :, 1]), (spin_1[i][:, :, 0] - spin_1[i][:, :, 1]), axes2[2], color_plot=color_plot_array[i], pos=pos_array_spin[i], text=text_array[i])
plt.tight_layout()
for i in range(len(database)):
    plt.savefig("{}/fit_1x3_folders_{}.png".format(model_spin[i], len(model_ener)), dpi=600)
    plt.savefig("{}/fit_1x3_folders_{}.png".format(model_ener[i], len(model_spin)), dpi=600)

# Plot spin moment training
plot_spin_time = True
# plot_spin_time = False
if plot_spin_time:
    axis_lim_y = [-0.02, 0.9]
    fig_spin_train, axes_spin_train = plt.subplots()
    plot_spin_time1((dft_s[0][:, :, 0] - dft_s[0][:, :, 1]), (spin_1[0][:, :, 0] - spin_1[0][:, :, 1]),
                    axes_spin_train, axis_lim_y, title="Energy, No-aparam")
    plt.tight_layout()
    plt.savefig("{}/spin_train.png".format(model_spin[0]), dpi=600)
    if zoom:
        plt.xlim(axis_lim_x_zoom[0], axis_lim_x_zoom[1])
    plt.savefig("{}/spin_train_zoom{}.png".format(model_spin[0], transition_time_plot), dpi=600)

# Plot spin moment validation
if plot_spin_time:
    axis_lim_y = [-0.02, 0.9]
    fig_spin_valid, axes_spin_valid = plt.subplots()
    plot_spin_time1((dft_s[1][:, :, 0] - dft_s[1][:, :, 1]), (spin_1[1][:, :, 0] - spin_1[1][:, :, 1]),
                    axes_spin_valid, axis_lim_y, title="Energy, No-aparam")
    plt.tight_layout()
    plt.savefig("{}/spin_valid.png".format(model_spin[0]), dpi=600)
    if zoom:
        plt.xlim(axis_lim_x_zoom[0], axis_lim_x_zoom[1])
    plt.savefig("{}/spin_valid_zoom{}.png".format(model_spin[0], transition_time_plot), dpi=600)

# Plot spin moment test
if plot_spin_time:
    axis_lim_y = [-0.02, 0.9]
    fig_spin_test, axes_spin_test = plt.subplots()
    plot_spin_time1((dft_s[2][:, :, 0] - dft_s[2][:, :, 1]), (spin_1[2][:, :, 0] - spin_1[2][:, :, 1]),
                    axes_spin_test, axis_lim_y, title="Energy, No-aparam")
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
