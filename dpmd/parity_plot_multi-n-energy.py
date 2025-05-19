import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import mean_absolute_error
from sklearn.metrics import root_mean_squared_error
from general import parameters as param
from scipy.stats import gaussian_kde

# Bulk hematite DFT-MD single-fit neutral
# model = ['single-fit-ener-se_e2_a-tf', 'single-fit-m-se_e2_a-tf']
# model = ['single-fit-ener-se_e2_a', 'single-fit-m-se_e2_a']
# model = ['single-fit-ener-se_e3', 'single-fit-m-se_e3']
# model = ['single-fit-ener-hybrid-se_e2_a-se_e3', 'single-fit-m-hybrid-se_e2_a-se_e3']
# model = ['single-fit-ener-dpa1-se_atten_v2-attn_layer-2', 'single-fit-m-dpa1-se_atten_v2-attn_layer-2']
# model = ['single-fit-ener-dpa1-se_atten_v2-attn_layer-3', 'single-fit-m-dpa1-se_atten_v2-attn_layer-3']
# model = ['single-fit-ener-dpa2-nlayers-6', 'single-fit-m-dpa2-nlayers-6']
# model = ['single-task-finetune-fe-o-only-rcut-6']
# model = ['single-task-finetune-fe-o-only']
# model = ['single-task-finetune-branch-H2O_H2O-PD-fe-o-only-use-pretrain-script']
# model = ['single-task-finetune-branch-H2O_H2O-PD-fe-o-only-rcut-6']
# model = ['single-task-finetune-branch-H2O_H2O-PD-fe-o-only']
# model = ['single-task-finetune-branch-Domains_Anode-fe-o-only-use-pretrain-script']
# model = ['single-task-finetune-branch-Domains_Anode-fe-o-only-rcut-6']
# model = ['single-task-finetune-branch-Domains_Anode-fe-o-only']
# model = ['multi-task-finetune-branch-Domains_Anode-fe-o-only']
# model = ['multi-task-finetune-branch-H2O_H2O-PD-fe-o-only']
# model_ener = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/400k-neutral-rcut-5/{}'.format(model[0]),
#               '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/400k-neutral-rcut-5/{}'.format(model[0])]
# database = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/400k-neutral-rcut-5/database_ener_force_train',
#             '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/400k-neutral-rcut-5/database_ener_force_test']
# model_ener = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/400k-neutral-rcut-5-new/{}'.format(model[0]),
#               '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/400k-neutral-rcut-5-new/{}'.format(model[0]),
#               '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/400k-neutral-rcut-5-new/{}'.format(model[0])]
# database = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/400k-neutral-rcut-5-new/database_ener_force_train',
#             '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/400k-neutral-rcut-5-new/database_ener_force_test/1',
#             '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/400k-neutral-rcut-5-new/database_ener_force_test/2']
# val = ['_0_', '_1_', '_2_']
# axis_lim_y = np.array([-4.05, -3.05])
# pos_array_energy = np.array(([0.38, 0.25], [0.38, 0.15], [0.38, 0.05]))
# pos_array_force = np.array(([0.6, 0.25], [0.6, 0.15], [0.6, 0.05]))
# pos_array_spin = np.array(([0.78, 0.25], [0.78, 0.15], [0.78, 0.05]))
# text_array = ['400 K train', '400 K valid', '400 K test']

# Bulk hematite DFT-MD single-fit neutral PBE
model = ['single-fit-ener-se_e2_a-tf', 'single-fit-m-se_e2_a-tf']
# model = ['single-fit-ener-se_e2_a', 'single-fit-m-se_e2_a']
# model = ['single-fit-ener-se_e3', 'single-fit-m-se_e3']
# model = ['single-fit-ener-hybrid-se_e2_a-se_e3', 'single-fit-m-hybrid-se_e2_a-se_e3']
# model = ['single-fit-ener-dpa1-se_atten_v2-attn_layer-2', 'single-fit-m-dpa1-se_atten_v2-attn_layer-2']
# model = ['single-fit-ener-dpa1-se_atten_v2-attn_layer-3', 'single-fit-m-dpa1-se_atten_v2-attn_layer-3']
# model = ['single-fit-ener-dpa2-nlayers-6', 'single-fit-m-dpa2-nlayers-6']
# model = ['single-fit-ener-dpa3-rcut-6', 'single-fit-ener-dpa3-rcut-6']
# model = ['single-fit-ener-dpa3-rcut-5', 'single-fit-ener-dpa3-rcut-5']
# model = ['single-fit-ener-dpa3-rcut-5-prefactor-20-60', 'single-fit-ener-dpa3-rcut-5-prefactor-20-60']
# model = ['single-fit-ener-dpa3-rcut-5-prefactor-20-60-asel-48', 'single-fit-ener-dpa3-rcut-5-prefactor-20-60-asel-48']
# model = ['single-task-finetune-fe-o-only-rcut-6']
# model = ['single-task-finetune-fe-o-only']
# model = ['single-task-finetune-branch-H2O_H2O-PD-fe-o-only-use-pretrain-script']
# model = ['single-task-finetune-branch-H2O_H2O-PD-fe-o-only-rcut-6']
# model = ['single-task-finetune-branch-H2O_H2O-PD-fe-o-only']
# model = ['single-task-finetune-branch-Domains_Anode-fe-o-only-use-pretrain-script']
# model = ['single-task-finetune-branch-Domains_Anode-fe-o-only-rcut-6']
# model = ['single-task-finetune-branch-Domains_Anode-fe-o-only']
# model = ['multi-task-finetune-branch-Domains_Anode-fe-o-only']
# model = ['multi-task-finetune-branch-H2O_H2O-PD-fe-o-only']
# model_ener = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liulab-hpc/bulk/400k-neutral-rcut-5-new/{}'.format(model[0]),
#               '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liulab-hpc/bulk/400k-neutral-rcut-5-new/{}'.format(model[0]),
#               '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liulab-hpc/bulk/400k-neutral-rcut-5-new/{}'.format(model[0])]
# database = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liulab-hpc/bulk/400k-neutral-rcut-5-new/database_ener_force_train',
#             '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liulab-hpc/bulk/400k-neutral-rcut-5-new/database_ener_force_test/1',
#             '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liulab-hpc/bulk/400k-neutral-rcut-5-new/database_ener_force_test/2']
# model_ener = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liulab-hpc/bulk/400k-neutral-rcut-5-new2/{}'.format(model[0]),
#               '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liulab-hpc/bulk/400k-neutral-rcut-5-new2/{}'.format(model[0]),
#               '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liulab-hpc/bulk/400k-neutral-rcut-5-new2/{}'.format(model[0])]
# database = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liulab-hpc/bulk/400k-neutral-rcut-5-new2/database_ener_force_train',
#             '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liulab-hpc/bulk/400k-neutral-rcut-5-new2/database_ener_force_test/1',
#             '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liulab-hpc/bulk/400k-neutral-rcut-5-new2/database_ener_force_test/2']
# model_ener = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liulab-hpc/bulk/400k-neutral-rcut-5-pbe-frames-2000/{}'.format(model[0]),
#               '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liulab-hpc/bulk/400k-neutral-rcut-5-pbe-frames-2000/{}'.format(model[0]),
#               '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liulab-hpc/bulk/400k-neutral-rcut-5-pbe-frames-2000/{}'.format(model[0])]
# database = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liulab-hpc/bulk/400k-neutral-rcut-5-pbe-frames-2000/database_ener_force_train',
#             '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liulab-hpc/bulk/400k-neutral-rcut-5-pbe-frames-2000/database_ener_force_test/1',
#             '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liulab-hpc/bulk/400k-neutral-rcut-5-pbe-frames-2000/database_ener_force_test/2']
# val = ['_0_', '_1_', '_2_']
# axis_lim_y = np.array([-4.05, -3.05])
# pos_array_energy = np.array(([0.38, 0.25], [0.38, 0.15], [0.38, 0.05]))
# pos_array_force = np.array(([0.6, 0.25], [0.6, 0.15], [0.6, 0.05]))
# pos_array_spin = np.array(([0.78, 0.25], [0.78, 0.15], [0.78, 0.05]))
# text_array = ['400 K train', '400 K valid', '400 K test']

# Bulk MgO 222 neutral
model = ['single-fit-ener-se_e2_a', 'single-fit-m-se_e2_a']
model = ['single-fit-ener-se_e3', 'single-fit-m-se_e3']
model = ['single-fit-ener-hybrid-se_e2_a-se_e3', 'single-fit-m-hybrid-se_e2_a-se_e3']
model = ['single-fit-ener-dpa1-se_atten_v2-attn_layer-2', 'single-fit-m-dpa1-se_atten_v2-attn_layer-2']
model = ['single-fit-ener-dpa1-se_atten_v2-attn_layer-3', 'single-fit-m-dpa1-se_atten_v2-attn_layer-3']
model = ['single-fit-ener-dpa2-nlayers-6', 'single-fit-m-dpa2-nlayers-6']
model = ['single-fit-ener-dpa3', 'single-fit-ener-dpa3']
model = ['single-fit-ener-dpa3-prefactor-20-60', 'single-fit-ener-dpa3-prefactor-20-60']
model = ['single-fit-ener-dpa3-prefactor-20-60-asel-48', 'single-fit-ener-dpa3-prefactor-20-60-asel-48']
model_ener = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/mgo/deepmd/cell-222/neutral/{}'.format(model[0]),
              '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/mgo/deepmd/cell-222/neutral/{}'.format(model[0]),
              '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/mgo/deepmd/cell-222/neutral/{}'.format(model[0])]
database = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/mgo/deepmd/cell-222/neutral/database_ener_force_train',
            '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/mgo/deepmd/cell-222/neutral/database_ener_force_test/1',
            '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/mgo/deepmd/cell-222/neutral/database_ener_force_test/2']
val = ['_0_', '_1_', '_2_']
axis_lim_y = np.array([-4.05, -3.05])
pos_array_energy = np.array(([0.38, 0.25], [0.38, 0.15], [0.38, 0.05]))
pos_array_force = np.array(([0.6, 0.25], [0.6, 0.15], [0.6, 0.05]))
pos_array_spin = np.array(([0.78, 0.25], [0.78, 0.15], [0.78, 0.05]))
text_array = ['400 K train', '400 K valid', '400 K test']

color_plot_array = ['r', 'g', 'b', 'm']
colorbar_plot_array = ['Reds', 'Greens', 'Blues', 'Oranges']
print('model', model)
axis_lim_y = np.array([-4, -3])
transition_time = np.array([460, 117.5, 315, 656, 721, 1283])
zoom = False
transition_time_plot = 5
axis_lim_x_zoom = np.array([transition_time[transition_time_plot]-40, transition_time[transition_time_plot]+40])

dft_e = []
dft_f = []
ener_0 = []
force_0 = []

for i in range(len(database)):
    dft_e.append(np.load("{}/set.000/energy.npy".format(database[i], val[i])))
    dft_f.append(np.load("{}/set.000/force.npy".format(database[i], val[i])))
    ener_0.append(np.load("{}/0{}ener.npy".format(model_ener[i], val[i])))
    force_0.append(np.load("{}/0{}force.npy".format(model_ener[i], val[i])))
    # ener_0.append(np.load("{}/1{}ener.npy".format(model_ener[i], val[i])))
    # force_0.append(np.load("{}/1{}force.npy".format(model_ener[i], val[i])))


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
    # else:
    #     ax.set_title("Energy")


def plot_force(dft, dp, ax, color_plot, pos, text, title=None):
    # dft = dft * param.hartree_per_bohr_to_ev_per_angstrom
    # dp = dp * param.hartree_per_bohr_to_ev_per_angstrom
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


def plot_force_with_density(dft, dp, ax, pos, text, title=None):
    dft = dft.flatten()
    dp = dp.flatten()

    mae = mean_absolute_error(dft, dp) * 1000  # Convert to meV/Å
    cmap = 'viridis'
    rmse = root_mean_squared_error(dft, dp) * 1000  # unit: meV/A

    print('plot_force_with_density')
    print('mean_absolute_error:', mae)
    print('root_mean_squared_error:', rmse)

    # Calculate density of points using a 2D histogram
    bins = 100

    counts, xedges, yedges = np.histogram2d(dft, dp, bins=bins)

    # Map the density back to the original points for scatter plot

    x_bin_indices = np.digitize(dft, xedges) - 1

    y_bin_indices = np.digitize(dp, yedges) - 1

    # Ensure indices are within bounds

    x_bin_indices = np.clip(x_bin_indices, 0, counts.shape[0] - 1)
    y_bin_indices = np.clip(y_bin_indices, 0, counts.shape[1] - 1)

    density = counts[x_bin_indices, y_bin_indices]

    # Create scatter plot with density-based coloring

    sc = ax.scatter(dft, dp, c=density, cmap=cmap, s=5, alpha=0.8)
    ax.set_xlabel("DFT force (eV/Å)")
    ax.set_ylabel("DP force (eV/Å)")

    # Add colorbar for density

    cbar = plt.colorbar(sc, ax=ax)
    cbar.ax.set_ylabel("Density")

    # Add RMSE annotation

    ax.text(pos[0], pos[1], f"{text} RMSE: {rmse:.2f} meV/Å", transform=ax.transAxes, color='black')

    # Add reference line y = x

    min_val, max_val = min(dft.min(), dp.min()), max(dft.max(), dp.max())
    ax.plot([min_val, max_val], [min_val, max_val], 'k--', linewidth=1)

    # Set title

    if title is not None:
        ax.set_title(title)
    else:
        ax.set_title("Force with Density")


def plot_multiple_force_with_density(datasets, ax, cmap='viridis', pos=(0.1, 0.9), text="", title=None, bins=100):
    """
    Plots multiple DFT vs DP datasets with a density-based shared colorbar.

    Parameters:
    datasets (list of dict): Each dict should have the structure:
                             {'dft': np.ndarray, 'dp': np.ndarray, 'label': str}
                             where 'label' is optional.
    ax (matplotlib.axes._axes.Axes): Axis object to plot on

    cmap (str): Colormap for the density of points

    pos (tuple): Position within the axis for RMSE text

    text (str): Prefix text for RMSE annotation

    title (str, optional): Title of the plot

    bins (int): Number of bins for the 2D histogram

    """

    # Combine all datasets into a single array for global density calculation

    all_dft = np.concatenate([data['dft'].flatten() for data in datasets])
    all_dp = np.concatenate([data['dp'].flatten() for data in datasets])

    # Compute a global 2D histogram

    counts, xedges, yedges = np.histogram2d(all_dft, all_dp, bins=bins)

    # Map density back to each dataset's points
    # colorbar_plot_array = ['Reds', 'Greens', 'Blues', 'Oranges']
    colorbar_plot_array = ['spring', 'summer', 'autumn', 'winter']

    densities = []
    for data in datasets:
        dft = data['dft'].flatten()
        dp = data['dp'].flatten()
        x_bin_indices = np.digitize(dft, xedges) - 1

        y_bin_indices = np.digitize(dp, yedges) - 1

        # Ensure indices are within bounds

        x_bin_indices = np.clip(x_bin_indices, 0, counts.shape[0] - 1)
        y_bin_indices = np.clip(y_bin_indices, 0, counts.shape[1] - 1)

        density = counts[x_bin_indices, y_bin_indices]
        densities.append(density)

    # Plot each dataset

    for i, data in enumerate(datasets):
        dft = data['dft'].flatten()
        dp = data['dp'].flatten()
        density = densities[i]
        label = data.get('label', f"Dataset {i+1}")
        # sc = ax.scatter(dft, dp, c=density, cmap=colorbar_plot_array[i], s=5, alpha=0.8, label=label)
        sc = ax.scatter(dft, dp, c=density, cmap='viridis', s=5, alpha=0.8, label=label)

    # Set axis labels

    ax.set_xlabel("DFT force (eV/Å)")
    ax.set_ylabel("DP force (eV/Å)")

    # Add colorbar for density (shared across all datasets)
    cbar = plt.colorbar(sc, ax=ax)
    cbar.ax.set_ylabel("Density")

    # Add reference line y = x

    min_val = min(all_dft.min(), all_dp.min())
    max_val = max(all_dft.max(), all_dp.max())
    ax.plot([min_val, max_val], [min_val, max_val], 'k--', linewidth=1)

    # Add RMSE annotations for each dataset

    for i, data in enumerate(datasets):
        dft = data['dft'].flatten()
        dp = data['dp'].flatten()
        mae = mean_absolute_error(dft, dp) * 1000  # Convert to meV/Å

        rmse = root_mean_squared_error(dft, dp) * 1000  # unit: meV/A

        label = data.get('label', f"Dataset {i+1}")
        ax.text(pos[0], pos[1] - i * 0.05, f"{label}: RMSE = {rmse:.2f} meV/Å",
                transform=ax.transAxes, color='black')

    # Set title

    if title is not None:
        ax.set_title(title)
    # else:
    #     ax.set_title("Force Comparison with Density")

    # Add legend

    # ax.legend()


# 2. Plot parity
fig2, axes2 = plt.subplots(1, 2, figsize=(10, 3))
for i in range(len(database)):
    plot_ener(dft_e[i], ener_0[i], axes2[0], color_plot=color_plot_array[i], pos=pos_array_energy[i], text=text_array[i])
    plot_force(dft_f[i], force_0[i], axes2[1], color_plot=color_plot_array[i], pos=pos_array_force[i], text=text_array[i])
    # plot_force_with_density(dft_f[i], force_0[i], axes2[1], pos=pos_array_energy[i], text=text_array[i])

# datasets = [
#     {'dft': dft_f[0], 'dp': force_0[0], 'label': "Dataset 1"},
#     {'dft': dft_f[1], 'dp': force_0[1], 'label': "Dataset 2"},
#     {'dft': dft_f[2], 'dp': force_0[2], 'label': "Dataset 3"}
# ]
# plot_multiple_force_with_density(datasets, axes2[1], cmap='plasma', pos=(0.1, 0.9), text="Test", bins=200)

plt.tight_layout()
for i in range(len(database)):
    plt.savefig("{}/fit_1x3_folders_{}.png".format(model_ener[i], len(model_ener)), dpi=600)

if __name__ == "__main__":
    print('Finished.')
    plt.show()
