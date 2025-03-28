import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import mean_absolute_error
from sklearn.metrics import root_mean_squared_error

# Bulk hematite DFT-MD single-fit neutral
model_ener = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/400k-neutral/test-set/single-fit-ener-se_e2_a-2-pt-train',
              '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/400k-neutral/single-fit-ener-se_e2_a-2-pt-train']
database = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/400k-neutral/test-set/database_ener_force_valid',
            '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/400k-neutral/database_ener_force_train']
axis_lim_y = np.array([-4.05, -3.05])

pos_array = np.array(([0.38, 0.15], [0.38, 0.25], [0.38, 0.05]))
color_plot_array = ['r', 'g', 'b']
num_atoms = 77

text_array = ['Validation', 'Train', 'Test 441']
# text_array = ['Test', 'Train', 'Test 441']
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
    dft_e.append(np.load("{}/set.000/energy.npy".format(database[i])))
    dft_f.append(np.load("{}/set.000/force.npy".format(database[i])))
    ener_0.append(np.load("{}/0_ener.npy".format(model_ener[i])))
    force_0.append(np.load("{}/0_force.npy".format(model_ener[i])))
    print(len(dft_e[i]))
    print(len(ener_0[i]))


def plot_ener(dft, dp, ax, color_plot, pos, text, title=None):
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
    ax.text(pos[0], pos[1], f"{text} RMSE: {rmse:.2f} meV/atom", transform=ax.transAxes)
    if title is not None:
        ax.set_title(title)
    # else:
    #     ax.set_title("Energy")


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
    # else:
    #     ax.set_title("Force")


# 2. Plot parity
fig2, axes2 = plt.subplots(1, 2, figsize=(10, 3))
for i in range(len(database)):
    plot_ener(dft_e[i], ener_0[i], axes2[0], color_plot=color_plot_array[i], pos=pos_array[i], text=text_array[i])
    plot_force(dft_f[i], force_0[i], axes2[1], color_plot=color_plot_array[i], pos=pos_array[i], text=text_array[i])
plt.tight_layout()
for i in range(len(database)):
    plt.savefig("{}/fit_1x3_folders_{}.png".format(model_ener[i], len(model_ener)), dpi=600)

if __name__ == "__main__":
    print('Finished.')
    plt.show()
