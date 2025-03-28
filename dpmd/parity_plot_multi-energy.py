import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import mean_absolute_error
from sklearn.metrics import root_mean_squared_error

# Bulk hematite DFT-MD single-fit neutral
model_ener = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/400k-neutral/single-fit-ener-se_e2_a-2-pt'
database = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/400k-neutral/database_ener_force'
axis_lim_y = np.array([-4.05, -3.05])


axis_lim_y = np.array([-4, -3])
transition_time = np.array([460, 117.5, 315, 656, 721, 1283])
zoom = False
transition_time_plot = 5
axis_lim_x_zoom = np.array([transition_time[transition_time_plot]-40, transition_time[transition_time_plot]+40])

dft_e = np.load("{}/set.000/energy.npy".format(database))
dft_f = np.load("{}/set.000/force.npy".format(database))
ener_0 = np.load("{}/0_ener.npy".format(model_ener))
force_0 = np.load("{}/0_force.npy".format(model_ener))


def plot_ener(dft, dp, ax, color_plot, pos, text, title=None):
    ax.plot(dft.flatten(), dp.flatten(), '.', color=color_plot)
    ax.set_xlabel("DFT energy (eV)")
    ax.set_ylabel("DP energy (eV)")
    mae = mean_absolute_error(dft.flatten(), dp.flatten()) / 77 * 1000  # unit: meV/atom
    rmse = root_mean_squared_error(dft.flatten(), dp.flatten()) / 77 * 1000  # unit: meV/atom
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

# 2. Plot 1x3
text_array = ['Test', 'Train']
pos_array = np.array(([0.5, 0.1], [0.5, 0.2]))
color_plot_array = ['r', 'g']
i = 1
fig2, axes2 = plt.subplots(1, 2, figsize=(10, 3))
plot_ener(dft_e, ener_0, axes2[0], color_plot=color_plot_array[i], pos=pos_array[i], text=text_array[i], title="Energy, Yes-aparam")
plot_force(dft_f, force_0, axes2[1], color_plot=color_plot_array[i], pos=pos_array[i], text=text_array[i], title="Force, Yes-aparam")
plt.tight_layout()
plt.savefig("{}/fit_1x3_folders_{}.png".format(model_ener, len(model_ener)), dpi=600)


if __name__ == "__main__":
    print('Finished.')
    plt.show()
