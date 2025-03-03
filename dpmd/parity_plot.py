import dpdata
import numpy as np
import matplotlib.pyplot as plt

# training_systems = dpdata.LabeledSystem(
#     "/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/shared/denan_li/chris/database", fmt="deepmd/npy"
# )

database = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/chain/hf/frozen-fe/single-task/database'
file = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/chain/hf/frozen-fe/single-task/single-fit-ener/graph.pb'

training_systems = dpdata.LabeledSystem(database, fmt="deepmd/npy")
predict = training_systems.predict(file)


plt.scatter(training_systems["energies"], predict["energies"])

x_range = np.linspace(plt.xlim()[0], plt.xlim()[1])

plt.plot(x_range, x_range, "r--", linewidth=0.25)
plt.xlabel("Energy of DFT")
plt.ylabel("Energy predicted by deep potential")
plt.plot()

# import numpy as np
# import matplotlib.pyplot as plt
# from sklearn.metrics import mean_absolute_error
#
# # 1. Load data
# model = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/shared/denan_li/multi-task'
# # model = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/shared/denan_li/testing/multi-task-copy'
# database = '{}/database_spin'.format(model)
#
# # model = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/shared/denan_li/delete/denan_li/multi-task'
# # database = '{}/database_spin'.format(model)
#
# # model = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/shared/denan_li/delete/denan_li/multi-task-copy'
# # database = '{}/database_spin'.format(model)
#
# # model = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/chain/hf/frozen-fe/multi-task'
# # database = '{}/database_spin'.format(model)
#
# # model = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/chain/hf/frozen-fe-h2o/multi-task'
# # database = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/chain/hf/frozen-fe-h2o/database'
#
# dft_e = np.load("{}/set.000/energy.npy".format(database))
# dft_f = np.load("{}/set.000/force.npy".format(database))
# dft_s = np.load("{}/set.000/atom_ener.npy".format(database))
#
# spin_0 = np.load("{}/0_spin.npy".format(model))
# ener_0 = np.load("{}/0_ener.npy".format(model))
# force_0 = np.load("{}/0_force.npy".format(model))
#
# spin_1 = np.load("{}/1_spin.npy".format(model))
# ener_1 = np.load("{}/1_ener.npy".format(model))
# force_1 = np.load("{}/1_force.npy".format(model))
#
#
# def plot_ener(dft, dp, ax, title=None):
#     ax.plot(dft.flatten(), dp.flatten(), 'o')
#     ax.set_xlabel("DFT energy (eV)")
#     ax.set_ylabel("DP energy (eV)")
#     mae = mean_absolute_error(dft.flatten(), dp.flatten()) / 77 * 1000 # unit: meV/atom
#     ax.text(0.5, 0.1, f"MAE: {mae:.2f} meV/atom", transform=ax.transAxes)
#     min_val, max_val = min(dft.min(), dp.min()), max(dft.max(), dp.max())
#     ax.plot([min_val, max_val], [min_val, max_val], 'k--')
#     if title is not None:
#         ax.set_title(title)
#     else:
#         ax.set_title("Energy")
#
# def plot_force(dft, dp, ax, title=None):
#     ax.plot(dft.flatten(), dp.flatten(), 'o')
#     ax.set_xlabel("DFT force (eV/$\AA$)")
#     ax.set_ylabel("DP force (eV/$\AA$)")
#     mae = mean_absolute_error(dft.flatten(), dp.flatten())
#     ax.text(0.5, 0.1, f"MAE: {mae:.2f} eV/$\AA$", transform=ax.transAxes)
#     min_val, max_val = min(dft.min(), dp.min()), max(dft.max(), dp.max())
#     ax.plot([min_val, max_val], [min_val, max_val], 'k--')
#     if title is not None:
#         ax.set_title(title)
#     else:
#         ax.set_title("Force")
#
# def plot_spin(dft, dp, ax, title=None):
#     ax.plot(dft.flatten(), dp.flatten(), 'o')
#     ax.set_xlabel("DFT spin")
#     ax.set_ylabel("DP spin")
#     min_val, max_val = min(dft.min(), dp.min()), max(dft.max(), dp.max())
#     mae = mean_absolute_error(dft.flatten(), dp.flatten())
#     ax.text(0.5, 0.1, f"MAE: {mae:.2f}", transform=ax.transAxes)
#     ax.plot([min_val, max_val], [min_val, max_val], 'k--')
#     if title is not None:
#         ax.set_title(title)
#     else:
#         ax.set_title("Spin")
#
# # 2. Plot
# fig, axes = plt.subplots(2, 3, figsize=(15, 6))
# plot_ener(dft_e, ener_0, axes[0, 0], title="Energy, No-aparam")
# plot_ener(dft_e, ener_1, axes[1, 0], title="Energy, Yes-aparam")
# plot_force(dft_f, force_0, axes[0, 1], title="Force, No-aparam")
# plot_force(dft_f, force_1, axes[1, 1], title="Force, Yes-aparam")
# plot_spin(dft_s, spin_0, axes[0, 2], title="Spin, No-aparam")
# plot_spin(dft_s, spin_1, axes[1, 2], title="Spin, Yes-aparam")
# plt.tight_layout()
# plt.savefig("fit.png", dpi=600)
#
if __name__ == "__main__":
    print('Finished.')
    plt.show()
