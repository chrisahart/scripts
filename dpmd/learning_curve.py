import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/chain/frozen-fe-h2o/single-task/single-fit-ener'
# folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/chain/hf/frozen-fe/single-task/single-fit-ener'
# legends = ["rmse_trn", "rmse_e_trn", "rmse_f_trn"]
# folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/chain/frozen-fe-h2o/single-task/single-fit-m'
# folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/chain/hf/frozen-fe/single-task/single-fit-m'
# legends = ["rmse_trn", "rmse_ae_trn"]

# folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/shared/denan_li/multi-task'
# folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/shared/denan_li/testing/multi-task-copy'

# HF with frozen Fe, H2O flies away
folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/chain/hf/frozen-fe/multi-task'

# HF with frozen Fe and H2O
folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/chain/hf/frozen-fe-h2o/multi-task'

# HSE06(35%) with frozen Fe and H2O, single polaron hop
# folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/chain/hse/frozen-fe-h2o/multi-task'

# Bulk hematite
# folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/bulk/hole/hops-0-species/multi-task'

legends = ["rmse_trn_ener_force", "rmse_e_trn_ener_force", "rmse_f_trn_ener_force", "rmse_trn_spin", "rmse_ae_trn_spin"]
set_ylim = [1e-6, 1e2]
set_xlim = [1e0, 1e6]

fig_lc, ax_lc = plt.subplots()
with open("{}/lcurve.out".format(folder)) as f:
    headers = f.readline().split()[1:]
lcurve = pd.DataFrame(
    np.loadtxt("{}/lcurve.out".format(folder)), columns=headers
)
for legend in legends:
    ax_lc.loglog(lcurve["step"], lcurve[legend], label=legend)
ax_lc.legend()
ax_lc.set_xlabel("Training steps")
ax_lc.set_ylabel("Loss")
ax_lc.set_xlim([set_xlim[0], set_xlim[1]])
ax_lc.set_ylim([set_ylim[0], set_ylim[1]])
fig_lc.tight_layout()

# raw = np.loadtxt("{}/lcurve.out".format(folder), skiprows=2)
# step, ener_force_rmse, ener_force_e_rmse, ener_force_f_rmse, spin_ae_rmse, learning_rate = raw[:, 0], raw[:, 1], raw[:, 2], raw[:, 3], raw[:, 5], raw[:, 6]
# legend = ["Network 1", "Network 1 Energy", "Network 1 Force", "Network 2 Spin"]
# data = [ener_force_rmse, ener_force_e_rmse, ener_force_f_rmse, spin_ae_rmse]
# fig, ax = plt.subplots(figsize=(4, 3))
# for i in range(4):
#     plt.loglog(step, data[i], label=legend[i])
# plt.legend()
# plt.grid()
# plt.tight_layout()
# plt.xlabel("Step")
# plt.ylabel("Loss")
# plt.savefig("loss.png", dpi=600)

if __name__ == "__main__":
    print('Finished.')
    plt.show()
