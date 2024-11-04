import numpy as np
from matplotlib import pyplot as plt
from ase.io.cube import read_cube_data
from general import parameters as param

""" Plotting of CP2K .cube files for Hartree potential and charge density """

params = {'axes.formatter.limits': [-4, 4],
          'axes.labelsize': 'x-large',
          'axes.titlesize': 'x-large',
          'legend.fontsize': 'large',
          'lines.markersize': '8',
          }
plt.rcParams.update(params)
plot_color = ['k', 'b', 'r', 'g', 'm', 'grey', 'orange', 'y', 'brown', 'cyan', 'pink']

# cp2k-smeagol-examples/examples/melamine
folder1 = [
    '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/melamine/hx1/cu_test/basis-test/eps-default-1e-20/cu-sz',
    '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/melamine/hx1/cu_test/basis-test/eps-default-1e-20/cu-dzvp',
    '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/melamine/hx1/cu_test/basis-test/eps-default-1e-20/au-sz',
    '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/melamine/hx1/cu_test/basis-test/eps-default-1e-20/au-dzvp'
    ]
file1 = ['2_dft_wfn-ELECTRON_DENSITY-1_0.cube'] * 5
labels = ['cu-sz', 'cu-dzvp', 'au-sz', 'au-dzvp']

# Read .cube using ASE
data_charge = []
for i in range(len(folder1)):
    print('reading .cube files', folder1[i])
    data_charge.append(read_cube_data('{}/{}'.format(folder1[i], file1[i])))
print('Finished reading .cube files')

# Calculate average along axis
axis = 2
average_charge = []
for i in range(len(folder1)):
    average_charge.append(np.zeros(data_charge[i][0].shape[axis]))
    for j in range(data_charge[i][0].shape[axis]):
        average_charge[i][j] = np.mean(data_charge[i][0][:, :, j])
        # average_charge[i][j] = np.max(data_charge[i][0][:, :, j])
        # if np.mean(data_charge[i][0][:, :, j]) < 1e-20:
        #     average_charge[i][j] = np.NaN
    print('calculating average of', folder1[i])
    # print(np.sum(average_charge[i]))
    # average_charge[i] = average_charge[i] * 1/np.sum(average_charge[i])
    # print(np.sum(average_charge[i]))
    # print(np.min(average_charge[i]))
    # print(average_charge[i])
# print(average_charge)

# Plot Hartree and charge .cube difference
print(np.min(average_charge[0]))
print(average_charge[0])
print(data_charge[0][0].shape)
fig_charge, ax_charge = plt.subplots()
energy_grid = np.linspace(start=-10, stop=10, num=data_charge[0][0].shape[0])
for i in range(len(folder1)):
    ax_charge.plot(energy_grid, average_charge[i], 'k-', color=plot_color[i], label=labels[i])
    # ax_charge.plot(energy_grid, data_charge[i][0][int(297/2), int(297/2), :], 'k-', color=plot_color[i], label=labels[i])
ax_charge.set_yscale('log')
ax_charge.legend(frameon=False)
ax_charge.set_xlabel(r'Position z / Å')
ax_charge.set_ylabel('Charge density z')
fig_charge.tight_layout()
for i in range(len(folder1)):
    fig_charge.savefig('{}/charge_z.png'.format(folder1[i]), dpi=300)


if __name__ == "__main__":
    print(folder1)
    print('Finished.')
    plt.show()
