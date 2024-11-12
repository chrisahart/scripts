import numpy as np
from matplotlib import pyplot as plt
from ase.io.cube import read_cube_data
from general import parameters as param

# CP2K
folder_save = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/other/basis_function_vacuum'
files = [
    '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/other/basis_function_vacuum/hf_EPS_RHO-1e-20/basis-ELECTRON_DENSITY-1_0.cube']
label = ['Hf']

# Read .cube using ASE
file = files[0]
data_1, atoms_1 = read_cube_data(file)

# Take z average of .cube
z_average_1 = np.zeros(data_1.shape[2])
for i in range(data_1.shape[2]):
    z_average_1[i] = np.mean(data_1[:, :, i])

# Set up plotting
energy_grid_1 = np.linspace(start=0, stop=atoms_1.get_cell()[2][2], num=data_1.shape[2])
cube_end = atoms_1.get_cell()[2][2]
xlim = [0, cube_end/2]

# Plot charge and hartree .cube
figcube_both, axcube_both = plt.subplots()
axcube_both.plot(energy_grid_1-cube_end/2, z_average_1, 'k-')
axcube_both.hlines(1e-10, 0, 100, 'r', alpha=0.5)
# axcube_both.hlines(1e-12, 0, 100, 'r', alpha=0.5)
axcube_both.set_xlim([xlim[0], xlim[1]])
axcube_both.set_yscale('log')
axcube_both.set_xlabel(r'Position / Å')
# axcube_both.set_ylabel('Value')
figcube_both.tight_layout()
# figcube_both.savefig('{}/plot_{}.png'.format(folder, filename), dpi=300)

# Read .cube using ASE
# cube = []
# z_average = []
# for i in range(len(files)):
#     cube.append(read_cube_data(files))
#     # z_average.append(np.zeros(cube[i][].shape[2]))

# print(cube)

if __name__ == "__main__":
    print('Finished.')
    plt.show()
