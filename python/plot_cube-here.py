import numpy as np
from matplotlib import pyplot as plt
from ase.io.cube import read_cube_data
import optparse
import os

# Find file
parser = optparse.OptionParser()
opts, args = parser.parse_args()
folder = os.getcwd()
filename = args[0]
file = '{}/{}'.format(folder, filename)
print('File to be plotted:', file)

# Read .cube using ASE
data_1, atoms_1 = read_cube_data(file)

# Take z average of .cube
z_average_1 = np.zeros(data_1.shape[2])
for i in range(data_1.shape[2]):
    z_average_1[i] = np.mean(data_1[:, :, i])

# Set up plotting
energy_grid_1 = np.linspace(start=0, stop=atoms_1.get_cell()[2][2], num=data_1.shape[2])
xlim = [0-1, atoms_1.get_cell()[2][2]+1]

# Plot charge and hartree .cube
figcube_both, axcube_both = plt.subplots()
axcube_both.plot(energy_grid_1, z_average_1, 'k-')
axcube_both.set_xlim([xlim[0], xlim[1]])
axcube_both.set_xlabel(r'Position / Å')
# axcube_both.set_ylabel('Value')
figcube_both.tight_layout()
figcube_both.savefig('{}/plot_{}.png'.format(folder, filename), dpi=300)
    
if __name__ == "__main__":
    print('Finished.')
    plt.show()
