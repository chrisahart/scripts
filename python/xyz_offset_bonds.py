from __future__ import absolute_import, division, print_function
import pandas as pd
import numpy as np
import glob
from general import load_coordinates
from general import print_xyz
import random

def calc_distance(x1, y1, z1, x2, y2, z2):
    distance = np.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2 + (z1 - z2) ** 2)
    return distance

""" .xyz offset bonds. 
    Modify bond lengths around specific atoms. 
    Used to encourage excess charge localisation by offsetting certain bond lengths """

# Experimental structure
# folder_out = '/media/chris/Elements/Backup/Archer-2/surfin/hematite/geo_opt/electron-offset-b/frozen-water-h-24hr-cg/'
# folder_in = '/media/chris/DATA1/Storage/University/PhD/Programming/dft_ml_md/output/fe_bulk/hematite/441_supercell_cdft/structures/offset'
# folder_out = '/media/chris/DATA1/Storage/University/PhD/Programming/dft_ml_md/output/fe_bulk/hematite/441_supercell_cdft/structures/offset'
# folder_out = '/media/chris/DATA1/Storage/University/PhD/Programming/dft_ml_md/output/fe_bulk/hematite/221_supercell_cdft/structures/offset/electron/one_site'
# folder_in = '/media/chris/DATA1/Storage/University/PhD/Programming/dft_ml_md/output/fe_bulk/hematite/221_supercell_cdft/structures/offset/electron/one_site'
# folder_out = '/media/chris/DATA1/Storage/University/PhD/Programming/dft_ml_md/output/fe_bulk/hematite/331_supercell_cdft/structures/electron/offset'
# folder_out = '/media/chris/DATA1/Storage/University/PhD/Programming/dft_ml_md/output/fe_bulk/hematite/441_supercell_cdft/structures/electron/offset'
# folder_in = folder_out
# filename_in = 'input-old.xyz'
# filename_in = 'input.xyz'
# filename_out = 'input.xyz'

folder_out = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/chris_phd/bulk/paper-2-jacs/221_supercell_cdft/hole/dft/neutral/m600_neutral-cubes2'
# folder_out = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/archer/bulk/hole/221_supercell/reftraj/geo_opt/atom-1/atom-1-c-print'
folder_in = folder_out
filename_in = 'input.xyz'
filename_out = 'input-offset-1-c.xyz'
cols = ['Species', 'X', 'Y', 'Z']

coordinates, coord_x, coord_y, coord_z, species, num_atoms, num_timesteps = load_coordinates.load_values_coord(
    folder_in, filename_in, cols)

# Atom number and bond target (hole)
# labels = np.array([156, 166, 167, 168, 170, 169, 171]) - 1
# bond_target = np.array([1.85, 1.85, 1.95, 2.15, 1.98, 1.98])
# bond_target = np.array([1.86, 1.96, 1.85, 2.00, 2.15, 1.97])
# bond_target = np.ones(np.shape(labels)[0]-1) * 2.4  # Hole

# Atom number and bond target (electron)
# labels = np.array([29, 118, 44, 55, 49, 47, 107]) - 1
# labels = np.array([133, 118, 44, 82, 146, 154, 140]) - 1

# labels = np.array([36, 90, 46, 51, 54, 73, 43]) - 1
# labels = np.array([127, 90, 46, 153, 117, 142, 139]) - 1

# labels = np.array([35, 122, 93, 53, 50, 26, 25]) - 1
# labels = np.array([96, 122, 93, 116, 152, 91, 125]) - 1


# atom 0
# labels = np.array([13, 50, 88, 108, 67, 61, 95]) - 1
# neutral = np.array([1.94, 1.94, 1.94, 2.12, 2.12, 2.12])
# a = np.array([1.84505, 1.96217, 1.85736, 1.96894, 2.14771, 1.99975])
# b = np.array([1.85736, 1.84505, 1.96217, 1.99975, 1.96894, 2.14771])
# c = np.array([1.96217, 1.85736, 1.84505, 2.14771, 1.99975, 1.96894])

# atom 1
# labels = np.array([38, 108, 104, 106, 95, 86, 109]) - 1
labels = np.array([38, 108, 104, 106, 95, 63, 109]) - 1
neutral = np.array([2.12, 2.12, 2.12, 1.94, 94, 1.94])
a = np.array([2.00392, 2.14936, 1.96812, 1.84155, 1.95954, 1.84966])
b = np.array([1.96812, 2.00392, 2.14936, 1.95954, 1.84966, 1.84155])
c = np.array([2.14936, 1.96812, 2.00392, 1.84155, 1.95954, 1.84966])


# bond_target = np.array([2.3]*6)
bond_target = np.array([1.8]*6)
bond_target = c

Fe_O = np.zeros(np.shape(labels)[0] - 1)
Fe_O_x = np.zeros(np.shape(labels)[0] - 1)
Fe_O_y = np.zeros(np.shape(labels)[0] - 1)
Fe_O_z = np.zeros(np.shape(labels)[0] - 1)

# Calculate bond lengths
for i in range(np.shape(labels)[0] - 1):
    Fe_O[i] = calc_distance(coord_x[0, labels[0]], coord_y[0, labels[0]], coord_z[0, labels[0]],
                                      coord_x[0, labels[1+i]], coord_y[0, labels[1+i]], coord_z[0, labels[1+i]])
    Fe_O_x[i] = coord_x[0, labels[0]] - coord_x[0, labels[1 + i]]
    Fe_O_y[i] = coord_y[0, labels[0]] - coord_y[0, labels[1 + i]]
    Fe_O_z[i] = coord_z[0, labels[0]] - coord_z[0, labels[1 + i]]

# Calculate change to coord
bond_change = -1*np.sign(bond_target - Fe_O) * np.sqrt((abs(bond_target - Fe_O) ** 2)/3)
bond_changex = bond_change
bond_changey = bond_change
bond_changez = bond_change

# Modify coord
for i in range(len(labels) - 1):

    coord_x[0, labels[1 + i]] = coord_x[0, labels[1 + i]] + \
                                (bond_changex[i] * np.sign(coord_x[0, labels[0]] - coord_x[0, labels[1 + i]]))
    coord_y[0, labels[1 + i]] = coord_y[0, labels[1 + i]] + \
                                (bond_changey[i] * np.sign(coord_y[0, labels[0]] - coord_y[0, labels[1 + i]]))
    coord_z[0, labels[1 + i]] = coord_z[0, labels[1 + i]] + \
                                (bond_changez[i] * np.sign(coord_z[0, labels[0]] - coord_z[0, labels[1 + i]]))

# Randomly modify positions of all atoms
# change = 0.01
# for i in range(num_atoms):
#
#     coord_x[0, i] = coord_x[0, i] + change * functions.random_sign()
#     coord_y[0, i] = coord_y[0, i] + change * functions.random_sign()
#     coord_z[0, i] = coord_z[0, i] + change * functions.random_sign()

# Create pandas dataframe from species and coordinates
coord = np.column_stack((coord_x.ravel(), coord_y.ravel(), coord_z.ravel()))
coord_xyz = pd.DataFrame(data=coord)
coord_xyz.insert(loc=0, column='A', value=pd.Series(species).values)

# Print pandas dataframe to file
print_xyz.print_from_pandas(coord_xyz, num_atoms, '{}/{}'.format(folder_out, filename_out))