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

# folder_out = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/archer/bulk/221_supercell/md/hole/cleaned/400k-neutral'
# folder_in = folder_out
# filename_in = 'topology-label-all.xyz'
# filename_out = 'topology-label-all2.xyz'
# cols = ['Species', 'X', 'Y', 'Z']
# mapping = {
#     'Fe_A': 0,
#     'Fe_B': 1,
#     'Fe_C': 2,
#     'Fe_D': 3,
#     'Fe_E': 4,
#     'Fe_F': 5,
#     'O': 6
# }

folder_out = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/structures'
folder_in = folder_out
filename_in = 'rutile_335.xyz'
filename_out = 'rutile_335-offset.xyz'

folder_out = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-336/md/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-csvr-timecon-1-COMVEL_TO-1e-10-nvt-u-3.5-electron-dz2-u-3.0-rs-print'
folder_in = folder_out
filename_in = 'last.xyz'
filename_out = 'last_offset.xyz'

cols = ['Species', 'X', 'Y', 'Z']
mapping = {
    'Ti': 0,
    'O': 1
}

# folder_out = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/mgo/structures'
# folder_in = folder_out
# filename_in = 'rocksalt_222.xyz'
# filename_out = 'rocksalt_222-offset.xyz'
# # filename_in = 'rocksalt_333.xyz'
# # filename_out = 'rocksalt_333-offset.xyz'
# cols = ['Species', 'X', 'Y', 'Z']
# mapping = {
#     'Mg': 0,
#     'O': 1
# }

coordinates, coord_x, coord_y, coord_z, species, num_atoms, num_timesteps = load_coordinates.load_values_coord(
    folder_in, filename_in, cols)


species_print = list(species)
my_string = "\n".join(species_print)
with open('{}/species.txt'.format(folder_out), 'w') as file:
    file.write(my_string)


converted_list = [mapping[item] for item in species_print]
converted_string = "\n".join(map(str, converted_list))
with open('{}/species_map.txt'.format(folder_out), 'w') as file:
    file.write(converted_string)


# atom 0
# labels = np.array([13, 50, 88, 108, 67, 61, 95]) - 1
# neutral = np.array([1.94, 1.94, 1.94, 2.12, 2.12, 2.12])
# a = np.array([1.84505, 1.96217, 1.85736, 1.96894, 2.14771, 1.99975])
# b = np.array([1.85736, 1.84505, 1.96217, 1.99975, 1.96894, 2.14771])
# c = np.array([1.96217, 1.85736, 1.84505, 2.14771, 1.99975, 1.96894])

# atom 1
# labels = np.array([38, 108, 104, 106, 95, 86, 109]) - 1
# labels = np.array([38, 108, 104, 106, 95, 63, 109]) - 1
# neutral = np.array([2.12, 2.12, 2.12, 1.94, 94, 1.94])
# a = np.array([2.00392, 2.14936, 1.96812, 1.84155, 1.95954, 1.84966])
# b = np.array([1.96812, 2.00392, 2.14936, 1.95954, 1.84966, 1.84155])
# c = np.array([2.14936, 1.96812, 2.00392, 1.84155, 1.95954, 1.84966])

# atom 1
# labels = np.array([38, 108, 104, 106, 95, 86, 109]) - 1
labels = np.array([6, 54, 58, 88, 71, 67, 69]) - 1
neutral = np.array([2.12, 2.12, 2.12, 1.94, 94, 1.94])
a = np.array([2.00392, 2.14936, 1.96812, 1.84155, 1.95954, 1.84966])
b = np.array([1.96812, 2.00392, 2.14936, 1.95954, 1.84966, 1.84155])
c = np.array([2.14936, 1.96812, 2.00392, 1.84155, 1.95954, 1.84966])
# bond_target = np.array([2.3]*6)
# bond_target = np.array([1.8]*6)
bond_target = c

# atom tio2 334 supercell
# labels = np.array([69, 211, 197, 193, 216, 198, 199]) - 1
# neutral = np.array([1.95181, 1.97966, 1.95181, 1.95180, 1.97966, 1.95180])
# bond_target = np.array([2.2, 2.2, 2.2, 2.2, 2.2, 2.2])

# atom tio2 335 supercell
# labels = np.array([83, 243, 239, 247, 261, 269, 245]) - 1
# neutral = np.array([1.95181, 1.97966, 1.95181, 1.95180, 1.97966, 1.95180])
# bond_target = np.array([2.2, 2.2, 2.2, 2.2, 2.2, 2.2])

# atom tio2 336 supercell
labels = np.array([105, 297, 293, 105, 299, 319, 298]) - 1
neutral = np.array([1.95181, 1.97966, 1.95181, 1.95180, 1.97966, 1.95180])
bond_target = np.array([2.4, 2.4, 2.4, 2.4, 2.4, 2.4])

# atom tio2 445 supercell
# labels = np.array([159, 455, 459, 454, 456, 397, 453]) - 1
# neutral = np.array([1.95181, 1.97966, 1.95181, 1.95180, 1.97966, 1.95180])
# bond_target = np.array([2.2, 2.2, 2.2, 2.2, 2.2, 2.2])

# atom tio2 446 supercell
# labels = np.array([189, 543, 547, 550, 541, 473, 549]) - 1
# neutral = np.array([1.95181, 1.97966, 1.95181, 1.95180, 1.97966, 1.95180])
# bond_target = np.array([2.2, 2.2, 2.2, 2.2, 2.2, 2.2])

# atom mgo 222 supercell
# labels = np.array([39, 8, 15, 19, 23, 13, 7]) - 1
# neutral = np.array([2.09700, 2.09700, 2.09700, 2.09700, 2.09700, 2.09700])
# bond_target = np.array([1.8, 1.8, 1.8, 1.8, 1.8, 1.8])

# atom mgo 333 supercell
# labels = np.array([150, 41, 18, 96, 87, 15, 42]) - 1
# neutral = np.array([2.09700, 2.09700, 2.09700, 2.09700, 2.09700, 2.09700])
# bond_target = np.array([1.8, 1.8, 1.8, 1.8, 1.8, 1.8])

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