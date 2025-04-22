import pandas as pd
import numpy as np
import glob
from matplotlib import pyplot as plt
from general import parameters as param
from matplotlib.colors import LogNorm
import dpdata
import cp2kdata
from cp2kdata import Cp2kOutput

def load_file_coord(folder, filename, cols, del_rows=None):
    """
        Return CP2K MD .XYZ coordinate file as Pandas database.
    """

    # Search for all files with path "data/*coordinates.xyz"
    # files = []
    # for file in glob.glob('{}{}{}'.format(folder, '/', filename)):
    #     files.append(file)
    #
    # if not files:
    #     print('\n No files were found, causing program to crash. \n')

    # Single file
    files = ['{}/{}'.format(folder, filename)]

    # Assign column identities
    # cols = ['Species', 'X', 'Y', 'Z']

    # Read as csv file with whitespace delimiter
    file_coord = pd.read_csv(files[0], names=cols, delim_whitespace=True, on_bad_lines='skip')
    # print(file_coord)

    # Determine number of atoms
    num_atoms = int(float(file_coord['Species'][0]))
    if del_rows: file_coord = file_coord.drop(del_rows)
    file_coord = file_coord.reset_index(drop=True)

    # print(file_coord)

    # Save species as separate Series
    species = file_coord['Species'][1:num_atoms + 1]
    if del_rows: species = file_coord['Species'][:num_atoms + 2]
    species = species.reset_index(drop=True)

    # print(species)

    # Force database to numeric, assigning any non-numeric as NaN
    # file_coord = file_coord.drop([0, 1])
    file_coord = file_coord.apply(pd.to_numeric, errors='coerce')

    # print(file_coord)

    # Filter rows with two or more NaN and columns with one of more NaN, leaving only coordinate data
    file_coord = file_coord.dropna(axis='rows', thresh=2)
    file_coord = file_coord.dropna(axis='columns', thresh=1)
    file_coord = file_coord.reset_index(drop=True)

    return file_coord, num_atoms, species


def load_values_coord(folder, filename, cols):
    """
        Return CP2K MD .XYZ file as Numpy array.
    """

    # Load coordinate data from Pandas database
    db_coord, num_atoms, species = load_file_coord(folder, filename, cols)
    coord_pandas_x = db_coord['X'].values
    coord_pandas_y = db_coord['Y'].values
    coord_pandas_z = db_coord['Z'].values

    # Assign variables
    num_timesteps = int(coord_pandas_x.shape[0] / num_atoms)

    # Initialise arrays
    coord_x = np.zeros((num_timesteps, num_atoms))
    coord_y = np.zeros((num_timesteps, num_atoms))
    coord_z = np.zeros((num_timesteps, num_atoms))
    coord = np.zeros((num_timesteps, 3, num_atoms))

    # Loop over each timestep and atoms
    for timestep in range(num_timesteps):
        for atom in range(num_atoms):

            # Re-structure coordinate arrays
            coord_x[timestep, atom] = coord_pandas_x[atom + timestep * num_atoms]
            coord_y[timestep, atom] = coord_pandas_y[atom + timestep * num_atoms]
            coord_z[timestep, atom] = coord_pandas_z[atom + timestep * num_atoms]

            coord[timestep, 0, atom] = coord_pandas_x[atom + timestep * num_atoms]
            coord[timestep, 1, atom] = coord_pandas_y[atom + timestep * num_atoms]
            coord[timestep, 2, atom] = coord_pandas_z[atom + timestep * num_atoms]

    return coord, coord_x, coord_y, coord_z, species, num_atoms, num_timesteps


def read_energy(folder, filename):
    """
        Return CP2K MD .ener file as re-structured Numpy array.
    """

    files = ['{}/{}'.format(folder, filename)]
    cols = ['Step', 'Time', 'E_kin', 'Temp', 'E_pot', 'E_tot', 'Time_per_step']
    file_energy = pd.read_csv(files[0], delim_whitespace=True, names=cols, skiprows=[0])

    # Load energy data from Pandas database
    energy_kinetic = file_energy['E_kin'].values
    energy_potential = file_energy['E_pot'].values
    energy_total = file_energy['E_tot'].values
    temperature = file_energy['Temp'].values
    step = file_energy['Step'].values
    time = file_energy['Time'].values
    time_per_step = file_energy['Time_per_step'].values

    return file_energy, energy_kinetic, energy_potential, energy_total, temperature, time, time_per_step, step


directory = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/archer/bulk/221_supercell/md/hole/cleaned-checked/400k-neutral'
files = ['hematite-1.ener', '/hematite-charges-1-clean.hirshfeld', 'hematite-pos-1.xyz', 'hematite-frc-1.xyz']
files_do = [True, False, True, True]
box = np.array([[10.071, 0, 0, -5.035, 8.722, 0, 0, 0, 13.747]])
num_atoms = 120

size_exclude_start = 400
size_exclude_end = 0
size_test = 200
size_validation = 800

# For energy just .flatten() potential energy column converted to eV Shape (200,1) -> (timesteps, 1)
_, _, energy, _, _, time_val, time_per_step, step = read_energy(directory, files[0])
num_timesteps = np.shape(time_val)[0]
energy = np.reshape(energy, (num_timesteps, 1)) * param.hartree_to_ev

# Shape (200,360) -> (timesteps, num_atoms*3) where it is arranged [x1, y1, z1, x2, y2, z2 ... xn, yn, zn]
coordinates, _, _, _, _, _, _ = load_values_coord(directory, files[2], ['Species', 'X', 'Y', 'Z'])
coordinates = np.transpose(coordinates, axes=(0, 2, 1))
coordinates = coordinates.reshape(num_timesteps, num_atoms*3)
# np.save('/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/archer/bulk/221_supercell/md/hole/cleaned-checked/400k-neutral/database_ener_force_test/2/set.000/coord2.npy', coordinates)
# coordinates_load = np.load('/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/archer/bulk/221_supercell/md/hole/cleaned-checked/400k-neutral/database_ener_force_test/2/set.000/coord2.npy')
# print(coordinates_load.shape)
# print(coordinates_load[0, 0:3])

# Shape (200,360) -> (timesteps, num_atoms*3) where it is arranged [x1, y1, z1, x2, y2, z2 ... xn, yn, zn]
forces, _, _, _, _, _, _ = load_values_coord(directory, files[3], ['Species', 'X', 'Y', 'Z'])
forces = np.transpose(forces, axes=(0, 2, 1))
forces = forces.reshape(num_timesteps, num_atoms*3)
# np.save('/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/archer/bulk/221_supercell/md/hole/cleaned-checked/400k-neutral/database_ener_force_test/2/set.000/force2.npy', forces)
# forces_load = np.load('/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/archer/bulk/221_supercell/md/hole/cleaned-checked/400k-neutral/database_ener_force_test/2/set.000/force2.npy')
# print(forces.shape)
# print(forces[0, 0:3])

# Cell information
box_array = np.zeros((num_timesteps, 9))
for i in range(box_array.shape[0]):
    box_array[i, :] = box
# np.save('/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/archer/bulk/221_supercell/md/hole/cleaned-checked/400k-neutral/database_ener_force_test/2/set.000/box2.npy', box_array)
# box_load = np.load('/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/archer/bulk/221_supercell/md/hole/cleaned-checked/400k-neutral/database_ener_force_test/2/set.000/box2.npy')
# print(box_load.shape)
# print(box_load[0])

print("# the data contains %d frames" % energy.shape[0])

energy = energy[range(size_exclude_start, num_timesteps - size_exclude_end)]
coordinates = coordinates[range(size_exclude_start, num_timesteps - size_exclude_end)]
forces = forces[range(size_exclude_start, num_timesteps - size_exclude_end)]
box_array = box_array[range(size_exclude_start, num_timesteps - size_exclude_end)]
print("# the data contains %d frames after excluding the first and last frames" % energy.shape[0])

energy_test = energy[-size_test:]
coordinates_test = coordinates[-size_test:]
forces_test = forces[-size_test:]
box_array_test = box_array[-size_test:]
print("# the test data contains %d frames" % energy_test.shape[0])

energy = energy[:-size_test]
coordinates = coordinates[:-size_test]
forces = forces[:-size_test]
box_array = box_array[:-size_test]

print("# the data contains %d frames after excluding the first and last frames and test set" % energy.shape[0])

rng = np.random.default_rng()
index_validation = rng.choice(energy.shape[0], size=size_validation, replace=False)
index_training = list(set(range(energy.shape[0])) - set(index_validation))

energy_train = energy[index_training]
coordinates_train = coordinates[index_training]
forces_train = forces[index_training]
box_array_train = box_array[index_training]

energy_valid = energy[index_validation]
coordinates_valid = coordinates[index_validation]
forces_valid = forces[index_validation]
box_array_valid = box_array[index_validation]

print("# the training data contains %d frames" % energy_train.shape[0])
print("# the validation data contains %d frames" % energy_valid.shape[0])

np.save('/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/archer/bulk/221_supercell/md/hole/cleaned-checked/400k-neutral/database_ener_force_test/2/set.000/energy2.npy', energy)
# energy_load = np.load('/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/archer/bulk/221_supercell/md/hole/cleaned-checked/400k-neutral/database_ener_force_test/2/set.000/energy2.npy')
# print(energy_load.shape)
# print(energy_load[0])

print('Finished')
