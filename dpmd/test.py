import pandas as pd
import numpy as np
import glob
from matplotlib import pyplot as plt
from general import parameters as param
from matplotlib.colors import LogNorm
import dpdata
import cp2kdata
from cp2kdata import Cp2kOutput


def read_hirsh(folder, filename, num_atoms):
    """
    Read Hirshfeld analysis from CP2K output file
    """

    # Read number of atoms and labels from .xyz file
    files = ['{}/{}'.format(folder, filename)]
    cols = ['Atom', 'Element', 'Kind', 'Ref_Charge', 'Pop_1', 'Pop_2', 'Spin', 'Charge', 'A', 'B']
    file_spec1 = pd.read_csv(files[0], names=cols, delim_whitespace=True, skiprows=5)

    file_spec1 = file_spec1.drop(columns=['A'])
    file_spec1 = file_spec1.drop(columns=['B'])

    species = file_spec1['Element']
    file_spec1 = file_spec1.drop(columns=['Element'])

    file_spec1 = file_spec1.apply(pd.to_numeric, errors='coerce')
    file_spec1 = file_spec1.dropna(axis='rows', thresh=2)
    file_spec1 = file_spec1.dropna(axis='columns', thresh=1)
    file_spec1 = file_spec1.reset_index(drop=True)
    cols_new = list(file_spec1.columns)

    # Loop over each timestep and atoms
    num_timesteps = int(file_spec1.shape[0]/num_atoms)
    hirsh_data = np.zeros((num_timesteps, len(cols_new), num_atoms))

    for timestep in range(num_timesteps):
        for atom in range(num_atoms):
            for i in range(len(cols_new)):
                hirsh_data[timestep, i, atom] = file_spec1[cols_new[i]].values[atom + timestep * num_atoms]

    return file_spec1, hirsh_data, species


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


# mgo 222
# pbc set {9.2 9.2 8.88 90 90 90} -all ; pbc box
# pbc wrap -all
files = ['mgo-1-cleaned.ener', 'mgo-charges-1-clean-cleaned.hirshfeld', 'mgo-pos-1-cleaned.xyz', 'mgo-frc-1-cleaned.xyz']
box = np.array([[8.38, 0, 0, 0, 8.38, 0, 0, 0, 8.38]])
num_atoms = 64


directory = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/mgo/archer/mgo/cell-222/md/pbe-u-8/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-600k-csvr-timecon-1-COMVEL_TO-1e-10-nvt-hole-u-6'
charged = True
# 1500 to 4000 fs
# Discard first 3000 frames

size_exclude_start = int(1.5e3 * 2)
size_exclude_end = int(0e3 * 2)
size_test = int(0.1e3 * 2)
size_validation = int(0.3e3 * 2)

# the data contains 8103 frames
# the data contains 5103 frames after excluding the first and last frames
# the test data contains 200 frames
# the data contains 4903 frames after excluding the first and last frames and test set
# the training data contains 4303 frames
# the validation data contains 600 frames


# For energy just .flatten() potential energy column converted to eV Shape (200,1) -> (timesteps, 1)
print('Loading energy')
_, _, energy, _, _, time_val, time_per_step, step = read_energy(directory, files[0])
num_timesteps = np.shape(time_val)[0]
energy = np.reshape(energy, (num_timesteps, 1)) * param.hartree_to_ev

# Spin moment as atom_ener and charge state as aparam
print('Loading Hirshfeld')
_, hirshfeld, _ = read_hirsh(directory, files[1], num_atoms)
population_alpha = hirshfeld[:, 3]
population_beta = hirshfeld[:, 4]
hirshfeld = hirshfeld[:, 5]
aparam = np.zeros((num_timesteps, num_atoms))
if charged:
    for i in range(num_timesteps-1):
        spin_moment = hirshfeld[i]
        value = np.max(np.abs(spin_moment))
        index = np.argmax(np.abs(spin_moment))
        aparam[i + 1, index] = 1.0
        # print(i, index+1, value)

coordinates, _, _, _, _, _, _ = load_values_coord(directory, files[2], ['Species', 'X', 'Y', 'Z'])
coordinates = np.transpose(coordinates, axes=(0, 2, 1))
coordinates = coordinates.reshape(num_timesteps, num_atoms*3)

print('Loading forces')
forces, _, _, _, _, _, _ = load_values_coord(directory, files[3], ['Species', 'X', 'Y', 'Z'])
forces = forces * param.hartree_per_bohr_to_ev_per_angstrom
forces = np.transpose(forces, axes=(0, 2, 1))
forces = forces.reshape(num_timesteps, num_atoms*3)

# Cell information
box_array = np.zeros((num_timesteps, 9))
for i in range(box_array.shape[0]):
    box_array[i, :] = box

print("# the data contains %d frames" % energy.shape[0])

energy = energy[range(size_exclude_start, num_timesteps - size_exclude_end)]
population_alpha = population_alpha[range(size_exclude_start, num_timesteps - size_exclude_end)]
population_beta = population_beta[range(size_exclude_start, num_timesteps - size_exclude_end)]
hirshfeld = hirshfeld[range(size_exclude_start, num_timesteps - size_exclude_end)]
aparam = aparam[range(size_exclude_start, num_timesteps - size_exclude_end)]
coordinates = coordinates[range(size_exclude_start, num_timesteps - size_exclude_end)]
forces = forces[range(size_exclude_start, num_timesteps - size_exclude_end)]
box_array = box_array[range(size_exclude_start, num_timesteps - size_exclude_end)]
print("# the data contains %d frames after excluding the first and last frames" % energy.shape[0])

energy_test = energy[-size_test:]
population_alpha_test = population_alpha[-size_test:]
population_beta_test = population_beta[-size_test:]
hirshfeld_test = hirshfeld[-size_test:]
aparam_test = aparam[-size_test:]
coordinates_test = coordinates[-size_test:]
forces_test = forces[-size_test:]
box_array_test = box_array[-size_test:]
print("# the test data contains %d frames" % energy_test.shape[0])

energy = energy[:-size_test]
population_alpha = population_alpha[:-size_test]
population_beta = population_beta[:-size_test]
hirshfeld = hirshfeld[:-size_test]
aparam = aparam[:-size_test]
coordinates = coordinates[:-size_test]
forces = forces[:-size_test]
box_array = box_array[:-size_test]

print("# the data contains %d frames after excluding the first and last frames and test set" % energy.shape[0])

# Load index from file
index_training = np.loadtxt('{}/database_ener_force_train/index_training.raw'.format(directory)).astype(np.int64)
index_validation = np.loadtxt('{}/database_ener_force_test/1/index_validation.raw'.format(directory)).astype(np.int64)

# Set random index
# rng = np.random.default_rng()
# index_validation = rng.choice(energy.shape[0], size=size_validation, replace=False)
# index_validation = np.sort(index_validation)
# index_training = list(set(range(energy.shape[0])) - set(index_validation))

energy_train = energy[index_training]
population_alpha_train = population_alpha[index_training]
population_beta_train = population_beta[index_training]
hirshfeld_train = hirshfeld[index_training]
aparam_train = aparam[index_training]
coordinates_train = coordinates[index_training]
forces_train = forces[index_training]
box_array_train = box_array[index_training]

energy_valid = energy[index_validation]
population_alpha_valid = hirshfeld[index_validation]
population_beta_valid = population_beta[index_validation]
hirshfeld_valid = hirshfeld[index_validation]
aparam_valid = aparam[index_validation]
coordinates_valid = coordinates[index_validation]
forces_valid = forces[index_validation]
box_array_valid = box_array[index_validation]

print("# the training data contains %d frames" % energy_train.shape[0])
print("# the validation data contains %d frames" % energy_valid.shape[0])

# Index
np.savetxt('{}/database_ener_force_train/index_training.raw'.format(directory), index_training)
np.savetxt('{}/database_ener_force_test/1/index_validation.raw'.format(directory), index_validation)
np.savetxt('{}/database_spin_train/index_training.raw'.format(directory), index_training)
np.savetxt('{}/database_spin_test/1/index_validation.raw'.format(directory), index_validation)
np.savetxt('{}/database_population_train/index_training.raw'.format(directory), index_training)
np.savetxt('{}/database_population_test/1/index_validation.raw'.format(directory), index_validation)

# Energy
np.save('{}/database_ener_force_train/set.000/energy.npy'.format(directory), energy_train)
np.save('{}/database_ener_force_test/1/set.000/energy.npy'.format(directory), energy_valid)
np.save('{}/database_ener_force_test/2/set.000/energy.npy'.format(directory), energy_test)
np.save('{}/database_spin_train/set.000/energy.npy'.format(directory), energy_train)
np.save('{}/database_spin_test/1/set.000/energy.npy'.format(directory), energy_valid)
np.save('{}/database_spin_test/2/set.000/energy.npy'.format(directory), energy_test)
np.save('{}/database_population_train/set.000/energy.npy'.format(directory), energy_train)
np.save('{}/database_population_test/1/set.000/energy.npy'.format(directory), energy_valid)
np.save('{}/database_population_test/2/set.000/energy.npy'.format(directory), energy_test)

# Hirshfeld
np.savetxt('{}/database_spin_train/atom_ener.raw'.format(directory), hirshfeld_train.flatten(), delimiter=' ')
np.save('{}/database_spin_train/set.000/atom_ener.npy'.format(directory), hirshfeld_train.flatten())
np.savetxt('{}/database_spin_test/1/atom_ener.raw'.format(directory), hirshfeld_valid.flatten(), delimiter=' ')
np.save('{}/database_spin_test/1/set.000/atom_ener.npy'.format(directory), hirshfeld_valid.flatten())
np.savetxt('{}/database_spin_test/2/atom_ener.raw'.format(directory), hirshfeld_test.flatten(), delimiter=' ')
np.save('{}/database_spin_test/2/set.000/atom_ener.npy'.format(directory), hirshfeld_test.flatten())

# Population
np.save('{}/database_population_train/set.000/atomic_spin.npy'.format(directory), np.concatenate((population_alpha_train, population_beta_train), axis=1))
np.save('{}/database_population_test/1/set.000/atomic_spin.npy'.format(directory), np.concatenate((population_alpha_valid, population_beta_valid), axis=1))
np.save('{}/database_population_test/2/set.000/atomic_spin.npy'.format(directory), np.concatenate((population_alpha_test, population_beta_test), axis=1))

# Aparam
np.savetxt('{}/database_ener_force_train/aparam.raw'.format(directory), aparam_train.flatten(), delimiter=' ')
np.save('{}/database_ener_force_train/set.000/aparam.npy'.format(directory), aparam_train.flatten())
np.savetxt('{}/database_ener_force_test/1/aparam.raw'.format(directory), aparam_valid.flatten(), delimiter=' ')
np.save('{}/database_ener_force_test/1/set.000/aparam.npy'.format(directory), aparam_valid.flatten())
np.savetxt('{}/database_ener_force_test/2/aparam.raw'.format(directory), aparam_test.flatten(), delimiter=' ')
np.save('{}/database_ener_force_test/2/set.000/aparam.npy'.format(directory), aparam_test.flatten())
np.savetxt('{}/database_spin_train/aparam.raw'.format(directory), aparam_train.flatten(), delimiter=' ')
np.save('{}/database_spin_train/set.000/aparam.npy'.format(directory), aparam_train.flatten())
np.savetxt('{}/database_spin_test/1/aparam.raw'.format(directory), aparam_valid.flatten(), delimiter=' ')
np.save('{}/database_spin_test/1/set.000/aparam.npy'.format(directory), aparam_valid.flatten())
np.savetxt('{}/database_spin_test/2/aparam.raw'.format(directory), aparam_test.flatten(), delimiter=' ')
np.save('{}/database_spin_test/2/set.000/aparam.npy'.format(directory), aparam_test.flatten())
np.savetxt('{}/database_population_train/aparam.raw'.format(directory), aparam_train.flatten(), delimiter=' ')
np.save('{}/database_population_train/set.000/aparam.npy'.format(directory), aparam_train.flatten())
np.savetxt('{}/database_population_test/1/aparam.raw'.format(directory), aparam_valid.flatten(), delimiter=' ')
np.save('{}/database_population_test/1/set.000/aparam.npy'.format(directory), aparam_valid.flatten())
np.savetxt('{}/database_population_test/2/aparam.raw'.format(directory), aparam_test.flatten(), delimiter=' ')
np.save('{}/database_population_test/2/set.000/aparam.npy'.format(directory), aparam_test.flatten())

# Coordinates
np.save('{}/database_ener_force_train/set.000/coord.npy'.format(directory), coordinates_train)
np.save('{}/database_ener_force_test/1/set.000/coord.npy'.format(directory), coordinates_valid)
np.save('{}/database_ener_force_test/2/set.000/coord.npy'.format(directory), coordinates_test)
np.save('{}/database_spin_train/set.000/coord.npy'.format(directory), coordinates_train)
np.save('{}/database_spin_test/1/set.000/coord.npy'.format(directory), coordinates_valid)
np.save('{}/database_spin_test/2/set.000/coord.npy'.format(directory), coordinates_test)
np.save('{}/database_population_train/set.000/coord.npy'.format(directory), coordinates_train)
np.save('{}/database_population_test/1/set.000/coord.npy'.format(directory), coordinates_valid)
np.save('{}/database_population_test/2/set.000/coord.npy'.format(directory), coordinates_test)

# Forces
np.save('{}/database_ener_force_train/set.000/force.npy'.format(directory), forces_train)
np.save('{}/database_ener_force_test/1/set.000/force.npy'.format(directory), forces_valid)
np.save('{}/database_ener_force_test/2/set.000/force.npy'.format(directory), forces_test)
np.save('{}/database_spin_train/set.000/force.npy'.format(directory), forces_train)
np.save('{}/database_spin_test/1/set.000/force.npy'.format(directory), forces_valid)
np.save('{}/database_spin_test/2/set.000/force.npy'.format(directory), forces_test)
np.save('{}/database_population_train/set.000/force.npy'.format(directory), forces_train)
np.save('{}/database_population_test/1/set.000/force.npy'.format(directory), forces_valid)
np.save('{}/database_population_test/2/set.000/force.npy'.format(directory), forces_test)

# Box
np.save('{}/database_ener_force_train/set.000/box.npy'.format(directory), box_array_train)
np.save('{}/database_ener_force_test/1/set.000/box.npy'.format(directory), box_array_valid)
np.save('{}/database_ener_force_test/2/set.000/box.npy'.format(directory), box_array_test)
np.save('{}/database_spin_train/set.000/box.npy'.format(directory), box_array_train)
np.save('{}/database_spin_test/1/set.000/box.npy'.format(directory), box_array_valid)
np.save('{}/database_spin_test/2/set.000/box.npy'.format(directory), box_array_test)
np.save('{}/database_population_train/set.000/box.npy'.format(directory), box_array_train)
np.save('{}/database_population_test/1/set.000/box.npy'.format(directory), box_array_valid)
np.save('{}/database_population_test/2/set.000/box.npy'.format(directory), box_array_test)

print('Finished')
