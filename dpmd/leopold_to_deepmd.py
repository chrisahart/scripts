import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import MDAnalysis as mda
from MDAnalysis.analysis import distances
from general import parameters as param
from ase.io import read, write
# from general import load_coordinates
# from general import print_xyz
from collections import OrderedDict
import csv

"""
    Plot energy and forces for bulk hematite
"""


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


def remove_duplicates(data):
    unique_frames = []
    unique_indices = []
    seen_frames = set()

    for idx, frame in enumerate(data):
        if frame not in seen_frames:
            unique_frames.append(frame)
            unique_indices.append(idx)
            seen_frames.add(frame)

    return unique_frames, unique_indices


def detect_missing_values(data):

    full_range = set(range(data[-1]))
    step_1_set = set(data)
    missing_values = full_range - step_1_set

    return list(missing_values)


def write_xyz(filename, coordinates, species, num_atoms, index, energy_clean):

    num_timesteps = np.shape(coordinates)[0]
    # print(num_timesteps)
    # print(coordinates.shape)
    # print('index', index[-1])

    with open(filename, 'w') as f:
        for timestep in range(num_timesteps):
            f.write(f"{num_atoms}\n")

            # Write a comment line (optional, here using timestep index)
            # f.write(f"i = {index[timestep]}, time = {index[timestep]/2}, E = 0\n")
            f.write(f"i = {index[timestep]}, time = {index[timestep]/2}, E = {energy_clean[timestep]}\n")

            # Write each atom's species and coordinates

            for atom in range(num_atoms):
                x, y, z = coordinates[timestep, :, atom]
                f.write(f"{species[atom]} {x:.6f} {y:.6f} {z:.6f}\n")


def write_hirshfeld(filename, species, hirshfeld_data, hirshfeld_index_no_duplicates):

    num_timesteps = np.shape(hirshfeld_data)[0]
    # print(num_timesteps)
    # print(hirshfeld_data.shape)
    # print('hirshfeld_index_no_duplicates', hirshfeld_index_no_duplicates[-1])

    with open(filename, 'w') as f:
        for timestep in range(num_timesteps):
            f.write(f"{num_atoms}\n")

            # Write a comment line (optional, here using timestep index)
            f.write(f"i = {hirshfeld_index_no_duplicates[timestep]}, time = {hirshfeld_index_no_duplicates[timestep]/2}, E = {energy_clean['E_pot'].values[timestep]}\n")

            # Write each atom's species and coordinates
            cols = ['Atom', 'Element', 'Kind', 'Ref_Charge', 'Pop_1', 'Pop_2', 'Spin', 'Charge']
            f.write(' '.join(cols) + '\n')

            for atom in range(num_atoms):
                Atom, Kind, Ref_Charge, Pop_1, Pop_2, Spin, Charge = hirshfeld_data[timestep, :, atom]
                f.write(f"{Atom} {species[atom]} {Kind} {Ref_Charge} {Pop_1} {Pop_2} {Spin} {Charge} \n")


def read_leopold(folder, filename, num_atoms):
    """
    Read lepold data
    """

    # Read number of atoms and labels from .xyz file
    files = ['{}/{}'.format(folder, filename)]
    cols = ['Element', 'coord_x', 'coord_y', 'coord_z', 'pol_state', 'force_x', 'force_y', 'force_z', 'magmom_1',
                'magmom_2', 'magmom_3', 'magmom_4', 'toccup_1', 'toccup_2']
    file_spec1 = pd.read_csv(files[0], names=cols, delim_whitespace=True, skiprows=0)

    file_spec1 = file_spec1.dropna()
    file_spec1 = file_spec1.reset_index(drop=True)
    species = file_spec1['Element']
    file_spec1 = file_spec1.drop(columns=['Element'])

    file_spec1 = file_spec1.apply(pd.to_numeric, errors='coerce')
    file_spec1 = file_spec1.dropna()
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


folder = '/Volumes/Elements/Data/Postdoc2/Data/Work/resources/other/leopold/data/TiO2/leopold_analysis'
name = 'train'
# name = 'valid'
# name = 'test'
folder_data1 = 'database_{}'.format(name)
data = 'data/{}_clean.xyz'.format(name)
num_atoms = 288
box = np.array([[12.9730319977, 0, 0, 0, 12.9730319977, 0, 0, 0, 17.7244205475]])
topology_file = '{}/data/system.xyz'.format(folder)

# Read leopold file
leopold_1_df, leopold_1_np, species = read_leopold(folder, data, num_atoms)
energy_clean = np.loadtxt('{}/data/{}_energy.txt'.format(folder, name))

num_timesteps = np.shape(leopold_1_np)[0]
leopold_index = np.linspace(start=0, stop=num_timesteps-1, num=num_timesteps)

# Print leopold to CP2K file coordinate .xyz file
coord = leopold_1_np[:, :3, :]
forces = leopold_1_np[:, 4:7, :]
trajectory_file = '{}/data/{}-pos-1-cleaned.xyz'.format(folder, name)
write_xyz(trajectory_file, coord, species, num_atoms, leopold_index, energy_clean)

# Spin moment as atom_ener and charge state as aparam
population_alpha = leopold_1_np[:, -2, :]
population_beta = leopold_1_np[:, -1, :]
spin_moment = population_alpha-population_beta
aparam = leopold_1_np[:, 3, :]

# Convert
coord = np.transpose(coord, axes=(0, 2, 1))
forces = np.transpose(forces, axes=(0, 2, 1))
coord = coord.reshape(num_timesteps, num_atoms*3)
forces = forces.reshape(num_timesteps, num_atoms*3)
energy = np.reshape(energy_clean, (num_timesteps, 1))
box_array = np.zeros((num_timesteps, 9))
for i in range(box_array.shape[0]):
    box_array[i, :] = box
population = np.zeros((np.shape(population_alpha)[0], num_atoms, 2))
for i in range(np.shape(population_alpha)[0]):
    for j in range(num_atoms):
        population[i, j, 0] = population_alpha[i, j]
        population[i, j, 1] = population_beta[i, j]

# Training
np.save('{}/{}/set.000/energy.npy'.format(folder, folder_data1), energy)
np.save('{}/{}/set.000/coord.npy'.format(folder, folder_data1), coord)
np.save('{}/{}/set.000/force.npy'.format(folder, folder_data1), forces)
np.save('{}/{}/set.000/box.npy'.format(folder, folder_data1), box_array)
np.save('{}/{}/set.000/aparam.npy'.format(folder, folder_data1), aparam)
np.save('{}/{}/set.000/atomic_population.npy'.format(folder, folder_data1), population)

if __name__ == "__main__":
    print('Finished.')
    plt.show()
