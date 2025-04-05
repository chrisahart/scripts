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
    print(file_coord)

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


def read_hirsh(folder, filename):
    """
    Read Hirshfeld analysis from CP2K output file
    """

    # Read number of atoms and labels from .xyz file
    files = ['{}/{}'.format(folder, filename)]
    cols = ['Atom', 'Element', 'Kind', 'Ref Charge', 'Pop 1', 'Pop 2', 'Spin', 'Charge', 'A', 'B']
    file_spec1 = pd.read_csv(files[0], names=cols, delim_whitespace=True, skiprows=5)

    file_spec1 = file_spec1.drop(columns=['A'])
    file_spec1 = file_spec1.drop(columns=['B'])
    file_spec1 = file_spec1.drop(columns=['Element'])

    file_spec1 = file_spec1.apply(pd.to_numeric, errors='coerce')
    file_spec1 = file_spec1.dropna(axis='rows', thresh=2)
    file_spec1 = file_spec1.dropna(axis='columns', thresh=1)
    file_spec1 = file_spec1.reset_index(drop=True)

    return file_spec1


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


def print_from_pandas3(coord_xyz, num_atoms, filename_output, save_dp='%.3f'):
    """ Print xyz from pandas dataframe of species and coordinates, adding header of number of atoms. """

    # Add number of atoms to header with blank line
    num_atoms = int(num_atoms)
    coord_xyz.loc[-1] = [None, None, None, None, None, None]
    coord_xyz.index = coord_xyz.index + 1
    coord_xyz = coord_xyz.sort_index()
    coord_xyz.loc[-1] = [str(num_atoms), None, None, None, None, None]  # Use string as Pandas would convert int to float
    coord_xyz.index = coord_xyz.index + 1
    coord_xyz = coord_xyz.sort_index()
    coord_xyz.to_csv(filename_output, index=False,header=False, quoting=csv.QUOTE_NONE, sep=" ", float_format=save_dp)


def write_xyz(filename, coordinates, species, num_atoms, index):

    num_timesteps = np.shape(coordinates)[0]

    with open(filename, 'w') as f:
        for timestep in range(num_timesteps):
            f.write(f"{num_atoms}\n")

            # Write a comment line (optional, here using timestep index)
            f.write(f"i =  {index[timestep]}\n")

            # Write each atom's species and coordinates

            for atom in range(num_atoms):
                x, y, z = coordinates[timestep, :, atom]
                f.write(f"{species[atom]} {x:.6f} {y:.6f} {z:.6f}\n")

# Data
folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/archer/bulk/221_supercell/md/hole/cleaned/400k-neutral-3'

# Energy
file_energy_1, energy_kinetic_1, energy_potential_1, energy_total_1, temperature_1, time_val_1, time_per_step_1, step_1 = read_energy(folder_1, 'hematite-1.ener')
file_energy_1_no_duplicates = file_energy_1.drop_duplicates(subset=['Step'], keep='first')
energy_step_1_no_duplicates = remove_duplicates(step_1)
energy_step_1_missing = detect_missing_values(energy_step_1_no_duplicates)
print('Energy step last value', energy_step_1_no_duplicates[-1])
print('Energy total number steps', energy_step_1_no_duplicates[-1] + 1)
print('Energy length', len(energy_step_1_no_duplicates))
print('Energy missing values', energy_step_1_missing)
print('Energy length + missing values', len(energy_step_1_no_duplicates)+len(energy_step_1_missing))
# plt.plot(energy_step_1_no_duplicates, 'k.', markersize=1)

# Hirshfeld
hirshfeld_1 = read_hirsh(folder_1, '/hematite-charges-1-clean.hirshfeld')
hirshfeld_index = list(np.loadtxt('{}/hirshfeld/index.hirshfeld'.format(folder_1), dtype=int))
# plt.plot(hirshfeld_index, 'k.', markersize=1)
hirshfeld_step_1_no_duplicates = remove_duplicates(list(hirshfeld_index))
hirshfeld_step_1_missing = detect_missing_values(list(hirshfeld_index))
print('\nHirshfeld step last value', hirshfeld_step_1_no_duplicates[-1])
print('Hirshfeld total number steps', hirshfeld_step_1_no_duplicates[-1] + 1)
print('Hirshfeld length', len(hirshfeld_step_1_no_duplicates))
print('Hirshfeld missing values', hirshfeld_step_1_missing)
print('Hirshfeld length + missing values', len(hirshfeld_step_1_no_duplicates)+len(hirshfeld_step_1_missing))
# plt.plot(hirshfeld_step_1_no_duplicates, 'k.', markersize=1)

# Pos
coordinates, coord_x, coord_y, coord_z, species, num_atoms, num_timesteps = load_values_coord(folder_1, 'hematite-pos-1.xyz', ['Species', 'X', 'Y', 'Z'])
print('coordinates[-1]', coordinates[-1])

pos_index = np.loadtxt('{}/index_pos.txt'.format(folder_1), dtype='str')
pos_index = [s.replace(',', '') for s in pos_index]
pos_index = [int(s) for s in pos_index]
pos_index_unique_frames, pos_index_unique_indices = remove_duplicates(list(pos_index))
# print(len(pos_index))
# print(len(pos_index)-len(pos_index_no_duplicates))
# print(pos_index_no_duplicates)
# print(coordinates.shape)
#
# coordinates_no_duplicates = coordinates[pos_index_no_duplicates]

# Since pos_index contains the correct frame numbers,
# we need to find the unique ones while preserving the order of their first occurrence.



# Step 2: Filter coordinates using unique frames

# Ensure that unique_frames are within the bounds of the coordinates array

# pos_index_no_duplicates = [frame for frame in unique_frames if frame < coordinates.shape[0]]
pos_index_no_duplicates = unique_indices
coordinates_no_duplicates = coordinates[pos_index_no_duplicates]


print('coordinates_no_duplicates[-1]', coordinates_no_duplicates[-1])

pos_step_1_missing = detect_missing_values(list(pos_index))
print('\nPosition step last value', pos_index_no_duplicates[-1])
print('Position total number steps', pos_index_no_duplicates[-1] + 1)
print('Position length', len(pos_index_no_duplicates))
print('Position missing values', pos_step_1_missing)
print('Position length + missing values', len(pos_index_no_duplicates)+len(pos_step_1_missing))
print('Position df total number steps', np.shape(coordinates_no_duplicates))

# Force
forces, force_x, force_y, force_z, _, _, _ = load_values_coord(folder_1, 'hematite-frc-1.xyz', ['Species', 'X', 'Y', 'Z'])
frc_index = np.loadtxt('{}/index_frc.txt'.format(folder_1), dtype='str')
frc_index = [s.replace(',', '') for s in frc_index]
frc_index = [int(s) for s in frc_index]
frc_step_1_no_duplicates = remove_duplicates(list(frc_index))
frc_no_duplicates = forces[frc_step_1_no_duplicates]
frc_step_1_missing = detect_missing_values(list(frc_index))
print('\nForce step last value', frc_step_1_no_duplicates[-1])
print('Force total number steps', frc_step_1_no_duplicates[-1] + 1)
print('Force length', len(frc_step_1_no_duplicates))
print('Force missing values', frc_step_1_missing)
print('Force length + missing values', len(frc_step_1_no_duplicates)+len(frc_step_1_missing))
# plt.plot(frc_step_1_no_duplicates, 'k.', markersize=1)
print(num_timesteps)

print('\nEnergy missing values', energy_step_1_missing)
print('Hirshfeld missing values', hirshfeld_step_1_missing)
print('Position missing values', pos_step_1_missing)
print('Force missing values', frc_step_1_missing)

missing_all = energy_step_1_missing + hirshfeld_step_1_missing + pos_step_1_missing + frc_step_1_missing
missing_all = (list(OrderedDict.fromkeys(missing_all)))
print('All missing values', missing_all)
print('\nWe will have this many values:',  frc_step_1_no_duplicates[-1] + 1 - len(missing_all))

# Clean energy
energy_step_1_no_duplicates = [x for x in energy_step_1_no_duplicates if x not in missing_all]
energy_clean = file_energy_1_no_duplicates[file_energy_1_no_duplicates['Step'].isin(energy_step_1_no_duplicates)]
print('Energy dataframe has this many values: ', np.shape(energy_clean))

# Clean position
# pos_index_no_duplicates = [x for x in pos_index_no_duplicates if x not in missing_all]
# missing_set = set(pos_step_1_missing)
# adjusted_index_map = []
# adjusted_index = 0
# for original_index in range(max(pos_index_no_duplicates) + 1):
#     if original_index not in missing_set:
#         adjusted_index_map.append(adjusted_index)
#         adjusted_index += 1
#     else:
#         adjusted_index_map.append(None)
# adjusted_indices = [adjusted_index_map[idx] for idx in pos_index_no_duplicates if adjusted_index_map[idx] is not None]
# coordinates_no_duplicates = coordinates[adjusted_indices]
# print('Coordinates has this many values: ', coordinates_no_duplicates.shape)
# print(pos_index_no_duplicates[-1])
# print(adjusted_indices[-1])
# Clean forces
# frc_index_no_duplicates = [x for x in frc_step_1_no_duplicates if x not in missing_all]
# missing_set = set(frc_step_1_missing)
# adjusted_index_map = []
# adjusted_index = 0
# for original_index in range(max(frc_index_no_duplicates) + 1):
#     if original_index not in missing_set:
#         adjusted_index_map.append(adjusted_index)
#         adjusted_index += 1
#     else:
#         adjusted_index_map.append(None)
# adjusted_indices = [adjusted_index_map[idx] for idx in frc_index_no_duplicates if adjusted_index_map[idx] is not None]
# frc_no_duplicates = forces[adjusted_indices]
# print('Forces has this many values: ', frc_no_duplicates.shape)

# Save energy
energy_clean.to_csv('{}/hematite-1-cleaned.ener'.format(folder_1), index=False, header=False, quoting=csv.QUOTE_NONE, sep=" ",)

# Save position
# system = pd.concat([file_coord_1, file_coord_2], ignore_index=True, sort=False)
# system = system.reset_index(drop=True)
# system.X = system.X.round(dp)
# system.Y = system.Y.round(dp)
# system.Z = system.Z.round(dp)
# num_atoms_1 = system.shape[0]
# print(system.shape[0])
#
# # Print to file CP2K normal
# cp2k = system.copy()
# species = 5 * 18 * ['Au'] + 3 * 18 * ['Au'] + 4 * 18 * ['Au']
# cp2k.insert(loc=0, column='A', value=species)
# print_xyz.print_from_pandas(cp2k, num_atoms_1, '{}/{}'.format(folder, output_filename_1))

write_xyz('{}/hematite-pos-1-cleaned.xyz'.format(folder_1), coordinates_no_duplicates, species, num_atoms, pos_index_no_duplicates)
print(coordinates_no_duplicates.shape)
print(coordinates_no_duplicates[-1])

# Save force

if __name__ == "__main__":
    print('Finished.')
    plt.show()
