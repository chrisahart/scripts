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
            f.write(f"i = {index[timestep]}, time = {index[timestep]/2}, E = {energy_clean['E_pot'].values[timestep]}\n")

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


# Data
# folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/archer/bulk/221_supercell/md/hole/cleaned-checked/400k-neutral'
# files = ['hematite-1.ener', '/hematite-charges-1-clean.hirshfeld', 'hematite-pos-1.xyz', 'hematite-frc-1.xyz']

# folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/archer/bulk/221_supercell/md/hole/cleaned-checked/400k-charged-philipp'
# files = ['hematite-1.ener', 'hirshfeld.xyz', 'hematite-pos-1.xyz', 'hematite-frc-1.xyz']

# folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/archer/bulk/221_supercell/md/hole/cleaned-checked/400k-charged'
# files = ['hematite-1.ener', '/hematite-charges-1-clean.hirshfeld', 'hematite-pos-1.xyz', 'hematite-frc-1.xyz']

# folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/archer/bulk/221_supercell/md/hole/cleaned-checked/330k-charged-philipp'
# files = ['hematite-1.ener', 'hirshfeld.xyz', 'hematite-pos-1.xyz', 'hematite-frc-1.xyz']

# folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/archer/bulk/221_supercell/md/hole/cleaned-checked/400k-neutral-pbe'
# files = ['hematite-1.ener', '/hematite-charges-1-clean.hirshfeld', 'hematite-pos-1-wrapped.xyz', 'hematite-frc-1.xyz']

# folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/group/huobei'
# files = ['hematite-1.ener', '/hematite-charges-1-clean.hirshfeld', 'hematite-pos-1.xyz', 'hematite-frc-1.xyz']
# folder_2 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/archer/bulk/221_supercell/md/hole/cleaned-checked/330k-charged-philipp'
# files2 = ['hematite-1.ener', '/hematite-charges-1-clean.hirshfeld', 'input.xyz', 'hematite-frc-1.xyz']

# folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-223/md/pbe-u-4.1/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-csvr-timecon-1-COMVEL_TO-1e-10'
# folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-223/md/pbe-u-4.1/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-electron-csvr-timecon-1-COMVEL_TO-1e-10'
# folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-223/md/pbe-u-4.1/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-electron-u-3-csvr-timecon-1-COMVEL_TO-1e-10'
# folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-223/md/pbe-u-4.1/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-csvr-timecon-1-COMVEL_TO-1e-10-nvt'
# folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-223/md/pbe-u-4.1/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-electron-csvr-timecon-1-COMVEL_TO-1e-10-nvt'
# folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-223/md/pbe-u-4.1/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-electron-u-3-csvr-timecon-1-COMVEL_TO-1e-10-nvt'
# files = ['tio2-1.ener', 'tio2-charges-1-clean.hirshfeld', 'tio2-pos-1.xyz', 'tio2-frc-1.xyz']
# name = 'tio2'
# num_atoms = 72

# folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/mgo/archer/mgo/cell-222/md/pbe-u-8/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-600k-csvr-timecon-1-COMVEL_TO-1e-10-nvt'
# folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/mgo/archer/mgo/cell-222/md/pbe-u-8/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-600k-csvr-timecon-1-COMVEL_TO-1e-10-nvt-hole-u-6'
# folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/mgo/archer/mgo/cell-222/md/pbe-u-8/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-600k-csvr-timecon-1-COMVEL_TO-1e-10-nvt-hole-u-8'
# folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/mgo/archer/mgo/cell-222/md/pbe-u-8/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-600k-csvr-timecon-1-COMVEL_TO-1e-10-nvt-hole-u-8-rs'
# files = ['mgo-1.ener', 'mgo-charges-1-clean.hirshfeld', 'mgo-pos-1.xyz', 'mgo-frc-1.xyz']
# name = 'mgo'
# num_atoms = 64

# folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/mgo/archer/mgo/cell-333/md/pbe-u-8/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-600k-csvr-timecon-1-COMVEL_TO-1e-10-nvt-hole-u-8'
# files = ['mgo-1.ener', 'mgo-charges-1-clean.hirshfeld', 'mgo-pos-1.xyz', 'mgo-frc-1.xyz']
# name = 'mgo'
# num_atoms = 216

# folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-336/md'
# folder_1 = '{}/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-csvr-timecon-1-COMVEL_TO-1e-10-nvt-u-3.5-electron-dz2-u-3.0-rs-print-hse-25-schwarz-e-1e-4-f-1e-6-cfit11-cpfit3'.format(folder)  # Polaron
# files = ['tio2-1.ener', 'tio2-charges-1-clean.hirshfeld', 'tio2-pos-1.xyz', 'tio2-frc-1.xyz']
# name = 'tio2'
# num_atoms = 324

# folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-336/md-cell-opt'
# folder_1 = '{}/electron-u-ti-3.0-300k-rs-hse-25-rs-3ps-nose'.format(folder)
# files = ['tio2-1.ener', 'tio2-charges-1-clean.hirshfeld', 'tio2-pos-1.xyz', 'tio2-frc-1.xyz']
# files = ['tio2-1.ener', 'tio2-charges-1-clean.hirshfeld', 'tio2-pos-1-vmd-wrap.xyz', 'tio2-frc-1.xyz']
# name = 'tio2'
# num_atoms = 324

# folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-336/md-cell-opt/backup/4'
# folder_1 = '{}/electron-u-ti-3.0-300k-rs-hse-25-rs-3ps-nose-rs-22-2'.format(folder)
# files = ['tio2-1.ener', 'tio2-charges-1-clean.hirshfeld', 'tio2-pos-1.xyz', 'tio2-frc-1.xyz']
# files = ['tio2-1.ener', 'tio2-charges-1-clean.hirshfeld', 'tio2-pos-1-vmd-wrap.xyz', 'tio2-frc-1.xyz']
# name = 'tio2'
# num_atoms = 324

# folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-336/md-cell-opt'
# folder_1 = '{}/electron-u-ti-3.0-300k-rs'.format(folder)
# files = ['tio2-1.ener', 'tio2-charges-1-clean.hirshfeld', 'tio2-pos-1.xyz', 'tio2-frc-1.xyz']
# files = ['tio2-1.ener', 'tio2-charges-1-clean.hirshfeld', 'tio2-pos-1-vmd-wrap.xyz', 'tio2-frc-1.xyz']
# name = 'tio2'
# num_atoms = 324

folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-336/md'
folder_1 = '{}/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k'.format(folder)
files = ['tio2-1.ener', 'tio2-charges-1-clean.hirshfeld', 'tio2-pos-1-vmd-wrap.xyz', 'tio2-frc-1.xyz']
name = 'tio2'
num_atoms = 324

# Energy Hirshfeld Position Force
files_do = [True, False, True, True]
# files_do = [True, True, True, True]
# files_do = [False, False, True, False]
# energy_step_1_missing = []
# hirshfeld_step_1_missing = []
# pos_step_1_missing = []
# frc_step_1_missing = []

# Energy
if files_do[0] is True:
    file_energy_1, energy_kinetic_1, energy_potential_1, energy_total_1, temperature_1, time_val_1, time_per_step_1, step_1 = read_energy(folder_1, files[0])
    file_energy_1_no_duplicates = file_energy_1.drop_duplicates(subset=['Step'], keep='first')
    energy_unique_frames, energy_unique_indices = remove_duplicates(step_1)
    energy_step_1_missing = detect_missing_values(energy_unique_frames)
    print('Energy step last value', energy_unique_frames[-1])
    print('Energy total number steps', energy_unique_frames[-1] + 1)
    print('Energy length', len(energy_unique_frames))
    # print('Energy missing values', energy_step_1_missing)
    print('Energy length + missing values', len(energy_unique_frames)+len(energy_step_1_missing))
    # plt.plot(energy_step_1_no_duplicates, 'k.', markersize=1)

# Hirshfeld
if files_do[1] is True:
    # will break if Hirshfeld has printed but energy has not, be careful
    hirshfeld_1_df, hirshfeld_1_np, _ = read_hirsh(folder_1, files[1], num_atoms)
    hirshfeld_index = list(np.loadtxt('{}/index.hirshfeld'.format(folder_1), dtype=int))
    print(hirshfeld_1_np.shape)
    print(len(hirshfeld_index))
    # plt.plot(hirshfeld_index, 'k.', markersize=1)
    hirshfeld_index_unique_frames, hirshfeld_index_unique_indices = remove_duplicates(list(hirshfeld_index))
    hirshfeld_1_np_no_duplicates = hirshfeld_1_np[hirshfeld_index_unique_indices]
    hirshfeld_step_1_missing = detect_missing_values(hirshfeld_index_unique_frames)
    print('\nHirshfeld step last value', hirshfeld_index_unique_frames[-1])
    print('Hirshfeld total number steps', hirshfeld_index_unique_frames[-1] + 1)
    print('Hirshfeld length', len(hirshfeld_index_unique_frames))
    # print('Hirshfeld missing values', hirshfeld_step_1_missing)
    print('Hirshfeld length + missing values', len(hirshfeld_index_unique_frames)+len(hirshfeld_step_1_missing))
    print('Hirshfeld df total number steps', np.shape(hirshfeld_1_np_no_duplicates))
    # plt.plot(hirshfeld_step_1_no_duplicates, 'k.', markersize=1)
else:
    hirshfeld_step_1_missing = []
    hirshfeld_index_no_duplicates = []
    hirshfeld_1_np_no_duplicates = np.NaN

# Position
if files_do[2] is True:
    # coordinates2, coord_x2, coord_y2, coord_z2, species2, num_atoms2, num_timesteps2 = load_values_coord(folder_2, files2[2], ['Species', 'X', 'Y', 'Z'])
    coordinates, coord_x, coord_y, coord_z, species, num_atoms, num_timesteps = load_values_coord(folder_1, files[2], ['Species', 'X', 'Y', 'Z'])
    # species = species2
    # pos_index = np.linspace(start=1, stop=50, num=50, dtype=int)
    # print(pos_index)
    pos_index = np.loadtxt('{}/index_pos.txt'.format(folder_1), dtype='str')
    pos_index = [s.replace(',', '') for s in pos_index]
    pos_index = [int(s) for s in pos_index]
    pos_index_unique_frames, pos_index_unique_indices = remove_duplicates(list(pos_index))
    coordinates_no_duplicates = coordinates[pos_index_unique_indices]
    pos_step_1_missing = detect_missing_values(pos_index_unique_frames)
    print('\nPosition step last value', pos_index_unique_frames[-1])
    print('Position total number steps', pos_index_unique_frames[-1] + 1)
    print('Position length', len(pos_index_unique_frames))
    # print('Position missing values', pos_step_1_missing)
    print('Position length + missing values', len(pos_index_unique_frames)+len(pos_step_1_missing))
    print('Position df total number steps', np.shape(coordinates_no_duplicates))
    # plt.plot(coordinates_no_duplicates, 'k.', markersize=1)

# Force
if files_do[3] is True:
    forces, force_x, force_y, force_z, _, _, _ = load_values_coord(folder_1, files[3], ['Species', 'X', 'Y', 'Z'])
    frc_index = np.loadtxt('{}/index_frc.txt'.format(folder_1), dtype='str')
    frc_index = [s.replace(',', '') for s in frc_index]
    frc_index = [int(s) for s in frc_index]
    frc_index_unique_frames, frc_index_unique_indices = remove_duplicates(list(frc_index))
    frc_no_duplicates = forces[frc_index_unique_indices]
    frc_step_1_missing = detect_missing_values(frc_index_unique_frames)
    print('\nForce step last value', frc_index_unique_frames[-1])
    print('Force total number steps', frc_index_unique_frames[-1] + 1)
    print('Force length', len(frc_index_unique_frames))
    # print('Force missing values', frc_step_1_missing)
    print('Force length + missing values', len(frc_index_unique_frames)+len(frc_step_1_missing))
    # plt.plot(frc_step_1_no_duplicates, 'k.', markersize=1)

# if files_do[0] is True: print('\nEnergy missing values', energy_step_1_missing)
# if files_do[1] is True: print('Hirshfeld missing values', hirshfeld_step_1_missing)
# if files_do[2] is True: print('Position missing values', pos_step_1_missing)
# if files_do[3] is True: print('Force missing values', frc_step_1_missing)

missing_all = energy_step_1_missing + hirshfeld_step_1_missing + pos_step_1_missing + frc_step_1_missing
missing_all = (list(OrderedDict.fromkeys(missing_all)))
# print('All missing values', missing_all)
# print('\nWe will have this many values:',  frc_index_unique_frames[-1] + 1 - len(missing_all))

# Clean energy
if files_do[0] is True:
    energy_step_1_no_duplicates = [x for x in energy_unique_frames if x not in missing_all]
    energy_clean = file_energy_1_no_duplicates[file_energy_1_no_duplicates['Step'].isin(energy_step_1_no_duplicates)]
    print('Energy dataframe has this many values: ', np.shape(energy_clean))

# Clean hirshfeld
if files_do[1] is True:
    hirshfeld_index_no_duplicates = [x for x in hirshfeld_index_unique_frames if x not in missing_all]
    missing_set = set(hirshfeld_step_1_missing)
    adjusted_index_map = []
    adjusted_index = 0
    for original_index in range(max(hirshfeld_index_no_duplicates) + 1):
        if original_index not in missing_set:
            adjusted_index_map.append(adjusted_index)
            adjusted_index += 1
        else:
            adjusted_index_map.append(None)
    adjusted_indices = [adjusted_index_map[idx] for idx in hirshfeld_index_no_duplicates if adjusted_index_map[idx] is not None]
    hirshfeld_1_np_no_duplicates = hirshfeld_1_np_no_duplicates[adjusted_indices]
    print('Hirshfeld index has this many values: ', len(hirshfeld_index_no_duplicates))
    print('Hirshfeld has this many values: ', hirshfeld_1_np_no_duplicates.shape)

# Clean position
if files_do[2] is True:
    pos_index_no_duplicates = [x for x in pos_index_unique_frames if x not in missing_all]
    missing_set = set(pos_step_1_missing)
    adjusted_index_map = []
    adjusted_index = 0
    for original_index in range(max(pos_index_no_duplicates) + 1):
        if original_index not in missing_set:
            adjusted_index_map.append(adjusted_index)
            adjusted_index += 1
        else:
            adjusted_index_map.append(None)
    adjusted_indices = [adjusted_index_map[idx] for idx in pos_index_no_duplicates if adjusted_index_map[idx] is not None]
    coordinates_no_duplicates = coordinates_no_duplicates[adjusted_indices]
    print('Coordinates index has this many values: ', len(pos_index_no_duplicates))
    print('Coordinates has this many values: ', coordinates_no_duplicates.shape)

# Clean forces
if files_do[3] is True:
    frc_index_no_duplicates = [x for x in frc_index_unique_frames if x not in missing_all]
    missing_set = set(frc_step_1_missing)
    adjusted_index_map = []
    adjusted_index = 0
    for original_index in range(max(frc_index_no_duplicates) + 1):
        if original_index not in missing_set:
            adjusted_index_map.append(adjusted_index)
            adjusted_index += 1
        else:
            adjusted_index_map.append(None)
    adjusted_indices = [adjusted_index_map[idx] for idx in frc_index_no_duplicates if adjusted_index_map[idx] is not None]
    frc_no_duplicates = frc_no_duplicates[adjusted_indices]
    print('Forces index has this many values: ', len(frc_index_no_duplicates))
    print('Forces has this many values: ', frc_no_duplicates.shape)


# Truncate to min length
index_array = np.array([np.shape(energy_clean)[0], len(hirshfeld_index_no_duplicates), len(pos_index_no_duplicates), len(frc_index_no_duplicates)])
truncate_length = np.min(index_array)
print('Truncate to min length', truncate_length)
energy_clean = energy_clean[:truncate_length]
if files_do[1] is True:
    hirshfeld_1_np_no_duplicates = hirshfeld_1_np_no_duplicates[:truncate_length]
    hirshfeld_index_no_duplicates = hirshfeld_index_no_duplicates[:truncate_length]
    print('Hirshfeld index has this many values: ', len(hirshfeld_index_no_duplicates))
    print('Hirshfeld has this many values: ', hirshfeld_1_np_no_duplicates.shape)

pos_index_no_duplicates = pos_index_no_duplicates[:truncate_length]
coordinates_no_duplicates = coordinates_no_duplicates[:truncate_length]
frc_no_duplicates = frc_no_duplicates[:truncate_length]
frc_index_no_duplicates = frc_index_no_duplicates[:truncate_length]
print('\nEnergy dataframe has this many values: ', np.shape(energy_clean))

print('Coordinates index has this many values: ', len(pos_index_no_duplicates))
print('Coordinates has this many values: ', coordinates_no_duplicates.shape)

print('Forces index has this many values: ', len(frc_index_no_duplicates))
print('Forces has this many values: ', frc_no_duplicates.shape)

# Save energy
if files_do[0] is True:
    energy_clean.to_csv('{}/{}-1-cleaned.ener'.format(folder_1, name), index=False, header=False, quoting=csv.QUOTE_NONE, sep=" ",)

# Save Hirshfeld
if files_do[1] is True:
    write_hirshfeld('{}/{}-charges-1-clean-cleaned.hirshfeld'.format(folder_1, name), species, hirshfeld_1_np_no_duplicates, hirshfeld_index_no_duplicates)

# Save position
if files_do[2] is True:
    # energy_clean = np.NaN
    write_xyz('{}/{}-pos-1-cleaned.xyz'.format(folder_1, name), coordinates_no_duplicates, species, num_atoms, pos_index_no_duplicates, energy_clean)

# Save force
if files_do[3] is True:
    write_xyz('{}/{}-frc-1-cleaned.xyz'.format(folder_1, name), frc_no_duplicates, species, num_atoms, frc_index_no_duplicates, energy_clean)

if __name__ == "__main__":
    print('Finished.')
    plt.show()
