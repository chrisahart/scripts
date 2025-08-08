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


folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/resources/other/leopold/data/TiO2/dataset'
# name = 'train'
name = 'valid'
# name = 'test'
data = '{}_clean.xyz'.format(name)
num_atoms = 288
box_size = np.array([12.9730319977, 12.9730319977, 17.7244205475, 90, 90, 90])
plot_zoom = True
topology_file = '{}/system.xyz'.format(folder)

# Read leopold file
# grep 'energy=' train.xyz  | awk -F 'energy=' '{print $2}'  | awk '{print $1}' > train_energy.txt
# sed '/Lattice=/d' train.xyz > train_clean.xyz
leopold_1_df, leopold_1_np, species = read_leopold(folder, data, num_atoms)
energy_clean = np.loadtxt('{}/{}_energy.txt'.format(folder, name))
print('leopold_1_np.shape', leopold_1_np.shape)

num_timesteps = np.shape(leopold_1_np)[0]
leopold_index = np.linspace(start=0, stop=num_timesteps-1, num=num_timesteps)
time_array = np.linspace(start=0, stop=num_timesteps-1, num=num_timesteps)
local_bonds = 6
timestep = 1

# Print leopold to CP2K file coordinate .xyz file
coordinates_no_duplicates = leopold_1_np[:, :3, :]
trajectory_file = '{}/{}-pos-1-cleaned.xyz'.format(folder, name)
write_xyz(trajectory_file, coordinates_no_duplicates, species, num_atoms, leopold_index, energy_clean)

population_alpha = leopold_1_np[:, -2, :]
population_beta = leopold_1_np[:, -1, :]
spin_moment = population_alpha-population_beta
print(population_alpha)

# Plot spin of all atoms 1
xlim_1 = [0, time_array[-1]]
ylim_1_spin = [0, 1.1]
ylim_1_bonds = [1.90, 2.11]
# fig_spin1, ax_spin1 = plt.subplots()
fig_spin1, ax_spin1 = plt.subplots(figsize=(10, 4))
temp = np.zeros(num_timesteps)
for j in range(num_atoms):
    ax_spin1.plot(time_array, spin_moment[:, j], '-', label='{}'.format(j+1))
ax_spin1.set_xlabel('Frame')
ax_spin1.set_ylabel('Spin moment')
ax_spin1.set_xlim(xlim_1)
ax_spin1.set_ylim(ylim_1_spin)
fig_spin1.tight_layout()
fig_spin1.savefig('{}/spin_all_{}.png'.format(folder, name), dpi=300)

# Setup md analysis environment
universe = mda.Universe(topology_file, trajectory_file)
num_timesteps2 = len(universe.trajectory)
time_val_1 = np.linspace(start=0, stop=len(universe.trajectory) * timestep, num=len(universe.trajectory))
atoms_ti = universe.select_atoms('name Ti')
num_atoms_ti = len(atoms_ti)
atoms_o = universe.select_atoms('name O')
num_atoms_o = len(atoms_o)
dist_arr = distances.distance_array(atoms_ti.positions, atoms_o.positions, box=box_size)
bond_lengths_time = np.zeros((num_timesteps2, num_atoms_ti, num_atoms_o))
for ts in universe.trajectory:
    frame = universe.trajectory.frame
    bond_lengths_time[frame] = distances.distance_array(atoms_ti.positions, atoms_o.positions, box=box_size)
bond_lengths_time_sorted = np.zeros((num_timesteps2, num_atoms_ti, local_bonds))
bond_lengths_time_sorted_mean = np.zeros((num_timesteps2, num_atoms_ti))
for i in range(num_atoms_ti):
    for j in range(num_timesteps2):
        bond_lengths_time_sorted[j, i] = np.sort(bond_lengths_time[j, i])[0:local_bonds]
        # bond_lengths_time_sorted_mean[j, i] = np.mean(bond_lengths_time_sorted[j, i])
bond_lengths_time_sorted_mean = np.mean(bond_lengths_time_sorted, axis=2)

# Plot average of 6 Ti-O bonds 1
metric = np.zeros((num_atoms_ti, num_timesteps2))
# fig_bonds_1, ax_bonds_1 = plt.subplots()
fig_bonds_1, ax_bonds_1 = plt.subplots(figsize=(10, 4))
for i in range(num_atoms_ti):
    ax_bonds_1.plot(time_val_1 - time_val_1[0], bond_lengths_time_sorted_mean[:, i], '-', label='Fe {}'.format(i + 1))
ax_bonds_1.set_xlabel('Frame')
ax_bonds_1.set_ylabel('Average of {} Ti-O bond lengths / A'.format(local_bonds))
ax_bonds_1.set_xlim(xlim_1)
ax_bonds_1.set_ylim(ylim_1_bonds)
fig_bonds_1.savefig('{}/bonds_{}.png'.format(folder, name), dpi=300)
fig_bonds_1.tight_layout()

# Calculate polaron atom
polaron_atom_time = np.zeros(num_timesteps, dtype=int)
for j in range(num_timesteps):
    polaron_atom_time[j] = int(np.argmax(leopold_1_np[j, -2, :]-leopold_1_np[j, -1, :]))
# polaron_atoms = np.unique(polaron_atom_time)
polaron_atoms = polaron_atom_time[np.insert(polaron_atom_time[:-1] != polaron_atom_time[1:], 0, True)]
print('polaron_atoms', polaron_atoms+1)

# Calculate distance between current timestep polaron atom and next timestep
# Then get all non-zero answer
polaron_distances = np.zeros(num_timesteps2)
for j in range(num_timesteps - 1):
    polaron_distances[j] = distances.distance_array(universe.select_atoms('index {}'.format(polaron_atom_time[j])).positions,
                                                    universe.select_atoms('index {}'.format(polaron_atom_time[j+1])).positions,
                                                    box=box_size)
polaron_distances_hop = polaron_distances[np.nonzero(polaron_distances)]
print('polaron_distance', polaron_distances)
print('polaron_distances_hop', polaron_distances_hop)

# Plot polaron distances
metric = np.zeros((num_atoms, num_timesteps2))
# fig_bonds_2, ax_bonds_2 = plt.subplots()
fig_bonds_2, ax_bonds_2 = plt.subplots(figsize=(10, 4))
# ax_bonds_2.plot(time_val_1 - time_val_1[0] - offset, polaron_distances[:-1], 'kx-')
ax_bonds_2.plot(time_val_1 - time_val_1[0], polaron_distances, 'kx-')
ax_bonds_2.set_xlabel('Time / fs')
# ax_bonds_2.set_xlabel('Timestep')
ax_bonds_2.set_ylabel('Polaron hopping distance / A')
# if draw_legend: ax_bonds_2.legend(frameon=True)
# # ax_bonds_2.set_xlim([0, len(universe.trajectory)])
ax_bonds_2.set_xlim(xlim_1)
# ax_bonds_2.set_xlim([0, len(universe.trajectory) * timestep])
# ax_bonds_2.set_ylim([0.06, -0.10])
# fig_bonds_2.savefig('{}/polaron_hopping_distance.png'.format(folder_save), dpi=300)
fig_bonds_2.tight_layout()


if plot_zoom:

    # Plot spin of all atoms 3 diabatic
    zoom_ts = 20
    # xlim_1 = [1123-zoom_ts, 1123+zoom_ts]
    # xlim_1 = [1381-zoom_ts, 1381+zoom_ts]
    xlim_1 = [1372-zoom_ts, 1372+zoom_ts]
    fig_spin3, ax_spin3 = plt.subplots()
    temp = np.zeros(num_timesteps)
    for j in range(num_atoms):
        ax_spin3.plot(time_array, leopold_1_np[:, -2, j]-leopold_1_np[:, -1, j], '-', label='{}'.format(j+1))
    ax_spin3.set_xlabel('Frame')
    ax_spin3.set_ylabel('Spin moment')
    ax_spin3.set_xlim(xlim_1)
    ax_spin3.set_ylim(ylim_1_spin)
    fig_spin3.tight_layout()
    fig_spin3.savefig('{}/spin_all3_{}.png'.format(folder, name), dpi=300)

    # Plot average of 6 Ti-O bonds 3
    fig_bonds_3, ax_bonds_3 = plt.subplots()
    for i in range(num_atoms_ti):
        ax_bonds_3.plot(time_val_1 - time_val_1[0], bond_lengths_time_sorted_mean[:, i], '-',
                        label='Fe {}'.format(i + 1))
    ax_bonds_3.set_xlabel('Frame')
    ax_bonds_3.set_ylabel('Average of {} Ti-O bond lengths / A'.format(local_bonds))
    ax_bonds_3.set_xlim(xlim_1)
    ax_bonds_3.set_ylim(ylim_1_bonds)
    fig_bonds_3.savefig('{}/bonds3_{}.png'.format(folder, name), dpi=300)
    fig_bonds_3.tight_layout()

    # Plot spin of all atoms 4 adiabatic
    # xlim_1 = [2056, 2056+40]
    xlim_1 = [532-20, 532+20]
    # ylim_1_bonds = [1.90, 2.05]
    fig_spin4, ax_spin4 = plt.subplots()
    temp = np.zeros(num_timesteps)
    for j in range(num_atoms):
        ax_spin4.plot(time_array, leopold_1_np[:, -2, j]-leopold_1_np[:, -1, j], '-', label='{}'.format(j+1))
    ax_spin4.set_xlabel('Frame')
    ax_spin4.set_ylabel('Spin moment')
    ax_spin4.set_xlim(xlim_1)
    ax_spin4.set_ylim(ylim_1_spin)
    fig_spin4.tight_layout()
    fig_spin4.savefig('{}/spin_all4_{}.png'.format(folder, name), dpi=300)

    # Plot average of 6 Ti-O bonds 4
    fig_bonds_4, ax_bonds_4 = plt.subplots()
    for i in range(num_atoms_ti):
        ax_bonds_4.plot(time_val_1 - time_val_1[0], bond_lengths_time_sorted_mean[:, i], '-',
                        label='Fe {}'.format(i + 1))
    ax_bonds_4.set_xlabel('Frame')
    ax_bonds_4.set_ylabel('Average of {} Ti-O bond lengths / A'.format(local_bonds))
    ax_bonds_4.set_xlim(xlim_1)
    ax_bonds_4.set_ylim(ylim_1_bonds)
    fig_bonds_4.savefig('{}/bonds4_{}.png'.format(folder, name), dpi=300)
    fig_bonds_4.tight_layout()

    # # Plot spin of all atoms 5
    # xlim_1 = [2296, 2296+40]
    # ylim_1_bonds = [1.90, 2.05]
    # fig_spin5, ax_spin5 = plt.subplots()
    # temp = np.zeros(num_timesteps)
    # for j in range(num_atoms):
    #     ax_spin5.plot(time_array, leopold_1_np[:, -2, j]-leopold_1_np[:, -1, j], '-', label='{}'.format(j+1))
    # ax_spin5.set_xlabel('Frame')
    # ax_spin5.set_ylabel('Spin moment')
    # ax_spin5.set_xlim(xlim_1)
    # ax_spin5.set_ylim(ylim_1_spin)
    # fig_spin5.tight_layout()
    # fig_spin5.savefig('{}/spin_all5_{}.png'.format(folder, name), dpi=300)
    #
    # # Plot average of 6 Ti-O bonds 5
    # fig_bonds_5, ax_bonds_5 = plt.subplots()
    # for i in range(num_atoms_ti):
    #     ax_bonds_5.plot(time_val_1 - time_val_1[0], bond_lengths_time_sorted_mean[:, i], '-',
    #                     label='Fe {}'.format(i + 1))
    # ax_bonds_5.set_xlabel('Frame')
    # ax_bonds_5.set_ylabel('Average of {} Ti-O bond lengths / A'.format(local_bonds))
    # ax_bonds_5.set_xlim(xlim_1)
    # ax_bonds_5.set_ylim(ylim_1_bonds)
    # fig_bonds_5.savefig('{}/bonds5_{}.png'.format(folder, name), dpi=300)
    # fig_bonds_5.tight_layout()

if __name__ == "__main__":
    print('Finished.')
    plt.show()
