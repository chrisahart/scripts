import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import MDAnalysis as mda
from MDAnalysis.analysis import distances


"""
    Plot energy and forces for bulk hematite
"""


def load_file_coord(folder, filename, cols, del_rows=None):
    """
        Return CP2K MD .XYZ coordinate file as Pandas database.
    """
    files = ['{}/{}'.format(folder, filename)]
    file_coord = pd.read_csv(files[0], names=cols, delim_whitespace=True, on_bad_lines='skip')

    # Determine number of atoms
    num_atoms = int(float(file_coord['Species'][0]))
    if del_rows: file_coord = file_coord.drop(del_rows)
    file_coord = file_coord.reset_index(drop=True)

    # Save species as separate Series
    species = file_coord['Species'][1:num_atoms + 1]
    if del_rows: species = file_coord['Species'][:num_atoms + 2]
    species = species.reset_index(drop=True)

    # Force database to numeric, assigning any non-numeric as NaN
    # file_coord = file_coord.drop([0, 1])
    file_coord = file_coord.apply(pd.to_numeric, errors='coerce')

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


folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/resources/wiki/wiki_cp2k/tio2/anatase/cell-221/md'
folder_1 = '{}/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-dz-400k-500k-csvr-timecon-1-COMVEL_TO-1e-10-nvt-u-3.5-hole-timestep-1.0-o-3.5-ti-0'.format(folder)
files = ['tio2-1.ener', 'tio2-charges-1-clean.hirshfeld', 'tio2-pos-1.xyz', 'tio2-frc-1.xyz']
num_atoms = 48
box_size = [7.56, 7.56, 9.62, 90, 90, 90]

file_energy_1, energy_kinetic_1, energy_potential_1, energy_total_1, temperature_1, time_array2, time_per_step_1, step_1 = read_energy(folder_1, files[0])
hirshfeld_1_df, hirshfeld_1_np, _ = read_hirsh(folder_1, files[1], num_atoms)
topology_file = '{}/system.xyz'.format(folder_1)
trajectory_file = '{}/tio2-pos-1.xyz'.format(folder_1)

local_bonds = 3
timestep = 1
frames_skip = 0

num_timesteps = np.shape(hirshfeld_1_np)[0]
time_array = np.linspace(start=0, stop=num_timesteps*timestep,num=num_timesteps)
xlim_1 = [0, time_array[-1]]
xlim_1 = [0, 400]
ylim_1 = [0, 1.0]
draw_legend = False
folder_save = folder_1
plotting_colors = ['r', 'g', 'b'] * 100

# Setup md analysis environment
universe = mda.Universe(topology_file, trajectory_file)
num_timesteps2 = len(universe.trajectory)
time_array2 = np.linspace(start=0, stop=len(universe.trajectory) * timestep, num=len(universe.trajectory))
atoms_ti = universe.select_atoms('name Ti')
num_atoms_ti = len(atoms_ti)
atoms_o = universe.select_atoms('name O')
num_atoms_o = len(atoms_o)
bond_lengths_time = np.zeros((num_timesteps2, num_atoms_o, num_atoms_ti))
for ts in universe.trajectory:
    frame = universe.trajectory.frame
    bond_lengths_time[frame] = distances.distance_array(atoms_o.positions, atoms_ti.positions, box=box_size)
bond_lengths_time_sorted = np.zeros((num_timesteps2, num_atoms_o, local_bonds))
bond_lengths_time_sorted_mean = np.zeros((num_timesteps2, num_atoms_o))
for i in range(num_atoms_o):
    for j in range(num_timesteps2):
        bond_lengths_time_sorted[j, i] = np.sort(bond_lengths_time[j, i])[0:local_bonds]
bond_lengths_time_sorted_mean = np.mean(bond_lengths_time_sorted, axis=2)

# Calculate polaron atom
polaron_atom_time = np.zeros(num_timesteps, dtype=int)
for j in range(num_timesteps):
    polaron_atom_time[j] = int(np.argmax(hirshfeld_1_np[frames_skip+j, 5, :]))
polaron_atoms = polaron_atom_time[np.insert(polaron_atom_time[:-1] != polaron_atom_time[1:], 0, True)]
print('polaron_atoms', polaron_atoms+1)

# Calculate distance between current timestep polaron atom and next timestep
polaron_distances = np.zeros(num_timesteps2)
for j in range(num_timesteps - 1):
    polaron_distances[j] = distances.distance_array(universe.select_atoms('index {}'.format(polaron_atom_time[j])).positions,
                                                    universe.select_atoms('index {}'.format(polaron_atom_time[j+1])).positions,
                                                    box=box_size)
polaron_distances_hop = polaron_distances[np.nonzero(polaron_distances)]
print('polaron_distances_hop', polaron_distances_hop)

# Plot Hirshfeld spin of all atoms
fig_spin1, ax_spin1 = plt.subplots()
temp = np.zeros(num_timesteps)
for j in range(num_atoms):
    ax_spin1.plot(time_array, hirshfeld_1_np[:, 5, j], '-', label='{}'.format(j+1))
for j in range(3):
    ax_spin1.plot(time_array, hirshfeld_1_np[:, 5, polaron_atoms[j]], '-', label='{}'.format(j+1), color=plotting_colors[j])
if draw_legend: ax_spin1.legend(frameon=True)
ax_spin1.set_xlabel('Time / fs')
ax_spin1.set_ylabel('Spin moment')
ax_spin1.set_xlim(xlim_1)
ax_spin1.set_ylim(ylim_1)
fig_spin1.tight_layout()
fig_spin1.savefig('{}/hirshfeld_spin_all.png'.format(folder_save), dpi=300)

# Plot average of 6 O-Ti bonds
fig_bonds_1, ax_bonds_1 = plt.subplots()
for i in range(num_atoms_o):
    ax_bonds_1.plot(time_array2 - time_array2[0], bond_lengths_time_sorted_mean[:, i], '-', label='Fe {}'.format(i + 1))
for j in range(3):
    ax_bonds_1.plot(time_array2 - time_array2[0], bond_lengths_time_sorted_mean[:, polaron_atoms[j]-num_atoms_ti], '-', label='Fe {}'.format(j + 1), color=plotting_colors[j])
ax_bonds_1.set_xlabel('Time / fs')
ax_bonds_1.set_ylabel('Average of {} O-Ti bond lengths / A'.format(local_bonds))
ax_bonds_1.set_xlim(xlim_1)
fig_bonds_1.savefig('{}/bond_lengths_average.png'.format(folder_save), dpi=300)
fig_bonds_1.tight_layout()

if __name__ == "__main__":
    print('Finished.')
    plt.show()
