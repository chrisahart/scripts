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
                

folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-223/md/pbe-u-4.1/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-electron-csvr-timecon-1-COMVEL_TO-1e-10-nvt'
# folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-223/md/pbe-u-4.1/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-electron-u-4.0-csvr-timecon-1-COMVEL_TO-1e-10-nvt'
# folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-223/md/pbe-u-4.1/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-electron-u-3.9-csvr-timecon-1-COMVEL_TO-1e-10-nvt'
# folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-223/md/pbe-u-4.1/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-electron-u-3.8-csvr-timecon-1-COMVEL_TO-1e-10-nvt'
# folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-223/md/pbe-u-4.1/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-electron-u-3.7-csvr-timecon-1-COMVEL_TO-1e-10-nvt'
# folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-223/md/pbe-u-4.1/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-electron-u-3.6-csvr-timecon-1-COMVEL_TO-1e-10-nvt'
# folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-223/md/pbe-u-4.1/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-electron-u-3.5-csvr-timecon-1-COMVEL_TO-1e-10-nvt'
# folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-223/md/pbe-u-4.1/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-electron-u-3-csvr-timecon-1-COMVEL_TO-1e-10-nvt'
files = ['tio2-1.ener', 'tio2-charges-1-clean.hirshfeld', 'tio2-pos-1.xyz', 'tio2-frc-1.xyz']
num_atoms = 72
box_size = [9.30263, 9.30263, 9.027579, 90, 90, 90]
# u_values = np.array([4.1, 4.0, 3.9, 3.8, 3.7, 3.6, 3.5, 3.0])
# No correct hopping for any value of U from 3.0 to 4.1

# folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-334/md/pbe-u-4.1/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-electron-csvr-timecon-1-COMVEL_TO-1e-10-nvt'
# folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-334/md/pbe-u-4.1/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-electron-u-4.0-csvr-timecon-1-COMVEL_TO-1e-10-nvt'
# folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-334/md/pbe-u-4.1/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-electron-u-3.9-csvr-timecon-1-COMVEL_TO-1e-10-nvt'
# folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-334/md/pbe-u-4.1/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-electron-u-3.8-csvr-timecon-1-COMVEL_TO-1e-10-nvt'
# folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-334/md/pbe-u-4.1/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-electron-u-3.7-csvr-timecon-1-COMVEL_TO-1e-10-nvt'
# folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-334/md/pbe-u-4.1/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-electron-u-3.6-csvr-timecon-1-COMVEL_TO-1e-10-nvt'
# folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-334/md/pbe-u-4.1/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-electron-u-3.5-csvr-timecon-1-COMVEL_TO-1e-10-nvt'
# folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-334/md/pbe-u-4.1/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-electron-u-3.4-csvr-timecon-1-COMVEL_TO-1e-10-nvt'
# folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-334/md/pbe-u-4.1/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-electron-u-3.3-csvr-timecon-1-COMVEL_TO-1e-10-nvt'
# folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-334/md/pbe-u-4.1/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-electron-u-3.2-csvr-timecon-1-COMVEL_TO-1e-10-nvt'
# folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-334/md/pbe-u-4.1/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-electron-u-3.1-csvr-timecon-1-COMVEL_TO-1e-10-nvt'
# folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-334/md/pbe-u-4.1/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-electron-u-3-csvr-timecon-1-COMVEL_TO-1e-10-nvt'
# folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-334/md/pbe-u-4.1/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-csvr-timecon-1-COMVEL_TO-1e-10-nvt-dz-electron-u-3.0'
# folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-334/md/pbe-u-4.1/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-csvr-timecon-1-COMVEL_TO-1e-10-nvt-dz-electron-u-3.5'
# folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-334/md/pbe-u-4.1/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-csvr-timecon-1-COMVEL_TO-1e-10-nvt-dz-electron-u-4.1'
# folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-334/md/pbe-u-4.1/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-csvr-timecon-1-COMVEL_TO-1e-10-nvt-tz-electron-u-3.5'
# folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-334/md/pbe-u-4.1/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-csvr-timecon-1-COMVEL_TO-1e-10-nvt-tz2p-electron-u-3.5'
# folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-334/md/pbe-u-4.1/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-csvr-timecon-1-COMVEL_TO-1e-10-nvt-dz-electron-u-3.5-lowdin'
# folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-334/md/pbe-u-4.1/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-csvr-timecon-1-COMVEL_TO-1e-10-nvt-dz-electron-u-4.1-lowdin'
# folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-334/md/pbe-u-4.1/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-csvr-timecon-1-COMVEL_TO-1e-10-nvt-dz-electron-u-4.1-lowdin-COMMENSURATE'
# folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-334/md/pbe-u-4.1/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-csvr-timecon-1-COMVEL_TO-1e-10-nvt-dz-electron-u-6-lowdin'
# folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-334/md/pbe-u-4.1/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-csvr-timecon-1-COMVEL_TO-1e-10-nvt-dz-electron-u-8-lowdin'
# folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-334/md/pbe-u-4.1/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-electron-u-3.5-csvr-timecon-1-COMVEL_TO-1e-10-nvt-rs-dz-lowdin'
# folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-334/md/pbe-u-4.1/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-electron-u-3.5-csvr-timecon-1-COMVEL_TO-1e-10-nvt-rs-tz-lowdin'
# folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-334/md/pbe-u-4.1/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-electron-u-3.5-csvr-timecon-1-COMVEL_TO-1e-10-nvt-rs-tz2p-lowdin'
# files = ['tio2-1.ener', 'tio2-charges-1-clean.hirshfeld', 'tio2-pos-1.xyz', 'tio2-frc-1.xyz']
# num_atoms = 216
# box_size = [13.8, 13.8, 11.84, 90, 90, 90]
# u_values = np.array([4.1, 4.0, 3.9, 3.8, 3.7, 3.6, 3.5, 3.0])
# time_localise = np.array([0, 0, ?, ?, ?, ?, 10, 0])

# folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-336/md/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-csvr-timecon-1-COMVEL_TO-1e-10-nvt-u-3.5-dz-electron'
# Failed
# folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-336/md/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-csvr-timecon-1-COMVEL_TO-1e-10-nvt-u-3.5-dz-electron-lowdin'
# folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-336/md/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-csvr-timecon-1-COMVEL_TO-1e-10-nvt-u-3.5-electron'
# polaron_atoms [31]
# polaron_distance []
# folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-336/md/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-csvr-timecon-1-COMVEL_TO-1e-10-nvt-u-3.5-electron-dz'
# polaron_atoms [31]
# polaron_distance []
# folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-336/md/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-csvr-timecon-1-COMVEL_TO-1e-10-nvt-u-3.5-electron-dz-lowdin'
# files = ['tio2-1.ener', 'tio2-charges-1-clean.hirshfeld', 'tio2-pos-1.xyz', 'tio2-frc-1.xyz']
# num_atoms = 324
# box_size = [13.8, 13.8, 17.76, 90, 90, 90]



print('folder_1', folder_1)
file_energy_1, energy_kinetic_1, energy_potential_1, energy_total_1, temperature_1, time_val_1, time_per_step_1, step_1 = read_energy(folder_1, files[0])
hirshfeld_1_df, hirshfeld_1_np, _ = read_hirsh(folder_1, files[1], num_atoms)
topology_file = '{}/system.xyz'.format(folder_1)
trajectory_file = '{}/tio2-pos-1.xyz'.format(folder_1)

timestep = 0.5
frames_skip = 100
frames_skip = 0
num_timesteps = np.shape(hirshfeld_1_np)[0]
time_array = np.linspace(start=0, stop=num_timesteps/2,num=num_timesteps)
# xlim_1 = [0, 10e3]
xlim_1 = [0, 1e3]
# xlim_1 = [0, time_array[-1]]
ylim_1 = [0, 1.0]
draw_legend = False
folder_save = folder_1

print('Polaron:', np.max(hirshfeld_1_np[-1, 5, :]), np.argmax(hirshfeld_1_np[-1, 5, :]))
polaron_atom = np.argmax(hirshfeld_1_np[-1, 5, :])

# Plot time taken
# fig_time_md, ax_time_md = plt.subplots()
# ax_time_md.plot(time_array, time_val_1, 'k-')
# ax_time_md.set_xlabel('Time / fs')
# ax_time_md.set_ylabel('Time per MD step')
# ax_time_md.set_xlim(xlim_1)
# fig_time_md.tight_layout()
# fig_time_md.savefig('{}/time_taken_md_step.png'.format(folder_save), dpi=300)

# Plot spin of all atoms
fig_spin1, ax_spin1 = plt.subplots()
temp = np.zeros(num_timesteps)
for j in range(num_atoms):
    ax_spin1.plot(time_array, hirshfeld_1_np[:, 5, j], '-', label='{}'.format(j+1))
# ax_spin1.plot(time_array, hirshfeld_1_np[:, 5, polaron_atom], 'k-',)
if draw_legend: ax_spin1.legend(frameon=True)
# ax_spin1.set_xlabel('Timestep')
ax_spin1.set_xlabel('Time / fs')
ax_spin1.set_ylabel('Spin moment')
ax_spin1.set_xlim(xlim_1)
ax_spin1.set_ylim(ylim_1)
fig_spin1.tight_layout()
fig_spin1.savefig('{}/spin_all.png'.format(folder_save), dpi=300)

# Plot spin of all atoms in 1D chain
# fig_spin2, ax_spin2 = plt.subplots()
# temp = np.zeros(num_data_1)
# for j in range(np.shape(atoms_draw)[0]):
#     for n in range(num_data_1):
#         temp[n] = hirshfeld_1_df.loc[num_atoms * n + atoms_draw[j], 'Spin']
#     ax_spin2.plot(time_val_1 - time_val_1[0], temp, label='{}'.format(atoms_draw[j+1]))
# if draw_legend: ax_spin2.legend(frameon=True)
# # ax_spin2.set_xlabel('Timestep')
# ax_spin2.set_xlabel('Time / fs')
# ax_spin2.set_ylabel('Spin moment')
# ax_spin2.set_ylim(ylim_1)
# ax_spin2.set_xlim([0, x_end])
# fig_spin2.tight_layout()
# fig_spin2.savefig('{}/spin_atoms.png'.format(folder_save), dpi=300)

# Setup md analysis environment
local_bonds = 6
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

# Plot average of 6 Ti-O bonds
metric = np.zeros((num_atoms_ti, num_timesteps2))
fig_bonds_1, ax_bonds_1 = plt.subplots()
for i in range(num_atoms_ti):
    ax_bonds_1.plot(time_val_1 - time_val_1[0], bond_lengths_time_sorted_mean[:, i], '-', label='Fe {}'.format(i + 1))
# ax_bonds_1.plot(time_val_1 - time_val_1[0], bond_lengths_time_sorted_mean[:, polaron_atom], 'k-', label='Fe {}'.format(polaron_atom + 1))
ax_bonds_1.set_xlabel('Time / fs')
# ax_bonds_1.set_xlabel('Timestep')
ax_bonds_1.set_ylabel('Average of 6 Ti-O bond lengths / A')
# # ax_bonds_1.set_xlim([0, len(universe.trajectory)])
ax_bonds_1.set_xlim(xlim_1)
# ax_bonds_1.set_xlim([0, len(universe.trajectory) * timestep])
# ax_bonds_1.set_ylim([0.06, -0.10])
fig_bonds_1.savefig('{}/bond_lengths_average.png'.format(folder_save), dpi=300)
fig_bonds_1.tight_layout()

# Plot Ti-O bonds
# metric = np.zeros((num_atoms_ti, num_timesteps))
# fig_bonds_2, ax_bonds_2 = plt.subplots()
# for i in range(num_atoms_ti):
#     for j in range(6):
#         ax_bonds_2.plot(time_val_1 - time_val_1[0], bond_lengths_time_sorted[:, i, j], '-', label='Fe {}'.format(i + 1))
# for j in range(local_bonds):
#     ax_bonds_2.plot(time_val_1 - time_val_1[0], bond_lengths_time_sorted[:, polaron_atom, j], 'k-', label='Fe {}'.format(polaron_atom + 1))
# ax_bonds_2.set_xlabel('Time / fs')
# # ax_bonds_2.set_xlabel('Timestep')
# ax_bonds_2.set_ylabel('Ti-O bond length / A')
# # # ax_bonds_2.set_xlim([0, len(universe.trajectory)])
# ax_bonds_2.set_xlim([0, len(universe.trajectory) * timestep])
# # ax_bonds_2.set_ylim([0.06, -0.10])
# fig_bonds_2.savefig('{}/bond_lengths_average.png'.format(folder_save), dpi=300)
# fig_bonds_2.tight_layout()

# Calculate polaron atom
polaron_atom_time = np.zeros(num_timesteps-frames_skip, dtype=int)
for j in range(num_timesteps-frames_skip):
    polaron_atom_time[j] = int(np.argmax(hirshfeld_1_np[frames_skip+j, 5, :]))
polaron_atoms = np.unique(polaron_atom_time)
polaron_distances = np.zeros(np.shape(polaron_atoms)[0]-1)
for i in range(polaron_atoms.shape[0]-1):
    polaron_distances[i] = distances.distance_array(universe.select_atoms('index {}'.format(polaron_atoms[i])).positions,
                                                universe.select_atoms('index {}'.format(polaron_atoms[i+1])).positions,
                                                box=box_size)
print('polaron_atoms', polaron_atoms)
print('polaron_distance', polaron_distances)

if __name__ == "__main__":
    print('Finished.')
    plt.show()
