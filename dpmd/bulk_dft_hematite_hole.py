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
from MDAnalysis.analysis import rdf
from scipy.optimize import curve_fit

"""
    Plot energy and forces for bulk hematite
"""

# Standard
# params = {'axes.formatter.limits': [-4, 4],
#           'axes.labelsize': '14',
#           'axes.titlesize': '16',
#           'lines.markersize': '8',
#           }

# TOC graphic
params = {'axes.formatter.limits': [-4, 4],
          'axes.labelsize': '14',
          'axes.titlesize': '16',
          'lines.markersize': '8',
          }
plt.rcParams.update(params)


# Define the function to fit: y = mx
def linear_func(x, m):
    return m * x


# Define the function to fit: y = mx + c
def linear_func2(x, m, c):
    return m * x + c


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


offset = 0
xlim_auto = True
plot_msd = False

# Hematite 441 hole yungyu
# folder_1 = '/Volumes/Elements/Data/Postdoc2/Data/Work/calculations/hematite/fe2o3/441/md/hole/300k-hse-25-csvr-1-cfit10-nodes-2-screen-t'  # 1.6 ps no polaron
# folder_1 = '/Volumes/Elements/Data/Postdoc2/Data/Work/calculations/hematite/fe2o3/441/md/hole/300k-hse-25-csvr-1-cfit10-nodes-2-screen-f'  # 35 fs no polaron
# folder_1 = '/Volumes/Elements/Data/Postdoc2/Data/Work/calculations/hematite/fe2o3/441/md/hole/300k-hse-25-csvr-1-cfit11-nodes-2-screen-t'  # 90 fs no polaron
# folder_1 = '/Volumes/Elements/Data/Postdoc2/Data/Work/calculations/hematite/fe2o3/441/md/hole/300k-hse-25-csvr-1-cfit11-nodes-2-screen-f'  # 500 fs no polaron
# folder_1 = '/Volumes/Elements/Data/Postdoc2/Data/Work/calculations/hematite/fe2o3/441/md/hole/300k-hse-25-csvr-1-cfit10-nodes-2-screen-f-hse-12'  # 22 fs no polaron ? confusing
# files = ['hematite-1.ener', 'hematite-charges-1-clean.hirshfeld', 'hematite-pos-1.xyz', 'hematite-frc-1.xyz']
# num_atoms = 480
# temperature_set = 300

# Hematite 441 hole archer
# folder_1 = '/Volumes/Elements/Data/Postdoc2/Data/Work/calculations/hematite/fe2o3/archer/441/md/hole/300k-hse-25-csvr-1-cfit10-nodes-2-screen-t'  # 77 fs no polaron
# folder_1 = '/Volumes/Elements/Data/Postdoc2/Data/Work/calculations/hematite/fe2o3/archer/441/md/hole/300k-hse-25-csvr-1-cfit10-nodes-2-screen-f'  # 159 fs no polaron
# folder_1 = '/Volumes/Elements/Data/Postdoc2/Data/Work/calculations/hematite/fe2o3/archer/441/md/hole/300k-hse-25-csvr-1-cfit11-nodes-2-screen-t'  # 76 fs no polaron
# folder_1 = '/Volumes/Elements/Data/Postdoc2/Data/Work/calculations/hematite/fe2o3/archer/441/md/hole/300k-hse-25-csvr-1-cfit11-nodes-2-screen-f'  # 60 fs no polaron
# folder_1 = '/Volumes/Elements/Data/Postdoc2/Data/Work/calculations/hematite/fe2o3/archer/441/md/hole/300k-hse-25-csvr-1-cfit10-nodes-2-screen-f-hse-12'  # 50 fs no polaron
# folder_1 = '/Volumes/Elements/Data/Postdoc2/Data/Work/calculations/hematite/fe2o3/archer/441/md/hole/300k-hse-25-csvr-1-cfit10-nodes-2-screen-f-hse-12-geo-opt'  # 158 fs no polaron
# folder_1 = '/Volumes/Elements/Data/Postdoc2/Data/Work/calculations/hematite/fe2o3/archer/441/md/hole/300k-hse-25-csvr-1-cfit10-nodes-2-screen-f-hse-12-nooffset'  # 72 fs no polaron
# folder_1 = '/Volumes/Elements/Data/Postdoc2/Data/Work/calculations/hematite/fe2o3/archer/441/md/hole/300k-hse-25-csvr-1-cfit10-nodes-2-screen-f-hse-12-nooffset-novelrs'  # 50 fs no polaron
# folder_1 = '/Volumes/Elements/Data/Postdoc2/Data/Work/calculations/hematite/fe2o3/archer/441/md/hole/300k-hse-25-csvr-1-cfit10-nodes-2-screen-f-hse-12-novelrs'  # 50 fs no polaron
# folder_1 = '/Volumes/Elements/Data/Postdoc2/Data/Work/calculations/hematite/fe2o3/archer/441/md/hole/300k-hse-25-csvr-1-cfit10-nodes-2-screen-f-hse-12-nooffset-novelrs-eps-1e-6'  # works **
# files = ['hematite-1.ener', 'hematite-charges-1-clean.hirshfeld', 'hematite-pos-1.xyz', 'hematite-frc-1.xyz']
# num_atoms = 480
# temperature_set = 300

# Hematite 441 electron yungyu
# folder_1 = '/Volumes/Elements/Data/Postdoc2/Data/Work/calculations/hematite/fe2o3/441/md/electron/300k-hse-25-csvr-1-cfit10-nodes-2-screen-t'  # polaron, but yungyu failed time x10 too slow
# folder_1 = '/Volumes/Elements/Data/Postdoc2/Data/Work/calculations/hematite/fe2o3/441/md/electron/300k-hse-25-csvr-1-cfit10-nodes-2-screen-t-2'  # polaron, but yungyu failed time x10 too slow
# folder_1 = '/Volumes/Elements/Data/Postdoc2/Data/Work/calculations/hematite/fe2o3/441/md/electron/300k-hse-25-csvr-1-cfit10-nodes-2-screen-f'  # polaron, but some sort of convergence issue
# folder_1 = '/Volumes/Elements/Data/Postdoc2/Data/Work/calculations/hematite/fe2o3/441/md/electron/300k-hse-25-csvr-1-cfit11-nodes-2-screen-t'  # polaron, but some sort of convergence issue
# folder_1 = '/Volumes/Elements/Data/Postdoc2/Data/Work/calculations/hematite/fe2o3/441/md/electron/300k-hse-25-csvr-1-cfit11-nodes-2-screen-f'  # polaron, best so far but still convergence issues
# files = ['hematite-1.ener', 'hematite-charges-1-clean.hirshfeld', 'hematite-pos-1.xyz', 'hematite-frc-1.xyz']
# num_atoms = 480
# temperature_set = 300

# Hematite 441 electron archer
# folder_1 = '/Volumes/Elements/Data/Postdoc2/Data/Work/calculations/hematite/fe2o3/archer/441/md/electron/300k-hse-25-csvr-1-cfit10-nodes-2-screen-f'  # polaron, but some sort of convergence issue
# folder_1 = '/Volumes/Elements/Data/Postdoc2/Data/Work/calculations/hematite/fe2o3/archer/441/md/electron/300k-hse-25-csvr-1-cfit10-nodes-2-screen-t'  # polaron, best so far though some convergence issues
# folder_1 = '/Volumes/Elements/Data/Postdoc2/Data/Work/calculations/hematite/fe2o3/archer/441/md/electron/300k-hse-25-csvr-1-cfit10-nodes-2-screen-t-rs'  # polaron, good **
# folder_1 = '/Volumes/Elements/Data/Postdoc2/Data/Work/calculations/hematite/fe2o3/archer/441/md/electron/300k-hse-25-csvr-1-cfit11-nodes-2-screen-f'  # polaron, but some sort of convergence issue
# folder_1 = '/Volumes/Elements/Data/Postdoc2/Data/Work/calculations/hematite/fe2o3/archer/441/md/electron/300k-hse-25-csvr-1-cfit11-nodes-2-screen-t'  # polaron, but some sort of convergence issue
# files = ['hematite-1.ener', 'hematite-charges-1-clean.hirshfeld', 'hematite-pos-1.xyz', 'hematite-frc-1.xyz']
# num_atoms = 480
# temperature_set = 300

# Hematite 221 hole yungyu
# folder_1 = '/Volumes/Elements/Data/Postdoc2/Data/Work/calculations/hematite/fe2o3/221/md/hole/300k-hse-25-csvr-1-cfit10-screen-scf-1e-5'  # 10 ps HSE(25%), 9 hops, 8 adiabatic
# folder_1 = '/Volumes/Elements/Data/Postdoc2/Data/Work/calculations/hematite/fe2o3/221/md/hole/300k-hse-25-csvr-1-cfit10-screen-scf-1e-5-rs-20'  # 10 HSE(20%) 4 hops, 2 adiabatic
# folder_1 = '/Volumes/Elements/Data/Postdoc2/Data/Work/calculations/hematite/fe2o3/221/md/hole/300k-hse-25-csvr-1-cfit10-screen-scf-1e-5-rs-20-rs-tight-nose-nodes-2'  # 10 ps HSE(20%) 8 hops, 6 adiabatic
# folder_1 = '/Volumes/Elements/Data/Postdoc2/Data/Work/calculations/hematite/fe2o3/221/md/hole/300k-hse-12-from-hse-20-step-5000'  # 3 ps HSE(12%) from HSE(25%) step 5000, 7 hops, 1 adiabatic **
# folder_1 = '/Volumes/Elements/Data/Postdoc2/Data/Work/calculations/hematite/fe2o3/221/md/hole/300k-hse-12-from-hse-20-step-8000'  # 3 ps HSE(12%) from HSE(25%) step 5000, 7 hops, 3 adiabatic, bad cancel
# folder_1 = '/Volumes/Elements/Data/Postdoc2/Data/Work/calculations/hematite/fe2o3/221/md/hole/300k-hse-15-from-hse-20-step-5000'  # 3 ps HSE(12%) from HSE(25%) step 5000, 3 hops, 1 adiabatic **
# files = ['hematite-1.ener', 'hematite-charges-1-clean.hirshfeld', 'hematite-pos-1.xyz', 'hematite-frc-1.xyz']
# num_atoms = 120
# temperature_set = 300

# Hematite 441 hole yungyu
# folder_1 = '/Volumes/Elements/Data/Postdoc2/Data/Work/calculations/hematite/fe2o3/441/md/hole/300k-hse-25-csvr-1-cfit10-nodes-2-screen-f-eps-1e-6'  # 2.2 ps HSE(25%), 0 hops
# folder_1 = '/Volumes/Elements/Data/Postdoc2/Data/Work/calculations/hematite/fe2o3/441/md/hole/300k-hse-25-csvr-1-cfit10-nodes-2-screen-t-eps-1e-6'  # 2.2 ps HSE(25%), 1 hops **
# files = ['hematite-1.ener', 'hematite-charges-1-clean.hirshfeld', 'hematite-pos-1.xyz', 'hematite-frc-1.xyz']
# num_atoms = 480
# temperature_set = 300

# Hematite 441 electron yungyu
# folder_1 = '/Volumes/Elements/Data/Postdoc2/Data/Work/calculations/hematite/fe2o3/441/md/electron/300k-hse-25-csvr-1-cfit10-nodes-2-screen-f-eps-1e-6'  # 1.2 ps HSE(25%), 5 hops, 5 adiabatic
# folder_1 = '/Volumes/Elements/Data/Postdoc2/Data/Work/calculations/hematite/fe2o3/441/md/electron/300k-hse-25-csvr-1-cfit10-nodes-2-screen-t-3-eps-1e-6'  # 2.6 ps HSE(25%), 6 hops, 6 adiabatic **
# files = ['hematite-1.ener', 'hematite-charges-1-clean.hirshfeld', 'hematite-pos-1.xyz', 'hematite-frc-1.xyz']
# num_atoms = 480
# temperature_set = 300

ylim_1_spin_o = [0, 1]
ylim_1_spin_fe = [-4.25, -3]
ylim_1_time = [10, 500]
draw_legend = False
folder_save = folder_1

print('folder_1', folder_1)
file_energy_1, energy_kinetic_1, energy_potential_1, energy_total_1, temperature_1, time_val_1, time_per_step_1, step_1 = read_energy(folder_1, files[0])
time_val_energy = time_val_1
hirshfeld_1_df, hirshfeld_1_np, _ = read_hirsh(folder_1, files[1], num_atoms)
print(hirshfeld_1_df)
calc_distance = False
calc_distance = True
save_fig = False
save_fig = True
topology_file = '{}/system.xyz'.format(folder_1)
trajectory_file = '{}/{}'.format(folder_1, files[2])

local_bonds = 6
timestep = 1
frames_skip = 0

num_timesteps = np.shape(hirshfeld_1_np)[0]
time_array = np.linspace(start=0, stop=num_timesteps*timestep,num=num_timesteps)
num_timesteps2 = np.shape(time_per_step_1)[0]
time_array2 = np.linspace(start=0, stop=num_timesteps*timestep, num=num_timesteps2)
print('num_timesteps np.shape(hirshfeld_1_np)[0]', num_timesteps)
print('num_timesteps2 np.shape(time_per_step_1)[0]', num_timesteps2)
if xlim_auto: xlim_1 = np.array([0, np.shape(hirshfeld_1_np)[0]])

# Plot time taken
fig_time_md, ax_time_md = plt.subplots(figsize=(10, 4))
ax_time_md.plot(time_array2, time_per_step_1, 'k-')
ax_time_md.set_xlabel('Time / fs')
ax_time_md.set_ylabel('Time per MD step')
ax_time_md.set_xlim(xlim_1)
ax_time_md.set_ylim(ylim_1_time)
fig_time_md.tight_layout()
if save_fig: fig_time_md.savefig('{}/time_taken_md_step.png'.format(folder_save), dpi=param.save_dpi)

# Plot energy
fig_energy_potential, ax_energy_potential = plt.subplots(figsize=(10, 4))
ax_energy_potential.plot(time_array2, energy_potential_1, 'k-')
ax_energy_potential.set_xlabel('Time / fs')
ax_energy_potential.set_ylabel('Energy / au')
ax_energy_potential.set_xlim(xlim_1)
fig_energy_potential.tight_layout()
if save_fig: fig_energy_potential.savefig('{}/energy_potential.png'.format(folder_save), dpi=param.save_dpi)

# Plot Hirshfeld spin o
fig_spin_o, ax_spin_o = plt.subplots(figsize=(10, 4))
temp = np.zeros(num_timesteps)
for j in range(num_atoms):
    ax_spin_o.plot(time_array/1e3, hirshfeld_1_np[:, 5, j], '-', label='{}'.format(j+1))
if draw_legend: ax_spin_o.legend(frameon=True)
ax_spin_o.set_xlabel('Time / ps')
ax_spin_o.set_ylabel('Spin moment')
ax_spin_o.set_xlim(np.array(xlim_1)/1e3)
ax_spin_o.set_ylim(ylim_1_spin_o)
fig_spin_o.tight_layout()
if save_fig: fig_spin_o.savefig('{}/hirshfeld_spin_o.png'.format(folder_save), dpi=param.save_dpi)

# Plot Hirshfeld spin Fe
fig_spin_fe, ax_spin_fe = plt.subplots(figsize=(10, 4))
temp = np.zeros(num_timesteps)
for j in range(num_atoms):
    ax_spin_fe.plot(time_array/1e3, hirshfeld_1_np[:, 5, j], '-', label='{}'.format(j+1))
if draw_legend: ax_spin_fe.legend(frameon=True)
ax_spin_fe.set_xlabel('Time / ps')
ax_spin_fe.set_ylabel('Spin moment')
ax_spin_fe.set_xlim(np.array(xlim_1)/1e3)
ax_spin_fe.set_ylim(ylim_1_spin_fe)
fig_spin_fe.tight_layout()
if save_fig: fig_spin_fe.savefig('{}/hirshfeld_spin_fe.png'.format(folder_save), dpi=param.save_dpi)

if __name__ == "__main__":
    print('Finished.')
    plt.show()
