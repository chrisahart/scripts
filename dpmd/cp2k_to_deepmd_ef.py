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


# directory = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/archer/bulk/221_supercell/md/hole/cleaned-checked/400k-neutral'
# size_exclude_start = 400
# size_exclude_end = 0
# size_test = 200
# size_validation = 800
# the data contains 2946 frames
# the data contains 2546 frames after excluding the first and last frames
# the test data contains 200 frames
# the data contains 2346 frames after excluding the first and last frames and test set
# the training data contains 1546 frames
# the validation data contains 800 frames

# directory = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/archer/bulk/221_supercell/md/hole/cleaned-checked/400k-neutral-pbe'
# 39154.000000-1674.500000 = 37,479.5 = 37 ps ~ 74,000 frames

# database-train-frames-200
# size_exclude_start = int(10e3*2)
# size_exclude_end = int(23.5e3*2+(1960-200))
# size_test = int(1.5e3*2)
# size_validation = int(1.5e3*2)
# the data contains 74960 frames
# the data contains 6200 frames after excluding the first and last frames
# the test data contains 3000 frames
# the data contains 3200 frames after excluding the first and last frames and test set
# the training data contains 200 frames
# the validation data contains 3000 frames

# database-train-frames-2000
# size_exclude_start = int(10e3*2)
# size_exclude_end = int(23.5e3*2)
# size_test = int(1.5e3*2)
# size_validation = int(1.5e3*2)
# the data contains 74960 frames
# the data contains 7960 frames after excluding the first and last frames
# the test data contains 3000 frames
# the data contains 4960 frames after excluding the first and last frames and test set
# the training data contains 1960 frames
# the validation data contains 3000 frames

# database-train-frames-22000
# size_exclude_start = int(10e3*2)
# size_exclude_end = int(13.5e3*2)
# size_test = int(1.5e3*2)
# size_validation = int(1.5e3*2)
# the data contains 74960 frames
# the data contains 27960 frames after excluding the first and last frames
# the test data contains 3000 frames
# the data contains 24960 frames after excluding the first and last frames and test set
# the training data contains 21960 frames
# the validation data contains 3000 frames

# files = ['hematite-1.ener', '/hematite-charges-1-clean.hirshfeld', 'hematite-pos-1.xyz', 'hematite-frc-1.xyz']
# pbc set {10.071 10.071 13.747 90 90 120} -all
# pbc box
# pbc wrap -all
# box = np.array([[10.071, 0, 0, -5.035, 8.722, 0, 0, 0, 13.747]])
# num_atoms = 120

# tio2
# pbc set {9.2 9.2 8.88 90 90 90} -all ; pbc box
# pbc wrap -all
# files = ['tio2-1.ener', 'tio2-charges-1-clean.hirshfeld', 'tio2-pos-1.xyz', 'tio2-frc-1.xyz']
# files = ['tio2-1-cleaned.ener', 'tio2-charges-1-clean-cleaned.hirshfeld', 'tio2-pos-1-cleaned.xyz', 'tio2-frc-1-cleaned.xyz']
# box = np.array([[9.2, 0, 0, 0, 9.2, 0, 0, 0, 8.88]])
# num_atoms = 72

# total frames 15,882 ~ 8 ps NVE 400 K neutral 223
# directory = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-223/md/pbe-u-4.1/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-csvr-timecon-1-COMVEL_TO-1e-10'
# total frames 55,637 ~ 22 ps NVT 400 K neutral 223
# directory = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-223/md/pbe-u-4.1/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-csvr-timecon-1-COMVEL_TO-1e-10-nvt'

# database-train-frames-2000
# size_exclude_start = int(5e3*2)
# size_exclude_end = int(0e3*2)
# size_test = int(1e3*2)
# size_validation = int(1e3*2)
# the data contains 15882 frames
# the data contains 5882 frames after excluding the first and last frames
# the test data contains 2000 frames
# the data contains 3882 frames after excluding the first and last frames and test set
# the training data contains 1882 frames
# the validation data contains 2000 frames

# total frames 14,400 ~ 7 ps
# directory = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-223/md/pbe-u-4.1/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-electron-csvr-timecon-1-COMVEL_TO-1e-10'
# total frames 14,400 ~ 7 ps
# directory = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-223/md/pbe-u-4.1/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-electron-csvr-timecon-1-COMVEL_TO-1e-10-nvt'

# database-train-frames-2000
# size_exclude_start = int(4.2e3*2)
# size_exclude_end = int(0e3*2)
# size_test = int(1e3*2)
# size_validation = int(1e3*2)
# the data contains 14400 frames
# the data contains 6400 frames after excluding the first and last frames
# the test data contains 2000 frames
# the data contains 4400 frames after excluding the first and last frames and test set
# the training data contains 2400 frames
# the validation data contains 2000 frames

# directory = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-223/md/pbe-u-4.1/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-electron-u-3-csvr-timecon-1-COMVEL_TO-1e-10'
# total frames 9060 ~ 4 ps

# database-train-frames-2000
# size_exclude_start = int(1.5e3*2)
# size_exclude_end = int(0e3*2)
# size_test = int(1e3*2)
# size_validation = int(1e3*2)
# the data contains 14400 frames
# the data contains 6400 frames after excluding the first and last frames
# the test data contains 2000 frames
# the data contains 4400 frames after excluding the first and last frames and test set
# the training data contains 2400 frames
# the validation data contains 2000 frames

# mgo 222
# pbc set {9.2 9.2 8.88 90 90 90} -all ; pbc box
# pbc wrap -all
files = ['mgo-1-cleaned.ener', 'mgo-charges-1-clean-cleaned.hirshfeld', 'mgo-pos-1-cleaned.xyz', 'mgo-frc-1-cleaned.xyz']
box = np.array([[8.38, 0, 0, 0, 8.38, 0, 0, 0, 8.38]])
num_atoms = 64

# directory = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/mgo/archer/mgo/cell-222/md/pbe-u-8/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-600k-csvr-timecon-1-COMVEL_TO-1e-10-nvt'
# charged = False
# 500 to 2000 fs
# Discard first 1000 frames

# size_exclude_start = int(0.5e3*2)
# size_exclude_end = int(0e3*2)
# size_test = int(0.1e3*2)
# size_validation = int(0.3e3*2)

# the data contains 3830 frames
# the data contains 2830 frames after excluding the first and last frames
# the test data contains 200 frames
# the data contains 2630 frames after excluding the first and last frames and test set
# the training data contains 2030 frames
# the validation data contains 600 frames

# directory = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/mgo/archer/mgo/cell-222/md/pbe-u-8/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-600k-csvr-timecon-1-COMVEL_TO-1e-10-nvt-hole-u-6'
# charged = True
# 1500 to 4000 fs
# Discard first 3000 frames

# size_exclude_start = int(1.5e3 * 2)
# size_exclude_end = int(0e3 * 2)
# size_test = int(0.1e3 * 2)
# size_validation = int(0.3e3 * 2)

# the data contains 8103 frames
# the data contains 5103 frames after excluding the first and last frames
# the test data contains 200 frames
# the data contains 4903 frames after excluding the first and last frames and test set
# the training data contains 4303 frames
# the validation data contains 600 frames

# directory = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/mgo/archer/mgo/cell-222/md/pbe-u-8/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-600k-csvr-timecon-1-COMVEL_TO-1e-10-nvt-hole-u-8'
# charged = True
# 2000 to 3000 fs
# Discard first 4000 frames

# size_exclude_start = int(2e3*2)
# size_exclude_end = int(0e3*2)
# size_test = int(0.1e3*2)
# size_validation = int(0.3e3*2)

# the data contains 6283 frames
# the data contains 2283 frames after excluding the first and last frames
# the test data contains 200 frames
# the data contains 2083 frames after excluding the first and last frames and test set
# the training data contains 1483 frames
# the validation data contains 600 frames

# directory = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/mgo/archer/mgo/cell-222/md/pbe-u-8/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-600k-csvr-timecon-1-COMVEL_TO-1e-10-nvt-hole-u-8-rs'
# charged = True
# 0 to 1800 fs
# Discard first 100 frames and last 1800-1250*2 frames

# size_exclude_start = int(0.05e3*2)
# size_exclude_end = int((1800-1250)*2)
# size_test = int(0.1e3*2)
# size_validation = int(0.3e3*2)

# the data contains 3904 frames
# the data contains 2704 frames after excluding the first and last frames
# the test data contains 200 frames
# the data contains 2504 frames after excluding the first and last frames and test set
# the training data contains 1904 frames
# the validation data contains 600 frames


# mgo 333
# pbc set {9.2 9.2 8.88 90 90 90} -all ; pbc box
# pbc wrap -all
# files = ['mgo-1-cleaned.ener', 'mgo-charges-1-clean-cleaned.hirshfeld', 'mgo-pos-1-cleaned.xyz', 'mgo-frc-1-cleaned.xyz']
# box = np.array([[12.57, 0, 0, 0, 12.57, 0, 0, 0, 12.57]])
# num_atoms = 216

# directory = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/mgo/archer/mgo/cell-333/md/pbe-u-8/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-600k-csvr-timecon-1-COMVEL_TO-1e-10-nvt-hole-u-8'
# charged = True
# 0 to 1400 fs
# Discard first 100 frames and last 1400-1300*2 frames

# size_exclude_start = int(0.05e3*2)
# size_exclude_end = int((1400-1300)*2)
# size_test = int(0.1e3*2)
# size_validation = int(0.3e3*2)

# the data contains 2800 frames
# the data contains 2500 frames after excluding the first and last frames
# the test data contains 200 frames
# the data contains 2300 frames after excluding the first and last frames and test set
# the training data contains 1700 frames
# the validation data contains 600 frames


# tio2 336
# pbc set {13.8 13.8 17.76 90 90 90} -all ; pbc box
# pbc wrap -all
# files = ['tio2-1.ener', 'tio2-charges-1-clean.hirshfeld', 'tio2-pos-1.xyz', 'tio2-frc-1.xyz']
files = ['tio2-1-cleaned.ener', 'tio2-charges-1-clean-cleaned.hirshfeld', 'tio2-pos-1-cleaned.xyz', 'tio2-frc-1-cleaned.xyz']
box = np.array([[13.8, 0, 0, 0, 13.8, 0, 0, 0, 17.76]])
num_atoms = 324

folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-336/md'
directory = '{}/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-csvr-timecon-1-COMVEL_TO-1e-10-nvt-u-3.5-electron-dz2-u-3.0-rs-print-hse-25-schwarz-e-1e-4-f-1e-6-cfit11-cpfit3'.format(folder)
charged = True
# 0 to 480 fs
# Discard first 390 fs

size_exclude_start = int(390*2)
size_exclude_end = int(0*2)
size_test = int(30*2)
size_validation = int(50)

# the data contains 975 frames
# the data contains 195 frames after excluding the first and last frames
# the test data contains 60 frames
# the data contains 135 frames after excluding the first and last frames and test set
# the training data contains 85 frames
# the validation data contains 50 frames




# For energy just .flatten() potential energy column converted to eV Shape (200,1) -> (timesteps, 1)
print('Loading energy')
_, _, energy, _, _, time_val, time_per_step, step = read_energy(directory, files[0])
num_timesteps = np.shape(time_val)[0]
energy = np.reshape(energy, (num_timesteps, 1)) * param.hartree_to_ev
# print(energy[-1]/param.hartree_to_ev)
# energy_load = np.load('/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/archer/bulk/221_supercell/md/hole/cleaned-checked/400k-neutral/database/dpdata/database_ener_force_test/2/set.000/energy.npy')
# energy_load = np.load('/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liulab-hpc/bulk/400k-neutral-rcut-5-new/database_ener_force_test/2/set.000/energy.npy')
# energy_load = np.load('/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/archer/bulk/221_supercell/md/hole/cleaned-checked/400k-neutral/database_ener_force_test/2/set.000/energy.npy')
# print(energy_load[-1]/param.hartree_to_ev)

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

# Shape (200,360) -> (timesteps, num_atoms*3) where it is arranged [x1, y1, z1, x2, y2, z2 ... xn, yn, zn]
print('Loading coordinates')
coordinates, _, _, _, _, _, _ = load_values_coord(directory, files[2], ['Species', 'X', 'Y', 'Z'])
print(coordinates.shape)
coordinates = np.transpose(coordinates, axes=(0, 2, 1))
coordinates = coordinates.reshape(num_timesteps, num_atoms*3)
# np.save('/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/archer/bulk/221_supercell/md/hole/cleaned-checked/400k-neutral/database_ener_force_test/2/set.000/coord2.npy', coordinates)
# coordinates_load = np.load('/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/archer/bulk/221_supercell/md/hole/cleaned-checked/400k-neutral/database_ener_force_test/2/set.000/coord2.npy')
# print(coordinates_load.shape)
# print(coordinates_load[0, 0:3])
# print(coordinates[-1, 0:3])
# coordinates_load = np.load('/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/archer/bulk/221_supercell/md/hole/cleaned-checked/400k-neutral/database/dpdata/database_ener_force_test/2/set.000/coord.npy')
# coordinates_load = np.load('/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liulab-hpc/bulk/400k-neutral-rcut-5-new/database_ener_force_test/2/set.000/coord.npy')
# coordinates_load = np.load('/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/archer/bulk/221_supercell/md/hole/cleaned-checked/400k-neutral/database_ener_force_test/2/set.000/coord.npy')
# print(coordinates_load[-1, 0:3])

# Shape (200,360) -> (timesteps, num_atoms*3) where it is arranged [x1, y1, z1, x2, y2, z2 ... xn, yn, zn]
print('Loading forces')
forces, _, _, _, _, _, _ = load_values_coord(directory, files[3], ['Species', 'X', 'Y', 'Z'])
forces = forces * param.hartree_per_bohr_to_ev_per_angstrom
forces = np.transpose(forces, axes=(0, 2, 1))
forces = forces.reshape(num_timesteps, num_atoms*3)
# np.save('/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/archer/bulk/221_supercell/md/hole/cleaned-checked/400k-neutral/database_ener_force_test/2/set.000/force2.npy', forces)
# forces_load = np.load('/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/archer/bulk/221_supercell/md/hole/cleaned-checked/400k-neutral/database_ener_force_test/2/set.000/force2.npy')
# print(forces.shape)
# print(forces[-1, 0:3]/param.hartree_per_bohr_to_ev_per_angstrom)
# forces_load = np.load('/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liulab-hpc/bulk/400k-neutral-rcut-5-new/database_ener_force_test/2/set.000/force.npy')
# forces_load = np.load('/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/archer/bulk/221_supercell/md/hole/cleaned-checked/400k-neutral/database_ener_force_test/2/set.000/force.npy')
# print(forces_load[-1, 0:3]/param.hartree_per_bohr_to_ev_per_angstrom)

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
# index_training = np.loadtxt('{}/database_ener_force_train/index_training.raw'.format(directory)).astype(np.int64)
# index_validation = np.loadtxt('{}/database_ener_force_test/1/index_validation.raw'.format(directory)).astype(np.int64)

# Set random index
rng = np.random.default_rng()
index_validation = rng.choice(energy.shape[0], size=size_validation, replace=False)
index_validation = np.sort(index_validation)
index_training = list(set(range(energy.shape[0])) - set(index_validation))

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
population_train = np.zeros((np.shape(population_alpha_train)[0], num_atoms, 2))
for i in range(np.shape(population_alpha_train)[0]):
    for j in range(num_atoms):
        population_train[i, j, 0] = population_alpha_train[i, j]
        population_train[i, j, 1] = population_beta_train[i, j]
population_valid = np.zeros((np.shape(population_alpha_valid)[0], num_atoms, 2))
for i in range(np.shape(population_alpha_valid)[0]):
    for j in range(num_atoms):
        population_valid[i, j, 0] = population_alpha_valid[i, j]
        population_valid[i, j, 1] = population_beta_valid[i, j]
population_test = np.zeros((np.shape(population_alpha_test)[0], num_atoms, 2))
for i in range(np.shape(population_alpha_test)[0]):
    for j in range(num_atoms):
        population_test[i, j, 0] = population_alpha_test[i, j]
        population_test[i, j, 1] = population_beta_test[i, j]
np.save('{}/database_population_train/set.000/atomic_spin.npy'.format(directory), population_train)
np.save('{}/database_population_test/1/set.000/atomic_spin.npy'.format(directory), population_valid)
np.save('{}/database_population_test/2/set.000/atomic_spin.npy'.format(directory), population_test)

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
