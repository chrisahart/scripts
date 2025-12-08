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

"""
    Plot energy and forces for bulk hematite
"""

params = {'axes.formatter.limits': [-4, 4],
          'axes.labelsize': '14',
          'axes.titlesize': '16',
          'lines.markersize': '8',
          }
plt.rcParams.update(params)


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


def read_mulliken(folder, filename, num_atoms):
    """
    Read Mulliken analysis from CP2K output file
    """

    # Read number of atoms and labels from .xyz file
    files = ['{}/{}'.format(folder, filename)]
    cols = ['Atom', 'Element', 'Kind', 'Pop_1', 'Pop_2', 'Charge', 'Spin', 'A', 'B', 'C', 'D']
    file_spec1 = pd.read_csv(files[0], names=cols, delim_whitespace=True, skiprows=5)

    file_spec1 = file_spec1.drop(columns=['A'])
    file_spec1 = file_spec1.drop(columns=['B'])
    file_spec1 = file_spec1.drop(columns=['C'])
    file_spec1 = file_spec1.drop(columns=['D'])

    species = file_spec1['Element']
    file_spec1 = file_spec1.drop(columns=['Element'])

    file_spec1 = file_spec1.apply(pd.to_numeric, errors='coerce')
    file_spec1 = file_spec1.dropna()
    # file_spec1 = file_spec1.dropna(axis='rows', thresh=0)
    # file_spec1 = file_spec1.dropna(axis='columns', thresh=1)
    file_spec1 = file_spec1.reset_index(drop=True)
    cols_new = list(file_spec1.columns)

    # Loop over each timestep and atoms
    num_timesteps = int(file_spec1.shape[0]/num_atoms)
    mulliken_data = np.zeros((num_timesteps, len(cols_new), num_atoms))

    for timestep in range(num_timesteps):
        for atom in range(num_atoms):
            for i in range(len(cols_new)):
                mulliken_data[timestep, i, atom] = file_spec1[cols_new[i]].values[atom + timestep * num_atoms]

    return file_spec1, mulliken_data, species


def read_hubbard(folder, filename, num_atoms):
    """
    Read Hubbard analysis from CP2K output file
    """

    # Read number of atoms and labels from .xyz file
    files = ['{}/{}'.format(folder, filename)]
    cols = ['Atom', 'Shell', 'py', 'pz', 'px', 'Trace', 'A', 'B', 'C', 'D', 'E', 'F', 'G']
    file_spec1 = pd.read_csv(files[0], names=cols, delim_whitespace=True, skiprows=2)
    print(file_spec1)
    file_spec1 = file_spec1.drop(columns=['A'])
    file_spec1 = file_spec1.drop(columns=['B'])
    file_spec1 = file_spec1.drop(columns=['C'])
    file_spec1 = file_spec1.drop(columns=['D'])
    file_spec1 = file_spec1.drop(columns=['E'])
    file_spec1 = file_spec1.drop(columns=['F'])
    file_spec1 = file_spec1.drop(columns=['G'])
    print(file_spec1)
    file_spec1 = file_spec1.apply(pd.to_numeric, errors='coerce')
    file_spec1 = file_spec1.dropna(axis='rows', thresh=2)
    file_spec1 = file_spec1.dropna(axis='columns', thresh=1)
    file_spec1 = file_spec1.reset_index(drop=True)
    cols_new = list(file_spec1.columns)

    # Loop over each timestep and atoms
    num_timesteps = int(file_spec1.shape[0] / num_atoms)
    hubbard_data = np.zeros((num_timesteps, len(cols_new), num_atoms))

    for timestep in range(num_timesteps):
        for atom in range(num_atoms):
            for i in range(len(cols_new)):
                hubbard_data[timestep, i, atom] = file_spec1[cols_new[i]].values[atom + timestep * num_atoms]

    return file_spec1, hubbard_data

offset = 0
xlim_auto = False

# Anatase
# folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/anatase/archer/anatase/cell-221/md'
# # folder_1 = '{}/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-dz-400k-500k-csvr-timecon-1-COMVEL_TO-1e-10-nvt-u-3.5-electron-timestep-1.0'.format(folder)  # 16 ps: Polaron hops 9 times, some with intermediate delocalisation
# # folder_1 = '{}/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-dz-400k-500k-csvr-timecon-1-COMVEL_TO-1e-10-nvt-u-3.5-electron-timestep-1.0-o-3.5'.format(folder)  # 16 ps: Similar to above, 10 hops, some with intermediate delocalisation
# # folder_1 = '{}/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-dz-400k-500k-csvr-timecon-1-COMVEL_TO-1e-10-nvt-u-3.5-hole-timestep-1.0'.format(folder)  # 1 ps: No polaron, fully delocalised
# folder_1 = '{}/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-dz-400k-500k-csvr-timecon-1-COMVEL_TO-1e-10-nvt-u-3.5-hole-timestep-1.0-o-3.5'.format(folder)  # 16 ps: Dynamics are good, polaron is stable whole time and does not delocalise ***
# folder_1 = '{}/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-dz-400k-500k-csvr-timecon-1-COMVEL_TO-1e-10-nvt-u-3.5-hole-timestep-1.0-o-3.5-ti-0'.format(folder)  # 14 ps: Dynamics seem worse than above, more likely for hop to fail.
# folder_1 = '{}/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-dz-400k-500k-csvr-timecon-1-COMVEL_TO-1e-10-nvt-u-3.5-hole-timestep-1.0-o-6-ti-0'.format(folder)  # 10 ps: bad dynamics
# files = ['tio2-1.ener', 'tio2-charges-1-clean.hirshfeld', 'tio2-pos-1.xyz', 'tio2-frc-1.xyz']
# num_atoms = 48
# box_size = [7.56, 7.56, 9.62, 90, 90, 90]

# folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/anatase/archer/anatase/cell-331/md'
# # folder_1 = '{}/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-dz-400k-500k-csvr-timecon-1-COMVEL_TO-1e-10-nvt-u-3.5-electron-timestep-1.0'.format(folder)  # 10 ps: No polaron hopping
# # folder_1 = '{}/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-dz-400k-500k-csvr-timecon-1-COMVEL_TO-1e-10-nvt-u-3.5-electron-timestep-1.0-o-3.5'.format(folder)  # 400 fs: No polaron hopping 400 fs, then tries to hop and fails SCF
# # folder_1 = '{}/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-dz-400k-500k-csvr-timecon-1-COMVEL_TO-1e-10-nvt-u-3.5-hole-timestep-1.0'.format(folder)  # 600 fs: No polaron
# folder_1 = '{}/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-dz-400k-500k-csvr-timecon-1-COMVEL_TO-1e-10-nvt-u-3.5-hole-timestep-1.0-o-3.5'.format(folder)  # 10 ps: 5 polaron hops, some discontinous. Not as good as 221. Larger cell polaron less mobile and not as well behaved?
# folder_1 = '{}/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-dz-400k-500k-csvr-timecon-1-COMVEL_TO-1e-10-nvt-u-3.5-hole-timestep-1.0-o-3.5-ti-0'.format(folder)  # 7 ps: no polaron hopping
# # folder_1 = '{}/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-dz-400k-500k-csvr-timecon-1-COMVEL_TO-1e-10-nvt-u-3.5-hole-timestep-1.0-o-6-ti-0'.format(folder)  # 5 ps: Bad dynamics
# # files = ['tio2-1.ener', 'tio2-charges-1-clean.hirshfeld', 'tio2-pos-1.xyz', 'tio2-frc-1.xyz']
# num_atoms = 108
# box_size = [11.34, 11.34, 9.62, 90, 90, 90]

# folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/anatase/archer/anatase/cell-441/md'
# folder_1 = '{}/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-csvr-timecon-1-COMVEL_TO-1e-10-nvt-dz-electron-timestep-1.0-o-0-ti-3.5'.format(folder)  # 3.5 ps: No polaron hopping. dt = 0.5
# xlim_1 = [0, 7101]
# folder_1 = '{}/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-csvr-timecon-1-COMVEL_TO-1e-10-nvt-dz-hole-timestep-1.0-o-3.5-ti-0'.format(folder)  # 3.5 ps: 2 hops, good. dt = 0.5
# xlim_1 = [0, 7120]
# temperature_set = 500

# folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/anatase/archer/anatase/cell-441/md-cell-opt-hse-20'
# folder = '/Volumes/Samsung/Data/Postdoc2/Data/Work/calculations/tio2/anatase/archer/anatase/cell-441/md-cell-opt-hse-20'
# folder_1 = '{}/hole-u-o-5-ti-0-600k'.format(folder)  # 10 ps lifetime fs 2500
# xlim_1 = [0, 10906]
# xlim_1 = [0, 10000]
# folder_1 = '{}/hole-u-o-5-ti-0-600k-hse-25'.format(folder)  # CSVR no hoppng
# xlim_1 = [0, 2192]
# folder_1 = '{}/hole-u-o-5-ti-0-600k-hse-25-center'.format(folder)  # CSVR no hoppng
# xlim_1 = [0, 1457]
# folder_1 = '{}/hole-u-o-5-ti-0-600k-hse-25-center-rs-19'.format(folder)  # lifetime fs 355
# xlim_1 = [0, 2922]
# folder_1 = '{}/hole-u-o-5-ti-0-600k-hse-25-center-rs-22'.format(folder)  # CSVR 3 ps total
# xlim_1 = [0, 1504]
# folder_1 = '{}/hole-u-o-5-ti-0-600k-hse-25-center-rs-22-yungu-rs'.format(folder)  # NOSE 10 ps total
# xlim_1 = [0, 10013-1504]
# folder_1 = '{}/hole-u-o-5-ti-0-600k-hse-25-center-rs-22-yungu-rs-rs-19'.format(folder)  # NOSE 10 ps total
# xlim_1 = [0, 10673]
# folder_1 = '{}/hole-u-o-5-ti-0-600k-hse-25-center-rs-22-yungu-rs-rs-19-rs'.format(folder)  # NOSE 10 ps total
# xlim_1 = [0, 5000]
# xlim_1 = [10674, 10674+5000]
# xlim_1 = [0, 11764]
# xlim_1 = [0, 5000]

# folder = '/Volumes/Samsung/Data/Postdoc2/Data/Work/calculations/tio2/anatase/archer/anatase/cell-441/md-cell-opt-hse-20/hse-19-complete'
# folder_1 = '{}/combined'.format(folder)
# folder_1 = '{}/trajectory_mdanalysis'.format(folder)
# xlim_1 = [0, 20000]
# xlim_1 = [14200, 14270]
# xlim_1 = [0, 10000]
# xlim_1 = [10000, 20000]
# offset = 10000
# offset = 0
# files = ['tio2-1-cleaned.ener', 'tio2-charges-1-clean-cleaned.hirshfeld', 'tio2-pos-1-cleaned.xyz', 'tio2-frc-1-cleaned.xyz']

folder = '/Volumes/Samsung/Data/Postdoc2/Data/Work/calculations/tio2/anatase/archer/anatase/cell-441/reftraj/trajectory_mdanalysis/md'
folder_1 = '{}/combined'.format(folder)
files = ['tio2-1-cleaned.ener', 'tio2-charges-1-clean-cleaned.hirshfeld', 'tio2-pos-1-cleaned.xyz', 'tio2-frc-1-cleaned.xyz']
xlim_auto = True

# folder = '/Volumes/Samsung/Data/Postdoc2/Data/Work/calculations/tio2/anatase/archer/anatase/cell-441/reftraj/trajectory_mdanalysis/ts/good'
# folder_1 = '{}/combined'.format(folder)
# xlim_auto = True
# files = ['tio2-1.ener', 'tio2-charges-1-clean.hirshfeld', 'tio2-pos-1.xyz', 'tio2-frc-1.xyz']

temperature_set = 600
num_atoms = 192
# box_size = [15.12, 15.12, 9.62, 90, 90, 90]
box_size = [15.08, 15.08, 9.68, 90, 90, 90]

# folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/anatase/archer/anatase/cell-442/md'
# # folder_1 = '{}/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-csvr-timecon-1-COMVEL_TO-1e-10-nvt-dz-electron-timestep-1.0-o-0-ti-3.5'.format(folder)  # 1.4 ps:
# folder_1 = '{}/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-csvr-timecon-1-COMVEL_TO-1e-10-nvt-dz-hole-timestep-1.0-o-3.5-ti-0'.format(folder)  # 1.4 ps: No hopping
# files = ['tio2-1.ener', 'tio2-charges-1-clean.hirshfeld', 'tio2-charges-1-clean.mulliken', 'tio2-charges-1-clean.hubbard', 'tio2-pos-1.xyz', 'tio2-frc-1.xyz']
# num_atoms = 384
# box_size = [15.12, 15.12, 19.23, 90, 90, 90]

# folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/anatase/archer/anatase/cell-442/md-cell-opt'
# folder_1 = '{}/hole-u-o-3.5-ti-0-500k'.format(folder)  # 3 ps : No polaron hopping
# folder_1 = '{}/hole-u-o-3.5-ti-0-500k-rs-300k'.format(folder)  # 3 ps : No polaron hopping
# folder_1 = '{}/hole-u-o-3.5-ti-0-500k-rs-800k'.format(folder)  # 2.5 ps: Polaron hopping not stable, restart at U = 5?
# folder_1 = '{}/hole-u-o-3.5-ti-0-500k-rs-1200k'.format(folder)  # 2.5 ps: Polaron hopping not stable, restart at U = 5?

# folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/anatase/archer/anatase/cell-442/md-cell-opt-hse-22'
# folder_1 = '{}/hole-u-o-5-ti-0-500k-atom-282'.format(folder)  # 1.75 ps: 1 bad hop at 250 fs
# folder_1 = '{}/hole-u-o-5-ti-0-500k-atom-282-rs-step-200'.format(folder)  # 11 ps : 4 polaron hops, not adiabatic but not smooth
# xlim_1 = [0, 11413]

# folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/anatase/archer/anatase/cell-442/md-cell-opt-hse-20'
# folder = '/Volumes/Samsung/Data/Postdoc2/Data/Work/calculations/tio2/anatase/archer/anatase/cell-442/md-cell-opt-hse-20'
# folder = '/Volumes/Elements/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/anatase/archer/anatase/cell-442/md-cell-opt-hse-20'
# folder_1 = '{}/hole-u-o-5-ti-0-300k'.format(folder)  # 10 ps : No polaron hopping
# folder_1 = '{}/hole-u-o-5-ti-0-300k-rs-hse-25'.format(folder)  # 2 ps : No polaron hopping
# folder_1 = '{}/hole-u-o-5-ti-0-800k'.format(folder)  # 2 ps : 4 polaron hops, bad and unstable. CANCEL
# folder_1 = '{}/hole-u-o-8-ti-0-800k'.format(folder)  # SCF run NOT converged CANCEL
# folder_1 = '{}/hole-u-o-10-ti-0-800k'.format(folder)  # SCF run NOT converged CANCEL
# folder_1 = '{}/hole-u-o-5-ti-0-1200k'.format(folder)  # 2 ps : 5 polaron hops, bad and unstable. CANCEL
# folder_1 = '{}/hole-u-o-8-ti-0-1200k'.format(folder)  # SCF run NOT converged CANCEL
# folder_1 = '{}/hole-u-o-10-ti-0-1200k'.format(folder)  # SCF run NOT converged CANCEL
# folder_1 = '{}/hole-u-o-5-ti-0-600k'.format(folder)  # CSVR 4 hops lifetime fs 2500
# xlim_1 = [0, 10000]
# xlim_1 = [0, 10710]
# folder_1 = '{}/hole-u-o-5-ti-0-600k-rs-hse-25'.format(folder)  # CSVR 1 hop lifetime fs 3126
# xlim_1 = [0, 3126]
# folder_1 = '{}/hole-u-o-5-ti-0-600k-rs-hse-25-rs-19'.format(folder)  # lifetime fs 295
# xlim_1 = [0, 5020]
# folder_1 = '{}/hole-u-o-5-ti-0-600k-rs-hse-25-rs-22'.format(folder)  # lifetime fs inf
# xlim_1 = [0, 1671]
# files = ['tio2-1.ener', 'tio2-charges-1-clean.hirshfeld', 'tio2-charges-1-clean.mulliken', 'tio2-charges-1-clean.hubbard', 'tio2-pos-1.xyz', 'tio2-frc-1.xyz']
# num_atoms = 384
# box_size = [15.12, 15.12, 19.23, 90, 90, 90]
# temperature_set = 600

print('folder_1', folder_1)
file_energy_1, energy_kinetic_1, energy_potential_1, energy_total_1, temperature_1, time_val_1, time_per_step_1, step_1 = read_energy(folder_1, files[0])
time_val_energy =  time_val_1
hirshfeld_1_df, hirshfeld_1_np, _ = read_hirsh(folder_1, files[1], num_atoms)
print(hirshfeld_1_df)
plot_mulliken = True
plot_mulliken = False
if plot_mulliken:
    mulliken_1_df, mulliken_1_np, _ = read_mulliken(folder_1, files[2], num_atoms)
    print(mulliken_1_df)
# plot_hubbard = True
plot_hubbard = False
calc_distance = False
# calc_distance = True
save_fig = False
save_fig = True
# calc_distance = False
if plot_hubbard:
    atoms_hubbard = num_atoms/3 * 2
    hubbard_1_df, hubbard_1_np = read_hubbard(folder_1, files[3], atoms_hubbard*2)
    print(hubbard_1_df)
    print(hubbard_1_np.shape)
topology_file = '{}/system.xyz'.format(folder_1)
trajectory_file = '{}/{}'.format(folder_1, files[2])

local_bonds = 3
# timestep = 0.5
timestep = 1
frames_skip = 100
frames_skip = 0

num_timesteps = np.shape(hirshfeld_1_np)[0]
time_array = np.linspace(start=0, stop=num_timesteps*timestep,num=num_timesteps)
num_timesteps2 = np.shape(time_per_step_1)[0]
time_array2 = np.linspace(start=0, stop=num_timesteps*timestep, num=num_timesteps2)
print('num_timesteps np.shape(hirshfeld_1_np)[0]', num_timesteps)
print('num_timesteps2 np.shape(time_per_step_1)[0]', num_timesteps2)
if xlim_auto: xlim_1 = [0, np.shape(hirshfeld_1_np)[0]]

if plot_mulliken:
    num_timesteps3 = np.shape(mulliken_1_np)[0]
    time_array3 = np.linspace(start=0, stop=num_timesteps3*timestep,num=num_timesteps3)
if plot_hubbard:
    num_timesteps4 = np.shape(hubbard_1_np)[0]
    time_array4 = np.linspace(start=0, stop=num_timesteps4*timestep,num=num_timesteps4)

# xlim_1 = [0, 10e3]
# xlim_1 = [0, 60]
# xlim_1 = [0, 1e3]
# xlim_1 = [0, 1e4]
# xlim_1 = [0, time_array[-1]]
ylim_1 = [0, 0.9]
draw_legend = False
folder_save = folder_1

print('Polaron:', np.max(hirshfeld_1_np[-1, 5, :]), np.argmax(hirshfeld_1_np[-1, 5, :]))
polaron_atom = np.argmax(hirshfeld_1_np[-1, 5, :])

# Plot time taken
# fig_time_md, ax_time_md = plt.subplots()
ylim_1_time = [10, 500]
fig_time_md, ax_time_md = plt.subplots(figsize=(10, 4))
ax_time_md.plot(time_array2, time_per_step_1, 'k-')
ax_time_md.set_xlabel('Time / fs')
ax_time_md.set_ylabel('Time per MD step')
ax_time_md.set_xlim(xlim_1)
ax_time_md.set_ylim(ylim_1_time)
fig_time_md.tight_layout()
if save_fig: fig_time_md.savefig('{}/time_taken_md_step.png'.format(folder_save), dpi=300)

# Plot energy
# fig_time_md, ax_time_md = plt.subplots()
ylim_1_time = [10, 500]
fig_energy_potential, ax_energy_potential = plt.subplots(figsize=(10, 4))
ax_energy_potential.plot(time_array2, energy_potential_1, 'k-')
ax_energy_potential.set_xlabel('Time / fs')
ax_energy_potential.set_ylabel('Energy / au')
ax_energy_potential.set_xlim(xlim_1)
fig_energy_potential.tight_layout()
if save_fig: fig_energy_potential.savefig('{}/energy_potential.png'.format(folder_save), dpi=300)

# Plot Hirshfeld spin of all atoms
fig_spin1, ax_spin1 = plt.subplots(figsize=(10, 4))
# fig_spin1, ax_spin1 = plt.subplots()
temp = np.zeros(num_timesteps)
for j in range(num_atoms):
    # ax_spin1.plot(time_array-1, hirshfeld_1_np[:, 5, j], '.-', label='{}'.format(j+1))
    ax_spin1.plot(time_array, hirshfeld_1_np[:, 5, j], '-', label='{}'.format(j+1))
# ax_spin1.plot(time_array, hirshfeld_1_np[:, 5, polaron_atom], 'k-',)
if draw_legend: ax_spin1.legend(frameon=True)
# ax_spin1.set_xlabel('Timestep')
ax_spin1.set_xlabel('Time / fs')
ax_spin1.set_ylabel('Spin moment')
ax_spin1.set_xlim(xlim_1)
ax_spin1.set_ylim(ylim_1)
fig_spin1.tight_layout()
if save_fig: fig_spin1.savefig('{}/hirshfeld_spin_all.png'.format(folder_save), dpi=300)

# Plot Hirshfeld spin of all atoms
# fig_spin1_square, ax_spin1_square = plt.subplots(figsize=(7, 5))
fig_spin1_square, ax_spin1_square = plt.subplots()
temp = np.zeros(num_timesteps)
for j in range(num_atoms):
    # ax_spin1_square.plot(time_array-1, hirshfeld_1_np[:, 5, j], '.-', label='{}'.format(j+1))
    # ax_spin1_square.plot(time_array, hirshfeld_1_np[:, 5, j], '.-', label='{}'.format(j+1))
    ax_spin1_square.plot(time_array, hirshfeld_1_np[:, 5, j], '-', label='{}'.format(j+1))
# ax_spin1_square.plot(time_array, hirshfeld_1_np[:, 5, polaron_atom], 'k-',)
if draw_legend: ax_spin1_square.legend(frameon=True)
# ax_spin1_square.set_xlabel('Timestep')
ax_spin1_square.set_xlabel('Time / fs')
ax_spin1_square.set_ylabel('Spin moment')
ax_spin1_square.set_xlim(xlim_1)
ax_spin1_square.set_ylim(ylim_1)
fig_spin1_square.tight_layout()
if save_fig: fig_spin1_square.savefig('{}/hirshfeld_spin_all2.png'.format(folder_save), dpi=300)

# Plot Mulliken spin of all atoms
if plot_mulliken:
    ylim_1 = [0, 1.1]
    # fig_spin2, ax_spin2 = plt.subplots()
    fig_spin2, ax_spin2 = plt.subplots(figsize=(10, 4))
    temp = np.zeros(num_timesteps3)
    for j in range(num_atoms):
        ax_spin2.plot(time_array3, mulliken_1_np[:, 5, j], '-', label='{}'.format(j+1))
    if draw_legend: ax_spin2.legend(frameon=True)
    # ax_spin2.set_xlabel('Timestep')
    ax_spin2.set_xlabel('Time / fs')
    ax_spin2.set_ylabel('Spin moment')
    ax_spin2.set_xlim(xlim_1)
    ax_spin2.set_ylim(ylim_1)
    fig_spin2.tight_layout()
    if save_fig: fig_spin2.savefig('{}/mulliken_spin_all.png'.format(folder_save), dpi=300)

# Plot Hubbard spin of all atoms
if plot_hubbard:
    ylim_1 = [0, 1.4]
    fig_spin3, ax_spin3 = plt.subplots(figsize=(10, 4))
    # fig_spin3, ax_spin3 = plt.subplots()
    temp = np.zeros(num_timesteps3)
    for j in range(atoms_hubbard):
        ax_spin3.plot(time_array4, hubbard_1_np[:, 5, j]-hubbard_1_np[:, 5, j+atoms_hubbard], '-', label='{}'.format(j+1))
    if draw_legend: ax_spin3.legend(frameon=True)
    # ax_spin3.set_xlabel('Timestep')
    ax_spin3.set_xlabel('Time / fs')
    ax_spin3.set_ylabel('Spin moment')
    ax_spin3.set_xlim(xlim_1)
    ax_spin3.set_ylim(ylim_1)
    fig_spin3.tight_layout()
    if save_fig: fig_spin3.savefig('{}/hubbard_spin_all.png'.format(folder_save), dpi=300)

# Setup md analysis environment
if calc_distance:
    universe = mda.Universe(topology_file, trajectory_file)
    universe.dimensions = box_size
    num_timesteps2 = len(universe.trajectory)
    time_val_1 = np.linspace(start=0, stop=len(universe.trajectory) * timestep, num=len(universe.trajectory))
    atoms_ti = universe.select_atoms('name Ti')
    num_atoms_ti = len(atoms_ti)
    atoms_o = universe.select_atoms('name O')
    num_atoms_o = len(atoms_o)
    dist_arr = distances.distance_array(atoms_o.positions, atoms_ti.positions, box=box_size)
    # dist_arr = distances.distance_array(atoms_ti.positions, atoms_o.positions, box=box_size)
    bond_lengths_time = np.zeros((num_timesteps2, num_atoms_o, num_atoms_ti))
    for ts in universe.trajectory:
        frame = universe.trajectory.frame
        bond_lengths_time[frame] = distances.distance_array(atoms_o.positions, atoms_ti.positions, box=box_size)
    bond_lengths_time_sorted = np.zeros((num_timesteps2, num_atoms_o, local_bonds))
    bond_lengths_time_sorted_mean = np.zeros((num_timesteps2, num_atoms_o))
    for i in range(num_atoms_o):
        for j in range(num_timesteps2):
            bond_lengths_time_sorted[j, i] = np.sort(bond_lengths_time[j, i])[0:local_bonds]
            # bond_lengths_time_sorted_mean[j, i] = np.mean(bond_lengths_time_sorted[j, i])
    bond_lengths_time_sorted_mean = np.mean(bond_lengths_time_sorted, axis=2)

# Plot average of 3 O-Ti bonds
if calc_distance:
    metric = np.zeros((num_atoms_o, num_timesteps2))
    fig_bonds_1, ax_bonds_1 = plt.subplots(figsize=(10, 4))
    # fig_bonds_1, ax_bonds_1 = plt.subplots()
    for i in range(num_atoms_o):
        # ax_bonds_1.plot(time_val_1 - time_val_1[0], np.sum(bond_lengths_time_sorted, axis=2)[:, i], '-', label='Fe {}'.format(i + 1))
        ax_bonds_1.plot(time_val_1 - time_val_1[0], bond_lengths_time_sorted_mean[:, i], '-', label='Fe {}'.format(i + 1))
    # ax_bonds_1.plot(time_val_1 - time_val_1[0], bond_lengths_time_sorted_mean[:, polaron_atom], 'k-', label='Fe {}'.format(polaron_atom + 1))
    ax_bonds_1.set_xlabel('Time / fs')
    # ax_bonds_1.set_xlabel('Timestep')
    ax_bonds_1.set_ylabel('Average of {} O-Ti bond lengths / A'.format(local_bonds))
    # # ax_bonds_1.set_xlim([0, len(universe.trajectory)])
    ax_bonds_1.set_xlim(xlim_1)
    # ax_bonds_1.set_xlim([0, len(universe.trajectory) * timestep])
    # ax_bonds_1.set_ylim([0.06, -0.10])
    if save_fig: fig_bonds_1.savefig('{}/bond_lengths_average.png'.format(folder_save), dpi=300)
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
# if save_fig: fig_bonds_2.savefig('{}/bond_lengths_average.png'.format(folder_save), dpi=300)
# fig_bonds_2.tight_layout()

# Calculate polaron atom
hirshfeld_mobility = hirshfeld_1_np[int(xlim_1[0]):int(xlim_1[1])]
num_timesteps_mobility = int(xlim_1[1]-xlim_1[0])
time_val_1 = time_val_energy[int(xlim_1[0]):int(xlim_1[1])]
# time_val_1 = np.linspace(start=xlim_1[0], stop=int(xlim_1[1]), num=num_timesteps_mobility)

if calc_distance:
    polaron_atom_time = np.zeros(num_timesteps, dtype=int)
    for j in range(num_timesteps):
        polaron_atom_time[j] = int(np.argmax(hirshfeld_1_np[j, 5, :]))
    # polaron_atoms = np.unique(polaron_atom_time)
    polaron_atoms = polaron_atom_time[np.insert(polaron_atom_time[:-1] != polaron_atom_time[1:], 0, True)]
    print('polaron_atoms', polaron_atoms+1)

    # Calculate distance between current timestep polaron atom and next timestep
    # Then get all non-zero answer
    polaron_distances = np.zeros(num_timesteps)
    for j in range(num_timesteps - 1):
        polaron_distances[j] = distances.distance_array(universe.select_atoms('index {}'.format(polaron_atom_time[j])).positions,
                                                        universe.select_atoms('index {}'.format(polaron_atom_time[j+1])).positions,
                                                        box=box_size)
    polaron_distances = polaron_distances[xlim_1[0]:int(xlim_1[1])]
    polaron_distances_hop = polaron_distances[np.nonzero(polaron_distances)]
    polaron_indices = np.nonzero(polaron_distances)[0]
    # print('polaron_distance', polaron_distances)
    print('polaron_distances_hop', polaron_distances_hop)
    print('np.shape(polaron_distances_hop)[0]', np.shape(polaron_distances_hop)[0])
    print('polaron hop index', polaron_indices)
    print('polaron hop time', time_val_energy[polaron_indices])

    # Plot polaron distances
    metric = np.zeros((num_atoms_o, num_timesteps_mobility))
    # fig_bonds_2, ax_bonds_2 = plt.subplots()
    fig_bonds_2, ax_bonds_2 = plt.subplots(figsize=(10, 2))
    # ax_bonds_2.plot(time_val_1 - time_val_1[0] - offset, polaron_distances[:-1], 'kx-')
    ax_bonds_2.plot(time_val_1[:int(xlim_1[1])] - time_val_1[0] - offset, polaron_distances, 'kx-')
    ax_bonds_2.set_xlabel('Time / fs')
    # ax_bonds_2.set_xlabel('Timestep')
    ax_bonds_2.set_ylabel('Distance / A')
    # if draw_legend: ax_bonds_2.legend(frameon=True)
    # # ax_bonds_2.set_xlim([0, len(universe.trajectory)])
    ax_bonds_2.set_xlim(xlim_1)
    # ax_bonds_2.set_xlim([0, len(universe.trajectory) * timestep])
    # ax_bonds_2.set_ylim([0.06, -0.10])
    if save_fig: fig_bonds_2.savefig('{}/polaron_hopping_distance.png'.format(folder_save), dpi=300)
    fig_bonds_2.tight_layout()

    # Calculate mobility using xlim_1
    # hops_distance = np.array([2.99547136, 2.85545014, 3.01149688]) * 1e-8  # Angstrom to cm
    # hops_time = (5103+200)/2 * 1e-12  # fs to s
    hops_distance = polaron_distances_hop * 1e-8  # Angstrom to cm
    hops_time = (xlim_1[1] - xlim_1[0]) * 1e-15  # ps 1e-12 fs 1e-15
    print('hops per ps ', np.shape(hops_distance)[0]/hops_time*1e-15*1e3)

    rate_constant = np.shape(hops_distance)[0] / hops_time
    print('rate_constant', rate_constant)
    print('rate_constant / 1e12', rate_constant/1e12)
    if np.shape(hops_distance)[0] > 2: print('lifetime fs', 1/rate_constant * 1e15)

    mean_distance = np.mean(hops_distance)
    print('mean_distance', mean_distance)

    site_multiplicity = 1
    diffusion_constant_analytical = (np.mean(hops_distance) ** 2 * site_multiplicity * rate_constant) / 2
    mobility = (1.60217662e-19 * diffusion_constant_analytical) / (1.380649e-23 * temperature_set)
    print('mobility analytical', mobility, 'cm^2/(V·s)')

    # print('lifetime hematite', 1/1.2e12 * 1e15)
    # test = 1.60217662e-19 * (3e-8**2*3*1.2e12)/2 / (1.380649e-23 * 600)
    # print('mobility test hematite', test)

    # mobility = 1.60217662e-19 * (3e-8**2*3*rate_constant)/2 / (1.380649e-23 * 600)
    # print('mobility analytical 2', mobility, 'cm^2/(V·s)')

    mean_square_displacement = np.sum(hops_distance ** 2) / hops_time
    diffusion_constant_numerical = mean_square_displacement / 2
    mobility = (1.60217662e-19 * diffusion_constant_numerical) / (1.380649e-23 * temperature_set)
    print('mobility numerical', mobility, 'cm^2/(V·s)')

    k_el = 1
    temp = temperature_set
    kb_t_au = 8.617333262145E-5 * temperature_set  # KbT in eV
    kb_t = 1.38e-23 * temperature_set  # KbT in SI units
    # vn = 2.4e13  # 0.10 eV to s^-1 Deskins Dupuis TiO2 rutile (optic-mode phonon frequencies)
    vn = 2.66e13  # 0.11 eV to s^-1 Deskins Dupuis TiO2 anatase (optic-mode phonon frequencies)
    # rate_constant = vn * k_el * np.exp(-energy/kb_t)
    activation_energy = -np.log(rate_constant / (vn * k_el)) * kb_t_au
    print('activation_energy / meV', activation_energy*1e3)

    # hirshfeld and distance subplot
    rows, cols = 2, 1
    # fig_plot_all, ax_plot_all = plt.subplots(rows, cols,sharex='col', sharey='row',
    #                             figsize=(10, 6), gridspec_kw={'height_ratios': [2, 1],  'hspace': 0.05})
    fig_plot_all, ax_plot_all = plt.subplots(rows, cols,sharex='col', sharey='row',
                                figsize=(18, 6), gridspec_kw={'height_ratios': [2, 1],  'hspace': 0.05})
    temp = np.zeros(num_timesteps)
    for j in range(num_atoms):
        ax_plot_all[0].plot((time_val_1 - offset)/1e3, hirshfeld_mobility[:int(xlim_1[1]), 5, j], '-', label='{}'.format(j + 1))
    if draw_legend: ax_plot_all[0].legend(frameon=True)
    # ax_plot_all[0].set_xlabel('Time / fs')
    ax_plot_all[0].set_ylabel('Spin moment')
    ax_plot_all[0].set_xlim((np.array(xlim_1)-offset)/1000)
    # ax_plot_all[0].set_ylim(ylim_1)
    ax_plot_all[0].set_ylim([0, 0.9])
    ax_plot_all[1].plot((time_val_1 - offset)/1e3, polaron_distances, 'kx-')
    ax_plot_all[1].set_xlabel('Time / ps')
    ax_plot_all[1].set_ylabel('Distance / A')
    ax_plot_all[1].set_xlim((np.array(xlim_1)-offset)/1000)
    ax_plot_all[1].set_ylim([0, 3.3])
    fig_plot_all.tight_layout()
    fig_plot_all.subplots_adjust(hspace=0.05)
    if save_fig: fig_plot_all.savefig('{}/polaron_subplot.png'.format(folder_save), dpi=300)
    
    # hirshfeld and energy subplot
    # fig_hirshfeld_energy, ax_hirshfeld_energy = plt.subplots(rows, cols,sharex='col', sharey='row',
    #                             figsize=(18, 6), gridspec_kw={'height_ratios': [1, 1],  'hspace': 0.05})
    # temp = np.zeros(num_timesteps)
    # for j in range(num_atoms):
    #     ax_hirshfeld_energy[0].plot((time_val_1 - offset)/1e3, hirshfeld_mobility[:int(xlim_1[1]), 5, j], '-', label='{}'.format(j + 1))
    # if draw_legend: ax_hirshfeld_energy[0].legend(frameon=True)
    # # ax_hirshfeld_energy[0].set_xlabel('Time / fs')
    # ax_hirshfeld_energy[0].set_ylabel('Spin moment')
    # ax_hirshfeld_energy[0].set_xlim((np.array(xlim_1)-offset)/1000)
    # # ax_hirshfeld_energy[0].set_ylim(ylim_1)
    # ax_hirshfeld_energy[0].set_ylim([0, 0.8])
    # potential_energy_plot = energy_potential_1[int(xlim_1[0]):int(xlim_1[1])]-np.mean(energy_potential_1[int(xlim_1[0]):int(xlim_1[1])])
    # potential_energy_plot = potential_energy_plot * param.hartree_to_ev
    # ax_hirshfeld_energy[1].plot((time_val_1 - offset)/1e3, potential_energy_plot, 'k-')
    # ax_hirshfeld_energy[1].set_xlabel('Time / ps')
    # ax_hirshfeld_energy[1].set_ylabel('Energy / eV')
    # ax_hirshfeld_energy[1].set_xlim((np.array(xlim_1)-offset)/1000)
    # fig_hirshfeld_energy.tight_layout()
    # fig_hirshfeld_energy.subplots_adjust(hspace=0.05)
    # if save_fig: fig_hirshfeld_energy.savefig('{}/polaron_energy_subplot.png'.format(folder_save), dpi=300)

    plot_ts = True
    if plot_ts:
        fig_spin1_square_ts, ax_spin1_square_ts = plt.subplots()
        temp = np.zeros(num_timesteps)
        # ts_plot_index = np.shape(polaron_distances_hop)[0] - 27 - 4
        ts_plot_index = np.shape(polaron_distances_hop)[0] - 48
        ts_plot_time = 20
        print('time value', time_val_energy[polaron_indices[ts_plot_index]])
        print('time index', polaron_indices[ts_plot_index])
        print('FIRST_SNAPSHOT', int(polaron_indices[ts_plot_index]-ts_plot_time))
        print('LAST_SNAPSHOT', int(polaron_indices[ts_plot_index]+ts_plot_time))

        for j in range(num_atoms):
            ax_spin1_square_ts.plot((time_val_1 - offset), hirshfeld_mobility[:int(xlim_1[1]), 5, j], '.-', label='{}'.format(j + 1))
            # ax_spin1_square_ts.plot((time_val_energy[:int(xlim_1[1])] - offset), hirshfeld_mobility[:int(xlim_1[1]), 5, j], '.-', label='{}'.format(j + 1))

        # hirshfeld_1_np_ts = hirshfeld_1_np[:, 5, j][polaron_indices[ts_plot_index] - ts_plot_time:polaron_indices[ts_plot_index] + ts_plot_time]
        # for j in range(num_atoms):
            # ax_spin1_square_ts.plot(time_val_1_plot[polaron_indices[ts_plot_index] - ts_plot_time:polaron_indices[ts_plot_index] + ts_plot_time],
            #                         hirshfeld_mobility[:int(xlim_1[1]), 5, j][polaron_indices[ts_plot_index] - ts_plot_time:polaron_indices[ts_plot_index] + ts_plot_time], '-', label='{}'.format(j + 1))
            # ax_spin1_square_ts.plot(time_array_ts, hirshfeld_1_np[:int(xlim_1[1]), 5, j], '-.', label='{}'.format(j + 1))
        # ax_spin1_square_ts.plot(time_array, hirshfeld_1_np[:, 5, polaron_atom], 'k-',)
        if draw_legend: ax_spin1_square_ts.legend(frameon=True)
        # ax_spin1_square_ts.set_xlabel('Timestep')
        ax_spin1_square_ts.set_xlabel('Time / fs')
        ax_spin1_square_ts.set_ylabel('Spin moment')
        ax_spin1_square_ts.set_xlim(np.array([time_val_energy[polaron_indices[ts_plot_index]] - ts_plot_time, time_val_energy[polaron_indices[ts_plot_index]] + ts_plot_time]))
        ax_spin1_square_ts.set_ylim(ylim_1)
        fig_spin1_square_ts.tight_layout()
        if save_fig: fig_spin1_square_ts.savefig('{}/hirshfeld_spin_ts_{}.png'.format(folder_save, time_val_energy[polaron_indices[ts_plot_index]]), dpi=300)

# Plot RDF for Ti - O
nbins = 300
rdf_xlim = 13.77 / 2

# if calc_distance:
#     ti_polaron = universe.select_atoms("bynum 78")
#     atoms_o = universe.select_atoms("element O")
#     rdf_ti_o_polaron = rdf.InterRDF(ti_polaron, atoms_o, range=(0, rdf_xlim), nbins=nbins)
#     rdf_ti_o_polaron.run()
#
#     ti_not_polaron = universe.select_atoms("element Ti and not bynum 78")
#     rdf_ti_o_not_polaron = rdf.InterRDF(ti_not_polaron, atoms_o, range=(0, rdf_xlim), nbins=nbins)
#     rdf_ti_o_not_polaron.run()
#
#     fig_rdf_ti_o, ax_rdf_ti_o = plt.subplots(figsize=(10, 4))
#     ax_rdf_ti_o.plot(rdf_ti_o_polaron.bins, rdf_ti_o_polaron.rdf, 'r-', label='RDF Ti (polaron) - O')
#     ax_rdf_ti_o.plot(rdf_ti_o_not_polaron.bins, rdf_ti_o_not_polaron.rdf, 'g-', label='RDF Ti - O')
#     ax_rdf_ti_o.set_xlabel("Radial distance / Å")
#     ax_rdf_ti_o.set_ylabel("RDF (arb. units)")
#     ax_rdf_ti_o.legend(frameon=False)
#     ax_rdf_ti_o.set_xlim([0, rdf_xlim])
#     ax_rdf_ti_o.tick_params(axis='y', which='both', left=False, right=False, labelleft=False)
#     fig_rdf_ti_o.tight_layout()
#     if save_fig: fig_rdf_ti_o.savefig("{}/rdf.png".format(folder_save), dpi=600)
#
# # Plot RDF for Ti - Ti
# if calc_distance:
#     ti_polaron = universe.select_atoms("bynum 78")
#     rdf_ti_ti_polaron = rdf.InterRDF(ti_polaron, ti_not_polaron, range=(0, rdf_xlim), nbins=nbins)
#     rdf_ti_ti_polaron.run()
#
#     ti_not_polaron = universe.select_atoms("element Ti and not bynum 78")
#     rdf_ti_ti_not_polaron = rdf.InterRDF(ti_not_polaron, ti_not_polaron, exclusion_block=(1, 1), range=(0, rdf_xlim), nbins=nbins)
#     rdf_ti_ti_not_polaron.run()
#
#     fig_rdf_ti_ti, ax_rdf_ti_ti = plt.subplots(figsize=(10, 4))
#     ax_rdf_ti_ti.plot(rdf_ti_ti_polaron.bins, rdf_ti_ti_polaron.rdf, 'r-', label='RDF Ti (polaron) - Ti')
#     ax_rdf_ti_ti.plot(rdf_ti_ti_not_polaron.bins, rdf_ti_ti_not_polaron.rdf, 'g-', label='RDF Ti - Ti')
#     ax_rdf_ti_ti.set_xlabel("Radial distance / Å")
#     ax_rdf_ti_ti.set_ylabel("RDF (arb. units)")
#     ax_rdf_ti_ti.legend(frameon=False)
#     ax_rdf_ti_ti.set_xlim([0, rdf_xlim])
#     ax_rdf_ti_ti.tick_params(axis='y', which='both', left=False, right=False, labelleft=False)
#     fig_rdf_ti_ti.tight_layout()
#     if save_fig: fig_rdf_ti_ti.savefig("{}/rdf.png".format(folder_save), dpi=600)

if __name__ == "__main__":
    print('Finished.')
    plt.show()
