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
                

folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/mgo/archer/mgo/cell-222/md/pbe-u-8'
# folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/mgo/archer/mgo/cell-222/md/pbe-u-8/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-600k-csvr-timecon-1-COMVEL_TO-1e-10-nvt'
# folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/mgo/archer/mgo/cell-222/md/pbe-u-8/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-600k-csvr-timecon-1-COMVEL_TO-1e-10-nvt-rs'
# folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/mgo/archer/mgo/cell-222/md/pbe-u-8/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-600k-csvr-timecon-1-COMVEL_TO-1e-10-nvt-hole-u-6'
# polaron_atoms [35 40 45 53 54 58 59 61 62 63 64]
# polaron_distance [2.98178364 5.87473399 3.00040799 5.1273932  3.047131   2.99547136
#  2.99547136 2.85545014 3.01149688 5.83736131 2.89373298 2.89373298
#  2.89373298 2.89373298 3.0580113  5.7055927 ]
# folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/mgo/archer/mgo/cell-222/md/pbe-u-8/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-600k-csvr-timecon-1-COMVEL_TO-1e-10-nvt-hole-u-6-rs'

# folder_1 = '{}/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-600k-csvr-timecon-1-COMVEL_TO-1e-10-nvt-hole-u-6-rs-rs-extrapolation-aspc-2'.format(folder)  # CANCEL, first 4 crossings are good but then loses stability and crossings become discontinuous
# folder_1 = '{}/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-600k-csvr-timecon-1-COMVEL_TO-1e-10-nvt-hole-u-6-rs-rs-extrapolation-aspc-3'.format(folder)  # DECENT CONTINUE as aspc-3-print
folder_1 = '{}/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-600k-csvr-timecon-1-COMVEL_TO-1e-10-nvt-hole-u-6-rs-rs-extrapolation-aspc-3-print'.format(folder)  # DECENT *** CONTINUE
# folder_1 = '{}/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-600k-csvr-timecon-1-COMVEL_TO-1e-10-nvt-hole-u-6-rs-rs-extrapolation-aspc-3-print-mg-tz'.format(folder)
# folder_1 = '{}/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-600k-csvr-timecon-1-COMVEL_TO-1e-10-nvt-hole-u-6-rs-rs-extrapolation-aspc-3-print-mg-tz2p'.format(folder)
folder_1 = '{}/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-600k-csvr-timecon-1-COMVEL_TO-1e-10-nvt-hole-u-6-rs-rs-extrapolation-aspc-3-print-timestep-1.0'.format(folder)  # DECENT *** CONTINUE
# folder_1 = '{}/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-600k-csvr-timecon-1-COMVEL_TO-1e-10-nvt-hole-u-6-rs-rs-extrapolation-aspc-3-print-timestep-1.5'.format(folder)  # DECENT, not as good as timestep 1.0
# folder_1 = '{}/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-600k-csvr-timecon-1-COMVEL_TO-1e-10-nvt-hole-u-6-rs-rs-extrapolation-aspc-3-print-FULL_KINETIC'.format(folder)  # DECENT cancelled
# folder_1 = '{}/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-600k-csvr-timecon-1-COMVEL_TO-1e-10-nvt-hole-u-6-rs-rs-extrapolation-aspc-3-print-broyden'.format(folder)  # DECENT, same as DIIS?
# folder_1 = '{}/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-600k-csvr-timecon-1-COMVEL_TO-1e-10-nvt-hole-u-6-rs-rs-extrapolation-aspc-4'.format(folder)  # CANCEL, worse than aspc-3. Crossings show large jumps
# folder_1 = '{}/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-600k-csvr-timecon-1-COMVEL_TO-1e-10-nvt-hole-u-6-rs-rs-extrapolation-aspc-3-cg'.format(folder)  # DECENT 1.5x more expensive than DIIS with around same accuracy
# folder_1 = '{}/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-600k-csvr-timecon-1-COMVEL_TO-1e-10-nvt-hole-u-6-rs-rs-extrapolation-aspc-3-cg-FULL_ALL'.format(folder)  # FAILED
# folder_1 = '{}/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-600k-csvr-timecon-1-COMVEL_TO-1e-10-nvt-hole-u-6-rs-rs-extrapolation-aspc-3-cg-FULL_ALL-ENERGY_GAP-0.001'.format(folder)  # FAILED
# folder_1 = '{}/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-600k-csvr-timecon-1-COMVEL_TO-1e-10-nvt-hole-u-6-rs-rs-extrapolation-aspc-3-cg-FULL_ALL-ENERGY_GAP-0.001-3PNT'.format(folder)  # FAILED
# folder_1 = '{}/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-600k-csvr-timecon-1-COMVEL_TO-1e-10-nvt-hole-u-6-rs-rs-extrapolation-aspc-3-cg-FULL_ALL-rs-wfn'.format(folder)  # BAD, crossings show large jumps
# folder_1 = '{}/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-600k-csvr-timecon-1-COMVEL_TO-1e-10-nvt-hole-u-6-rs-rs-extrapolation-aspc-3-eps-scf-1e-5'.format(folder)  # BAD, crossings are discontinuous
# folder_1 = '{}/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-600k-csvr-timecon-1-COMVEL_TO-1e-10-nvt-hole-u-6-rs-rs-extrapolation-aspc-3-eps-scf-5e-6'.format(folder)  # BAD, crossings are discontinuous
# folder_1 = '{}/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-600k-csvr-timecon-1-COMVEL_TO-1e-10-nvt-hole-u-6-rs-rs-extrapolation-aspc-3-eps-scf-5e-7'.format(folder)  # BAD, more unstable than 1e-6. Interesting
# folder_1 = '{}/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-600k-csvr-timecon-1-COMVEL_TO-1e-10-nvt-hole-u-6-rs-rs-extrapolation-aspc-3-FULL_ALL'.format(folder)  # FAILED
# folder_1 = '{}/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-600k-csvr-timecon-1-COMVEL_TO-1e-10-nvt-hole-u-6-rs-rs-extrapolation-aspc-3-outer-scf-2'.format(folder)  # BAD, crossings are discontinuous
# folder_1 = '{}/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-600k-csvr-timecon-1-COMVEL_TO-1e-10-nvt-hole-u-6-rs-rs-extrapolation-aspc-3-outer-scf-10-IGNORE_CONVERGENCE_FAILURE'.format(folder)  # CANCEL
# folder_1 = '{}/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-600k-csvr-timecon-1-COMVEL_TO-1e-10-nvt-hole-u-6-rs-rs-extrapolation-aspc-3-timestep-0.75-outer-scf-10-IGNORE_CONVERGENCE_FAILURE'.format(folder)  # CANCEL
# folder_1 = '{}/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-600k-csvr-timecon-1-COMVEL_TO-1e-10-nvt-hole-u-6-rs-rs-extrapolation-aspc-3-timestep-1-outer-scf-10-IGNORE_CONVERGENCE_FAILURE'.format(folder)  # CANCEL
# folder_1 = '{}/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-600k-csvr-timecon-1-COMVEL_TO-1e-10-nvt-hole-u-6-rs-rs-extrapolation-aspc-3-timestep-1.5-outer-scf-10-IGNORE_CONVERGENCE_FAILURE'.format(folder)  # CANCEL
# folder_1 = '{}/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-600k-csvr-timecon-1-COMVEL_TO-1e-10-nvt-hole-u-6-rs-rs-extrapolation-aspc-3-timestep-2-outer-scf-10-IGNORE_CONVERGENCE_FAILURE'.format(folder)
# folder_1 = '{}/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-600k-csvr-timecon-1-COMVEL_TO-1e-10-nvt-hole-u-6-rs-rs-extrapolation-aspc-3-timecon-100'.format(folder)  # DECENT, no difference to TIMECON 1
# folder_1 = '{}/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-600k-csvr-timecon-1-COMVEL_TO-1e-10-nvt-hole-u-6-rs-rs-extrapolation-aspc-3-timestep-0.25'.format(folder)  # DECENT, no difference to TIMESTEP 0.5
# folder_1 = '{}/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-600k-csvr-timecon-1-COMVEL_TO-1e-10-nvt-hole-u-6-rs-rs-extrapolation-aspc-3-timestep-1'.format(folder)  # DECENT, might be better than TIMESTEP 0.5? *** CONTINUE
# folder_1 = '{}/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-600k-csvr-timecon-1-COMVEL_TO-1e-10-nvt-hole-u-6-rs-rs-extrapolation-aspc-4'.format(folder)  # BAD, crossings show large jumps
# folder_1 = '{}/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-600k-csvr-timecon-1-COMVEL_TO-1e-10-nvt-hole-u-6-rs-rs-extrapolation-ps-2'.format(folder)  # BAD, crossings are discontinuous
# folder_1 = '{}/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-600k-csvr-timecon-1-COMVEL_TO-1e-10-nvt-hole-u-6-rs-rs-extrapolation-ps-3'.format(folder)  # DECENT, although twice as expensive as ASPC CANCEL
# folder_1 = '{}/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-600k-csvr-timecon-1-COMVEL_TO-1e-10-nvt-hole-u-6-rs-rs-extrapolation-ps-3-cg'.format(folder)  # DECENT, although 4x cost ASPC CANCEL
# folder_1 = '{}/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-600k-csvr-timecon-1-COMVEL_TO-1e-10-nvt-hole-u-6-rs-rs-extrapolation-ps-4'.format(folder)  # CANCEL, worse than PS-3
# folder_1 = '{}/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-600k-csvr-timecon-1-COMVEL_TO-1e-10-nvt-hole-u-6-rs-rs-extrapolation-USE_GUESS'.format(folder)  # BAD
# folder_1 = '{}/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-600k-csvr-timecon-1-COMVEL_TO-1e-10-nvt-hole-u-6-rs-rs-extrapolation-USE_PREV_P'.format(folder)  # BAD, crossings are discontinuous
# folder_1 = '{}/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-600k-csvr-timecon-1-COMVEL_TO-1e-10-nvt-hole-u-6-rs-rs-extrapolation-USE_PREV_WF'.format(folder)  # BAD, crossings are discontinuous. Good example of MD time vs polaron hopping

# folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/mgo/archer/mgo/cell-222/md/pbe-u-8/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-600k-csvr-timecon-1-COMVEL_TO-1e-10-nvt-hole-u-7'
# polaron_atoms [34 60]
# polaron_distance [3.03819753]
# folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/mgo/archer/mgo/cell-222/md/pbe-u-8/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-600k-csvr-timecon-1-COMVEL_TO-1e-10-nvt-hole-u-8'
# polaron_atoms [36 47 52 55 57 60 64]
# polaron_distance [2.89983298 3.0466969  3.01701384 3.01701384 3.01701384 2.92059895
#  2.82335941 3.06215134 3.06215134]
# folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/mgo/archer/mgo/cell-222/md/pbe-u-8/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-600k-csvr-timecon-1-COMVEL_TO-1e-10-nvt-hole-u-8-rs'
# polaron_atoms [37 48 52 54]
# polaron_distance [3.05700085 2.79967143 2.88156667 2.88156667 2.79967143]
# folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/mgo/archer/mgo/cell-222/md/pbe-u-8/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-600k-csvr-timecon-1-COMVEL_TO-1e-10-nvt-hole-u-9'
# polaron_atoms [35 49 59 60]
# polaron_distance [2.964933   3.06289796 2.89168198]
# folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/mgo/archer/mgo/cell-222/md/pbe-u-8/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-600k-csvr-timecon-1-COMVEL_TO-1e-10-nvt-hole-u-6-lowdin'
# folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/mgo/archer/mgo/cell-222/md/pbe-u-8/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-600k-csvr-timecon-1-COMVEL_TO-1e-10-nvt-hole-u-7-lowdin'
# polaron_atoms [33 36 41 49 52 53 54 58 61 64]
# polaron_distance [2.98182276 2.92679186 2.95280503 2.94787832 2.93923825 2.98919492  2.82341608 5.85930257 2.91022249 5.84971569]
# folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/mgo/archer/mgo/cell-222/md/pbe-u-8/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-600k-csvr-timecon-1-COMVEL_TO-1e-10-nvt-hole-u-8-lowdin'
# polaron_atoms [54]
# polaron_distance []
# folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/mgo/archer/mgo/cell-222/md/pbe-u-8/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-600k-csvr-timecon-1-COMVEL_TO-1e-10-nvt-hole-u-9-lowdin'
# polaron_atoms [54]
# polaron_distance []
# folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/mgo/archer/mgo/cell-222/md/pbe-u-8/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-600k-csvr-timecon-1-COMVEL_TO-1e-10-nvt-hole-u-10-lowdin'
# polaron_atoms [54]
# polaron_distance []
num_atoms = 64
box_size = [8.38, 8.38, 8.38, 90, 90, 90]

# folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/mgo/archer/mgo/cell-333/md/pbe-u-8/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-600k-csvr-timecon-1-COMVEL_TO-1e-10-nvt'
# folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/mgo/archer/mgo/cell-333/md/pbe-u-8/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-600k-csvr-timecon-1-COMVEL_TO-1e-10-nvt-hole-u-6'
# polaron_atoms [154 210]
# polaron_distance [2.88922998]
# folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/mgo/archer/mgo/cell-333/md/pbe-u-8/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-600k-csvr-timecon-1-COMVEL_TO-1e-10-nvt-hole-u-7'
# polaron_atoms [111 116 144 210]
# polaron_distance [2.81204124 3.08964831 3.04692436]
# folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/mgo/archer/mgo/cell-333/md/pbe-u-8/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-600k-csvr-timecon-1-COMVEL_TO-1e-10-nvt-hole-u-8-lowdin'
# folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/mgo/archer/mgo/cell-333/md/pbe-u-8/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-600k-csvr-timecon-1-COMVEL_TO-1e-10-nvt-hole-u-8'
# polaron_atoms [128 209]
# polaron_distance [2.90048786]
# folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/mgo/archer/mgo/cell-333/md/pbe-u-8/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-600k-csvr-timecon-1-COMVEL_TO-1e-10-nvt-hole-u-9'
# polaron_atoms [129 210]
# polaron_distance [2.90050665]
# folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/mgo/archer/mgo/cell-333/md/pbe-u-8/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-600k-csvr-timecon-1-COMVEL_TO-1e-10-nvt-hole-u-10'
# folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/mgo/archer/mgo/cell-333/md/pbe-u-8/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-600k-csvr-timecon-1-COMVEL_TO-1e-10-nvt-hole-u-6-rs'  # Large number unstable hops
# folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/mgo/archer/mgo/cell-333/md/pbe-u-8/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-600k-csvr-timecon-1-COMVEL_TO-1e-10-nvt-hole-u-8-rs'  # Large number unstable hops
# folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/mgo/archer/mgo/cell-333/md/pbe-u-8/geo-opt-hole-u-8-400k'
# folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/mgo/archer/mgo/cell-333/md/pbe-u-8/geo-opt-hole-u-8-600k'

# folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/mgo/archer/mgo/cell-333/md/pbe-u-8/'
# folder_1 = '{}/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-600k-csvr-timecon-1-COMVEL_TO-1e-10-nvt-hole-u-6-rs-print'.format(folder)  # Large number unstable hops
# folder_1 = '{}/neutral-4hours-100k-COMVEL_TO-1e-10-TEMPTOL-10-200k-300k-400k-500k-600k-csvr-timecon-1-COMVEL_TO-1e-10-nvt-hole-u-6-rs-print-timestep-1.0'.format(folder)  # Large number unstable hops
# num_atoms = 216
# box_size = [12.57, 12.57, 12.57, 90, 90, 90]

files = ['mgo-1.ener', 'mgo-charges-1-clean.hirshfeld', 'mgo-pos-1.xyz', 'mgo-frc-1.xyz']
timestep = 0.5
# timestep = 0.25
# timestep = 1.0
local_bonds = 6

print('folder_1', folder_1)
file_energy_1, energy_kinetic_1, energy_potential_1, energy_total_1, temperature_1, time_val_1, time_per_step_1, step_1 = read_energy(folder_1, files[0])
hirshfeld_1_df, hirshfeld_1_np, _ = read_hirsh(folder_1, files[1], num_atoms)
topology_file = '{}/system.xyz'.format(folder_1)
trajectory_file = '{}/mgo-pos-1.xyz'.format(folder_1)
plotting_colors = ['r', 'g', 'b', 'm', 'grey', 'orange', 'brown', 'hotpink'] * 100
plotting_colors = ['b', 'g', 'r', 'm', 'grey', 'orange', 'brown', 'hotpink'] * 100

num_timesteps = np.shape(hirshfeld_1_np)[0]
time_array = np.linspace(start=0, stop=num_timesteps/2,num=num_timesteps)
num_timesteps2 = np.shape(time_per_step_1)[0]
time_array2 = np.linspace(start=0, stop=num_timesteps/2,num=num_timesteps2)
print(num_timesteps, num_timesteps2)
num_timesteps_min = np.min((num_timesteps, num_timesteps2))
print(num_timesteps_min)
# num_timesteps = num_timesteps2
# num_timesteps2 = num_timesteps
num_timesteps = num_timesteps_min
num_timesteps2 = num_timesteps_min

# xlim_1 = [0, 10e3]
# xlim_1 = [0, 1e3]
# xlim_1 = [0, 2e3]
# xlim_1 = [2e3, 4e3]
# xlim_1 = [1700, 1740]
xlim_1 = [0, time_array[-1]]
# xlim_1 = [1.5e3, 1.5e3+(5103+200)/2]

offset = 1500
offset = 0
xlim_1 = np.array(xlim_1) - offset

ylim_1 = [-0.02, 0.9]
draw_legend = False
folder_save = folder_1

print('Polaron:', np.max(hirshfeld_1_np[-1, 5, :]), np.argmax(hirshfeld_1_np[-1, 5, :]))
polaron_atom = np.argmax(hirshfeld_1_np[-1, 5, :])

# Setup md analysis environment
print(num_timesteps2)
universe = mda.Universe(topology_file, trajectory_file)
num_timesteps2 = len(universe.trajectory)
time_val_1 = np.linspace(start=0, stop=len(universe.trajectory) * timestep, num=len(universe.trajectory))
print(time_val_1.shape)
atoms_mg = universe.select_atoms('name Mg')
num_atoms_mg = len(atoms_mg)
atoms_o = universe.select_atoms('name O')
num_atoms_o = len(atoms_o)
dist_arr = distances.distance_array(atoms_o.positions, atoms_mg.positions, box=box_size)
bond_lengths_time = np.zeros((num_timesteps2, num_atoms_o, num_atoms_mg))
for ts in universe.trajectory:
    frame = universe.trajectory.frame
    bond_lengths_time[frame] = distances.distance_array(atoms_o.positions, atoms_mg.positions, box=box_size)
bond_lengths_time_sorted = np.zeros((num_timesteps2, num_atoms_o, local_bonds))
bond_lengths_time_sorted_mean = np.zeros((num_timesteps2, num_atoms_o))
for i in range(num_atoms_mg):
    for j in range(num_timesteps2):
        bond_lengths_time_sorted[j, i] = np.sort(bond_lengths_time[j, i])[0:local_bonds]
        # bond_lengths_time_sorted_mean[j, i] = np.mean(bond_lengths_time_sorted[j, i])
bond_lengths_time_sorted_mean = np.mean(bond_lengths_time_sorted, axis=2)

print('num_timesteps', num_timesteps)
print('num_timesteps2', num_timesteps2)

# Calculate polaron atom
polaron_atom_time = np.zeros(num_timesteps, dtype=int)
# polaron_atom_time = np.zeros(num_timesteps2, dtype=int)
# for j in range(num_timesteps):
#     polaron_atom_time[j] = int(np.argmax(hirshfeld_1_np[j, 5, :]))
# polaron_atoms = np.unique(polaron_atom_time)
# for j in range(num_timesteps2):
for j in range(num_timesteps):
    polaron_atom_time[j] = int(np.argmax(hirshfeld_1_np[j, 5, :]))
polaron_atoms = np.unique(polaron_atom_time)

# polaron_distances = np.zeros(np.shape(polaron_atoms)[0]-1)
# for i in range(polaron_atoms.shape[0]-1):
#     polaron_distances[i] = distances.distance_array(universe.select_atoms('index {}'.format(polaron_atoms[i])).positions,
#                                                 universe.select_atoms('index {}'.format(polaron_atoms[i+1])).positions,
#                                                 box=box_size)
# print('polaron_distance', polaron_distances)

# Calculate distance between current timestep polaron atom and next timestep
# Then get all non-zero answer
# polaron_distances = np.zeros(num_timesteps)
polaron_distances = np.zeros(num_timesteps2)
for j in range(num_timesteps - 1):
# for j in range(num_timesteps2-1):
    polaron_distances[j] = distances.distance_array(universe.select_atoms('index {}'.format(polaron_atom_time[j])).positions,
                                                    universe.select_atoms('index {}'.format(polaron_atom_time[j+1])).positions,
                                                    box=box_size)
print('polaron_distance', polaron_distances)
print('polaron_atoms', polaron_atoms+1)
polaron_distances_hop = polaron_distances[np.nonzero(polaron_distances)]
print('polaron_distances_hop', polaron_distances_hop)

# Plot polaron distances
metric = np.zeros((num_atoms_mg, num_timesteps2))
# fig_bonds_2, ax_bonds_2 = plt.subplots()
fig_bonds_2, ax_bonds_2 = plt.subplots(figsize=(10, 4))
# ax_bonds_2.plot(time_val_1 - time_val_1[0] - offset, polaron_distances[:-1], 'kx-')
ax_bonds_2.plot(time_val_1 - time_val_1[0] - offset, polaron_distances, 'kx-')
ax_bonds_2.set_xlabel('Time / fs')
# ax_bonds_2.set_xlabel('Timestep')
ax_bonds_2.set_ylabel('Polaron hopping distance / A')
# if draw_legend: ax_bonds_2.legend(frameon=True)
# # ax_bonds_2.set_xlim([0, len(universe.trajectory)])
ax_bonds_2.set_xlim(xlim_1)
# ax_bonds_2.set_xlim([0, len(universe.trajectory) * timestep])
# ax_bonds_2.set_ylim([0.06, -0.10])
fig_bonds_2.savefig('{}/polaron_hopping_distance.png'.format(folder_save), dpi=300)
fig_bonds_2.tight_layout()

# Calculate mobility
# hops_distance = np.array([2.99547136, 2.85545014, 3.01149688]) * 1e-8  # Angstrom to cm
# hops_time = (5103+200)/2 * 1e-12  # fs to s
hops_distance = polaron_distances_hop * 1e-8  # Angstrom to cm
hops_time = (time_val_1[-1] - time_val_1[0]) * 1e-15  # ps 1e-12 fs 1e-15
print('hops per ps ', np.shape(hops_distance)[0]/hops_time*1e-15*1e3)

rate_constant = np.shape(hops_distance)[0] / hops_time
print('rate_constant', rate_constant)
print('lifetime fs', 1/rate_constant * 1e15)

mean_distance = np.mean(hops_distance)
print('mean_distance', mean_distance)

site_multiplicity = 1
diffusion_constant_analytical = (np.mean(hops_distance) ** 2 * site_multiplicity * rate_constant) / 2
mobility = (1.60217662e-19 * diffusion_constant_analytical) / (1.380649e-23 * 600)
print('mobility analytical', mobility, 'cm^2/(V·s)')

# print('lifetime hematite', 1/1.2e12 * 1e15)
# test = 1.60217662e-19 * (3e-8**2*3*1.2e12)/2 / (1.380649e-23 * 600)
# print('mobility test hematite', test)

# mobility = 1.60217662e-19 * (3e-8**2*3*rate_constant)/2 / (1.380649e-23 * 600)
# print('mobility analytical 2', mobility, 'cm^2/(V·s)')

mean_square_displacement = np.sum(hops_distance ** 2) / hops_time
diffusion_constant_numerical = mean_square_displacement / (2 * 3)
mobility = (1.60217662e-19 * diffusion_constant_numerical) / (1.380649e-23 * 600)
print('mobility numerical', mobility, 'cm^2/(V·s)')

# Plot average of 6 O-Mg bonds
metric = np.zeros((num_atoms_mg, num_timesteps2))
# fig_bonds_1, ax_bonds_1 = plt.subplots()
fig_bonds_1, ax_bonds_1 = plt.subplots(figsize=(10, 4))
for i in range(num_atoms_o):
    ax_bonds_1.plot(time_val_1 - time_val_1[0] - offset, bond_lengths_time_sorted_mean[:, i], '-')
for j in range(polaron_atoms.shape[0]):
    ax_bonds_1.plot(time_val_1 - time_val_1[0] - offset, bond_lengths_time_sorted_mean[:, polaron_atoms[j]-num_atoms_mg], '-', color=plotting_colors[j], label='{}'.format(polaron_atoms[j]+1))
ax_bonds_1.set_xlabel('Time / fs')
# ax_bonds_1.set_xlabel('Timestep')
ax_bonds_1.set_ylabel('Average of 6 O-Mg bond lengths / A')
if draw_legend: ax_bonds_1.legend(frameon=True)
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

# Plot total energy
# fig_time_md, ax_time_md = plt.subplots()
ylim_1_time = [10, 500]
fig_energy, ax_energy = plt.subplots(figsize=(10, 4))
ax_energy.plot(time_array2, energy_total_1, 'k-')
ax_energy.set_xlabel('Time / fs')
ax_energy.set_ylabel('Energy total / H')
ax_energy.set_xlim(xlim_1)
# ax_energy.set_ylim(ylim_1_time)
fig_energy.tight_layout()
fig_energy.savefig('{}/energy_total.png'.format(folder_save), dpi=300)

# Plot kinetic energy
# fig_time_md, ax_time_md = plt.subplots()
ylim_1_time = [10, 500]
fig_energy_kin, ax_energy_kin = plt.subplots(figsize=(10, 4))
ax_energy_kin.plot(time_array2, energy_kinetic_1, 'k-')
ax_energy_kin.set_xlabel('Time / fs')
ax_energy_kin.set_ylabel('Energy kinetic / H')
ax_energy_kin.set_xlim(xlim_1)
# ax_energy_kin.set_ylim(ylim_1_time)
fig_energy_kin.tight_layout()
fig_energy_kin.savefig('{}/energy_kinetic.png'.format(folder_save), dpi=300)

# Plot potential energy
# fig_time_md, ax_time_md = plt.subplots()
ylim_1_time = [10, 500]
fig_energy_pot, ax_energy_pot = plt.subplots(figsize=(10, 4))
ax_energy_pot.plot(time_array2, energy_potential_1, 'k-')
ax_energy_pot.set_xlabel('Time / fs')
ax_energy_pot.set_ylabel('Energy potential / H')
ax_energy_pot.set_xlim(xlim_1)
# ax_energy_pot.set_ylim(ylim_1_time)
fig_energy_pot.tight_layout()
fig_energy_pot.savefig('{}/energy_potential.png'.format(folder_save), dpi=300)

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
fig_time_md.savefig('{}/time_taken_md_step.png'.format(folder_save), dpi=300)

# Plot spin of all atoms
# polaron_atoms = np.array([53, 62, 63]) - 1
fig_spin1, ax_spin1 = plt.subplots(figsize=(10, 4))
# fig_spin1, ax_spin1 = plt.subplots()
temp = np.zeros(num_timesteps)
print(num_timesteps)
for j in range(num_atoms):
    ax_spin1.plot(time_array-offset, hirshfeld_1_np[:, 5, j], '-')
for j in range(polaron_atoms.shape[0]):
    ax_spin1.plot(time_array-offset, hirshfeld_1_np[:, 5, polaron_atoms[j]], '-', color=plotting_colors[j], label='{}'.format(polaron_atoms[j]+1))
if draw_legend: ax_spin1.legend(frameon=True)
# ax_spin1.set_xlabel('Timestep')
ax_spin1.set_xlabel('Time / fs')
ax_spin1.set_ylabel('Spin moment')
ax_spin1.set_xlim(xlim_1)
ax_spin1.set_ylim(ylim_1)
fig_spin1.tight_layout()
fig_spin1.savefig('{}/spin_all.png'.format(folder_save), dpi=300)

# Plot spin alpha of all atoms
ylim_1 = [3.2, 3.52]
# polaron_atoms = np.array([53, 62, 63]) - 1
fig_pop1, ax_pop1 = plt.subplots(figsize=(10, 4))
# fig_pop1, ax_pop1 = plt.subplots()
temp = np.zeros(num_timesteps)
print(num_timesteps)
for j in range(num_atoms):
    ax_pop1.plot(time_array-offset, hirshfeld_1_np[:, 3, j], '-')
for j in range(polaron_atoms.shape[0]):
    ax_pop1.plot(time_array-offset, hirshfeld_1_np[:, 3, polaron_atoms[j]], '-', color=plotting_colors[j], label='{}'.format(polaron_atoms[j]+1))
if draw_legend: ax_pop1.legend(frameon=True)
# ax_pop1.set_xlabel('Timestep')
ax_pop1.set_xlabel('Time / fs')
ax_pop1.set_ylabel('Alpha spin population')
ax_pop1.set_xlim(xlim_1)
ax_pop1.set_ylim(ylim_1)
fig_pop1.tight_layout()
fig_pop1.savefig('{}/population_alpha_all.png'.format(folder_save), dpi=300)

# Plot pop beta of all atoms
ylim_1 = [3.25, 2.65]
# polaron_atoms = np.array([53, 62, 63]) - 1
fig_pop2, ax_pop2 = plt.subplots(figsize=(10, 4))
# fig_pop2, ax_pop2 = plt.subplots()
temp = np.zeros(num_timesteps)
print(num_timesteps)
for j in range(num_atoms):
    ax_pop2.plot(time_array-offset, hirshfeld_1_np[:, 4, j], '-')
for j in range(polaron_atoms.shape[0]):
    ax_pop2.plot(time_array-offset, hirshfeld_1_np[:, 4, polaron_atoms[j]], '-', color=plotting_colors[j], label='{}'.format(polaron_atoms[j]+1))
if draw_legend: ax_pop2.legend(frameon=True)
# ax_pop2.set_xlabel('Timestep')
ax_pop2.set_xlabel('Time / fs')
ax_pop2.set_ylabel('Beta spin population')
ax_pop2.set_xlim(xlim_1)
ax_pop2.set_ylim(ylim_1)
fig_pop2.tight_layout()
fig_pop2.savefig('{}/population_beta_all.png'.format(folder_save), dpi=300)

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

if __name__ == "__main__":
    print('Finished.')
    plt.show()
