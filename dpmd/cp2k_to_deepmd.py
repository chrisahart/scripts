import pandas as pd
import numpy as np
import glob
from matplotlib import pyplot as plt
from general import parameters as param
from matplotlib.colors import LogNorm
import dpdata
import cp2kdata
from cp2kdata import Cp2kOutput
# from hirshfeld_to_raw import read_hirsh


def read_hirsh(folder, filename):
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
    num_atoms = int(120)
    num_timesteps = int(file_spec1.shape[0]/num_atoms)
    hirsh_data = np.zeros((num_timesteps, len(cols_new), num_atoms))

    for timestep in range(num_timesteps):
        for atom in range(num_atoms):
            for i in range(len(cols_new)):
                hirsh_data[timestep, i, atom] = file_spec1[cols_new[i]].values[atom + timestep * num_atoms]

    return file_spec1, hirsh_data, species

# directory = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/archer/chain/units-7/hf/md/frozen-fe-h2o/csvr-1/'
# directory = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/archer/chain/units-7/hse/frozen-fe-h2o/md/hf-35-short-auto-eps-1e-5-full-all-xyz-15A-MIN_PAIR_LIST_RADIUS'
# polaron_atom = 2
# num_atoms = 77

# 400k-f
# directory = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/chris_phd/bulk/hole/hse/400k-f'
# polaron_atom = 0

# 400k-b
# directory = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/chris_phd/bulk/hole/hse/400k-b'
# polaron_atom = 13

# REFTRAJ
# directory = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/archer/bulk/221_supercell/reftraj/neb/concat/all'
# directory = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/archer/bulk/221_supercell/reftraj/neb/concat/train'
# directory = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/archer/bulk/221_supercell/reftraj/neb/concat/test'
# polaron_atom = 12

# REFTRAJ + 400k-b
# directory = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/archer/bulk/221_supercell/reftraj/neb/concat/400k-b-geo-opt-all'
# polaron_atom = 12

# REFTRAJ 0K 221
# directory = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/archer/bulk/221_supercell/reftraj/neb/atom-0-3/database'
# polaron_atom = 12
# num_fe = 47
# num_atoms = 120

# REFTRAJ 0K 441
# directory = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/archer/bulk/441_supercell/reftraj/neb/atom-0-1/database/geo-opt-441'
# polaron_atom = 156
# num_fe = 192
# num_atoms = 480
# num_fe = num_atoms

# MD 400K neutral 221
# directory = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/archer/bulk/221_supercell/md/hole/400K/input-file-EPS_SCHWARZ_FORCES-neutral'
# directory = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/archer/bulk/221_supercell/md/hole/400K/input-file-EPS_SCHWARZ_FORCES-neutral2'
# directory = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/archer/bulk/221_supercell/md/hole/cleaned-checked/400k-neutral'
directory = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/archer/bulk/221_supercell/md/hole/cleaned-checked/400k-charged'

polaron_atom = 13
num_atoms = 120

size_exclude_start = 3600
size_exclude_end = 8820 - 5600
size_test = 250
size_validation = 250

# Energy and forces
data = dpdata.LabeledSystem(directory, cp2k_output_name='cp2k_log.log', fmt="cp2kdata/md", type_map=['Fe_a', 'Fe_b', 'O'])

data = data.sub_system(range(size_exclude_start, len(data) - size_exclude_end))
print("# the data contains %d frames after excluding the first and last frames" % len(data))

data_test = data.sub_system(range(len(data) - size_test, len(data)))
data = data.sub_system(range(len(data) - size_test))
rng = np.random.default_rng()
index_validation = rng.choice(len(data), size=size_validation, replace=False)
index_training = list(set(range(len(data))) - set(index_validation))
data_training = data.sub_system(index_training)
data_validation = data.sub_system(index_validation)

print("# the training data contains %d frames" % len(data_training))
print("# the validation data contains %d frames" % len(data_validation))
print("# the test data contains %d frames" % len(data_test))

data_training.to_deepmd_npy("{}/database_ener_force_train".format(directory))
data_validation.to_deepmd_npy("{}/database_ener_force_test/1".format(directory))
data_test.to_deepmd_npy("{}/database_ener_force_test/2".format(directory))

# Hirshfeld spin moment / magnetization m (atom_ener.npy)
hirshfeld_1_df, hirshfeld_1_np, species = read_hirsh(directory, 'hematite-charges-1-clean.hirshfeld')

species_list = list(species[2:num_atoms+2].values)
spin_moment = hirshfeld_1_np[size_exclude_start:-size_exclude_end, 5, :]
print(spin_moment.shape)

num_timesteps = np.shape(spin_moment)[0]
print(species_list)
print(num_timesteps)

# Polaron species identifier / charge state c (aparam.npy)
aparam = np.zeros((num_timesteps, num_atoms))
aparam[:, polaron_atom] = 1.0

species_indices = [j for j, s in enumerate(species_list[0:num_atoms+1]) if s == 'Fe']
num_fe = len(species_indices)
print(num_fe)

for i in range(num_timesteps-1):
    spin_moment_fe = spin_moment[i, species_indices]
    value = np.min(np.abs(spin_moment_fe))
    index = np.argmin(np.abs(spin_moment_fe))
    original_index = species_indices[index]
    aparam[i + 1, original_index] = 1.0

    # print(i, original_index, value)
    # print(species[original_index])

spin_moment_test = spin_moment[-size_test:]
aparam_test = aparam[-size_test:]
print(spin_moment_test.shape)
np.savetxt('{}/database_ener_force_test/2/atom_ener.raw'.format(directory), spin_moment_test.flatten(), delimiter=' ')
np.save('{}/database_ener_force_test/2/set.000/atom_ener.npy'.format(directory), spin_moment_test.flatten())
np.savetxt('{}/database_ener_force_test/2/aparam.raw'.format(directory), aparam_test.flatten(), delimiter=' ')
np.save('{}/database_ener_force_test/2/set.000/aparam.npy'.format(directory), aparam_test.flatten())

spin_moment_validation = spin_moment[index_validation]
aparam_validation = aparam[index_validation]
print(spin_moment_validation.shape)
np.savetxt('{}/database_ener_force_test/1/atom_ener.raw'.format(directory), spin_moment_validation.flatten(), delimiter=' ')
np.save('{}/database_ener_force_test/1/set.000/atom_ener.npy'.format(directory), spin_moment_validation.flatten())
np.savetxt('{}/database_ener_force_test/1/aparam.raw'.format(directory), aparam_validation.flatten(), delimiter=' ')
np.save('{}/database_ener_force_test/1/set.000/aparam.npy'.format(directory), aparam_validation.flatten())

spin_moment_train = spin_moment[index_training]
aparam_train = aparam[index_training]
print(spin_moment_train.shape)
np.savetxt('{}/database_ener_force_train/atom_ener.raw'.format(directory), spin_moment_train.flatten(), delimiter=' ')
np.save('{}/database_ener_force_train/set.000/atom_ener.npy'.format(directory), spin_moment_train.flatten())
np.savetxt('{}/database_ener_force_train/aparam.raw'.format(directory), aparam_train.flatten(), delimiter=' ')
np.save('{}/database_ener_force_train/set.000/aparam.npy'.format(directory), aparam_train.flatten())

# # for i in range(size_exclude, num_timesteps):
# #     for j in range(num_atoms):
# #         # data[i, j] = np.abs(test['Spin moment'][(i*num_atoms) + j])
# #         data[i, j] = test['Spin moment'][(i*num_atoms) + j]
# # print('data[0, :]', data[0, :])
# data_flat = data.flatten()
#
# np.savetxt('{}/{}/atom_ener.raw'.format(directory, directory_out), data_flat, delimiter=' ')
# np.save('{}/{}/set.000/atom_ener.npy'.format(directory, directory_out), data_flat)
#

# Polaron species identifier / charge state c (aparam.npy)
# aparam[:, polaron_atom] = 1.0
# aparam[:, polaron_atom] = 1.0
# aparam[:, polaron_atom] = 1.0
# aparam = np.zeros((data.shape[0], num_atoms))
# aparam[0, polaron_atom] = 1.0
# # aparam[1, polaron_atom] = 1.0
# # aparam[2, polaron_atom] = 1.0
# timestep_start = 0
#
# species_indices = [j for j, s in enumerate(species[0:num_atoms+1]) if s == 'Fe']
# num_fe = len(species_indices)
# print(num_fe)
#
# for i in range(timestep_start, data.shape[0] - 1):
#     data_fe = data[i, species_indices]
#     value = np.min(np.abs(data_fe))
#     index = np.argmin(np.abs(data_fe))
#     original_index = species_indices[index]
#     aparam[i + 1, original_index] = 1.0
#
#     print(i, original_index, value)
#     print(species[original_index])
#
# print(aparam[:10, polaron_atom])
# print(aparam[-10:, original_index])
#
# aparam_flat = aparam.flatten()
# np.savetxt('{}/{}/aparam.raw'.format(directory, directory_out), aparam_flat, delimiter=' ')
# np.save('{}/{}/set.000/aparam.npy'.format(directory, directory_out), aparam_flat)
# test = np.load('{}/{}/set.000/aparam.npy'.format(directory, directory_out))
# # print(test)

print('Finished')
