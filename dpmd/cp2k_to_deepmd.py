import pandas as pd
import numpy as np
import glob
from matplotlib import pyplot as plt
from general import parameters as param
from matplotlib.colors import LogNorm
import dpdata
import cp2kdata
from cp2kdata import Cp2kOutput
from hirshfeld_to_raw import read_hirsh

# directory = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/archer/chain/units-7/hf/md/frozen-fe-h2o/csvr-1/'
# directory = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/archer/chain/units-7/hse/frozen-fe-h2o/md/hf-35-short-auto-eps-1e-5-full-all-xyz-15A-MIN_PAIR_LIST_RADIUS'
# polaron_atom = 2
# num_atoms = 77

directory = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/chris_phd/bulk/hole/hse/400k-f'
polaron_atom = 0

directory = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/chris_phd/bulk/hole/hse/400k-b'
polaron_atom = 13

num_fe = 47
num_atoms = 120

directory_out = 'database'

# file_output = Cp2kOutput('cp2k_log.log', path_prefix=directory, run_type="md")
# print(file_output)

# Energy and forces
dp = dpdata.LabeledSystem(directory, cp2k_output_name='cp2k_log.log', fmt="cp2kdata/md", type_map=['Fe_a', 'Fe_b', 'O'])
# dp = dpdata.LabeledSystem(directory, cp2k_output_name='test', fmt="cp2kdata/md", type_map=['Fe', 'O', 'H'])
dp.to_deepmd_npy('{}/{}'.format(directory, directory_out))
print(dp)

# Hirshfeld spin moment / magnetization m (atom_ener.npy)
test, species = read_hirsh('{}/hematite-charges-1-clean.hirshfeld'.format(directory))
num_timesteps = int(len(test)/num_atoms)
data = np.zeros((num_timesteps, num_atoms))
for i in range(num_timesteps):
    for j in range(num_atoms):
        data[i, j] = np.abs(test['Spin moment'][(i*num_atoms) + j])
        # data[i, j] = test['Spin moment'][(i*num_atoms) + j]

# Remove data
# start_timestep = 0
# end_timestep = 1107
# data = data[start_timestep:end_timestep, :]
# print(data.shape)
# print(end_timestep*77)

data_flat = data.flatten()
# print(data_flat[-120:])

np.savetxt('{}/{}/atom_ener.raw'.format(directory, directory_out), data_flat, delimiter=' ')
np.save('{}/{}/set.000/atom_ener.npy'.format(directory, directory_out), data_flat)

# Polaron species identifier / charge state c (aparam.npy)
# aparam = np.zeros((num_timesteps, num_atoms))
# aparam[:, polaron_atom] = 1.0

aparam = np.zeros((data.shape[0], num_atoms))
aparam[0, polaron_atom] = 1.0
aparam[1, polaron_atom] = 1.0
aparam[2, polaron_atom] = 1.0
timestep_start = 3
for i in range(timestep_start, data.shape[0]-1):
        value = np.min(np.abs(data[i, 0:num_fe+1]))
        index = np.argmin(np.abs(data[i, 0:num_fe+1]))
        print(i, index, value)
        aparam[i+1, index] = 1.0

aparam_flat = aparam.flatten()
np.savetxt('{}/{}/aparam.raw'.format(directory, directory_out), aparam_flat, delimiter=' ')
np.save('{}/{}/set.000/aparam.npy'.format(directory, directory_out), aparam_flat)
test = np.load('{}/{}/set.000/aparam.npy'.format(directory, directory_out))
print(test)