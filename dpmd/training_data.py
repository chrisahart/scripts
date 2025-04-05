import pandas as pd
import numpy as np
import glob
from matplotlib import pyplot as plt
from general import parameters as param
from matplotlib.colors import LogNorm
import dpdata
import cp2kdata
from cp2kdata import Cp2kOutput

# load data of cp2k/md format
# directory = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/archer/bulk/221_supercell/md/hole/400K/input-file-EPS_SCHWARZ_FORCES-neutral'
# directory = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/archer/bulk/221_supercell/md/hole/cleaned/400k-neutral'
directory = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/archer/bulk/221_supercell/md/hole/cleaned/400k-neutral-2'

# Load the data from the specified directory
data = dpdata.LabeledSystem(directory, cp2k_output_name='cp2k_log.log', fmt="cp2kdata/md", type_map=['Fe_a', 'Fe_b', 'O'])

# Print the total number of frames in the data
print("# the data contains %d frames" % len(data))

# Exclude the first 200 frames
data = data.sub_system(range(200, len(data)))

# Print the number of frames after excluding the first 200 frames
print("# the data contains %d frames after excluding the first 200 frames" % len(data))

# Create the test set from the last 100 frames
data_test = data.sub_system(range(len(data) - 100, len(data)))

# Remove the test set frames from the remaining data
data = data.sub_system(range(len(data) - 100))

# Print the number of frames after creating the test set
print("# the data contains %d frames after creating the test set" % len(data))

# Randomly choose indices for the validation set
rng = np.random.default_rng()
index_validation = rng.choice(len(data), size=400, replace=False)

# Get the indices for the training set
index_training = list(set(range(len(data))) - set(index_validation))

# Create the training and validation datasets
data_training = data.sub_system(index_training)
data_validation = data.sub_system(index_validation)

# Save the training data to the specified directory
data_training.to_deepmd_npy("{}/database_ener_force_train".format(directory))

# Save the validation data to the specified directory
data_validation.to_deepmd_npy("{}/database_ener_force_valid".format(directory))

# Save the test data to the specified directory
data_test.to_deepmd_npy("{}/database_ener_force_test".format(directory))

# Print the number of frames in the training, validation, and test sets
print("# the training data contains %d frames" % len(data_training))
print("# the validation data contains %d frames" % len(data_validation))
print("# the test data contains %d frames" % len(data_test))