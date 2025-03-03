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
directory = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/chain/units-7/ml/hf/frozen-fe/equil/csvr-1-auto-rs/00.data'
data = dpdata.LabeledSystem('{}/cp2k_md'.format(directory), cp2k_output_name='cp2k_log.log', fmt="cp2kdata/md", type_map=['Fe', 'O', 'H'])
print("# the data contains %d frames" % len(data))

# random choose 40 index for validation_data
rng = np.random.default_rng()
index_validation = rng.choice(len(data), size=100, replace=False)

# other indexes are training_data
index_training = list(set(range(len(data))) - set(index_validation))
data_training = data.sub_system(index_training)
data_validation = data.sub_system(index_validation)

# all training data put into directory:"training_data"
data_training.to_deepmd_npy("{}/training_data".format(directory))

# all validation data put into directory:"validation_data"
data_validation.to_deepmd_npy("{}/validation_data".format(directory))

print("# the training data contains %d frames" % len(data_training))
print("# the validation data contains %d frames" % len(data_validation))
