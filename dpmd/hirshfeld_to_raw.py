import pandas as pd
import numpy as np
import glob
from matplotlib import pyplot as plt
from general import parameters as param
from matplotlib.colors import LogNorm
import dpdata
import cp2kdata
from cp2kdata import Cp2kOutput


def read_hirsh(filename):
    """
    Read Hirshfeld analysis from CP2K output file
    """

    # Read number of atoms and labels from .xyz file
    cols = ['Atom', 'Element', 'Kind', 'Ref Charge', 'Pop 1', 'Pop 2', 'Spin moment', 'Charge', 'A', 'B']
    file_spec1 = pd.read_csv(filename, names=cols, delim_whitespace=True, skiprows=5)
    species = file_spec1['Element']
    # print(file_spec1)

    file_spec1 = file_spec1.drop(columns=['A'])
    file_spec1 = file_spec1.drop(columns=['B'])
    file_spec1 = file_spec1.drop(columns=['Element'])

    file_spec1 = file_spec1.apply(pd.to_numeric, errors='coerce')
    file_spec1 = file_spec1.dropna(axis='rows', thresh=2)
    # file_spec1 = file_spec1.dropna(axis='columns', thresh=1)
    file_spec1 = file_spec1.reset_index(drop=True)
    print(file_spec1)

    return file_spec1, species


# directory = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/archer/chain/units-7/hf/md/frozen-fe-h2o/csvr-1/'
# directory = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/archer/chain/units-7/hse/frozen-fe-h2o/md/hf-50-short-auto-rs/'
# directory = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/archer/chain/units-7/hse/frozen-fe-h2o/geo_opt/hf-50-schwartz-1e-6'
# directory = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/archer/chain/units-7/hse/frozen-fe-h2o/geo_opt/hf-25-schwartz-1e-6'
# filename_hirsh_0 = '{}/fe_chain-hirshfeld-1.hirshfeld'.format(directory)

# directory = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/archer/chain/units-7/hse/frozen-fe-h2o/md/hf-50-short-auto-rs'
# directory = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/archer/chain/units-7/hse/frozen-fe-h2o/md/hf-35-short-auto-rs'
# directory = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/archer/chain/units-7/hse/frozen-fe-h2o/md/hf-25-short-auto-rs'
# filename_hirsh_0 = '{}/test'.format(directory)

# directory = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/archer/chain/units-7/hse/frozen-fe-h2o/md/hf-12-short-auto-rs'
# directory = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/archer/chain/units-7/hse/frozen-fe-h2o/md/hf-20-short-auto-rs'
# directory = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/archer/chain/units-7/hse/frozen-fe-h2o/md/hf-25-short-auto-rs'
# directory = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/archer/chain/units-7/hse/frozen-fe-h2o/md/hf-30-short-auto-rs'
# directory = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/archer/chain/units-7/hse/frozen-fe-h2o/md/hf-30-short-auto-eps-1e-5-full-all-xyz-15A-MIN_PAIR_LIST_RADIUS'
# directory = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/archer/chain/units-7/hse/frozen-fe-h2o/md/hf-35-short-auto-eps-1e-5-full-all'
# directory = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/archer/chain/units-7/hse/frozen-fe-h2o/md/hf-35-short-auto-eps-1e-5-full-all-xyz-15A-MIN_PAIR_LIST_RADIUS'
# directory = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/archer/chain/units-7/hse/frozen-fe-h2o/md/hf-35-short-auto-eps-1e-5-full-all-xyz-15A-MIN_PAIR_LIST_RADIUS-rs-24hours-temp-800'
# directory = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/archer/chain/units-7/hse/frozen-fe-h2o/md/hf-35-short-auto-eps-1e-5-full-all-xyz-15A-MIN_PAIR_LIST_RADIUS-rs-24hours'
# directory = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/archer/chain/units-7/hse/frozen-fe-h2o/md/hf-35-short-auto-eps-1e-5-full-all-xyz-15A-MIN_PAIR_LIST_RADIUS-rs-24hours-tight'
# directory = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/archer/chain/units-7/hse/frozen-fe-h2o/md/hf-35-short-auto-eps-1e-6-full-all-xyz-15A-MIN_PAIR_LIST_RADIUS-IGNORE'
# directory = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/archer/chain/units-7/hse/frozen-fe-h2o/md/hf-35-short-auto-eps-1e-5-full-all-temp-800'
# directory = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/archer/chain/units-7/hse/frozen-fe-h2o/md/hf-35-eps-1e-6-temp-400'
# directory = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/archer/chain/units-7/hse/frozen-fe-h2o/md/hf-50-short-auto-rs'
# filename_hirsh_0 = '{}/hirshfeld/combined.hirshfeld'.format(directory)
# filename_hirsh_0 = '{}/hematite-charges-1-clean.hirshfeld'.format(directory)

# directory = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/chris_phd/bulk/hole/hse/hops-1'
# filename_hirsh_0 = '{}/hematite-charges-1-clean.hirshfeld'.format(directory)

directory = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/chris_phd/bulk/hole/hse/400k'
directory = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/chris_phd/bulk/hole/hse/400k-f'
filename_hirsh_0 = '{}/hematite-charges-1-clean.hirshfeld'.format(directory)

# directory = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/chris_phd/bulk/hole/philipp_schienbein/analysis/hirshfeld'
# filename_hirsh_0 = '{}/330K.out'.format(directory)

directory_out = 'database'
test, species = read_hirsh(filename_hirsh_0)

# num_atoms = 77
num_atoms = 120
num_timesteps = int(len(test)/num_atoms)
print(num_timesteps)
time_array = np.linspace(0, int(num_timesteps/2), num=num_timesteps)

data_raw = np.zeros((num_timesteps, num_atoms))
for i in range(num_timesteps):
    for j in range(num_atoms):
        data_raw[i, j] = test['Spin moment'][(i*num_atoms) + j]

fig_plot_1, ax_plot_1 = plt.subplots()
# num_fe = 7
num_fe = 120
for i in range(num_fe):
    # ax_plot_1.plot(data_raw[:, i], '-', label='Fe {}'.format(i))
    ax_plot_1.plot(time_array, data_raw[:, i], '-', label='Fe {}'.format(i))
ax_plot_1.legend(frameon=False)
ax_plot_1.set_xlabel('Time / fs')
# ax_plot_1.set_xlabel('Timestep')
ax_plot_1.set_ylabel('Fe spin moment')
fig_plot_1.tight_layout()
fig_plot_1.savefig('{}/hirshfeld.png'.format(directory), dpi=300)

# np.savetxt('{}/{}/spin.raw'.format(directory, directory_out), data_raw, delimiter=' ')
# data = np.loadtxt('{}/{}/spin.raw'.format(directory, directory_out), ndmin=2)
# data = data.astype(np.float64)
# np.save('{}/{}/spin'.format(directory, directory_out), data)

if __name__ == "__main__":
    print('Finished.')
    plt.show()

