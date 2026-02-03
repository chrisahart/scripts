import numpy as np
import pandas as pd


"""
    Convert leopold data to deepmd-kit files
"""


def read_leopold(folder, filename, num_atoms):
    """
    Read lepold data
    """

    file = '{}/{}'.format(folder, filename)
    cols = ['Element', 'coord_x', 'coord_y', 'coord_z', 'pol_state', 'force_x', 'force_y', 'force_z', 'magmom_1',
                'magmom_2', 'magmom_3', 'magmom_4', 'toccup_1', 'toccup_2']
    data_pd = pd.read_csv(file, names=cols, delim_whitespace=True, skiprows=0)

    data_pd = data_pd.dropna()
    data_pd = data_pd.reset_index(drop=True)
    species_pd = data_pd['Element']
    data_pd = data_pd.drop(columns=['Element'])

    data_pd = data_pd.apply(pd.to_numeric, errors='coerce')
    data_pd = data_pd.dropna()
    data_pd = data_pd.dropna(axis='rows', thresh=2)
    data_pd = data_pd.dropna(axis='columns', thresh=1)
    data_pd = data_pd.reset_index(drop=True)
    cols_new = list(data_pd.columns)

    # Loop over each timestep and atoms
    num_timesteps_pd = int(data_pd.shape[0]/num_atoms)
    data_np = np.zeros((num_timesteps_pd, len(cols_new), num_atoms))

    for timestep in range(num_timesteps_pd):
        for atom in range(num_atoms):
            for i in range(len(cols_new)):
                data_np[timestep, i, atom] = data_pd[cols_new[i]].values[atom + timestep * num_atoms]

    return data_pd, data_np, species_pd


folder = '/Volumes/Elements/Data/Postdoc2/Data/Work/resources/other/leopold/data/TiO2/leopold_analysis'
# name_leopold = 'train'
# name_deepmd = 'train'
# name_leopold = 'valid'
# name_deepmd = 'test/1'
# name_leopold = 'test'
# name_deepmd = 'test/2'
name_leopold = 'test'
name_deepmd = 'example'
num_atoms = 288
box = np.array([[12.9730319977, 0, 0, 0, 12.9730319977, 0, 0, 0, 17.7244205475]])
topology_file = '{}/data/system.xyz'.format(folder)

# Read leopold files
data = 'data/{}_clean.xyz'.format(name_leopold)
folder_data1 = 'database_{}'.format(name_deepmd)
energy_clean = np.loadtxt('{}/data/{}_energy.txt'.format(folder, name_leopold))
_, leopold_1_np, _ = read_leopold(folder, data, num_atoms)
num_timesteps = np.shape(leopold_1_np)[0]
coord = leopold_1_np[:, :3, :]
forces = leopold_1_np[:, 4:7, :]
aparam = leopold_1_np[:, 3, :]
population_alpha = leopold_1_np[:, -2, :]
population_beta = leopold_1_np[:, -1, :]

# Convert
coord = np.transpose(coord, axes=(0, 2, 1))
forces = np.transpose(forces, axes=(0, 2, 1))
coord = coord.reshape(num_timesteps, num_atoms*3)
forces = forces.reshape(num_timesteps, num_atoms*3)
energy = np.reshape(energy_clean, (num_timesteps, 1))
box_array = np.zeros((num_timesteps, 9))
for i in range(box_array.shape[0]):
    box_array[i, :] = box
population = np.zeros((np.shape(population_alpha)[0], num_atoms, 2))
for i in range(np.shape(population_alpha)[0]):
    for j in range(num_atoms):
        population[i, j, 0] = population_alpha[i, j]
        population[i, j, 1] = population_beta[i, j]

# Cut data
num_frames = 50
energy = energy[:num_frames]
coord = coord[:num_frames]
forces = forces[:num_frames]
box_array = box_array[:num_frames]
aparam = aparam[:num_frames]
population = population[:num_frames]

# Training
np.save('{}/{}/set.000/energy.npy'.format(folder, folder_data1), energy)
np.save('{}/{}/set.000/coord.npy'.format(folder, folder_data1), coord)
np.save('{}/{}/set.000/force.npy'.format(folder, folder_data1), forces)
np.save('{}/{}/set.000/box.npy'.format(folder, folder_data1), box_array)
np.save('{}/{}/set.000/aparam.npy'.format(folder, folder_data1), aparam)
np.save('{}/{}/set.000/atomic_population.npy'.format(folder, folder_data1), population)

if __name__ == "__main__":
    print('Finished.')
