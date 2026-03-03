import numpy as np
import pandas as pd


"""
    Convert leopold data to deepmd-kit files
"""


def read_n2p2(folder, filename, num_atoms):
    """
    Read lepold data
    """

    file = '{}/{}'.format(folder, filename)
    cols = ['atom', 'coord_x', 'coord_y', 'coord_z', 'species', 'pop_alpha', 'pop_beta', 'force_x', 'force_y', 'force_z']
    data_pd = pd.read_csv(file, names=cols, delim_whitespace=True, skiprows=0)

    data_pd = data_pd.dropna()
    data_pd = data_pd.reset_index(drop=True)

    species_pd = data_pd['species']
    species_pd = species_pd.reset_index(drop=True)
    print(species_pd)

    data_pd = data_pd.drop(columns=['atom'])
    data_pd = data_pd.drop(columns=['species'])
    data_pd = data_pd.drop(columns=['pop_alpha'])
    data_pd = data_pd.drop(columns=['pop_beta'])
    data_pd = data_pd.drop(columns=['force_x'])
    data_pd = data_pd.drop(columns=['force_y'])
    data_pd = data_pd.drop(columns=['force_z'])
    data_pd = data_pd.apply(pd.to_numeric, errors='coerce')
    data_pd = data_pd.dropna()
    data_pd = data_pd.dropna(axis='rows', thresh=2)
    data_pd = data_pd.dropna(axis='columns', thresh=1)
    data_pd = data_pd.reset_index(drop=True)
    cols_new = list(data_pd.columns)

    # Loop over each timestep and atoms
    num_timesteps_pd = int(data_pd.shape[0]/num_atoms)
    species = species_pd[:num_atoms]
    print('data_pd.shape[0], num_timesteps_pd, num_atoms',
          data_pd.shape[0], num_timesteps_pd, num_atoms)
    data_np = np.zeros((num_timesteps_pd, len(cols_new), num_atoms))
    species_array = np.empty((num_timesteps_pd, num_atoms), dtype='U2')
    for timestep in range(num_timesteps_pd):
        for atom in range(num_atoms):
            species_array[timestep, atom] = species_pd.values[atom + timestep * num_atoms]
            for i in range(len(cols_new)):
                data_np[timestep, i, atom] = data_pd[cols_new[i]].values[atom + timestep * num_atoms]

    return data_np, species, species_array


def write_xyz(filename, coordinates, species, num_atoms):

    num_timesteps = np.shape(coordinates)[0]
    with open(filename, 'w') as f:
        for timestep in range(num_timesteps):
            f.write(f"{num_atoms}\n")
            f.write(f"\n")
            # f.write(f"i = \n")

            for atom in range(num_atoms):
                x, y, z = coordinates[timestep, :, atom]
                f.write(f"{species[timestep, atom]} {x:.10f} {y:.10f} {z:.10f}\n")


def write_xyz_all(base_filename, coordinates, species, num_atoms):

    num_timesteps = np.shape(coordinates)[0]
    for timestep in range(num_timesteps):
        filename = f"{base_filename}_{timestep}.xyz"

        with open(filename, 'w') as f:
            f.write(f"{num_atoms}\n")
            f.write(f"\n")

            # Write atomic coordinates
            for atom in range(num_atoms):
                x = coordinates[timestep, 0, atom]
                y = coordinates[timestep, 1, atom]
                z = coordinates[timestep, 2, atom]
                f.write(f"{species[timestep, atom]:2s} {x:12.6f} {y:12.6f} {z:12.6f}\n")


folder = '/Users/chris/Documents/Storage/other/tio2_water/data/anatase/vacuum'

# selected-200-cp2k-md-anatase-001-v-p0.02
filename_in = 'selected-200-cp2k-md-anatase-001-v-p0.02.n2p2'
filename_out = 'selected-200-cp2k-md-anatase-001-v-p0.02.xyz'
filename_out2 = 'all/selected-200-cp2k-md-anatase-001-v-p0.02'
num_atoms = 288
box = np.array([[21.126760, 0, 0, 0, 21.126760, 0, 0, 0, 78.332548]])

data_np, species, species_array = read_n2p2(folder, filename_in, num_atoms)
data_np = data_np * 0.529
num_timesteps = np.shape(data_np)[0]

write_xyz('{}/{}'.format(folder, filename_out), data_np, species_array, num_atoms)
write_xyz_all('{}/{}'.format(folder, filename_out2), data_np, species_array, num_atoms)
