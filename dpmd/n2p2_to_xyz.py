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
    # print(species_pd.shape)
    # print(species_pd[350:])

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
    # print(data_pd.shape)

    # Loop over each timestep and atoms
    print('data_pd.shape[0]/num_atoms', data_pd.shape[0]/num_atoms)
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


# rutile 100 (100 structures instead of 20 as for selected-20 in /cp2k-aimd-bulk)
# much higher water density than file in /cp2k-aimd-bulk
# folder = '/Users/chris/Documents/Storage/other/tio2_water/cp2k-aimd/rutile/bulk/100'
# filename_in = 'selected-100-cp2k-md-rutile-100-water64.n2p2'
# filename_out = 'selected-100-cp2k-md-rutile-100-water64.xyz'
# filename_out2 = 'all/selected-100-cp2k-md-rutile-100-water64'
# box = np.array([[26.380198, 16.832924, 51.226884]])
# box = box * 0.529
# print(box)
# num_atoms = 354

# rutile 100
# folder = '/Users/chris/Documents/Storage/other/tio2_water/cp2k-aimd-bulk/rutile/bulk/100'
# filename_in = 'selected-20-cp2k-md-rutile-100.n2p2'
# filename_out = 'selected-20-cp2k-md-rutile-100.xyz'
# filename_out2 = 'all/selected-20-cp2k-md-rutile-100'
# box = np.array([[26.380576, 16.837460, 59.241590]])
# box = box * 0.529
# print(box)
# num_atoms = 354

# rutile 110
# folder = '/Users/chris/Documents/Storage/other/tio2_water/cp2k-aimd-bulk/rutile/bulk/110'
# filename_in = 'selected-20-cp2k-md-rutile-110.n2p2'
# filename_out = 'selected-20-cp2k-md-rutile-110.xyz'
# filename_out2 = 'all/selected-20-cp2k-md-rutile-110'
# box = np.array([[24.868795, 16.837460, 59.337400]])
# box = box * 0.529
# print(box)
# num_atoms = 336

# rutile 101 == 011
# folder = '/Users/chris/Documents/Storage/other/tio2_water/cp2k-aimd-bulk/rutile/bulk/101'
# filename_in = 'selected-20-cp2k-md-rutile-101.n2p2'
# filename_out = 'selected-20-cp2k-md-rutile-101.xyz'
# filename_out2 = 'all/selected-20-cp2k-md-rutile-101'
# num_atoms = 300
# box = np.array([[20.862198, 26.380198, 41.434890]])
# box = box * 0.529
# print(box)

# rutile 001
# folder = '/Users/chris/Documents/Storage/other/tio2_water/cp2k-aimd-bulk/rutile/bulk/001'
# filename_in = 'selected-20-cp2k-md-rutile-001.n2p2'
# filename_out = 'selected-20-cp2k-md-rutile-001.xyz'
# filename_out2 = 'all/selected-20-cp2k-md-rutile-001'
# num_atoms = 354
# box = np.array([[26.380576, 26.380576, 39.169298]])
# box = box * 0.529
# print(box)

# anatase 100
folder = '/Users/chris/Documents/Storage/other/tio2_water/cp2k-aimd-bulk/anatase/bulk/001'
filename_in = 'selected-20-cp2k-md-anatase-001.n2p2'
filename_out = 'selected-20-cp2k-md-anatase-001.xyz'
filename_out2 = 'all/selected-20-cp2k-md-anatase-001'
num_atoms = 300
box = np.array([[21.557995, 21.557995, 49.986657]])
box = box * 0.529
print(box)

# anatase 101
# folder = '/Users/chris/Documents/Storage/other/tio2_water/cp2k-aimd-bulk/anatase/bulk/101'
# filename_in = 'selected-20-cp2k-md-anatase-101.n2p2'
# filename_out = 'selected-20-cp2k-md-anatase-101.xyz'
# filename_out2 = 'all/selected-20-cp2k-md-anatase-101'
# num_atoms = 288
# box = np.array([[19.772582, 14.372123, 76.001949]])
# box = box * 0.529
# print(box)

# anatase 110 (100 structures instead of 20 as for selected-20 in /cp2k-aimd-bulk)
# no file in /cp2k-aimd-bulk
# folder = '/Users/chris/Documents/Storage/other/tio2_water/cp2k-aimd/anatase/bulk/110'
# filename_in = 'selected-100-cp2k-md-anatase-110-water64.n2p2'
# filename_out = 'selected-100-cp2k-md-anatase-110-water64.xyz'
# filename_out2 = 'all/selected-100-cp2k-md-anatase-110-water64'
# num_atoms = 288
# box = np.array([[20.230652, 35.978873, 27.042169]])
# box = box * 0.529
# print(box)

data_np, species, species_array = read_n2p2(folder, filename_in, num_atoms)
data_np = data_np * 0.529
num_timesteps = np.shape(data_np)[0]

# print(str(species_array[0].flatten()))
# print(np.shape(species_array))
np.savetxt('{}/species'.format(folder, filename_out), [str(species_array[0].flatten())], fmt='%s')
write_xyz('{}/{}'.format(folder, filename_out), data_np, species_array, num_atoms)
write_xyz_all('{}/{}'.format(folder, filename_out2), data_np, species_array, num_atoms)
