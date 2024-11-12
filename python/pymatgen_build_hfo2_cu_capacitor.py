from ase import Atoms
from ase.visualize import view
from ase.io import read, write
from ase.build import fcc100, fcc111
import numpy as np
import pandas as pd
from general import print_xyz

folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/structures/pymatgen/5-5-4-layers'
# folder = r'C:\Users\storm\Desktop\old'
polymorph = 't-hfo2'
cu_lattice = 3.61
hfo2_lattice_a = 5.1
cu_lattice_spacing_xy = cu_lattice*np.sqrt(2)/2
cu_o = 5.57326 - 3.73217
layers_xy = 5
layers_z = 4

hafnia_supercell = read("{}/{}_supercell.xyz".format(folder, polymorph), format='xyz')
cu_slab = fcc100('Cu', size=(6, 6, 11), a=cu_lattice, periodic=True)
junction = Atoms(pbc=True)

cu_max = np.max(cu_slab.positions[:, 2])
cu_xy_unique = np.unique(cu_slab.positions[:, 0])
cu_z_unique = np.unique(cu_slab.positions[:, 2])

hafnia_xy_unique = np.unique(hafnia_supercell.positions[:, 0])
cut_index = layers_xy * 2 - len(hafnia_xy_unique)
cut_index_cu = cut_index
hafnia_z_unique = np.unique(hafnia_supercell.positions[:, 2])
z_cut_hafnia = hafnia_z_unique[2]

print(cu_xy_unique)
print(cu_xy_unique[cut_index_cu])

for atom in cu_slab:
    if (atom.position[0] < cu_xy_unique[cut_index_cu]
            and atom.position[1] < cu_xy_unique[cut_index_cu]):
        junction.append(atom)

z_cut = np.max(cu_slab.positions[:, 2]) - 1
junction_max = hfo2_lattice_a/2 * layers_z + np.max(cu_slab.positions[:, 2]) + cu_o

for atom in cu_slab:
    if (atom.position[2] < z_cut and
            atom.position[0] < cu_xy_unique[cut_index_cu]
            and atom.position[1] < cu_xy_unique[cut_index_cu]):
        atom.position[2] = atom.position[2] + junction_max + cu_o
        junction.append(atom)

# label au atoms
unique_z = np.unique(junction.positions[:, 2])
# species = []
# j = 1
# for i in range(unique_z.shape[0]):
#     for atom in junction:
#         if atom.symbol == 'Cu':
#             if atom.position[2] == unique_z[0]:
#                 species.append('Cu_{}'.format(str(j)))
#             j=j+1

for atom in hafnia_supercell:
    if (atom.position[2] > z_cut_hafnia
            and atom.position[0] < hafnia_xy_unique[cut_index]
            and atom.position[1] < hafnia_xy_unique[cut_index]):
        atom.position[2] = atom.position[2] - hafnia_z_unique[3] + cu_max + cu_o
        atom.position[0] = atom.position[0] - 1.27500
        atom.position[1] = atom.position[1] + (2.60000 - 1.27633)
        junction.append(atom)

# Remove Hf and O to form capacitor
# mask = np.isin(junction.symbols, ['Hf', 'O'], invert=True)
# junction = Atoms(symbols=junction.symbols[mask],
#                  positions=junction.positions[mask])

# cell_size = np.array([10.21, 10.21, 49.982])
cell_size = np.array([cu_lattice_spacing_xy*layers_xy,
                      cu_lattice_spacing_xy*layers_xy,
                      np.max(junction.positions[:, 2])+cu_lattice/2])
junction.set_cell(cell_size)

# Wrap the positions within the new cell and sort
junction.wrap()
sorted_indices = np.lexsort((junction.positions[:, 1],
                             junction.positions[:, 0],
                             junction.positions[:, 2]))
junction = junction[sorted_indices]

# Round all atom positions to 3 dp
for atom in junction:
    atom.position = np.round(atom.position, 3)
print('cell_size', np.round(cell_size, 3))

# Adjust species to "Cu_bulk" based on z-coordinate condition
species = list(junction.symbols)
print(species)
for i, atom in enumerate(junction):
    if atom.symbol == 'Cu' and (atom.position[2] <= 10.830 or atom.position[2] >= 39.15200):
        species[i] = 'Au_bulk'
    elif atom.symbol == 'Cu' and (atom.position[2] <= 14.44000 or atom.position[2] >= 35.5420):
        species[i] = 'Au_1'
    elif atom.symbol == 'Cu':
        species[i] = 'Au_2'

# Convert to Pandas dataframe so we can print custom species
df1 = pd.DataFrame({"Species":species,
                    "X":junction.positions[:, 0],
                    "Y":junction.positions[:, 1],
                    "Z":junction.positions[:, 2]})
num_atoms = len(species)
filename_output = "{}/{}_junction.xyz".format(folder, polymorph)
print_xyz.print_from_pandas(df1, num_atoms, filename_output, save_dp='%.3f')
view(junction)



