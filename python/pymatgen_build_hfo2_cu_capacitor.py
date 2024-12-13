from ase import Atoms
from ase.visualize import view
from ase.io import read, write
from ase.build import fcc100, fcc111
import numpy as np
import pandas as pd
from general import print_xyz

folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/structures/pymatgen/cu/supercell-1-1-4-cu-1.86'
# folder = r'C:\Users\storm\Desktop\old'
polymorph = 't-hfo2'
hfo2_lattice_a = 5.105
hfo2_lattice_a_primitive = np.round(hfo2_lattice_a/np.sqrt(2), 3)
cu_lattice = hfo2_lattice_a_primitive
cu_lattice_spacing_xy = cu_lattice*np.sqrt(2)/2
print('cu_lattice = hfo2_lattice_a/np.sqrt(2)', hfo2_lattice_a_primitive)
print('DFT optimised cu lattice / DFT optimised hafnia lattice * 100 - 100 = ', 3.650042/hfo2_lattice_a_primitive * 100 - 100)

cu_o = 1.84  # 5.57326 - 3.73217
cu_o = 2.00 - 0.156  # 5.57326 - 3.73217
print(cu_o)

print(23.50000-21.65900)

layers_xy = 2
layers_z = 13

hafnia_supercell = read("{}/{}_supercell.xyz".format(folder, polymorph), format='xyz')
cu_slab = fcc100('Cu', size=(6, 6, layers_z), a=cu_lattice, periodic=True)
junction = Atoms(pbc=True)

cu_max = np.max(cu_slab.positions[:, 2])
cu_xy_unique = np.unique(cu_slab.positions[:, 0])
cu_z_unique = np.unique(cu_slab.positions[:, 2])

hafnia_xy_unique = np.unique(hafnia_supercell.positions[:, 0])
cut_index = layers_xy * 2 - len(hafnia_xy_unique)
cut_index_cu = cut_index
hafnia_z_unique = np.unique(hafnia_supercell.positions[:, 2])
print('hafnia_z_unique', hafnia_z_unique)
z_cut_hafnia = hafnia_z_unique[2]

print(cu_xy_unique)
print(cu_xy_unique[cut_index_cu])

for atom in cu_slab:
    # junction.append(atom)
    if (atom.position[0] < cu_xy_unique[cut_index_cu]
            and atom.position[1] < cu_xy_unique[cut_index_cu]):
        junction.append(atom)

z_cut = np.max(cu_slab.positions[:, 2]) - 1
junction_max = hfo2_lattice_a/2 * layers_z + np.max(cu_slab.positions[:, 2]) + cu_o
print('junction_max', junction_max)
junction_max = 68.95000
junction_max = 48.75000



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

z_upper_hafnia = hafnia_z_unique[-5]
for atom in hafnia_supercell:
    # atom.position[2] = atom.position[2] - hafnia_z_unique[3] + cu_max + cu_o
    # atom.position[0] = atom.position[0] - 1.27500
    # atom.position[1] = atom.position[1] + (2.60000 - 1.27633)
    # junction.append(atom)
    # if atom.position[2] > z_cut_hafnia:
    #     atom.position[2] = atom.position[2] - hafnia_z_unique[3] + cu_max + cu_o
    #     atom.position[0] = atom.position[0] - 1.27500
    #     atom.position[1] = atom.position[1] + (2.60000 - 1.27633)
    #     junction.append(atom)
    if (atom.position[2] > z_cut_hafnia and atom.position[2] <= z_upper_hafnia
            and atom.position[0] < hafnia_xy_unique[cut_index]
            and atom.position[1] < hafnia_xy_unique[cut_index]):
        atom.position[2] = atom.position[2] - hafnia_z_unique[3] + cu_max + cu_o
        atom.position[0] = atom.position[0] - (3.78800 - 2.55200)
        atom.position[1] = atom.position[1] - (2.59500 -1.27600)
        junction.append(atom)

print(np.max(junction.positions[:, 2]))
test = np.max(junction.positions[:, 2])

for atom in cu_slab:
    # if atom.position[2] < z_cut:
    #     atom.position[2] = atom.position[2] + junction_max + cu_o
    #     junction.append(atom)
    if (atom.position[2] < z_cut and
            atom.position[0] < cu_xy_unique[cut_index_cu]
            and atom.position[1] < cu_xy_unique[cut_index_cu]):
        atom.position[2] = atom.position[2] + test + cu_o
        junction.append(atom)

for atom in junction:
     atom.position[1] = atom.position[1] + 0.25000

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
left_bulk = 10.830
left_au_1 = 14.44000
right_bulk = 54.45200
right_au_1 = 50.84200
# right_bulk = 44.25200
# right_au_1 =40.64200

left_bulk = np.unique(junction.positions[:, 2])[8]
print('left_bulk', left_bulk)
left_au_1 = np.unique(junction.positions[:, 2])[10]
print('left_au_1', left_au_1)
right_bulk = 58.08500
right_bulk = np.unique(junction.positions[:, 2])[-8]
print('right_bulk', right_bulk)
right_au_1 = 54.47600
right_au_1 = np.unique(junction.positions[:, 2])[-10]
print('right_au_1', right_au_1)
for i, atom in enumerate(junction):
    if atom.symbol == 'Cu' and (atom.position[2] <= left_bulk or atom.position[2] >= right_bulk):
        species[i] = 'Cu_bulk'
    elif atom.symbol == 'Cu' and (atom.position[2] <= left_au_1 or atom.position[2] >= right_au_1):
        species[i] = 'Cu_1'
    elif atom.symbol == 'Cu':
        species[i] = 'Cu_2'

# Convert to Pandas dataframe so we can print custom species
df1 = pd.DataFrame({"Species":species,
                    "X":junction.positions[:, 0],
                    "Y":junction.positions[:, 1],
                    "Z":junction.positions[:, 2]})
num_atoms = len(species)
filename_output = "{}/{}_junction.xyz".format(folder, polymorph)
print_xyz.print_from_pandas(df1, num_atoms, filename_output, save_dp='%.3f')
view(junction)



