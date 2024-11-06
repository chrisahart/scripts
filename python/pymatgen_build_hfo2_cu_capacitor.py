from ase import Atoms
from ase.visualize import view
from ase.io import read, write
from ase.build import fcc100, fcc111
import numpy as np
import pandas as pd
from general import print_xyz

folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/structures/pymatgen'
polymorph = 't-hfo2'
cu_lattice = 3.61
print(cu_lattice*np.sqrt(2))
print(cu_lattice*np.sqrt(2)*2)
print(5.10500000*2)
# cu_lattice = 4.17971
cu_o = 5.57326 - 3.73217
z_cut = 3
hafnia_min = 3.82500

hafnia_supercell = read("{}/{}_supercell.xyz".format(folder, polymorph), format='xyz')
cu_slab = fcc100('Cu', size=(4, 4, 11), a=cu_lattice, periodic=True)
# cu_slab = fcc111('Cu', size=(4, 4, 7), a=cu_lattice, periodic=True)
print(cu_slab.cell)
junction = Atoms(pbc=True)

cu_max = np.max(cu_slab.positions[:, 2])

for atom in cu_slab:
    atom.position[0] = atom.position[0] + 1.27500
    atom.position[1] = atom.position[1] + (2.60000 - 1.27633)
    junction.append(atom)

for atom in hafnia_supercell:
    if atom.position[2] > z_cut:
        atom.position[2] = atom.position[2] - hafnia_min + cu_max + cu_o
        junction.append(atom)

z_cut = np.max(cu_slab.positions[:, 2]) - 1
junction_max = np.max(junction.positions[:, 2])

for atom in cu_slab:
    if atom.position[2] < z_cut:
        atom.position[2] = atom.position[2] + junction_max + cu_o
        junction.append(atom)

cell_size = np.array([2.55000 + 7.65000,
                      2.60000 + 7.80000,
                      np.max(junction.positions[:, 2]) + cu_lattice / 2])
cell_size = np.array([cu_lattice*3,
                      cu_lattice*3,
                      np.max(junction.positions[:, 2]) + cu_lattice / 2])
junction.set_cell(cell_size)
print(cell_size)

# Wrap the positions within the new cell and sort
junction.wrap()
sorted_indices = np.lexsort((junction.positions[:, 1], junction.positions[:, 0], junction.positions[:, 2]))
junction = junction[sorted_indices]

# Label each according to z index

# Remove Hf and O to form capacitor
# mask = np.isin(junction.symbols, ['Hf', 'O'], invert=True)
# junction = Atoms(symbols=junction.symbols[mask],
#                  positions=junction.positions[mask])

# Round all atom positions to 3 dp
for atom in junction:
    atom.position = np.round(atom.position, 3)

# Round all atom positions to 3 dp
# unique_z = np.unique(junction.positions[:, 2])
# species = []
# j = 1
# for i in range(unique_z.shape[0]):
#     for atom in junction:
#         if atom.symbol == 'Cu':
#             if atom.position[2] == unique_z[0]:
#                 species.append('Cu_{}'.format(str(j)))
#             j=j+1
#

# Need copy symbols from atoms object and then change species au to au_1 etc.
# Clean up file and

# Convert to Pandas dataframe so we can print custom species
# df1 = pd.DataFrame({"X":junction.positions[:, 0],
#                     "Y":junction.positions[:, 1],
#                     "Z":junction.positions[:, 2]})
# df1.insert(loc=0, column='Species', value=pd.Series(species).values)
# print(df1)
# num_atoms = len(species)
# filename_output = "{}/{}_junction_10layers.xyz".format(folder, polymorph)
# # filename_output = "{}/{}_junction_capacitor_10layers.xyz".format(folder, polymorph)
# print_xyz.print_from_pandas(df1, num_atoms, filename_output, save_dp='%.8f')

write("{}/{}_junction.xyz".format(folder, polymorph), junction)
# write("{}/{}_junction_capacitor_10layers.xyz".format(folder, polymorph), junction)
view(junction)



