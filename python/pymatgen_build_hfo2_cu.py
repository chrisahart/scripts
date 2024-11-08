from ase import Atoms
from ase.visualize import view
from ase.io import read, write
from ase.build import fcc100, fcc111
import numpy as np

folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/structures/pymatgen'
polymorph = 't-hfo2'
cu_lattice = 3.61
# cu_lattice = 4.17
cu_o = 5.57326 - 3.73217
z_cut = 3
hafnia_min = 3.82500

hafnia_supercell = read("{}/{}_supercell.xyz".format(folder, polymorph), format='xyz')
cu_slab = fcc100('Cu', size=(4, 4, 7), a=cu_lattice, periodic=True)
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

# Wrap the positions within the new cell
junction.wrap()

sorted_indices = np.lexsort((junction.positions[:, 1], junction.positions[:, 0], junction.positions[:, 2]))
junction = junction[sorted_indices]

print(10.25796-(2.60000-(10.25796-7.70531 )))

print(cu_lattice*3)
print('x', 7.65797+(7.65797-5.10531))
print('y', 7.65797+(7.65797-5.10531))

print(cu_lattice/2)
print(3.61000/2)

print(7.65797 +2.55266)

print(33.73718+ 1.80500)

# Remove Hf and O to form capacitor
# Remove all atoms with elements Hf or O
# mask = np.isin(junction.symbols, ['Hf', 'O'], invert=True)
# junction = Atoms(symbols=junction.symbols[mask],
#                  positions=junction.positions[mask])

write("{}/{}_junction_au_111.xyz".format(folder, polymorph), junction)
# write("{}/{}_junction_capacitor.xyz".format(folder, polymorph), junction)
view(junction)



