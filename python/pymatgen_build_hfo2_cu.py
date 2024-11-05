from ase import Atoms
from ase.visualize import view
from ase.io import read, write
from ase.build import fcc100, fcc111
import numpy as np

folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/structures/pymatgen'
polymorph = 't-hfo2'
hafnia_supercell = read("{}/{}_supercell.xyz".format(folder, polymorph), format='xyz')
print(hafnia_supercell)

cu_slab = fcc100('Cu', size=(4, 4, 3), a=3.61)
print(cu_slab)

cu_o = 5.57326 - 3.73217
print(cu_o)

cell_size = np.array([2.55000 + 7.65000, 2.60000+7.80000, 21.10218])
z_cut = 3

junction = Atoms(cell=cell_size)

for atom in cu_slab:
    atom.position[0] = atom.position[0] + 1.27500
    atom.position[1] = atom.position[1] + ( 2.60000 -  1.27633)
    junction.append(atom)

# Move Cu slab so that x coordinate of outer Cu is same as O
hafnia_min = 3.82500
hafnia_max = 15.65109
cu_max = 3.61000
for atom in hafnia_supercell:
    if atom.position[2] > z_cut:
        atom.position[2] = atom.position[2] - hafnia_min + cu_max + cu_o
        junction.append(atom)

z_cut = 20
for atom in cu_slab:
    # atom.position[0] = atom.position[0] + 1.27500
    atom.position[2] = atom.position[2] + hafnia_max + cu_o
    if atom.position[2] < z_cut:
        junction.append(atom)

# Save the supercell to an XYZ file
write("{}/{}_junction.xyz".format(folder, polymorph), junction)
view(junction)

print(cell_size)

# # Lattice constant for copper
# a = 3.61
#
# # Create the bulk copper structure with FCC lattice
# cu_bulk = Structure(
#     Lattice.cubic(a),
#     ["Cu", "Cu", "Cu", "Cu"],
#     [[0, 0, 0], [0, 0.5, 0.5],
#      [0.5, 0, 0.5], [0.5, 0.5, 0]],
# )
# #
# # Use SlabGenerator to create a (100) slab
# slabgen = SlabGenerator(
#     initial_structure=cu_bulk,
#     miller_index=(1, 0, 0),
#     min_slab_size=2 * a,   # This is 1.5a for 3 layers
#     min_vacuum_size=10,           # Vacuum layer thickness
#     center_slab=True
# )
#
# # Generate the slab
# slab = slabgen.get_slabs(bonds=None, max_broken_bonds=0, symmetrize=True)[0]
# slab.make_supercell([4, 4, 1])  # Expands slab to 4x4 in a and b, and 1 in c
#
# # Convert Pymatgen structure to ASE Atoms object
# ase_atoms = Atoms(
#     symbols=[str(s) for s in slab.species],  # Convert Element objects to strings
#     positions=slab.cart_coords,
#     cell=slab.lattice.matrix,
#     pbc=True
# )
#
# # Visualize the structure using ASE
# view(ase_atoms)

# # Use Materials Project API to get HfO2 structure
# api_key = np.loadtxt('mp_api_key', dtype=str)
# print(api_key)
# mpr = MPRester(api_key)
#
# # Retrieve the HfO2 structure
# structure = mpr.get_structure_by_material_id("mp-1018721")  # mp-2023 is a known ID for tetragonal HfO2
# print(structure)
#
# # Convert Pymatgen structure to ASE Atoms object
# ase_atoms = Atoms(
#     symbols=[str(s) for s in structure.species],  # Convert Element objects to strings
#     positions=structure.cart_coords,
#     cell=structure.lattice.matrix,
#     pbc=True
# )
#
# # Save the unit cell to an XYZ file
# xyz_file_path = "HfO2_unit_cell.xyz"  # Specify your desired file name
# write(xyz_file_path, ase_atoms)
#
# # Optionally, visualize the structure using ASE
# view(ase_atoms)

# Define lattice parameters for tetragonal HfO2
# a = 5.1
# b = a
# c = 5.2
#
# # Create the lattice
# lattice = Lattice.from_parameters(a, b, c, 90, 90, 90)  # Tetragonal cell
#
# # Define the atomic positions for HfO2 in P42/nmc
# sites = [
#     ("Hf", [0, 0, 0]),              # Hafnium
#     ("O", [0.25, 0.25, 0.25]),      # Oxygen
#     ("O", [0.75, 0.75, 0.75]),      # Equivalent position of O
#     ("O", [0.25, 0.75, 0.75]),      # Equivalent position of O
#     ("O", [0.75, 0.25, 0.25]),      # Equivalent position of O
# ]
#
# # Create the structure
# structure = Structure(lattice, [site[0] for site in sites], [site[1] for site in sites])
#
# # Convert Pymatgen structure to ASE Atoms object
# ase_atoms = Atoms(
#     symbols=[str(s) for s in structure.species],  # Convert Element objects to strings
#     positions=structure.cart_coords,
#     cell=structure.lattice.matrix,
#     pbc=True
# )
#
# # Save the unit cell to an XYZ file
# # xyz_file_path = "HfO2_unit_cell_xyz"  # Specify your desired file name
# # write(xyz_file_path, ase_atoms)
#
# # Optionally, visualize the structure using ASE
# view(ase_atoms)



# # Generate a 2x2x2 supercell
# supercell = structure.copy()  # Create a copy of the structure
# supercell.make_supercell([2, 2, 2])  # Expand to 2x2x2
#
# # Convert Pymatgen structure to ASE Atoms object
# ase_atoms = Atoms(
#     symbols=[str(s) for s in supercell.species],  # Convert Element objects to strings
#     positions=supercell.cart_coords,
#     cell=supercell.lattice.matrix,
#     pbc=True
# )
#
# # Save the supercell to an XYZ file
# xyz_file_path = "HfO2_supercell.xyz"  # Specify your desired file name
# write(xyz_file_path, ase_atoms)
#
# # Visualize the structure using ASE
# view(ase_atoms)
