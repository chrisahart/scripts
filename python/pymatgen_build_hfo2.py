from pymatgen.core.structure import Structure
from pymatgen.core.lattice import Lattice
from pymatgen.core.surface import SlabGenerator
from ase.build import bulk
from ase import Atoms
from ase.visualize import view
from mp_api.client import MPRester
from ase.io import write
import numpy as np
from pymatgen.io.ase import AseAtomsAdaptor

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
#
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

polymorph = 't-hfo2'
polymorph_id = 'mp-1018721'

# # Use Materials Project API to get HfO2 structure
# api_key = str(np.loadtxt('mp_api_key', dtype=str))
# mpr = MPRester(api_key)
# unit_cell = mpr.get_structure_by_material_id(polymorph_id)
# print(unit_cell)

# Primitive cell for monoclinic has 12 atoms, for tetragonal has 6 so use a sqrt(2) doubled cell
a = 5.1
b = a
c = 5.2
lattice = Lattice.from_parameters(a, b, c, 90, 90, 90)  # Tetragonal cell
sites = [
    ("Hf", [0.000, 0.000, 0.000]),
    ("Hf", [0.000, 0.500, 0.500]),
    ("Hf", [0.500, 0.000, 0.500]),
    ("Hf", [0.500, 0.500, 0.000]),
    ("O", [0.250, 0.250, 0.206]),
    ("O", [0.750, 0.750, 0.206]),
    ("O", [0.750, 0.250, 0.794]),
    ("O", [0.250, 0.750, 0.794]),
    ("O", [0.250, 0.250, 0.706]),
    ("O", [0.750, 0.750, 0.706]),
    ("O", [0.750, 0.250, 0.294]),
    ("O", [0.250, 0.750, 0.294]),
]
unit_cell = Structure(lattice, [site[0] for site in sites], [site[1] for site in sites])

supercell = unit_cell.copy()
supercell.make_supercell([3, 3, 3])

# Switch the y and z axes
rotation_matrix = np.array([[1, 0, 0],  # x remains the same
                             [0, 0, 1],  # y becomes z
                             [0, 1, 0]])  # z becomes y

# Apply rotation to all atomic coordinates in the supercell
for site in supercell:
    rotated_coords = rotation_matrix @ site.coords
    site.coords = rotated_coords

# Create the (111) surface cut
# Adjust min_slab_size and min_vacuum_size as needed
min_slab_size = 3  # in Å, thickness of the slab
min_vacuum_size = 15  # in Å, thickness of the vacuum layer
slabgen = SlabGenerator(unit_cell, miller_index=(1, 0, 0), min_slab_size=1,
                        min_vacuum_size=0)
slabs = slabgen.get_slabs()
if slabs:
    surface_slab = slabs[0]
else:
    print("Failed to generate (111) slab.")
    exit()
surface_slab = surface_slab * (2,2,2)

ase_unit_cell = AseAtomsAdaptor.get_atoms(unit_cell)
ase_supercell = AseAtomsAdaptor.get_atoms(supercell)
ase_supercell_111 = AseAtomsAdaptor.get_atoms(surface_slab)

# Save the unit cell to an XYZ file
folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/structures/pymatgen'
write("{}/{}_unit_cell.xyz".format(folder, polymorph), ase_unit_cell)
write("{}/{}_supercell.xyz".format(folder, polymorph), ase_supercell)
write("{}/{}_supercell_111.xyz".format(folder, polymorph), ase_supercell_111)

view(ase_unit_cell)
view(ase_supercell)
view(ase_supercell_111)

# Define lattice parameters for tetragonal HfO2 P42/nmc
# a = 5.1
# b = a
# c = 5.2
# lattice = Lattice.from_parameters(a, b, c, 90, 90, 90)  # Tetragonal cell
# sites = [
#     ("Hf", [0, 0, 0]),              # Hafnium
#     ("O", [0.25, 0.25, 0.25]),      # Oxygen
#     ("O", [0.75, 0.75, 0.75]),      # Equivalent position of O
#     ("O", [0.25, 0.75, 0.75]),      # Equivalent position of O
#     ("O", [0.75, 0.25, 0.25]),      # Equivalent position of O
# ]
#
# structure = Structure(lattice, [site[0] for site in sites], [site[1] for site in sites])
# ase_atoms = Atoms(
#     symbols=[str(s) for s in structure.species],  # Convert Element objects to strings
#     positions=structure.cart_coords,
#     cell=structure.lattice.matrix,
#     pbc=True
# )

# Save the unit cell to an XYZ file
# xyz_file_path = "HfO2_unit_cell_xyz"  # Specify your desired file name
# write(xyz_file_path, ase_atoms)

# Optionally, visualize the structure using ASE
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

