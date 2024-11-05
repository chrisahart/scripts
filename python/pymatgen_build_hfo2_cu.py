from pymatgen.core.structure import Structure
from pymatgen.core.lattice import Lattice
from pymatgen.core.surface import SlabGenerator
from ase.build import bulk
from ase import Atoms
from ase.visualize import view
from mp_api.client import MPRester
from ase.io import write

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

# # Use Materials Project API to get HfO2 structure
# api_key = "YklFNMvJH9AjJIx00MJcrFEsl4Etc70n"  # Replace with your actual API key
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
a = 5.1  # Lattice parameter a (Å)
c = 5.2  # Lattice parameter c (Å)

# Create the lattice
lattice = Lattice.from_parameters(a, a, c, 90, 90, 90)  # Tetragonal cell

# Define the atomic positions for HfO2 in P42/nmc
# Hf at (0, 0, 0)
# O at (0.25, 0.25, 0.25) and its symmetry equivalents
sites = [
    ("Hf", [0, 0, 0]),              # Hafnium
    ("O", [0.25, 0.25, 0.25]),      # Oxygen
    ("O", [0.75, 0.75, 0.75]),      # Equivalent position of O
    ("O", [0.25, 0.75, 0.75]),      # Equivalent position of O
    ("O", [0.75, 0.25, 0.25]),      # Equivalent position of O
]

# Create the structure
structure = Structure(lattice, [site[0] for site in sites], [site[1] for site in sites])

# Convert Pymatgen structure to ASE Atoms object
ase_atoms = Atoms(
    symbols=[str(s) for s in structure.species],  # Convert Element objects to strings
    positions=structure.cart_coords,
    cell=structure.lattice.matrix,
    pbc=True
)

# Save the unit cell to an XYZ file
# xyz_file_path = "HfO2_unit_cell_xyz"  # Specify your desired file name
# write(xyz_file_path, ase_atoms)

# Optionally, visualize the structure using ASE
view(ase_atoms)



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

