from pymatgen.core import Lattice, Structure, Element
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.io.ase import AseAtomsAdaptor
from ase.visualize import view

# Define lattice parameters for monoclinic HfO2
a, b, c = 5.09, 5.15, 5.278
alpha, beta, gamma = 90, 99.56, 90

# Create the lattice
lattice = Lattice.from_parameters(a, b, c, alpha, beta, gamma)

# Define atomic positions and species for the conventional cell
# These positions are based on the 4e Wyckoff positions (4 for Hf and 4 for O)
atomic_positions = [
    [0.723864, 0.541935, 0.292093],  # Hf at 4e
    [0.552002, 0.257474, 0.021422],  # O at 4e
    [0.931639, 0.669241, 0.653747],  # O at 4e
]

# The complete conventional unit cell for P2₁/c has specific positions
# based on the symmetry operations. For HfO2 with 12 atoms:
wyckoff_positions = [
    # Hf Wyckoff position (4e) replicated based on symmetry
    (0.723864, 0.541935, 0.292093),
    (0.276136, 0.458065, 0.207907),  # Equivalent by symmetry
    (0.223864, 0.541935, 0.292093),  # Equivalent by symmetry
    (0.776136, 0.458065, 0.207907),  # Equivalent by symmetry

    # O Wyckoff positions (4e)
    (0.552002, 0.257474, 0.021422),
    (0.447998, 0.742526, 0.021422),  # Equivalent by symmetry
    (0.931639, 0.669241, 0.653747),
    (0.068361, 0.330759, 0.346253),  # Equivalent by symmetry
]

# Define species corresponding to the positions
species = [Element("Hf"), Element("Hf"), Element("Hf"), Element("Hf"),
           Element("O"), Element("O"), Element("O"), Element("O")]

# Create the initial structure
initial_structure = Structure(lattice, species, wyckoff_positions)

# Analyze the space group to get the conventional unit cell
spacegroup_analyzer = SpacegroupAnalyzer(initial_structure, symprec=0.01)
conventional_structure = spacegroup_analyzer.get_conventional_standard_structure()

# Convert to ASE Atoms object for visualization
ase_atoms = AseAtomsAdaptor.get_atoms(conventional_structure)

# Visualize the conventional unit cell
view(ase_atoms, repeat=(1, 1, 1))

# Print out atomic sites to verify there are 12
print(f"Number of atomic sites in the conventional unit cell: {len(conventional_structure.sites)}")
for site in conventional_structure.sites:
    print(f"{site.species_string} at {site.frac_coords}")
