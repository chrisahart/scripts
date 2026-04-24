import parmed as pmd

folder = '/Users/chris/Documents/Storage/calculations/h2o/packmol/cubic/64/testing/xyz'

# Load PDB
struct = pmd.load_file('{}/system_precise.pdb'.format(folder))

# Set box
struct.box = [12.4, 12.4, 12.4, 90.0, 90.0, 90.0]

# Write PSF
struct.save('{}/system.psf'.format(folder), overwrite=True)

print(f"Atoms  : {len(struct.atoms)}")
print(f"Bonds  : {len(struct.bonds)}")
print(f"Waters : {len(struct.residues)}")