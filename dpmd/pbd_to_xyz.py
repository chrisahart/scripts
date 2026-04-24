from ase.io import read, write
from ase import Atoms
import numpy as np
import MDAnalysis as mda

folder = '/Users/chris/Documents/Storage/calculations/h2o/packmol/cubic/64/testing/xyz'

u_pdb = mda.Universe('{}/system.pdb'.format(folder))
u_xyz = mda.Universe('{}/system.xyz'.format(folder))

u_pdb.atoms.positions = u_xyz.atoms.positions
u_pdb.dimensions = [12.4, 12.4, 12.4, 90.0, 90.0, 90.0]

# Fix residue names
for res in u_pdb.residues:
    res.resname = 'WAT'

# Fix atom names
for res in u_pdb.residues:
    idx = res.atoms.indices
    res.atoms[0].name = 'O'
    res.atoms[1].name = 'H1'
    res.atoms[2].name = 'H2'

# Manual bonds
bond_list = []
for res in u_pdb.residues:
    idx = res.atoms.indices
    bond_list.append((idx[0], idx[1]))
    bond_list.append((idx[0], idx[2]))

u_pdb.add_TopologyAttr('bonds', bond_list)
u_pdb.atoms.write('{}/system_precise.pdb'.format(folder))

# Verify
print(f"Atoms  : {len(u_pdb.atoms)}")
print(f"Bonds  : {len(u_pdb.bonds)}")
print(f"Waters : {len(u_pdb.residues)}")

# Verify bond lengths
lengths = np.array([
    np.linalg.norm(u_pdb.atoms[b[0]].position - u_pdb.atoms[b[1]].position)
    for b in bond_list
])
print(f"Min O-H: {lengths.min():.6f} Å")
print(f"Max O-H: {lengths.max():.6f} Å")
