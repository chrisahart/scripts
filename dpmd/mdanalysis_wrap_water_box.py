from ase.io import read, write
import numpy as np

folder = '/Volumes/Elements/Data/Postdoc2/Data/Work/calculations/h2o/packmol/tio2_water/110/md/equil/temp-300'
input_file = f"{folder}/last.xyz"
output_file = f"{folder}/last_center_fixed2.xyz"

# folder = '/Volumes/Elements/Data/Postdoc2/Data/Work/calculations/h2o/packmol/tio2_water/110/md/equil/pbe/from-tip3p-1-ns'
# input_file = f"{folder}/last.xyz"
# output_file = f"{folder}/last_center_fixed.xyz"

# Box dimensions
Lx, Ly, Lz = 12.9824805, 17.76, 16.6100492

# Load water
water = read(input_file)

n_atoms = len(water)
assert n_atoms % 3 == 0, "Number of atoms not divisible by 3. Check water box."
n_molecules = n_atoms // 3

# Wrap each H2O molecule by its center-of-mass
for i in range(n_molecules):
    indices = [3*i, 3*i+1, 3*i+2]  # O,H,H
    coords = water.positions[indices]

    # Compute molecule COM
    com = coords.mean(axis=0)

    # Shift COM into [0,L) for each axis
    shift = np.array([0.0,0.0,0.0])
    shift[0] = -np.floor(com[0]/Lx)*Lx
    shift[1] = -np.floor(com[1]/Ly)*Ly
    shift[2] = -np.floor(com[2]/Lz)*Lz

    # Apply shift to all atoms in molecule
    water.positions[indices] = coords + shift

# Save wrapped water
write(output_file, water)
print("Finished: water molecules wrapped along x, y, z and intact.")