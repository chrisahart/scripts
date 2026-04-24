from ase.io import read, write
from ase import Atoms
import numpy as np
import MDAnalysis as mda

# folder = '/Users/chris/Documents/Storage/calculations/h2o/packmol/cubic/64/testing/xyz'
# input_file = '{}/system.xyz'.format(folder)
# output_file = '{}/system_h2o.xyz'.format(folder)

# folder = '/Users/chris/Documents/Storage/calculations/h2o/packmol/cubic/64/md/tip3p/testing/xyz2'
folder = '/Users/chris/Documents/Storage/calculations/h2o/packmol/cubic/64/md/tip3p/equil/temp-300'
input_file = '{}/last_center.xyz'.format(folder)
output_file = '{}/last_center_h2o.xyz'.format(folder)

# folder = '/Users/chris/Documents/Storage/calculations/h2o/packmol/cubic/128'
# folder = '/Users/chris/Documents/Storage/calculations/h2o/packmol/tio2_water/100'
# folder = '/Users/chris/Documents/Storage/calculations/h2o/packmol/tio2_water/001'
# folder = '/Users/chris/Documents/Storage/calculations/h2o/packmol/tio2_water/011'
# folder = '/Users/chris/Documents/Storage/calculations/h2o/packmol/tio2_water/110'
# input_file = '{}/system.xyz'.format(folder)
# output_file = '{}/system_h2o.xyz'.format(folder)

# Read the file
with open(input_file, 'r') as f:
    lines = f.readlines()

# Add H2O to the end of each coordinate line
with open(output_file, 'w') as f:
    f.write(lines[0])  # Number of atoms
    f.write(lines[1])  # Comment line
    for line in lines[2:]:
        f.write(line.strip() + " H2O\n")

print(f"Done! Modified file: {output_file}")