from ase.io import read, write
from ase.build import surface
from ase.spacegroup import crystal
import numpy as np
from ase import Atom
from ase.neighborlist import NeighborList
from ase.build import molecule

tio2_file = '/Volumes/Elements/Data/Postdoc2/Data/Work/calculations/tio2-h2o/ahart/110/monolayer/geo_opt/pbe-frozen-tio2/system.xyz'
water_file = '/Volumes/Elements/Data/Postdoc2/Data/Work/calculations/h2o/packmol/tio2_water/110/md/equil/temp-300/last_center_fixed.xyz'
folder = '/Volumes/Elements/Data/Postdoc2/Data/Work/calculations/tio2-h2o/ahart/110/monolayer/geo_opt/pbe-frozen-tio2'

# cell is water box + tio2 box + 2 * (3.5 A - Ti-OH2 distance)
# ti_o = 11.71669 - 9.73686
Lz_water = 16.6100492  # Å
# cell = np.array([12.9824805,  17.76,  (24.96230878-10)+Lz_water+(3.5-ti_o), 90, 90, 90])
cell = np.array([12.9824805,  17.76,  9.73686+Lz_water+2*3.5, 90, 90, 90])

# Load structure
slab = read(tio2_file)
water = read(water_file)

# Get position of atom index 0
symbols = np.array(slab.get_chemical_symbols())
z_positions = slab.positions[:, 2]   # <-- define here
top_ti_z = z_positions[symbols == 'Ti'].max()

# --- Water center (robust, not atom-dependent) ---
z_center_water = water.positions[:, 2].mean()

# --- Compute bottom of water slab ---
water_bottom_z = z_center_water - 0.5 * Lz_water

# --- Desired placement ---
target_gap = 3.5
shift_z = (top_ti_z + target_gap) - water_bottom_z

# Apply translation (only z-direction)
water.positions[:, 2] += shift_z

# --- Step 5: merge systems ---
combined = slab + water
combined.positions[:, 2] -= 5
# combined.positions[:, 2] += Lz_water/2
combined.set_cell(cell)
combined.set_pbc([False, False, True])
combined.wrap()
combined.positions[:, 2] += Lz_water/2
combined.wrap()

# Order atoms along z axis
# combined = combined[np.argsort(combined.get_positions()[:, 2])]
print("Slab cell (Å):", np.diag(combined.get_cell()))

# Save if needed
write('{}/shifted_system.xyz'.format(folder), slab, format='xyz')
write('{}/tio2_water.xyz'.format(folder), combined, format='xyz')
