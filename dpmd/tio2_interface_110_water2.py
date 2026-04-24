from ase.io import read, write
import numpy as np

folder = '/Volumes/Elements/Data/Postdoc2/Data/Work/calculations/tio2-h2o/ahart/110/monolayer/geo_opt/pbe-frozen-tio2'
tio2_file = '/Volumes/Elements/Data/Postdoc2/Data/Work/calculations/tio2-h2o/ahart/110/monolayer/geo_opt/pbe-frozen-tio2/system.xyz'
water_file = '/Volumes/Elements/Data/Postdoc2/Data/Work/calculations/h2o/packmol/tio2_water/110/md/equil/temp-300/last_center_fixed.xyz'
tio2_water_file = '/Volumes/Elements/Data/Postdoc2/Data/Work/calculations/tio2-h2o/ahart/110/monolayer/geo_opt/pbe-frozen-tio2/tio2_water2.xyz'
target_gap = 3.2

# load structures
slab = read(tio2_file)
water = read(water_file)

# ----------------------------
# 1. Identify slab top
# ----------------------------
symbols = np.array(slab.get_chemical_symbols())
z_positions = slab.positions[:, 2]

top_ti_z = z_positions[symbols == 'Ti'].max()

# ----------------------------
# 2. Water geometry reference
# ----------------------------
symbols_water = np.array(water.get_chemical_symbols())
water_z_min = water.positions[symbols_water == 'O', 2].min()

# ----------------------------
# 3. Compute rigid shift for water
#    (slab stays fixed, water moves once)
# ----------------------------
shift_z = (top_ti_z + target_gap) - water_z_min
water.positions[:, 2] += shift_z

# ----------------------------
# 4. Merge (STRICT ORDER: slab first, then water)
# ----------------------------
combined = slab + water

r0 = combined.positions[41-1].copy()
combined.translate(-r0)
combined.translate([0.5, 0.5, 0])

# ----------------------------
# 5. Define cell (no post-centering)
# ----------------------------
symbols_combined = np.array(combined.get_chemical_symbols())
combined_o_max = combined.positions[symbols_combined == 'O', 2].max()

Lx = 12.9824805
Ly = 17.76
Lz = combined_o_max + target_gap - (22.26844-21.00265)
# Lz = combined_o_max + 1

cell = np.array([Lx, Ly, Lz, 90, 90, 90])
print(cell)

combined.set_cell(cell)
combined.set_pbc([False, False, True])

# IMPORTANT: only ONE wrap, at the end (optional)
# combined.translate([0,0,10])
combined.wrap()

# ----------------------------
# 6. Save
# ----------------------------
write(tio2_water_file, combined, format='xyz')

print("Done: slab + water stacked (no re-centering, no density distortion).")