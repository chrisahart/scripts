from ase.build import surface
from ase.spacegroup import crystal
from ase.io import write
import numpy as np

folder = '/Users/chris/Documents/Storage/calculations/tio2-h2o/ahart'
a_hse = 4.59
c_hse = 2.96
u = 0.305

# Build rutile TiO2 manually (space group 136)
rutile = crystal(
    symbols=["Ti", "O"],
    basis=[(0, 0, 0), (0.3, 0.3, 0)],
    spacegroup=136,
    cellpar=[a_hse, a_hse, c_hse, 90, 90, 90]
)
slab = surface(rutile, (1, 0, 0), layers=7, vacuum=10)
slab = slab.repeat((3, 3, 1))

# Remove topmost Ti layer to expose O termination
sym = np.array(slab.get_chemical_symbols())
z = slab.get_positions()[:, 2]
z_top_Ti = z[sym == "Ti"].max()
del slab[np.where((sym == "Ti") & (z > z_top_Ti - 0.5))[0]]

write('{}/rutile_3x3x3_001_claude.xyz'.format(folder), slab)
