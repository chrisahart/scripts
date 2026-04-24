from ase.io import write
from ase.build import surface
from ase.spacegroup import crystal
import numpy as np
from ase import Atom
from ase.neighborlist import NeighborList
from ase.build import molecule

folder = '/Users/chris/Documents/Storage/calculations/tio2-h2o/ahart'
a_hse = 4.59
c_hse = 2.96
u = 0.305

# Build unit cell (space group 136: P42/mnm)
rutile = crystal(
    symbols=['Ti', 'O'],
    basis=[(0, 0, 0), (u, u, 0)],
    spacegroup=136,
    cellpar=[a_hse, a_hse, c_hse, 90, 90, 90]
)
# Create 3x3x3 supercell
supercell = rutile.repeat((3, 3, 3))

# Create (1 0 0) surface slab of 3x3x3 supercell
slab = surface(rutile, (1, 1, 0), layers=5, vacuum=5)
slab = slab.repeat((2, 6, 1))

# Detect Ti layers along y, remove atoms below second Ti layer
sym = np.array(slab.get_chemical_symbols())
z = slab.get_positions()[:, 2]
ti_z = np.sort(np.unique(np.round(z[sym == "Ti"], 1)))
ti_filer = ti_z[0]
slab = slab[z > ti_filer]

#
# # Remove bottom row of Ti atoms
# sym = np.array(slab.get_chemical_symbols())
# z = slab.get_positions()[:, 2]
# ti_z = np.sort(np.unique(np.round(z[sym == "Ti"], 1)))
# slab = slab[~((sym == "Ti") & (np.abs(z - ti_z[0]) < 0.5))]
# slab = slab[~((np.array(slab.get_chemical_symbols()) == "Ti") & (np.abs(slab.get_positions()[:, 2] - ti_z[-1]) < 0.5))]
#
# slab.rotate(180, "z", rotate_cell=False)
# slab.rotate(180, "y", rotate_cell=False)
# slab.rotate(180, "z", rotate_cell=False)
# slab.set_center_of_mass(slab.get_cell().sum(axis=0) / 2)
# # slab.positions += slab.get_cell().diagonal()[:3] * [1, 1, 0]  # shift x,y back into cell
# # slab.wrap()
#

slab_water = slab.copy()

oh_bond = 0.97  # Å
cutoffs = [1.5 if s == "Ti" else 1.0 for s in slab_water.get_chemical_symbols()]
nl = NeighborList(cutoffs, self_interaction=False, bothways=True)
nl.update(slab_water)

sym = np.array(slab_water.get_chemical_symbols())
z = slab_water.get_positions()[:, 2]
z_top = z[sym == "O"].max()
z_bot = z[sym == "O"].min()

for i in np.where(sym == "O")[0]:
    indices, _ = nl.get_neighbors(i)
    ti_neighbors = np.sum(sym[indices] == "Ti")
    if ti_neighbors == 1:
        if np.abs(z[i] - z_top) < 0.5:
            slab_water.append(Atom("H", slab_water.get_positions()[i] + [0.79, 0, 0.62]))
            slab_water.append(Atom("H", slab_water.get_positions()[i] + [-0.79, 0, 0.62]))
        elif np.abs(z[i] - z_bot) < 0.5:
            slab_water.append(Atom("H", slab_water.get_positions()[i] + [0.79, 0, -0.62]))
            slab_water.append(Atom("H", slab_water.get_positions()[i] + [-0.79, 0, -0.62]))

# Set center of mass to center of supercell
slab.set_center_of_mass(slab.get_cell().sum(axis=0) / 2)
slab_water.set_center_of_mass(slab_water.get_cell().sum(axis=0) / 2)

# Order atoms along z axis
slab = slab[np.argsort(slab.get_positions()[:, 2])]
slab_water = slab_water[np.argsort(slab_water.get_positions()[:, 2])]

# Save to .xyz file
print("Slab cell (Å):", np.diag(slab.get_cell()))
print("Slab monolayer cell (Å):", np.diag(slab_water.get_cell()))
write('{}/rutile_3x3x3.xyz'.format(folder), supercell)
write('{}/rutile_3x3x3_110.xyz'.format(folder), slab, format='xyz')
write('{}/rutile_3x3x3_110_monolayer.xyz'.format(folder), slab_water, format='xyz')