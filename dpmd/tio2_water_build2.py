from ase.spacegroup import crystal
from ase.build import surface
from ase.io import write
import numpy as np


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
# Create 3x3x6 supercell
supercell = rutile.repeat((3, 3, 3))

# Build a thick slab (no vacuum yet)
slab = surface(rutile, (1, 0, 0), layers=12, vacuum=0.0)

# --- Step 1: identify atomic planes along z ---
z = slab.positions[:, 2]
symbols = np.array(slab.get_chemical_symbols())

# Cluster atoms into planes (tolerance ~0.2 Å)
planes = []
used = np.zeros(len(z), dtype=bool)

for i in np.argsort(z):
    if used[i]:
        continue
    plane = np.where(abs(z - z[i]) < 0.2)[0]
    used[plane] = True
    planes.append(plane)

# Sort planes from bottom → top
planes = sorted(planes, key=lambda p: z[p].mean())

# --- Step 2: classify each plane ---
plane_types = []
for p in planes:
    elems = set(symbols[p])
    plane_types.append(elems)

# --- Step 3: find a window with O termination on BOTH ends ---
start, end = None, None

for i in range(len(planes)):
    if plane_types[i] == {'O'}:  # bottom O layer
        for j in range(i + 4, len(planes)):  # ensure thickness
            if plane_types[j] == {'O'}:  # top O layer
                start, end = i, j
                break
    if start is not None:
        break

if start is None:
    raise RuntimeError("Could not find O–O terminated slab window.")

# --- Step 4: extract slab ---
selected_indices = np.concatenate(planes[start:end+1])
slab = slab[selected_indices]

# --- Step 5: center + vacuum ---
slab.center(vacuum=15.0, axis=2)

# --- Optional: lateral expansion ---
slab = slab.repeat((3, 3, 1))

# Save to .xyz file
write('{}/rutile_3x3x3.xyz'.format(folder), supercell)
write('{}/rutile_3x3x3_001.xyz'.format(folder), slab)
