import MDAnalysis as mda
import numpy as np

# Hafnia
# box = [15.12, 15.12, 9.62, 90, 90, 90]
# folder = '/Volumes/Samsung/Data/Postdoc2/Data/Work/calculations/hfo2/polaron/mengyu/chris/333/cdft/structures/0-1'
# folder = '/Volumes/Samsung/Data/Postdoc2/Data/Work/calculations/hfo2/polaron/mengyu/chris/333/cdft/structures/0-2a'
# folder = '/Volumes/Samsung/Data/Postdoc2/Data/Work/calculations/hfo2/polaron/mengyu/chris/333/cdft/structures/0-2b'
# folder = '/Volumes/Samsung/Data/Postdoc2/Data/Work/calculations/hfo2/polaron/mengyu/chris/333/cdft/structures/0-3a'
# folder = '/Volumes/Samsung/Data/Postdoc2/Data/Work/calculations/hfo2/polaron/mengyu/chris/333/cdft/structures/0-3b'

# Rutile 336
box = [13.77, 13.77, 17.76, 90, 90, 90]
# folder = '/Volumes/Samsung/Data/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-336/cdft/atoms-78-79/hse-22-dz/structures2'
folder = '/Volumes/Samsung/Data/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-336/cdft/atoms-79-101/hse-22-dz/structures2'
# folder = '/Volumes/Samsung/Data/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-336/cdft/atoms-79-107/hse-22-dz/structures2'
# box = [13.77, 13.77, 23.68, 90, 90, 90]
# folder = '/Volumes/Samsung/Data/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-338/cdft/atom-106-107/hse-22-dz/structures2'
# box = [13.77, 13.77, 29.60, 90, 90, 90]
# folder = '/Volumes/Samsung/Data/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-3310/cdft/atom-134-135/hse-22-dz/structures2'
# box = [13.77, 13.77, 35.52, 90, 90, 90]
# folder = '/Volumes/Samsung/Data/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-3312/cdft/atom-162-163/hse-22-dz/structures2'

# Anatase 441
# box = [15.08, 15.08, 9.68, 90, 90, 90]
# folder = '/Volumes/Samsung/Data/Postdoc2/Data/Work/calculations/tio2/anatase/archer/anatase/cell-441/cdft/atoms-122-129/structures2'
# folder = '/Volumes/Samsung/Data/Postdoc2/Data/Work/calculations/tio2/anatase/archer/anatase/cell-441/cdft/atoms-122-147/structures2'

# Anatase 442
# box = [15.08, 15.08, 17.76, 90, 90, 90]
# box = [15.08, 15.08, 19.36, 90, 90, 90]
# folder = '/Volumes/Samsung/Data/Postdoc2/Data/Work/calculations/tio2/anatase/archer/anatase/cell-442/cdft-cell-opt-hse-20/atoms-282-241/hse-19-dz/structures2'
# folder = '/Volumes/Samsung/Data/Postdoc2/Data/Work/calculations/tio2/anatase/archer/anatase/cell-442/cdft-cell-opt-hse-20/atoms-282-305/hse-19-dz/structures2'
# folder = '/Volumes/Samsung/Data/Postdoc2/Data/Work/calculations/tio2/anatase/archer/anatase/cell-442/cdft-cell-opt-hse-20/atoms-282-310/hse-19-dz/structures2'
# folder = '/Volumes/Samsung/Data/Postdoc2/Data/Work/calculations/tio2/anatase/archer/anatase/cell-442/cdft-cell-opt-hse-20/atoms-282-310/hse-19-dz/structures2/test'
# folder = '/Volumes/Samsung/Data/Postdoc2/Data/Work/calculations/tio2/anatase/archer/anatase/cell-442/cdft-cell-opt-hse-20/atoms-282-340/hse-19-dz/structures2'
# folder = '/Volumes/Samsung/Data/Postdoc2/Data/Work/calculations/tio2/anatase/archer/anatase/cell-442/cdft-cell-opt-hse-20/atoms-282-361/hse-19-dz/structures2'
# folder = '/Volumes/Samsung/Data/Postdoc2/Data/Work/calculations/tio2/anatase/archer/anatase/cell-442/cdft-cell-opt-hse-20/atoms-282-367/hse-19-dz/structures2'

# Load the initial and final structures
initial = mda.Universe(f'{folder}/initial.xyz')
final = mda.Universe(f'{folder}/final.xyz')

# Set the box dimensions for both universes
initial.dimensions = box
final.dimensions = box

# Create a new universe for the interpolated structure
ts = mda.Merge(initial.atoms)

# Set the box dimensions for the interpolated structure
ts.dimensions = box

# Calculate the displacement vector under PBC
displacement = final.atoms.positions - initial.atoms.positions
displacement -= np.round(displacement / np.array(box[:3])) * np.array(box[:3])

# Interpolate the coordinates (alpha = 0.5 for midpoint)
alpha = 0.5
ts.atoms.positions = initial.atoms.positions + alpha * displacement

# Manually write the .xyz file without the MDAnalysis header
with open(f'{folder}/ts.xyz', 'w') as f:
    f.write(f"{len(ts.atoms)}\n")  # Number of atoms
    f.write("\n")  # Comment line (can be customized)
    for atom in ts.atoms:
        f.write(f"{atom.name:2s} {atom.position[0]:12.6f} {atom.position[1]:12.6f} {atom.position[2]:12.6f}\n")