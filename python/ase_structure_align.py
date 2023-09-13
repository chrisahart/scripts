import ase
import numpy as np
import matplotlib.pyplot as plt
import math
from ase.cluster.cubic import FaceCenteredCubic, SimpleCubic, BodyCenteredCubic
from ase.io import read, write
from ase import visualize
from ase.visualize import view
from ase.cluster.octahedron import Octahedron
from ase.cluster.decahedron import Decahedron
from ase.io import write
import ase
from ase.io import read, write

"""

"""

folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/theory_support/hannah/structures/nanoparticle-ase'
# atoms = ase.io.read("{}/struct-pt_atoms-1289_rotate.xyz".format(folder))
atoms = ase.io.read("{}/struct-pt_atoms-1289_aligned_grid_size-100.xyz".format(folder))

# Align atom positions manually
atoms_align_1 = np.array([598, 1283]) - 1  # Align along x and z
atoms_align_2 = np.array([1181, 492]) - 1  # Align along x and y
rotate_align = np.zeros(3)
converged = False

# Exhaustive 3D grid search
error = 0.1
grid_size = int(1e3)
rotate_x = np.linspace(start=-2, stop=2, num=grid_size)
rotate_y = np.linspace(start=-1, stop=1, num=grid_size)
rotate_z = np.linspace(start=-1, stop=1, num=grid_size)

# error 0.1, grid_size 1000,

# Change species of aligned atoms
atoms.symbols[atoms_align_1[0]] = 'Au'
atoms.symbols[atoms_align_1[1]] = 'Au'
atoms.symbols[atoms_align_2[0]] = 'C'
atoms.symbols[atoms_align_2[1]] = 'C'

print('error', error)
print('grid_size', grid_size)
print('atom pair 1', atoms.positions[atoms_align_1[0]], atoms.positions[atoms_align_1[1]])
print('atom pair 2', atoms.positions[atoms_align_2[0]], atoms.positions[atoms_align_2[1]])
for i in range(0, np.shape(rotate_x)[0]):
    print('iteration', i, 'of', np.shape(rotate_x)[0])
    print('angle_x', rotate_x[i])
    print('converged', converged)
    if converged: break
    for j in range(0, np.shape(rotate_y)[0]):
        if converged: break
        for k in range(0, np.shape(rotate_z)[0]):
            if converged: break

            angle_x = rotate_x[i]
            angle_y = rotate_y[j]
            angle_z = rotate_z[k]

            # Rotate structure
            temp = atoms.copy()
            temp.rotate(angle_x, 'x')
            temp.rotate(angle_y, 'y')
            temp.rotate(angle_z, 'z')

            # Get positions and their difference
            pos_diff_align_1 = temp.positions[atoms_align_1[0]] - temp.positions[atoms_align_1[1]]
            pos_diff_align_2 = temp.positions[atoms_align_2[0]] - temp.positions[atoms_align_2[1]]

            # Check if difference in x coordinate is zero and value of y coordinate is equal 0
            if np.abs(pos_diff_align_1[0]) < error and \
                    np.abs(pos_diff_align_1[2]) < error and \
                    np.abs(pos_diff_align_2[0]) < error and \
                    np.abs(pos_diff_align_2[1]) < error:
                rotate_align = np.array([angle_x, angle_y, angle_z])
                converged = True
                print('Position difference atom pair 1 x', pos_diff_align_1[0])
                print('Position difference atom pair 1 z', pos_diff_align_1[2])
                print('Position difference atom pair 2 x', pos_diff_align_2[0])
                print('Position difference atom pair 2 y', pos_diff_align_2[1])
                print('Angle = ', angle_x, angle_y, angle_z)

                # View in yz plane
                view(temp)
                write('{}/testing2.xyz'.format(folder), temp)
                # write('{}/struct-pt_atoms-{}_aligned_grid_size-{}.xyz'.format(folder, num_atoms, grid_size), temp)

if not converged:
    print('Loop did not converge.')

if __name__ == "__main__":
    # view(atoms)
    # plt.show()
    print('Finished.')
