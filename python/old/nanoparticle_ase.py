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

# https://wiki.fysik.dtu.dk/ase/ase/cluster/cluster.html

# Create nanoparticle with only (100) and (111) surfaces, inspired by Clotilde nanoparticle
folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/theory_support/hannah/structures/nanoparticle-ase'

# lc = 3.61000  # Cu
lc = 3.94  # Pt

surfaces = [(1, 0, 0),  (1, 1, 1)]
# layers = [4, 3] # Pt 201
layers = [8, 6]  # Pt 1289
atoms = ase.cluster.cubic.FaceCenteredCubic('Pt', surfaces, layers,  latticeconstant=lc)
num_atoms = atoms.get_global_number_of_atoms()
print('num_atoms', num_atoms)

# Save
write('{}/struct-pt_atoms-{}.xyz'.format(folder, num_atoms), atoms)

# Orientation can be fixed using align tool in Avogadro to align along y and z avis
# Unfortunately Avogadro crashes if too many atoms

# Align atom positions manually
atoms_align_1 = np.array([598, 1283]) - 1  # Align along x and z
atoms_align_2 = np.array([1181, 492]) - 1  # Align along x and y
rotate_align = np.zeros(3)
converged = False

# Exhaustive 3D grid search
error = 1
grid_size = int(1e2)
rotate_x = np.linspace(start=-360, stop=360, num=grid_size)
rotate_y = np.linspace(start=-360, stop=360, num=grid_size)
rotate_z = np.linspace(start=-360, stop=360, num=grid_size)
# error = 1.0, grid_size = 1e2, angle = -287.27272727272725 69.09090909090907 149.09090909090907
# error = 0.9, grid_size = 1e2, angle = -185.45454545454547 -170.9090909090909 163.63636363636363
# error = 0.8, grid_size = 1e2, angle = -185.45454545454547 -170.9090909090909 163.63636363636363
# error = 0.7, grid_size = 1e2, angle = -185.45454545454547 -170.9090909090909 163.63636363636363
# error = 0.6, grid_size = 1e2, angle = -83.63636363636363 185.45454545454538 -214.54545454545456
# error = 0.5, grid_size = 1e2, angle = 207.27272727272725 -83.63636363636363 -3.636363636363626

# error = 0.1, grid_size = 1e3, angle = -355.6756756756757 -121.44144144144144 261.9819819819819
# error = 0.1, grid_size = 1e4, angle = -360.0 -289.5049504950495 -11.053105310531066

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
        for k in range(0, np.shape(rotate_z)[0]):

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
                write('{}/struct-pt_atoms-{}_aligned_grid_size-{}.xyz'.format(folder, num_atoms, grid_size), temp)

if not converged:
    print('Loop did not converge.')

if __name__ == "__main__":
    # view(atoms)
    # plt.show()
    print('Finished.')
