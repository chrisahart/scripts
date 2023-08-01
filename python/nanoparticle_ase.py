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
layers = [8, 6]
atoms = ase.cluster.cubic.FaceCenteredCubic('Pt', surfaces, layers,  latticeconstant=lc)
print(atoms.get_global_number_of_atoms())

# Save
write('{}/struct-pt_larger.xyz'.format(folder), atoms)

# Attempt to rotate
# atoms.rotate(-130, 'x')
# atoms.rotate(-54, 'y')
# atoms.rotate(-86, 'z')
# atoms.rotate(-90, 'x')
# atoms.rotate(-90, 'y')
# atoms.rotate(-49.417141333078206, 'x')
# atoms.rotate(-53.65189183391107, 'y')
# atoms.rotate(-85.93881580223298, 'z')
# atoms.rotate(-90, 'x')

# Rotation can be fixed using align tool in Avogadro to align along y and z avis
# Unfortunately Avogadro crashes if too many atoms

# Align atom positions manually
atoms_align_z = np.array([598, 1283]) - 1
atoms_align_y = np.array([1181, 492]) - 1
atoms_align_x = np.array([598, 1283]) - 1
rotate_align = np.zeros(3)
converged = False

# Coarse search
error = 0.4
rotate_x = np.linspace(start=-360, stop=360, num=int(1e2))
rotate_y = np.linspace(start=-360, stop=360, num=int(1e2))
rotate_z = np.linspace(start=-360, stop=360, num=int(1e2))
# error = 1,  angle: -287.27272727272725 69.09090909090907 149.09090909090907
# error = 0.9,  angle: -185.45454545454547 -170.9090909090909 163.63636363636363
# error = 0.8,  angle: -185.45454545454547 -170.9090909090909 163.63636363636363
# error = 0.7,  angle: -185.45454545454547 -170.9090909090909 163.63636363636363
# error = 0.6,  angle: -83.63636363636363 185.45454545454538 -214.54545454545456
# error = 0.5,  angle: 207.27272727272725 -83.63636363636363 -3.636363636363626

# error = 0.1
# rotate_x = np.linspace(start=-360, stop=360, num=int(1e4))
# rotate_y = np.linspace(start=-360, stop=360, num=int(1e4))
# rotate_z = np.linspace(start=-360, stop=360, num=int(1e4))

# Fine search
# error = 0.1
# rotate_x = np.linspace(start=-186, stop=-185, num=int(1e2))
# rotate_y = np.linspace(start=-171, stop=170, num=int(1e2))
# rotate_z = np.linspace(start=163, stop=164, num=int(1e2))
# error = 0.5,  angle: -185.95959595959596 -112.44444444444444 163.2828282828283

error = 0.5
rotate_x = np.linspace(start=-84, stop=-83, num=int(1e2))
rotate_y = np.linspace(start=185, stop=186, num=int(1e2))
rotate_z = np.linspace(start=-215, stop=-214, num=int(1e2))
# 207.27272727272725 -83.63636363636363 -3.636363636363626
print('error', error)

# error = 0.1
# rotate_x = np.linspace(start=-360, stop=360, num=int(1e2))
# rotate_y = np.linspace(start=-360, stop=360, num=int(1e2))
# rotate_z = np.linspace(start=-360, stop=360, num=int(1e2))

print(atoms.positions[atoms_align_x[0]])
print(atoms.positions[atoms_align_x[1]])

# Change species of aligned atoms
atoms.symbols[atoms_align_x[0]] = 'Au'
atoms.symbols[atoms_align_x[1]] = 'Au'
atoms.symbols[atoms_align_y[0]] = 'C'
atoms.symbols[atoms_align_y[1]] = 'C'
# print('species atom 1', atoms.symbols[atoms_align_1[0]])
# print('species atom 2', atoms.symbols[atoms_align_1[1]])

for i in range(0, np.shape(rotate_x)[0]):
    if not converged:
        print('iteration', i, 'of', np.shape(rotate_x)[0])
        print('angle_x', rotate_x[i])

        for j in range(0, np.shape(rotate_y)[0]):
            for k in range(0, np.shape(rotate_z)[0]):
                angle_x = rotate_x[i]
                angle_y = rotate_y[j]
                angle_z = rotate_z[k]

                # Rotate structure
                atoms.rotate(angle_x, 'x')
                atoms.rotate(angle_y, 'y')
                atoms.rotate(angle_z, 'z')

                # Get positions and their difference
                pos_1_align_x = atoms.positions[atoms_align_x[0]]
                pos_2_align_x = atoms.positions[atoms_align_x[1]]
                pos_diff_align_x = pos_1_align_x - pos_2_align_x
                pos_1_align_y = atoms.positions[atoms_align_y[0]]
                pos_2_align_y = atoms.positions[atoms_align_y[1]]
                pos_diff_align_y = pos_1_align_y - pos_2_align_y
                pos_1_align_z = atoms.positions[atoms_align_z[0]]
                pos_2_align_z = atoms.positions[atoms_align_z[1]]
                pos_diff_align_z = pos_1_align_z - pos_2_align_z

                # Check if difference in x coordinate is zero and value of y coordinate is equal 0
                if np.abs(pos_diff_align_x[0]) < error and \
                        np.abs(pos_diff_align_y[0]) < error and \
                        np.abs(pos_diff_align_y[1]) < error and \
                        np.abs(pos_diff_align_z[2]) < error:

                    rotate_align = np.array([angle_x, angle_y, angle_z])
                    converged = True
                    print('pos_diff_align_x[0]', pos_diff_align_x[0])
                    print('pos_diff_align_y[0]', pos_diff_align_y[0])
                    print('pos_diff_align_y[0]', pos_diff_align_y[1])
                    print('pos_diff_align_y[0]', pos_diff_align_z[2])
                    print('angle:', angle_x, angle_y, angle_z)

                    # View in yz plane
                    view(atoms)
                    write('{}/struct-pt_larger_rotate.xyz'.format(folder), atoms)

                # Reset positions
                atoms.rotate(-angle_x, 'x')
                atoms.rotate(-angle_y, 'y')
                atoms.rotate(-angle_z, 'z')


# rotate_align = np.array([-360.0, -18.181818181818187, -170.9090909090909])
# print(rotate_align[0], rotate_align[1], rotate_align[2])
# atoms.rotate(rotate_align[0], 'x')
# atoms.rotate(rotate_align[1], 'y')
# atoms.rotate(rotate_align[2], 'z')
# print(atoms.positions[atoms_align_x[0]])
# print(atoms.positions[atoms_align_x[1]])

# Save
# write('{}/struct-pt_rotate.xyz'.format(folder), atoms)
# write('{}/struct-pt_larger_rotate.xyz'.format(folder), atoms)

if __name__ == "__main__":
    # view(atoms)
    # plt.show()
    print('Finished.')

