import numpy as np
from scipy import integrate
from matplotlib import pyplot as plt
from ase.io.cube import read_cube_data
from general import parameters as param

distance = np.array([0, 0.1, 0.2])
force = np.array([0.07134469, 0.06332366, 0.05514261])
energy = np.array([-0.908177903559681, -0.895657666594070, -0.884564823160404])

# /gpfs/home/cahart/transport/force_integral/h2

print('energy')
# print((energy[1] - energy[0]))
print((energy[1] - energy[0])*param.hartree_to_ev)

print('integrate.simpson', integrate.simpson(force[:-1], distance[:-1])*param.hartree_to_ev*2)

print('energy')
# print((energy[2] - energy[1]))
print((energy[2] - energy[1])*param.hartree_to_ev)

print('integrate.simpson', integrate.simpson(force[1:], distance[1:])*param.hartree_to_ev*2)