from __future__ import division, print_function
import pandas as pd
import numpy as np
from general import parameters as param
import matplotlib.pyplot as plt

"""
    HfO2. 3 atoms per formula unit, 4 formula units per cell
"""

# cp2k
# m = np.array([-326.278013863037472]) * param.hartree_to_ev
# po = np.array([-326.267945906153273]) * param.hartree_to_ev
# t = np.array([-326.257945121013336]) * param.hartree_to_ev

m = np.array([-326.278013892778802]) * param.hartree_to_ev
po = np.array([-326.267946120097974]) * param.hartree_to_ev
t = np.array([-326.257818085635222]) * param.hartree_to_ev
# print('PO eV per unit cell:', (m - po))
# print('PO meV per atom:', (m - po) * 1e3 / 12)
# print('T eV per unit cell:', (m-t))
# print('T meV per atom:', (m-t)*1e3/12)

print('PO meV per atom:', 0.317312*1e3/12)
print('T meV per atom:', 0.561317*1e3/12)


# qe Yudi
m = np.array([-717.64750235]) * param.rydberg_to_hartree * param.hartree_to_ev
po = np.array([-717.62479461]) * param.rydberg_to_hartree * param.hartree_to_ev
t = np.array([-717.60129133]) * param.rydberg_to_hartree * param.hartree_to_ev
print('PO eV per unit cell:', (m - po))
print('PO meV per atom:', (m - po) * 1e3 / 12)
print('T eV per unit cell:', (m-t))
print('T meV per atom:', (m-t)*1e3/12)

# [2] Fan
print('PO meV per atom:', (0.47*1e3/12))
print('T meV per atom:', (0.50*1e3/12))

# [3] Liu
print('PO eV per unit cell:', (128-75)*4/1e3)
print('PO meV per atom:', (128-75)/3)
print('T eV per unit cell:', 128*4/1e3)
print('T meV per atom:', 128/3)

# [3] Ma
print('PO eV per unit cell:', 84.3*4/1e3)
print('PO meV per atom:', 84.3/3)
print('T eV per unit cell:', 166.3*4/1e3)
print('T meV per atom:', 166.3/3)
