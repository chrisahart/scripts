import optparse

import numpy as np
from ase.visualize import view
from ase.data import covalent_radii
from ase.io.cube import read_cube_data
from ase.data.colors import cpk_colors
from ase.calculators.calculator import get_calculator_class


# eV/Ang^2
conv = 1e-8 #/ 1.6e-19
k1 = 4 * conv
k2 = 4 * conv
k = 7 * conv

f1 = (k1+k2)/2+k
f2 = np.sqrt((k1-k2)**2+4*(k**2))/2

w1 = np.sqrt(f1-f2)
w2 = np.sqrt(f1+f2)
print(w1)
print(w2)
print(w2/w1)
print(2230/1062)
