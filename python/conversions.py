from __future__ import division, print_function
from general import parameters as param

x = 0.707

print('Hartree/Angstrom to cm-1', x*param.rydberg_to_hartree)

# print('Ry to Hartree', x*param.rydberg_to_hartree)
# print('Ry/Bohr to eV/Angstrom', x*param.rydberg_to_hartree*param.hartree_per_bohr_to_ev_per_angstrom)
