import ase
from ase.io import read, write

"""
Convert VASP POSCAR to .xyz using ASE
"""

folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Files/mcc/YIKE/ORR/water/Pt3Ni/water_bilayer'
bec = ase.io.read("{}/POSCAR".format(folder))
bec_new = ase.io.write("{}/coord.xyz".format(folder), bec, format="xyz")
