import ase
from ase.io import read, write

"""
Convert VASP POSCAR to .xyz using ASE
"""

# folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Masters/2021-2022/Yike/Archive/YIKE/ORR/water/Pt3Ni/water_bilayer'
folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/masters/xiangzhuochun/pt3ni/yike/vacuum/Pt3Ni/O-f'
bec = ase.io.read("{}/POSCAR".format(folder))
bec_new = ase.io.write("{}/coord.xyz".format(folder), bec, format="xyz")
