import numpy as np
import matplotlib.pyplot as plt
from scripts.main import parameters as param
import ase
from mayavi import mlab

from ase.data import covalent_radii
from ase.io.cube import read_cube_data
from ase.data.colors import cpk_colors
from ase.calculators.calculator import get_calculator_class

""" Plotting of SMEAGOL output _TRC.agr by filename"""

fileobj = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/iv/li/cp2k/transmission/dev-chris/single-points/V-0_HLB-F_z-0-0/0V-ELECTRON_DENSITY-1_0.cube'
test = ase.io.cube.read_cube(fileobj, read_data=True, program=None, verbose=False)