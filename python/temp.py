import ase
from ase.spacegroup import crystal

"""
    Temp
"""

# Create an Ag(110)-Si(110) interface with three atomic layers
# on each side.
a_ag = 4.09
ag = crystal(['Ag'], basis=[(0,0,0)], spacegroup=225,
             cellpar=[a_ag, a_ag, a_ag, 90., 90., 90.])
ag110 = ase.build.cut(ag, (0, 0, 3), (-1.5, 1.5, 0), nlayers=3)

a_si = 5.43
si = crystal(['Si'], basis=[(0,0,0)], spacegroup=227,
             cellpar=[a_si, a_si, a_si, 90., 90., 90.])
si110 = ase.build.cut(si, (0, 0, 2), (-1, 1, 0), nlayers=3)