from ase.io import read
import numpy as np
from ase import Atoms
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md.npt import NPT
from ase.md.verlet import VelocityVerlet
from ase.md.langevin import Langevin
from ase import units

data = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/mgo/deepmd/cell-222/electron-u-6/database_spin_train/set.000'
coord = np.load('{}/coord.npy'.format(data))
box = np.load('{}/box.npy'.format(data))
aparam = np.load('{}/aparam.npy'.format(data))
vel = np.loadtxt('{}/vel-test.xyz'.format(data))

coord_use = np.reshape(coord[0].flatten(), (64, 3))
symbols = ['Mg'] * 32 + ['O'] * 32

test = Atoms(symbols=symbols, positions=coord_use)
cell = np.reshape(box[0], (3, 3))
test.set_cell(cell)
test.set_pbc(True)

print(test)
print(test.get_number_of_atoms)
test.set_velocities(np.reshape(vel.flatten(), (64, 3)))
print(test)

print(aparam[:64])
print(np.argmax(aparam[:64]))
