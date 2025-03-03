from deepmd.infer import DeepPot
import dpdata
import numpy as np

# 1. Load the database
# folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/chain/frozen-fe/multi-task'
# database = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/chain/frozen-fe/multi-task/database_ener_force'

# 1. Load the database
system = dpdata.LabeledSystem("../database", fmt='deepmd/npy')
coords = system['coords']
cells = system['cells']
atomic_types = system['atom_types']

polaron_atom = 1

nframes = coords.shape[0]
natoms = coords.shape[1]

# 2. Load the model
# model = 'model.ckpt.pt'
model = 'model.ckpt-400000.pt'
model_ener = DeepPot(model, head='ener_force')
model_spin = DeepPot(model, head='spin')

# 3. Test the model

# 3.1 aparam: all zeros
aparam = np.zeros((nframes, natoms, 1))
_, _, _, spin, _ = model_spin.eval(
        coords = coords,
        cells = cells,
        atom_types=atomic_types,
        atomic = True,
        aparam = aparam
    )
ener, force, _, _, _, = model_ener.eval(
        coords = coords,
        cells = cells,
        atom_types=atomic_types,
        atomic = True,
        aparam = aparam
    )
np.save('0_spin.npy', spin)
np.save('0_ener.npy', ener)
np.save('0_force.npy', force)


# 3.2 aparam: index 1 is 1, others are zeros
aparam = np.zeros((nframes, natoms, 1))
aparam[:, 0, 0] = polaron_atom
_, _, _, spin, _ = model_spin.eval(
        coords = coords,
        cells = cells,
        atom_types=atomic_types,
        atomic = True,
        aparam = aparam
    )
ener, force, _, _, _, = model_ener.eval(
        coords = coords,
        cells = cells,
        atom_types=atomic_types,
        atomic = True,
        aparam = aparam
    )
np.save('1_spin.npy', spin)
np.save('1_ener.npy', ener)
np.save('1_force.npy', force)