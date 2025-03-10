from DPCharge import DP
from ase.io import read
import numpy as np
from ase import Atoms
from ase.constraints import FixedLine, FixAtoms
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md.npt import NPT
from ase.md.verlet import VelocityVerlet
from ase.md.langevin import Langevin
from ase import units


def identify_polaron_site(atoms: Atoms, spin: np.ndarray) -> np.ndarray:
    spin = spin.flatten()
    charge_state = np.zeros(len(atoms))

    # identify the atom with the minimum spin
    allowed_element = ['Fe'] # only consider Fe atoms
    allowed_element_idx = [atom.index for atom in atoms if atom.symbol in allowed_element]
    spin_element = spin[allowed_element_idx]
    min_idx = np.argmin(spin_element)
    charge_state[allowed_element_idx[min_idx]] = 1

    return charge_state


folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/liushiLab/chain/hse/frozen-fe-h2o/hops-1/ase_md/chris'
atoms = read("{}/POSCAR".format(folder))

dt = 0.5*units.fs
nstep = int(1e4)
thermo_freq = 10
T = 300
md_type = 'nvt'
fixed_idx = list(range(7)) + list(range(35, 77))
c2 = FixAtoms(
    indices=fixed_idx
)
atoms.set_constraint(c2)
charge_state = np.zeros(len(atoms))
charge_state[2] = 1


calc = DP(
    model='model.ckpt.pt',
    force_head='ener_force',
    charge_head='spin',
    initial_charge_state=charge_state,
    identifier=identify_polaron_site
)
atoms.calc = calc
MaxwellBoltzmannDistribution(atoms, temperature_K=T)
if md_type == 'npt':
    md = NPT(
        atoms=atoms,
        timestep=dt,
        temperature_K=T,
        externalstress=0.0,
        ttime=100 * units.fs,
        pfactor=100 * units.fs,
        trajectory='md.traj',
        logfile='md.log',
        loginterval=thermo_freq,
        append_trajectory=False
    )
elif md_type == 'nve':
    md = VelocityVerlet(
        atoms=atoms,
        timestep=dt,
        trajectory='md.traj',
        logfile='md.log',
        loginterval=thermo_freq,
        append_trajectory=False
    )
elif md_type == 'nvt':
    md = Langevin(
        atoms=atoms,
        timestep=dt,
        temperature_K=T,
        friction=0.01 / units.fs,
        trajectory='md.traj',
        logfile='md.log',
        loginterval=thermo_freq,
        append_trajectory=False
    )

md.run(nstep)

# 3. save spin and charge_state history
spin_history = calc.spin_history.copy()[2:]
charge_state_history = calc.charge_state_history.copy()
np.save('spin_history.npy', spin_history)
np.save('charge_state_history.npy', charge_state_history)
