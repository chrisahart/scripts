from DPCharge_single_task import DP
from ase.io import read
import numpy as np
from ase import Atoms
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.md.npt import NPT
from ase.md.verlet import VelocityVerlet
from ase.md.langevin import Langevin
from ase import units


def identify_polaron_site(atoms: Atoms, spin: np.ndarray) -> np.ndarray:
    spin = spin.flatten()
    charge_state = np.zeros(len(atoms))

    # Identify the atom with the minimum spin

    allowed_element = ['O']  # only consider Fe atoms

    allowed_element_idx = [atom.index for atom in atoms if atom.symbol in allowed_element]
    spin_element = spin[allowed_element_idx]
    min_idx = np.argmax(spin_element)
    charge_state[allowed_element_idx[min_idx]] = 1

    return charge_state


# Read the atomic structure from a file
data = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/mgo/deepmd/cell-222/electron-u-6/database_spin_train/set.000'
coord = np.load('{}/coord.npy'.format(data))
box = np.load('{}/box.npy'.format(data))
vel = np.loadtxt('{}/vel.xyz'.format(data))
aparam = np.load('{}/aparam.npy'.format(data))

# MD simulation parameters
dt = 0.5 * units.fs
nstep = int(1e4)
thermo_freq = 10
T = 600
num_atoms = 64
symbols = ['Mg'] * 32 + ['O'] * 32
polaron_index_start = np.argmax(aparam[:num_atoms])

atoms = Atoms(symbols=symbols, positions=np.reshape(coord[0].flatten(), (num_atoms, 3)))

cell = np.reshape(box[0], (3, 3))
atoms.set_cell(cell)
atoms.set_pbc(True)

# Convert velocities from CP2K units (Å/fs) to ASE units (Å/ps)
vel_ase = np.reshape(vel.flatten(), (num_atoms, 3)) * units.Angstrom / units.fs
atoms.set_velocities(vel_ase)

# CP2K velocity angstrom*fs^-1
# atoms.set_velocities(np.reshape(vel.flatten(), (num_atoms, 3)))

md_type = 'nvt'
charge_state = np.zeros(len(atoms))
charge_state[polaron_index_start] = 1

# Create the DP calculator with two separate models

calc = DP(
    model_ener_force='model_ener.pth',
    model_spin='model_spin.pth',
    initial_charge_state=charge_state,
    identifier=identify_polaron_site

)
atoms.calc = calc

# Initialize velocities
# MaxwellBoltzmannDistribution(atoms, temperature_K=T)

# Set up and run the MD simulation

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

# Run the MD simulation

md.run(nstep)

# Save spin and charge state history

spin_history = calc.spin_history.copy()[2:]
charge_state_history = calc.charge_state_history.copy()
np.save('spin_history.npy', spin_history)
np.save('charge_state_history.npy', charge_state_history)