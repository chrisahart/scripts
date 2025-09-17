import numpy as np
import matplotlib.pyplot as plt
from ase.io import read, write, Trajectory
import MDAnalysis as mda
from MDAnalysis.analysis import distances
from MDAnalysis.analysis import rdf
# from MDAnalysis.analysis import vacf
from ase import Atoms
import ase.io
from general import parameters as param


def compute_vacf(vels):
    # Flatten atom and axis
    v = vels.reshape(vels.shape[0], -1)
    v -= v.mean(axis=0)
    nframes = v.shape[0]
    vacf = np.zeros(nframes)
    for t in range(nframes):
        vacf[t] = np.mean(np.sum(v[:nframes-t] * v[t:], axis=1))
    return vacf


# MgO
# folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/mgo/deepmd/md/cell-222/electron-u-6/temp-600-velocity-dft/single-task-dpa3'
# folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/mgo/deepmd/md/cell-222/electron-u-6/temp-600-velocity-dft-population2/single-task-dpa3-new'
# num_atoms = 64
# box_size = [8.38, 8.38, 8.38, 90, 90, 90]

# Rutile 336 22%
# folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/deepmd/rutile/336/md-cell-opt/deepmd-md/hse-22-v2/single-fit-ener-dpa3-nlayers-6-new2'
# folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/deepmd/rutile/336/md-cell-opt/deepmd-md/hse-22-v2/single-fit-ener-dpa3-nlayers-6-new2'
# num_atoms = 324
# box_size = [13.77, 13.77, 17.76, 90, 90, 90]

# Rutile 336 PBE + U
topology_file = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/deepmd/rutile/336/md-cell-opt/system.xyz'
# folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/deepmd/rutile/336/md-cell-opt/deepmd-md/pbe-u-4-ps-v2/single-fit-ener-se_e2_a-official-v3.1.0-sel-90'
# folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/deepmd/rutile/336/md-cell-opt/deepmd-md/pbe-u-4-ps-v2/single-fit-ener-dpa3-nlayers-6-new2-v3.1.0'
# folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/deepmd/rutile/336/md-cell-opt/deepmd-md/pbe-u-4-ps-v2/old/single-fit-m-dpa3-nlayers-6-new2-spin-v3'
# folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/deepmd/rutile/336/md-cell-opt/deepmd-md/pbe-u-4-ps-v2/single-fit-spin-se_e2_a-official-v3.1.0'
# folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/deepmd/rutile/336/md-cell-opt/deepmd-md/pbe-u-4-ps-v2/single-fit-spin-dpa3-nlayers-6-new2-v3.1.0-test-2-md'
folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/deepmd/rutile/336/md-cell-opt/deepmd-md/pbe-u-4-ps-v2/single-fit-m-dpa3-lr-single-fit-ener-dpa3-md'

num_atoms = 324
box_size = [13.77, 13.77, 17.76, 90, 90, 90]

print(folder)
spin = np.load("{}/spin_history.npy".format(folder))
# print(spin)
charge_state = np.load("{}/charge_state_history.npy".format(folder))
# print(charge_state[0])
# print(charge_state[-1])
calc_distance = True
# calc_distance = False

with Trajectory('{}/md.traj'.format(folder)) as traj:
    nframes = len(traj)
    for i, atoms in enumerate(traj):
        write('{}/tio2-pos-1.xyz'.format(folder), atoms, format='xyz', append=i > 0)
        v_atoms = Atoms(symbols=atoms.get_chemical_symbols(), positions=atoms.get_velocities())
        write('{}/tio2-vel-1.xyz'.format(folder), v_atoms, format='xyz', append=i > 0)

trajectory_file_pos = '{}/tio2-pos-1.xyz'.format(folder)

# all_atoms = list(ase.io.read('{}/tio2-vel-1.xyz'.format(folder), index=':'))
# n_atoms = 324
# n_frames = len(all_atoms)
# velocities = np.zeros((n_frames, n_atoms, 3))
# for i, atoms in enumerate(all_atoms):
#     velocities[i] = atoms.positions.copy()
# universe_velocity = mda.Universe(v_atoms)

num_timesteps = spin.shape[0]
print(num_timesteps)
time_array = np.linspace(0, int(num_timesteps), num=num_timesteps)
timestep = 1
local_bonds = 6

draw_legend = False
plotting_colors = ['r', 'g', 'b', 'm', 'grey', 'orange', 'brown', 'hotpink'] * 100
# plotting_colors = ['b', 'g', 'r', 'm', 'grey', 'orange', 'brown', 'hotpink'] * 100
# xlim_1 = [0, (5103+200)/2]
xlim_1 = [0, time_array[-1]]
# ylim_1 = [0, 0.9]
ylim_1 = [-0.02, 0.9]

# Setup md analysis environment
if calc_distance:
    universe = mda.Universe(topology_file, trajectory_file_pos)
    universe.dimensions = box_size
    num_timesteps2 = len(universe.trajectory)
    print(num_timesteps2)
    xlim_2 = [0, num_timesteps2]
    time_val_1 = np.linspace(start=0, stop=len(universe.trajectory) * timestep, num=len(universe.trajectory))
    atoms_ti = universe.select_atoms('name Ti')
    num_atoms_ti = len(atoms_ti)
    atoms_o = universe.select_atoms('name O')
    num_atoms_o = len(atoms_o)
    dist_arr = distances.distance_array(atoms_ti.positions, atoms_o.positions)
    bond_lengths_time = np.zeros((num_timesteps2, num_atoms_ti, num_atoms_o))
    for ts in universe.trajectory:
        frame = universe.trajectory.frame
        bond_lengths_time[frame] = distances.distance_array(atoms_ti.positions, atoms_o.positions, box=box_size)
    bond_lengths_time_sorted = np.zeros((num_timesteps2, num_atoms_ti, local_bonds))
    bond_lengths_time_sorted_mean = np.zeros((num_timesteps2, num_atoms_ti))
    for i in range(num_atoms_ti):
        for j in range(num_timesteps2):
            bond_lengths_time_sorted[j, i] = np.sort(bond_lengths_time[j, i])[0:local_bonds]
            # bond_lengths_time_sorted_mean[j, i] = np.mean(bond_lengths_time_sorted[j, i])
    bond_lengths_time_sorted_mean = np.mean(bond_lengths_time_sorted, axis=2)

if calc_distance:
    metric = np.zeros((num_atoms_ti, num_timesteps2))
    fig_bonds_1, ax_bonds_1 = plt.subplots(figsize=(10, 4))
    # fig_bonds_1, ax_bonds_1 = plt.subplots()
    for i in range(num_atoms_ti):
        # ax_bonds_1.plot(time_val_1 - time_val_1[0], np.sum(bond_lengths_time_sorted, axis=2)[:, i], '-', label='Fe {}'.format(i + 1))
        ax_bonds_1.plot(time_val_1 - time_val_1[0], bond_lengths_time_sorted_mean[:, i], '-', label='Fe {}'.format(i + 1))
    # ax_bonds_1.plot(time_val_1 - time_val_1[0], bond_lengths_time_sorted_mean[:, polaron_atom], 'k-', label='Fe {}'.format(polaron_atom + 1))
    ax_bonds_1.set_xlabel('Time / fs')
    # ax_bonds_1.set_xlabel('Timestep')
    ax_bonds_1.set_ylabel('Average of {} Ti-O bond lengths / A'.format(local_bonds))
    # # ax_bonds_1.set_xlim([0, len(universe.trajectory)])
    ax_bonds_1.set_xlim(xlim_2)
    # ax_bonds_1.set_xlim([0, len(universe.trajectory) * timestep])
    # ax_bonds_1.set_ylim([0.06, -0.10])
    # fig_bonds_1.savefig('{}/bond_lengths_average.png'.format(folder_save), dpi=300)
    # fig_bonds_1.tight_layout()

    print('np.mean(bond_lengths_time_sorted_mean)', np.mean(bond_lengths_time_sorted_mean))

# Calculate polaron atom
polaron_atom_time = np.zeros(num_timesteps, dtype=int)
for j in range(num_timesteps):
    polaron_atom_time[j] = int(np.argmax(spin[j, :]))
polaron_atoms = np.unique(polaron_atom_time)
# print('polaron_atoms', polaron_atoms)

# Plot total spin
fig_spin2, ax_spin2 = plt.subplots(figsize=(10, 4))
# fig_spin2, ax_spin2 = plt.subplots()
ax_spin2.plot(time_array, np.sum(spin, axis=1), 'k-')
ax_spin2.set_xlim(0, time_array[-1])
ax_spin2.set_xlabel("Time / fs")
ax_spin2.set_ylabel("Spin")
ax_spin2.set_xlim(xlim_2)
fig_spin2.tight_layout()
fig_spin2.savefig("{}/dp_md_spin_sum.png".format(folder), dpi=600)
ax_spin2.set_xlim(xlim_1)
fig_spin2.savefig("{}/dp_md_spin_sum_lim.png".format(folder), dpi=600)
fig_spin2.tight_layout()

# Plot all spin
fig_spin1, ax_spin1 = plt.subplots(figsize=(10, 4))
# fig_spin1, ax_spin1 = plt.subplots()
for i in range(num_atoms):
    ax_spin1.plot(time_array, spin[:, i], '-')
# for j in range(polaron_atoms.shape[0]):
#     print('(polaron_atoms[j]', polaron_atoms[j])
#     ax_spin1.plot(time_array, spin[:, polaron_atoms[j]], '-', color=plotting_colors[j], label='{}'.format(polaron_atoms[j]+1))
# ax_spin1.plot(time_array, spin[:, 52], '-', color=plotting_colors[0], label='{}'.format(52))
# ax_spin1.plot(time_array, np.sum(spin, axis=1), 'k-', alpha=0.2)
if draw_legend: ax_spin1.legend(frameon=True)
ax_spin1.set_xlim(0, time_array[-1])
ax_spin1.set_xlabel("Time / fs")
ax_spin1.set_ylabel("Spin")
ax_spin1.set_ylim(ylim_1)
ax_spin1.set_xlim(xlim_2)
fig_spin1.tight_layout()
fig_spin1.savefig("{}/dp_md_spin.png".format(folder), dpi=600)
ax_spin1.set_xlim(xlim_1)
fig_spin1.savefig("{}/dp_md_spin_lim.png".format(folder), dpi=600)
fig_spin1.tight_layout()

# Plot RDF for Ti - O
nbins = 300
rdf_xlim = 13.77 / 2

if calc_distance:
    ti_polaron = universe.select_atoms("bynum 78")
    atoms_o = universe.select_atoms("element O")
    rdf_ti_o_polaron = rdf.InterRDF(ti_polaron, atoms_o, range=(0, rdf_xlim), nbins=nbins)
    rdf_ti_o_polaron.run()

    ti_not_polaron = universe.select_atoms("element Ti and not bynum 78")
    rdf_ti_o_not_polaron = rdf.InterRDF(ti_not_polaron, atoms_o, range=(0, rdf_xlim), nbins=nbins)
    rdf_ti_o_not_polaron.run()

    fig_rdf_ti_o, ax_rdf_ti_o = plt.subplots(figsize=(10, 4))
    ax_rdf_ti_o.plot(rdf_ti_o_polaron.bins, rdf_ti_o_polaron.rdf, 'r-', label='RDF Ti (polaron) - O')
    ax_rdf_ti_o.plot(rdf_ti_o_not_polaron.bins, rdf_ti_o_not_polaron.rdf, 'g-', label='RDF Ti - O')
    ax_rdf_ti_o.set_xlabel("Radial distance / Å")
    ax_rdf_ti_o.set_ylabel("RDF (arb. units)")
    ax_rdf_ti_o.legend(frameon=False)
    ax_rdf_ti_o.set_xlim([0, rdf_xlim])
    ax_rdf_ti_o.tick_params(axis='y', which='both', left=False, right=False, labelleft=False)
    fig_rdf_ti_o.tight_layout()
    fig_rdf_ti_o.savefig("{}/rdf.png".format(folder), dpi=600)

# Plot RDF for Ti - Ti
if calc_distance:
    ti_polaron = universe.select_atoms("bynum 78")
    rdf_ti_ti_polaron = rdf.InterRDF(ti_polaron, ti_not_polaron, range=(0, rdf_xlim), nbins=nbins)
    rdf_ti_ti_polaron.run()

    ti_not_polaron = universe.select_atoms("element Ti and not bynum 78")
    rdf_ti_ti_not_polaron = rdf.InterRDF(ti_not_polaron, ti_not_polaron, exclusion_block=(1, 1), range=(0, rdf_xlim), nbins=nbins)
    rdf_ti_ti_not_polaron.run()

    fig_rdf_ti_ti, ax_rdf_ti_ti = plt.subplots(figsize=(10, 4))
    ax_rdf_ti_ti.plot(rdf_ti_ti_polaron.bins, rdf_ti_ti_polaron.rdf, 'r-', label='RDF Ti (polaron) - Ti')
    ax_rdf_ti_ti.plot(rdf_ti_ti_not_polaron.bins, rdf_ti_ti_not_polaron.rdf, 'g-', label='RDF Ti - Ti')
    ax_rdf_ti_ti.set_xlabel("Radial distance / Å")
    ax_rdf_ti_ti.set_ylabel("RDF (arb. units)")
    ax_rdf_ti_ti.legend(frameon=False)
    ax_rdf_ti_ti.set_xlim([0, rdf_xlim])
    ax_rdf_ti_ti.tick_params(axis='y', which='both', left=False, right=False, labelleft=False)
    fig_rdf_ti_ti.tight_layout()
    fig_rdf_ti_ti.savefig("{}/rdf.png".format(folder), dpi=600)


if calc_distance:
    all_atoms = list(ase.io.read('{}/tio2-vel-1.xyz'.format(folder), index=':'))  # Each frame
    velocities = np.array([atoms.positions for atoms in all_atoms])  # shape: (n_frames, n_atoms, 3)
    vacf_all = compute_vacf(velocities)

    dt = 1.0
    N = len(vacf_all)
    freqs = np.fft.fftfreq(N, dt*1e-15)[:N//2]
    vdos_all = np.abs(np.fft.fft(vacf_all))[:N//2]
    c = 2.99792458e10
    freqs_cm1 = freqs / c
    freqs_thz = freqs * 1e-12

    # fig_vdos1, ax_vdos1 = plt.subplots(figsize=(10, 4))
    # ax_vdos1.plot(freqs_cm1, vdos_all, 'k-')
    # ax_vdos1.set_xlim([0, 1000])
    # ax_vdos1.set_xlabel('Frequency (1/cm)')
    # ax_vdos1.set_ylabel('VDOS (arb. units)')
    # fig_vdos1.tight_layout()
    # fig_vdos1.savefig("{}/vdos_cm.png".format(folder), dpi=600)

    fig_vdos2, ax_vdos2 = plt.subplots(figsize=(10, 4))
    ax_vdos2.plot(freqs_thz, vdos_all, 'k-')
    ax_vdos2.set_xlim([0, 32])
    ax_vdos2.set_xlabel('Frequency (THz)')
    ax_vdos2.set_ylabel('VDOS (arb. units)')
    fig_vdos2.tight_layout()
    # fig_vdos2.savefig("{}/vdos_thz.png".format(folder), dpi=600)

if __name__ == "__main__":
    print('Finished.')
    plt.show()
