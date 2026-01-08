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
import os
from scipy.optimize import curve_fit


def compute_vacf(vels):
    # Flatten atom and axis
    v = vels.reshape(vels.shape[0], -1)
    v -= v.mean(axis=0)
    nframes = v.shape[0]
    vacf = np.zeros(nframes)
    for t in range(nframes):
        vacf[t] = np.mean(np.sum(v[:nframes-t] * v[t:], axis=1))
    return vacf


# Define the function to fit: y = mx
def linear_func(x, m):
    return m * x


# Define the function to fit: y = mx + c
def linear_func2(x, m, c):
    return m * x + c


# Bulk anatase hse-19-ts-md-9500-9900-removed
topology_file = '/Volumes/Samsung/Data/Postdoc2/Data/Work/calculations/tio2/anatase/archer/anatase/cell-441/md-cell-opt-hse-20/hse-19-complete/combined/system.xyz'
folder_1 = '/Volumes/Samsung/Data/Postdoc2/Data/Work/calculations/tio2/anatase/deepmd/anatase/441/deepmd/hse-19-ts-md-9500-9900-removed'
folder_energy = 'single-fit-ener-dpa3-nlayers-6-official-v3.1.0-start_pref-0.02-1000_limit_pref-1-1-rcut-4.5-twostep-lr-1e-5-1e-8'
folder_spin = 'single-fit-pop-dpa3-nlayers-6-official-v3.1.0-dev-polaron-loss-mae-pref-1-pref_pop-1000-1'
# folder_md = '600k-vel-17-ps'
# folder_md = '600k-vel-18-ps' #
# folder_md = '600k-vel-18-ps-50-ps'
# folder_md = '600k-vel-19-ps'
# temperature_set = 600

# folder_md = '200k-vel-18-ps-10-ps'
# folder_md = '200k-vel-18-ps-50-ps'
# folder_md = '200k-vel-18-ps-1000-ps'
# temperature_set = 200

# folder_md = '300k-vel-18-ps-10-ps'
# folder_md = '300k-vel-18-ps-50-ps'
# folder_md = '300k-vel-18-ps-1000-ps'
# temperature_set = 300

# folder_md = '400k-vel-18-ps-10-ps'
# folder_md = '400k-vel-18-ps-50-ps'
# folder_md = '400k-vel-18-ps-1000-ps'
# temperature_set = 400

# folder_md = '500k-vel-18-ps-10-ps'
# folder_md = '500k-vel-18-ps-50-ps'
# folder_md = '500k-vel-18-ps-1000-ps'
# temperature_set = 500

folder_md = '600k-vel-18-ps-10-ps'
# folder_md = '600k-vel-18-ps-50-ps'
# folder_md = '600k-vel-18-ps-1000-ps'
temperature_set = 600
# fit_start = 20000
fit_start = 5000
# fit_start = 0

folder = '{}/{}/{}/{}'.format(folder_1, folder_energy, folder_spin, folder_md)

num_atoms = 192
box_size = [15.08, 15.08, 9.68, 90, 90, 90]
save_fig = True
folder_save = folder

print(folder)
print(folder_md)
spin = np.load("{}/spin_history.npy".format(folder))
# print(spin)
charge_state = np.load("{}/charge_state_history.npy".format(folder))
# print(charge_state[0])
# print(charge_state[-1])
calc_distance = True
# calc_distance = False
plot_rdf = False

pos_file = '{}/tio2-pos-1.xyz'.format(folder)
vel_file = '{}/tio2-vel-1.xyz'.format(folder)

# Delete files if they exist
for file in [pos_file, vel_file]:
    if os.path.exists(file):
        os.remove(file)

# Write trajectory data
with Trajectory('{}/md.traj'.format(folder)) as traj:
    for i, atoms in enumerate(traj):
        if i == 0:
            continue
        write(pos_file, atoms, format='xyz', append=(i > 1))
        v_atoms = Atoms(symbols=atoms.get_chemical_symbols(), positions=atoms.get_velocities())
        write(vel_file, v_atoms, format='xyz', append=(i > 1))

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
# xlim_1 = [0, 1000]
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

if calc_distance:
    polaron_atom_time = np.zeros(num_timesteps, dtype=int)
    for j in range(num_timesteps):
        polaron_atom_time[j] = int(np.argmax(spin[j, :]))
    # polaron_atoms = np.unique(polaron_atom_time)
    polaron_atoms = polaron_atom_time[np.insert(polaron_atom_time[:-1] != polaron_atom_time[1:], 0, True)]
    print('polaron_atoms', polaron_atoms+1)

    # Calculate distance between current timestep polaron atom and next timestep
    # Then get all non-zero answer
    polaron_distances = np.zeros(num_timesteps)
    for j in range(num_timesteps - 1):
        polaron_distances[j] = distances.distance_array(universe.select_atoms('index {}'.format(polaron_atom_time[j])).positions,
                                                        universe.select_atoms('index {}'.format(polaron_atom_time[j+1])).positions,
                                                        box=box_size)
    polaron_distances = polaron_distances[xlim_1[0]:int(xlim_1[1])]
    # print(polaron_distances)

    # Remove polaron hops that occur within 50 fs of another hop
    min_residence = 50
    hops = polaron_distances > 0
    hop_indices = np.where(hops)[0]
    mask = np.ones_like(hop_indices, dtype=bool)
    for i in range(len(hop_indices)):
        if i < len(hop_indices) - 1:
            if hop_indices[i + 1] - hop_indices[i] <= min_residence:
                mask[i] = False
        else:
            if len(polaron_distances) - hop_indices[i] <= min_residence:
                mask[i] = False
    filtered_polaron_distances = polaron_distances.copy()
    filtered_polaron_distances[hop_indices[~mask]] = 0
    polaron_distances = filtered_polaron_distances

    polaron_distances_hop = polaron_distances[np.nonzero(polaron_distances)]
    polaron_indices = np.nonzero(polaron_distances)[0]
    # print('polaron_distance', polaron_distances)
    print('polaron_distances_hop', polaron_distances_hop)
    print('np.shape(polaron_distances_hop)[0]', np.shape(polaron_distances_hop)[0])
    print('polaron hop index', polaron_indices)

    print(np.shape(atoms_ti.positions))

    # print('polaron hop time', time_val_energy[polaron_indices])

    # Plot polaron distances
    metric = np.zeros((num_atoms_ti, num_timesteps))
    offset = 0
    # fig_bonds_2, ax_bonds_2 = plt.subplots()
    fig_bonds_2, ax_bonds_2 = plt.subplots(figsize=(10, 2))
    # ax_bonds_2.plot(time_val_1 - time_val_1[0] - offset, polaron_distances[:-1], 'kx-')
    ax_bonds_2.plot(time_val_1[:int(xlim_1[1])] - time_val_1[0] - offset, polaron_distances, 'kx-')
    ax_bonds_2.set_xlabel('Time / fs')
    # ax_bonds_2.set_xlabel('Timestep')
    ax_bonds_2.set_ylabel('Distance / A')
    # if draw_legend: ax_bonds_2.legend(frameon=True)
    # # ax_bonds_2.set_xlim([0, len(universe.trajectory)])
    ax_bonds_2.set_xlim(xlim_1)
    # ax_bonds_2.set_xlim([0, len(universe.trajectory) * timestep])
    # ax_bonds_2.set_ylim([0.06, -0.10])
    if save_fig: fig_bonds_2.savefig('{}/polaron_hopping_distance.png'.format(folder_save), dpi=300)
    fig_bonds_2.tight_layout()

    # Calculate mobility using xlim_1
    # hops_distance = np.array([2.99547136, 2.85545014, 3.01149688]) * 1e-8  # Angstrom to cm
    # hops_time = (5103+200)/2 * 1e-12  # fs to s
    hops_distance = polaron_distances_hop * 1e-8  # Angstrom to cm
    hops_time = (xlim_1[1] - xlim_1[0]) * 1e-15  # ps 1e-12 fs 1e-15
    print('hops per ps ', np.shape(hops_distance)[0]/hops_time*1e-15*1e3)

    rate_constant = np.shape(hops_distance)[0] / hops_time
    print('rate_constant', rate_constant)
    print('rate_constant / 1e12', rate_constant/1e12)
    if np.shape(hops_distance)[0] > 2: print('lifetime fs', 1/rate_constant * 1e15)

    mean_distance = np.mean(hops_distance)
    print('mean_distance', mean_distance)

    site_multiplicity = 1
    diffusion_constant_analytical = (np.mean(hops_distance) ** 2 * site_multiplicity * rate_constant) / 2
    print('diffusion_constant_analytical', diffusion_constant_analytical)
    mobility = (1.60217662e-19 * diffusion_constant_analytical) / (1.380649e-23 * temperature_set)
    print('mobility analytical', mobility, 'cm^2/(V·s)')

    # print('lifetime hematite', 1/1.2e12 * 1e15)
    # test = 1.60217662e-19 * (3e-8**2*3*1.2e12)/2 / (1.380649e-23 * 600)
    # print('mobility test hematite', test)

    # mobility = 1.60217662e-19 * (3e-8**2*3*rate_constant)/2 / (1.380649e-23 * 600)
    # print('mobility analytical 2', mobility, 'cm^2/(V·s)')

    mean_square_displacement = np.sum(hops_distance ** 2) / hops_time
    diffusion_constant_numerical = mean_square_displacement / 2
    mobility = (1.60217662e-19 * diffusion_constant_numerical) / (1.380649e-23 * temperature_set)
    print('diffusion_constant_numerical', diffusion_constant_numerical)
    print('mobility numerical', mobility, 'cm^2/(V·s)')

    k_el = 1
    temp = temperature_set
    kb_t_au = 8.617333262145E-5 * temperature_set  # KbT in eV
    kb_t = 1.38e-23 * temperature_set  # KbT in SI units
    # vn = 2.4e13  # 0.10 eV to s^-1 Deskins Dupuis TiO2 rutile (optic-mode phonon frequencies)
    vn = 2.66e13  # 0.11 eV to s^-1 Deskins Dupuis TiO2 anatase (optic-mode phonon frequencies)
    # rate_constant = vn * k_el * np.exp(-energy/kb_t)
    activation_energy = -np.log(rate_constant / (vn * k_el)) * kb_t_au
    print('activation_energy / meV', activation_energy*1e3)

    # hirshfeld and distance subplot
    rows, cols = 2, 1
    # fig_plot_all, ax_plot_all = plt.subplots(rows, cols,sharex='col', sharey='row',
    #                             figsize=(10, 6), gridspec_kw={'height_ratios': [2, 1],  'hspace': 0.05})
    fig_plot_all, ax_plot_all = plt.subplots(rows, cols,sharex='col', sharey='row',
                                figsize=(18, 6), gridspec_kw={'height_ratios': [2, 1],  'hspace': 0.05})
    temp = np.zeros(num_timesteps)
    for j in range(num_atoms):
        ax_plot_all[0].plot((time_val_1[:int(xlim_1[1])] - offset)/1e3, spin[:int(xlim_1[1]), j], '-', label='{}'.format(j + 1))
    if draw_legend: ax_plot_all[0].legend(frameon=True)
    # ax_plot_all[0].set_xlabel('Time / fs')
    ax_plot_all[0].set_ylabel('Spin moment')
    ax_plot_all[0].set_xlim((np.array(xlim_1)-offset)/1000)
    ax_plot_all[0].set_ylim(ylim_1)
    # ax_plot_all[0].set_ylim([0, 0.8])
    ax_plot_all[1].plot((time_val_1[:int(xlim_1[1])] - offset)/1e3, polaron_distances, 'kx-')
    ax_plot_all[1].set_xlabel('Time / ps')
    ax_plot_all[1].set_ylabel(r'Distance / $\mathrm{\AA}$')

    ax_plot_all[1].set_xlim((np.array(xlim_1)-offset)/1000)
    ax_plot_all[1].set_ylim([0, 3.3])
    fig_plot_all.tight_layout()
    fig_plot_all.subplots_adjust(hspace=0.05)
    if save_fig: fig_plot_all.savefig('{}/polaron_subplot.png'.format(folder_save), dpi=300)

# Plot total spin
fig_spin2, ax_spin2 = plt.subplots(figsize=(10, 4))
# fig_spin2, ax_spin2 = plt.subplots()
ax_spin2.plot(time_array[:int(xlim_1[1])], np.sum(spin[:int(xlim_1[1])], axis=1), 'k-')
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
    ax_spin1.plot(time_array[:int(xlim_1[1])], spin[:int(xlim_1[1]), i], '-')
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

# Plot MSD = cumulative polaron hopping distance**2
cumulative_sum = np.cumsum(polaron_distances**2)
# cumulative_sum_m, cumulative_sum_c = np.polyfit(time_val_1, cumulative_sum, 1)
# fitted_line = cumulative_sum_m * time_val_1 + cumulative_sum_c
# cumulative_sum_m, _, _, _ = np.linalg.lstsq(time_val_1, cumulative_sum, rcond=None)

# y = mx
# cumulative_sum_m, _ = curve_fit(linear_func, time_val_1[fit_start:int(xlim_1[1])], cumulative_sum[fit_start:])
# fitted_line = linear_func(time_val_1[fit_start:int(xlim_1[1])], cumulative_sum_m)

# y = mx + c
cumulative_sum_fit, _ = curve_fit(linear_func2, time_val_1[fit_start:int(xlim_1[1])], cumulative_sum[fit_start:])
cumulative_sum_m = cumulative_sum_fit[0]
cumulative_sum_c = cumulative_sum_fit[1]
fitted_line = linear_func2(time_val_1[fit_start:int(xlim_1[1])], cumulative_sum_m, cumulative_sum_c)

print('diffusion coefficient from gradient msd (units A**2 / fs)', cumulative_sum_m/2)
print('diffusion coefficient from gradient msd (units cm**2 / s)', 0.1*cumulative_sum_m/2)
diffusion_constant_msd = 0.1*cumulative_sum_m/2
mobility = (1.60217662e-19 * diffusion_constant_msd) / (1.380649e-23 * temperature_set)
rate_constant = 2 * diffusion_constant_msd / np.mean(hops_distance)**2
print('mobility from msd (units cm**2 / s)', mobility)
print('rate constant from msd (units / s)', rate_constant)
print('rate constant from msd (units e12 / s)', rate_constant/1e12)
activation_energy = -np.log(rate_constant / (vn * k_el)) * kb_t_au
print('activation_energy from msd (units meV)', activation_energy*1e3)

fig_msd, ax_msd = plt.subplots(figsize=(4, 4))
ax_msd.plot((time_val_1[:int(xlim_1[1])] - offset)/1e3, cumulative_sum, 'k-')
ax_msd.plot((time_val_1[fit_start:int(xlim_1[1])] - offset)/1e3, fitted_line, '--', color='grey')
ax_msd.set_xlim(0, time_array[-1])
ax_msd.set_xlabel("Time / ps")
ax_msd.set_ylabel(r"MSD / $\mathrm{\AA}^2$")
ax_msd.set_xlim(np.array(xlim_1)/1e3)
ax_msd.set_ylim([0, np.max(fitted_line)*1.02])
fig_msd.tight_layout()
fig_msd.savefig("{}/msd_cumulative.png".format(folder_save), dpi=600)
fig_msd.tight_layout()

# Plot RDF for Ti - O
nbins = 300
rdf_xlim = 13.77 / 2

# if calc_distance and plot_rdf:
#     ti_polaron = universe.select_atoms("bynum 78")
#     atoms_o = universe.select_atoms("element O")
#     rdf_ti_o_polaron = rdf.InterRDF(ti_polaron, atoms_o, range=(0, rdf_xlim), nbins=nbins)
#     rdf_ti_o_polaron.run()
#
#     ti_not_polaron = universe.select_atoms("element Ti and not bynum 78")
#     rdf_ti_o_not_polaron = rdf.InterRDF(ti_not_polaron, atoms_o, range=(0, rdf_xlim), nbins=nbins)
#     rdf_ti_o_not_polaron.run()
#
#     fig_rdf_ti_o, ax_rdf_ti_o = plt.subplots(figsize=(10, 4))
#     ax_rdf_ti_o.plot(rdf_ti_o_polaron.bins, rdf_ti_o_polaron.rdf, 'r-', label='RDF Ti (polaron) - O')
#     ax_rdf_ti_o.plot(rdf_ti_o_not_polaron.bins, rdf_ti_o_not_polaron.rdf, 'g-', label='RDF Ti - O')
#     ax_rdf_ti_o.set_xlabel("Radial distance / Å")
#     ax_rdf_ti_o.set_ylabel("RDF (arb. units)")
#     ax_rdf_ti_o.legend(frameon=False)
#     ax_rdf_ti_o.set_xlim([0, rdf_xlim])
#     ax_rdf_ti_o.tick_params(axis='y', which='both', left=False, right=False, labelleft=False)
#     fig_rdf_ti_o.tight_layout()
#     fig_rdf_ti_o.savefig("{}/rdf.png".format(folder), dpi=600)

# Plot RDF for Ti - Ti
# if calc_distance and plot_rdf:
#     ti_polaron = universe.select_atoms("bynum 78")
#     rdf_ti_ti_polaron = rdf.InterRDF(ti_polaron, ti_not_polaron, range=(0, rdf_xlim), nbins=nbins)
#     rdf_ti_ti_polaron.run()
#
#     ti_not_polaron = universe.select_atoms("element Ti and not bynum 78")
#     rdf_ti_ti_not_polaron = rdf.InterRDF(ti_not_polaron, ti_not_polaron, exclusion_block=(1, 1), range=(0, rdf_xlim), nbins=nbins)
#     rdf_ti_ti_not_polaron.run()
#
#     fig_rdf_ti_ti, ax_rdf_ti_ti = plt.subplots(figsize=(10, 4))
#     ax_rdf_ti_ti.plot(rdf_ti_ti_polaron.bins, rdf_ti_ti_polaron.rdf, 'r-', label='RDF Ti (polaron) - Ti')
#     ax_rdf_ti_ti.plot(rdf_ti_ti_not_polaron.bins, rdf_ti_ti_not_polaron.rdf, 'g-', label='RDF Ti - Ti')
#     ax_rdf_ti_ti.set_xlabel("Radial distance / Å")
#     ax_rdf_ti_ti.set_ylabel("RDF (arb. units)")
#     ax_rdf_ti_ti.legend(frameon=False)
#     ax_rdf_ti_ti.set_xlim([0, rdf_xlim])
#     ax_rdf_ti_ti.tick_params(axis='y', which='both', left=False, right=False, labelleft=False)
#     fig_rdf_ti_ti.tight_layout()
#     fig_rdf_ti_ti.savefig("{}/rdf.png".format(folder), dpi=600)


# if calc_distance and plot_rdf:
#     all_atoms = list(ase.io.read('{}/tio2-vel-1.xyz'.format(folder), index=':'))  # Each frame
#     velocities = np.array([atoms.positions for atoms in all_atoms])  # shape: (n_frames, n_atoms, 3)
#     vacf_all = compute_vacf(velocities)

    # dt = 1.0
    # N = len(vacf_all)
    # freqs = np.fft.fftfreq(N, dt*1e-15)[:N//2]
    # vdos_all = np.abs(np.fft.fft(vacf_all))[:N//2]
    # c = 2.99792458e10
    # freqs_cm1 = freqs / c
    # freqs_thz = freqs * 1e-12

    # fig_vdos1, ax_vdos1 = plt.subplots(figsize=(10, 4))
    # ax_vdos1.plot(freqs_cm1, vdos_all, 'k-')
    # ax_vdos1.set_xlim([0, 1000])
    # ax_vdos1.set_xlabel('Frequency (1/cm)')
    # ax_vdos1.set_ylabel('VDOS (arb. units)')
    # fig_vdos1.tight_layout()
    # fig_vdos1.savefig("{}/vdos_cm.png".format(folder), dpi=600)

    # fig_vdos2, ax_vdos2 = plt.subplots(figsize=(10, 4))
    # ax_vdos2.plot(freqs_thz, vdos_all, 'k-')
    # ax_vdos2.set_xlim([0, 32])
    # ax_vdos2.set_xlabel('Frequency (THz)')
    # ax_vdos2.set_ylabel('VDOS (arb. units)')
    # fig_vdos2.tight_layout()
    # fig_vdos2.savefig("{}/vdos_thz.png".format(folder), dpi=600)

if __name__ == "__main__":
    print('Finished.')
    plt.show()