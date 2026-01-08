import numpy as np
import matplotlib.pyplot as plt
from ase.io import read, write, Trajectory
import MDAnalysis as mda
from MDAnalysis.analysis import distances
from scipy.optimize import curve_fit
import os
from ase import Atoms
from general import parameters as param


def compute_vacf(vels):
    v = vels.reshape(vels.shape[0], -1)
    v -= v.mean(axis=0)
    nframes = v.shape[0]
    vacf = np.zeros(nframes, dtype=np.float32)
    for t in range(nframes):
        vacf[t] = np.mean(np.sum(v[:nframes-t] * v[t:], axis=1))
    return vacf

def linear_func(x, m):
    return m * x

def linear_func2(x, m, c):
    return m * x + c

# --- Parameters ---
topology_file = '/Volumes/Samsung/Data/Postdoc2/Data/Work/calculations/tio2/rutile/deepmd/rutile/336/md-cell-opt/system.xyz'
folder_1 = '/Volumes/Samsung/Data/Postdoc2/Data/Work/calculations/tio2/rutile/deepmd/rutile/336/md-cell-opt/deepmd-md/hse-22-ts-md2'
folder_energy = 'single-fit-ener-dpa3-nlayers-6-official-v3.1.0-start_pref-0.02-1000_limit_pref-1-1-twostep-lr-1e-5-1e-8'
folder_spin = 'single-fit-pop-dpa3-nlayers-6-official-v3.1.0-dev-polaron-loss-mae-pref-1-pref_pop-1000-1'

# folder_md = '50k-1000-ps-vel-9.3-ps-friction-0.001'  # 1 hops in 1000 ps, no polaron mobility
# fit_start = int(400*1e3)
# temperature_set = 50

# folder_md = '100k-10-ps-vel-9.3-ps-friction-0.001'  # 4 hops in 10 ps
# folder_md = '100k-50-ps-vel-9.3-ps-friction-0.001'  # 9 hops in 50 ps
# fit_start = 5000
# temperature_set = 100

# folder_md = '100k-1000-ps-vel-9.3-ps-friction-0.001'  # 109 hops in 1000 ps
# fit_start = int(400*1e3)
# temperature_set = 100

folder = '{}/{}/{}/{}'.format(folder_1, folder_energy, folder_spin, folder_md)

num_atoms = 324
box_size = [13.77, 13.77, 17.76, 90, 90, 90]
local_bonds = 6
timestep = 1
save_fig = True
calc_distance = True
plot_rdf = False

# --- Load Data ---
spin = np.load(f"{folder}/spin_history.npy", mmap_mode='r')  # Memory-map large arrays
charge_state = np.load(f"{folder}/charge_state_history.npy", mmap_mode='r')
num_timesteps = spin.shape[0]
time_array = np.linspace(0, num_timesteps, num=num_timesteps, dtype=np.float32)

# --- Convert Trajectory to XYZ (Chunked) ---
pos_file = f'{folder}/tio2-pos-1.xyz'
vel_file = f'{folder}/tio2-vel-1.xyz'

# Delete old files
for file in [pos_file, vel_file]:
    if os.path.exists(file):
        os.remove(file)

# Write positions and velocities in chunks
with Trajectory(f'{folder}/md.traj') as traj:
    for i, atoms in enumerate(traj):
        if i == 0:
            continue
        write(pos_file, atoms, format='xyz', append=(i > 1))
        v_atoms = Atoms(symbols=atoms.get_chemical_symbols(), positions=atoms.get_velocities())
        write(vel_file, v_atoms, format='xyz', append=(i > 1))
        if i % 100 == 0:  # Free memory periodically
            del atoms, v_atoms

# --- MDAnalysis Setup ---
universe = mda.Universe(topology_file, pos_file)
universe.dimensions = box_size

# Select atoms
atoms_ti = universe.select_atoms('name Ti')
atoms_o = universe.select_atoms('name O')
num_atoms_ti = len(atoms_ti)
num_atoms_o = len(atoms_o)

# --- Bond Length Analysis (Vectorized) ---
bond_lengths_time = np.zeros((num_timesteps, num_atoms_ti, local_bonds), dtype=np.float32)
for ts in universe.trajectory:
    frame = ts.frame
    dist_matrix = distances.distance_array(atoms_ti.positions, atoms_o.positions, box=box_size)
    bond_lengths_time[frame] = np.partition(dist_matrix, local_bonds, axis=1)[:, :local_bonds]
    del dist_matrix  # Free memory

bond_lengths_time_sorted_mean = np.mean(bond_lengths_time, axis=2)

# --- Polaron Analysis ---
polaron_atom_time = np.argmax(spin, axis=1)
polaron_atoms = np.unique(polaron_atom_time)

# Filter rapid hops
min_residence = 50
hop_mask = np.ones_like(polaron_atom_time, dtype=bool)
hop_mask[:-1] = (polaron_atom_time[:-1] != polaron_atom_time[1:])
hop_mask[-1] = True
polaron_hops = polaron_atom_time[hop_mask]

# Calculate hop distances
polaron_distances = np.zeros(num_timesteps, dtype=np.float32)
for i in range(num_timesteps - 1):
    if hop_mask[i]:
        polaron_distances[i] = distances.distance_array(
            universe.select_atoms(f'index {polaron_atom_time[i]}').positions,
            universe.select_atoms(f'index {polaron_atom_time[i+1]}').positions,
            box=box_size
        )

# Filter hops within `min_residence` timesteps
hop_indices = np.where(polaron_distances > 0)[0]
mask = np.ones_like(hop_indices, dtype=bool)
for j in range(len(hop_indices) - 1):
    if hop_indices[j+1] - hop_indices[j] <= min_residence:
        mask[j] = False
polaron_distances[hop_indices[~mask]] = 0
polaron_distances_hop = polaron_distances[np.nonzero(polaron_distances)]
polaron_indices = np.nonzero(polaron_distances)[0]
print('polaron_distances_hop', polaron_distances_hop)
print('np.shape(polaron_distances_hop)[0]', np.shape(polaron_distances_hop)[0])
print('polaron hop index', polaron_indices)

# --- Mobility and Diffusion ---
hops_distance = polaron_distances[polaron_distances > 0] * 1e-8  # Angstrom to cm
hops_time = num_timesteps * timestep * 1e-15  # Total time in seconds
rate_constant = len(hops_distance) / hops_time
diffusion_constant = (np.mean(hops_distance)**2 * rate_constant) / 2
mobility = (1.60217662e-19 * diffusion_constant) / (1.380649e-23 * temperature_set)

# --- MSD Analysis ---
draw_legend = False
k_el = 1
temp = temperature_set
kb_t_au = 8.617333262145E-5 * temperature_set  # KbT in eV
kb_t = 1.38e-23 * temperature_set  # KbT in SI units
vn = 2.4e13  # 0.10 eV to s^-1 Deskins Dupuis TiO2 rutile (optic-mode phonon frequencies)

xlim_1 = [0, time_array[-1]]
offset = 0
cumulative_sum = np.cumsum(polaron_distances**2)

cumulative_sum_fit, _ = curve_fit(linear_func2, time_array[fit_start:int(xlim_1[1])], cumulative_sum[fit_start:])
cumulative_sum_m = cumulative_sum_fit[0]
cumulative_sum_c = cumulative_sum_fit[1]
fitted_line = linear_func2(time_array[fit_start:int(xlim_1[1])], cumulative_sum_m, cumulative_sum_c)

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
ax_msd.plot((time_array[:int(xlim_1[1])] - offset)/1e3, cumulative_sum, 'k-')
ax_msd.plot((time_array[fit_start:int(xlim_1[1])] - offset)/1e3, fitted_line, '--', color='grey')
ax_msd.set_xlim(0, time_array[-1])
ax_msd.set_xlabel("Time / ps")
ax_msd.set_ylabel(r"MSD / $\mathrm{\AA}^2$")
ax_msd.set_xlim(np.array(xlim_1)/1e3)
ax_msd.set_ylim([0, np.max(fitted_line)*1.02])
fig_msd.tight_layout()
fig_msd.savefig("{}/msd_cumulative.png".format(folder), dpi=600)
fig_msd.tight_layout()

# hirshfeld and distance subplot
ylim_1 = [-0.02, 0.8]
rows, cols = 2, 1
fig_plot_all, ax_plot_all = plt.subplots(rows, cols, sharex='col', sharey='row',
                                         figsize=(18, 6), gridspec_kw={'height_ratios': [2, 1], 'hspace': 0.05})
temp = np.zeros(num_timesteps)
for j in range(num_atoms):
    ax_plot_all[0].plot((time_array[:int(xlim_1[1])] - offset) / 1e3, spin[:int(xlim_1[1]), j], '-',
                        label='{}'.format(j + 1))
if draw_legend: ax_plot_all[0].legend(frameon=True)
# ax_plot_all[0].set_xlabel('Time / fs')
ax_plot_all[0].set_ylabel('Spin moment')
ax_plot_all[0].set_xlim((np.array(xlim_1) - offset) / 1000)
ax_plot_all[0].set_ylim(ylim_1)
# ax_plot_all[0].set_ylim([0, 0.8])
ax_plot_all[1].plot((time_array[:int(xlim_1[1])] - offset) / 1e3, polaron_distances, 'kx-')
ax_plot_all[1].set_xlabel('Time / ps')
ax_plot_all[1].set_ylabel(r'Distance / $\mathrm{\AA}$')
ax_plot_all[1].set_xlim((np.array(xlim_1) - offset) / 1000)
ax_plot_all[1].set_ylim([0, 3.3])
fig_plot_all.tight_layout()
fig_plot_all.subplots_adjust(hspace=0.05)
if save_fig: fig_plot_all.savefig('{}/polaron_subplot.png'.format(folder), dpi=300)


# --- Cleanup ---
del bond_lengths_time, bond_lengths_time_sorted_mean, polaron_distances, spin

if __name__ == "__main__":
    print('Finished.')
    plt.show()

