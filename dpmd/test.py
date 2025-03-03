import numpy as np
import ase
from ase.io import read
import matplotlib.pyplot as plt

folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/chris_phd/bulk/hole/hse/hops-0'
trajectory_position = read(f'{folder}/hematite-pos-1.xyz', index=':')
trajectory_force = read(f'{folder}/hematite-frc-1.xyz', index=':')

fe_indices = [i for i, atom in enumerate(trajectory_position[0]) if atom.symbol == 'Fe']
o_indices = [i for i, atom in enumerate(trajectory_position[0]) if atom.symbol == 'O']

# Define a cutoff distance for nearest neighbors

cutoff = 3

# Function to find nearest neighbors within a cutoff distance

def find_nearest_neighbors(fe_indices, o_indices, positions, cutoff):
    neighbors = {}
    for fe_idx in fe_indices:
        fe_pos = positions[fe_idx]
        nearest = []
        for o_idx in o_indices:
            o_pos = positions[o_idx]
            dist = np.linalg.norm(fe_pos - o_pos)
            if dist < cutoff:
                nearest.append((o_idx, dist, fe_pos - o_pos))
        nearest.sort(key=lambda x: x[1])
        neighbors[fe_idx] = nearest[:6]
    return neighbors

# Function to project forces along the direction of the bond
def project_forces_along_bond(bond_direction, force):
    bond_direction_normalized = bond_direction / np.linalg.norm(bond_direction)
    projection = np.dot(force, bond_direction_normalized) * bond_direction_normalized

    return projection

# Function to calculate bond lengths and forces for a single frame

def calculate_bond_lengths_and_forces(frame, forces, fe_indices, o_indices, cutoff):
    bond_lengths = []
    bond_forces = []
    bond_forces_magnitude = []
    neighbors = find_nearest_neighbors(fe_indices, o_indices, frame.positions, cutoff)
    for fe_idx, nearest in neighbors.items():
        bond_lengths_for_fe = []
        bond_forces_magnitude_for_fe = []
        for o_idx, dist, bond_direction in nearest:
            bond_lengths_for_fe.append(dist)
            force = forces[o_idx]
            projected_force = project_forces_along_bond(bond_direction, force)
            bond_forces_magnitude_for_fe.append(np.linalg.norm(projected_force))
        bond_lengths.append(np.mean(bond_lengths_for_fe))
        bond_forces_magnitude.append(np.mean(bond_forces_magnitude_for_fe))
    return bond_lengths, bond_forces_magnitude

# Initialize lists to store bond lengths and time steps

time_steps = []
all_bond_lengths = []
all_bond_forces_magnitude = []

# Loop through each frame in the trajectory

for i, frame in enumerate(trajectory_position):
    forces = trajectory_force[i].positions
    # forces = trajectory_force[i].get_forces()
    bond_lengths, bond_forces_magnitude = calculate_bond_lengths_and_forces(frame, forces, fe_indices, o_indices, cutoff)
    time_steps.extend([i] * len(bond_lengths))  # Time step for each bond length

    all_bond_lengths.extend(bond_lengths)
    all_bond_forces_magnitude.extend(bond_forces_magnitude)

# Plot bond lengths against time

plt.figure(figsize=(12, 6))
plt.plot(time_steps, all_bond_lengths, 'k-', markersize=2, alpha=0.5)
ylim_2 = [2.16, 1.89]
plt.ylim(ylim_2)
plt.title('Fe-O Bond Lengths Over Time')
plt.xlabel('Time Step')
plt.ylabel('Bond Length (Å)')
plt.show()

# Plot bond forces magnitude against time

plt.figure(figsize=(12, 6))
plt.plot(time_steps, all_bond_forces_magnitude, 'o', markersize=2, alpha=0.5)
plt.title('Fe-O Bond Forces Magnitude Over Time')
plt.xlabel('Time Step')
plt.ylabel('Bond Force Magnitude (eV/Å)')
plt.show()