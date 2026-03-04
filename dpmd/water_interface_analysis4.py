import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt

# Load the XYZ file
u = mda.Universe('/Users/chris/Documents/Storage/other/tio2_water/cp2k-aimd-bulk/rutile/bulk/100/selected-20-cp2k-md-rutile-100.xyz')

# Select oxygen and hydrogen atoms
oxygen = u.select_atoms('name O')
hydrogen = u.select_atoms('name H')

# Select Ti atoms
ti = u.select_atoms('name Ti')

# Get the maximum z-coordinate of Ti atoms (outmost Ti layer)
max_ti_z = max(ti.positions[:, 2])

# Group oxygen and hydrogen atoms into water molecules
# This is a simplified approach; for a real system, use a proper bonding analysis
water_coms = []
for o_atom in oxygen:
    o_pos = o_atom.position
    # Find the two closest hydrogens to this oxygen
    distances = np.linalg.norm(hydrogen.positions - o_pos, axis=1)
    closest_h_indices = np.argsort(distances)[:2]
    closest_h_positions = hydrogen.positions[closest_h_indices]
    # Calculate the center of mass (COM) of the water molecule
    com = (o_pos + np.sum(closest_h_positions, axis=0)) / 3
    water_coms.append(com)

water_coms = np.array(water_coms)

# Calculate the water density profile
bin_width = 0.2  # in Angstrom
bins = np.arange(0, 30, bin_width)
water_z = water_coms[:, 2]
water_density, _ = np.histogram(water_z, bins=bins, weights=np.ones_like(water_z))

# Normalize by bin width to get density (molecules/Å)
water_density = water_density / bin_width

# Calculate the height relative to the outmost Ti layer
bin_centers = (bins[:-1] + bins[1:]) / 2
relative_height = bin_centers - max_ti_z

# Plot the water density profile
plt.figure(figsize=(10, 6))
plt.plot(relative_height, water_density, label='Water Density Profile')
plt.axvline(x=0, color='red', linestyle='--', label='Outmost Ti Layer')
plt.xlabel('Height from Outmost Ti Layer (Å)')
plt.ylabel('Water Density (molecules/Å)')
plt.title('Water Density Profile at TiO₂/Water Interface')
plt.legend()
plt.grid(True)
plt.savefig('water_density_profile.png', format='png', dpi=200, bbox_inches='tight')
plt.show()