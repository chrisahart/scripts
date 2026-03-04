import numpy as np
import matplotlib.pyplot as plt
import MDAnalysis as mda
from MDAnalysis.lib.distances import distance_array

# ================= USER INPUT =================
filename = "/Users/chris/Documents/Storage/other/tio2_water/cp2k-aimd-bulk/rutile/bulk/100/selected-20-cp2k-md-rutile-100.xyz"
box = np.array([13.9553247, 8.90701634, 31.33880111])  # Å
z_min = -20
z_max = 20
nbins = 400
top_ti_average = 4
OH_cutoff = 1.25  # Å
# ==============================================

# ---- load trajectory ----
u = mda.Universe(filename, format="XYZ")

# Set box dimensions manually (XYZ does not contain box)
for ts in u.trajectory:
    ts.dimensions = np.array([box[0], box[1], box[2], 90, 90, 90])

Lx, Ly, Lz = box
area = Lx * Ly

z_bins = np.linspace(z_min, z_max, nbins)
density_counts = np.zeros(len(z_bins) - 1)

# =====================================================
#  Loop over trajectory
# =====================================================

for ts in u.trajectory:

    positions = u.atoms.positions
    symbols = u.atoms.names  # for XYZ this stores element names
    z = positions[:, 2]

    # ---- Identify Ti surfaces ----
    ti = u.select_atoms("name Ti")
    z_ti_sorted = np.sort(ti.positions[:, 2])

    z_bottom = np.mean(z_ti_sorted[:top_ti_average])
    z_top    = np.mean(z_ti_sorted[-top_ti_average:])

    slab_mid = 0.5 * (z_bottom + z_top)

    # ---- Bulk midpoint via PBC ----
    bulk_mid = (slab_mid + Lz / 2.0) % Lz

    # =====================================================
    #  Identify water oxygen via O–H bonding
    # =====================================================

    O_atoms = u.select_atoms("name O")
    H_atoms = u.select_atoms("name H")

    # distance matrix with PBC
    dists = distance_array(O_atoms.positions,
                           H_atoms.positions,
                           box=ts.dimensions)

    water_O_indices = []

    for i, dist_row in enumerate(dists):
        bonded_H = np.where(dist_row < OH_cutoff)[0]
        if len(bonded_H) == 2:
            water_O_indices.append(O_atoms.indices[i])

    water_O_indices = np.array(water_O_indices)

    z_water = positions[water_O_indices, 2]

    # ---- relative coordinate (minimum image) ----
    z_rel = z_water - bulk_mid
    z_rel = (z_rel + Lz/2) % Lz - Lz/2

    hist, _ = np.histogram(z_rel, bins=z_bins)
    density_counts += hist

# =====================================================
#  Normalize density
# =====================================================

nframes = len(u.trajectory)
dz = z_bins[1] - z_bins[0]

number_density = density_counts / (nframes * area * dz)

# convert to g/cm³
rho_g_cm3 = number_density * 18.01528 / 6.022e23 * 1e24

z_centers = 0.5 * (z_bins[:-1] + z_bins[1:])

# =====================================================
#  Plot
# =====================================================

plt.figure(figsize=(6,4))
plt.plot(z_centers, rho_g_cm3)
plt.axvline(0, linestyle="--")
plt.xlabel("z (Å)  [0 = midpoint between Ti surfaces]")
plt.ylabel("Water density (g/cm³)")
plt.title("Water Density Profile (MDAnalysis)")
plt.tight_layout()
plt.show()