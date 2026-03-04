import numpy as np
import matplotlib.pyplot as plt
from ase.io import read
from ase.neighborlist import neighbor_list

# ================= USER INPUT =================
filename = "/Users/chris/Documents/Storage/other/tio2_water/cp2k-aimd-bulk/rutile/bulk/100/selected-20-cp2k-md-rutile-100.xyz"
box = [13.9553247, 8.90701634, 31.33880111]  # Å
z_min = -20
z_max = 20
nbins = 400
top_ti_average = 4
OH_cutoff = 1.25  # Å (water O–H bond cutoff)
# ==============================================

atoms_list = read(filename, index=":")

Lx, Ly, Lz = box
area = Lx * Ly

z_bins = np.linspace(z_min, z_max, nbins)
density_counts = np.zeros(len(z_bins) - 1)

# =====================================================
#  Loop over trajectory
# =====================================================

for atoms in atoms_list:

    atoms.set_cell(box)
    atoms.set_pbc([True, True, True])

    positions = atoms.get_positions()
    symbols = np.array(atoms.get_chemical_symbols())
    z = positions[:, 2]

    # ---- Identify slab surfaces ----
    ti_mask = (symbols == "Ti")
    z_ti_sorted = np.sort(z[ti_mask])
    z_bottom = np.mean(z_ti_sorted[:top_ti_average])
    z_top    = np.mean(z_ti_sorted[-top_ti_average:])
    slab_mid = 0.5 * (z_bottom + z_top)

    bulk_mid = (slab_mid + Lz/2.0) % Lz

    # =====================================================
    #  Identify water oxygen via bonding
    # =====================================================

    # Get neighbor list for O-H pairs
    cutoffs = [OH_cutoff if s in ["O","H"] else 0.0 for s in symbols]
    i_indices, j_indices, distances = neighbor_list(
        "ijd", atoms, cutoffs
    )

    water_O_indices = []

    for i in range(len(atoms)):
        if symbols[i] != "O":
            continue

        # find bonded H within cutoff
        bonded_H = [
            j for ii, j, d in zip(i_indices, j_indices, distances)
            if ii == i and symbols[j] == "H"
        ]

        if len(bonded_H) == 2:
            water_O_indices.append(i)

    water_O_indices = np.array(water_O_indices)

    z_water = z[water_O_indices]

    # ---- relative to bulk midpoint ----
    z_rel = z_water - bulk_mid
    z_rel = (z_rel + Lz/2) % Lz - Lz/2

    hist, _ = np.histogram(z_rel, bins=z_bins)
    density_counts += hist

# =====================================================
#  Normalize density
# =====================================================

nframes = len(atoms_list)
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
plt.title("Water Density Profile (Water Only)")
plt.tight_layout()
plt.show()