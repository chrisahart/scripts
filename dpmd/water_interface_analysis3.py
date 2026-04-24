import numpy as np
import matplotlib.pyplot as plt
from ase.io import read

# ================= USER INPUT =================
filename = "/Users/chris/Documents/Storage/calculations/tio2-h2o/zeng/rutile/cp2k-aimd-bulk/100/md/zeng_input/optb88/tio2-pos-1.xyz"
box = [13.9553247, 8.90701634, 31.33880111]  # Å
# filename = "/Volumes/Elements/Data/Postdoc2/Data/Work/calculations/tio2-h2o/ahart/110/slab/geo_opt/pbe-frozen-tio2/tio2_water.xyz"
# box = [12.9824805, 17.76, 33.3469092]  # Å
z_min = -10
z_max = 10
nbins = 400
top_ti_average = 9
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

    ti_mask = (symbols == "Ti")

    # --- determine bottom and top Ti surfaces ---
    z_ti_sorted = np.sort(z[ti_mask])

    z_bottom = np.mean(z_ti_sorted[:top_ti_average])
    z_top    = np.mean(z_ti_sorted[-top_ti_average:])

    # --- slab midpoint ---
    slab_mid = 0.5 * (z_bottom + z_top)

    # --- due to PBC, bulk midpoint is opposite slab midpoint ---
    # shift by half box length
    bulk_mid = slab_mid + Lz / 2.0

    # wrap into box
    bulk_mid = bulk_mid % Lz

    # --- water oxygen ---
    o_mask = (symbols == "O")
    water_o_mask = o_mask & (~ti_mask)

    z_rel = z[water_o_mask] - bulk_mid

    # apply minimum image convention along z
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
plt.title("Water Density Profile (PBC Midpoint Reference)")
plt.tight_layout()
plt.show()