import numpy as np
import matplotlib.pyplot as plt
from ase.io import read

# -------- USER INPUT --------
# filename = read("/Users/chris/Documents/Storage/other/tio2_water/cp2k-aimd-bulk/rutile/bulk/100/selected-20-cp2k-md-rutile-100.xyz", index=":")
filename = "/Users/chris/Documents/Storage/other/tio2_water/cp2k-aimd-bulk/rutile/bulk/100/selected-20-cp2k-md-rutile-100.xyz"
box = [13.9553247, 8.90701634, 31.33880111]  # Å
z_min = 0
z_max = 20
nbins = 300
# ----------------------------

atoms_list = read(filename, index=":")

Lx, Ly, Lz = box
area = Lx * Ly
z_bins = np.linspace(z_min, z_max, nbins)
density_counts = np.zeros(len(z_bins) - 1)

for atoms in atoms_list:

    atoms.set_cell(box)
    atoms.set_pbc([True, True, True])

    positions = atoms.get_positions()
    symbols = np.array(atoms.get_chemical_symbols())

    z = positions[:, 2]

    # ---- surface reference: outermost Ti ----
    ti_mask = (symbols == "Ti")
    z_surface = z[ti_mask].max()

    # ---- water oxygen identification ----
    # assume water O above slab
    o_mask = (symbols == "O")
    water_o_mask = o_mask & (z > z_surface - 1.0)

    z_rel = z[water_o_mask] - z_surface

    hist, _ = np.histogram(z_rel, bins=z_bins)
    density_counts += hist

# ---- normalization ----
nframes = len(atoms_list)
dz = z_bins[1] - z_bins[0]

number_density = density_counts / (nframes * area * dz)

# Optional: convert to g/cm³
rho_g_cm3 = number_density * 18.01528 / 6.022e23 * 1e24

z_centers = 0.5 * (z_bins[:-1] + z_bins[1:])

# ---- plot ----
plt.figure()
plt.plot(z_centers, rho_g_cm3)
plt.xlabel("z - z_Ti_surface (Å)")
plt.ylabel("Water density (g/cm³)")
plt.title("Water Density Profile at TiO2 Interface")
plt.show()