import numpy as np
import matplotlib.pyplot as plt
from ase.io import read

# ── Settings ──────────────────────────────────────────────────────────────────
FILE = "/Users/chris/Documents/Storage/other/tio2_water/cp2k-aimd-bulk/rutile/bulk/100/selected-20-cp2k-md-rutile-100.xyz"

BIN_SIZE = 0.1               # Å
Z_MAX    = 15.0              # Å above surface to plot

# ── Load structure ────────────────────────────────────────────────────────────
atoms   = read(FILE)
symbols = np.array(atoms.get_chemical_symbols())
pos     = atoms.get_positions()

# ── Outermost Ti layer z-reference ───────────────────────────────────────────
ti_z  = pos[symbols == "Ti", 2]
z_ref = ti_z[ti_z >= ti_z.max() - 0.5].mean()
print(f"Outermost Ti layer z = {z_ref:.3f} Å")

# ── Find water O atoms (have H neighbour < 1.3 Å, no Ti neighbour < 2.5 Å) ──
water_o_z = []
for i in np.where(symbols == "O")[0]:
    d_h  = np.linalg.norm(pos[symbols == "H"]  - pos[i], axis=1)
    d_ti = np.linalg.norm(pos[symbols == "Ti"] - pos[i], axis=1)
    if d_h.min() < 1.3 and d_ti.min() >= 2.5:
        water_o_z.append(pos[i, 2])

water_o_z = np.array(water_o_z) - z_ref
water_o_z = water_o_z[water_o_z >= 0]
print(f"Water molecules found: {len(water_o_z)}")

# ── Density profile ───────────────────────────────────────────────────────────
Lx, Ly    = atoms.get_cell()[0, 0], atoms.get_cell()[1, 1]
bins      = np.arange(0, Z_MAX + BIN_SIZE, BIN_SIZE)
counts, edges = np.histogram(water_o_z, bins=bins)
z_bins    = 0.5 * (edges[:-1] + edges[1:])
density   = counts / (Lx * Ly * BIN_SIZE)   # Å⁻³

# ── Plot ──────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(8, 4))
ax.plot(z_bins, density, color="steelblue", linewidth=1.5)
ax.fill_between(z_bins, density, alpha=0.25, color="steelblue")
ax.axhline(0.0334, color="firebrick", linestyle="--", linewidth=1, label="Bulk water (0.0334 Å⁻³)")
ax.set_xlabel("Distance from outermost Ti layer (Å)", fontsize=12)
ax.set_ylabel("Water O number density (Å⁻³)", fontsize=12)
ax.set_title("Water density profile above TiO₂ surface", fontsize=13)
ax.set_xlim(0, Z_MAX)
ax.set_ylim(bottom=0)
ax.legend()
plt.tight_layout()
plt.savefig("water_density_profile.png", dpi=150)
plt.show()