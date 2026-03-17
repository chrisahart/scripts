from ase.io import read
import numpy as np
import matplotlib.pyplot as plt
from general import parameters as param

folder = "/Users/chris/Documents/Storage/calculations/tio2-h2o/zeng/rutile/cp2k-aimd-bulk/100/md/zeng_input/optb88"
filename = "/Users/chris/Documents/Storage/calculations/tio2-h2o/zeng/rutile/cp2k-aimd-bulk/100/md/zeng_input/optb88/tio2-pos-1.xyz"
box = [13.9553247, 8.90701634, 31.33880111]  # Å
z_min = 0
z_max = 20
nbins = 300
box = [13.9553247, 8.90701634, 31.33880111]
Lx, Ly, Lz = box
area = Lx * Ly
z_bins = np.linspace(z_min, z_max, nbins)
density_counts = np.zeros(len(z_bins)-1)
atoms_list = read(filename, index=":")

for atoms in atoms_list:

    # print(atoms)

    positions = atoms.get_positions()

    atoms.set_cell(box)
    atoms.set_pbc([True, True, True])

    positions = atoms.get_positions()
    symbols = np.array(atoms.get_chemical_symbols())

    z = positions[:, 2]
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

fig_trapping1, ax_trapping1 = plt.subplots(figsize=(5, 3))
ax_trapping1.plot(0.5*(z_bins[1:]+z_bins[:-1]), number_density)
ax_trapping1.set_ylabel('ρ(O)')
ax_trapping1.set_xlabel('h / Å')
ax_trapping1.set_xlim([0, 8])
ax_trapping1.set_ylim([0, 0.32])
fig_trapping1.tight_layout()
fig_trapping1.savefig('{}/density_o.png'.format(folder), dpi=600)

plt.show()
