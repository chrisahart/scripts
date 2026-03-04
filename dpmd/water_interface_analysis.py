from ase.io import read
import numpy as np
import matplotlib.pyplot as plt

atoms_list = read("/Users/chris/Documents/Storage/other/tio2_water/cp2k-aimd-bulk/rutile/bulk/100/all/selected-20-cp2k-md-rutile-100_0.xyz", index=":")
box = [13.9553247, 8.90701634, 31.33880111]
z_bins = np.linspace(-5, 25, 300)
density = np.zeros(len(z_bins)-1)

for atoms in atoms_list:
    atoms.set_cell(box)
    atoms.set_pbc([True, True, True])
    positions = atoms.get_positions()
    symbols = atoms.get_chemical_symbols()
    z = positions[:,2]

    ti_mask = np.array(symbols) == "Ti"
    o_mask = np.array(symbols) == "O"
    z_ti_max = z[ti_mask].max()
    print(z_ti_max)
    print(14.29389-z_ti_max)

    z_rel = z[o_mask] - z_ti_max

    hist, _ = np.histogram(z_rel, bins=z_bins)
    density += hist

density /= len(atoms_list)

cell = atoms_list[0].get_cell()
area = cell[0,0] * cell[1,1]
dz = z_bins[1] - z_bins[0]

number_density = density / (area * dz)

plt.plot(0.5*(z_bins[1:]+z_bins[:-1]), number_density)
plt.show()
