import MDAnalysis as mda

import numpy as np


# folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/archer/bulk/221_supercell/md/hole/400K/input-file-EPS_SCHWARZ_FORCES-neutral2-rs-pbe'
# filename_input = "{}/hematite-pos-1.xyz".format(folder)
# filename_output = "{}/hematite-pos-1-wrapped.xyz".format(folder)
# cell = [10.071, 10.071, 13.747, 90, 90, 120]
# print('cell')

folder = '/Volumes/Samsung/Data/Postdoc2/Data/Work/temp/files/trajectory_mdanalysis'
filename_input = "{}/tio2-pos-1-orig.xyz".format(folder)
filename_output = "{}/tio2-pos-1.xyz".format(folder)
cell = [13.77, 13.77, 17.76, 90, 90, 90]

u = mda.Universe(filename_input, format="XYZ")

# Define the simulation cell dimensions (triclinic box)
u.dimensions = np.array(cell)  # [x_length, y_length, z_length, alpha, beta, gamma]
coordinates = u.atoms.positions

wrapped_coordinates = u.atoms.pack_into_box()

with open(filename_output, 'w') as f:
    # Loop through all frames in the trajectory

    for ts in u.trajectory:
        wrapped_coordinates = u.atoms.pack_into_box()

        f.write(f"{len(u.atoms)}\n")
        f.write(f"i = {ts.frame}, time = {(ts.frame)/2}, E = 0\n")

        for atom, coord in zip(u.atoms.types, wrapped_coordinates):
            f.write(f"{atom} {coord[0]:.6f} {coord[1]:.6f} {coord[2]:.6f}\n")

if __name__ == "__main__":
    print('Finished.')
    # plt.show()

