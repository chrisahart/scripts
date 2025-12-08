import MDAnalysis as mda

import numpy as np


# folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/liu_group/archer/bulk/221_supercell/md/hole/400K/input-file-EPS_SCHWARZ_FORCES-neutral2-rs-pbe'
# filename_input = "{}/hematite-pos-1.xyz".format(folder)
# filename_output = "{}/hematite-pos-1-wrapped.xyz".format(folder)
# cell = [10.071, 10.071, 13.747, 90, 90, 120]
# print('cell')

# folder = '/Volumes/Samsung/Data/Postdoc2/Data/Work/temp/files/trajectory_mdanalysis_wrap'
# filename_input = "{}/tio2-pos-1-orig.xyz".format(folder)
# filename_output = "{}/tio2-pos-1.xyz".format(folder)
# cell = [13.77, 13.77, 17.76, 90, 90, 90]
# desired_pos = np.array([0.00000, 0.00000, 1.47961])

folder = '/Volumes/Samsung/Data/Postdoc2/Data/Work/calculations/tio2/anatase/archer/anatase/cell-441/md-cell-opt-hse-20/hse-19-complete/trajectory_mdanalysis_cell_opt'
filename_input = "{}/tio2-pos-1-orig.xyz".format(folder)
filename_output = "{}/tio2-pos-1.xyz".format(folder)
# cell = [15.12, 15.12, 9.62, 90, 90, 90]
cell = [15.08, 15.08, 9.68, 90, 90, 90]
desired_pos = np.array([0.00000, 0.00000, 0.00000])

u = mda.Universe(filename_input, format="XYZ")
u.dimensions = np.array(cell)

# Get the position of atom 1 in the first frame
atom1_pos = u.atoms.positions[0]

# Get the position of atom 1 in the first frame after centering the COM
first_frame = u.trajectory[0]
com = u.atoms.center_of_mass()
centered_coords = u.atoms.positions - com
atom1_pos = centered_coords[0]

# Calculate the translation vector to move atom 1 to the desired position
translation = desired_pos - atom1_pos

with open(filename_output, 'w') as f:
    for ts in u.trajectory:
        # Center the COM at the origin
        com = u.atoms.center_of_mass()
        centered_coords = u.atoms.positions - com

        # Translate all atoms to move atom 1 to the desired position
        translated_coords = centered_coords + translation
        u.atoms.positions = translated_coords

        # Wrap all atoms into the simulation box
        # wrapped_coords = u.atoms.wrap(compound='atoms')
        wrapped_coords = translated_coords

        f.write(f"{len(u.atoms)}\n")
        f.write(f"i = {ts.frame}, time = {ts.frame/2}, E = 0\n")
        for atom, coord in zip(u.atoms.names, wrapped_coords):
            f.write(f"{atom} {coord[0]:.6f} {coord[1]:.6f} {coord[2]:.6f}\n")

print('Finished.')
if __name__ == "__main__":
    print('Finished.')
    # plt.show()

