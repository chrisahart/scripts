import MDAnalysis as mda

import numpy as np


# folder = '/Users/chris/Documents/Storage/calculations/h2o/packmol/cubic/64/md/tip3p/testing/xyz2'
# folder = '/Users/chris/Documents/Storage/calculations/h2o/packmol/cubic/64/md/tip3p/equil/temp-300'
# filename_input = "{}/last.xyz".format(folder)
# filename_output = "{}/last_center.xyz".format(folder)
# cell = [12.4, 12.4, 12.4, 90, 90, 90]

# folder = '/Volumes/Elements/Data/Postdoc2/Data/Work/calculations/h2o/packmol/cubic/128/md/equil/temp-300'
# cell = [15.65, 15.65, 15.65, 90, 90, 90]
# folder = '/Volumes/Elements/Data/Postdoc2/Data/Work/calculations/h2o/packmol/tio2_water/110/md/equil/temp-300'
# cell = [12.9824805,  17.76,  16.6100492, 90, 90, 90]
# filename_input = "{}/last.xyz".format(folder)
# filename_output = "{}/last_center.xyz".format(folder)

folder = '/Volumes/Elements/Data/Postdoc2/Data/Work/calculations/tio2-h2o/ahart/110/water/md/pbe-frozen-tio2'
cell = [12.9824805, 17.76, 33.3469092, 90.0, 90.0, 90.0]
filename_input = "{}/last.xyz".format(folder)
filename_output = "{}/last_center.xyz".format(folder)

folder = '/Users/chris/Documents/Storage/calculations/tio2-h2o/zeng/rutile/cp2k-aimd-bulk/110/md/pbe-d-d3-300k-temptol-30-nose-100'
cell = [13.155592560000001, 8.9070163400000002, 31.389484599999999, 90.0, 90.0, 90.0]
filename_input = "{}/last.xyz".format(folder)
filename_output = "{}/last_center.xyz".format(folder)
filename_input = "{}/tio2-pos-1.xyz".format(folder)
filename_output = "{}/tio2-pos-1_center.xyz".format(folder)

folder = '/Volumes/Elements/Data/Postdoc2/Data/Work/calculations/h2o/packmol/tio2_water/110/md/equil/pbe/from-tip3p-1-ns'
input_file = f"{folder}/last.xyz"
output_file = f"{folder}/last_center.xyz"
cell = [12.9824805, 17.76, 16.6100492, 90.0, 90.0, 90.0]

u = mda.Universe(filename_input, format="XYZ")
u.dimensions = np.array(cell)

with open(filename_output, 'w') as f:
    for ts in u.trajectory:

        # Wrap atoms into box
        u.atoms.wrap(compound='atoms')

        # Center system in box (optional but implied by filename)
        u.atoms.translate(-u.atoms.center_of_geometry())
        u.atoms.translate(u.dimensions[:3] / 2)

        coords = u.atoms.positions

        f.write(f"{len(u.atoms)}\n")
        f.write("\n")
        for atom, coord in zip(u.atoms.names, coords):
            f.write(f"{atom} {coord[0]:.10f} {coord[1]:.10f} {coord[2]:.10f}\n")

print('Finished.')
