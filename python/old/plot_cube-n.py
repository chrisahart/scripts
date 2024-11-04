import numpy as np
from matplotlib import pyplot as plt
from ase.io.cube import read_cube_data
from scripts.general import parameters as param

# CP2K Li chain
folder_save = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/iv/li/cp2k/transmission/dev-chris/single-points/V-0_HLB-F_z-0-0'
files = [
    '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/iv/li/cp2k/transmission/dev-chris/single-points/V-0_HLB-F_z-0-0/0V-ELECTRON_DENSITY-1_0.cube',
    '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/iv/li/cp2k/transmission/dev-chris/single-points/V-0_HLB-F_z-0-0/0V-v_hartree-1_0.cube',
    '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/iv/li/cp2k/transmission/dev-chris/single-points/V-0_HLB-F_z-0-0/dft_wfn-ELECTRON_DENSITY-1_0.cube',
    '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/iv/li/cp2k/transmission/dev-chris/single-points/V-0_HLB-F_z-0-0/dft_wfn-v_hartree-1_0.cube',
    '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/iv/li/cp2k/transmission/dev-chris/single-points/V-0_HLB-F_z-0-0/bulk-ELECTRON_DENSITY-1_0.cube',
    '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/iv/li/cp2k/transmission/dev-chris/single-points/V-0_HLB-F_z-0-0/bulk-v_hartree-1_0.cube']
label = ['1', '2', '3', '4', '5', '6']

# Read .cube using ASE
cube = []
z_average = []
for i in range(len(files)):
    cube.append(read_cube_data(files))
    # z_average.append(np.zeros(cube[i][].shape[2]))

print(cube)

if __name__ == "__main__":
    print('Finished.')
    plt.show()
