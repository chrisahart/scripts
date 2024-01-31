import numpy as np
from matplotlib import pyplot as plt
from general import parameters as param
import pandas as pd
from general import load_coordinates
import csv

""" Plotting of CP2K .cube files for Hartree potential and charge density """

files_bader = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-wire-ismael/jad/calculations/single-points/snapshots/from-md-0/bader/ACF.dat'
filename_out_bader = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-wire-ismael/jad/calculations/single-points/snapshots/from-md-0/bader/bader_charges.out'
filename_save = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-wire-ismael/jad/calculations/single-points/snapshots/from-md-0/bader/bader_charges.png'
cols = ['#', 'X', 'Y', 'Z', 'CHARGE', 'MIN DIST', 'ATOMIC VOL']
file_bader_1 = pd.read_csv(files_bader, names=cols, delim_whitespace=True, skiprows=2)
file_bader_1 = file_bader_1.apply(pd.to_numeric, errors='coerce')
file_bader_1 = file_bader_1.reset_index(drop=True)

files_structure = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-wire-ismael/jad/calculations/md/restart-auto-long_full/md_V-1.5_long'
cols = ['Species', 'X', 'Y', 'Z']
file_coord, num_atoms_1, species_1 = load_coordinates.load_file_coord(files_structure, 'em.xyz', cols)
file_coord = file_coord.reset_index(drop=True)

# Read number of atoms and labels from .xyz file
filename = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-wire-ismael/jad/calculations/single-points/snapshots/md-structure-0/dft/2_dft_wfn-charges-1.hirshfeld'
cols_hirsh = ['Atom', 'Element', 'Kind', 'Ref Charge', 'Population', 'Net charge']
file_hirsh_1 = pd.read_csv(filename, names=cols_hirsh, delim_whitespace=True, skiprows=5)
file_hirsh_1 = file_hirsh_1.apply(pd.to_numeric, errors='coerce')
file_hirsh_1 = file_hirsh_1.reset_index(drop=True)
hirsh_charges = file_hirsh_1['Net charge'][0:num_atoms_1]
print('Max Hirshfeld charge', np.max(hirsh_charges))
print('Min Hirshfeld charge', np.min(hirsh_charges))

# print(hirsh_charges)
print('Number of atoms', num_atoms_1)

marker_array = np.zeros(num_atoms_1)
marker_size = 150
for i in range(0, num_atoms_1):
    if species_1[i] == "Au_lead":
        marker_array[i] = marker_size
    elif species_1[i] == "Au_screen":
        marker_array[i] = marker_size
    elif species_1[i] == "Au_al":
        marker_array[i] = marker_size
    elif species_1[i] == "Au_bl":
        marker_array[i] = marker_size
    elif species_1[i] == "Au_cl":
        marker_array[i] = marker_size
    elif species_1[i] == "H":
        marker_array[i] = 0
    elif species_1[i] == "O":
        marker_array[i] = 0
    elif species_1[i] == "Au_wire":
        marker_array[i] = marker_size
    elif species_1[i] == "Au_br":
        marker_array[i] = marker_size
    elif species_1[i] == "Au_cr":
        marker_array[i] = marker_size
    else:
        print('Error', species_1[i])

ref_charge = np.zeros(num_atoms_1)
for i in range(0, num_atoms_1):
    if species_1[i] == "Au_lead":
        ref_charge[i] = np.NaN
    elif species_1[i] == "Au_screen":
        ref_charge[i] = np.NaN
    elif species_1[i] == "Au_al":
        ref_charge[i] = np.NaN
    elif species_1[i] == "Au_bl":
        ref_charge[i] = np.NaN
    elif species_1[i] == "Au_cl":
        ref_charge[i] = 11
    elif species_1[i] == "H":
        ref_charge[i] = np.NaN
    elif species_1[i] == "O":
        ref_charge[i] = np.NaN
    elif species_1[i] == "Au_wire":
        ref_charge[i] = 11
    elif species_1[i] == "Au_br":
        ref_charge[i] = 11
    elif species_1[i] == "Au_cr":
        ref_charge[i] = np.NaN
    else:
        print('Error', species_1[i])

# Setup pandas dataframe
bader_charges = ref_charge-file_bader_1['CHARGE'][0:num_atoms_1]
bader_charges_pandas = pd.DataFrame([], columns=['Species', 'Charge'])
bader_charges_pandas['Species'] = species_1
bader_charges_pandas['Charge'] = bader_charges
print('Max Bader charge', np.max(bader_charges))
print('Min Bader charge', np.min(bader_charges))
lims = [0.1, -0.3]
bader_charges_pandas.to_csv(filename_out_bader, index=False, header=False, quoting=csv.QUOTE_NONE, sep=" ")

print(bader_charges)
print(hirsh_charges)

# Create a 3D plot of Bader charges
fig_bader = plt.figure(figsize=(4, 8))
ax_bader = fig_bader.add_subplot(111, projection='3d')
sc = ax_bader.scatter(file_bader_1['X'][0:num_atoms_1] / param.angstrom_to_bohr,
                      file_bader_1['Y'][0:num_atoms_1] / param.angstrom_to_bohr,
                      file_bader_1['Z'][0:num_atoms_1] / param.angstrom_to_bohr,
                      c=bader_charges,
                      cmap='viridis',
                      marker='o',
                      s=marker_array,
                      vmin=lims[1], vmax=lims[0])
for i in range(0, num_atoms_1):
    if species_1[i] == "Au_wire":
        if abs(bader_charges[i]) > 0.1:
            ax_bader.text((file_bader_1['X'][i] / param.angstrom_to_bohr) + 0.0,
                          (file_bader_1['Y'][i] / param.angstrom_to_bohr) + -2.0,
                          (file_bader_1['Z'][i] / param.angstrom_to_bohr) + -1.1,
                          '{0:.2f}'.format(bader_charges[i]))
ax_bader.set_xlabel('X axis')
ax_bader.set_ylabel('Y axis')
ax_bader.set_zlabel('Z axis')
ax_bader.set_proj_type('ortho')
ax_bader.axis("off")
# cbar = fig_bader.colorbar(sc, ax=ax_bader, fraction=0.046, pad=0.04)
# cbar.set_label('Charge |e|')
fig_bader.tight_layout()
ax_bader.view_init(elev=0, azim=0, roll=-90)
ax_bader.set_box_aspect((1, 1, 3))
fig_bader.tight_layout()
fig_bader.savefig(filename_save, dpi=600)

# Create a 3D plot of Hirshfeld charges
fig_hirsh = plt.figure(figsize=(4, 8))
ax_hirsh = fig_hirsh.add_subplot(111, projection='3d')
sc = ax_hirsh.scatter(file_bader_1['X'][0:num_atoms_1] / param.angstrom_to_bohr,
                      file_bader_1['Y'][0:num_atoms_1] / param.angstrom_to_bohr,
                      file_bader_1['Z'][0:num_atoms_1] / param.angstrom_to_bohr,
                      c=hirsh_charges,
                      cmap='viridis',
                      marker='o',
                      s=marker_array,
                      vmin=lims[1], vmax=lims[0])
for i in range(0, num_atoms_1):
    if species_1[i] == "Au_wire":
        if abs(bader_charges[i]) > 0.1:
            ax_hirsh.text((file_bader_1['X'][i] / param.angstrom_to_bohr) + 0.0,
                          (file_bader_1['Y'][i] / param.angstrom_to_bohr) + -2.0,
                          (file_bader_1['Z'][i] / param.angstrom_to_bohr) + -1.1,
                          '{0:.2f}'.format(hirsh_charges[i]))
ax_hirsh.set_xlabel('X axis')
ax_hirsh.set_ylabel('Y axis')
ax_hirsh.set_zlabel('Z axis')
ax_hirsh.set_proj_type('ortho')
ax_hirsh.axis("off")
# cbar = fig_hirsh.colorbar(sc, ax=ax_hirsh, fraction=0.046, pad=0.04)
# cbar.set_label('Charge |e|')
fig_hirsh.tight_layout()
ax_hirsh.view_init(elev=0, azim=0, roll=-90)
ax_hirsh.set_box_aspect((1, 1, 3))
fig_hirsh.tight_layout()
fig_hirsh.savefig(filename_save, dpi=600)

if __name__ == "__main__":
    print('Finished.')
    plt.show()