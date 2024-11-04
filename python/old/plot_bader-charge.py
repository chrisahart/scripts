import numpy as np
from matplotlib import pyplot as plt
from general import parameters as param
import pandas as pd
from general import load_coordinates
import csv

""" Plotting of CP2K .cube files for Hartree potential and charge density """

# Plotting variables
save_dpi = 600
plot_colorbar = False
colorbar_limit = [0.1, -0.3]
marker_size = 150
plot_charge_text = True
species_to_marker = {"Au_lead": np.NaN, "Au_screen": np.NaN, "Au_al": np.NaN, "Au_bl": np.NaN,
                     "Au_cl": marker_size, "H": np.NaN, "O": np.NaN, "Au_wire": marker_size, "Au_br": marker_size,
                     "Au_cr": np.NaN}
# colorbar_limit = [1, -1]
# marker_size = 10
# plot_charge_text = False
# species_to_marker = {"Au_lead": marker_size, "Au_screen": marker_size, "Au_al": marker_size, "Au_bl": marker_size,
#                      "Au_cl": marker_size,   "H": np.NaN, "O": np.NaN, "Au_wire": marker_size,  "Au_br": marker_size,
#                      "Au_cr": marker_size}

# Structure
files_structure = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-wire-ismael/jad/calculations/md/restart-auto-long_full/md_V-1.5_long'
cols = ['Species', 'X', 'Y', 'Z']
file_coord, num_atoms_1, species_1 = load_coordinates.load_file_coord(files_structure, 'em.xyz', cols)
file_coord = file_coord.reset_index(drop=True)

# Bader
array_time = ['0', '50', '100', '180']
plot_time = array_time[0]
array_bias = ['2_dft_wfn', '0.0V', '0.1V', '1.0V']
plot_bias = array_bias[1]
plot_au_species = 'Au_wire'
# plot_au = [30, 44]  # Au wire and first atom of tip
plot_au = [25, 30]  # Au tip left
# plot_au = [44, 46]  # Au tip right
# plot_au_species = 'Au_cl'
# plot_au = [20, 30]  # Au slab left
# plot_au_species = 'Au_br'
# plot_au = [40, 50]  # Au slab right
folder_bader = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-wire-ismael/jad/calculations/single-points/snapshots/from-md-cube/md_{}-fs/charge/{}'.format(plot_time, plot_bias)
files_bader = '{}/ACF.dat'.format(folder_bader)
filename_out_bader = '{}/bader_charges.out'.format(folder_bader)
filename_save_bader = '{}/bader_charges_1.png'.format(folder_bader)
cols = ['#', 'X', 'Y', 'Z', 'CHARGE', 'MIN DIST', 'ATOMIC VOL']
file_bader_1 = pd.read_csv(files_bader, names=cols, delim_whitespace=True, skiprows=2)
file_bader_1 = file_bader_1.apply(pd.to_numeric, errors='coerce')
file_bader_1 = file_bader_1.reset_index(drop=True)
print(file_bader_1)

# Hirshfeld
plot_hirsh = True
skip_rows = 2
folder_hirsh = folder_bader
files_hirsh = '{}/{}-charges-1.hirshfeld'.format(folder_hirsh, plot_bias)
filename_save_hirsh = '{}/hirshfeld_charges_1.png'.format(folder_hirsh)
cols_hirsh = ['Atom', 'Element', 'Kind', 'Ref Charge', 'Population', 'Net charge']
file_hirsh_1 = pd.read_csv(files_hirsh, names=cols_hirsh, delim_whitespace=True, on_bad_lines='skip', skiprows=skip_rows)
file_hirsh_1 = file_hirsh_1.apply(pd.to_numeric, errors='coerce')
file_hirsh_1 = file_hirsh_1[file_hirsh_1['Net charge'].notna()]
file_hirsh_1 = file_hirsh_1.reset_index(drop=True)
hirsh_charges = file_hirsh_1['Net charge'][0:num_atoms_1]

# Plotting marker size
marker_array = np.zeros(num_atoms_1)
for i in range(num_atoms_1):
    species = species_1[i]
    if species in species_to_marker:
        marker_array[i] = species_to_marker[species]
    else:
        print('Error', species)

# Bader reference charges
species_to_charge = {"Au_lead": 1, "Au_screen": 1, "Au_al": 11, "Au_bl": 11, "Au_cl": 11,
                     "H": 1, "O": 6, "Au_wire": 11,
                     "Au_br": 11, "Au_cr": 11}
ref_charge = np.zeros(num_atoms_1)
for i in range(num_atoms_1):
    species = species_1[i]
    if species in species_to_charge:
        ref_charge[i] = species_to_charge[species]
    else:
        print('Error', species)

# Setup pandas dataframe
bader_charges = ref_charge-file_bader_1['CHARGE'][0:num_atoms_1]
bader_charges_pandas = pd.DataFrame([], columns=['Species', 'Charge'])
bader_charges_pandas['Species'] = species_1
bader_charges_pandas['Charge'] = bader_charges
bader_charges_pandas.to_csv(filename_out_bader, index=False, header=False, quoting=csv.QUOTE_NONE, sep=" ")

# Printing
print('Number of atoms', num_atoms_1)
print('Max Bader charge', np.max(bader_charges))
print('Min Bader charge', np.min(bader_charges))
print('Max Hirshfeld charge', np.max(hirsh_charges))
print('Min Hirshfeld charge', np.min(hirsh_charges))

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
                      vmin=colorbar_limit[1], vmax=colorbar_limit[0])
store_charges = []
if plot_charge_text:
    for i in range(0, num_atoms_1):
        if species_1[i] == plot_au_species and plot_au[0] < (file_bader_1['Z'][i] / param.angstrom_to_bohr) < plot_au[1]:
            store_charges.append(bader_charges[i])
            ax_bader.text((file_bader_1['X'][i] / param.angstrom_to_bohr) + 0.0,
                          (file_bader_1['Y'][i] / param.angstrom_to_bohr) + -2.0,
                          (file_bader_1['Z'][i] / param.angstrom_to_bohr) + -1.1,
                          '{0:.2f}'.format(bader_charges[i]))
store_charges = np.array(store_charges)
print('Charges:', store_charges)
print('Average charge:', np.sum(store_charges)/np.shape(store_charges)[0])
ax_bader.set_xlabel('X axis')
ax_bader.set_ylabel('Y axis')
ax_bader.set_zlabel('Z axis')
ax_bader.set_proj_type('ortho')
ax_bader.axis("off")
if plot_colorbar: cbar = fig_bader.colorbar(sc, ax=ax_bader, fraction=0.046, pad=0.04)
if plot_colorbar: cbar.set_label('Charge |e|')
fig_bader.tight_layout()
ax_bader.view_init(elev=0, azim=0, roll=-90)
ax_bader.set_box_aspect((1, 1, 3))
fig_bader.tight_layout()
fig_bader.savefig(filename_save_bader, dpi=save_dpi)

# Create a 3D plot of Hirshfeld charges
if plot_hirsh:
    fig_hirsh = plt.figure(figsize=(4, 8))
    ax_hirsh = fig_hirsh.add_subplot(111, projection='3d')
    sc = ax_hirsh.scatter(file_bader_1['X'][0:num_atoms_1] / param.angstrom_to_bohr,
                          file_bader_1['Y'][0:num_atoms_1] / param.angstrom_to_bohr,
                          file_bader_1['Z'][0:num_atoms_1] / param.angstrom_to_bohr,
                          c=hirsh_charges,
                          cmap='viridis',
                          marker='o',
                          s=marker_array,
                          vmin=colorbar_limit[1], vmax=colorbar_limit[0])
    if plot_charge_text:
        for i in range(0, num_atoms_1):
            if species_1[i] == "Au_wire" and 30 < (file_bader_1['Z'][i] / param.angstrom_to_bohr) < 44:
                ax_hirsh.text((file_bader_1['X'][i] / param.angstrom_to_bohr) + 0.0,
                              (file_bader_1['Y'][i] / param.angstrom_to_bohr) + -2.0,
                              (file_bader_1['Z'][i] / param.angstrom_to_bohr) + -1.1,
                              '{0:.2f}'.format(hirsh_charges[i]))
    ax_hirsh.set_xlabel('X axis')
    ax_hirsh.set_ylabel('Y axis')
    ax_hirsh.set_zlabel('Z axis')
    ax_hirsh.set_proj_type('ortho')
    ax_hirsh.axis("off")
    if plot_colorbar: cbar = fig_hirsh.colorbar(sc, ax=ax_hirsh, fraction=0.046, pad=0.04)
    if plot_colorbar: cbar.set_label('Charge |e|')
    fig_hirsh.tight_layout()
    ax_hirsh.view_init(elev=0, azim=0, roll=-90)
    ax_hirsh.set_box_aspect((1, 1, 3))
    fig_hirsh.tight_layout()
    fig_hirsh.savefig(filename_save_hirsh, dpi=save_dpi)

if __name__ == "__main__":
    print('Finished.')
    plt.show()