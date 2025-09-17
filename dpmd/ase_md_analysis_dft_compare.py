import numpy as np
import matplotlib.pyplot as plt
from ase.io import read, write, Trajectory
import MDAnalysis as mda
from MDAnalysis.analysis import distances
from MDAnalysis.analysis import rdf
from general import parameters as param


# Rutile 336 PBE + U
topology_file = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/deepmd/rutile/336/md-cell-opt/system.xyz'
folder_all = ['/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-336/md-cell-opt/electron-u-ti-3.0-300k-rs',
              '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/deepmd/rutile/336/md-cell-opt/deepmd-md/pbe-u-4-ps-v2/old/single-fit-m-dpa3-nlayers-6-new2-spin-v3']
trajectory_file_all = ['{}/tio2-pos-1.xyz'.format(folder_all[0])], '{}/md.xyz'.format(folder_all[1]),
start_frame = [6000, 0]
num_atoms = 324
box_size = [13.77, 13.77, 17.76, 90, 90, 90]

calc_distance = True
draw_legend = False
colors_1 = ['r', 'g']
colors_2 = ['b', 'm']
labels_1 = ['DFT Ti (polaron) - O', 'ML Ti (polaron) - O']
labels_2 = ['DFT Ti - O', 'ML Ti - O']

# Plot RDF for Ti - O
nbins = 200
rdf_xlim = 13.77 / 2
fig_rdf_ti_o, ax_rdf_ti_o = plt.subplots(figsize=(10, 4))

for i in range(len(folder_all)):
    print(i)
    folder = folder_all[i]
    print(folder)
    trajectory_file = trajectory_file_all[i]

    universe = mda.Universe(topology_file, trajectory_file)
    universe.dimensions = box_size
    atoms_ti = universe.select_atoms('name Ti')
    num_atoms_ti = len(atoms_ti)
    atoms_o = universe.select_atoms('name O')
    num_atoms_o = len(atoms_o)

    ti_polaron = universe.select_atoms("bynum 78")
    atoms_o = universe.select_atoms("element O")
    rdf_ti_o_polaron = rdf.InterRDF(ti_polaron, atoms_o, range=(0, rdf_xlim), nbins=nbins, start=start_frame[i])
    rdf_ti_o_polaron.run()

    ax_rdf_ti_o.plot(rdf_ti_o_polaron.bins, rdf_ti_o_polaron.rdf, '-', color=colors_1[i], label=labels_1[i])

for i in range(len(folder_all)):
    print(i)
    folder = folder_all[i]
    print(folder)
    trajectory_file = trajectory_file_all[i]

    universe = mda.Universe(topology_file, trajectory_file)
    universe.dimensions = box_size
    atoms_ti = universe.select_atoms('name Ti')
    num_atoms_ti = len(atoms_ti)
    atoms_o = universe.select_atoms('name O')
    num_atoms_o = len(atoms_o)

    ti_polaron = universe.select_atoms("bynum 78")
    atoms_o = universe.select_atoms("element O")
    ti_not_polaron = universe.select_atoms("element Ti and not bynum 78")
    rdf_ti_o_not_polaron = rdf.InterRDF(ti_not_polaron, atoms_o, range=(0, rdf_xlim), nbins=nbins, start=start_frame[i])
    rdf_ti_o_not_polaron.run()

    ax_rdf_ti_o.plot(rdf_ti_o_not_polaron.bins, rdf_ti_o_not_polaron.rdf, '-', color=colors_2[i], label=labels_2[i])

ax_rdf_ti_o.set_xlabel("Radial distance / Å")
ax_rdf_ti_o.set_ylabel("RDF (arb. units)")
ax_rdf_ti_o.legend(frameon=False)
ax_rdf_ti_o.set_xlim([0, rdf_xlim])
ax_rdf_ti_o.tick_params(axis='y', which='both', left=False, right=False, labelleft=False)
fig_rdf_ti_o.tight_layout()

for i in range(len(folder_all)):
    print(i)
    folder = folder_all[i]
    fig_rdf_ti_o.savefig("{}/rdf.png".format(folder), dpi=600)

if __name__ == "__main__":
    print('Finished.')
    plt.show()
