from __future__ import division, print_function
import pandas as pd
import numpy as np
import glob
from general import load_coordinates
from general import print_xyz
import matplotlib.pyplot as plt
from general import parameters as param

"""
    Plot SMEAGOL file IV.SystemLabel_TRC_Channels_K.dat with k-dependent number of channels
"""

# SIESTA+SMEAGOL Au capacitor
folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/siesta-smeagol/capacitor/testing/other/kpoints-2-2-20_hlb-auto_cores-64_au_TransmissionOverk_EM.TRCChannels'
num_kpoints = 4
reshape = [2, 2]
# folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/siesta-smeagol/capacitor/testing/other/kpoints-4-4-20_hlb-auto_cores-64_au_TransmissionOverk_EM.TRCChannels'
folder_save = folder
filename = '0.transport_TRC_Channels_K'
# filename = '0.V_TRC_Channels_K'
# num_kpoints = 12
# reshape = [3, 4]

# CP2K+SMEAGOL Au capacitor
# folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-capacitor/layers-1-2-3-4/cp2k/weightrho/weightrho/testing/kpoints-4-4-20_hlb-auto_NEnergReal-64_TransmissionOverk_EM.TRCChannels'
# folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-capacitor/layers-1-2-3-4/cp2k/weightrho/weightrho/testing/kpoints-4-4-20_hlb-auto_NEnergReal-64_TransmissionOverk_EM.TRCChannels_chris'
# num_kpoints = 16
# reshape = [4, 4]
# folder = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-capacitor/layers-1-2-3-4/cp2k/weightrho/weightrho/testing/kpoints-2-2-20_hlb-auto_NEnergReal-64_TransmissionOverk_EM.TRCChannels_chris'
# num_kpoints = 4
# reshape = [2, 2]
# filename = '0.0V_TRC_Channels_K'
# filename = '0.1V_TRC_Channels_K'

# Au capacitor
input_filename = '{}.dat'.format(filename)
folder_save = folder
ylim = [0, 12]

# Colorbar for line graphs
cmap = plt.get_cmap("tab10")
cmap_colors = np.zeros((10, 4))
for i in range(0, 10):
    cmap_colors[i, :] =cmap(i)
cmap_colors = list(cmap_colors) * 10

# Colors for 2D plots
colorbar_color= 'jet'

# Read number of atoms and labels from .xyz file
cols = ['Energy', 'kx', 'ky', 'wk', 'Ttotal', 'T1', 'T2', 'T3', 'T4', 'T5', 'NLk', 'NRk']
file_spec = pd.read_csv('{}/{}'.format(folder, input_filename), names=cols, delim_whitespace=True, skiprows=1)
file_spec = file_spec.apply(pd.to_numeric, errors='coerce')
file_spec = file_spec.dropna(axis='rows', thresh=2)
file_spec = file_spec.reset_index(drop=True)
data_length = int(file_spec.shape[0]/num_kpoints)

# Find Fermi level
# kpoints_list = np.arange()
index_fermi = np.zeros(num_kpoints)
energy_fermi = np.zeros(num_kpoints)
leads_left_fermi = np.zeros(num_kpoints)
leads_right_fermi = np.zeros(num_kpoints)
for i in range(0, num_kpoints):
    index_fermi[i] = np.argmin(np.abs(file_spec['Energy'][data_length*i:data_length*(i+1)]))+(data_length*i)
    energy_fermi[i] = file_spec['Energy'][index_fermi[i]]
    leads_left_fermi[i] = file_spec['NLk'][index_fermi[i]]
    leads_right_fermi[i] = file_spec['NRk'][index_fermi[i]]

# Plot energy against number of channels (left lead)
fig_plot_1, ax_plot_1 = plt.subplots()
for i in range(0, num_kpoints):
    ax_plot_1.plot(file_spec['Energy'][data_length*i:data_length*(i+1)],
                   file_spec['NLk'][data_length*i:data_length*(i+1)],
                   label=i+1)
ax_plot_1.set_xlabel(r'E-E$_{\mathrm{F}}$ (eV)')
ax_plot_1.set_ylabel('Number of channels (left lead)')
ax_plot_1.legend(frameon=True)
fig_plot_1.tight_layout()
fig_plot_1.savefig('{}/{}_energy_channels_left.png'.format(folder_save, filename), dpi=300)

# Plot energy against number of channels (right lead)
fig_plot_2, ax_plot_2 = plt.subplots()
for i in range(0, num_kpoints):
    ax_plot_2.plot(file_spec['Energy'][data_length*i:data_length*(i+1)],
                   file_spec['NRk'][data_length*i:data_length*(i+1)],
                   label=i+1)
ax_plot_2.set_xlabel(r'E-E$_{\mathrm{F}}$ (eV)')
ax_plot_2.set_ylabel('Number of channels (right lead)')
ax_plot_2.legend(frameon=True)
fig_plot_2.tight_layout()
fig_plot_2.savefig('{}/{}_energy_channels_right.png'.format(folder_save, filename), dpi=300)

# Plot energy against number of channels (both leads)
fig_plot_5, ax_plot_5 = plt.subplots()
for i in range(0, num_kpoints):
    ax_plot_5.plot(file_spec['Energy'][data_length * i:data_length * (i + 1)],
                   file_spec['NLk'][data_length * i:data_length * (i + 1)],  '-',
                   label=i + 1, color=tuple(cmap_colors[i]))
for i in range(0, num_kpoints):
    ax_plot_5.plot(file_spec['Energy'][data_length*i:data_length*(i+1)],
                   file_spec['NRk'][data_length*i:data_length*(i+1)], '--',
                   color=tuple(cmap_colors[i]))
ax_plot_5.set_xlabel(r'E-E$_{\mathrm{F}}$ (eV)')
ax_plot_5.set_ylabel('Number of channels')
ax_plot_5.legend(frameon=True)
fig_plot_5.tight_layout()
fig_plot_5.savefig('{}/{}_energy_channels_left_right.png'.format(folder_save, filename), dpi=300)

# Plot number of channels for each kpoint (left lead)
fig_plot_3, ax_plot_3 = plt.subplots()
c = plt.imshow(leads_left_fermi.reshape(reshape[0], reshape[1]), vmin=ylim[0], vmax=ylim[1], cmap=colorbar_color)
fig_plot_3.colorbar(c)
ax_plot_3.tick_params(top=False, bottom=False, left=False, right=False, labelleft=False, labelbottom=False)
fig_plot_3.tight_layout()
fig_plot_3.savefig('{}/{}_kpoints_channels_left.png'.format(folder_save, filename), dpi=300)

# Plot number of channels for each kpoint (left lead)
fig_plot_4, ax_plot_4 = plt.subplots()
c = plt.imshow(leads_right_fermi.reshape(reshape[0], reshape[1]), vmin=ylim[0], vmax=ylim[1], cmap=colorbar_color)
fig_plot_4.colorbar(c)
ax_plot_4.tick_params(top=False, bottom=False, left=False, right=False, labelleft=False, labelbottom=False)
fig_plot_4.tight_layout()
fig_plot_4.savefig('{}/{}_kpoints_channels_right.png'.format(folder_save, filename), dpi=300)

if __name__ == "__main__":
    print('Finished.')
    plt.show()
