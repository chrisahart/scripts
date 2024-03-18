import numpy as np
import matplotlib.pyplot as plt
from general import parameters as param
from general import load_energy as load_energy

""" Scaling """

folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/AuH2/transport/cp2k-smeagol/archer/archer_pbe/md/hx1/au-6s-only/md/production/'
energy_kinetic_1, energy_potential_1, energy_total_1, temperature_1, time_val_1, time_per_step_1 = load_energy.load_values_energy(folder_1, 'V-0.0/0.0V-1.ener')
energy_kinetic_2, energy_potential_2, energy_total_2, temperature_2, time_val_2, time_per_step_2 = load_energy.load_values_energy(folder_1, 'V-0.1/0.1V-1.ener')
energy_kinetic_3, energy_potential_3, energy_total_3, temperature_3, time_val_3, time_per_step_3 = load_energy.load_values_energy(folder_1, 'V-0.1_dc/0.1V-1.ener')
energy_kinetic_4, energy_potential_4, energy_total_4, temperature_4, time_val_4, time_per_step_4 = load_energy.load_values_energy(folder_1, 'V-0.5/0.5V-1.ener')
energy_kinetic_5, energy_potential_5, energy_total_5, temperature_5, time_val_5, time_per_step_5 = load_energy.load_values_energy(folder_1, 'V-0.5_dc/0.5V-1.ener')
energy_kinetic_6, energy_potential_6, energy_total_6, temperature_6, time_val_6, time_per_step_6 = load_energy.load_values_energy(folder_1, 'V-1/1.0V-1.ener')
energy_kinetic_7, energy_potential_7, energy_total_7, temperature_7, time_val_7, time_per_step_7 = load_energy.load_values_energy(folder_1, 'V-1_dc/1.0V-1.ener')
energy_kinetic_8, energy_potential_8, energy_total_8, temperature_8, time_val_8, time_per_step_8 = load_energy.load_values_energy(folder_1, 'V-2/2.0V-1.ener')
energy_kinetic_9, energy_potential_9, energy_total_9, temperature_9, time_val_9, time_per_step_9 = load_energy.load_values_energy(folder_1, 'V-2_dc/2.0V-1.ener')
energy_kinetic_10, energy_potential_10, energy_total_10, temperature_10, time_val_10, time_per_step_10 = load_energy.load_values_energy(folder_1, 'V-4/4.0V-1.ener')
energy_kinetic_11, energy_potential_11, energy_total_11, temperature_11, time_val_11, time_per_step_11 = load_energy.load_values_energy(folder_1, 'V-4_dc/4.0V-1.ener')
energy_kinetic_12, energy_potential_12, energy_total_12, temperature_12, time_val_12, time_per_step_12 = load_energy.load_values_energy(folder_1, 'V-6/6.0V-1.ener')
energy_kinetic_13, energy_potential_13, energy_total_13, temperature_13, time_val_13, time_per_step_13 = load_energy.load_values_energy(folder_1, 'V-6_dc/6.0V-1.ener')
energy_kinetic_14, energy_potential_14, energy_total_14, temperature_14, time_val_14, time_per_step_14 = load_energy.load_values_energy(folder_1, 'V-8/8.0V-1.ener')
energy_kinetic_15, energy_potential_15, energy_total_15, temperature_15, time_val_15, time_per_step_15 = load_energy.load_values_energy(folder_1, 'V-8_dc/8.0V-1.ener')
energy_kinetic_16, energy_potential_16, energy_total_16, temperature_16, time_val_16, time_per_step_16 = load_energy.load_values_energy(folder_1, 'V-10/10.0V-1.ener')
energy_kinetic_17, energy_potential_17, energy_total_17, temperature_17, time_val_17, time_per_step_17 = load_energy.load_values_energy(folder_1, 'V-10_dc/10.0V-1.ener')
energy_kinetic_18, energy_potential_18, energy_total_18, temperature_18, time_val_18, time_per_step_18 = load_energy.load_values_energy(folder_1, 'dft/dft_wfn-1.ener')


labels = ['V=0.0 SC', 'V=0.1 SC', 'V=0.1 DC', 'V=0.5 SC', 'V=0.5 DC', 'V=1.0 SC', 'V=1.0 DC', 'V=2.0 SC',
          'V=2.0 DC', 'V=4.0 SC', 'V=4.0 DC', 'V=6.0 SC', 'V=6.0 DC', 'V=8.0 SC', 'V=8.0 DC', 'V=10.0 SC',
          'V=10.0 DC', 'DFT']

plotting_colors = ['k', 'b', 'r', 'g', 'm', 'grey', 'orange', 'y', 'brown', 'cyan', 'pink'] * 5
# plotting_colors = ['r', 'g', 'b', 'm', 'grey', 'orange', 'brown']

# Plot total energy all
atoms = 113
scale_x = 1
scale_y = 1/atoms * 1e6
x_lim = [0, 1000]
y_lim = np.array([-3e-6, 1.2e-5]) * 1e6
y_lim = np.array([-200, 500])
fig_plot_energy_all, ax_plot_energy_all = plt.subplots(figsize=(16, 6))
ax_plot_energy_all.plot((time_val_18-time_val_18[0])*scale_x, (energy_total_18-energy_total_18[0])*scale_y, '-', color=plotting_colors[0], label=labels[17])
ax_plot_energy_all.plot((time_val_1-time_val_1[0])*scale_x, (energy_total_1-energy_total_1[0])*scale_y, '-', color=plotting_colors[17], label=labels[0])
ax_plot_energy_all.plot((time_val_2-time_val_2[0])*scale_x, (energy_total_2-energy_total_2[0])*scale_y, '-', color=plotting_colors[1], label=labels[1])
ax_plot_energy_all.plot((time_val_3-time_val_3[0])*scale_x, (energy_total_3-energy_total_3[0])*scale_y, '-', color=plotting_colors[2], label=labels[2])
ax_plot_energy_all.plot((time_val_4-time_val_4[0])*scale_x, (energy_total_4-energy_total_4[0])*scale_y, '-', color=plotting_colors[3], label=labels[3])
ax_plot_energy_all.plot((time_val_5-time_val_5[0])*scale_x, (energy_total_5-energy_total_5[0])*scale_y, '-', color=plotting_colors[4], label=labels[4])
ax_plot_energy_all.plot((time_val_6-time_val_6[0])*scale_x, (energy_total_6-energy_total_6[0])*scale_y, '-', color=plotting_colors[5], label=labels[5])
ax_plot_energy_all.plot((time_val_7-time_val_7[0])*scale_x, (energy_total_7-energy_total_7[0])*scale_y, '-', color=plotting_colors[6], label=labels[6])
ax_plot_energy_all.plot((time_val_8-time_val_8[0])*scale_x, (energy_total_8-energy_total_8[0])*scale_y, '-', color=plotting_colors[7], label=labels[7])
ax_plot_energy_all.plot((time_val_9-time_val_9[0])*scale_x, (energy_total_9-energy_total_9[0])*scale_y, '-', color=plotting_colors[8], label=labels[8])
ax_plot_energy_all.plot((time_val_10-time_val_10[0])*scale_x, (energy_total_10-energy_total_10[0])*scale_y, '-', color=plotting_colors[9], label=labels[9])
ax_plot_energy_all.plot((time_val_11-time_val_11[0])*scale_x, (energy_total_11-energy_total_11[0])*scale_y, '-', color=plotting_colors[10], label=labels[10])
# ax_plot_energy_all.plot((time_val_12-time_val_12[0])*scale_x, (energy_total_12-energy_total_12[0])*scale_y, '-', color=plotting_colors[11], label=labels[11])
ax_plot_energy_all.plot((time_val_13-time_val_13[0])*scale_x, (energy_total_13-energy_total_13[0])*scale_y, '-', color=plotting_colors[12], label=labels[12])
# ax_plot_energy_all.plot((time_val_14-time_val_14[0])*scale_x, (energy_total_14-energy_total_14[0])*scale_y, '-', color=plotting_colors[13], label=labels[13])
# ax_plot_energy_all.plot((time_val_15-time_val_15[0])*scale_x, (energy_total_15-energy_total_15[0])*scale_y, '-', color=plotting_colors[14], label=labels[14])
# ax_plot_energy_all.plot((time_val_16-time_val_16[0])*scale_x, (energy_total_16-energy_total_16[0])*scale_y, '-', color=plotting_colors[15], label=labels[15])
# ax_plot_energy_all.plot((time_val_17-time_val_17[0])*scale_x, (energy_total_17-energy_total_17[0])*scale_y, '-', color=plotting_colors[16], label=labels[16])
ax_plot_energy_all.legend(frameon=False)
ax_plot_energy_all.set_xlim([x_lim[0], x_lim[1]])
ax_plot_energy_all.set_ylim([y_lim[0], y_lim[1]])
ax_plot_energy_all.set_xlabel('Time / fs')
ax_plot_energy_all.set_ylabel('Energy drift per atom / µHa')
fig_plot_energy_all.tight_layout()
fig_plot_energy_all.savefig('{}/energy_all.png'.format(folder_1), dpi=param.save_dpi)

# Plot total energy all
atoms = 113
scale_x = 1
scale_y = 1/atoms * 1e6
x_lim = [0, 1000]
y_lim = np.array([-3e-6, 1.2e-5]) * 1e6
y_lim = np.array([-4, 4])
fig_plot_energy_all_2, ax_plot_energy_all_2 = plt.subplots(figsize=(16, 6))
ax_plot_energy_all_2.plot((time_val_18-time_val_18[0])*scale_x, (energy_total_18-energy_total_18[0])*scale_y, '-', color=plotting_colors[0], label=labels[17])
ax_plot_energy_all_2.plot((time_val_1-time_val_1[0])*scale_x, (energy_total_1-energy_total_1[0])*scale_y, '-', color=plotting_colors[17], label=labels[0])
ax_plot_energy_all_2.plot((time_val_2-time_val_2[0])*scale_x, (energy_total_2-energy_total_2[0])*scale_y, '-', color=plotting_colors[1], label=labels[1])
# ax_plot_energy_all_2.plot((time_val_3-time_val_3[0])*scale_x, (energy_total_3-energy_total_3[0])*scale_y, '-', color=plotting_colors[2], label=labels[2])
ax_plot_energy_all_2.plot((time_val_4-time_val_4[0])*scale_x, (energy_total_4-energy_total_4[0])*scale_y, '-', color=plotting_colors[3], label=labels[3])
# ax_plot_energy_all_2.plot((time_val_5-time_val_5[0])*scale_x, (energy_total_5-energy_total_5[0])*scale_y, '-', color=plotting_colors[4], label=labels[4])
ax_plot_energy_all_2.plot((time_val_6-time_val_6[0])*scale_x, (energy_total_6-energy_total_6[0])*scale_y, '-', color=plotting_colors[5], label=labels[5])
# ax_plot_energy_all_2.plot((time_val_7-time_val_7[0])*scale_x, (energy_total_7-energy_total_7[0])*scale_y, '-', color=plotting_colors[6], label=labels[6])
# ax_plot_energy_all_2.plot((time_val_8-time_val_8[0])*scale_x, (energy_total_8-energy_total_8[0])*scale_y, '-', color=plotting_colors[7], label=labels[7])
# ax_plot_energy_all_2.plot((time_val_9-time_val_9[0])*scale_x, (energy_total_9-energy_total_9[0])*scale_y, '-', color=plotting_colors[8], label=labels[8])
# ax_plot_energy_all_2.plot((time_val_10-time_val_10[0])*scale_x, (energy_total_10-energy_total_10[0])*scale_y, '-', color=plotting_colors[9], label=labels[9])
# ax_plot_energy_all_2.plot((time_val_11-time_val_11[0])*scale_x, (energy_total_11-energy_total_11[0])*scale_y, '-', color=plotting_colors[10], label=labels[10])
# ax_plot_energy_all_2.plot((time_val_12-time_val_12[0])*scale_x, (energy_total_12-energy_total_12[0])*scale_y, '-', color=plotting_colors[11], label=labels[11])
# ax_plot_energy_all_2.plot((time_val_13-time_val_13[0])*scale_x, (energy_total_13-energy_total_13[0])*scale_y, '-', color=plotting_colors[12], label=labels[12])
# ax_plot_energy_all_2.plot((time_val_14-time_val_14[0])*scale_x, (energy_total_14-energy_total_14[0])*scale_y, '-', color=plotting_colors[13], label=labels[13])
# ax_plot_energy_all_2.plot((time_val_15-time_val_15[0])*scale_x, (energy_total_15-energy_total_15[0])*scale_y, '-', color=plotting_colors[14], label=labels[14])
# ax_plot_energy_all_2.plot((time_val_16-time_val_16[0])*scale_x, (energy_total_16-energy_total_16[0])*scale_y, '-', color=plotting_colors[15], label=labels[15])
# ax_plot_energy_all_2.plot((time_val_17-time_val_17[0])*scale_x, (energy_total_17-energy_total_17[0])*scale_y, '-', color=plotting_colors[16], label=labels[16])
ax_plot_energy_all_2.legend(frameon=False)
ax_plot_energy_all_2.set_xlim([x_lim[0], x_lim[1]])
ax_plot_energy_all_2.set_ylim([y_lim[0], y_lim[1]])
ax_plot_energy_all_2.set_xlabel('Time / fs')
ax_plot_energy_all_2.set_ylabel('Energy drift per atom / µHa')
fig_plot_energy_all_2.tight_layout()
fig_plot_energy_all_2.savefig('{}/energy_all_2.png'.format(folder_1), dpi=param.save_dpi)

if __name__ == "__main__":
    print('Finished.')
    plt.show()
