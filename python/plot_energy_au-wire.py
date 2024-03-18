import numpy as np
import matplotlib.pyplot as plt
from general import parameters as param
from general import load_energy as load_energy

""" Scaling """

folder_1 = '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-wire-ismael/jad/calculations/md/restart-auto-long_au-frozen/'
energy_kinetic_1, energy_potential_1, energy_total_1, temperature_1, time_val_1, time_per_step_1 = load_energy.load_values_energy(folder_1, 'md_dft_long_LINEAR_P/2_dft_wfn-1.ener')
energy_kinetic_2, energy_potential_2, energy_total_2, temperature_2, time_val_2, time_per_step_2 = load_energy.load_values_energy(folder_1, 'md_dft_long_LINEAR_P_NVE/2_dft_wfn-1.ener')
energy_kinetic_3, energy_potential_3, energy_total_3, temperature_3, time_val_3, time_per_step_3 = load_energy.load_values_energy(folder_1, 'md_V-0.0_long_USE_GUESS-failed-2237/0.0V-1.ener')
energy_kinetic_4, energy_potential_4, energy_total_4, temperature_4, time_val_4, time_per_step_4 = load_energy.load_values_energy(folder_1, 'md_V-0.0_long_USE_GUESS_NVE-failed-2238/0.0V-1.ener')
energy_kinetic_5, energy_potential_5, energy_total_5, temperature_5, time_val_5, time_per_step_5 = load_energy.load_values_energy(folder_1, 'md_V-0.1_long_USE_GUESS/0.1V-1.ener')
energy_kinetic_6, energy_potential_6, energy_total_6, temperature_6, time_val_6, time_per_step_6 = load_energy.load_values_energy(folder_1, 'md_V-1.0_long_USE_GUESS-failed-2358/1.0V-1.ener')
energy_kinetic_7, energy_potential_7, energy_total_7, temperature_7, time_val_7, time_per_step_7 = load_energy.load_values_energy(folder_1, 'md_V-1.0_long_USE_GUESS_NVE/1.0V-1.ener')

# plot_color = ['k', 'b', 'r', 'g', 'm', 'grey', 'orange', 'y', 'brown', 'cyan', 'pink']
plotting_colors = ['k', 'b', 'r', 'g', 'm', 'grey', 'orange', 'y', 'brown', 'cyan', 'pink'] * 5
# plotting_colors = ['r', 'g', 'b', 'm', 'grey', 'orange', 'brown']

# Plot total energy all
scale_x = 1
scale_y = 1/1338 * 1e6
x_lim = [0, 100]
y_lim = np.array([-3e-6, 1.2e-5]) * 1e6
fig_plot_energy_all, ax_plot_energy_all = plt.subplots()
# ax_plot_energy_all.plot((time_val_1-time_val_1[0])*scale_x, (energy_total_1-energy_total_1[0])*scale_y, '-', color=plotting_colors[0], label='CP2K NVT')
# ax_plot_energy_all.plot((time_val_2-time_val_2[0])*scale_x, (energy_total_2-energy_total_2[0])*scale_y, '-', color=plotting_colors[1], label='CP2K NVE')
# ax_plot_energy_all.plot((time_val_3-time_val_3[0])*scale_x, (energy_total_3-energy_total_3[0])*scale_y, '-', color=plotting_colors[2], label='CP2K+SMEAGOL NVT V=0')
# ax_plot_energy_all.plot((time_val_4-time_val_4[0])*scale_x, (energy_total_4-energy_total_4[0])*scale_y, '-', color=plotting_colors[3], label='CP2K+SMEAGOL NVE V=0')
# ax_plot_energy_all.plot((time_val_5-time_val_5[0])*scale_x, (energy_total_5-energy_total_5[0])*scale_y, '-', color=plotting_colors[4], label='CP2K+SMEAGOL NVT V=0.1')
# ax_plot_energy_all.plot((time_val_6-time_val_6[0])*scale_x, (energy_total_6-energy_total_6[0])*scale_y, '-', color=plotting_colors[5], label='CP2K+SMEAGOL NVT V=1')
# ax_plot_energy_all.plot((time_val_7-time_val_7[0])*scale_x, (energy_total_7-energy_total_7[0])*scale_y, '-', color=plotting_colors[6], label='CP2K+SMEAGOL NVE V=1')
ax_plot_energy_all.plot((time_val_2-time_val_2[0])*scale_x, (energy_total_2-energy_total_2[0])*scale_y, '-', color=plotting_colors[0], label='CP2K')
ax_plot_energy_all.plot((time_val_4-time_val_4[0])*scale_x, (energy_total_4-energy_total_4[0])*scale_y, '-', color=plotting_colors[3], label='CP2K+SMEAGOL V=0.0')
ax_plot_energy_all.plot((time_val_5-time_val_5[0])*scale_x, (energy_total_5-energy_total_5[0])*scale_y, '-', color=plotting_colors[4], label='CP2K+SMEAGOL V=0.1')
ax_plot_energy_all.plot((time_val_7-time_val_7[0])*scale_x, (energy_total_7-energy_total_7[0])*scale_y, '-', color=plotting_colors[6], label='CP2K+SMEAGOL V=1.0')
ax_plot_energy_all.legend(frameon=False)
ax_plot_energy_all.set_xlim([x_lim[0], x_lim[1]])
ax_plot_energy_all.set_ylim([y_lim[0], y_lim[1]])
ax_plot_energy_all.set_xlabel('Time / fs')
ax_plot_energy_all.set_ylabel('Energy drift per atom / µHa')
fig_plot_energy_all.tight_layout()
fig_plot_energy_all.savefig('{}/energy_all.png'.format(folder_1), dpi=param.save_dpi)

# Plot temperature all
scale_x = 1
scale_y = 1
x_lim_2 = x_lim
y_lim_2 = [280, 340]
fig_plot_temperature_all, ax_plot_temperature_all = plt.subplots()
# ax_plot_temperature_all.plot((time_val_1-time_val_1[0])*scale_x, temperature_1*scale_y, '-', color=plotting_colors[0], label='CP2K NVT')
# ax_plot_temperature_all.plot((time_val_2-time_val_2[0])*scale_x, temperature_2*scale_y, '-', color=plotting_colors[1], label='CP2K NVE')
# ax_plot_temperature_all.plot((time_val_3-time_val_3[0])*scale_x, temperature_3*scale_y, '-', color=plotting_colors[2], label='CP2K+SMEAGOL NVT V=0')
# ax_plot_temperature_all.plot((time_val_4-time_val_4[0])*scale_x, temperature_4*scale_y, '-', color=plotting_colors[3], label='CP2K+SMEAGOL NVE V=0')
# ax_plot_temperature_all.plot((time_val_5-time_val_5[0])*scale_x, temperature_5*scale_y, '-', color=plotting_colors[4], label='CP2K+SMEAGOL NVT V=0.1')
# ax_plot_temperature_all.plot((time_val_6-time_val_6[0])*scale_x, temperature_6*scale_y, '-', color=plotting_colors[5], label='CP2K+SMEAGOL NVT V=1')
# ax_plot_temperature_all.plot((time_val_7-time_val_7[0])*scale_x, temperature_7*scale_y, '-', color=plotting_colors[6], label='CP2K+SMEAGOL NVE V=1')
ax_plot_temperature_all.plot((time_val_2-time_val_2[0])*scale_x, temperature_2*scale_y, '-', color=plotting_colors[0], label='CP2K')
ax_plot_temperature_all.plot((time_val_4-time_val_4[0])*scale_x, temperature_4*scale_y, '-', color=plotting_colors[3], label='CP2K+SMEAGOL V=0.0')
ax_plot_temperature_all.plot((time_val_5-time_val_5[0])*scale_x, temperature_5*scale_y, '-', color=plotting_colors[4], label='CP2K+SMEAGOL V=0.1')
ax_plot_temperature_all.plot((time_val_7-time_val_7[0])*scale_x, temperature_7*scale_y, '-', color=plotting_colors[6], label='CP2K+SMEAGOL V=1.0')
ax_plot_temperature_all.legend(frameon=False)
ax_plot_temperature_all.set_xlim([x_lim_2[0], x_lim_2[1]])
ax_plot_temperature_all.set_ylim([y_lim_2[0], y_lim_2[1]])
ax_plot_temperature_all.set_xlabel('Time / fs')
ax_plot_temperature_all.set_ylabel('Temperature / K')
fig_plot_temperature_all.tight_layout()
fig_plot_temperature_all.savefig('{}/temperature_all.png'.format(folder_1), dpi=param.save_dpi)

if __name__ == "__main__":
    print('Finished.')
    plt.show()
