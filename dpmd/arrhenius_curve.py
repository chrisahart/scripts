import numpy as np
import matplotlib.pyplot as plt
from general import parameters as param

folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/tio2/rutile/deepmd/rutile/336/md-cell-opt/deepmd-md/hse-22-v2'
rate_constant = np.array([2899999999999.9995, 2200000000000.0, 1100000000000.0])
mobility_analytical = np.array([0.049421577096343476, 0.05594334565322291, 0.03711935362092547])
temperature = np.array([300, 200, 150])

# Plot rate constant vs temperature
fig_rate_constant, ax_rate_constant = plt.subplots(figsize=(6, 6))
ax_rate_constant.plot(temperature, rate_constant, 'kx-')
ax_rate_constant.set_xlabel("Temperature / K")
ax_rate_constant.set_ylabel("Rate constant / s$^{-1}$")
ax_rate_constant.set_xlim([150-10, 300+10])
fig_rate_constant.tight_layout()
fig_rate_constant.savefig("{}/rate_constant.png".format(folder), dpi=600)

# Arrhenius plot
inv_T = 1 / temperature
ln_k = np.log(rate_constant)
coeffs = np.polyfit(inv_T, ln_k, 1)
ln_k_fit = np.polyval(coeffs, inv_T)
k_au = 8.617333262145E-5  # Kb in eV
Ea = -coeffs[0] * k_au * 1000  # Ea in meV
text_str = r'$\mathrm{{E_A}} = {0:.0f}\ \mathrm{{meV}}$'.format(Ea)

fig_rate_constant_log, ax_rate_constant_log = plt.subplots(figsize=(6, 6))
ax_rate_constant_log.plot(1/temperature, np.log(rate_constant), 'kx')
ax_rate_constant_log.plot(inv_T, ln_k_fit, 'k-', alpha=0.7, label=f'Fit: slope = {coeffs[0]:.2f}')
ax_rate_constant_log.text(0.95, 0.95, text_str, transform=ax_rate_constant_log.transAxes,
                          fontsize=12, color='black', weight='bold', ha='right', va='top')
ax_rate_constant_log.set_xlabel("1 / Temperature (1 / K)")
ax_rate_constant_log.set_ylabel(r"ln [Rate constant (s$^{-1}$)] ")
fig_rate_constant_log.tight_layout()
fig_rate_constant_log.savefig("{}/arrhenius.png".format(folder), dpi=600)

# Plot mobility vs temperature
fig_mobility, ax_mobility = plt.subplots(figsize=(6, 6))
ax_mobility.plot(temperature, mobility_analytical, 'kx-')
ax_mobility.set_xlabel("Temperature / K")
ax_mobility.set_ylabel(r"Mobility / cm$^2$ V$^{-1}$ s$^{-1}$")
fig_mobility.tight_layout()
fig_mobility.savefig("{}/mobility.png".format(folder), dpi=600)

# Plot log(mobility) vs temperature
fig_mobility_log, ax_mobility_log = plt.subplots(figsize=(6, 6))
ax_mobility_log.plot(temperature, np.log10(mobility_analytical), 'kx-')
ax_mobility_log.set_xlabel("Temperature (K)")
ax_mobility_log.set_ylabel(r"Log [mobility (cm$^2$/Vs)]")
fig_mobility_log.tight_layout()
fig_mobility_log.savefig("{}/mobility_log.png".format(folder), dpi=600)

if __name__ == "__main__":
    print('Finished.')
    plt.show()
