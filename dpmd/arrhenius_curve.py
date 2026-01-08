import numpy as np
import matplotlib.pyplot as plt
from general import parameters as param

folder = '/Volumes/Samsung/Data/Postdoc2/Data/Work/calculations/tio2/rutile/deepmd/rutile/336/md-cell-opt/deepmd-md/hse-22-v2'
rate_constant = np.array([2899999999999.9995, 2200000000000.0, 1100000000000.0])
mobility = np.array([0.049421577096343476, 0.05594334565322291, 0.03711935362092547])
temperature = np.array([300, 200, 150])

folder2 = '/Volumes/Samsung/Data/Postdoc2/Data/Work/calculations/tio2/rutile/deepmd/rutile/336/md-cell-opt/deepmd-md/hse-22-ts-md2'
# rate_constant2 = np.array([2699909319544.7314, 3162885773028.014, 1607279517009.998, 1043797660327.0316])
# mobility2 = np.array([0.04613585262396951, 0.06410966464081157, 0.040988928201939036, 0.035473713688457635])
# temperature2 = np.array([300, 250, 200, 150])
rate_constant2 = np.array([2699909319544.7314,1607279517009.998, 1043797660327.0316, 86527635758.60518])
mobility2 = np.array([0.04613585262396951, 0.040988928201939036, 0.035473713688457635, 0.004410082150441339])
temperature2 = np.array([300, 200, 150, 100])

folder = folder2
rate_constant = rate_constant2
mobility = mobility2
temperature = temperature2

folder3 = '/Volumes/Samsung/Data/Postdoc2/Data/Work/calculations/tio2/anatase/deepmd/anatase/441/deepmd/hse-19-ts-md-9500-9900-removed'

rate_constant3 = np.array([1463332478658.7344, 672956769482.8204, 241988029193.7231, 43152745244.30003])
mobility3 = np.array([0.011144128372621817, 0.006315264815752413, 0.0028502005634936617, 0.0006294358382421668])
temperature3 = np.array([600, 500, 400, 300])

rate_constant3 = np.array([1422749145062.692, 833346616448.9613, 403902936462.9768, 94455838920.67673])
mobility3 = np.array([0.010993409695600897, 0.007573700430471499, 0.004588796525929469, 0.001414400073473747])
temperature3 = np.array([600, 500, 400, 300])


# folder = folder3
# rate_constant = rate_constant3
# mobility = mobility3
# temperature = temperature3

# Plot rate constant vs temperature
# fig_rate_constant, ax_rate_constant = plt.subplots(figsize=(6, 6))
# ax_rate_constant.plot(temperature, rate_constant, 'kx-')
# ax_rate_constant.plot(temperature2, rate_constant2, 'gx-')
# ax_rate_constant.set_xlabel("Temperature / K")
# ax_rate_constant.set_ylabel("Rate constant / s$^{-1}$")
# ax_rate_constant.set_xlim([150-10, 300+10])
# fig_rate_constant.tight_layout()
# fig_rate_constant.savefig("{}/rate_constant.png".format(folder), dpi=600)

# Arrhenius plot
inv_T = 1 / temperature
ln_k = np.log(rate_constant)
coeffs = np.polyfit(inv_T, ln_k, 1)
ln_k_fit = np.polyval(coeffs, inv_T)
k_au = 8.617333262145E-5  # Kb in eV
Ea = -coeffs[0] * k_au * 1000  # Ea in meV
text_str = r'$\mathrm{{E_A}} = {0:.0f}\ \mathrm{{meV}}$'.format(Ea)

inv_T3 = 1 / temperature3
ln_k3 = np.log(rate_constant3)
coeffs3 = np.polyfit(inv_T3, ln_k3, 1)
ln_k_fit3 = np.polyval(coeffs3, inv_T3)
Ea3 = -coeffs3[0] * k_au * 1000  # Ea in meV
text_str3 = r'$\mathrm{{E_A}} = {0:.0f}\ \mathrm{{meV}}$'.format(Ea3)

# Calculate pre-exponential factor
rutile_freq = rate_constant / np.exp(-(Ea/1000)/(k_au * temperature))
anatase_freq = rate_constant3 / np.exp(-(Ea3/1000)/(k_au * temperature3))
print('rutile_freq', rutile_freq/1e13, np.mean(rutile_freq)/1e13)
print('anatase_freq', anatase_freq/1e13, np.mean(anatase_freq)/1e13)

fig_rate_constant_log, ax_rate_constant_log = plt.subplots(figsize=(6, 6))
ax_rate_constant_log.plot(1/temperature, np.log(rate_constant), 'bx')
ax_rate_constant_log.plot(1/temperature3, np.log(rate_constant3), 'rx')
ax_rate_constant_log.plot(inv_T, ln_k_fit, 'b-', alpha=0.7, label=f'Fit: slope = {coeffs[0]:.2f}')
ax_rate_constant_log.plot(inv_T3, ln_k_fit3, 'r-', alpha=0.7, label=f'Fit: slope = {coeffs3[0]:.2f}')
ax_rate_constant_log.text(0.90, 0.90, text_str, transform=ax_rate_constant_log.transAxes,
                          fontsize=12, color='black', weight='bold', ha='right', va='top')
ax_rate_constant_log.text(0.4, 0.65, text_str3, transform=ax_rate_constant_log.transAxes,
                          fontsize=12, color='black', weight='bold', ha='right', va='top')
ax_rate_constant_log.set_xlabel("1 / Temperature (1 / K)")
ax_rate_constant_log.set_ylabel(r"ln [Rate constant (s$^{-1}$)] ")
fig_rate_constant_log.tight_layout()
fig_rate_constant_log.savefig("{}/arrhenius.png".format(folder), dpi=600)

# Plot mobility vs temperature
# fig_mobility, ax_mobility = plt.subplots(figsize=(6, 6))
# ax_mobility.plot(temperature, mobility, 'bx-')
# ax_mobility.plot(temperature3, mobility3, 'rx-')
# ax_mobility.set_xlabel("Temperature / K")
# ax_mobility.set_ylabel(r"Mobility / cm$^2$ V$^{-1}$ s$^{-1}$")
# fig_mobility.tight_layout()
# fig_mobility.savefig("{}/mobility.png".format(folder), dpi=600)

# Plot log(mobility) vs temperature
# fig_mobility_log, ax_mobility_log = plt.subplots(figsize=(6, 6))
# ax_mobility_log.plot(temperature, np.log10(mobility), 'bx-')
# ax_mobility_log.plot(temperature3, np.log10(mobility3), 'rx-')
# ax_mobility_log.set_xlabel("Temperature (K)")
# ax_mobility_log.set_ylabel(r"Log [mobility (cm$^2$/Vs)]")
# fig_mobility_log.tight_layout()
# fig_mobility_log.savefig("{}/mobility_log.png".format(folder), dpi=600)

if __name__ == "__main__":
    print('Finished.')
    plt.show()
