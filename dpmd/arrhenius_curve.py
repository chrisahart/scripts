import numpy as np
import matplotlib.pyplot as plt
from general import parameters as param


def calc_rate(vn, kb_t, k_el, energy):
    """ Calculate electron transfer rate constant [Spencer 2016]. """
    return vn * k_el * np.exp(-energy/kb_t)


def calc_diffusion(multiplicity, r, k):
    """ Calculate diffusion coefficient. """
    return (r**2 * multiplicity * k) / 2


def calc_mobility(diffusion, kb_t):
    """ Calculate mobility. """
    return (1.6e-19 * diffusion) / kb_t


folder_rutile = '/Volumes/Samsung/Data/Postdoc2/Data/Work/calculations/tio2/rutile/deepmd/rutile/336/md-cell-opt/deepmd-md/hse-22-ts-md2'

# All hops
# rate_constant_rutile = np.array([2466806150140.934, 1922901242452.6802, 1193889820010.756, 563188642737.2736, 95850623481.80663])
# mobility_rutile = np.array([0.04392254697643786, 0.03992144978775102, 0.03061677378314711, 0.019304924800510884, 0.004902278090620809])
# temperature_rutile = np.array([300, 250, 200, 150, 100])

# Hops greater than 4 A removed
rate_constant_rutile = np.array([2275753014439.75, 1868394642018.7478, 1188002204103.6323, 553257682674.263, 95850623481.80663])
mobility_rutile = np.array([0.04052076352971199, 0.03878983550393645, 0.03046578849009147, 0.018964512330043627, 0.004902278090620809])
temperature_rutile = np.array([300, 250, 200, 150, 100])

rate_constant_rutile = rate_constant_rutile[:-1]
mobility_rutile = mobility_rutile[:-1]
temperature_rutile = temperature_rutile[:-1]

folder3 = '/Volumes/Samsung/Data/Postdoc2/Data/Work/calculations/tio2/anatase/deepmd/anatase/441/deepmd/hse-19-ts-md-9500-9900-removed'

# folder_spin = 'single-fit-pop-dpa3-nlayers-6-official-v3.1.0-dev-polaron-loss-mae-pref-1-pref_pop-1000-1'
# rate_constant3 = np.array([1422749145062.692, 833346616448.9613, 403902936462.9768, 94455838920.67673])
# mobility3 = np.array([0.010993409695600897, 0.007573700430471499, 0.004588796525929469, 0.001414400073473747])
# temperature3 = np.array([600, 500, 400, 300])

# folder_spin = 'single-fit-pop-dpa3-nlayers-6-official-v3.1.0-dev-polaron-loss-mae-pref-1-pref_pop-1000-1-rcut-4.5'
# rate_constant3 = np.array([1531240835791.8442, 931905377310.2635, 455628998959.0827, 112557518910.90729])
# mobility3 = np.array([0.011672085531292287, 0.008508652656509129, 0.005206972090464953, 0.0016916012658779792])
# temperature3 = np.array([600, 500, 400, 300])

# folder_spin = 'single-fit-pop-dpa3-nlayers-6-official-v3.1.0-dev-polaron-loss-mae-pref-1-pref_pop-1000-1-rcut-4.5-twostep-lr-1e-5-1e-8'
# rate_constant3 = np.array([1503624403816.1074, 1214594101485.2908, 896123298214.007, 682763115402.9434, 425681041881.4069, 227440979090.41116, 104679543405.30695])
# mobility3 = np.array([0.0114328287435331,  0.010114301896832892, 0.008219924281108409, 0.0068804432778992106, 0.00483449432398798, 0.002923648595162479, 0.0015799973199492336])
# temperature3 = np.array([600, 550, 500, 450, 400, 350, 300])
rate_constant3 = np.array([1503624403816.1074, 1214594101485.2908, 896123298214.007, 682763115402.9434, 425681041881.4069, 227440979090.41116])
mobility3 = np.array([0.0114328287435331,  0.010114301896832892, 0.008219924281108409, 0.0068804432778992106, 0.00483449432398798, 0.002923648595162479])
temperature3 = np.array([600, 550, 500, 450, 400, 350])

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
inv_T = 1 / temperature_rutile
ln_k = np.log(rate_constant_rutile)
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
rutile_freq = rate_constant_rutile / np.exp(-(Ea/1000)/(k_au * temperature_rutile))
anatase_freq = rate_constant3 / np.exp(-(Ea3/1000)/(k_au * temperature3))
print('rutile_freq', rutile_freq/1e13, np.mean(rutile_freq)/1e13)
print('anatase_freq', anatase_freq/1e13, np.mean(anatase_freq)/1e13)

fig_rate_constant_log, ax_rate_constant_log = plt.subplots(figsize=(6, 6))
ax_rate_constant_log.plot(1/temperature_rutile, np.log(rate_constant_rutile), 'bx')
ax_rate_constant_log.plot(1/temperature3, np.log(rate_constant3), 'rx')
ax_rate_constant_log.plot(inv_T, ln_k_fit, 'b-', alpha=0.7, label=f'Fit: slope = {coeffs[0]:.2f}')
ax_rate_constant_log.plot(inv_T3, ln_k_fit3, 'r-', alpha=0.7, label=f'Fit: slope = {coeffs3[0]:.2f}')
ax_rate_constant_log.text(0.90, 0.90, text_str, transform=ax_rate_constant_log.transAxes,
                          fontsize=12, color='black', weight='bold', ha='right', va='top')
ax_rate_constant_log.text(0.4, 0.65, text_str3, transform=ax_rate_constant_log.transAxes,
                          fontsize=12, color='black', weight='bold', ha='right', va='top')
ax_rate_constant_log.set_xlabel("1 / Temperature (1 / K)")
ax_rate_constant_log.set_ylabel(r"Log [Rate constant (s$^{-1}$)] ")
fig_rate_constant_log.tight_layout()
fig_rate_constant_log.savefig("{}/arrhenius.png".format(folder_rutile), dpi=600)

# Plot mobility vs temperature
# fig_mobility, ax_mobility = plt.subplots(figsize=(6, 6))
# ax_mobility.plot(temperature, mobility, 'bx-')
# ax_mobility.plot(temperature3, mobility3, 'rx-')
# ax_mobility.set_xlabel("Temperature / K")
# ax_mobility.set_ylabel(r"Mobility / cm$^2$ V$^{-1}$ s$^{-1}$")
# fig_mobility.tight_layout()
# fig_mobility.savefig("{}/mobility.png".format(folder), dpi=600)

# Plot log(mobility) vs temperature
fig_mobility_log, ax_mobility_log = plt.subplots(figsize=(6, 6))
ax_mobility_log.plot(temperature_rutile, np.log10(mobility_rutile), 'bx-')
ax_mobility_log.plot(temperature3, np.log10(mobility3), 'rx-')
ax_mobility_log.set_xlabel("Temperature (K)")
ax_mobility_log.set_ylabel(r"Log [mobility (cm$^2$/Vs)]")
fig_mobility_log.tight_layout()
fig_mobility_log.savefig("{}/mobility_log.png".format(folder_rutile), dpi=600)

ev_to_joules = 1.60218e-19

# Rutile analytical
# Ea = 35
temperature_array_rutile = np.linspace(10, 1000, 1000)
kb_t = 1.38e-23 * temperature_array_rutile  # KbT in SI units
rutile_rate_analytical = calc_rate(np.mean(rutile_freq), kb_t, 1.0, Ea/1e3*ev_to_joules)
rutile_diff_analytical = calc_diffusion(1.0, 2.96e-10, rutile_rate_analytical)
rutile_mobility_analytical = calc_mobility(rutile_diff_analytical, kb_t)

# Anatase analytical
temperature_array_anatase = np.linspace(10, 1000, 1000)
kb_t = 1.38e-23 * temperature_array_anatase  # KbT in SI units
anatase_rate_analytical = calc_rate(np.mean(anatase_freq), kb_t, 1.0, Ea3/1e3*ev_to_joules)
anatase_diff_analytical = calc_diffusion(1.0, 2.81388e-10, anatase_rate_analytical)
anatase_mobility_analytical = calc_mobility(anatase_diff_analytical, kb_t)

# literature
literature_rutile_deskins = 5.2E-02
literature_rutile_morita = 5.3E-02
literature_rutile_dai = 7.6e-2
# literature_rutile_birschitzky = 1.3
literature_rutile_austin = 1e-2
# literature_rutile_austin = 1
literature_rutile_yagi = 1e-2
literature_rutile_tamaki = 1e-1

# Plot subplot
fig2, axes2 = plt.subplots(1, 2, figsize=(12, 6))

axes2[0].plot(1/temperature_rutile, np.log(rate_constant_rutile), 'b*', label=r'e$^-$ This work', fillstyle='full')
axes2[0].plot(1/temperature3, np.log(rate_constant3), 'r*',  label=r'h$^+$ This work', fillstyle='full')
axes2[0].plot(inv_T, ln_k_fit, 'b-', alpha=0.7, label=f'Fit: slope = {coeffs[0]:.2f}')
axes2[0].plot(inv_T3, ln_k_fit3, 'r-', alpha=0.7, label=f'Fit: slope = {coeffs3[0]:.2f}')
axes2[0].text(0.90, 0.90, text_str, transform=axes2[0].transAxes,
                          fontsize=12, color='black', weight='bold', ha='right', va='top')
axes2[0].text(0.4, 0.65, text_str3, transform=axes2[0].transAxes,
                          fontsize=12, color='black', weight='bold', ha='right', va='top')
# axes2[0].plot(1/temperature_array_rutile, np.log(rutile_rate_analytical), 'k-')
# axes2[0].plot(1/temperature3, np.log(anatase_rate_analytical), 'k-')
axes2[0].set_xlabel("1 / Temperature (1 / K)")
axes2[0].set_ylabel(r"Log [Rate constant (s$^{-1}$)] ")

axes2[1].plot(temperature_rutile, np.log10(mobility_rutile), 'b*', label=r'e$^-$ This work', fillstyle='full')
axes2[1].plot(temperature3, np.log10(mobility3), 'r*',  label=r'h$^+$ This work', fillstyle='full')
axes2[1].plot(temperature_array_rutile, np.log10(rutile_mobility_analytical*100**2), 'b-')
axes2[1].plot(temperature_array_anatase, np.log10(anatase_mobility_analytical*100**2), 'r-')

axes2[1].plot(300, np.log10(5.2E-02), 'bo', label=r'e$^-$ Deskins', fillstyle='full')
axes2[1].plot(300, np.log10(5.3E-02), 'bs', label=r'e$^-$ Morita', fillstyle='full')
axes2[1].plot(300, np.log10(7.6e-2), 'bv', label=r'e$^-$ Dai', fillstyle='full')
# axes2[1].plot(300, np.log10(1.3), 'bD', label=r'e$^-$ Birschitzky', fillstyle='full')

axes2[1].plot(300, np.log10(1.6e-3), 'ro', label=r'h$^+$ Deskins', fillstyle='full')
axes2[1].plot(300, np.log10(1.3e-3), 'rs', label=r'h$^+$ Carey', fillstyle='full')

# axes2[1].plot(300, np.log10(1), 'b^', label=r'e$^-$ Austin, Hendry', fillstyle='none')
axes2[1].plot(300, np.log10(1E-02), 'bX', label=r'e$^-$ Yagi', fillstyle='none')
axes2[1].plot(300, np.log10(1E-01), 'bP', label=r'e$^-$ Tamaki', fillstyle='none')

# axes2[1].legend(frameon=True)
# axes2[1].set_ylim([-4.4, 0.25])

# axes2[1].legend(frameon=True, loc='upper right')
# axes2[1].set_ylim([-3.3, 0.5])

axes2[1].legend(frameon=True, loc='lower right')
axes2[1].set_ylim([-5.2, 0.25])
axes2[1].set_ylim([-3.2, -0.9])

axes2[1].set_xlim([100-20, 600+20])
axes2[1].set_xlabel("Temperature (K)")
axes2[1].set_ylabel(r"Log [mobility (cm$^2$/Vs)]")

fig2.tight_layout()
fig2.savefig("{}/subplot.png".format(folder_rutile), dpi=600)

if __name__ == "__main__":
    print('Finished.')
    plt.show()
