import numpy as np
from general import parameters as param
import scipy.constants as const
from matplotlib import pyplot as plt

""" Calculate hematite mobility. All units are au unless otherwise specified """


def calc_adiabaticity(vn, kb_t, planck, lambda_tot, v_ab):
    """ Calculate adiabaticity parameter (2 pi gamma). """
    return (np.pi**(3/2)*v_ab**2) / (planck*vn*(lambda_tot*kb_t)**(1/2))


def calc_probability(vn, kb_t, planck, lambda_tot, v_ab):
    """ Calculate Landau-Zener transition probability (P_LZ). """
    return 1-np.exp(-calc_adiabaticity(vn, kb_t, planck, lambda_tot, v_ab))


def calc_transmission(vn, kb_t, planck, lambda_tot, v_ab):
    """ Calculate electronic transmission coefficient (k_el). """
    p_lz = calc_probability(vn, kb_t, planck, lambda_tot, v_ab)
    return (2 * p_lz) / (1 + p_lz)


def calc_energy_rosso(lambda_tot, v_ab):
    """ Calculate activation energy using Rosso method [Iordanova 2005]. """
    return -lambda_tot/4 + (lambda_tot**2+4*v_ab**2)**(1/2)/2 - v_ab


def calc_energy_spencer(lambda_tot, v_ab):
    """ Calculate activation energy using Jochen preferred method [Spencer 2016]. """
    return lambda_tot / 4 - (v_ab - 1/lambda_tot * v_ab**2)


def calc_energy_na_spencer(lambda_tot):
    """ Calculate activation energy using Marcus method. Non-adiabatic limit v_ab = 0 """
    """ Calculate non-adiabatic activation energy using Jochen preferred method [Spencer 2016]. """
    return lambda_tot / 4


def calc_energy_weak_coupling(lambda_tot, v_ab):
    """ Calculate activation energy using Walsh method. Weak coupling """
    return lambda_tot / 4 - v_ab


def calc_marcus_factor(lambda_tot, v_ab, kb_t):
    """ Marcus pre-factor. """
    return ((2*np.pi)/6.63e-34) * ((v_ab **2) / np.sqrt(4*np.pi*lambda_tot*kb_t))


def calc_rate(vn, kb_t, k_el, energy):
    """ Calculate electron transfer rate constant [Spencer 2016]. """
    return vn * k_el * np.exp(-energy/kb_t)


def calc_rate_ad(vn, kb_t, v_ab, lambda_tot, energy):
    """ Calculate non-adiabatic electron transfer rate constant [Spencer 2016]. """
    factor = 2*np.pi * v_ab**2 * (4*np.pi*lambda_tot*kb_t)**(-1/2)
    return vn * factor * np.exp(-energy/kb_t)


def calc_diffusion(multiplicity, r, k):
    """ Calculate diffusion coefficient. """
    return ((r*angstrom_to_cm)**2 * multiplicity * k) / 2


def calc_mobility(diffusion, kb_t):
    """ Calculate mobility. """
    return diffusion / kb_t


# Constants
planck = 6.63e-34  # Planck constant in SI units
planck_au = 2 * np.pi  # Planck constant in SI units
angstrom_to_cm = 1e-8
ev_to_joules = 1.60218e-19

# Parameters
temp = 300  # K
multiplicity = 1  # Site multiplicity
vn = 1.13e13  # rutile NNP-MD
kb_t_au = 8.617333262145E-5 * temp  # KbT in eV
kb_t = 1.38e-23 * temp  # KbT in SI units
r_hop = np.array([2.96])
folder = '/Users/chris/Documents/Storage/images'

energy_array = np.linspace(start=0.01, stop=150, num=int(1e3)) / 1e3
mobility_spencer = np.zeros_like(energy_array)
for i in range(len(energy_array)):
    rate_spencer = calc_rate(vn, kb_t_au, 1.0, energy_array[i])
    diffusion_spencer = calc_diffusion(multiplicity, r_hop, rate_spencer)
    mobility_spencer[i] = calc_mobility(diffusion_spencer, kb_t_au)

fig_energy_mobility, ax_energy_mobility = plt.subplots()
ax_energy_mobility.plot(energy_array*1e3, mobility_spencer, 'k-')

ax_energy_mobility.plot(39, 4.4E-02, 'b*', label=r'e$^-$ This work', fillstyle='full', markersize=10)
ax_energy_mobility.plot(39, 4.24E-02, 'b*', label=r'e$^-$ From DP ∆A$^{‡}$', fillstyle='none', markersize=10)
ax_energy_mobility.plot(123, 1.7e-3, 'bo', label=r'e$^-$ Deskins', fillstyle='full')
ax_energy_mobility.plot(56, 2.19E-02, 'bs', label=r'e$^-$ Morita', fillstyle='full')
ax_energy_mobility.plot(13, 1.16E-01, 'bv', label=r'e$^-$ Dai', fillstyle='full')

ax_energy_mobility.set_xlabel('Activation energy / meV')
ax_energy_mobility.set_ylabel('Mobility')
ax_energy_mobility.legend(frameon=True, loc='upper right')
fig_energy_mobility.tight_layout()
fig_energy_mobility.savefig("{}/rutile_energy_mobility.png".format(folder), dpi=600)

fig_energy_mobility_log, ax_energy_mobility_log = plt.subplots()
ax_energy_mobility_log.plot(energy_array*1e3, mobility_spencer, 'k-')

ax_energy_mobility_log.plot(39, 4.4E-02, 'b*', label=r'e$^-$ This work', fillstyle='full', markersize=10)
ax_energy_mobility_log.plot(39, 4.24E-02, 'b*', label=r'e$^-$ From DP ∆A$^{‡}$', fillstyle='none', markersize=10)
ax_energy_mobility_log.plot(123, 1.7e-3, 'bo', label=r'e$^-$ Deskins', fillstyle='full')
ax_energy_mobility_log.plot(56, 2.19E-02, 'bs', label=r'e$^-$ Morita', fillstyle='full')
ax_energy_mobility_log.plot(13, 1.16E-01, 'bv', label=r'e$^-$ Dai', fillstyle='full')

ax_energy_mobility_log.set_xlabel('Activation energy / meV')
ax_energy_mobility_log.set_ylabel('Mobility')
ax_energy_mobility_log.set_yscale('log')
ax_energy_mobility_log.legend(frameon=True, loc='upper right')
fig_energy_mobility_log.tight_layout()
fig_energy_mobility_log.savefig("{}/rutile_energy_mobility_log.png".format(folder), dpi=600)

if __name__ == "__main__":
    print('Finished.')
    plt.show()