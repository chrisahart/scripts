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
# temp = 600  # K
multiplicity = 1  # Site multiplicity
# vn = 1.85e13  # Effective nuclear frequency Fe-O
# vn_ev = 98.8/1e3 # Dai et al from phonon spectra
# vn_ev = 0.10  # 0.10 eV to s^-1 Deskins Dupuis TiO2 rutile (optic-mode phonon frequencies)
# vn_ev = 0.11  # 0.11 eV to s^-1 Deskins Dupuis TiO2 anatase (optic-mode phonon frequencies)
# vn_s = vn_ev * ev_to_joules / planck
# vn = vn_s
# print('Effective nuclear frequency e13 s-1', vn_s/1e13)
# print('Effective nuclear frequency fs', 1/vn_s * 1e15)
# vn = 1.13e13  # rutile NNP-MD
vn = 2.66e13  # anatase NNP-MD
# vn = 2.42e13  # 0.10 eV to s^-1 Deskins Dupuis TiO2 rutile (optic-mode phonon frequencies)
# Effective nuclear frequency e13 s-1 2.4165610859728504
# Effective nuclear frequency fs 41.38111822641651
# vn = 2.66e13  # 0.11 eV to s^-1 Deskins Dupuis TiO2 anatase (optic-mode phonon frequencies)
kb_t_au = 8.617333262145E-5 * temp  # KbT in eV
kb_t = 1.38e-23 * temp  # KbT in SI units
r_hop = np.array([2.81388])
folder = '/Users/chris/Documents/Storage/images'

energy_array = np.linspace(start=0.01, stop=150, num=int(1e3)) / 1e3
mobility_spencer = np.zeros_like(energy_array)
for i in range(len(energy_array)):
    rate_spencer = calc_rate(vn, kb_t_au, 1.0, energy_array[i])
    diffusion_spencer = calc_diffusion(multiplicity, r_hop, rate_spencer)
    mobility_spencer[i] = calc_mobility(diffusion_spencer, kb_t_au)

fig_energy_mobility, ax_energy_mobility = plt.subplots()
ax_energy_mobility.plot(energy_array*1e3, mobility_spencer, 'k-')

ax_energy_mobility.plot(136, 1.7e-3, 'r*', label=r'e$^-$ This work', fillstyle='full', markersize=10)
ax_energy_mobility.plot(136, 2.11e-3, 'r*', label=r'e$^-$ From DP ∆A$^{‡}$', fillstyle='none', markersize=10)
ax_energy_mobility.plot(143, 1.6e-3, 'ro', label=r'h$^+$ Deskins', fillstyle='full')
ax_energy_mobility.plot(133, 2.38E-03, 'rs', label=r'h$^+$ Carey', fillstyle='full')

ax_energy_mobility.set_xlabel('Activation energy / meV')
ax_energy_mobility.set_ylabel('Mobility')
ax_energy_mobility.legend(frameon=True, loc='upper right')
fig_energy_mobility.tight_layout()
fig_energy_mobility.savefig("{}/anatase_energy_mobility.png".format(folder), dpi=600)

fig_energy_mobility_log, ax_energy_mobility_log = plt.subplots()
ax_energy_mobility_log.plot(energy_array*1e3, mobility_spencer, 'k-')

ax_energy_mobility_log.plot(136, 1.7e-3, 'r*', label=r'e$^-$ This work', fillstyle='full', markersize=10)
ax_energy_mobility_log.plot(136, 2.11e-3, 'r*', label=r'e$^-$ From DP ∆A$^{‡}$', fillstyle='none', markersize=10)
ax_energy_mobility_log.plot(143, 1.6e-3, 'ro', label=r'h$^+$ Deskins', fillstyle='full')
ax_energy_mobility_log.plot(133, 2.38E-03, 'rs', label=r'h$^+$ Carey', fillstyle='full')

ax_energy_mobility_log.set_xlabel('Activation energy / meV')
ax_energy_mobility_log.set_ylabel('Mobility')
ax_energy_mobility_log.set_yscale('log')
ax_energy_mobility_log.legend(frameon=True, loc='upper right')
fig_energy_mobility_log.tight_layout()
fig_energy_mobility_log.savefig("{}/anatase_energy_mobility_log.png".format(folder), dpi=600)

if __name__ == "__main__":
    print('Finished.')
    plt.show()