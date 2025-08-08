import numpy as np
from general import parameters as param

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
    """ Calculate non-adiabatic activation energy using Jochen preferred method [Spencer 2016]. """
    return lambda_tot / 4


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


# Parameters
temp = 300  # K
multiplicity = 2  # Site multiplicity
# vn = 1.85e13  # Effective nuclear frequency Fe-O
vn = 2.4e13  # 0.10 eV to s^-1 Deskins Dupuis TiO2 rutile (optic-mode phonon frequencies)

# Constants
kb_t_au = 8.617333262145E-5 * temp  # KbT in eV
kb_t = 1.38e-23 * temp  # KbT in SI units
planck = 6.63e-34  # Planck constant in SI units
planck_au = 2 * np.pi  # Planck constant in SI units
angstrom_to_cm = 1e-8
ev_to_joules = 1.60218e-19

# Hematite Hole Boltzmann factor at 300 K
# r_hop = np.array([2.97])
# coupling = np.array([147]) / 1e3
# reorg = np.array([752]) / 1e3

# TiO2 336 22% HFX 1st nearest neighbour
# r_hop = np.array([2.96])
# coupling = np.array([14.446622491241]) * param.hartree_to_ev / 1e3
# energy_cdft_ts = np.array([-9762.745898065810252])
# energy_dft_gs = np.array([-9762.753627197094829])
# reorg = (energy_cdft_ts - energy_dft_gs) * param.hartree_to_ev * 4

# TiO2 336 23% HFX 1st nearest neighbour  ***** energy does not make sense
# r_hop = np.array([2.96])
# coupling = np.array([11.448954242300]) * param.hartree_to_ev / 1e3
# energy_cdft_ts = np.array([-9761.974060083877703])
# energy_dft_gs = np.array([-9762.494274678130751])
# reorg = (energy_cdft_ts - energy_dft_gs) * param.hartree_to_ev * 4

# TiO2 336 24% HFX 1st nearest neighbour ***** energy does not make sense
# r_hop = np.array([2.96])
# coupling = np.array([11.449932666881]) * param.hartree_to_ev / 1e3
# energy_cdft_ts = np.array([-9761.973854499945446])
# energy_dft_gs = np.array([-9762.237305633563665])
# reorg = (energy_cdft_ts - energy_dft_gs) * param.hartree_to_ev * 4

# TiO2 336 25% HFX 1st nearest neighbour
# r_hop = np.array([2.96])
# coupling = np.array([11.587773387758]) * param.hartree_to_ev / 1e3
# energy_cdft_ts = np.array([-9761.973670271818264])
# energy_dft_gs = np.array([-9761.98275975738943])
# reorg = (energy_cdft_ts - energy_dft_gs) * param.hartree_to_ev * 4

# TiO2 336 25% HFX 1st nearest neighbour
r_hop = np.array([2.96])
coupling = np.array([14.408699460995]) * param.hartree_to_ev / 1e3
energy_cdft_ts = np.array([-13017.114620704207482])
energy_dft_gs = np.array([-13017.120554204866494])
# energy_cdft_ts = np.array([-13017.114949657321631])
# energy_dft_gs = np.array([-13017.12058427568445])
reorg = (energy_cdft_ts - energy_dft_gs) * param.hartree_to_ev * 4

# TiO2 446 25% HFX 1st nearest neighbour ***** energy does not make sense
# r_hop = np.array([2.96])
# coupling = np.array([11.053469188353]) * param.hartree_to_ev / 1e3
# energy_cdft_ts = np.array([-17354.882408187037072])
# energy_dft_gs = np.array([-17354.898964491683728])
# reorg = (energy_cdft_ts - energy_dft_gs) * param.hartree_to_ev * 4

# TiO2 336 25% HFX 2nd nearest neighbour
# r_hop = np.array([3.57])
# coupling = np.array([1.687770049697]) * param.hartree_to_ev / 1e3
# energy_cdft_ts = np.array([-9762.740718427978209])
# energy_dft_gs = np.array([-9762.753294884205388])
# reorg = (energy_cdft_ts - energy_dft_gs) * param.hartree_to_ev * 4

# TiO2 336 25% HFX 3rd nearest neighbour
# r_hop = np.array([4.59])
# coupling = np.array([0.225872628396]) * param.hartree_to_ev / 1e3
# energy_cdft_ts = np.array([-9762.74181232870068])
# energy_dft_gs = np.array([-9762.753847150284855 ])
# reorg = (energy_cdft_ts - energy_dft_gs) * param.hartree_to_ev * 4

for i in range(0, np.shape(coupling)[0]):

    adiabaticity_parameter = calc_adiabaticity(vn, kb_t, planck, reorg[i]*ev_to_joules, coupling[i]*ev_to_joules)
    lz_probability = calc_probability(vn, kb_t, planck, reorg[i]*ev_to_joules, coupling[i]*ev_to_joules)
    transmission_coefficient = calc_transmission(vn, kb_t, planck, reorg[i]*ev_to_joules, coupling[i]*ev_to_joules)

    energy_rosso = calc_energy_rosso(reorg[i], coupling[i])
    rate_rosso = calc_rate(vn, kb_t_au, 1, energy_rosso)
    diffusion_rosso = calc_diffusion(multiplicity, r_hop, rate_rosso)
    mobility_rosso = calc_mobility(diffusion_rosso, kb_t_au)

    print("\nRosso Activation energy (delta G*): {0:.2} eV".format(energy_rosso))
    print("Rosso Activation energy (delta G*): {} meV".format(energy_rosso*1e3))
    print("Rosso Electron transfer rate constant (k_et): {0:.2E} s-1".format(rate_rosso))
    print("Rosso Electron transfer rate constant (i k_et): {0:.2E} s-1".format(multiplicity * rate_rosso))
    print("Rosso Mobility: {0:.2E} cm2/V".format(float(mobility_rosso)))

    energy_spencer = calc_energy_spencer(reorg[i], coupling[i])
    rate_spencer_ad = calc_rate(vn, kb_t_au, 1, energy_spencer)
    diffusion_spencer_ad = calc_diffusion(multiplicity, r_hop[i], rate_spencer_ad)
    mobility_spencer_ad = calc_mobility(diffusion_spencer_ad, kb_t_au)

    energy_na_spencer = calc_energy_na_spencer(reorg[i])
    rate_spencer_na = calc_rate_ad(vn, kb_t_au, coupling[i], reorg[i], energy_na_spencer)
    diffusion_spencer_na = calc_diffusion(multiplicity, r_hop[i], rate_spencer_na)
    mobility_spencer_na = calc_mobility(diffusion_spencer_na, kb_t_au)

    rate_spencer = calc_rate(vn, kb_t_au, transmission_coefficient, energy_spencer)
    diffusion_spencer = calc_diffusion(multiplicity, r_hop[i], rate_spencer)
    mobility_spencer = calc_mobility(diffusion_spencer, kb_t_au)

    print("-------------------------------------------")
    print("Distance: {} A".format(r_hop[i]))
    print("Coupling (v_ab): {} meV".format(int(np.round(coupling[i]*1e3))))
    print("Reorganisation energy energy (lambda): {} meV".format(int(np.round(reorg[i]*1e3))))

    print("\nAdiabaticity parameter (2 pi gamma): {}".format(adiabaticity_parameter))
    # print("Landau-Zener transition probability (P_LZ): {0:.2}".format(lz_probability))
    print("Electronic transmission coefficient (k_el): {0:.1}".format(transmission_coefficient))

    print("\nActivation energy (delta G*): {} meV".format(int(np.round(energy_spencer*1e3))))
    print("Electron transfer rate constant (k_et): {0:.1E} s-1".format(rate_spencer))
    print("1/Electron transfer rate constant (1/k_et): {} fs".format((1/rate_spencer)*1e15))
    print("1/Electron transfer rate constant (1/k_et) / 1e6: {} fs".format((1/rate_spencer)*1e15/1e6))
    # print("Mobility: {0:.2} cm2/V".format(mobility_spencer))
    print("Mobility: {0:.1E} cm2/V".format(mobility_spencer))

    # print("\nAdiabatic electron transfer rate constant (k_et): {0:.1E} s-1".format(rate_spencer_ad))
    # print("Adiabatic mobility: {0:.2} cm2/V".format(mobility_spencer_ad))

    # print("\nNon-adiabatic activation energy (sigma A dagger na): {} meV".format(int(np.round(energy_na_spencer*1e3))))
    # print("Non-adiabatic electron transfer rate constant (k_et): {0:.1E} s-1".format(rate_spencer_na))
    # print("Non-adiabatic mobility: {0:.2} cm2/V".format(mobility_spencer_na))