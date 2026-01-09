import numpy as np
import matplotlib.pyplot as plt
from general import parameters as param
from scipy.optimize import curve_fit


def pes_fit(x, a, b, c):
    """ PES function """
    return a * x**2 + b * x + c
    # return a * np.exp(b * x) + c


def calc_energy_marcus(lambda_tot, v_ab):
    """ Calculate activation energy using Marcus method. Non-adiabatic limit v_ab = 0 """
    return lambda_tot / 4


def calc_energy_walsh(lambda_tot, v_ab):
    """ Calculate activation energy using Walsh method. Weak coupling """
    return lambda_tot / 4 - v_ab


def calc_energy_spencer(lambda_tot, v_ab):
    """ Calculate activation energy using Jochen preferred method [Spencer 2016]. Strong coupling """
    return lambda_tot / 4 - (v_ab - 1/lambda_tot * v_ab**2)


replica = np.linspace(start=1, stop=9, num=9)
folder_save = '/Users/chris/Library/Mobile Documents/com~apple~CloudDocs/Work/Postdoc2/images/other'

# Rutile 336
# folder = '/Volumes/Samsung/Data/Postdoc2/Data/Work/calculations/tio2/rutile/archer/rutile/cell-336/cdft-yungyu/cdft/atoms-78-79/hse-22-dz/neb/opt-loose'
# folder_save = folder
# datapoints_fit = 3
# energy_linear = np.array(
#     [-9762.753577167, -9762.752918532, -9762.751333435, -9762.749775396, -9762.749141520, -9762.749667441,
#      -9762.751047882, -9762.752538172, -9762.753258105])
# distances_linear = np.array([0, 0.139290, 0.139290, 0.139290, 0.139290, 0.139290, 0.139290, 0.139290, 0.139290])
# energy_neb_au = np.array(
#     [-9762.753609068, -9762.753334480, -9762.752622105, -9762.751595710, -9762.750578583, -9762.751644531,
#      -9762.752508052, -9762.753055152, -9762.753280049])
# distances_neb = np.array([0, 0.153036, 0.158806, 0.160883, 0.159180, 0.160263, 0.165753, 0.161765, 0.155795])

# Anatase 442 nn-1
# folder = '/Volumes/Samsung/Data/Postdoc2/Data/Work/calculations/tio2/anatase/archer/anatase/cell-442/cdft-yungyu/cdft/atoms-282-305/neb/opt-loose'
# folder_save = folder
# datapoints_fit = 4
# energy_linear = np.array(
#     [-11572.502128082, -11572.501140159, -11572.498538581, -11572.494347994, -11572.490913239, -11572.494324440,
#      -11572.498551182, -11572.501235707, -11572.502215435])
# distances_linear = np.array([0, 0.138637, 0.138637, 0.138637, 0.138637, 0.138637, 0.138637, 0.138637, 0.138637])
# energy_neb_au = np.array(
#     [-11572.502122635, -11572.501932193, -11572.501205352, -11572.498275480, -11572.492368385, -11572.498609926,
#      -11572.501212763, -11572.502021656, -11572.502231501])
# distances_neb = np.array([0, 0.156266, 0.155167, 0.156760, 0.159355, 0.201935, 0.180115, 0.170900, 0.159059])

# Anatase 442 nn-2
folder = '/Volumes/Samsung/Data/Postdoc2/Data/Work/calculations/tio2/anatase/archer/anatase/cell-442/cdft-yungyu/cdft/atoms-282-340/neb/opt-loose'
folder_save = folder
datapoints_fit = 4
energy_linear = np.array(
    [-11572.502138208, -11572.501274767, -11572.499138112, -11572.496351278, -11572.494680470, -11572.496323467,
     -11572.499091583, -11572.501194533, -11572.502155607])
distances_linear = np.array([0, 0.134995, 0.134995, 0.134995, 0.134995, 0.134995, 0.134995, 0.134995, 0.134995])
energy_neb_au = np.array(
    [-11572.502145486, -11572.501993647, -11572.501428642, -11572.499120237, -11572.496346942, -11572.499376534,
     -11572.501499010, -11572.502003350, -11572.502163805])
distances_neb = np.array([0, 0.172638, 0.175487, 0.181760, 0.180325, 0.184482, 0.178338, 0.172304, 0.172153])

# -----------------

energy_linear_ev = (energy_linear-np.min(energy_linear)) * param.hartree_to_ev
cumulative_distances_linear = np.cumsum(distances_linear)
cumulative_distances_linear_norm = cumulative_distances_linear / cumulative_distances_linear[-1]

energy_neb_ev = (energy_neb_au-np.min(energy_neb_au)) * param.hartree_to_ev
cumulative_distances_neb = np.cumsum(distances_neb)
cumulative_distances_neb_norm = cumulative_distances_neb / cumulative_distances_neb[-1]

# coefficients_linear_1 = np.polyfit(cumulative_distances_linear_norm[:datapoints_fit], energy_linear_ev[:datapoints_fit], 2)
# coefficients_linear_2 = np.polyfit(cumulative_distances_linear_norm[-datapoints_fit:], energy_linear_ev[-datapoints_fit:], 2)
# parabola_linear_1 = np.poly1d(coefficients_linear_1)
# parabola_linear_2 = np.poly1d(coefficients_linear_2)
popt_linear_1, _ = curve_fit(pes_fit, cumulative_distances_linear_norm[:datapoints_fit], energy_linear_ev[:datapoints_fit])
popt_linear_2, _ = curve_fit(pes_fit, cumulative_distances_linear_norm[-datapoints_fit:], energy_linear_ev[-datapoints_fit:])
parabola_linear_1 = lambda x: pes_fit(x, *popt_linear_1)
parabola_linear_2 = lambda x: pes_fit(x, *popt_linear_2)

coupling_linear_1 = parabola_linear_1(cumulative_distances_linear_norm[4])-energy_linear_ev[4]
coupling_linear_2 = parabola_linear_2(cumulative_distances_linear_norm[4])-energy_linear_ev[4]
coupling_linear_12 = np.mean(np.array([coupling_linear_1, coupling_linear_2]))
reorg_linear_1 = parabola_linear_2(cumulative_distances_linear_norm[0])
reorg_linear_2 = parabola_linear_1(cumulative_distances_linear_norm[-1])
reorg_linear_12 = np.mean(np.array([reorg_linear_1, reorg_linear_2]))
linear_activation_marcus = calc_energy_marcus(reorg_linear_12, coupling_linear_12)
linear_activation_walsh = calc_energy_walsh(reorg_linear_12, coupling_linear_12)
linear_activation_jacob = calc_energy_spencer(reorg_linear_12, coupling_linear_12)
print('\nLinear Reorganisation energy:', reorg_linear_1*1e3, reorg_linear_2*1e3, reorg_linear_12*1e3)
print('Linear Electronic coupling:', coupling_linear_1*1e3, coupling_linear_2*1e3, coupling_linear_12*1e3)
print('Linear adiabatic barrier height:', energy_linear_ev[4]*1e3)
print('Linear diabatic barrier height:', parabola_linear_1(cumulative_distances_linear_norm[4])*1e3)
print('Linear Marcus activation barrier:', linear_activation_marcus*1e3)
print('Linear Walsh activation barrier:', linear_activation_walsh*1e3)
print('Linear Jacob activation barrier:', linear_activation_jacob*1e3)

# coefficients_neb_1 = np.polyfit(cumulative_distances_neb_norm[:datapoints_fit], energy_neb_ev[:datapoints_fit], 2)
# coefficients_neb_2 = np.polyfit(cumulative_distances_neb_norm[-datapoints_fit:], energy_neb_ev[-datapoints_fit:], 2)
# parabola_neb_1 = np.poly1d(coefficients_neb_1)
# parabola_neb_2 = np.poly1d(coefficients_neb_2)
popt_neb_1, _ = curve_fit(pes_fit, cumulative_distances_neb_norm[:datapoints_fit], energy_neb_ev[:datapoints_fit])
popt_neb_2, _ = curve_fit(pes_fit, cumulative_distances_neb_norm[-datapoints_fit:], energy_neb_ev[-datapoints_fit:])
parabola_neb_1 = lambda x: pes_fit(x, *popt_neb_1)
parabola_neb_2 = lambda x: pes_fit(x, *popt_neb_2)

coupling_neb_1 = parabola_neb_1(cumulative_distances_neb_norm[4])-energy_neb_ev[4]
coupling_neb_2 = parabola_neb_2(cumulative_distances_neb_norm[4])-energy_neb_ev[4]
coupling_neb_12 = np.mean(np.array([coupling_neb_1, coupling_neb_2]))
reorg_neb_1 = parabola_neb_2(cumulative_distances_neb_norm[0])
reorg_neb_2 = parabola_neb_1(cumulative_distances_neb_norm[-1])
reorg_neb_12 = np.mean(np.array([reorg_neb_1, reorg_neb_2]))
neb_activation_marcus = calc_energy_marcus(reorg_neb_12, coupling_neb_12)
neb_activation_walsh = calc_energy_walsh(reorg_neb_12, coupling_neb_12)
neb_activation_jacob = calc_energy_spencer(reorg_neb_12, coupling_neb_12)
print('\nneb Reorganisation energy:', reorg_neb_1*1e3, reorg_neb_2*1e3, reorg_neb_12*1e3)
print('neb Electronic coupling:', coupling_neb_1*1e3, coupling_neb_2*1e3, coupling_neb_12*1e3)
print('neb adiabatic barrier height:', energy_neb_ev[4]*1e3)
print('neb diabatic barrier height:', parabola_neb_1(cumulative_distances_neb_norm[4])*1e3)
print('neb Marcus activation barrier:', neb_activation_marcus*1e3)
print('neb Walsh activation barrier:', neb_activation_walsh*1e3)
print('neb Jacob activation barrier:', neb_activation_jacob*1e3)

fig_spin_valid, axes_spin_valid = plt.subplots()
# axes_spin_valid.plot(5, (-9762.7491174583-np.min(energy_linear))*param.hartree_to_ev*1e3, 'bx')
# axes_spin_valid.plot(5, (-9762.7492265083-np.min(energy_linear))*param.hartree_to_ev*1e3, 'yx')
axes_spin_valid.plot(cumulative_distances_linear_norm, energy_linear_ev, 'rx-', label='Linear')
axes_spin_valid.plot(cumulative_distances_linear_norm, parabola_linear_1(cumulative_distances_linear_norm), 'r-')
axes_spin_valid.plot(cumulative_distances_linear_norm, parabola_linear_2(cumulative_distances_linear_norm), 'r-')

axes_spin_valid.plot(cumulative_distances_neb_norm, energy_neb_ev, 'gx-', label='NEB')
axes_spin_valid.plot(cumulative_distances_neb_norm, parabola_neb_1(cumulative_distances_neb_norm), 'g-')
axes_spin_valid.plot(cumulative_distances_neb_norm, parabola_neb_2(cumulative_distances_neb_norm), 'g-')

axes_spin_valid.set_xlabel("Fractional reaction coordinate")
axes_spin_valid.set_ylabel("Relative energy / eV")
axes_spin_valid.legend(frameon=False)
plt.tight_layout()
fig_spin_valid.savefig('{}/pes.png'.format(folder_save), dpi=param.save_dpi)

if __name__ == "__main__":
    print('Finished.')
    plt.show()
