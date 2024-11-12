import numpy as np
from matplotlib import pyplot as plt
from ase.io.cube import read_cube_data
from general import parameters as param

kau_to_cu = 1 / 1.5156
cu_to_corehour = 128
corehour_se_rmb = 0.06
corehour_archer_pound = 0.2/128
rmb_to_pound = 1 / 9.28

# print('corehour_archer_pound', corehour_archer_pound)
# print('corehour_archer_pound rmb_to_pound', corehour_archer_pound / rmb_to_pound)
# print('SE cluster is x amount for expensive than ARCHER2:',
#       0.06 / (corehour_archer_pound / rmb_to_pound))

# 1.5156 kAU : 4.21 ARCHER node hour : 1 ARCHER2 node hour : 1 CU
# 9.28 Chinese Yuan = 1.00 Pound Sterling

# Per 6 month period
kau = np.array([16133, 14191, 13140, 13243, 10364, 11826, 12000, 9700, 9700, 15200, 15200])

# Per 12 month period
# kau = np.array([16133+14191, 13140+13243, 10364+11826, 12000, 9700+9700, 15200+15200])

cu = kau * kau_to_cu
core_hours = cu * cu_to_corehour
archer2_cost_pounds = core_hours * corehour_archer_pound
se_cost_rmb = core_hours * corehour_se_rmb
se_cost_pounds = core_hours * corehour_se_rmb * rmb_to_pound

print('kAU', kau)
print('CU', str(np.round(cu)))
print('Core hours', str(np.round(core_hours)))
print('Cost archer £', str(np.round(archer2_cost_pounds)))
print('Cost SE RMB', str(np.round(se_cost_rmb)))
print('Cost SE £', str(np.round(se_cost_pounds)))
#
# print('Average ARCHER usage 6 months:', str(np.round(np.average(cu))))
# print('Average ARCHER usage cost 6 months(£):', str(np.round(np.average(archer2_cost_pounds))))
# print('Average SE usage cost 6 months(£):', str(np.round(np.average(se_cost_pounds))))
# print('Average SE usage cost 6 months(RMB):', str(np.round(np.average(se_cost_rmb))))

print('Average core hours usage :', str(np.round(2*np.average(core_hours))))
print('Average SE usage cost per year (£):', str(np.round(2*np.average(se_cost_pounds))))
print('Average SE usage cost per year (RMB):', str(np.round(2*np.average(se_cost_rmb))))

# 12 month analysis
# print('Core hours usage :', str(np.round(core_hours)))
# print('SE usage :', str(np.round(se_cost_rmb)))

# print('Average core hours usage :', str(np.round(np.average(core_hours))))
# print('Average SE usage cost per year (£):', str(np.round(np.average(se_cost_pounds))))
# print('Average SE usage cost per year (RMB):', str(np.round(np.average(se_cost_rmb))))
