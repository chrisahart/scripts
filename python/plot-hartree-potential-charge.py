import numpy as np
import matplotlib.pyplot as plt
from general import parameters as param

""" Plotting of SMEAGOL output _TRC.agr by filename"""

plot_markers = False
plot_diff = False

# SIESTA Lithium 1D wire
# bulk_charge = '0.Li.leads-RHO_AV.dat'
# bulk_hartree = '0.Li.leads-VH_AV.dat'
# em_hartree = '0.Liwire-VH_AV.dat'
# em_charge = '0.Liwire-RHO_AV.dat'
# labels_em = ['EM']
# labels_bulk = ['Bulk']
# ylim = [-3, 2]
# folder = [
#     '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/iv/li/siesta/single-points/v-0_bottom-0',
#     '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/iv/li/siesta/single-points/v-0_bottom-1']

# CP2K Lithium 1D wire 3e
# bulk_charge = 'bulk-RHO_z.dat'
# bulk_hartree = 'bulk-VH_AV.dat'
# em_hartree = '0V-VH_z.dat'
# em_charge = '0V-RHO_z.dat'
# labels_em = ['EM']
# labels_bulk = ['Bulk']
# ylim = [-3, 2]
# folder = [
#     '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/iv/li/cp2k/transmission/dev-chris/single-points/V-0_HLB-F_z-0-0/',
#     '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/iv/li/cp2k/transmission/dev-chris/single-points/V-0_HLB-1_z-0-0']

# CP2K Lithium 1D wire 1e
# bulk_charge = 'bulk-RHO_z.dat'
# bulk_hartree = 'bulk-VH_AV.dat'
# em_hartree = '0V-VH_z.dat'
# em_charge = '0V-RHO_z.dat'
# labels_em = ['EM']
# labels_bulk = ['Bulk']
# ylim = [-3, 2]
# folder = [
#     '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/iv/li/cp2k/transmission/dev-chris/single-points/V-0_HLB-F_z-0-0-pbe-DZVP-MOLOPT-SR-GTH-q1',
#     '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/iv/li/cp2k/transmission/dev-chris/single-points/V-0_HLB-1_z-0-0-pbe-DZVP-MOLOPT-SR-GTH-q1']

# CP2K Au 1D wire
# bulk_charge = 'bulk-RHO_z.dat'
# bulk_hartree = 'bulk-VH_AV.dat'
# em_hartree = '0V-VH_z.dat'
# em_charge = '0V-RHO_z.dat'
# labels_em = ['EM']
# labels_bulk = ['Bulk']
# ylim = [-3, 2]
# folder = [
#     '/Volumes/Storage/Data/delete/Work/Postdoc/Work/calculations/transport/iv/au/cp2k/dev/single-points/V-0_HLB-0_z-0-0',
#     '/Volumes/Storage/Data/delete/Work/Postdoc/Work/calculations/transport/iv/au/cp2k/dev/single-points/V-0_HLB-1_z-0-0']

# SIESTA Au(s,d) 1D wire
# bulk_charge = '0.Au_wire_leads-RHO_AV.dat'
# bulk_hartree = '0.Au_wire_leads-VH_AV.dat'
# em_hartree = '0.Au-VH_AV.dat'
# em_charge = '0.Au-RHO_AV.dat'
# labels_em = ['EM']
# labels_bulk = ['Bulk']
# ylim = [-3, 2]
# folder = [
#     '/Volumes/Storage/Data/delete/Work/Postdoc/Work/calculations/transport/iv/au/siesta/V-0_HLB-0_z-0-0',
#     '/Volumes/Storage/Data/delete/Work/Postdoc/Work/calculations/transport/iv/au/siesta/V-0_HLB-1_z-0-0']

# SIESTA Au(s,p,d) 1D wire
# bulk_charge = '0.Au_wire_leads-RHO_AV.dat'
# bulk_hartree = '0.Au_wire_leads-VH_AV.dat'
# em_hartree = '0.Au-VH_AV.dat'
# em_charge = '0.Au-RHO_AV.dat'
# labels_em = ['EM']
# labels_bulk = ['Bulk']
# ylim = [-3, 2]
# folder = [
#     '/Volumes/Storage/Data/delete/Work/Postdoc/Work/calculations/transport/iv/au/siesta/V-0_HLB-0_z-0-0-sdp',
#     '/Volumes/Storage/Data/delete/Work/Postdoc/Work/calculations/transport/iv/au/siesta/V-0_HLB-1_z-0-0-sdp']

# CP2K Lithium 1D wire 1e
# bulk_charge = 'bulk-RHO_z.dat'
# bulk_hartree = 'bulk-VH_AV.dat'
# em_hartree = '0V-VH_z_ORIG-0001.dat'
# em_charge = 'RHO_AO-0001.txt'
# labels_em = ['EM']
# labels_bulk = ['Bulk']
# ylim = [-3, 2]
# folder = [
#     '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/au-capacitor-100/struct_5-4/cp2k-smeagol/single-points/old/V-0_HLB-F_z-0-0_hirshfeld-gaussian_ordered']

# SIESTA Au-BDT
# bulk_charge = '0.bulk-RHO_AV.dat'
# bulk_hartree = '0.bulk-VH_AV.dat'
# em_hartree = '0.transport-VH_AV.dat'
# em_charge = '0.transport-RHO_AV.dat'
# labels_em = ['EM']
# labels_bulk = ['Bulk']
# ylim = [-3, 2]
# x = 33.834 - 8.292
# folder = ['/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/siesta-smeagol/transmission/exp/bulk-4-4-100-em-4-4-1_hlb-15.2496']
# folder = ['/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/siesta-smeagol/transmission/opt-cg/bulk-4-4-100-em-4-4-1_hlb-15.412-0-0']
# folder_save = []
# for i in range(len(folder)):
#     folder_save.append('{}/output'.format(folder[i]))
# print('Saving to', folder_save)

# # SIESTA Au-H2-Au
# bulk_charge = '0.lead-RHO_AV.dat'
# bulk_hartree = '0.lead-VH_AV.dat'
# em_hartree = '0.mx-VH_AV.dat'
# em_charge = '0.mx-RHO_AV.dat'
# labels_em = ['EM']
# labels_bulk = ['Bulk']
# ylim = [-3, 2]
# x = 38.72 - 8.36
# # folder = ['/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/AuH2/transport/siesta-smeagol/transmission/HLB-F']
# folder = ['/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/AuH2/transport/siesta-smeagol/transmission/HLB-T-12.19948578']
# folder_save = []
# for i in range(len(folder)):
#     folder_save.append('{}/output'.format(folder[i]))
# print('Saving to', folder_save)

# Au capacitor
# bulk_charge = '0.bulk-RHO_AV.dat'
# bulk_hartree = '0.bulk-VH_AV.dat'
# em_hartree = '0.transport-VH_AV.dat'
# em_charge = '0.transport-RHO_AV.dat'
# labels_em = ['EM']
# labels_bulk = ['Bulk']
# ylim = [-3, 2]
# plot_markers = True
# markers = np.array([2.08400*0, 2.08400*1, 2.08400*2, 2.08400*3, 2.08400*4, 2.08400*5, 21.42200+2.08400*0,
#                     21.42200+2.08400*1, 21.42200+2.08400*2, 21.42200+2.08400*3, 21.42200+2.08400*4, 21.42200+2.08400*5])
# x = 33.926 - 8.336
# folder = ['/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/siesta-smeagol/capacitor/bulk-4-4-100-em-4-4-1_hlb-15.2496_0-0']
# folder = ['/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/siesta-smeagol/capacitor/V-1_bulk-4-4-100-em-4-4-1_hlb-15.2496_0-0']
# folder_save = []
# for i in range(len(folder)):
#     folder_save.append('{}/output'.format(folder[i]))
# print('Saving to', folder_save)

# Au wire
# bulk_charge = '0.Au_wire_leads-RHO_AV.dat'
# bulk_hartree = '0.Au_wire_leads-VH_AV.dat'
# em_hartree = '0.Au-VH_AV.dat'
# em_charge = '0.Au-RHO_AV.dat'
# labels_em = ['EM']
# labels_bulk = ['Bulk']
# ylim = [-3, 2]
# plot_dft = False
# plot_markers = False
# plot_diff = False
# x = 78.4 - 11.2
# folder = [
#     '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-chain/siesta/transmission/V-0_HLB-0_z-0-0']
# folder_save = []
# for i in range(len(folder)):
#     folder_save.append('{}/output'.format(folder[i]))
# print('Saving to', folder_save)

# # Au wire
# bulk_charge = '0.bulk-RHO_AV.dat'
# bulk_hartree = '0.bulk-VH_AV.dat'
# dft_hartree = '0.dft-VH_AV.dat'
# dft_charge = '0.dft-RHO_AV.dat'
# em_hartree = '0.transport-VH_AV.dat'
# em_charge = '0.transport-RHO_AV.dat'
# labels_em = ['EM']
# labels_bulk = ['Bulk']
# ylim = [-3, 2]
# plot_dft = True
# plot_markers = False
# plot_diff = True
# x = 78.4 - 11.2
# # folder = ['/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-chain/siesta/transmission/testing/single-points/V-0_HLB-F_kpoints-1-1-20']
# # folder = ['/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-chain/siesta/transmission/testing/single-points/V-0_HLB-auto_kpoints-1-1-20']
# # folder = ['/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-chain/siesta/transmission/testing/single-points/V-1_HLB-auto_kpoints-1-1-20']
# folder = ['/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-chain/siesta/transmission/testing/single-points/V-0_HLB-auto_kpoints-1-1-20',
#           '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-chain/siesta/transmission/testing/single-points/V-1_HLB-auto_kpoints-1-1-20']
# folder_save = []
# for i in range(len(folder)):
#     folder_save.append('{}/output'.format(folder[i]))
# print('Saving to', folder_save)

# Au capacitor bias dependence (bound states?)
bulk_charge = '0.bulk-RHO_AV.dat'
bulk_hartree = '0.bulk-VH_AV.dat'
# dft_hartree = '0.dft-VH_AV.dat'
# dft_charge = '0.dft-RHO_AV.dat'
dft_hartree = '0.0V-VH_AV.dat'
dft_charge = '0.0V-RHO_AV.dat'
em_hartree = '0.4V-VH_AV.dat'
em_charge = '0.4V-RHO_AV.dat'
labels = ['V=0', 'V=1', 'Bulk']
ylim = [-3, 2]
plot_dft = True
plot_markers = False
plot_diff = False
plot_diff_multiple = True
# diff_label = ['V=4 - V=0 delta=1e-2', 'V=4 - V=0 delta=1e-4', 'V=4 - V=0 delta=0', 'V=4 - V=0 BSCS delta=1e-1', 'V=4 - V=0 BSCS delta=1e-2', 'V=4 - V=0 BSCS delta=1e-4', 'V=4 - V=0 BSCS delta=0']
diff_label = ['V=1 - V=0 EM.WeightRho 0.5', 'V=4 - V=0 EM.WeightRho 0.5', 'V=4 - V=0 weighted double contour', 'V=4 - V=0 weighted double contour 1']
# diff_label = ['V=4 - V=0 WeightRho 0', 'V=4 - V=0 WeightRho 0.5', 'V=4 - V=0 WeightRho 1', 'V=4 - V=0 double contour 1', 'V=4 - V=0 double contour 2']
# diff_label = ['V=1 - V=0', 'V=4 - V=0', ]
# diff_label = ['V=1 - V=0', 'V=2 - V=0', 'V=3 - V=0', 'V=4 - V=0', 'V=5 - V=0', 'V=10 - V=0']
x = 78.4 - 11.2
folder = ['/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/au-capacitor/siesta1-smeagol/kpoints-2-2-20_V-1_WeightRho-0.5',
          '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/au-capacitor/siesta1-smeagol/kpoints-2-2-20_V-4_WeightRho-0.5']
# folder = ['/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/au-capacitor/siesta1-smeagol/kpoints-2-2-20_V-1_WeightRho-0.5',
#           '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/au-capacitor/siesta1-smeagol/kpoints-2-2-20_V-4_WeightRho-0.5',
#           '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/cp2k-smeagol-examples/examples/au-capacitor/siesta1-smeagol/kpoints-4-4-20_V-4_double-contour']
# folder = ['/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-capacitor/layers-1-2-3-4/siesta/archer/li/1V/kpoints-4-4-20_hlb-auto_cores-20_li_WeightRho-0.5_delta-1e-4_ivan',
#           '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-capacitor/layers-1-2-3-4/siesta/archer/li/4V/kpoints-4-4-20_hlb-auto_cores-20_li_4V_WeightRho-0.5_delta-1e-4_ivan',
#           '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-capacitor/layers-1-2-3-4/siesta/archer/li/4V/weight/kpoints-4-4-20_hlb-auto_cores-20_li_4V_delta-1e-4_weight-3-squared']
# folder = ['/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-capacitor/layers-1-2-3-4/siesta/archer/li/4V/kpoints-4-4-20_hlb-auto_cores-20_li_4V_WeightRho-0_delta-1e-4_ivan',
#           '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-capacitor/layers-1-2-3-4/siesta/archer/li/4V/kpoints-4-4-20_hlb-auto_cores-20_li_4V_WeightRho-0.5_delta-1e-4_ivan',
#           '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-capacitor/layers-1-2-3-4/siesta/archer/li/4V/kpoints-4-4-20_hlb-auto_cores-20_li_4V_WeightRho-1_delta-1e-4_ivan',
#           '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-capacitor/layers-1-2-3-4/siesta/archer/li/4V/kpoints-4-4-20_hlb-auto_cores-20_li_4V_delta-1e-4_chris2-weight',
#           '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-capacitor/layers-1-2-3-4/siesta/archer/li/4V/kpoints-4-4-20_hlb-auto_cores-20_li_4V_delta-1e-4_chris2-weight2']
# folder = ['/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-capacitor/layers-1-2-3-4/siesta/archer/li/kpoints-4-4-20_hlb-auto_cores-20_li_4V_WeightRho-0.5_delta-1e-4-chris',
#           '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-capacitor/layers-1-2-3-4/siesta/archer/li/kpoints-4-4-20_hlb-auto_cores-20_li_4V_WeightRho-0.5_delta-1e-4-alex']
# folder = ['/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/siesta-smeagol/capacitor/testing/other/kpoints-2-2-20_hlb-auto_cores-64_au_TransmissionOverk_EM.TRCChannels_3V',
#           '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/siesta-smeagol/capacitor/testing/other/kpoints-4-4-20_hlb-auto_cores-64_au_TransmissionOverk_EM.TRCChannels_3V',
#           ]
# folder = ['/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/siesta-smeagol/capacitor/testing/other/kpoints-2-2-20_hlb-auto_cores-20_li_4V_WeightRho-0.5',
#           '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/siesta-smeagol/capacitor/testing/other/kpoints-3-3-20_hlb-auto_cores-20_li_4V_WeightRho-0.5',
#           '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/siesta-smeagol/capacitor/testing/other/kpoints-4-4-20_hlb-auto_cores-20_li_4V_WeightRho-0.5',
#           ]
# folder = ['/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/siesta-smeagol/capacitor/testing/archer2/bound_states/li/kpoints-4-4-20_hlb-auto_cores-20_li_4V_WeightRho-0.5_delta-1e-2',
#           '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/siesta-smeagol/capacitor/testing/archer2/bound_states/li/kpoints-4-4-20_hlb-auto_cores-20_li_4V_WeightRho-0.5_delta-1e-4',
#           '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/siesta-smeagol/capacitor/testing/archer2/bound_states/li/kpoints-4-4-20_hlb-auto_cores-20_li_4V_WeightRho-0.5_delta-0',
#           '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/siesta-smeagol/capacitor/testing/archer2/bound_states/li/kpoints-4-4-20_hlb-auto_cores-20_li_4V_WeightRho-0.5_BS.MiddleOrbital-0_BS.Tolerance-1e-5_delta-1e-1',
#           '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/siesta-smeagol/capacitor/testing/archer2/bound_states/li/kpoints-4-4-20_hlb-auto_cores-20_li_4V_WeightRho-0.5_BS.MiddleOrbital-0_BS.Tolerance-1e-5_delta-1e-2',
#           '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/siesta-smeagol/capacitor/testing/archer2/bound_states/li/kpoints-4-4-20_hlb-auto_cores-20_li_4V_WeightRho-0.5_BS.MiddleOrbital-0_BS.Tolerance-1e-5_delta-1e-4',
#           '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/siesta-smeagol/capacitor/testing/archer2/bound_states/li/kpoints-4-4-20_hlb-auto_cores-20_li_4V_WeightRho-0.5_BS.MiddleOrbital-0_BS.Tolerance-1e-5_delta-0',
#           ]
# folder = [
#     '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/siesta-smeagol/capacitor/testing/other/kpoints-4-4-20_hlb-auto_cores-64_au_TransmissionOverk_EM.TRCChannels',
#     ]
# folder = ['/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/siesta-smeagol/capacitor/testing/other/kpoints-4-4-20_hlb-auto_cores-64_au_TransmissionOverk_EM.TRCChannels',
#           '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/siesta-smeagol/capacitor/testing/other/kpoints-4-4-20_hlb-auto_cores-64_au_TransmissionOverk_EM.TRCChannels_4V',
#           ]
# folder = ['/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/siesta-smeagol/capacitor/testing/other/kpoints-4-4-20_hlb-auto_cores-20_li',
#           '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/siesta-smeagol/capacitor/testing/other/kpoints-4-4-20_hlb-auto_cores-20_li_2V',
#           '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/siesta-smeagol/capacitor/testing/other/kpoints-4-4-20_hlb-auto_cores-20_li_3V',
#           '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/siesta-smeagol/capacitor/testing/other/kpoints-4-4-20_hlb-auto_cores-20_li_4V_WeightRho-1',
#           ]
# folder = ['/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/siesta-smeagol/capacitor/testing/other/kpoints-4-4-20_hlb-auto_cores-20_li',
#           '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/siesta-smeagol/capacitor/testing/other/kpoints-4-4-20_hlb-auto_cores-20_li_2V',
#           '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/siesta-smeagol/capacitor/testing/other/kpoints-4-4-20_hlb-auto_cores-20_li_3V',
#           '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/siesta-smeagol/capacitor/testing/other/kpoints-4-4-20_hlb-auto_cores-20_li_4V',
#           '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/siesta-smeagol/capacitor/testing/other/kpoints-4-4-20_hlb-auto_cores-20_li_5V',
#           '/Volumes/ELEMENTS/Storage/Postdoc/Data/Work/Postdoc/Work/calculations/transport/2023/au-bdt/siesta-smeagol/capacitor/testing/other/kpoints-4-4-20_hlb-auto_cores-20_li_10V'
#           ]
folder_save = []
for i in range(len(folder)):
    folder_save.append('{}/output'.format(folder[i]))
print('Saving to', folder_save)

data_bulk_charge = []
data_bulk_hartree = []
data_em_hartree = []
data_dft_charge = []
data_dft_hartree = []
data_em_charge = []
data_cp2k_hartree = []
data_cp2k_charge = []
for i in range(len(folder)):
    data_bulk_charge.append(np.genfromtxt('{}/{}'.format(folder[i], bulk_charge), skip_header=0, skip_footer=0))
    data_bulk_hartree.append(np.genfromtxt('{}/{}'.format(folder[i], bulk_hartree), skip_header=0, skip_footer=0))
    if plot_dft: data_dft_hartree.append(np.genfromtxt('{}/{}'.format(folder[i], dft_hartree), skip_header=0, skip_footer=0))
    if plot_dft: data_dft_charge.append(np.genfromtxt('{}/{}'.format(folder[i], dft_charge), skip_header=0, skip_footer=0))
    data_em_hartree.append(np.genfromtxt('{}/{}'.format(folder[i], em_hartree), skip_header=0, skip_footer=0))
    data_em_charge.append(np.genfromtxt('{}/{}'.format(folder[i], em_charge), skip_header=0, skip_footer=0))


# Hartree potential
# fig_plot_1, ax_plot_1 = plt.subplots()
# ax_plot_1.plot(data_em_hartree[0][:, 0], data_em_hartree[0][:, 1], 'r-', label='HLB 0')
# ax_plot_1.plot(data_em_hartree[1][:, 0], data_em_hartree[1][:, 1], 'g-', label='HLB 1 eV')
# ax_plot_1.legend(frameon=False)
# ax_plot_1.set_xlabel(r'Position / Å')
# ax_plot_1.set_ylabel('Hartree potential / eV')
# fig_plot_1.tight_layout()
# for i in range(len(folder)):
#     fig_plot_1.savefig('{}/hartree_potential.png'.format(folder_save[i]), dpi=param.save_dpi)

# Plot charge density
# fig_plot_2, ax_plot_2 = plt.subplots()
# ax_plot_2.plot(data_em_charge[0][:, 0], data_em_charge[0][:, 1], 'r-', label='HLB 0')
# ax_plot_2.plot(data_em_charge[1][:, 0], data_em_charge[1][:, 1], 'g-', label='HLB 1 eV')
# print(np.sum(np.abs(data_em_charge[0][:, 1])))
# ax_plot_2.legend(frameon=False)
# ax_plot_2.set_xlabel(r'Position / Å')
# ax_plot_2.set_ylabel('Charge density')
# fig_plot_2.tight_layout()
# for i in range(len(folder)):
#     fig_plot_2.savefig('{}/charge_density.png'.format(folder_save[i]), dpi=param.save_dpi)

# Plot charge density difference
# fig_plot_3, ax_plot_3 = plt.subplots()
# ax_plot_3.plot(data_em_charge[0][:, 0], data_em_charge[0][:, 1]-data_em_charge[1][:, 1], 'k-')
# ax_plot_3.set_xlabel(r'Position / Å')
# ax_plot_3.set_ylabel('Charge density difference (0 - 1)')
# fig_plot_3.tight_layout()
# for i in range(len(folder)):
#     fig_plot_3.savefig('{}/charge_density_difference.png'.format(folder_save[i]), dpi=param.save_dpi)

# Plot Hartree potential and charge density
# rows, cols = 2, 1
# fig_plot_both, ax_plot_both = plt.subplots(rows, cols, sharex='col', sharey='row', figsize=(6, 8))
# ax_plot_both[0].plot(data_em_hartree[0][:, 0], data_em_hartree[0][:, 1], 'g-', label='EM HLB 0=F')
# ax_plot_both[0].plot(data_em_hartree[1][:, 0], data_em_hartree[1][:, 1], 'r-', label='EM HLB 1 eV')
# ax_plot_both[0].plot(data_bulk_hartree[0][:, 0], data_bulk_hartree[0][:, 1], 'k-', label='Bulk')
# ax_plot_both[0].legend(frameon=False)
# ax_plot_both[0].set_ylabel('Hartree potential / eV')
# ax_plot_both[1].plot(data_em_charge[0][:, 0], data_em_charge[0][:, 1], 'g-', label='EM HLB 0=F')
# ax_plot_both[1].plot(data_em_charge[1][:, 0], data_em_charge[1][:, 1], 'r-', label='EM HLB 1 eV')
# ax_plot_both[1].plot(data_bulk_charge[0][:, 0], data_bulk_charge[0][:, 1], 'k-', label='Bulk')
# ax_plot_both[1].legend(frameon=False)
# ax_plot_both[1].set_xlabel(r'Position / Å')
# ax_plot_both[1].set_ylabel('Charge density')
# fig_plot_both.tight_layout()
# for i in range(len(folder)):
#     fig_plot_both.savefig('{}/hartree_potential_charge.png'.format(folder_save[i]), dpi=param.save_dpi)

# Plot Hartree potential and charge density
rows, cols = 2, 1
fig_plot_both, ax_plot_both = plt.subplots(rows, cols, sharex='col', sharey='row', figsize=(6, 8))
if plot_dft: ax_plot_both[0].plot(data_dft_hartree[0][:, 0], data_dft_hartree[0][:, 1], 'r-', label=labels[0])
ax_plot_both[0].plot(data_em_hartree[0][:, 0], data_em_hartree[0][:, 1], 'g-', label=labels[1])
ax_plot_both[0].plot(data_bulk_hartree[0][:, 0], data_bulk_hartree[0][:, 1], 'k-', label=labels[2])
# ax_plot_both[0].plot(x+data_bulk_hartree[0][:, 0], data_bulk_hartree[0][:, 1], 'k-')
if plot_markers: ax_plot_both[0].plot(markers, markers*0, 'o', color='orange', fillstyle='none')
ax_plot_both[0].legend(frameon=False)
ax_plot_both[0].set_ylabel('Hartree potential / eV')
if plot_dft: ax_plot_both[1].plot(data_dft_charge[0][:, 0], data_dft_charge[0][:, 1], 'r-', label=labels[0])
ax_plot_both[1].plot(data_em_charge[0][:, 0], data_em_charge[0][:, 1], 'g-', label=labels[1])
ax_plot_both[1].plot(data_bulk_charge[0][:, 0], data_bulk_charge[0][:, 1], 'k-', label=labels[2])
# ax_plot_both[1].plot(x+data_bulk_charge[0][:, 0], data_bulk_charge[0][:, 1], 'k-')
if plot_markers: ax_plot_both[1].plot(markers, markers*0, 'o', color='orange', fillstyle='none')
ax_plot_both[1].legend(frameon=False)
ax_plot_both[1].set_xlabel(r'Position / Å')
ax_plot_both[1].set_ylabel('Charge density')
fig_plot_both.tight_layout()
for i in range(len(folder)):
    fig_plot_both.savefig('{}/hartree_potential_charge.png'.format(folder_save[i]), dpi=param.save_dpi)

# Plot Hartree potential and charge density difference
if plot_diff:
    rows, cols = 2, 1
    fig_plot_both_diff, ax_plot_both_diff = plt.subplots(rows, cols, sharex='col', sharey='row', figsize=(6, 8))
    # ax_plot_both_diff[0].plot(data_em_hartree[0][:, 0], data_em_hartree[1][:, 1]-data_em_hartree[0][:, 1], 'k-', label='EM')
    ax_plot_both_diff[0].plot(data_em_hartree[0][:, 0], data_em_hartree[0][:, 1]-data_dft_hartree[0][:, 1], 'k-', label='EM')
    if plot_markers: ax_plot_both_diff[0].plot(markers, markers*0, 'o', color='orange', fillstyle='none')
    # ax_plot_both_diff[0].legend(frameon=False)
    ax_plot_both_diff[0].set_ylabel('Hartree potential / eV')
    # ax_plot_both_diff[1].plot(data_em_charge[0][:, 0], data_em_charge[1][:, 1]-data_em_charge[0][:, 1], 'k-', label='EM')
    ax_plot_both_diff[1].plot(data_em_charge[0][:, 0], data_em_charge[0][:, 1]-data_dft_charge[0][:, 1], 'k-', label='EM')
    if plot_markers: ax_plot_both_diff[1].plot(markers, markers*0, 'o', color='orange', fillstyle='none')
    ax_plot_both_diff[1].set_xlabel(r'Position / Å')
    ax_plot_both_diff[1].set_ylabel('Charge density')
    fig_plot_both_diff.tight_layout()
    for i in range(len(folder)):
        fig_plot_both_diff.savefig('{}/hartree_potential_charge_diff.png'.format(folder_save[i]), dpi=param.save_dpi)

# Plot Hartree potential and charge density difference
rows, cols = 2, 1
fig_plot_both_diff_multiple, ax_plot_both_diff_multiple = plt.subplots(rows, cols, sharex='col', sharey='row', figsize=(6, 8))
if plot_diff_multiple:
    for i in range(len(folder)):
        ax_plot_both_diff_multiple[0].plot(data_em_hartree[i][:, 0], data_em_hartree[i][:, 1]-data_dft_hartree[i][:, 1], '-', label=diff_label[i])
        ax_plot_both_diff_multiple[1].plot(data_em_charge[i][:, 0], data_em_charge[i][:, 1]-data_dft_charge[i][:, 1], '-', label=diff_label[i])
if plot_markers: ax_plot_both_diff_multiple[0].plot(markers, markers * 0, 'o', color='orange', fillstyle='none')
if plot_markers: ax_plot_both_diff_multiple[1].plot(markers, markers*0, 'o', color='orange', fillstyle='none')
ax_plot_both_diff_multiple[0].legend(frameon=False)
ax_plot_both_diff_multiple[1].set_xlabel(r'Position / Å')
ax_plot_both_diff_multiple[1].set_ylabel('Charge density')
ax_plot_both_diff_multiple[0].set_ylabel('Hartree potential / eV')
fig_plot_both_diff_multiple.tight_layout()
for i in range(len(folder)):
    fig_plot_both_diff_multiple.savefig('{}/hartree_potential_charge_diff_multiple.png'.format(folder_save[i]), dpi=param.save_dpi)
        
print('Value(z=Min) of Hartree potential z EM .DAT', data_em_hartree[0][0, 1])
print('Value(z=Min) of Hartree potential z bulk .DAT', data_bulk_hartree[0][0, 1])

print('Value(z=Max) of Hartree potential z EM .DAT', data_em_hartree[0][-1, 1])
print('Value(z=Max) of Hartree potential z bulk .DAT', data_bulk_hartree[0][-1, 1])

print('Average of z=min, z=max of Hartree potential z EM .DAT', np.average([data_em_hartree[0][0, 1], data_em_hartree[0][-1, 1]]))
print('Average of z=min, z=max of Hartree potential z bulk .DAT', np.average([data_bulk_hartree[0][0, 1], data_bulk_hartree[0][-1, 1]]))

if __name__ == "__main__":
    print('Finished.')
    plt.show()
