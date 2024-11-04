#!/home/yangyudi/software/anaconda3/bin/python
import numpy as np
import re
import matplotlib.pyplot as plt


##############################################################################
# This program is used to plot the band structures only when spin = 1.
# This program needs to read the input and output of the band calculations in
# QE. Please remember to change the name of the the input and output files.
##############################################################################

def readfile(filename):
    '''
    This function is used to find 6 groups of data:
    lattice: lattice parameter
    index: number of k point
    valence: number of valence bands
    k_matrix: 3*index matrix to represent the coordinates of k points
    e_matrix: energy value for different k points.
    '''

    file = open(filename, 'r')
    read = file.read()  # Save a file as a string.
    file.close()
    files = open(filename, 'r')  # Save a file as a file

    search1 = 'lattice parameter'  # Find the lattice parameter
    search2 = 'number of Kohn-Sham states'  # Find the number of bands
    search3 = 'number of k points'
    search4 = 'number of electrons'
    search5 = 'reciprocal axes'
    search6 = 'k(    1) = ('  # Find the value of energy for corresponding k point (reciprocal)
    search7 = 'number of atoms/cell'
    search8 = 'number of atomic types'

    p = re.compile(
        r"k[\w ]*=([\- ]*\d\.\d{4})([\- ]*\d\.\d{4})([\- ]*\d\.\d{4})([\w\d\(\)\s\:]*)([\d\s\.\-]*)",
        re.M)
    matches = p.findall(read)  # Find all the data for the energy with the corresponding k points
    b_matrix = np.zeros((3, 3))  # The matrix of the

    # initialize the count numbers
    b_count = 0  # count the elenent in the b matrix
    reciprocal = 0  # count the line of the b matrix
    cart = 0  # start to count the line of k points matrix in cart. coordinates

    for line in files:
        if search7 in line:
            N_atom = choosefloat(line)[0]  # number of atoms
        if search8 in line:
            N_type = choosefloat(line)[0]  # number of types of atoms
        if search1 in line:
            lattice = choosefloat(line)[0]  # lattice parameter
        if search2 in line:
            index = int(choosefloat(line)[0])  # number of bands
        if search4 in line:
            electron = int(choosefloat(line)[0])  # for spin=1, occupation = N of electron/2
            occupation = electron / 2
        if search3 in line:
            numberk = int(choosefloat(line)[0])  # number of k points
            cryst_matrix = np.zeros((numberk, 3))  # The matrix to record the position of k points in reciprocal coord
            cart_matrix = np.zeros((numberk, 3))  # The matrix to record the position of k points in cart. coord
            e_matrix = np.zeros((index, numberk))  # create the matrix of energy for each k-point
        if search6 in line:
            cart += 1
        if cart > 1 and cart < numberk + 2:
            k_line = choosefloat(line)
            cryst_matrix[cart - 2][0] = k_line[1]
            cryst_matrix[cart - 2][1] = k_line[2]
            cryst_matrix[cart - 2][2] = k_line[3]
            cart += 1
        if search5 in line:
            reciprocal += 1  # start to read the k points matrix for reciprocal coordinates
        elif reciprocal > 0 and reciprocal < 4:
            b_number = choosefloat(line)
            b_matrix[b_count][0] = b_number[1]
            b_matrix[b_count][1] = b_number[2]
            b_matrix[b_count][2] = b_number[3]
            b_count += 1
            reciprocal += 1
    # print(str(matches))
    for i in range(numberk):

        the_string = str(matches[i])

        the_float = choosefloat(the_string)
        cart_matrix[i][0] = the_float[0]
        cart_matrix[i][1] = the_float[1]
        cart_matrix[i][2] = the_float[2]
        for j in range(index):
            num = j + 4
            e_matrix[j][i] = the_float[num]

    return lattice, index, occupation, cryst_matrix, e_matrix, b_matrix, cart_matrix, N_atom, N_type


def choosefloat(string):
    '''
    This function could find all numbers in a given string.
    '''

    matches = re.findall("[+-]?\d+\.{0,1}\d*", string)

    ary = np.zeros(len(matches))
    for i in range(len(matches)):
        ary[i] = float(matches[i])

    return ary


def fermisurface(energy, valence):
    '''
    This function is used to adjust the 0 point on the y-axis aline to the
    fermi energy.
    '''

    valence_bands = energy[int(valence) - 1]
    valence_max = max(valence_bands)
    new_energy = energy - valence_max

    return valence_max, new_energy


def specialk(input_file):
    '''
    This function is used to find the indexs of the special k points and
    the corresponding names.
    '''

    file = open(input_file, 'r')

    count = 0  # Find the line for K_POINTS
    line_n = 0  # Line for specialk_mareix
    search = 'K_POINTS {crystal_b}'  # Find the data of special k points
    number_kpt = 0  # count the total number of k points
    for line in file:
        if search in line:
            count += 1
        elif count == 1:
            new_line = line.split(' ')
            k_total = int(choosefloat(new_line[0]))

            count += 1
            specialk_matrix = np.zeros((k_total, 3))
            specialk_name = []
        elif count == 2:
            if line_n < k_total:
                nameline = line.split(' ')
                points = choosefloat(line)
                specialk_matrix[line_n][0] = points[0]
                specialk_matrix[line_n][1] = points[1]
                specialk_matrix[line_n][2] = points[2]
                number_kpt += int(points[3])
                if int(points[3]) == 0:
                    number_kpt += 1
                if nameline[-1] == '\n':
                    specialk_name.append(nameline[-2] + nameline[-1])
                else:
                    specialk_name.append(nameline[-1])
                line_n += 1

    return specialk_matrix, specialk_name, number_kpt


def renameLatex(specialk_name):
    '''
    This function is used to assign the label on x-axis.
    '''

    name = []
    for i in range(len(specialk_name)):
        if specialk_name[i].replace('\n', '') == 'Gamma':
            name.append(r'$\Gamma$')
        else:
            new_name = '$' + specialk_name[i].replace('\n', '') + '$'
            name.append(new_name)

    return name


def find(kpt_matrix, special_kpt):
    '''
    This fuction is used to find the index of the special k points.
    '''

    k_total = len(special_kpt)
    specialk_index = np.zeros(k_total)
    count_spk = 0
    delta = 1 / 10000
    for i in range(len(kpt_matrix)):
        if abs(kpt_matrix[i][0] - special_kpt[count_spk][0]) < delta:
            if abs(kpt_matrix[i][1] - special_kpt[count_spk][1]) < delta:
                if abs(kpt_matrix[i][2] - special_kpt[count_spk][2]) < delta:
                    specialk_index[count_spk] = i
                    count_spk += 1

    return specialk_index


def rescale(b_matrix, specialk_matrix, kpt_matrix, specialk_index):
    '''
    This function is used to rescale the x-axis. Let the distance between
    different special k points be different.
    '''

    new_array = np.zeros(len(kpt_matrix))
    distance1 = 0
    for i in range(len(specialk_matrix) - 1):
        new_x = specialk_matrix[i + 1][0] * b_matrix[0][0] + specialk_matrix[i + 1][1] * b_matrix[1][0] + \
                specialk_matrix[i + 1][2] * b_matrix[2][0]
        new_y = specialk_matrix[i + 1][0] * b_matrix[0][1] + specialk_matrix[i + 1][1] * b_matrix[1][1] + \
                specialk_matrix[i + 1][2] * b_matrix[2][1]
        new_z = specialk_matrix[i + 1][0] * b_matrix[0][2] + specialk_matrix[i + 1][1] * b_matrix[1][2] + \
                specialk_matrix[i + 1][2] * b_matrix[2][2]
        old_x = specialk_matrix[i][0] * b_matrix[0][0] + specialk_matrix[i][1] * b_matrix[1][0] + specialk_matrix[i][
            2] * b_matrix[2][0]
        old_y = specialk_matrix[i][0] * b_matrix[0][1] + specialk_matrix[i][1] * b_matrix[1][1] + specialk_matrix[i][
            2] * b_matrix[2][1]
        old_z = specialk_matrix[i][0] * b_matrix[0][2] + specialk_matrix[i][1] * b_matrix[1][2] + specialk_matrix[i][
            2] * b_matrix[2][2]
        diff = np.sqrt((new_x - old_x) ** 2 + (new_y - old_y) ** 2 + (new_z - old_z) ** 2)
        distance2 = distance1 + diff
        index1 = int(specialk_index[i])
        index2 = int(specialk_index[i + 1])
        positions = abs(index2 - index1)
        new_array[index1: index2] = np.linspace(distance1, distance2, positions) * 10
        distance1 = distance2

    new_array[-1] = new_array[-2] + (new_array[-2] - new_array[-3])

    return new_array


def writing(filename, valence, energy, coord, fermi, k_index, k_name):
    '''
    This function is used to write the data into a datafile.
    '''

    file = open(filename, 'x')
    file.write('YYD is so fucking smart!!! \n')
    file.write('YYD is so fucking smart!!! \n')
    file.write('YYD is so fucking smart!!! \n\n')
    # Write in the Fermi energy
    file.write('Fermi energy = ' + str(fermi) + '\n\n')

    # Write the number of band for VBM and CBM
    file.write('CBM : ' + str(int(valence)) + '\n')
    file.write('VBM : ' + str(int(valence) - 1) + '\n\n')

    # Write in the Energy data and special k points info
    # Energy data
    for i in range(len(energy)):
        # file.write('(' + str(coord[i][0]) + ',' + str(coord[i][1]) + ',' + str(coord[i][2]) + ') \n')
        file.write(str(i + 1) + '  ')
        for j in range(len(energy[0])):
            file.write(str(energy[i][j]) + ' ')
        file.write('\n\n')

    file.write('\n\n Special k points \n\n')
    nameinx = 0
    for n in k_index:
        file.write(k_name[nameinx] + '\n')
        for m in range(len(energy)):
            file.write(str(energy[m][int(n)]) + ' ')
        nameinx += 1
        file.write('\n\n')

    file.write('\n\n Special k points \n\n')
    for k in range(2):
        if k == 0:
            file.write('--- Spin Up --- \n\n')
        if k == 1:
            file.write('--- Spin Down --- \n\n')
        nameinx = 0
        for n in k_index:
            file.write(k_name[nameinx] + '\n')
            for m in range(len(energy[0])):
                file.write(str(energy[k][m][int(n)]) + ' ')
            nameinx += 1
            file.write('\n\n')


def plot_band():
    '''
    Plotting part.
    '''

    length = len(cart_matrix)
    zero = np.ones(len(new_xaxis)) * valence_max
    zero_0 = np.zeros(len(new_xaxis))
    zero_fermi = np.ones(len(new_xaxis)) * energy_fermi
    fig = plt.figure()

    label_index = []
    label = []

    for i in range(len(namelist)):
        index = int(special_index[i])
        if i > 0 and i < len(namelist) - 1:
            ind_diff = index - int(special_index[i - 1])
            if ind_diff < 2:
                posi_diff = new_xaxis[index] - new_xaxis[index - 1]
                new_xaxis[index:] = new_xaxis[index:] - posi_diff
                newname = namelist[i - 1] + '/' + namelist[i]
                namelist[i] = newname
                del label[-1]
                del label_index[-1]
        label.append(namelist[i])
        label_index.append(new_xaxis[index])

    sn = 0
    for i in special_index:
        if int(i - sn) != 1:
            plt.axvline(new_xaxis[int(i)], linestyle='--', color='blue', alpha=0.2)
        sn = i
    # print(np.shape(energy_matrix[0]))
    VBM = max(energy_matrix[vbm_idx])
    CBM = min(energy_matrix[cbm_idx])
    Gap = CBM - VBM

    for i in range(len(energy_matrix)):
        plt.plot(new_xaxis, energy_matrix[i]-VBM, zorder=1)
        if i == cbm_idx:
            plt.plot(new_xaxis, energy_matrix[i]-VBM, zorder=1)
        if i == vbm_idx:
            plt.plot(new_xaxis, energy_matrix[i]-VBM, zorder=1)

    # print(len(zero_fermi))
    # plt.plot(new_xaxis, zero_fermi*0, linestyle = '--', color = 'gray')

    plt.axhspan(0, Gap, alpha=0.2)
    print(VBM)
    print(CBM)

    plt.xticks(label_index, label)
    # plt.legend()

    # plt.ylim(min(energy_matrix[vbm_idx-1]-1), max(energy_matrix[cbm_idx]+2))
    plt.xlim(new_xaxis[0], new_xaxis[-1])
    plt.ylim(-2, 6)
    # plt.show()


###############################################################################

params = {'axes.formatter.limits': [-4, 4],
          'axes.labelsize': 'x-large',
          'xtick.labelsize' : '10',
          'axes.titlesize': 'x-large',
          'legend.fontsize': 'medium',
          'lines.markersize': '8',
          }
plt.rcParams.update(params)

# Setting part
# folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/bulk/monoclinic/calculations/band_structure/setyawan-kpoints/yudi-pbe-qe/band'
# folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/bulk/monoclinic/calculations/band_structure/setyawan-kpoints/mp-352-pbe-qe/band'

# folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/bulk/tetragonal/calculations/band_structure/setyawan-kpoints/yudi-pbe-qe/band'
# folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/bulk/tetragonal/calculations/band_structure/setyawan-kpoints/mp-1018721-qe/band'

folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/bulk/po/calculations/band_structure/setyawan-kpoints/yudi-pbe-qe/band'
folder = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hfo2/archer/bulk/po/calculations/band_structure/setyawan-kpoints/mp-685097-qe/band'

input_filename = "input.inp"
output_filename = "qe_log.log"
# input_filename = "input"
# output_filename = "output"

input_file = '{}/{}'.format(folder, input_filename)
output_file = '{}/{}'.format(folder, output_filename)

lattice_parameter, N_bands, valence_num, cryst_matrix, energy_matrix, b_matrix, cart_matrix, N_atom, N_type = readfile(output_file)
valence_max, energy_fermi = fermisurface(energy_matrix, valence_num)
special_kpt, special_name, N_kpt = specialk(input_file)
special_index = find(cryst_matrix, special_kpt)
namelist = renameLatex(special_name)
new_xaxis = rescale(b_matrix, special_kpt, cart_matrix, special_index)

# Minus the Fermi energy, if do not need, comment this line
# energy_matrix = energy_fermi
energy_fermi = 14.0466
energy_matrix = energy_matrix - energy_fermi

# BOHR2ANG = 0.529177
# factor = 2 * np.pi / (lattice_parameter * BOHR2ANG)
# kpoints_matrix_factor = cryst_matrix*factor

# These two lines will be labled as red
vbm_idx = int(valence_num - 1)
# chg_idx = int(valence_num-1)
cbm_idx = int(valence_num)

plot_band()
CBM = min(energy_matrix[cbm_idx])
charge = max(energy_matrix[cbm_idx - 1])
VBM = max(energy_matrix[vbm_idx])
Gap = CBM - VBM
Core = min(energy_matrix[0])
print('CBM = ', CBM)
print('Charge state = ', charge)
print('VBM = ', VBM)
print('Band Gap = ', Gap)
print('Core =', Core)
plt.ylabel(r'E - E$_\mathrm{f}$ (eV)')

# plt.title("PBE+U")
# plt.ylim(-10, 4)
plt.tight_layout()
plt.savefig('{}/bands.png'.format(folder), dpi=500)
plt.show()
