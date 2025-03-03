import pandas as pd
import numpy as np
import glob
from matplotlib import pyplot as plt
from general import parameters as param
from matplotlib.colors import LogNorm
import dpdata
import cp2kdata
from cp2kdata import Cp2kOutput

file_in = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/chris_phd/bulk/hole/hse/400k/hematite-charges-1-clean.hirshfeld'
file_out = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/chris_phd/bulk/hole/hse/400k-f/hematite-charges-1-clean.hirshfeld'
file_out = '/Volumes/ELEMENTS/Storage/Postdoc2/Data/Work/calculations/hematite/chris_phd/bulk/hole/hse/400k-b/hematite-charges-1-clean.hirshfeld'

# line_start = 2
# line_end = 3  # inclusive

# step number starts at 0, therefore line_start = 1 is equal to step 0
# step number to line number is +1
line_start = 1699 + 1
line_end = 2899 + 1  # inclusive

line_start = 3300 + 1
line_end = 6000 + 1  # inclusive

with open(file_in, 'r') as input_file:

    with open(file_out, 'w') as output_file:

        lines = input_file.readlines()

        count = 0

        for i, line in enumerate(lines):

            if 'Hirshfeld Charges' in line:

                count += 1

                if count >= line_start and count <= line_end:

                    output_file.write(line)

                elif count == line_end+1:

                    break

            elif count >= line_start and count <= line_end:

                output_file.write(line)
