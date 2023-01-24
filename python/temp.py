from __future__ import division, print_function
import pandas as pd
import numpy as np
import glob
from general import load_coordinates
from general import print_xyz
import matplotlib.pyplot as plt

"""
    Temp
"""

# Ni_1 = 2
# Pt_1 = 0
# Pt_2 = 0
# Ni_3 = 2
# Pt_3 = 0
# Pt_4 = 0

Ni_1 = 1
Pt_1 = 0.5  # 4/8
Pt_2 = 0.25  # 2/8
Ni_3 = 1
Pt_3 = 0.5  # 4/8
Pt_4 = 0.25  # 2/8

multiplicity = (2*0.5*Ni_1*4+2*0.5*Pt_1*4+2*0.5*Pt_2*8+2*0.5*Pt_3*4+2*0.5*Ni_3*4+2*0.5*Pt_4*8)+1
print(multiplicity)
