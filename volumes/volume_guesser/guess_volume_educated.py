from pymatgen import Element, Composition
import numpy as np
import sys
import re

l = float(sys.argv[1])  # Starting cube side length
margin = float(sys.argv[2])  # Fractional Margin
n_tests = int(sys.argv[3])  # The number of points to test

vol = l**3
vol_min = (1-margin)*vol
vol_max = (1+margin)*vol

l_min = vol_min**(1/3)
l_max = vol_max**(1/3)

print('Test Lengths [Angstroms]: '+str(np.linspace(l_min, l_max, n_tests)))
