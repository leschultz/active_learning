from pymatgen import Element, Composition
import numpy as np
import sys
import re

comp = sys.argv[1]  # Composition
margin = float(sys.argv[2])  # Fractional Margin
n_tests = int(sys.argv[3])  # The number of points to test

comp = re.split('(\d+)', comp)
comp = [i for i in comp if i != '']

n = [int(i) for i in comp if i.isdigit()]
e = [i for i in comp if not i.isdigit()]
r = []

for i in e:
    el = Element(i)
    rad = []

    # Gather applicable radii
    try:
        rad.append(el.atomic_radius)
    except Exception:
        pass
    try:
        rad.append(el.average_anionic_radius)
    except Exception:
        pass
    try:
        rad.append(el.average_cationic_radius)
    except Exception:
        pass
    try:
        rad.append(el.average_ionic_radius)
    except Exception:
        pass
    try:
        rad.append(el.metallic_radius)
    except Exception:
        pass

    r.append(max(rad))  # Maximum radius

r = np.array(r)
v = 4/3*np.pi*r**3

vol = 0.0
for i, j in zip(n, v):
    vol += i*j

vol_max = (1+margin)*vol

l = vol**(1/3)
l_max = vol_max**(1/3)

print('Test Lengths [Angstroms]: '+str(np.linspace(l, l_max, n_tests)))
