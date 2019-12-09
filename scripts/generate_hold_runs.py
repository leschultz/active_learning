from ast import literal_eval
from shutil import copyfile

import numpy as np
import functions
import sys
import os

sys.argv[1:] = [literal_eval(i) for i in sys.argv[1:]]  # Convert types

# Input paramters
a = sys.argv[1]  # Cubic lattice constant
elements = sys.argv[2]  # Elements
numbers = sys.argv[3]  # Number corresponding to each element
coords = sys.argv[4]  # Starting coordinates
save = sys.argv[5]  # Location to save POSCAR
runs = sys.argv[6]  # The number of runs to generate
incar = sys.argv[7]  # The VASP input file
kpoints = sys.argv[8]  # The VASP kpoints file
potcar = sys.argv[9]  # The VASP potential file
submit = sys.argv[10]  # The cluster submit file

coords = np.loadtxt(coords)  # Load starting coordinates
coords = coords/coords.max(axis=0)  # Make fractional
np.random.shuffle(coords)  # Randomize coordinates

# Generates structure
structure = functions.gen_cubic(
                                elements,
                                numbers,
                                coords,
                                a,
                                )

# Random seeds
numbers = list(range(runs))

# Generate runs
count = 1
total = str(runs)
for number in numbers:
    run = os.path.join(save, str(number))

    # Find matching POTCARS
    pots = []
    for i in elements:

        # Missing Zr potential
        if i == 'Zr':
            i = 'Zr_sv'

        pots.append(os.path.join(*[potcar, i, 'POTCAR']))

    # Write run
    functions.create_dir(run)

    # Concat POTCARs and save
    with open(os.path.join(run, 'POTCAR'), 'wb') as outfile:
        for i in pots:
            with open(i, 'rb') as potfile:
                outfile.write(potfile.read())

    copyfile(incar, os.path.join(run, 'INCAR'))  # Save INCAR
    copyfile(kpoints, os.path.join(run, 'KPOINTS'))  # Save KPOINTS
    copyfile(submit, os.path.join(run, 'paralle.sh'))  # Save KPOINTS

    # Save POSCAR
    structure.to(fmt='poscar', filename=os.path.join(run, 'POSCAR'))

    # Status update
    print('Generated ('+str(count)+'/'+total+')'+run)

    count += 1
