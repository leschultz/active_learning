from ast import literal_eval
from shutil import copyfile

import numpy as np
import functions
import sys
import os

# Recomended potentials
recomended_sv = {
                 'Zr',
                 'Li',
                 'Ca',
                 'Nb',
                 'K',
                 'Sc',
                 'Ti',
                 'V',
                 'Rb',
                 'Sr',
                 'Y',
                 'Mo',
                 'Cs',
                 'Ba',
                 'Fr',
                 'Ra'
                 }
recomended_pv = {
                 'Na',
                 'Cr',
                 'Tc',
                 'Ru',
                 'Mn',
                 'Rh',
                 'Hf',
                 'Ta',
                 'W'
                 }
recomended_d = {
                'Sn',
                'Ga',
                'Ge',
                'In',
                'Tl',
                'Pb',
                'Bi',
                'Po',
                'At'
                }

sys.argv[1:] = [literal_eval(i) for i in sys.argv[1:]]  # Convert types

# Input paramters
elements = sys.argv[1]  # Elements
numbers = sys.argv[2]  # Number corresponding to each element
save = sys.argv[3]  # Location to save POSCAR
runs = sys.argv[4]  # List of cubic lattice constants
incar = sys.argv[5]  # The VASP input file
kpoints = sys.argv[6]  # The VASP kpoints file
potcar = sys.argv[7]  # The VASP potential file
submit = sys.argv[8]  # The cluster submit file

# Generate runs
count = 1
total = str(len(runs))
for number in runs:
    run = os.path.join(save, str(number))

    # Generates structure
    structure = functions.gen_cubic(
                                    elements,
                                    numbers,
                                    number,
                                    )

    # Find matching POTCARS
    pots = []
    for i in elements:

        if i in recomended_sv:
            i += '_sv'

        if i in recomended_pv:
            i += '_pv'

        if i in recomended_d:
            i += '_d'

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
    copyfile(submit, os.path.join(run, 'parallel.sh'))  # Save KPOINTS

    # Save POSCAR
    structure.to(fmt='poscar', filename=os.path.join(run, 'POSCAR'))

    # Status update
    print('Generated ('+str(count)+'/'+total+')'+run)

    count += 1
