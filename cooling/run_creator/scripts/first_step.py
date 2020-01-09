from shutil import copyfile
from os.path import join
import pandas as pd
import numpy as np
import functions
import re

# Input paramters
coords = '../templates/poscar/365.txt'  # Starting coordinates
save = '../../runs'  # Location to save POSCAR
incar = '../templates/incar/hold'  # The VASP input file
kpoints = '../templates/kpoints/M'  # The VASP kpoints file
potcar = '/home/leschultz/work/POTCARs/paw/LDA/5.4'  # The VASP potential file
submit = '../templates/submit/bardeen_morgan.q'  # The cluster submit file
fits = '../../../volumes/tv_curves/data/data.csv'  # Data for linear fits
start_temp = 2000.0  # The starting temperature

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

start_temp_str = str(start_temp)

coords = np.loadtxt(coords)  # Load starting coordinates
coords = coords/coords.max(axis=0)  # Make fractional

# Load the linear fits
fits = pd.read_csv(fits)

# Open and read template
incar = open(incar)
incar_contents = incar.read()
incar.close()

groups = fits.groupby(['composition'])

count = 1
total = str(fits.shape[0])
for group, values in groups:

    i = re.split('(\d+)', group)
    i = [j for j in i if j != '']

    for j in range(len(i)):
        try:
            i[j] = int(i[j])
        except Exception:
            pass

    numbers = [j for j in i if isinstance(j, int)]
    elements = [j for j in i if isinstance(j, str)]

    new_coords = coords[:sum(numbers), :]  # The first n atoms
    np.random.shuffle(new_coords)  # Randomize coordinates

    m = values['slope'].values[0]
    b = values['intercept'].values[0]

    length = (m*start_temp+b)**(1/3)

    run = join(*[save, group, start_temp_str])

    # Generates structure
    structure = functions.gen_cubic(
                                    elements,
                                    numbers,
                                    length,
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

        pots.append(join(*[potcar, i, 'POTCAR']))

    # Write run
    functions.create_dir(run)

    # Concat POTCARs and save
    with open(join(run, 'POTCAR'), 'wb') as outfile:
        for i in pots:
            with open(i, 'rb') as potfile:
                outfile.write(potfile.read())

    # Write INCAR file
    incar = incar_contents
    incar = incar.replace('{temp}', start_temp_str)
    file_out = open(join(run, 'INCAR'), 'w')
    file_out.write(incar)
    file_out.close()

    copyfile(kpoints, join(run, 'KPOINTS'))  # Save KPOINTS
    copyfile(submit, join(run, 'parallel.sh'))  # Save KPOINTS

    # Save POSCAR
    structure.to(fmt='poscar', filename=join(run, 'POSCAR'))

    # Status update
    print('Generated ('+str(count)+'/'+total+')'+run)

    count += 1
