from pymatgen import Lattice, Structure
from ast import literal_eval
from shutil import copyfile
import pymatgen as mg
import pandas as pd
import subprocess
import functions
import sys
import os

sys.argv[1:] = [literal_eval(i) for i in sys.argv[1:]]  # Convert types

# Input parameters
contcar = 'CONTCAR'  # The CONTCAR file
kpoints = 'KPOINTS'  # The KPOINTS file
potcar = 'POTCAR'  # The POTCAR file
incar = '../../../run_creator/templates/incar/hold'  # The VASP input file
fits = '../../../../volumes/tv_curves/data/data.csv'  # Data for linear fits
temp = sys.argv[1]  # The temperature

fits = pd.read_csv(fits)  # Load TV curves

dir_name = os.path.join('../', str(temp))

# Read the composition from contcar and find the volume
count = 0
with open(contcar) as f:
    for line in f:
        if count == 5:
            elements = line.strip().split(' ')
            elements = [i for i in elements if i != '']
        if count == 6:
            numbers = line.strip().split(' ')
            numbers = [i for i in numbers if i != '']
        count += 1

composition = ''
for i, j in zip(elements, numbers):
    composition += i+j

fit = fits.loc[fits['composition'] == composition].values[0]
m = fit[1]
b = fit[2]

# Change volume
structure = Structure.from_file(contcar)
structure.lattice = Lattice.cubic((m*temp+b)**(1/3))

# Open and read template
incar = open(incar)
incar_contents = incar.read()
incar.close()

# Write INCAR file
incar_contents = incar_contents.replace('{temp}', str(temp))
file_out = open(os.path.join(dir_name, 'INCAR'), 'w')
file_out.write(incar_contents)
file_out.close()

# Save input files
structure.to(fmt='poscar', filename=os.path.join(dir_name, 'POSCAR'))
copyfile(kpoints, os.path.join(dir_name, 'KPOINTS'))
copyfile(potcar, os.path.join(dir_name, 'POTCAR'))
copyfile(submit, os.path.join(dir_name, 'parallel.sh'))
