from ast import literal_eval
from shutil import copyfile
import subprocess
import functions
import sys
import os

sys.argv[1:] = [literal_eval(i) for i in sys.argv[1:]]  # Convert types

# Input parameters
contcar = 'CONTCAR'  # The CONTCAR file
kpoints = 'KPOINTS'  # The KPOINTS file
potcar = 'POTCAR'  # The POTCAR file
incar = (
         '/home/leschultz/work'
         '/vasp_runs/volumes/run_creator'
         '/cooling_holds/templates/incar/hold'
         )  # The VASP input file
submit = (
         '/home/leschultz/work'
         '/vasp_runs/volumes/run_creator'
         '/cooling_holds'
         '/templates/submit/bardeen_morgan.q'
         )  # The cluster submit file
fits = '../../../tv_curves/data/data.csv'  # Data for linear fits
dT = sys.argv[1]  # Change in temperature
min_temp = sys.argv[2]  # The minimum allowable temperature hold

cwd = os.getcwd()
temp = float(os.path.basename(cwd))-dT

if temp >= min_temp:

    dir_name = os.path.join('../', str(temp))
    functions.create_dir(dir_name)

    # Open and read template
    incar = open(incar)
    incar_contents = incar.read()
    incar.close()

    # Write INCAR file
    incar_contents = incar_contents.replace('{temp}', str(temp))
    file_out = open(os.path.join(dir_name, 'INCAR'), 'w')
    file_out.write(incar_contents)
    file_out.close()

    copyfile(contcar, os.path.join(dir_name, 'POSCAR'))  # Save POSCAR
    copyfile(kpoints, os.path.join(dir_name, 'KPOINTS'))  # Save KPOINTS
    copyfile(potcar, os.path.join(dir_name, 'POTCAR'))  # Save POSCAR
    copyfile(submit, os.path.join(dir_name, 'parallel.sh'))  # Save KPOINTS

    os.chdir(dir_name)  # Change working directory

    # Submit new hold
    subprocess.run(
                   'qsub parallel.sh', cwd=dir_name,
                   stdout=open(os.devnull, 'wb')
                   )
